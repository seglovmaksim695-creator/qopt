"""c10: НЕЗАВИСИМАЯ проверка больших детерминированных схем Clifford+T с измерениями (30-70 кубитов), без вектора состояний.

Метод: sum-over-paths (Amy) + точное исключение переменных (bucket elimination на тензорах, ширина ~10-15).
Для входного базисного состояния x и ВЫБРАННЫХ исходов измерений m (ветвь) состояние данных равно
   2^{-nH/2} * sum_y w^{F(y)} |W y + e>       (см. pathsum.py)
Сравнивается с эталоном (схема из .qc, lower_std, тоже путь-интеграл): в каждой ветви должно быть
   A_det(d) * 2^{nm/2} = c_m * A_ref(d)  для всех d, |c_m| = 1, c_m не зависит от x  (=> U_branch = c_m U).
Условия XOR/AND в if поддерживаются (pathsum.parse_cond).
Режимы: 'enum' (n<=8: все d) и 'classical' (эталонный выход d* считается классически, подтверждается путь-интегралом
эталона |A_ref(d*)|=1; плюс проверка нулевых амплитуд в нескольких других точках образа).
Ветви выбираются равномерно (не по весам) — для корректных детерминированных протоколов исходы равновероятны;
нулевая ветвь помечается отдельно (zero), а не считается ошибкой молча.
"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, math, random, cmath, time
import numpy as np
sys.path[:0] = [_QPKG, (_QPKG + '/portfolio'), (_QPKG + '/contrib/c10')]
import pathsum as P

OM = [cmath.exp(1j * math.pi * k / 4) for k in range(8)]


def _bits(m):
    r = []; i = 0
    while m:
        if m & 1: r.append(i)
        m >>= 1; i += 1
    return r


def _par_table(w):
    """чётность по w битам -> массив формы (2,)*w"""
    idx = np.indices((2,) * w).reshape(w, -1).sum(0) & 1
    return idx.reshape((2,) * w)


_PT = {}


def par(w):
    if w not in _PT: _PT[w] = _par_table(w)
    return _PT[w]


def amplitude_sum(ps, targets, maxw=24):
    """S = sum_{y: W y + e = targets} w^{F(y)}  (без нормировки 2^{-nH/2}). targets: {wire: bit} для живых проводов."""
    scal = OM[ps.gph % 8]
    factors = []   # (scope list, table)
    for m, c in ps.lin.items():
        vs = _bits(m)
        assert len(vs) <= maxw, "слишком тяжёлая чётность %d" % len(vs)
        t = np.where(par(len(vs)) == 1, OM[c % 8], 1.0 + 0j)
        factors.append((vs, t))
    for (ma, mb) in ps.quad:
        va, vb = _bits(ma), _bits(mb); vs = sorted(set(va) | set(vb))
        assert len(vs) <= maxw
        pos = {v: i for i, v in enumerate(vs)}
        ia = np.zeros(2 ** len(vs), dtype=np.int64)
        grid = np.indices((2,) * len(vs)).reshape(len(vs), -1)
        pa = grid[[pos[v] for v in va]].sum(0) & 1; pb = grid[[pos[v] for v in vb]].sum(0) & 1
        factors.append((vs, np.where((pa & pb) == 1, -1.0 + 0j, 1.0 + 0j).reshape((2,) * len(vs))))
    if not ps.consistent: return 0j
    for m, bit in ps.cons:
        vs = _bits(m); assert len(vs) <= maxw
        factors.append((vs, (par(len(vs)) == bit).astype(complex)))
    for q, bit in targets.items():
        m, c = ps.w[q]
        if m == 0:
            if c != bit: return 0j
            continue
        vs = _bits(m); assert len(vs) <= maxw
        want = bit ^ c
        factors.append((vs, (par(len(vs)) == want).astype(complex)))
    # исключение переменных
    var_f = {}
    alive = set()
    fl = {}
    for i, (vs, t) in enumerate(factors):
        fl[i] = (vs, t)
        for v in vs: var_f.setdefault(v, set()).add(i); alive.add(v)
    nxt = len(factors)
    free_vars = set(range(ps.k)) - alive   # переменные без единого множителя: сумма = 2
    scal *= 2 ** len(free_vars)
    while alive:
        best = None
        for v in alive:
            sc = set()
            for i in var_f[v]: sc |= set(fl[i][0])
            key = (len(sc), len(var_f[v]))
            if best is None or key < best[0]: best = (key, v, sc)
        _, v, sc = best
        fi = sorted(var_f[v]); scope = sorted(sc)
        lab = {u: j for j, u in enumerate(scope)}
        out = [lab[u] for u in scope if u != v]
        args = []
        for i in fi:
            args += [fl[i][1], [lab[u] for u in fl[i][0]]]
        r = np.einsum(*args, out, optimize=False) if len(scope) <= 50 else None
        newscope = [u for u in scope if u != v]
        for i in fi:
            for u in fl[i][0]: var_f[u].discard(i)
            del fl[i]
        alive.discard(v); del var_f[v]
        if newscope:
            fl[nxt] = (newscope, r)
            for u in newscope: var_f.setdefault(u, set()).add(nxt)
            nxt += 1
        else:
            scal *= complex(r)
        # переменные, потерявшие все множители, уже вне alive (были в scope => получили новый множитель)
        alive = {u for u in alive if var_f.get(u)}
    return scal


# ------------------------------------------------------------ эталон
def ref_ps(gates, n, x):
    ps = P.PS(n, x)
    for g, q in gates:
        if g == "h": ps.hadamard(q[0])
        else: ps.gate(g, q)
    return ps


def classical_target(hl, n, x):
    """Эталонный выход для классических схем (h ccz h = ccx). Бросает ValueError, если схема не такая."""
    b = [x.get(q, 0) for q in range(n)]
    frame = set()
    for g, q in hl:
        if g == "h":
            frame ^= {q[0]}
            continue
        fr = [u for u in q if u in frame]
        if g == "x":
            if fr: continue   # X в H-рамке = Z: на классические биты не влияет
            b[q[0]] ^= 1
        elif g == "cx":
            if fr: raise ValueError("cx на H-рамке")
            if b[q[0]]: b[q[1]] ^= 1
        elif g in ("z", "cz", "ccz"):
            if not fr: continue
            if len(fr) != 1: raise ValueError("две H-рамки в одном CCZ")
            t = fr[0]; ctl = [u for u in q if u != t]
            if all(b[u] for u in ctl): b[t] ^= 1
        else:
            raise ValueError("не классический вентиль " + g)
    if frame: raise ValueError("незакрытая H-рамка")
    return b


def image_points(ps, wires, cap, rng):
    """Точки образа аффинного отображения y -> (W y + e) на данных 0..n-1: список целых d (бит q = провод q)."""
    cols = {}
    for j, q in enumerate(wires):
        for v in _bits(ps.w[q][0]): cols[v] = cols.get(v, 0) | (1 << j)
    basis = []
    for c in cols.values():
        for b in basis: c = min(c, c ^ b)
        if c: basis.append(c)
    e = sum(ps.w[q][1] << j for j, q in enumerate(wires))
    r = len(basis)
    if 2 ** r <= cap:
        pts = [e]
        for b in basis: pts = pts + [p ^ b for p in pts]
        return pts, True
    pts = {e}
    while len(pts) < cap:
        p = e
        for b in basis:
            if rng.random() < 0.5: p ^= b
        pts.add(p)
    return sorted(pts), False


def check(qasm, qc_path, x_trials=4, m_trials=4, seed=1, verbose=True, cap=2048, tol=1e-7, max_zero_tries=40):
    """Главная проверка. Возвращает dict(ok, F, worst_dev, ...). F = 1 - наибольшее отклонение (грубая оценка, см. REPORT).
    Для каждой пары (x, ветвь m): амплитуды детерминированной схемы на ВСЕХ точках образа карты данных (2^rank точек,
    либо случайная подвыборка cap), эталон — на точках носителя; вектор ветви должен быть c_m * (эталонный выход),
    c_m не зависит от x (в т.ч. по модулю: вероятность ветви не должна зависеть от входа)."""
    import qoptct, evo_lib
    names, hl = qoptct.load_qc(qc_path); n = len(names); ref = evo_lib.lower_std(hl)
    prog = P.parse(qasm); N, nb, ops = prog
    rng = random.Random(seed)
    worst = 0.0; zero_branches = 0; cs = {}; t0 = time.time(); complete = True; maxsupp = 0
    xs = [{q: rng.randint(0, 1) for q in range(n)} for _ in range(x_trials)]
    ms = []
    while len(ms) < m_trials: ms.append([rng.randint(0, 1) for _ in range(nb)])
    rps = [ref_ps(ref, n, x) for x in xs]
    refcache = [dict() for _ in xs]

    def aref(xi, d):
        if d not in refcache[xi]:
            rp = rps[xi]
            refcache[xi][d] = amplitude_sum(rp, {q: (d >> q) & 1 for q in range(n)}) * 2 ** (-rp.nH / 2)
        return refcache[xi][d]
    mi = 0; tries = 0
    pending = list(ms)
    done_m = []
    while pending:
        m = pending.pop(0)
        vecs = []; branch_zero = False
        for xi, x in enumerate(xs):
            try: ps, c = P.build(prog, x, outcomes=m)
            except ValueError as e:
                return dict(ok=False, F=0.0, error=str(e))
            al = [q for q in range(N) if ps.alive[q]]
            assert al[:n] == list(range(n)), "данные 0..n-1 должны остаться живыми"
            scale = 2 ** (-ps.nH / 2) * 2 ** (ps.nmeas / 2)
            pts, full = image_points(ps, al, cap, rng); complete &= full
            amp = {}
            for d in pts:
                tg = {q: (d >> j) & 1 for j, q in enumerate(al)}
                amp[d] = amplitude_sum(ps, tg) * scale
            mx = max(abs(v) for v in amp.values())
            if mx < 1e-9:
                branch_zero = True; break
            supp = [d for d, v in amp.items() if abs(v) > 1e-9 * max(1.0, mx)]
            maxsupp = max(maxsupp, len(supp))
            # эталон на носителе: ||ref|| на носителе = 1
            # строгость: живые анциллы обязаны завершить в одном базисном состоянии (не запутаны с данными)
            hi = {d >> n for d in supp}
            if len(hi) > 1: worst = max(worst, 1.0)
            h0 = min(hi)
            ar = {d: (aref(xi, d & ((1 << n) - 1)) if (d >> n) == h0 else 0j) for d in supp}
            nrm = sum(abs(v) ** 2 for v in ar.values())
            worst = max(worst, abs(nrm - 1))
            d0 = max(supp, key=lambda d: abs(ar[d]))
            if abs(ar[d0]) < 1e-9:
                worst = max(worst, 1.0); vecs.append(complex(1.0)); continue
            cc = amp[d0] / ar[d0]
            dev = max(abs(amp[d] - cc * ar[d]) / abs(cc) for d in supp)
            worst = max(worst, dev)
            vecs.append(cc)
        if branch_zero:
            zero_branches += 1
            if zero_branches > max_zero_tries * m_trials: return dict(ok=False, F=0.0, error="слишком много нулевых ветвей")
            pending.append([rng.randint(0, 1) for _ in range(nb)]); continue
        for cc in vecs:
            worst = max(worst, abs(cc - vecs[0]) / abs(vecs[0]))
        done_m.append(abs(vecs[0]) ** 2 * 2 ** (-ps.nmeas))
        if verbose: print("  ветвь done t=%.1fs worst=%.2e |c|^2*2^-nb=%.3g supp<=%d" % (time.time() - t0, worst, done_m[-1], maxsupp), flush=True)
    return dict(ok=bool(worst < tol), F=max(0.0, 1 - worst), worst_dev=worst, zero_branches=zero_branches,
                complete_image=complete, branch_prob=done_m, n=n, N=N, measurements=nb, time=time.time() - t0)


if __name__ == "__main__":
    qasm, qc = sys.argv[1], sys.argv[2]
    xt = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    mt = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    print(check(qasm, qc, x_trials=xt, m_trials=mt))
