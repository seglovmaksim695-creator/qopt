"""Этап 4, проба А: T с измерениями (логическое И Гидни, arXiv:1709.06648).
Автоматически находит в схеме .qc:
  R1 «вычисление»: Тоффоли пишет в кубит, про который известно, что он |0>  -> логическое И, 4 T (вместо 7);
  R2 «отмена»:     Тоффоли стирает кубит, про который известно, что он = a·b (a, b не менялись) -> H, измерение,
                   при исходе 1 — CZ(a,b), сброс; 0 T.
Кубиты не из строки .i — вспомогательные, стартуют в |0>.
Проверка: симуляция векторов состояния с измерениями (случайные входы на данных, вспомогательные |0>,
случайные исходы по их вероятностям) против исходной схемы, допуск 1e-10. До 24 кубитов, иначе «не проверено».
python gidney.py bench_qc [bench_new ...]"""
import glob, json, os, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "portfolio"))
import qoptct, evo_lib

# логическое И Гидни (проверено на матрицах: |a,b,0> -> |a,b,ab> с общей фазой)
def and_gate(a, b, t):
    """Логическое И: |a,b,0> -> |a,b,ab> точно, 4 T и 3 CX (rccx Марголуса + S† на цели; находка инженера 3, 11 нояб.,
    перепроверено на матрицах). Прежний вариант Гидни — 6 CX: and_gate6."""
    return evo_lib.rccx(a, b, t) + [("sdg", (t,))]


def and_gate6(a, b, t):
    return [("h", (t,)), ("t", (t,)), ("cx", (a, t)), ("cx", (b, t)), ("cx", (t, a)), ("cx", (t, b)),
            ("tdg", (a,)), ("tdg", (b,)), ("t", (t,)), ("cx", (t, a)), ("cx", (t, b)), ("h", (t,)), ("s", (t,))]

DIAG = {"z", "s", "sdg", "t", "tdg", "cz", "ccz"}


def inputs_of(path):
    for line in open(path):
        if line.startswith(".i"):
            return line.split()[1:]
    return None


def transform(hl, n, anc):
    """-> список: обычные команды уровня Тоффоли, ('AND',(a,b,t)), ('MZ',(a,b,t))."""
    st = {q: ("zero" if q in anc else None) for q in range(n)}
    ver = [0] * n; out = []; r1 = r2 = 0
    for g, q in hl:
        q = tuple(q)
        if g == "ccx":
            a, b, t = q; s = st[t]
            if s == "zero":
                out.append(("AND", q)); st[t] = ("and", frozenset((a, b)), {a: ver[a], b: ver[b]}); r1 += 1; continue
            if isinstance(s, tuple) and s[1] == frozenset((a, b)) and s[2][a] == ver[a] and s[2][b] == ver[b]:
                out.append(("MZ", q)); st[t] = "zero"; r2 += 1; continue
            tg = (t,)
        elif g in DIAG:
            out.append((g, q)); continue
        elif g in ("cx",):
            tg = (q[1],)
        else:  # h, x и прочее — меняет значение кубита
            tg = q
        out.append((g, q))
        for x in tg:
            ver[x] += 1; st[x] = None
    return out, r1, r2


def anf_check(items, n, anc, cap=20000):
    """Точное доказательство для схем без H: значения кубитов — многочлены над GF(2) от входов.
    Требуем: в AND цель = 0, в MZ цель = a·b. Тогда замена точна (фазы — диагональные команды)."""
    if any(g == "h" for g, _ in items): return None
    v = [set() if q in anc else {frozenset([q])} for q in range(n)]
    def mul(p, r):
        out = set()
        for x in p:
            for y in r:
                out ^= {x | y}
        return out
    for g, q in items:
        if g == "x": v[q[0]] ^= {frozenset()}
        elif g == "cx": v[q[1]] = v[q[1]] ^ v[q[0]]
        elif g == "ccx": v[q[2]] = v[q[2]] ^ mul(v[q[0]], v[q[1]])
        elif g == "AND":
            if v[q[2]]: return False
            v[q[2]] = mul(v[q[0]], v[q[1]])
        elif g == "MZ":
            if v[q[2]] != mul(v[q[0]], v[q[1]]): return False
            v[q[2]] = set()
        elif g in DIAG: continue
        else: return None
        if g in ("cx", "ccx", "AND") and len(v[q[-1]]) > cap: return None
    return True


def segments(items):
    """Разбить на унитарные куски между измерениями; вернуть [(низкоуровневый кусок, MZ или None)]."""
    segs, cur = [], []
    for g, q in items:
        if g == "AND":
            cur += and_gate(*q)
        elif g == "MZ":
            cur.append(("h", (q[2],))); segs.append((cur, q)); cur = []
        else:
            cur += evo_lib.lower_gate(g, q)
    segs.append((cur, None))
    return segs


# --- симулятор векторов состояния ---
_H = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
_1Q = {"h": _H, "x": np.array([[0, 1], [1, 0]], complex), "z": np.diag([1, -1]).astype(complex),
       "s": np.diag([1, 1j]), "sdg": np.diag([1, -1j]),
       "t": np.diag([1, np.exp(1j * np.pi / 4)]), "tdg": np.diag([1, np.exp(-1j * np.pi / 4)])}


def apply(psi, n, g, q):
    psi = psi.reshape([2] * n)
    if g in _1Q:
        psi = np.moveaxis(np.tensordot(_1Q[g], psi, axes=([1], [q[0]])), 0, q[0])
    elif g in ("cx", "ccx"):
        *c, t = q; idx = [slice(None)] * n
        for x in c: idx[x] = 1
        sub = psi[tuple(idx)]; ax = t - sum(1 for x in c if x < t)
        psi = psi.copy(); psi[tuple(idx)] = np.flip(sub, axis=ax)
    elif g in ("cz", "ccz"):
        idx = [slice(None)] * n
        for x in q: idx[x] = 1
        psi = psi.copy(); psi[tuple(idx)] *= -1
    else:
        raise ValueError(g)
    return psi.reshape(-1)


def run_plain(hl, n, psi):
    for g, q in hl: psi = apply(psi, n, g, tuple(q))
    return psi


def run_meas(segs, n, psi, rng):
    for seg, mz in segs:
        for g, q in seg: psi = apply(psi, n, g, q)
        if mz:
            a, b, t = mz; p = psi.reshape([2] * n)
            idx1 = [slice(None)] * n; idx1[t] = 1
            p1 = float(np.sum(np.abs(p[tuple(idx1)]) ** 2))
            m = int(rng.random() < p1)
            idx0 = [slice(None)] * n; idx0[t] = 1 - m
            p = p.copy(); p[tuple(idx0)] = 0; psi = p.reshape(-1) / np.sqrt(p1 if m else 1 - p1)
            if m:
                psi = apply(psi, n, "cz", (a, b)); psi = apply(psi, n, "x", (t,))  # поправка и сброс
    return psi


def verify(hl, segs, n, anc, trials=4, seed=0):
    if n > 24: return None
    rng = np.random.default_rng(seed); data = [q for q in range(n) if q not in anc]
    for _ in range(trials):
        v = rng.normal(size=2 ** len(data)) + 1j * rng.normal(size=2 ** len(data)); v /= np.linalg.norm(v)
        psi = np.zeros([2] * n, complex)
        idx = [slice(None) if q in data else 0 for q in range(n)]
        psi[tuple(idx)] = v.reshape([2] * len(data)) if data else v[0]
        psi = psi.reshape(-1)
        ref = run_plain(hl, n, psi.copy()); got = run_meas(segs, n, psi.copy(), rng)
        if abs(abs(np.vdot(ref, got)) - 1) > 1e-10: return False
    return True


def tc(c): return sum(g in ("t", "tdg") for g, _ in c)


if __name__ == "__main__":
    out = {}
    for folder in sys.argv[1:]:
        for f in sorted(glob.glob(str(ROOT / folder / "*.qc"))):
            nm = os.path.basename(f)[:-3]
            names, hl0 = qoptct.load_qc(f); n = len(names)
            ins = inputs_of(f); anc = set() if ins is None else {i for i, x in enumerate(names) if x not in ins}
            hl = evo_lib.fold_h(hl0)
            items, r1, r2 = transform(hl, n, anc)
            segs = segments(items)
            low = [(evo_lib.fold_cancel(s, n), mz) for s, mz in segs]
            T = sum(tc(s) for s, _ in low); CX = sum(evo_lib.cx(s) for s, _ in low)
            ok = anf_check(items, n, anc) if (r1 or r2) else True
            how = "anf" if ok else None
            if ok is None:
                try:
                    ok = qoptct.limited(600, verify, hl0, low, n, anc); how = "sim"
                except qoptct._TO:
                    ok = None
            base = evo_lib.fold_cancel(evo_lib.lower_std(hl), n)
            r = dict(n=n, anc=len(anc), and_compute=r1, meas_uncompute=r2, T_meas=T, CX_meas=CX, cond_cz=r2,
                     T_nomeas_fold=tc(base), CX_nomeas_fold=evo_lib.cx(base), verified=ok, how=how)
            out[nm] = r; print(nm, json.dumps(r), flush=True)
        (ROOT / "results" / "gidney.json").write_text(json.dumps(out, indent=1))
