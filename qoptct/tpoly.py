"""Свой T-этап (TODD-подобный, Heyfron–Campbell arXiv:1712.01557; см. также arXiv:2407.08695).
Безадамаровы куски схемы -> фазовый полином (маска чётности -> вес mod 8) + аффинная линейная часть.
Нечётные столбцы P (T-count = их число). Шаг: для пары столбцов a,b берём z = p_a^p_b и ищем y в ядре системы
  [все единицы; строки P; для троек i<j<k: z_i(p_j∘p_k)+z_j(p_i∘p_k)+z_k(p_i∘p_j)] с y_a != y_b;
тогда P' = P + z·y^T имеет тот же кубический тензор, столбцы a и b совпадают -> T-count −2.
Разность фаз (Δ) — клиффордова; её вычисляем точно по коэффициентам полинома (≤3 переменных) и добавляем чётными весами.
Синтез обратно — graysynth (c0005), линейная часть — исключением Гаусса. Результат всей схемы доказывается qoptct.prove.
python tpoly.py схема.qc ..."""
import importlib.util, itertools, json, random, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "portfolio"))
import qoptct, evo_lib

_s = importlib.util.spec_from_file_location("pp", ROOT / "portfolio/programs/c0005_phasepoly_graysynth.py")
pp = importlib.util.module_from_spec(_s); _s.loader.exec_module(pp)
PH = {"t", "tdg", "s", "sdg", "z", "cx", "x"}


def chunks(c):
    """Разбить на максимальные куски из CX/X/фаз; прочее (H) — отдельно."""
    out, cur = [], []
    for g in c:
        if g[0] in PH: cur.append(g)
        else:
            if cur: out.append(("poly", cur)); cur = []
            out.append(("raw", [g]))
    if cur: out.append(("poly", cur))
    return out


def coeffs(terms):
    """Коэффициенты многочлена (мономы ≤3 переменных) mod 8 для суммы весов на чётностях."""
    C = {}
    for m, w in terms:
        bits = [i for i in range(m.bit_length()) if m >> i & 1]
        for r in (1, 2, 3):
            f = w * (-2) ** (r - 1)
            for S in itertools.combinations(bits, r):
                C[S] = (C.get(S, 0) + f) % 8
    return {k: v for k, v in C.items() if v}


def kernel(rows, m):
    """Базис ядра системы (строки — битовые маски по m столбцам)."""
    piv = {}
    for r in rows:
        for p, pr in piv.items():
            if r >> p & 1: r ^= pr
        if r:
            p = r.bit_length() - 1
            for q in list(piv):
                if piv[q] >> p & 1: piv[q] ^= r
            piv[p] = r
    free = [j for j in range(m) if j not in piv]
    basis = []
    for f in free:
        y = 1 << f
        for p, pr in piv.items():
            if pr >> f & 1: y |= 1 << p
        basis.append(y)
    return basis


def todd(cols, n, rng, max_iter=10000):
    """cols — список масок нечётных столбцов (вес 1). Возврат новый список (меньше или тот же)."""
    cols = list(cols)
    for _ in range(max_iter):
        m = len(cols)
        if m < 2: break
        rowsP = [sum((cols[j] >> i & 1) << j for j in range(m)) for i in range(n)]
        pairs = [(a, b) for a in range(m) for b in range(a + 1, m)]; rng.shuffle(pairs)
        done = False
        tried = set()
        for a, b in pairs:
            z = cols[a] ^ cols[b]
            if z in tried: continue
            tried.add(z)
            zb = [i for i in range(n) if z >> i & 1]
            rows = [(1 << m) - 1] + rowsP
            act = [i for i in range(n) if rowsP[i]]
            for i, j, k in itertools.combinations(act, 3):
                r = 0
                if z >> i & 1: r ^= rowsP[j] & rowsP[k]
                if z >> j & 1: r ^= rowsP[i] & rowsP[k]
                if z >> k & 1: r ^= rowsP[i] & rowsP[j]
                if r: rows.append(r)
            for y in kernel(rows, m):
                if (y >> a & 1) != (y >> b & 1):
                    new = [cols[j] ^ (z if y >> j & 1 else 0) for j in range(m)]
                    cnt = {}
                    for c in new:
                        if c: cnt[c] = cnt.get(c, 0) ^ 1
                    new2 = [c for c, v in cnt.items() if v]
                    if len(new2) < m:
                        cols = new2; done = True
                    break
            if done: break
        if not done: break
    return cols


def lin_ops(W, n):
    """CX-операции (c,t): W[t]^=W[c], приводящие W к единичной."""
    W = list(W); ops = []
    for i in range(n):
        if not W[i] >> i & 1:
            k = next(k for k in range(i + 1, n) if W[k] >> i & 1)
            W[i] ^= W[k]; ops.append((k, i))
        for k in range(n):
            if k != i and W[k] >> i & 1: W[k] ^= W[i]; ops.append((i, k))
    return ops


def resynth(Pold, Wend, Kend, n, rng):
    odd = [m for m, w in Pold.items() if w % 2]
    even = {m: (w - 1) % 8 for m, w in Pold.items() if w % 2}
    for m, w in Pold.items():
        if not w % 2: even[m] = w
    newcols = todd(odd, n, rng)
    if len(newcols) >= len(odd): return None
    old_terms = [(m, 1) for m in odd]; new_terms = [(m, 1) for m in newcols]
    Co, Cn = coeffs(old_terms), coeffs(new_terms)
    D = {k: (Co.get(k, 0) - Cn.get(k, 0)) % 8 for k in set(Co) | set(Cn)}
    P = {}
    for m in newcols: pp.add(P, m, 1)
    for m, w in even.items():
        if w: pp.add(P, m, w)
    for S, v in D.items():
        if not v: continue
        if len(S) == 3: assert v % 8 == 0, "кубическая часть не совпала"; continue
        if len(S) == 1: assert v % 2 == 0; pp.add(P, 1 << S[0], v)
        else:
            assert v % 4 == 0; a, b = S; h = v // 2
            pp.add(P, 1 << a, h); pp.add(P, 1 << b, h); pp.add(P, (1 << a) | (1 << b), -h)
    W = [1 << i for i in range(n)]
    out = pp.graysynth(P, W, n, rng)
    if out is None: return None
    ops = lin_ops(W, n) + list(reversed(lin_ops(Wend, n)))
    out += [("cx", o) for o in ops] + [("x", (q,)) for q in range(n) if Kend[q]]
    return out


def tstage(c, n, seed=0, tries=4):
    res = []
    for kind, g in chunks(c):
        if kind == "raw": res += g; continue
        W = [1 << i for i in range(n)]; K = [0] * n; P = {}
        pp.block_poly(g, W, K, P)
        best = g
        if sum(w % 2 for w in P.values()) >= 4:
            for t in range(tries):
                r = resynth(P, W, K, n, random.Random(seed + t))
                if r is not None and (evo_lib.tcount(r), evo_lib.cx(r)) < (evo_lib.tcount(best), evo_lib.cx(best)): best = r
        res += best
    return evo_lib.fold_cancel(res, n)


if __name__ == "__main__":
    for p in sys.argv[1:]:
        names, ref, best, how = qoptct.optimize(p); n = len(names)
        c2 = tstage(best, n)
        ok = qoptct.prove(c2, ref, n)
        print(Path(p).stem, json.dumps(dict(cx_before=evo_lib.cx(best), t_before=evo_lib.tcount(best),
                                            cx_after=evo_lib.cx(c2), t_after=evo_lib.tcount(c2), proven=ok)), flush=True)
