"""Совместный синтез фазовых полиномов блоков с отложенным линейным состоянием (lowering для gf2 и др.).
CCX -> H·CCZ·H; схема режется на блоки {cx,x,диаг.} с учётом коммутации (H на q задерживает только команды на q):
у gf2 все старшие CCZ (общие a_i, b_j, разные c_k) попадают в ОДИН блок. Фазы блока складываются по чётностям,
сеть CX строится жадно (каждый CX делает одну чётность готовой; штраф за уход от логического состояния) или
GraySynth, K вариантов на блок, оценка = сеть + починка. Линейное состояние не сбрасывается: перед H на q
чинятся только строка/столбец q (жадно + исключение), полный сброс — в конце. Чётности можно переносить между
блоками, где они выразимы (режимы assign)."""
NAME = "phasepoly_graysynth"; KIND = "lower"; HL_RAW = True
import random
from evo_lib import fold_cancel, key, HL_DIAG

SEQ = {0: [], 1: ["t"], 2: ["s"], 3: ["s", "t"], 4: ["s", "s"], 5: ["sdg", "tdg"], 6: ["sdg"], 7: ["tdg"]}
PHG = {"t": 1, "tdg": 7, "s": 2, "sdg": 6, "z": 4}
def comm(g, h):
    """Команды уровня {h,x,cx,диаг} коммутируют (достаточное условие)."""
    (n1, q1), (n2, q2) = g, h
    if not set(q1) & set(q2): return True
    if n1 == "h" or n2 == "h": return False
    d1, d2 = n1 in HL_DIAG, n2 in HL_DIAG
    if d1 and d2: return True
    if n1 == "x" and n2 == "x": return True
    if n1 == "cx" and n2 == "cx": return q1[0] != q2[1] and q2[0] != q1[1]
    if (d1 and n2 == "x") or (d2 and n1 == "x"): return False
    if n1 == "cx": n1, q1, n2, q2 = n2, q2, n1, q1
    if n1 == "x": return q1[0] == q2[1]
    return q2[1] not in q1  # диагональная и cx: цель cx не задета
def layers(g):
    """Чередование: блок (не-H, вытащенные вперёд) / слой H. Возврат список ('B', гейты) | ('H', кубиты)."""
    R = list(g); out = []
    while R:
        for want_h in (False, True):
            rest, took = [], []
            for x in R:
                if (x[0] == "h") == want_h and all(comm(x, y) for y in rest):
                    took.append(x)
                else:
                    rest.append(x)
            R = rest
            if took: out.append(("H", [x[1][0] for x in took]) if want_h else ("B", took))
    return out
def add(P, m, v):
    P[m] = (P.get(m, 0) + v) % 8
    if not P[m]: del P[m]
def block_poly(gates, W, K, P):
    """Прогнать блок по аффинному состоянию проводов (W — маски чётностей, K — константы), набрать фазы P."""
    def term(qs, v):
        m = 0; k = 0
        for q in qs: m ^= W[q]; k ^= K[q]
        add(P, m, -v if k else v)
    for nm, qs in gates:
        if nm == "cx":
            a, b = qs; W[b] ^= W[a]; K[b] ^= K[a]
        elif nm == "x": K[qs[0]] ^= 1
        elif nm in PHG: term(qs, PHG[nm])
        elif nm == "cz":
            a, b = qs; term((a,), 2); term((b,), 2); term((a, b), -2)
        elif nm in ("cs", "csdg"):
            w = 1 if nm == "cs" else -1; a, b = qs
            term((a,), w); term((b,), w); term((a, b), -w)
        elif nm == "ccz":
            a, b, c = qs
            for s, v in (((a,), 1), ((b,), 1), ((c,), 1), ((a, b), -1), ((a, c), -1), ((b, c), -1), ((a, b, c), 1)):
                term(s, v)
        else:
            raise ValueError(nm)
def solve(W, n, m):
    """Выразить маску m через строки W: вернуть маску строк или None."""
    basis = []; res = 0
    for i in range(n):
        v, c = W[i], 1 << i
        for bv, bc in basis:
            if v ^ bv < v: v ^= bv; c ^= bc
        if v: basis.append((v, c)); basis.sort(reverse=True)
    for bv, bc in basis:
        if m ^ bv < m: m ^= bv; res ^= bc
    return None if m else res
def graysynth(P, W, n, rng):
    """Синтез фаз P (маски над переменными) на проводах с состоянием W (меняется). Возврат схема."""
    out = []; cols = {m: solve(W, n, m) for m in P}; vals = dict(P)
    def emit_ready(ids):
        for m in list(ids):
            c = cols[m]
            if c and c & (c - 1) == 0:
                q = c.bit_length() - 1
                out.extend((s, (q,)) for s in SEQ[vals[m]])
                del cols[m]
    def cnot(c, t):
        out.append(("cx", (c, t))); W[t] ^= W[c]
        for m in cols:
            if cols[m] >> t & 1: cols[m] ^= 1 << c
        emit_ready(list(cols))
    emit_ready(list(cols))
    stack = [(set(cols), frozenset(range(n)), None)]
    while stack:
        S, I, i = stack.pop()
        S = {m for m in S if m in cols}
        if not S or not I: continue
        if i is not None:
            while True:
                S = {m for m in S if m in cols}
                if not S: break
                js = [j for j in range(n) if j != i and all(cols[m] >> j & 1 for m in S)]
                if not js: break
                cnot(rng.choice(js) if len(js) > 1 else js[0], i)
            S = {m for m in S if m in cols}
            if not S: continue
        best = -1; cand = []
        for j in I:
            o = sum(cols[m] >> j & 1 for m in S); sc = max(o, len(S) - o)
            if sc > best: best, cand = sc, [j]
            elif sc == best: cand.append(j)
        j = rng.choice(cand) if len(cand) > 1 else cand[0]
        S0 = {m for m in S if not cols[m] >> j & 1}; S1 = S - S0; I2 = I - {j}
        if i is None:
            stack.append((S1, I2, j)); stack.append((S0, I2, None))
        else:
            stack.append((S1, I2, i)); stack.append((S0, I2, i))
    return None if cols else out
def greedy(P, W, n, rng, Lend, lam, mu=0.0, noise=0.5):
    """Жадная сеть чётностей: каждый шаг — CX, делающий одну чётность готовой (popcount 2 -> 1), с минимумом
    Σ|столбцов| + lam·(удаление от логического состояния конца блока); иначе — шаг, уменьшающий Σ|столбцов|."""
    out = []; cols = {}; vals = dict(P)
    for m in P: cols[m] = solve(W, n, m)
    M = [solve(Lend, n, W[r]) for r in range(n)]
    pc = lambda x: bin(x).count("1")
    def emit_ready():
        for m in [m for m, c in cols.items() if c & (c - 1) == 0]:
            q = cols[m].bit_length() - 1
            out.extend((s, (q,)) for s in SEQ[vals[m]]); del cols[m]
    def cnot(c, t):
        out.append(("cx", (c, t))); W[t] ^= W[c]; M[t] ^= M[c]
        for m in cols:
            if cols[m] >> t & 1: cols[m] ^= 1 << c
        emit_ready()
    def delta(c, t):
        d = 0; d2 = 0
        for x in cols.values():
            if x >> t & 1:
                p = pc(x)
                if x >> c & 1: d -= 1; d2 += (p == 3) - (p == 2)
                else: d += 1; d2 += (p == 1) - (p == 2)
        return d + lam * (pc(M[t] ^ M[c]) - pc(M[t])) - mu * d2

    emit_ready(); lim = 6 * len(P) + 6 * n + 50
    while cols:
        lim -= 1
        if lim < 0: return None  # защита от зацикливания: вариант отбрасывается
        ops = set()
        lo = min(pc(y) for y in cols.values())
        src_ = [x for x in cols.values() if pc(x) <= max(lo, 2)]
        for x in src_:
            bits = [k for k in range(n) if x >> k & 1]
            for a in bits:
                for b in bits:
                    if a != b: ops.add((a, b))
        sc = [(delta(c, t) + rng.random() * noise, c, t) for c, t in ops]
        _, c, t = min(sc)
        cnot(c, t)
    return out
def lin_fix(M, n, Q, rng):
    """M — координаты физ. проводов в логическом базисе; CX до M[q]=e_q (q∈Q) и бита q нет в других строках."""
    out = []
    def cnot(c, t):
        out.append((c, t)); M[t] ^= M[c]
    Q = list(Q); rng.shuffle(Q); fixed = set()
    for q in Q:
        b = 1 << q
        if M[q] & b == 0:
            r = min((k for k in range(n) if M[k] & b), key=lambda k: (k in Q, bin(M[k]).count("1"), rng.random()))
            cnot(r, q)
        for r in range(n):
            if r != q and M[r] & b: cnot(q, r)
        extra = M[q] & ~b
        if extra:
            c = solve([M[k] if k != q and k not in fixed else 0 for k in range(n)], n, extra)
            for k in range(n):
                if c >> k & 1: cnot(k, q)
        fixed.add(q)
    return out
def lin_greedy(M, n, Q, rng):
    """Жадно: CX, сильнее всего уменьшающий нарушение (строки/столбцы Q) + 0.3·вес остальных строк."""
    out = []; Qs = set(Q); qm = sum(1 << q for q in Qs); pc = lambda x: bin(x).count("1")
    def cost(r, v):
        if r in Qs: return pc(v ^ (1 << r)) * 2
        return pc(v & qm) * 2 + 0.3 * (pc(v) - 1) * (len(Qs) == n)
    while True:
        best = 0; ops = []
        for t in range(n):
            ct = cost(t, M[t])
            for c in range(n):
                if c == t: continue
                d = cost(t, M[t] ^ M[c]) - ct
                if d < best - 1e-9: best, ops = d, [(c, t)]
                elif d < 0 and abs(d - best) < 1e-9: ops.append((c, t))
        if not ops: break
        c, t = rng.choice(ops); out.append((c, t)); M[t] ^= M[c]
    return out
def lin_best(M, n, Q, rng, tries=8):
    best = None
    for k in range(tries):
        M2 = list(M)
        ops = (lin_greedy(M2, n, Q, rng) if k % 2 else []) + lin_fix(M2, n, Q, rng)
        if best is None or len(ops) < len(best): best = ops
    for c, t in best: M[t] ^= M[c]
    return best
def prepare(hl, n):
    """Проход 1: слои, фазы каждого блока (маски над переменными) и базисы L на блоках."""
    g = []
    for nm, qs in hl:
        qs = tuple(qs)
        if nm == "ccx": g += [("h", (qs[2],)), ("ccz", qs), ("h", (qs[2],))]
        else: g.append((nm, qs))
    L = [1 << q for q in range(n)]; KL = [0] * n; nv = n; lay = []
    for kind, x in layers(g):
        if kind == "B":
            Lb = list(L); P = {}; block_poly(x, L, KL, P); lay.append(("B", P, Lb, list(L), list(KL)))
        else:
            lay.append(("H", x))
            for q in x: L[q] = 1 << nv; nv += 1
    return lay
def assign(lay, n, mode, rng):
    """Переносит каждую чётность в блок, где она выразима (span не меняется внутри блока): собирает
    общие чётности (a_i^b_j и т.п.) в один блок для совместного синтеза."""
    blocks = [k for k, x in enumerate(lay) if x[0] == "B"]
    if mode == 0 or len(blocks) < 2: return [lay[k][1] for k in blocks]
    terms = {}
    for bi, k in enumerate(blocks):
        for m, v in lay[k][1].items():
            terms.setdefault(m, []).append((bi, v))
    can = {}
    for m in terms:
        can[m] = [bi for bi, k in enumerate(blocks) if solve(lay[k][2], n, m) is not None]
    newP = [dict() for _ in blocks]; left = set(terms)
    while left:
        cnt = {}
        for m in left:
            w = 1 if mode == 1 else bin(m).count("1")
            for bi in can[m]: cnt[bi] = cnt.get(bi, 0) + w
        top = max(cnt.values()); bs = [b for b in cnt if cnt[b] == top]
        b = rng.choice(bs)
        for m in [m for m in left if b in can[m]]:
            for _, v in terms[m]: add(newP[b], m, v)
            left.discard(m)
    return newP
def synth(lay, Ps, n, seed, K=20):
    """Поблочно: K вариантов сети (жадная/GraySynth) + починка перед следующим слоем H, лучший по CX."""
    rng = random.Random(seed)
    L = [1 << q for q in range(n)]; KL = [0] * n; W = list(L); KW = [0] * n; nv = n; out = []
    def fixops(W2, Lt, Q):
        M = [solve(Lt, n, W2[r]) for r in range(n)]
        return lin_best(M, n, Q, rng, 4)
    def emit_fix(ops, Q, KLt):
        for c, t in ops:
            out.append(("cx", (c, t))); W[t] ^= W[c]; KW[t] ^= KW[c]
        for q in Q:
            if KW[q] != KLt[q]: out.append(("x", (q,))); KW[q] ^= 1

    bi = 0
    for k, x in enumerate(lay):
        if x[0] == "B":
            P = Ps[bi]; bi += 1
            Q = set(lay[k + 1][1]) if k + 1 < len(lay) else set(range(n))
            best = None
            for v in range(K if P else 1):
                W2 = list(W)
                if v == K - 1:
                    try: net = graysynth(P, W2, n, rng)
                    except Exception: net = None
                elif v < 5: net = greedy(P, W2, n, rng, x[3], (0.0, 0.3, 0.6, 1.0, 1.5)[v % 5])
                else: net = greedy(P, W2, n, rng, x[3], rng.uniform(0, 1.5), rng.uniform(0, 1.5), rng.uniform(0.1, 1.0))
                if net is None: continue
                fo = fixops(W2, x[3], Q)
                sc = sum(g[0] == "cx" for g in net) + len(fo)
                if best is None or sc < best[0]: best = (sc, net, fo)
            _, net, fo = best
            for nm, qs in net:
                if nm == "cx":
                    W[qs[1]] ^= W[qs[0]]; KW[qs[1]] ^= KW[qs[0]]; out.append((nm, qs))
                else:
                    q = qs[0]
                    out.append((INVP[nm], qs) if KW[q] else (nm, qs))
            L[:] = x[3]; KL[:] = x[4]
            emit_fix(fo, Q, KL)
        else:
            Q = set(x[1])
            emit_fix(fixops(W, L, Q), Q, KL)
            for q in x[1]:
                out.append(("h", (q,))); L[q] = W[q] = 1 << nv; nv += 1
    emit_fix(fixops(W, L, range(n)), range(n), KL)
    return out
INVP = {"t": "tdg", "tdg": "t", "s": "sdg", "sdg": "s"}
def run(hl, n):
    lay = prepare(hl, n); best = None
    for seed in range(6):
        rng = random.Random(1000 + seed)
        Ps = assign(lay, n, (1, 2, 0)[seed % 3], rng)
        c = fold_cancel(synth(lay, Ps, n, seed), n)
        if best is None or key(c) < key(best): best = c
    return best
