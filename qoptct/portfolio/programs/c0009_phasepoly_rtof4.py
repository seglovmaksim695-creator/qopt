"""Скрещивание phasepoly_graysynth + rtof4. rtof4: если уровень Тоффоли классический и управляющие всех CCX — аффинные
от входа (gf2), каждый CCX -> rccx4^{±1} (4 CX, D=CS^{∓1}(a,b)), а все поправки переносятся в начало (U — перестановка
с фазами) и синтезируются ОДНОЙ сетью: цепочка разностей a1, a1^a2, ..., каждая b_j проходит её (m^2+3m-2 CX).
Итого ~5m^2+3m против ~6m^2 — выгодно при больших m. phasepoly: добавочные зёрна по бюджету ~len·n^2; берётся лучшее."""
NAME = "phasepoly_rtof4"; KIND = "lower"; HL_RAW = True
import random
from evo_lib import fold_cancel, key, expand, fold_h_choice, rccx4, inv_circ, lower_gate, HL_DIAG, combo2
SEQ = {0: [], 1: ["t"], 2: ["s"], 3: ["s", "t"], 4: ["s", "s"], 5: ["sdg", "tdg"], 6: ["sdg"], 7: ["tdg"]}
PHG = {"t": 1, "tdg": 7, "s": 2, "sdg": 6, "z": 4}; INVP = {"t": "tdg", "tdg": "t", "s": "sdg", "sdg": "s"}; CAP = 4096
def comm(g, h):
    (n1, q1), (n2, q2) = g, h
    if not set(q1) & set(q2): return True
    if n1 == "h" or n2 == "h": return False
    d1, d2 = n1 in HL_DIAG, n2 in HL_DIAG
    if (d1 and d2) or (n1 == "x" and n2 == "x"): return True
    if n1 == "cx" and n2 == "cx": return q1[0] != q2[1] and q2[0] != q1[1]
    if (d1 and n2 == "x") or (d2 and n1 == "x"): return False
    if n1 == "cx": n1, q1, n2, q2 = n2, q2, n1, q1
    if n1 == "x": return q1[0] == q2[1]
    return q2[1] not in q1  # диагональная и cx: цель cx не задета
def layers(g):
    R = list(g); out = []
    while R:
        for want_h in (False, True):
            rest, took = [], []
            for x in R:
                if (x[0] == "h") == want_h and all(comm(x, y) for y in rest): took.append(x)
                else: rest.append(x)
            R = rest
            if took: out.append(("H", [x[1][0] for x in took]) if want_h else ("B", took))
    return out
def add(P, m, v):
    P[m] = (P.get(m, 0) + v) % 8; P[m] or P.pop(m)
def block_poly(gates, W, K, P):
    def term(qs, v):
        m = 0; k = 0
        for q in qs: m ^= W[q]; k ^= K[q]
        add(P, m, -v if k else v)
    for nm, qs in gates:
        if nm == "cx": a, b = qs; W[b] ^= W[a]; K[b] ^= K[a]
        elif nm == "x": K[qs[0]] ^= 1
        elif nm in PHG: term(qs, PHG[nm])
        elif nm == "cz": a, b = qs; term((a,), 2); term((b,), 2); term((a, b), -2)
        elif nm in ("cs", "csdg"):
            w = 1 if nm == "cs" else -1; a, b = qs; term((a,), w); term((b,), w); term((a, b), -w)
        elif nm == "ccz":
            a, b, c = qs
            for s, v in (((a,), 1), ((b,), 1), ((c,), 1), ((a, b), -1), ((a, c), -1), ((b, c), -1), ((a, b, c), 1)):
                term(s, v)
        else: raise ValueError(nm)
def solve(W, n, m):
    basis = []; res = 0
    for i in range(n):
        v, c = W[i], 1 << i
        for bv, bc in basis:
            if v ^ bv < v: v ^= bv; c ^= bc
        if v: basis.append((v, c)); basis.sort(reverse=True)
    for bv, bc in basis:
        if m ^ bv < m: m ^= bv; res ^= bc
    return None if m else res
def greedy(P, W, n, rng, Lend, lam, mu=0.0, noise=0.5):
    out = []; cols = {m: solve(W, n, m) for m in P}; vals = dict(P)
    M = [solve(Lend, n, W[r]) for r in range(n)]; pc = lambda x: bin(x).count("1")
    def emit_ready():
        for m in [m for m, c in cols.items() if c & (c - 1) == 0]:
            q = cols[m].bit_length() - 1; out.extend((s, (q,)) for s in SEQ[vals[m]]); del cols[m]
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
        if lim < 0: return None
        ops = set(); lo = min(pc(y) for y in cols.values()); src_ = [x for x in cols.values() if pc(x) <= max(lo, 2)]
        for x in src_:
            bits = [k for k in range(n) if x >> k & 1]; ops |= {(a, b) for a in bits for b in bits if a != b}
        sc = [(delta(c, t) + rng.random() * noise, c, t) for c, t in ops]; _, c, t = min(sc); cnot(c, t)
    return out
def lin_fix(M, n, Q, rng):
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
        M2 = list(M); ops = (lin_greedy(M2, n, Q, rng) if k % 2 else []) + lin_fix(M2, n, Q, rng)
        if best is None or len(ops) < len(best): best = ops
    for c, t in best: M[t] ^= M[c]
    return best
def prepare(hl, n):
    g = []
    for nm, qs in hl:
        qs = tuple(qs)
        if nm == "ccx": g += [("h", (qs[2],)), ("ccz", qs), ("h", (qs[2],))]
        else: g.append((nm, qs))
    L = [1 << q for q in range(n)]; KL = [0] * n; nv = n; lay = []
    for kind, x in layers(g):
        if kind == "B": Lb = list(L); P = {}; block_poly(x, L, KL, P); lay.append(("B", P, Lb, list(L), list(KL)))
        else:
            lay.append(("H", x))
            for q in x: L[q] = 1 << nv; nv += 1
    return lay
def assign(lay, n, mode, rng):
    blocks = [k for k, x in enumerate(lay) if x[0] == "B"]
    if mode == 0 or len(blocks) < 2: return [lay[k][1] for k in blocks]
    terms = {}
    for bi, k in enumerate(blocks):
        for m, v in lay[k][1].items(): terms.setdefault(m, []).append((bi, v))
    can = {m: [bi for bi, k in enumerate(blocks) if solve(lay[k][2], n, m) is not None] for m in terms}
    newP = [dict() for _ in blocks]; left = set(terms)
    while left:
        cnt = {}
        for m in left:
            for bi in can[m]: cnt[bi] = cnt.get(bi, 0) + (1 if mode == 1 else bin(m).count("1"))
        top = max(cnt.values()); bs = [b for b in cnt if cnt[b] == top]; b = rng.choice(bs)
        for m in [m for m in left if b in can[m]]:
            for _, v in terms[m]: add(newP[b], m, v)
            left.discard(m)
    return newP
def synth(lay, Ps, n, seed, K=20):
    rng = random.Random(seed); L = [1 << q for q in range(n)]; KL = [0] * n; W = list(L); KW = [0] * n; nv = n; out = []
    def fixops(W2, Lt, Q):
        M = [solve(Lt, n, W2[r]) for r in range(n)]
        return lin_best(M, n, Q, rng, 4)
    def emit_fix(ops, Q, KLt):
        for c, t in ops: out.append(("cx", (c, t))); W[t] ^= W[c]; KW[t] ^= KW[c]
        for q in Q:
            if KW[q] != KLt[q]: out.append(("x", (q,))); KW[q] ^= 1
    bi = 0
    for k, x in enumerate(lay):
        if x[0] == "B":
            P = Ps[bi]; bi += 1; Q = set(lay[k + 1][1]) if k + 1 < len(lay) else set(range(n)); best = None
            for v in range(K if P else 1):
                W2 = list(W)
                if v < 5: net = greedy(P, W2, n, rng, x[3], (0.0, 0.3, 0.6, 1.0, 1.5)[v % 5])
                else: net = greedy(P, W2, n, rng, x[3], rng.uniform(0, 1.5), rng.uniform(0, 1.5), rng.uniform(0.1, 1.0))
                if net is None: continue
                fo = fixops(W2, x[3], Q); sc = sum(g[0] == "cx" for g in net) + len(fo)
                if best is None or sc < best[0]: best = (sc, net, fo)
            _, net, fo = best
            for nm, qs in net:
                if nm == "cx": W[qs[1]] ^= W[qs[0]]; KW[qs[1]] ^= KW[qs[0]]; out.append((nm, qs))
                else: q = qs[0]; out.append((INVP[nm], qs) if KW[q] else (nm, qs))
            L[:] = x[3]; KL[:] = x[4]; emit_fix(fo, Q, KL)
        else:
            Q = set(x[1]); emit_fix(fixops(W, L, Q), Q, KL)
            for q in x[1]: out.append(("h", (q,))); L[q] = W[q] = 1 << nv; nv += 1
    emit_fix(fixops(W, L, range(n)), range(n), KL)
    return out
def chain_choice(g, base):
    """Цепочки h(t)..h(t), где всё на t между ними — CCZ/CZ с t: всем цель t (как в combo2_permpair)."""
    choice = dict(base); last = {}
    for i, (nm, qs) in enumerate(g):
        if nm == "h":
            t = qs[0]; s = last.get(t)
            if s is not None:
                seg = [k for k in range(s + 1, i) if t in g[k][1]]
                if seg and all(g[k][0] in ("ccz", "cz") and choice.get(k) in (None, t) for k in seg):
                    for k in seg: choice[k] = t
                    last[t] = None; continue
            last[t] = i
    return choice
def analyse(g, n):
    """Классическая схема -> [(индекс ccx, A, B)], A, B — аффинные (маска, конст.) значения управляющих; иначе None."""
    poly = {}; rec = []; P = lambda q: poly.get(q) or {1 << q}
    def aff(p):
        m = k = 0
        for x in p:
            if x & (x - 1): return None
            if x: m ^= x
            else: k ^= 1
        return (m, k) if m else None
    for i, (nm, qs) in enumerate(g):
        if nm in HL_DIAG: continue
        if nm == "x": poly[qs[0]] = P(qs[0]) ^ {0}
        elif nm == "cx": poly[qs[1]] = P(qs[1]) ^ P(qs[0])
        elif nm == "ccx":
            A, B = aff(P(qs[0])), aff(P(qs[1]))
            if A is None or B is None: return None
            rec.append((i, A, B)); r = set()
            for x in P(qs[0]):
                for y in P(qs[1]): r ^= {x | y}
            poly[qs[2]] = P(qs[2]) ^ r
            if len(r) > CAP or len(poly[qs[2]]) > CAP: return None
        else: return None
    return rec
def chain_corr(P, n, rng, tries=6):
    """Чётности веса 2 — двудольный граф U-V: цепочка на U, каждая v проходит v^u1..v^u_k, откат, v ^= u_k."""
    pc = lambda x: bin(x).count("1"); adj = {}; col = {}; best = None
    if any(pc(m) > 2 for m in P): return None
    for m in P:
        if pc(m) == 2:
            a, b = [q for q in range(n) if m >> q & 1]; adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    for s0 in adj:
        if s0 in col: continue
        col[s0] = 0; st = [s0]
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in col: col[y] = col[x] ^ 1; st.append(y)
                elif col[y] == col[x]: return None
    for it in range(tries):
        U = [q for q in adj if col[q] == it % 2]; V = [q for q in adj if col[q] != it % 2]
        if it < 2: U.sort(key=lambda q: -len(adj[q]))
        else: rng.shuffle(U)
        pos = {u: i for i, u in enumerate(U)}; kv = {v: max(pos[u] for u in adj[v]) for v in V}
        out = [(g, (m.bit_length() - 1,)) for m in P if pc(m) == 1 for g in SEQ[P[m]]]
        ch = [("cx", (U[i - 1], U[i])) for i in range(len(U) - 1, 0, -1)]; out += ch
        for v in V:
            for i in range(kv[v] + 1):
                out.append(("cx", (U[i], v)))
                if U[i] in adj[v]: out += [(g, (v,)) for g in SEQ[P[(1 << U[i]) | (1 << v)]]]
        out += ch[::-1] + [("cx", (U[kv[v]], v)) for v in V]; sc = sum(x[0] == "cx" for x in out)
        if best is None or sc < best[0]: best = (sc, out)
    return best
def synth_corr(P, n, rng, K=8):
    """Фазы P на тождественном состоянии + возврат к тождеству: явная цепочка или жадная сеть, минимум CX."""
    ident = [1 << q for q in range(n)]; best = chain_corr(P, n, rng)
    for v in range(K):
        W = list(ident); net = greedy(P, W, n, rng, ident, rng.uniform(0, 1.5), rng.uniform(0, 1.5), rng.uniform(0.1, 1.0))
        if net is None: continue
        c = net + [("cx", ct) for ct in lin_best([solve(ident, n, W[r]) for r in range(n)], n, range(n), rng, 6)]
        sc = sum(x[0] == "cx" for x in c)
        if best is None or sc < best[0]: best = (sc, c)
    return None if best is None else best[1]
def rtof4(hl, n, choice, rng, tries=8):
    """Все CCX -> rccx4 (=CS^†·Tof) или обращённый (=Tof·CS); поправки CS^{±1}(A,B) — одна сеть в начале."""
    g = expand(hl, choice); rec = analyse(g, n); best = None; L = len(rec or ())
    if not rec: return None
    sidx = {r[0]: k for k, r in enumerate(rec)}
    for it in range(tries):
        sg = [1] * L if it == 0 else [-1] * L if it == 1 else [rng.choice((1, -1)) for _ in range(L)]
        sw = [0] * L if it < 2 else [rng.randint(0, 1) for _ in range(L)]; P = {}
        for (i, A, B), s in zip(rec, sg):  # s·(π/2)·A·B = s·(π/4)(A + B - A^B)
            for (m, k), v in ((A, s), (B, s), ((A[0] ^ B[0], A[1] ^ B[1]), -s)):
                if m: add(P, m, -v if k else v)
        out = synth_corr(P, n, rng) if P else []
        if out is None: continue
        for i, (nm, qs) in enumerate(g):
            if nm != "ccx": out += lower_gate(nm, qs); continue
            k = sidx[i]; a, b, t = qs
            c = rccx4(b, a, t) if sw[k] else rccx4(a, b, t); out += c if sg[k] > 0 else inv_circ(c)
        c = fold_cancel(out, n)
        if best is None or key(c) < key(best): best = c
    return best
def run(hl, n):
    hl = [(x[0], tuple(x[1])) for x in hl]; best = None; rng = random.Random(7)
    for choice in (chain_choice(hl, fold_h_choice(hl)), chain_choice(hl, {})):
        c = rtof4(hl, n, choice, rng)
        if c is not None and (best is None or key(c) < key(best)): best = c
    lay = prepare(hl, n); ns = max(6, min(40, 60000 // max(1, len(hl) * n * n)))
    for seed in range(ns):
        try: c = fold_cancel(synth(lay, assign(lay, n, (1, 2, 0)[seed % 3], random.Random(1000 + seed)), n, seed), n)
        except Exception: continue
        if best is None or key(c) < key(best): best = c
    return best if best is not None else combo2(hl, n, rounds=1, max_evals=200)
