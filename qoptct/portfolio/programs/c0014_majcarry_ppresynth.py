"""combo2_majcarry + проход ppresynth: точный пересинтез фазово-полиномиальных блоков после разложения.
Блок — CX/X/фазы на ≤K кубитах (K=4, затем 5 под бюджетом), с переносом непересекающихся/независимых команд;
он задан аффинной картой и многочленом фаз. BFS с допустимой оценкой ищет кратчайшую цепочку CX, посещающую все
чётности с ненулевой фазой и дающую ту же линейную карту (X — в конец). T не меняется; замена при CX строго меньше."""
NAME = "majcarry_ppresynth"; KIND = "lower"; HL_RAW = True
from evo_lib import (fold_cancel, key, expand, fold_h_choice, diag_on, rccx, rccx4, inv_circ, ccx_std,
                     lower_gate, HL_DIAG)

CAP = 256  # предел числа одночленов в многочлене провода (иначе проверка блока прекращается)

def pmul(P, Q):
    r = set()
    for m in P:
        for k in Q:
            r ^= {m | k}
    return r

def pair_perm(g, use4=True):
    """{индекс: (вид, сопряжённый?, кубиты)} — как evo_lib.pair_toffolis, но с точной проверкой перестановки."""
    mark = {}
    for i, (nm, qs) in enumerate(g):
        if nm != "ccx" or i in mark:
            continue
        a, b, t = qs; C, S = {a, b}, {a, b, t}
        ok3 = okc = True          # старый критерий: всё между — диагонально на S / на C
        live = True               # новый: блок моном. (классика+диагональ), многочлены ведутся
        poly = {}

        def P(q):
            return poly.get(q) or {1 << q}

        def same(q):
            p = poly.get(q)
            return p is None or p == {1 << q}

        for j in range(i + 1, len(g)):
            h = g[j]; hn, hq = h[0], tuple(h[1])
            if hn == "ccx" and set(hq[:2]) == C and j not in mark:
                n3 = live and same(a) and same(b) and same(t)
                nc = live and same(a) and same(b)
                if hq[2] == t and (ok3 or n3):
                    mark[i], mark[j] = ("R", 0, tuple(qs)), ("R", 1, tuple(qs)); break
                if use4 and (okc or nc):
                    mark[i], mark[j] = ("R4", 0, tuple(qs)), ("R4", 1, (a, b, hq[2])); break
            ok3 = ok3 and all(diag_on(h, q) for q in S)
            okc = okc and all(diag_on(h, q) for q in C)
            if live:
                if hn in HL_DIAG:
                    pass
                elif hn == "x":
                    poly[hq[0]] = P(hq[0]) ^ {0}
                elif hn == "cx":
                    poly[hq[1]] = P(hq[1]) ^ P(hq[0])
                elif hn == "ccx":
                    poly[hq[2]] = P(hq[2]) ^ pmul(P(hq[0]), P(hq[1]))
                else:
                    live = False
                if live and any(len(poly.get(q, ())) > CAP for q in hq):
                    live = False
            if not live and (not okc or (not use4 and not ok3)):
                break
    return mark

def lower_marked(g, mark, var):
    out = []
    for i, (nm, qs) in enumerate(g):
        qs = tuple(qs)
        if nm == "ccx" and i in mark:
            kind, dg, qq = mark[i]
            c = rccx(*qq) if kind == "R" else rccx4(*qq)
            out += inv_circ(c) if dg else c
        elif nm == "ccx" and i in var:
            a, b, t = qs; sw, iv = var[i]
            if sw: a, b = b, a
            d = ccx_std(a, b, t); out += inv_circ(d) if iv else d
        else:
            out += lower_gate(nm, qs)
    return out

def chain_choice(g, base):
    """Цепочки h(t)..h(t) из CCZ/CZ с t: цель t сразу всем."""
    choice = dict(base); last = {}
    for i, (nm, qs) in enumerate(g):
        if nm == "h":
            t = qs[0]; s = last.get(t)
            if s is not None:
                seg = [k for k in range(s + 1, i) if t in g[k][1]]
                if seg and all(g[k][0] in ("ccz", "cz") for k in seg) and \
                        all(choice.get(k) in (None, t) for k in seg):
                    for k in seg:
                        choice[k] = t
                    last[t] = None; continue
            last[t] = i
    return choice

def combo(hl, n, rounds=3, max_evals=None, use4=True, start=None):
    g = [(x[0], tuple(x[1])) for x in hl]; var = {}; cnt = [0]

    def full(g2):
        return fold_cancel(lower_marked(g2, pair_perm(g2, use4), var), n)

    def sc(g2):
        cnt[0] += 1; return key(full(g2))

    def spent():
        return max_evals is not None and cnt[0] >= max_evals

    c0 = fold_h_choice(g); c1 = chain_choice(g, {}); c2 = chain_choice(g, c0)
    best, choice = min(((sc(expand(g, c)), j, c) for j, c in enumerate((c0, c1, c2)) if start in (None, j)))[::2]
    choice = dict(choice)
    idx = [i for i, x in enumerate(g) if x[0] in ("ccz", "cz")]
    for _ in range(rounds):
        imp = False
        for i in idx:
            if spent(): break
            cur = choice.get(i)
            for t in [None] + list(g[i][1]):
                if t == cur: continue
                choice[i] = t; s = sc(expand(g, choice))
                if s < best: best, cur, imp = s, t, True
            choice[i] = cur
        g2 = expand(g, choice); m = pair_perm(g2, use4)
        var = {k: v for k, v in var.items() if k < len(g2) and g2[k][0] == "ccx"}
        for i, x in enumerate(g2):
            if spent(): break
            if x[0] != "ccx" or i in m: continue
            cur = var.get(i, (0, 0))
            for v in [(0, 0), (1, 0), (0, 1), (1, 1)]:
                if v == cur: continue
                var[i] = v; s = sc(g2)
                if s < best: best, cur, imp = s, v, True
            var[i] = cur
        if not imp or spent(): break
    return full(expand(g, choice))

def csim(gs, x):
    for nm, qs in gs:
        if nm == "cx":
            x ^= ((x >> qs[0]) & 1) << qs[1]
        elif nm == "ccx":
            x ^= ((x >> qs[0]) & (x >> qs[1]) & 1) << qs[2]
    return x

def same_perm(g1, g2, wires):
    w = sorted(wires); loc = {q: i for i, q in enumerate(w)}
    r1 = [(nm, tuple(loc[q] for q in qs)) for nm, qs in g1]
    r2 = [(nm, tuple(loc[q] for q in qs)) for nm, qs in g2]
    return all(csim(r1, x) == csim(r2, x) for x in range(1 << len(w)))

def classical(seq, t):
    """H(t)·[CCZ с t, CX]·H(t) -> CCX с целью t и CX (точно: H·CCZ·H = CCX)."""
    return [("ccx", tuple(x for x in qs if x != t) + (t,)) if nm == "ccz" else (nm, qs)
            for nm, qs in seq if nm != "h"]

def maj_rewrite(hl):
    g = [(x[0], tuple(x[1])) for x in hl]; out = []; i = 0; hits = 0
    while i < len(g):
        rep = None
        if g[i][0] == "h":
            t = g[i][1][0]; W = {t}; seq = []; skip = []; j = i + 1
            want = ["ccz", "cx", "ccz", "h"]
            while j < len(g) and len(seq) < 4:
                nm, qs = g[j]
                if not set(qs) & W:
                    skip.append(g[j]); j += 1; continue
                if nm != want[len(seq)]:
                    break
                seq.append(g[j]); W |= set(qs); j += 1
            if len(seq) == 4 and seq[3] == g[i] and all(not set(s[1]) & W for s in skip):
                k1, (_, (p, q)), k2, _ = seq
                X, Y = set(k1[1]), set(k2[1])
                ok = t in X and t in Y and t not in (p, q) and len(X) == len(Y) == 3
                X -= {t}; Y -= {t}; T = [("h", (t,)), None, ("h", (t,))]
                if ok and X == {p, q} and q in Y and len(Y - {p, q}) == 1:
                    a, b = p, q; c = (Y - {b}).pop(); T[1] = ("ccz", (b, c, t))
                    rep = [("cx", (a, b)), ("cx", (a, c))] + T + [("cx", (a, t)), ("cx", (a, c))]
                elif ok and Y == {p, q} and q in X and len(X - {p, q}) == 1:
                    a, b = p, q; c = (X - {b}).pop(); T[1] = ("ccz", (b, c, t))
                    rep = [("cx", (a, c)), ("cx", (a, t))] + T + [("cx", (a, c)), ("cx", (a, b))]
                if rep is not None:
                    if same_perm(classical(seq, t), classical(rep, t), W):
                        out += skip + rep; i = j; hits += 1; continue
                    rep = None
        out.append(g[i]); i += 1
    return out, hits

# ---- проход ppresynth --------------------------------------------------------------------------------------
from evo_lib import SEQ, cx as ncx_of
PHZ = {"t": 1, "s": 2, "sdg": 6, "tdg": 7, "z": 4}
BUD = [0]; MEMO = {}
def collect(c, i, K, win=400):
    S = []; B = []; bad = set(); taint = set()
    for j in range(i, min(len(c), i + win)):
        nm, qs = c[j]; Q = set(qs)
        ok = (nm in ("cx", "x") or nm in PHZ) and not (Q & bad) and len(set(S) | Q) <= K
        if ok:
            B.append(j); taint |= Q; S += [q for q in qs if q not in S]
        elif Q & taint:
            bad |= Q; taint |= Q
            if set(S) <= bad: break
    return B, S
def polyof(c, B, S):
    loc = {q: k for k, q in enumerate(S)}; k = len(S)
    rows = [1 << a for a in range(k)]; cst = [0] * k; ph = {}; m = 0
    for j in B:
        nm, qs = c[j]
        if nm == "cx":
            a, b = loc[qs[0]], loc[qs[1]]; rows[b] ^= rows[a]; cst[b] ^= cst[a]; m += 1
        elif nm == "x": cst[loc[qs[0]]] ^= 1
        else:
            w = loc[qs[0]]; ph[rows[w]] = (ph.get(rows[w], 0) + PHZ[nm] * (-1 if cst[w] else 1)) % 8
    return rows, cst, {r: v for r, v in ph.items() if v}, m
def synth(k, rows, ph, maxd, cap):
    kk = (k, tuple(rows), frozenset(ph), maxd)
    if kk not in MEMO: MEMO[kk] = synth0(k, rows, ph, maxd, cap)
    return MEMO[kk]
def synth0(k, rows, ph, maxd, cap):
    need = list(ph); nid = {m: i for i, m in enumerate(need)}; full = (1 << len(need)) - 1
    st = tuple(1 << a for a in range(k)); goal = tuple(rows)
    v0 = 0
    for r in st:
        if r in nid: v0 |= 1 << nid[r]

    L = len(need)

    def h(r, v):
        return max(L - v.bit_count(), (r[0] != goal[0]) + (r[1] != goal[1]) + sum(x != y for x, y in zip(r[2:], goal[2:])))
    if h(st, v0) > maxd: return None
    start = (st, v0); par = {start: None}; fr = [start]; d = 0
    moves = [(a, b) for a in range(k) for b in range(k) if a != b]
    while fr and d < maxd:
        d += 1; nf = []
        for s in fr:
            r, v = s
            for a, b in moves:
                nr = list(r); nr[b] ^= nr[a]; nr = tuple(nr); nv = v
                if nr[b] in nid: nv |= 1 << nid[nr[b]]
                ns = (nr, nv)
                if ns in par or d + h(nr, nv) > maxd: continue
                par[ns] = (s, (a, b)); BUD[0] += 1
                if ns == (goal, full):
                    path = []
                    while par[ns] is not None:
                        ns, mv = par[ns]; path.append(mv)
                    return path[::-1]
                nf.append(ns)
            if len(par) > cap: return None
        fr = nf
    return None
def emit(k, path, ph, cst, S):
    rows = [1 << a for a in range(k)]; done = set(); out = []

    def put(w):
        m = rows[w]
        if m in ph and m not in done:
            done.add(m); out.extend((g, (S[w],)) for g in SEQ[ph[m]])
    for w in range(k): put(w)
    for a, b in path:
        rows[b] ^= rows[a]; out.append(("cx", (S[a], S[b]))); put(b)
    return out + [("x", (S[w],)) for w in range(k) if cst[w]]
def resynth(c, K, cap, budget):
    c = list(c); i = 0
    while i < len(c):
        if c[i][0] != "cx" or BUD[0] > budget: i += 1; continue
        B, S = collect(c, i, K)
        rows, cst, ph, m = polyof(c, B, S)
        p = synth(len(S), rows, ph, m - 1, cap) if m >= 2 else None
        if p is not None and len(p) < m:
            Bs = set(B); last = B[-1]; tnt = set(); pre = []; post = []
            for j in range(i, last + 1):
                if j in Bs: tnt |= set(c[j][1]); continue
                if set(c[j][1]) & tnt: tnt |= set(c[j][1]); post.append(c[j])
                else: pre.append(c[j])
            c = c[:i] + pre + emit(len(S), p, ph, cst, S) + post + c[last + 1:]
            continue
        i += 1
    return c
def ppass(c, n, K=4, cap=40000, budget=600000):
    best = fold_cancel(c, n)
    for _ in range(6):
        c2 = fold_cancel(resynth(best, K, cap, budget), n)
        if key(c2) >= key(best): break
        best = c2
    return best
def run(hl, n):
    """Варианты родителя (старт поиска целей c0/c1/c2, без R4) дают разные схемы; проход — к каждой, лучшая."""
    BUD[0] = 0; MEMO.clear()
    h2, hits = maj_rewrite(hl); hs = [hl, h2] if hits else [hl]
    cands = [combo(h, n, rounds=3, max_evals=1500) for h in hs]
    ev = max(30, min(1500, 100000 // (len(hl) * 4 * len(hs) + 1)))
    for h in hs:
        for st, u4 in ((0, True), (1, True), (2, True), (None, False)):
            cands.append(combo(h, n, rounds=3, max_evals=ev, start=st, use4=u4))
    k = len(hs); cands = cands[:k] + sorted((fold_cancel(c, n) for c in cands[k:]), key=key)[:2]
    outs = [ppass(c, n, 4, 10000, 300000) for c in cands]
    best = min(outs, key=key)
    b5 = ppass(best, n, 5, 10000, BUD[0] + 100000)
    return min([best, b5], key=key)
