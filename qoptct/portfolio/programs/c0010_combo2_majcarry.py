"""combo2_permpair + перенос вычисления переноса (VBE carry) в форму MAJ с одним Тоффоли.
Блок carry уровня Тоффоли H(t)·CCZ(a,b,t)·CX(a,b)·CCZ(c,b,t)·H(t) = CCX(a,b,t)·CX(a,b)·CCX(c,b,t):
t ^= ab ⊕ c(a⊕b) = MAJ(a,b,c), b -> a⊕b. Тождество MAJ = a ⊕ (a⊕b)(a⊕c) даёт тот же блок с ОДНИМ Тоффоли:
CX(a,b) CX(a,c) CCX(b',c',t) CX(a,t) CX(a,c); откат (uncarry) — зеркально. CX(a,t) ставится СНАРУЖИ интервала
пары (до Тоффоли в uncarry, после — в carry), иначе t не возвращается между парой и permpair теряет вид R (3 CX).
Каждая замена проверяется точно (таблица истинности: блоки классические, H·CCZ·H = CCX без фаз). Блоки ищутся по
проводам (команды на непересекающихся проводах пропускаются). Итог — лучшее из combo(исходник) и combo(замена)."""
NAME = "combo2_majcarry"; KIND = "lower"; HL_RAW = True
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
    """Цепочки: h(t) ... h(t), где каждая команда между ними, задевающая t, — CCZ/CZ с t. Тогда все они
    получают цель t сразу (поодиночке локальный поиск этого не находит: H между ними не сокращаются)."""
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


def combo(hl, n, rounds=3, max_evals=None, use4=True):
    g = [(x[0], tuple(x[1])) for x in hl]; var = {}; cnt = [0]

    def full(g2):
        return fold_cancel(lower_marked(g2, pair_perm(g2, use4), var), n)

    def sc(g2):
        cnt[0] += 1; return key(full(g2))

    def spent():
        return max_evals is not None and cnt[0] >= max_evals

    c0 = fold_h_choice(g); c1 = chain_choice(g, {}); c2 = chain_choice(g, c0)
    best, choice = min(((sc(expand(g, c)), j, c) for j, c in enumerate((c0, c1, c2))))[::2]
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


def run(hl, n):
    best = combo(hl, n, rounds=3, max_evals=1500)
    h2, hits = maj_rewrite(hl)
    if hits:
        c2 = combo(h2, n, rounds=3, max_evals=1500)
        if key(fold_cancel(c2, n)) < key(fold_cancel(best, n)):
            best = c2
    return best
