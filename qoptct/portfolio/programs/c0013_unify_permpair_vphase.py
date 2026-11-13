"""unify_permpair_vphase: один конвейер вместо c0010 majcarry + c0011 rtof3_vphase + c0008 cla_collapse, и
скрещивание: точная (АНФ) проверка пар CCX из combo2_permpair внутри разложения с отложенными фазами RTOF3/RTOF4.
Старты уровня Тоффоли: fold_h; multi_choice + свёртка пар CCX (cla); MAJ-перенос (vbe). Каждый — спуск по маске
переноса CX (c0011); общий бюджет оценок делится между стартами. Лучшее по (CX, T, длина)."""
NAME = "unify_permpair_vphase"; KIND = "lower"; HL_RAW = True
from collections import defaultdict
from evo_lib import (lower_gate, HL_DIAG, inv_circ, fold_cancel, key, PH, SEQ, expand, fold_h_choice, diag_on)
INVHL = {"t": "tdg", "tdg": "t", "s": "sdg", "sdg": "s", "cs": "csdg", "csdg": "cs"}
CAP = 256
# ---------- точные пары CCX (c0006/c0010) ----------
def pmul(P, Q):
    r = set()
    for m in P:
        for k in Q: r ^= {m | k}
    return r
def pair_perm(g):
    mark = {}
    for i, (nm, qs) in enumerate(g):
        if nm != "ccx" or i in mark: continue
        a, b, t = qs; C, S = {a, b}, {a, b, t}; ok3 = okc = live = True; poly = {}
        P = lambda q: poly.get(q) or {1 << q}
        same = lambda q: poly.get(q) is None or poly[q] == {1 << q}
        for j in range(i + 1, len(g)):
            hn, hq = g[j][0], tuple(g[j][1])
            if hn == "ccx" and set(hq[:2]) == C and j not in mark:
                nc = live and same(a) and same(b)
                if hq[2] == t and (ok3 or (nc and same(t))):
                    mark[i], mark[j] = ("R", 0, tuple(qs)), ("R", 1, tuple(qs)); break
                if okc or nc:
                    mark[i], mark[j] = ("R4", 0, tuple(qs)), ("R4", 1, (a, b, hq[2])); break
            ok3 = ok3 and all(diag_on(g[j], q) for q in S)
            okc = okc and all(diag_on(g[j], q) for q in C)
            if live:
                if hn in HL_DIAG: pass
                elif hn == "x": poly[hq[0]] = P(hq[0]) ^ {0}
                elif hn == "cx": poly[hq[1]] = P(hq[1]) ^ P(hq[0])
                elif hn == "ccx": poly[hq[2]] = P(hq[2]) ^ pmul(P(hq[0]), P(hq[1]))
                else: live = False
                if live and any(len(poly.get(q, ())) > CAP for q in hq): live = False
            if not live and not okc: break
    return mark
# ---------- MAJ-перенос (c0010) ----------
def csim(gs, x):
    for nm, qs in gs:
        if nm == "cx": x ^= ((x >> qs[0]) & 1) << qs[1]
        elif nm == "ccx": x ^= ((x >> qs[0]) & (x >> qs[1]) & 1) << qs[2]
    return x
def same_perm(g1, g2, wires):
    loc = {q: i for i, q in enumerate(sorted(wires))}
    r1 = [(nm, tuple(loc[q] for q in qs)) for nm, qs in g1]
    r2 = [(nm, tuple(loc[q] for q in qs)) for nm, qs in g2]
    return all(csim(r1, x) == csim(r2, x) for x in range(1 << len(loc)))
def classical(seq, t): return [("ccx", tuple(x for x in qs if x != t) + (t,)) if nm == "ccz" else (nm, qs)
            for nm, qs in seq if nm != "h"]
def maj_rewrite(g):
    out = []; i = 0; hits = 0
    while i < len(g):
        rep = None
        if g[i][0] == "h":
            t = g[i][1][0]; W = {t}; seq = []; skip = []; j = i + 1; want = ["ccz", "cx", "ccz", "h"]
            while j < len(g) and len(seq) < 4:
                nm, qs = g[j]
                if not set(qs) & W:
                    skip.append(g[j]); j += 1; continue
                if nm != want[len(seq)]: break
                seq.append(g[j]); W |= set(qs); j += 1
            if len(seq) == 4 and seq[3] == g[i] and all(not set(s[1]) & W for s in skip):
                k1, (_, (p, q)), k2, _ = seq; X, Y = set(k1[1]), set(k2[1])
                ok = t in X and t in Y and t not in (p, q) and len(X) == len(Y) == 3
                X -= {t}; Y -= {t}; T = [("h", (t,)), None, ("h", (t,))]
                if ok and X == {p, q} and q in Y and len(Y - {p, q}) == 1:
                    a, b = p, q; c = (Y - {b}).pop(); T[1] = ("ccz", (b, c, t))
                    rep = [("cx", (a, b)), ("cx", (a, c))] + T + [("cx", (a, t)), ("cx", (a, c))]
                elif ok and Y == {p, q} and q in X and len(X - {p, q}) == 1:
                    a, b = p, q; c = (X - {b}).pop(); T[1] = ("ccz", (b, c, t))
                    rep = [("cx", (a, c)), ("cx", (a, t))] + T + [("cx", (a, c)), ("cx", (a, b))]
                if rep is not None and same_perm(classical(seq, t), classical(rep, t), W):
                    out += skip + rep; i = j; hits += 1; continue
        out.append(g[i]); i += 1
    return out, hits
# ---------- свёртка пар CCX и цели-цепочки (c0008) ----------
def anf_linear(terms, nv):
    vs = sorted({b for p in terms for m, _ in p for b in range(nv) if m >> b & 1})
    if len(vs) > 12: return None
    def f(x):
        r = 0
        for p in terms:
            v = 1
            for m, c in p: v &= (bin(m & x).count("1") + c) & 1
            r ^= v
        return r
    c0 = f(0); mask = 0
    for b in vs:
        if f(1 << b) ^ c0: mask |= 1 << b
    for k in range(1 << len(vs)):
        x = sum(1 << b for i, b in enumerate(vs) if k >> i & 1)
        if f(x) != (bin(mask & x).count("1") + c0) & 1: return None
    return mask, c0
def collapse_once(g, n):
    for i, (nm, qs) in enumerate(g):
        if nm != "ccx": continue
        a, b, t = qs; form = {w: (1 << w, 0) for w in range(n)}; bad = set()
        for j in range(i + 1, len(g)):
            h, hq = g[j]
            if h == "ccx" and hq[2] == t:
                c, d = hq[0], hq[1]
                if c in bad or d in bad: break
                lin = anf_linear([[(1 << a, 0), (1 << b, 0)], [form[c], form[d]]], n)
                if lin is None: break
                mask, c0 = lin
                ins = [("cx", (w, t)) for w in range(n) if mask >> w & 1] + ([("x", (t,))] if c0 else [])
                return g[:i] + ins + g[i + 1:j] + g[j + 1:]
            if t in hq:
                if (h in ("cx", "ccx") and hq[-1] == t) or h == "x": continue
                break
            if h == "cx":
                p, q = hq
                if p in bad: bad.add(q)
                else: form[q] = (form[q][0] ^ form[p][0], form[q][1] ^ form[p][1])
            elif h == "x": form[hq[0]] = (form[hq[0]][0], form[hq[0]][1] ^ 1)
            elif h == "ccx": bad.add(hq[2])
            elif h not in ("t", "tdg", "s", "sdg", "z", "cz", "ccz", "cs", "csdg"): bad |= set(hq)
    return None
def collapse(g, n):
    while True:
        r = collapse_once(g, n)
        if r is None: return g
        g = r
def multi_choice(g):
    ch = fold_h_choice(g)
    for i, (nm, qs) in enumerate(g):
        if nm != "h": continue
        t = qs[0]; idx = []; j = i + 1
        while j < len(g):
            if t in g[j][1]:
                if g[j] == ("h", (t,)): break
                if g[j][0] in ("ccz", "cz"): idx.append(j)
                else: idx = None; break
            j += 1
        if j >= len(g) or not idx or any(k in ch and ch[k] != t for k in idx): continue
        for k in idx: ch[k] = t
    return ch
def rev_hl(g): return [(INVHL.get(nm, nm), qs) for nm, qs in reversed(g)]
# ---------- разложение с отложенными фазами, RTOF3/RTOF4 (c0011) ----------
def role(g, q):
    nm, qs = g
    if q not in qs: return None
    if nm in HL_DIAG: return "c"
    if nm in ("cx", "ccx", "cxf"): return "x" if q == qs[-1] else "c"
    return "x" if nm == "x" else "o"
def push_one(g, i):
    p, q = g[i][1]; head, rest = g[:i], g[i + 1:]; k = 0
    while k < len(rest):
        h = rest[k]
        if h[0] in ("cx", "cxf") and tuple(h[1]) == (p, q): break
        if role(h, p) in (None, "c") and role(h, q) in (None, "x"):
            k += 1; continue
        if h[0] == "ccx" and p in h[1][:2] and q in h[1][:2]:
            rest = rest[:k + 1] + [("cxf", (p, h[1][2]))] + rest[k + 1:]; k += 2; continue
        break
    return head + rest[:k] + [("cx", (p, q))] + rest[k:]
def rccx4(a, b, t):  # = CS^dag(a,b) * Tof
    return [("h", (t,)), ("t", (t,)), ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("t", (t,)),
            ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("h", (t,))]
def rccx3(a, b, t): return rccx4(a, b, t)[:-2] + [("h", (t,))]
def options(g):
    nm, qs = g
    if nm == "ccx":
        a, b, t = qs; o = ((a, b, t), (b, a, t))
        return [(m, x) for m in ("r3", "i3") for x in o] + [("r4", qs)] + \
            [(m, x) for m in ("r3z", "i3z", "r3c", "i3c") for x in o] + [("std", qs)]
    if nm == "ccz":
        a, b, c = qs; tr = ((a, b, c), (a, c, b), (b, c, a))
        return [("r3", x) for x in tr] + [("r3", (x[1], x[0], x[2])) for x in tr] + [("r4", x) for x in tr] + [("std", qs)]
    if nm in ("cz", "cs", "csdg"): return [("ph", qs), ("std", qs)]
    return [("std", qs)]
def lower(hl, choice, mark):
    out = []
    for i, (nm, qs) in enumerate(hl):
        if i in mark:
            kind, dg, qq = mark[i]; c = rccx3(*qq) if kind == "R" else rccx4(*qq)
            out += inv_circ(c) if dg else c; continue
        md, qq = options((nm, qs))[choice[i]] if nm in ("ccx", "ccz", "cz", "cs", "csdg") else ("std", qs)
        if nm == "cxf":
            p, t = qs
            out += [("h", (t,)), ("s", (p,)), ("s", (t,)), ("cx", (p, t)), ("sdg", (t,)), ("cx", (p, t)), ("h", (t,))]
        elif md == "std":
            out += lower_gate(nm, qs)
        elif md[1] == "3":
            a, b, t = qq; hh = [("h", (t,))] if nm == "ccz" else []; pre = md[0] == "i"
            cz = lower_gate("cz", (a, t)) if md[2:] == "z" else [("s", (a,)), ("s", (t,)), ("VPH", (a, t), 6, i)]
            if md[2:] == "c": cs = lower_gate("csdg" if pre else "cs", (a, b))
            elif pre: cs = [("tdg", (a,)), ("tdg", (b,)), ("VPH", (a, b), 1, i)]
            else: cs = [("t", (a,)), ("t", (b,)), ("VPH", (a, b), 7, i)]
            out += hh + (cz + cs + inv_circ(rccx3(a, b, t)) if pre else rccx3(a, b, t) + cz + cs) + hh
        elif md == "r4":
            a, b, t = qq; body = rccx4(a, b, t)
            if nm == "ccz": body = [("h", (t,))] + body + [("h", (t,))]
            out += body + [("t", (a,)), ("t", (b,)), ("VPH", (a, b), 7, i)]
        else:
            a, b = qq; s1, v = {"cz": ("s", 6), "cs": ("t", 7), "csdg": ("tdg", 1)}[nm]
            out += [(s1, (a,)), (s1, (b,)), ("VPH", (a, b), v, i)]
    return out
def vfold(c, n):
    PHZ = dict(PH, z=4); wire = [(1 << q, 0) for q in range(n)]; nv = n
    first = {1 << q: (-1, q) for q in range(n)}; acc = defaultdict(int); sk = []; vph = []
    for g in c:
        nm, qs = g[0], g[1]
        if nm in PHZ:
            m, k = wire[qs[0]]; acc[m] += PHZ[nm] * (-1 if k else 1); continue
        if nm == "VPH":
            a, b = qs; vph.append((wire[a][0] ^ wire[b][0], g[2] * (-1 if wire[a][1] ^ wire[b][1] else 1), g[3])); continue
        if nm == "h":
            wire[qs[0]] = (1 << nv, 0); nv += 1; first[wire[qs[0]][0]] = (len(sk), qs[0])
        elif nm == "x":
            m, k = wire[qs[0]]; wire[qs[0]] = (m, k ^ 1)
        elif nm == "cx":
            a, b = qs; wire[b] = (wire[b][0] ^ wire[a][0], wire[b][1] ^ wire[a][1])
            if wire[b][0] not in first: first[wire[b][0]] = (len(sk), b)
        sk.append((nm, tuple(qs)))
    fails = set()
    for m, v, tag in vph:
        if m == 0: continue
        if m in first: acc[m] += v
        else: fails.add(tag)
    if fails: return None, fails
    ins = defaultdict(list)
    for m, v in acc.items():
        if v % 8 and m: ins[first[m][0]].append((first[m][1], v % 8))
    wire = [(1 << q, 0) for q in range(n)]; nv = n; res = []
    def emit(pos):
        for q, v in ins.get(pos, ()): res.extend((s, (q,)) for s in SEQ[v * (-1 if wire[q][1] else 1) % 8])
    emit(-1)
    for pos, (nm, qs) in enumerate(sk):
        if nm == "h": wire[qs[0]] = (1 << nv, 0); nv += 1
        elif nm == "x": m, k = wire[qs[0]]; wire[qs[0]] = (m, k ^ 1)
        elif nm == "cx":
            a, b = qs; wire[b] = (wire[b][0] ^ wire[a][0], wire[b][1] ^ wire[a][1])
        res.append((nm, qs)); emit(pos)
    return res, set()
def solve(hl, n, use_pairs):
    mark = pair_perm(hl) if use_pairs else {}; choice = [0] * len(hl)
    for _ in range(4 * len(hl) + 4):
        r, fails = vfold(lower(hl, choice, mark), n)
        if r is not None: break
        for i in fails: choice[i] = min(choice[i] + 1, len(options(hl[i])) - 1)
    best = fold_cancel(r, n)
    for i in range(len(hl)):
        if choice[i] == 0 or i in mark or hl[i][0] not in ("ccx", "ccz", "cz", "cs", "csdg"): continue
        for o in range(choice[i]):
            ch = list(choice); ch[i] = o; r2, _ = vfold(lower(hl, ch, mark), n)
            if r2 is not None:
                r2 = fold_cancel(r2, n)
                if key(r2) < key(best): best, choice = r2, ch
                break
    return best
def push_sel(g, mask):
    pos = [i for i, x in enumerate(g) if x[0] == "cx"]
    for k in range(len(pos) - 1, -1, -1):
        if mask[k]: g = push_one(g, pos[k])
    return g
def best_solve(g, n):
    r = None
    for up in (False, True):
        try: c = solve(g, n, up)
        except Exception: continue
        if r is None or key(c) < key(r): r = c
    return r
def search(h0, n, budget):
    m = sum(x[0] == "cx" for x in h0); cnt = 0; best = bm = None
    for mask in ([0] * m, [1] * m):
        r = best_solve(push_sel(h0, mask), n); cnt += 1
        if r is not None and (best is None or key(r) < key(best)): best, bm = r, mask
    if best is None: return None
    for _ in range(2):
        imp = False
        for k in range(m):
            if cnt >= budget: return best
            mk = list(bm); mk[k] ^= 1; r = best_solve(push_sel(h0, mk), n); cnt += 1
            if r is not None and key(r) < key(best): best, bm, imp = r, mk, True
        if not imp: break
    return best
def run(hl, n):
    g = [(nm, tuple(qs)) for nm, qs in hl]
    h0 = expand(g, multi_choice(g)); starts = [expand(g, fold_h_choice(g)), collapse(h0, n),
                                              rev_hl(collapse(rev_hl(h0), n))]
    g2, hits = maj_rewrite(g)
    if hits: starts.append(expand(g2, fold_h_choice(g2)))
    uniq = []
    for s in starts:
        s = [(nm, tuple(qs)) for nm, qs in s]
        if s not in uniq: uniq.append(s)
    cands = [fold_cancel([x for nm, q in uniq[0] for x in lower_gate(nm, q)], n)]
    share = 1.5 / len(uniq)  # общий бюджет ~ c0011 x 1.5 (0.12 с на оценку при 60 командах)
    for s in uniq:
        budget = max(4, int(share * min(200, 25 * (60.0 / max(len(s), 1)) ** 2)))
        for h1, inv in ((s, False), (rev_hl(s), True)):
            r = search(h1, n, budget)
            if r is not None: cands.append(inv_circ(r) if inv else r)
    return min(cands, key=key)
