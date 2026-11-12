"""rtof3_vphase + цели CCZ по цепочкам H + ESOP-слияние команд с общей целью (уровень Тоффоли).
1) Цепочки: h(t) ... h(t), где всё между, задевающее t, — CCZ/CZ с t => у всех цель t (как в combo2_majcarry).
   Так в форме CCZ+H восстанавливаются подряд идущие Тоффоли с ОДНОЙ целью (fold_h их не видит: H сокращены).
2) ESOP по цели: команды t ^= f(управляющие) (CX, CCX, X) с общей целью коммутируют. Группа собирается вперёд
   от команды на t, пока t никто не читает (не управляющий, не H); провода между ведутся аффинно в значениях
   начала группы (CX/X — аффинно, CCX/H — провод «мёртв», член с мёртвым управляющим не берётся). f — многочлен
   Жегалкина степени ≤2; если квадратичная часть — одно произведение аффинных форм A·B, группа заменяется на
   A в провод p, B в провод q (CX; X бесплатно), CCX(p,q,t), откат, остаток — CX(v,t)/X(t). Примеры:
   ab⊕a = a·¬b (−1 CX), ab⊕ac = a(b⊕c) (2 CX вместо Тоффоли), ab⊕a⊕b = ¬(¬a¬b). Замена точна по построению
   (классические команды, та же функция). Член берётся в группу, только если f остаётся одним произведением.
3) Каждый вариант уровня Тоффоли (fold_h; цепочки+ESOP от fold_h; цепочки+ESOP с нуля) идёт в разложение
   rtof3_vphase (RTOF3/RTOF4 с поправками по чётностям, перенос CX, прямой и обращённый); бюджет оценок
   делится между вариантами (время как у rtof3_vphase); берётся лучший по (CX, T, длина)."""
NAME = "rtof3_esop"; KIND = "lower"; HL_RAW = True
from collections import defaultdict
from evo_lib import lower_gate, HL_DIAG, inv_circ, fold_cancel, key, PH, SEQ, fold_h, fold_h_choice, expand
INVHL = {"t": "tdg", "tdg": "t", "s": "sdg", "sdg": "s", "cs": "csdg", "csdg": "cs"}
def role(g, q):
    nm, qs = g
    if q not in qs: return None
    if nm in HL_DIAG: return "c"
    if nm in ("cx", "ccx", "cxf"): return "x" if q == qs[-1] else "c"
    if nm == "x": return "x"
    return "o"
def push_one(g, i):
    p, q = g[i][1]; head, rest = g[:i], g[i + 1:]; k = 0
    while k < len(rest):
        h = rest[k]
        if h[0] in ("cx", "cxf") and tuple(h[1]) == (p, q): break
        rp, rq = role(h, p), role(h, q)
        if rp in (None, "c") and rq in (None, "x"):
            k += 1; continue
        if h[0] == "ccx" and p in h[1][:2] and q in h[1][:2]:
            rest = rest[:k + 1] + [("cxf", (p, h[1][2]))] + rest[k + 1:]; k += 2; continue
        break
    return head + rest[:k] + [("cx", (p, q))] + rest[k:]
def rccx4(a, b, t):  # = CS^dag(a,b) * Tof
    return [("h", (t,)), ("t", (t,)), ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("t", (t,)),
            ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("h", (t,))]
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
def diag_on(g, q):
    nm, qs = g
    if q not in qs or nm in HL_DIAG: return True
    if nm in ("cx", "ccx", "cxf"): return q != qs[-1]
    return False
def rccx3(a, b, t):
    return rccx4(a, b, t)[:-2] + [("h", (t,))]
def pair(g):
    mark = {}
    for i, (nm, qs) in enumerate(g):
        if nm != "ccx" or i in mark: continue
        C, S = set(qs[:2]), set(qs); ok3 = True
        for j in range(i + 1, len(g)):
            h = g[j]
            if h[0] == "ccx" and set(h[1][:2]) == C and j not in mark:
                if ok3 and h[1][-1] == qs[-1]:
                    mark[i], mark[j] = ("R", 0, tuple(qs)), ("R", 1, tuple(qs)); break
                mark[i], mark[j] = ("R4", 0, tuple(qs)), ("R4", 1, tuple(qs[:2]) + (h[1][2],)); break
            ok3 = ok3 and all(diag_on(h, q) for q in S)
            if not all(diag_on(h, q) for q in C): break
    return mark
def lower(hl, choice, mark=None):
    out = []; mark = mark or {}
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
        elif md[1] == "3":  # rccx = CS^dag(a,b) CZ(a,t) Tof (поправка после) или Tof = rccx^dag ... (до)
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
            a, b = qq
            s1, v = {"cz": ("s", 6), "cs": ("t", 7), "csdg": ("tdg", 1)}[nm]
            out += [(s1, (a,)), (s1, (b,)), ("VPH", (a, b), v, i)]
    return out
def vfold(c, n):
    PHZ = dict(PH, z=4)
    wire = [(1 << q, 0) for q in range(n)]; nv = n
    first = {1 << q: (-1, q) for q in range(n)}
    acc = defaultdict(int); sk = []; vph = []
    for g in c:
        nm, qs = g[0], g[1]
        if nm in PHZ:
            m, k = wire[qs[0]]; acc[m] += PHZ[nm] * (-1 if k else 1); continue
        if nm == "VPH":
            a, b = qs; m = wire[a][0] ^ wire[b][0]; k = wire[a][1] ^ wire[b][1]
            vph.append((m, g[2] * (-1 if k else 1), g[3])); continue
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
        for q, v in ins.get(pos, ()):
            res.extend((s, (q,)) for s in SEQ[v * (-1 if wire[q][1] else 1) % 8])
    emit(-1)
    for pos, (nm, qs) in enumerate(sk):
        if nm == "h": wire[qs[0]] = (1 << nv, 0); nv += 1
        elif nm == "x": m, k = wire[qs[0]]; wire[qs[0]] = (m, k ^ 1)
        elif nm == "cx":
            a, b = qs; wire[b] = (wire[b][0] ^ wire[a][0], wire[b][1] ^ wire[a][1])
        res.append((nm, qs)); emit(pos)
    return res, set()
def solve(hl, n, use_pairs=False):
    mark = pair(hl) if use_pairs else {}
    choice = [0] * len(hl)
    for _ in range(4 * len(hl) + 4):
        r, fails = vfold(lower(hl, choice, mark), n)
        if r is not None: break
        for i in fails:
            choice[i] = min(choice[i] + 1, len(options(hl[i])) - 1)
    best = fold_cancel(r, n)
    # жадно вернуть дешёвые варианты тем, кому уже стало можно
    for i in range(len(hl)):
        if choice[i] == 0 or i in mark or hl[i][0] not in ("ccx", "ccz", "cz", "cs", "csdg"): continue
        for o in range(choice[i]):
            ch = list(choice); ch[i] = o
            r2, f2 = vfold(lower(hl, ch, mark), n)
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
        try:
            c = solve(g, n, up)
        except Exception:
            continue
        if r is None or key(c) < key(r): r = c
    return r
def search(h0, n, budget):
    m = sum(x[0] == "cx" for x in h0); cnt = 0; best = None; bm = None
    for mask in ([0] * m, [1] * m):
        r = best_solve(push_sel(h0, mask), n); cnt += 1
        if r is not None and (best is None or key(r) < key(best)): best, bm = r, mask
    if best is None: return None
    for _ in range(2):
        imp = False
        for k in range(m):
            if cnt >= budget: return best
            mk = list(bm); mk[k] ^= 1
            r = best_solve(push_sel(h0, mk), n); cnt += 1
            if r is not None and key(r) < key(best): best, bm, imp = r, mk, True
        if not imp: break
    return best
def lower_std(hl):
    out = []
    for nm, q in hl: out += lower_gate(nm, q)
    return out
def pmul(P, Q):
    r = set()
    for m in P:
        for k in Q:
            r ^= {m | k}
    return r
def chain_choice(g, base):
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
XOR = ("x", "cx", "ccx")
TC = 4          # оценка цены Тоффоли в CX при решении о замене группы
LIM = 120       # окно поиска членов группы
SCALE = 1.0     # доля бюджета rtof3_vphase на все варианты
SYN = {}
deg = lambda m: bin(m).count("1")
bits = lambda m: [v for v in range(m.bit_length()) if m >> v & 1]
lin = lambda M: {1 << v for v in bits(M)}
def synth(f):
    """f — frozenset одночленов (битмаски, 0 — константа) -> (цена, план) или None (не одно произведение)."""
    if f in SYN: return SYN[f]
    Q = {m for m in f if deg(m) == 2}; VQ = 0; best = None
    for m in Q: VQ |= m
    vs = bits(VQ); masks = [sum(1 << vs[k] for k in range(len(vs)) if s >> k & 1) for s in range(1, 1 << len(vs))]
    if any(deg(m) > 2 for m in f): pass
    elif not Q: best = (sum(1 for m in f if m), ("lin", f))
    else:
        for x, y in ((x, y) for x in range(len(masks)) for y in range(x + 1, len(masks))):
            A, B = masks[x], masks[y]
            if {m for m in pmul(lin(A), lin(B)) if deg(m) == 2} != Q: continue
            piv = next(((p, q) for p in bits(A) for q in bits(B) if p != q and (not B >> p & 1 or not A >> q & 1)), None)
            for ca, cb in ((0, 0), (0, 1), (1, 0), (1, 1)) if piv else ():
                R = set(f) ^ pmul(lin(A) ^ ({0} if ca else set()), lin(B) ^ ({0} if cb else set()))
                c = 2 * (deg(A) + deg(B) - 2) + TC + sum(1 for m in R if m)
                if all(deg(m) <= 1 for m in R) and (best is None or c < best[0]):
                    best = (c, ("prod", A, ca, B, cb, piv, frozenset(R)))
    SYN[f] = best
    return best
def emit(plan, t):
    out = []; f = plan[-1]
    if plan[0] == "prod":
        _, A, ca, B, cb, (p, q), f = plan; pre = []
        for M, c, w in (((A, ca, p), (B, cb, q)) if not B >> p & 1 else ((B, cb, q), (A, ca, p))):
            pre += [("cx", (v, w)) for v in bits(M) if v != w] + ([("x", (w,))] if c else [])
        out = pre + [("ccx", (p, q, t))] + pre[::-1]
    return out + [("cx", (v, t)) for m in sorted(f) if m for v in bits(m)] + ([("x", (t,))] if 0 in f else [])
def esop(g):
    """g — уровень Тоффоли с целями -> (новый g, число замен)."""
    g = [(x[0], tuple(x[1])) for x in g]; used = set(); out = []; hits = 0
    for i, G in enumerate(g):
        if i in used: continue
        if G[0] not in XOR: out.append(G); continue
        t = G[1][-1]; aff = {}; dead = set(); mem = [i]
        A = lambda q: aff.get(q, {1 << q})
        term = lambda K: pmul(A(K[1][0]), A(K[1][1])) if K[0] == "ccx" else set(A(K[1][0])) if K[0] == "cx" else {0}
        f = set(term(G))
        for j in range(i + 1, min(len(g), i + LIM)):
            if j in used: continue
            nm, qs = g[j]
            if nm in XOR and qs[-1] == t:
                if not set(qs[:-1]) & dead and synth(frozenset(f ^ term(g[j]))) is not None:
                    f ^= term(g[j]); mem.append(j)
                continue
            if t in qs: break
            if nm == "cx":
                if qs[0] in dead: dead.add(qs[1])
                elif qs[1] not in dead: aff[qs[1]] = A(qs[1]) ^ A(qs[0])
            elif nm == "x": aff[qs[0]] = A(qs[0]) ^ {0}
            elif nm == "ccx": dead.add(qs[2])
            elif nm not in HL_DIAG: dead |= set(qs)
        s = synth(frozenset(f))
        old = sum(1 if g[k][0] == "cx" else TC if g[k][0] == "ccx" else 0 for k in mem)
        if len(mem) > 1 and s is not None and s[0] < old:
            used |= set(mem); out += emit(s[1], t); hits += 1
        else:
            out.append(G)
    return out, hits
def inv_hl(g):
    return [(INVHL.get(nm, nm), qs) for nm, qs in reversed(g)]
def core(hl, n, scale):
    rev = inv_hl(hl); cands = [fold_cancel(lower_std(hl), n)]
    budget = max(3, min(200, int(scale * 25 * (60.0 / max(len(hl), 1)) ** 2)))
    for h0, inv in ((hl, False), (rev, True)):
        r = search(h0, n, budget)
        if r is not None: cands.append(inv_circ(r) if inv else r)
    return min(cands, key=key)
def run(hl, n):
    raw = [(nm, tuple(qs)) for nm, qs in hl]; c0 = fold_h_choice(raw)
    vs = [expand(raw, c0)]
    for ch in (chain_choice(raw, {}), chain_choice(raw, c0)):
        e = esop(expand(raw, ch))[0]
        if e not in vs: vs.append(e)
    sc = SCALE / len(vs)
    return min((core(v, n, sc) for v in vs), key=key)
