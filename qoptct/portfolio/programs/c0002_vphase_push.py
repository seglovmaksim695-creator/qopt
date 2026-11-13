"""Отложенная фаза управляющих (virtual control phase) + перенос CX вправо через Тоффоли.
CCX(a,b,t) = CS(a,b) * RTOF4(a,b,t): RTOF4 стоит 4 CX, а CS(a,b) — это фазы T на a, T на b и T^-1 на
чётности a^b. Фазу на a^b не синтезируем у гейта, а кладём (как phase_fold) в ПЕРВОЕ место схемы, где
какой-либо провод несёт ту же аффинную чётность (по сумме путей). Если такой чётности в схеме нет —
откат к стандартному разложению этого гейта (жадно, до неподвижной точки).
Чтобы чётности появлялись, CX(p,q) переносятся вправо: через команды, где p — управляющий/диагональ, а q —
X-цель, и через CCX(q,p,t) с превращением в CCX(q,p,t)·CX(p,t) (CX(p,t) рядом с целью t раскладывается как
H·CZ·H фазами — бесплатно). В csum_mux хвостовой CX(b,a) Фредкина уходит в конец и даёт чётность a^b для
CS-поправки звена цепочки: −2 CX на каждое звено. То же для CZ/CS (CS = 2 CX -> 0 при найденной чётности)."""
NAME = "vphase_push"; KIND = "lower"
from collections import defaultdict
from evo_lib import lower_gate, HL_DIAG, inv_circ, fold_cancel, key, PH, SEQ

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


def push_all(hl):
    g = list(hl)
    for i in range(len(g) - 1, -1, -1):
        if g[i][0] == "cx": g = push_one(g, i)
    return g


def rccx4(a, b, t):  # = CS^dag(a,b) * Tof
    return [("h", (t,)), ("t", (t,)), ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("t", (t,)),
            ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("h", (t,))]


def options(g):
    nm, qs = g
    if nm == "ccx": return [("r4", qs), ("std", qs)]
    if nm == "ccz":
        a, b, c = qs
        return [("r4", (a, b, c)), ("r4", (a, c, b)), ("r4", (b, c, a)), ("std", qs)]
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
    """rtof-пары (как у родителя): CCX с теми же управляющими, между ними опора диагональна."""
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
        elif md == "r4":
            a, b, t = qq; body = rccx4(a, b, t)
            if nm == "ccz": body = [("h", (t,))] + body + [("h", (t,))]
            out += body + [("t", (a,)), ("t", (b,)), ("VPH", (a, b), 7, i)]
        else:  # ph: cz = S a S b Sdg(a^b); cs = T a T b Tdg(a^b); csdg — сопряжённое
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


def run(hl, n):
    hl = [(nm, tuple(qs)) for nm, qs in hl]
    rev = [(INVHL.get(nm, nm), qs) for nm, qs in reversed(hl)]
    cands = []
    for h0, inv in ((hl, False), (rev, True)):
        for g in (h0, push_all(h0)):
            for up in (False, True):
                try:
                    r = solve(g, n, up)
                except Exception:
                    continue
                cands.append(inv_circ(r) if inv else r)
    cands.append(fold_cancel(lower_std(hl), n))
    return min(cands, key=key)


def lower_std(hl):
    out = []
    for nm, q in hl: out += lower_gate(nm, q)
    return out
