"""rtof3_vphase: Тоффоли с относительной фазой за 3 CX (Margolus), поправка которого берётся фазами по
чётностям, уже имеющимся в схеме (отложенная фаза управляющих vphase_push, обобщённая).
rccx(a,b,t) (3 CX) = CS^dag(a,b)·CZ(a,t)·Tof — проверено численно. Отсюда два разложения CCX:
  r3: Tof = rccx, затем CS(a,b)·CZ(a,t) — поправка ПОСЛЕ (на значениях после Тоффоли);
  i3: Tof = CS^dag(a,b)·CZ(a,t), затем rccx^dag — поправка ДО (на значениях до Тоффоли).
CZ(a,t) = S a · S t · S^dag(a^t), CS = T a · T b · T^dag(a^b): одиночные фазы бесплатны, а фаза на чётности
кладётся (как в phase_fold) туда, где какой-либо провод несёт эту аффинную чётность; если нет — эта часть
синтезируется явно (варианты r3z/i3z: CZ явно = 1 CX; r3c/i3c: CS явно = 2 CX), иначе RTOF4 (4 CX), иначе 7-T.
Фредкин CX(b,a)·CCX(s,a,b)·CX(b,a): в i3 поправка CZ(a,b) до Тоффоли — на проводе a в этот момент чётность
a^b, а a^b^b = a — входное значение провода a, оно в схеме есть всегда => CZ бесплатна; CS(s,a^b) явно (2 CX).
Фредкин = 1 + 3 + 2 + 1 = 7 CX (было 8). Хвостовой CX(b,a) переносится вправо через звено цепочки
CCX(a,b,t) и даёт чётность a^b для CS-поправки RTOF4 этого звена (4 CX). Ячейка csum_mux: 11 CX вместо 12.
Какие CX уровня Тоффоли переносить — координатный спуск по маске (перенос нужен звену цепочки, а иногда
выгоднее оставить CX на месте); прямая и обращённая схема; бюджет оценок обратно квадратичен по длине.
Шаблоны не завязаны на размер: работают для любого числа ячеек и для любых CCX, где чётность найдётся."""
NAME = "rtof3_vphase"; KIND = "lower"; HL_RAW = True
from collections import defaultdict
from evo_lib import lower_gate, HL_DIAG, inv_circ, fold_cancel, key, PH, SEQ, fold_h

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
    """Перенести вправо только CX уровня Тоффоли с mask[k]=1 (k — номер CX слева); с конца, индексы целы."""
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
    """Координатный спуск по маске переноса CX: перенос даёт чётность a^b для CS-поправки звена цепочки,
    а оставленный на месте хвостовой CX Фредкина даёт чётность a^b' для CZ-поправки RTOF3 (3 CX)."""
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


def run(hl, n):
    raw = [(nm, tuple(qs)) for nm, qs in hl]
    hl = [(nm, tuple(qs)) for nm, qs in fold_h(raw)]
    rev = [(INVHL.get(nm, nm), qs) for nm, qs in reversed(hl)]
    cands = [fold_cancel(lower_std(hl), n)]
    # одна оценка ~ квадратична по длине (0.12 с при 60 командах): держим поиск в нескольких секундах
    budget = max(4, min(200, int(25 * (60.0 / max(len(hl), 1)) ** 2)))
    for h0, inv in ((hl, False), (rev, True)):
        r = search(h0, n, budget)
        if r is not None: cands.append(inv_circ(r) if inv else r)
    return min(cands, key=key)


def lower_std(hl):
    out = []
    for nm, q in hl: out += lower_gate(nm, q)
    return out
