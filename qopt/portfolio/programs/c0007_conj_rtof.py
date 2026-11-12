"""conj_rtof: сопряжение S·B·S^{-1} -> S'·B·S'^{-1}, где S' — S с Тоффоли с относительной фазой.
Если S состоит только из мономиальных команд (перестановка × диагональ: ccx, cx, x, ccz, cz, z, cs, ...),
то любое разложение S' = E·S, E диагональна на supp(S). Если B диагонально на supp(S) (кубиты supp(S) у команд B —
только управляющие/фазовые), то E коммутирует с B и S'^{-1}·B·S' = S^{-1}·B·S точно. Значит, ВСЕ CCX внутри S
можно брать как rccx (3 CX), а S^{-1} — как точное обращение той же схемы. Для mct_dirty (Барэнко: T·L·T·L,
L палиндром на грязных) это даёт 2 точных Тоффоли + (4k-10) rccx = 12k-18 CX (рекорд barenco_tof_10: 102).
Остальное — пары CCX (pair_toffolis) и локальный поиск целей CCZ/вариантов, как в combo2, с этим же разложением."""
NAME = "conj_rtof"; KIND = "lower"; HL_RAW = True
from evo_lib import (fold_cancel, key, inv_circ, lower_gate, ccx_std, rccx, fold_h_choice,
                     expand, diag_on, HL_DIAG)

MONO = {"ccx", "cx", "x", "ccz", "cz", "z", "cs", "csdg", "t", "tdg", "s", "sdg"}
HINV = {"t": "tdg", "tdg": "t", "s": "sdg", "sdg": "s", "cs": "csdg", "csdg": "cs"}
SYM = {"ccz", "cz", "cs", "csdg"}


def inv_eq(a, b):
    """Команда b — обратная к a (уровень Тоффоли)."""
    if HINV.get(a[0], a[0]) != b[0]:
        return False
    if a[0] in SYM:
        return set(a[1]) == set(b[1])
    if a[0] == "ccx":
        return set(a[1][:2]) == set(b[1][:2]) and a[1][2] == b[1][2]
    return a[1] == b[1]


def find_conj(g, lo, hi):
    """Жадно: непересекающиеся тройки (i, L, e): g[i:i+L]=S, g[e-L+1:e+1]=S^{-1}, между — B, диагональное на supp(S).
    Выбор по числу CCX в S (больше — лучше)."""
    cands = []
    for i in range(lo, hi):
        if g[i][0] not in MONO:
            continue
        for e in range(hi - 1, i, -1):
            if not inv_eq(g[i], g[e]):
                continue
            best = None; sup = set(); L = 0
            while i + L < e - L and g[i + L][0] in MONO and inv_eq(g[i + L], g[e - L]):
                sup |= set(g[i + L][1]); L += 1
                ok = all(all(diag_on(h, q) for q in sup) for h in g[i + L:e - L + 1])
                if ok:
                    nc = sum(x[0] == "ccx" for x in g[i:i + L])
                    if nc and (best is None or nc > best[0]):
                        best = (nc, L)
            if best:
                cands.append((best[0], -(e - i), i, best[1], e))
    cands.sort(reverse=True)
    used = [False] * len(g); res = []
    for nc, _, i, L, e in cands:
        span = list(range(i, i + L)) + list(range(e - L + 1, e + 1))
        if any(used[k] for k in span):
            continue
        # пересечение с уже выбранными: только вложение целиком в B или непересечение
        bad = False
        for i2, L2, e2 in res:
            if e2 < i or e < i2: continue
            if i + L <= i2 and e2 <= e - L: continue
            if i2 + L2 <= i and e <= e2 - L2: continue
            bad = True; break
        if bad:
            continue
        for k in span:
            used[k] = True
        res.append((i, L, e))
    return res


def rel_var(qs, v):
    a, b, t = qs
    if v & 1: a, b = b, a
    c = rccx(a, b, t)
    return inv_circ(c) if v & 2 else c


def std_var(qs, v):
    a, b, t = qs
    if v & 1: a, b = b, a
    c = ccx_std(a, b, t)
    return inv_circ(c) if v & 2 else c


def pairs(g, conj):
    """Как evo_lib.pair_toffolis, но пара не должна пересекать границу B какого-либо сопряжения."""
    mark = {}

    def same_side(p, q):
        return all((i + L <= p <= e - L) == (i + L <= q <= e - L) for i, L, e in conj)

    for i, (nm, qs) in enumerate(g):
        if nm != "ccx" or i in mark:
            continue
        C, S = set(qs[:2]), set(qs); ok3 = True
        for j in range(i + 1, len(g)):
            h = g[j]
            if h[0] == "ccx" and set(h[1][:2]) == C and j not in mark and same_side(i, j):
                if ok3 and h[1][-1] == qs[-1]:
                    mark[i], mark[j] = ("R", 0, tuple(qs)), ("R", 1, tuple(qs)); break
                mark[i], mark[j] = ("R4", 0, tuple(qs)), ("R4", 1, tuple(qs[:2]) + (h[1][2],)); break
            ok3 = ok3 and all(diag_on(h, q) for q in S)
            if not all(diag_on(h, q) for q in C):
                break
    return mark


def lower_all(g, conj, var):
    """conj — список (i, L, e); var — {индекс: вариант 0..3}."""
    role = {}
    for i, L, e in conj:
        for r in range(L):
            role[i + r] = ("S", i + r)
            role[e - r] = ("Si", i + r)
    # пары с относительной фазой — только среди команд вне сопряжений
    g_p = [(("ccx_c", x[1]) if k in role else x) for k, x in enumerate(g)]
    mark = pairs(g_p, conj)
    pieces = {}
    for k in sorted(role):
        tag, src = role[k]
        if tag == "S":
            nm, qs = g[k]
            pieces[k] = rel_var(qs, var.get(k, 0)) if nm == "ccx" else lower_gate(nm, qs)
    out = []
    for k, (nm, qs) in enumerate(g):
        qs = tuple(qs)
        if k in role:
            tag, src = role[k]
            out += pieces[k] if tag == "S" else inv_circ(pieces[src])
        elif nm == "ccx" and k in mark:
            kind, dg, qq = mark[k]
            c = rccx(*qq) if kind == "R" else rccx(*qq)[:-1] + [("cx", (qq[0], qq[2])), ("h", (qq[2],))]
            out += inv_circ(c) if dg else c
        elif nm == "ccx":
            out += std_var(qs, var.get(k, 0))
        else:
            out += lower_gate(nm, qs)
    return out


def run(hl, n):
    g0 = [(x[0], tuple(x[1])) for x in hl]
    return _search(g0, n, max_evals=600 if len(g0) < 150 else 300)


def frame_choice(g):
    """Цель CCZ — единственный кубит, находящийся «в базисе H» (нечётное число H до команды на этом проводе)."""
    par = {}; ch = {}
    for i, (nm, qs) in enumerate(g):
        if nm == "h":
            par[qs[0]] = par.get(qs[0], 0) ^ 1
        elif nm in ("ccz", "cz"):
            odd = [q for q in qs if par.get(q, 0)]
            if len(odd) == 1:
                ch[i] = odd[0]
    return ch


def _search(g0, n, rounds=2, max_evals=600):
    cnt = [0]

    def build(choice, var):
        g = expand(g0, choice)
        conj = find_conj(g, 0, len(g))
        return g, conj

    def score(choice, var):
        cnt[0] += 1
        g, conj = build(choice, var)
        return key(fold_cancel(lower_all(g, conj, var), n))

    var = {}; fh = fold_h_choice(g0); fr = dict(fh); fr.update(frame_choice(g0))
    best, choice = None, None
    for ch in (fr, fh):
        s = score(ch, var)
        if best is None or s < best:
            best, choice = s, dict(ch)
    idx = [i for i, x in enumerate(g0) if x[0] in ("ccz", "cz")]
    for _ in range(rounds):
        imp = False
        for i in idx:
            if cnt[0] >= max_evals: break
            cur = choice.get(i)
            for t in [None] + list(g0[i][1]):
                if t == cur: continue
                choice[i] = t; s = score(choice, var)
                if s < best: best, cur, imp = s, t, True
            choice[i] = cur
        g, conj = build(choice, var)
        var = {k: v for k, v in var.items() if k < len(g) and g[k][0] == "ccx"}
        inS = set()
        for i, L, e in conj:
            inS |= set(range(i, i + L))
        for k, x in enumerate(g):
            if cnt[0] >= max_evals: break
            if x[0] != "ccx": continue
            cur = var.get(k, 0)
            for v in range(4):
                if v == cur: continue
                var[k] = v; s = score(choice, var)
                if s < best: best, cur, imp = s, v, True
            var[k] = cur
        if not imp or cnt[0] >= max_evals: break
    g, conj = build(choice, var)
    return fold_cancel(lower_all(g, conj, var), n)
