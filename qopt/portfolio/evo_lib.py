"""evo_lib — библиотека, доступная программам-кандидатам (единственный разрешённый импорт проекта).
Формат схемы: список (имя, (кубиты,)) с именами из GATES. Точные правила взяты из проекта
(block_synth.commutes / global_cancel.py / proto2/phase_fold.py) и продублированы, чтобы песочница не тянула torch."""
GATES = {"h", "x", "cx", "t", "tdg", "s", "sdg"}
INV = {"h": "h", "x": "x", "s": "sdg", "sdg": "s", "t": "tdg", "tdg": "t", "cx": "cx"}
DIAG = {"s", "sdg", "t", "tdg"}
PH = {"t": 1, "s": 2, "sdg": 6, "tdg": 7}
SEQ = {0: [], 1: ["t"], 2: ["s"], 3: ["s", "t"], 4: ["s", "s"], 5: ["sdg", "tdg"], 6: ["sdg"], 7: ["tdg"]}
SELF = {"cx", "x", "h"}


def cx(c):
    return sum(g[0] == "cx" for g in c)


def tcount(c):
    return sum(g[0] in ("t", "tdg") for g in c)


def key(c):
    return (cx(c), tcount(c), len(c))


def inv_circ(c):
    return [(INV[nm], tuple(qb)) for nm, qb in reversed(c)]


def commutes(g1, g2):
    (n1, q1), (n2, q2) = g1, g2
    if not set(q1) & set(q2):
        return True
    if n1 != "cx" and n2 != "cx":
        return (n1 in DIAG and n2 in DIAG) or (n1 == n2 == "x") or (n1 == n2 == "h")
    if n1 == "cx" and n2 == "cx":
        return q1[0] != q2[1] and q2[0] != q1[1]
    one, c = (g1, g2) if n2 == "cx" else (g2, g1)
    return (one[0] in DIAG and one[1][0] == c[1][0]) or (one[0] == "x" and one[1][0] == c[1][1])


def _sweep(c):
    c = list(c); i = 0; changed = False
    while i < len(c):
        g = c[i]; j = i + 1; done = False
        while j < len(c):
            h = c[j]
            if not (set(g[1]) & set(h[1])):
                j += 1; continue
            if g[0] in SELF and h == g:
                del c[j]; del c[i]; changed = done = True; break
            if g[0] in DIAG and h[0] in DIAG and g[1] == h[1]:
                new = [(nm, g[1]) for nm in SEQ[(PH[g[0]] + PH[h[0]]) % 8]]
                if len(new) < 2:
                    del c[j]; c[i:i + 1] = new; changed = done = True; break
            if commutes(g, h):
                j += 1; continue
            break
        if not done:
            i += 1
    return c, changed


def global_cancel(c):
    c = [(g[0], tuple(g[1])) for g in c]; changed = True
    while changed:
        c, changed = _sweep(c)
    return c


def phase_fold(c, n):
    """Свёртка фаз по аффинным чётностям (сумма по путям); точно до общей фазы."""
    PHZ = dict(PH, z=4)
    wire = [(1 << q, 0) for q in range(n)]; nv = n; first, acc, out = {}, {}, []
    for nm, qs in c:
        if nm in PHZ:
            m, k = wire[qs[0]]; v = PHZ[nm] * (-1 if k else 1) % 8
            if m in first: acc[m] = (acc[m] + v) % 8
            else: first[m] = len(out); acc[m] = v; out.append(("PH", m))
            continue
        if nm == "h": wire[qs[0]] = (1 << nv, 0); nv += 1
        elif nm == "x": m, k = wire[qs[0]]; wire[qs[0]] = (m, k ^ 1)
        elif nm == "cx":
            a, b = qs; wire[b] = (wire[b][0] ^ wire[a][0], wire[b][1] ^ wire[a][1])
        out.append((nm, qs))
    wire = [(1 << q, 0) for q in range(n)]; nv = n; res = []
    for nm, qs in out:
        if nm == "PH":
            q = next(i for i in range(n) if wire[i][0] == qs)
            res += [(s, (q,)) for s in SEQ[acc[qs] * (-1 if wire[q][1] else 1) % 8]]
            continue
        if nm == "h": wire[qs[0]] = (1 << nv, 0); nv += 1
        elif nm == "x": m, k = wire[qs[0]]; wire[qs[0]] = (m, k ^ 1)
        elif nm == "cx":
            a, b = qs; wire[b] = (wire[b][0] ^ wire[a][0], wire[b][1] ^ wire[a][1])
        res.append((nm, qs))
    return res


def fold_cancel(c, n, iters=4):
    best = global_cancel(c)
    for _ in range(iters):
        c2 = global_cancel(phase_fold(best, n))
        if key(c2) >= key(best):
            break
        best = c2
    return best


# ---- правила «шаблон -> короче» (формат для кандидатов вида rules) -----------------------
def take_block(c, i, support):
    """С позиции i собрать подряд идущие команды внутри support, перепрыгивая непересекающиеся;
    остановиться на команде, частично выходящей за support. Возвращает индексы."""
    S = set(support); idx = []; blocked = set()
    for j in range(i, len(c)):
        qs = set(c[j][1])
        if not qs & S:
            continue
        if qs <= S and not qs & blocked:
            idx.append(j)
        else:
            blocked |= qs & S
            if blocked == S:
                break
    return idx


def match_rule(c, i, lhs, k):
    """Сопоставить lhs с c начиная с i. Команды, не задевающие уже связанные кубиты, пропускаются
    (они не пересекаются с ранее сопоставленными => переносить можно). Команда на связанном кубите
    обязана совпасть со следующей командой lhs. Возвращает (позиции, отображение) или None."""
    mp, used, pos, t = {}, set(), [], 0
    for j in range(i, len(c)):
        if t == len(lhs):
            break
        nm, qs = c[j]; g_nm, g_qs = lhs[t]
        touches = bool(set(qs) & used)
        if nm == g_nm and len(qs) == len(g_qs):
            ok = True; m2 = dict(mp)
            for a, b in zip(g_qs, qs):
                if a in m2 and m2[a] != b: ok = False; break
                if a not in m2 and b in m2.values(): ok = False; break
                m2[a] = b
            if ok:
                mp = m2; used = set(mp.values()); pos.append(j); t += 1; continue
        if touches or t == 0:
            return None
    if t != len(lhs) or len(mp) != k:
        return None
    return pos, mp


def apply_rules(c, rules, max_rounds=3):
    """rules: [{"k":k, "lhs":[[g,[q..]],..], "rhs":[...]}], локальные кубиты 0..k-1; правило проверено точно
    (равенство lhs и rhs как операций на k кубитах). Замена ставится на место последней сопоставленной команды."""
    c = [(g[0], tuple(g[1])) for g in c]
    rules = [r for r in rules if sum(g[0] == "cx" for g in r["rhs"]) < sum(g[0] == "cx" for g in r["lhs"])]
    for _ in range(max_rounds):
        changed = False; i = 0
        while i < len(c):
            for r in rules:
                m = match_rule(c, i, r["lhs"], r["k"])
                if not m:
                    continue
                pos, mp = m; P = set(pos)
                new = [(g[0], tuple(mp[q] for q in g[1])) for g in r["rhs"]]
                c = c[:pos[0]] + [c[j] for j in range(pos[0], pos[-1] + 1) if j not in P] + new + c[pos[-1] + 1:]
                changed = True; break
            i += 1
        if not changed:
            break
    return c


# ---- уровень Тоффоли (для кандидатов вида 'lower'): имена ccx, ccz, cz, cs, csdg, z + GATES -----------
HL_DIAG = {"z", "cz", "ccz", "t", "tdg", "s", "sdg", "cs", "csdg"}


def ccx_std(a, b, t):
    return [("h", (t,)), ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("t", (t,)), ("cx", (b, t)),
            ("tdg", (t,)), ("cx", (a, t)), ("t", (b,)), ("t", (t,)), ("h", (t,)), ("cx", (a, b)),
            ("t", (a,)), ("tdg", (b,)), ("cx", (a, b))]


def lower_gate(nm, q):
    """Стандартное разложение одной команды уровня Тоффоли в Clifford+T."""
    q = tuple(q)
    if nm == "ccx": return ccx_std(*q)
    if nm == "ccz": return [("h", (q[2],))] + ccx_std(*q) + [("h", (q[2],))]
    if nm == "cz": return [("h", (q[1],)), ("cx", q), ("h", (q[1],))]
    if nm in ("cs", "csdg"):
        a, b = q; s, sd = ("t", "tdg") if nm == "cs" else ("tdg", "t")
        return [(s, (a,)), (s, (b,)), ("cx", (a, b)), (sd, (b,)), ("cx", (a, b))]
    if nm == "z": return [("s", q), ("s", q)]
    return [(nm, q)]


def lower_std(hl):
    out = []
    for nm, q in hl: out += lower_gate(nm, q)
    return out


# ---- выбор целей CCZ (fold_h / expand) и combo2 — перенос из proto2/rtof.py + run_combo2.py ----------------
# В реальных схемах (и в наборах val после v2) Тоффоли записаны как H·CCZ·H с сокращёнными H; цель CCX
# не задана. hl_target_default(hl) — эвристика fold_h (цель там, где на проводе по обе стороны H);
# combo2 — локальный поиск целей CCZ и вариантов непарных CCX по итоговому CX после fold_cancel.
def _wire_nb(g, i, t, d):
    k = i + d
    while 0 <= k < len(g):
        if g[k] and t in g[k][1]:
            return k
        k += d
    return None


def cancel_h(g):
    """H_t H_t, соседние на проводе t (между ними нет команд на t) -> убрать. Точно."""
    g = [(x[0], tuple(x[1])) for x in g]; last = {}
    for i, x in enumerate(g):
        if x is None:
            continue
        nm, qs = x
        if nm == "h":
            t = qs[0]
            if t in last and g[last[t]] == ("h", (t,)):
                g[last[t]] = g[i] = None; last.pop(t); continue
        for q in qs:
            last[q] = i
    return [x for x in g if x]


def fold_h_choice(g):
    """Эвристика выбора цели для каждого CCZ/CZ: {индекс: цель} (как rtof.fold_h)."""
    g = [(x[0], tuple(x[1])) for x in g]; choice = {}
    for i, (nm, qs) in enumerate(g):
        if nm not in ("ccz", "cz"):
            continue
        best, bs = None, (0, 0)
        for t in qs:
            sc = 0; dist = 0
            for d in (-1, 1):
                k = _wire_nb(g, i, t, d)
                while k is not None and g[k][0] in ("ccz", "cz") and t in g[k][1] and set(g[k][1]) == set(qs):
                    k = _wire_nb(g, k, t, d)
                if k is not None and g[k] == ("h", (t,)):
                    sc += 1; dist -= abs(k - i)
            if (sc, dist) > bs:
                best, bs = t, (sc, dist)
        if best is not None and bs[0] == 2:
            choice[i] = best
    return choice


def expand(g, choice):
    """CCZ(i) с целью choice[i] -> H·CCX·H (CZ -> H·CX·H), затем cancel_h. choice[i]=None — оставить CCZ."""
    out = []
    for i, x in enumerate(g):
        nm, qs = x[0], tuple(x[1]); t = choice.get(i)
        if t is None or nm not in ("ccz", "cz"):
            out.append((nm, qs))
        else:
            ctrl = tuple(q for q in qs if q != t)
            out += [("h", (t,)), ("ccx" if nm == "ccz" else "cx", ctrl + (t,)), ("h", (t,))]
    return cancel_h(out)


def fold_h(g):
    """Уровень Тоффоли с целями по эвристике (то, что получают кандидаты KIND='lower' без HL_RAW)."""
    return expand(g, fold_h_choice(g))


def diag_on(gate, q):
    """Команда уровня Тоффоли действует на кубит q диагонально?"""
    nm, qs = gate
    if q not in qs or nm in HL_DIAG:
        return True
    if nm in ("cx", "ccx"):
        return q != qs[-1]
    return False


def pair_toffolis(g, use4=True):
    """Пары CCX с общими управляющими -> Тоффоли с относительной фазой: 'R' (3 CX, та же цель, диагонально на
    {a,b,t} между ними) или 'R4' (4 CX, D=CS^†(a,b), цель любая). Возврат {индекс: (вид, сопряжённый?, кубиты)}."""
    mark = {}
    for i, (nm, qs) in enumerate(g):
        if nm != "ccx" or i in mark:
            continue
        C, S = set(qs[:2]), set(qs); ok3 = True
        for j in range(i + 1, len(g)):
            h = g[j]
            if h[0] == "ccx" and set(h[1][:2]) == C and j not in mark:
                if ok3 and h[1][-1] == qs[-1]:
                    mark[i], mark[j] = ("R", 0, tuple(qs)), ("R", 1, tuple(qs)); break
                if use4:
                    mark[i], mark[j] = ("R4", 0, tuple(qs)), ("R4", 1, tuple(qs[:2]) + (h[1][2],)); break
            ok3 = ok3 and all(diag_on(h, q) for q in S)
            if not all(diag_on(h, q) for q in C) or (not use4 and not ok3):
                break
    return mark


def rccx(a, b, t):
    """= D·Tof, D диагональна (3 CX)."""
    return [("h", (t,)), ("t", (t,)), ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("t", (t,)),
            ("cx", (b, t)), ("tdg", (t,)), ("h", (t,))]


def rccx4(a, b, t):
    """= CS^†(a,b)·Tof (4 CX), фаза не зависит от цели."""
    return rccx(a, b, t)[:-1] + [("cx", (a, t)), ("h", (t,))]


def lower_marked(g, mark=None, var=None):
    """Разложение уровня Тоффоли: пары из mark -> rccx/rccx4 (и сопряжённые); непарные CCX — вариант
    var[i]=(swap управляющих, обращённое разложение); остальное — lower_gate."""
    mark = mark or {}; var = var or {}; out = []
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


def combo2(hl, n, rounds=3, max_evals=None, use4=True):
    """Сильнейший известный метод (proto2/run_combo2.py): hl в форме CCZ+H (цели не выбраны) ->
    локальный поиск целей CCZ/CZ (fold_h + перебор {нет, q1, q2, q3}) и вариантов непарных CCX
    (порядок управляющих × обращённое разложение) по ключу (CX, T, длина) после пар с относительной фазой
    и fold_cancel. max_evals — потолок числа оценок (детерминированная замена тайм-аута)."""
    g = [(x[0], tuple(x[1])) for x in hl]; var = {}; cnt = [0]

    def full(g2):
        return fold_cancel(lower_marked(g2, pair_toffolis(g2, use4), var), n)

    def sc(g2):
        cnt[0] += 1; return key(full(g2))

    def spent():
        return max_evals is not None and cnt[0] >= max_evals

    choice = fold_h_choice(g)
    idx = [i for i, x in enumerate(g) if x[0] in ("ccz", "cz")]
    best = sc(expand(g, choice))
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
        g2 = expand(g, choice); m = pair_toffolis(g2, use4)
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
