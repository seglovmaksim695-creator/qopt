"""Ресинтез CX-фазовой области: фазовый полином + аффинное отображение -> схема с минимальным (эвристически) числом CX."""
import random, numpy as np
PH = {"t": 1, "s": 2, "sdg": 6, "tdg": 7}
SEQ = {0: [], 1: ["t"], 2: ["s"], 3: ["s", "t"], 4: ["s", "s"], 5: ["sdg", "tdg"], 6: ["sdg"], 7: ["tdg"]}

def poly_of(gates, qlist):
    """gates на глобальных кубитах -> (m, terms{mask:angle}, wires[(mask,const)])"""
    idx = {q: i for i, q in enumerate(qlist)}; m = len(qlist)
    wire = [(1 << i, 0) for i in range(m)]; terms = {}
    for nm, qs in gates:
        if nm == "cx":
            a, b = idx[qs[0]], idx[qs[1]]; wire[b] = (wire[b][0] ^ wire[a][0], wire[b][1] ^ wire[a][1])
        elif nm == "x":
            i = idx[qs[0]]; wire[i] = (wire[i][0], wire[i][1] ^ 1)
        else:
            i = idx[qs[0]]; mk, k = wire[i]; terms[mk] = (terms.get(mk, 0) + PH[nm] * (-1 if k else 1)) % 8
    terms = {k: v for k, v in terms.items() if v}
    return m, terms, wire

def apply_circ(local, m):
    """локальная схема (имена, локальные индексы) -> (terms, wires) для сверки"""
    wire = [(1 << i, 0) for i in range(m)]; terms = {}
    for nm, qs in local:
        if nm == "cx":
            a, b = qs; wire[b] = (wire[b][0] ^ wire[a][0], wire[b][1] ^ wire[a][1])
        elif nm == "x":
            i = qs[0]; wire[i] = (wire[i][0], wire[i][1] ^ 1)
        else:
            i = qs[0]; mk, k = wire[i]; terms[mk] = (terms.get(mk, 0) + PH[nm] * (-1 if k else 1)) % 8
    return {k: v for k, v in terms.items() if v}, wire

# ---------- линейный синтез: перевести состояние S (строки-маски) в целевое Tg операциями row_b ^= row_a
def _inv(rows, m):
    A = list(rows); I = [1 << i for i in range(m)]
    for c in range(m):
        p = next(r for r in range(c, m) if A[r] >> c & 1)
        A[c], A[p] = A[p], A[c]; I[c], I[p] = I[p], I[c]
        for r in range(m):
            if r != c and A[r] >> c & 1: A[r] ^= A[c]; I[r] ^= I[c]
    return I  # обратная матрица (строки)

def _mul(Aa, Bb, m):  # (A*B) строки: row_i(A*B)=xor of rows of B selected by A_i
    out = []
    for r in Aa:
        v = 0; j = 0
        while r:
            if r & 1: v ^= Bb[j]
            r >>= 1; j += 1
        out.append(v)
    return out

def _gauss(Q, m, rng):
    """Q -> I рядными операциями; список (a,b) (row_b ^= row_a) в порядке применения к Q."""
    A = list(Q); ops = []
    cols = list(range(m)); 
    if rng: rng.shuffle(cols)
    # приводим к единичной с произвольным порядком столбцов нельзя (нужна диагональ) -> порядок только пивотов
    for c in (cols if rng else range(m)):
        # нужен единичный бит на диагонали c; ищем строку r>=? среди неиспользованных строк. Упрощение: полная обработка по диагонали
        pass
    return None

def lin_synth(S, Tg, m, rng=None, tries=1):
    """Список CX (a,b) последовательно, переводящий S в Tg."""
    Q = _mul(Tg, _inv(S, m), m)  # E*S=Tg => E=Tg*S^-1
    best = None
    for t in range(tries):
        ops = (_elim(Q, m, rng if t else None) if t % 3 == 0 else _pmh(Q, m, rng) if t % 3 == 1 else _greedy(Q, m, rng))
        if best is None or len(ops) < len(best): best = ops
    # ops приводят Q к I: L_j..L_1 Q = I => порядок схемы: reversed(ops)
    return list(reversed(best))

def _elim(Q, m, rng):
    A = list(Q); ops = []
    order = list(range(m))
    if rng: rng.shuffle(order)
    # гауссово исключение: для каждого столбца c (в порядке) ставим пивот в строку c
    for c in order:
        if not (A[c] >> c & 1):
            cand = [r for r in range(m) if r != c and A[r] >> c & 1 and (r not in order[:order.index(c)])]
            if not cand: cand = [r for r in range(m) if r != c and A[r] >> c & 1]
            r = rng.choice(cand) if rng else cand[0]
            A[c] ^= A[r]; ops.append((r, c))
        for r in range(m):
            if r != c and A[r] >> c & 1: A[r] ^= A[c]; ops.append((c, r))
    # проверка: A должна быть I
    return ops if A == [1 << i for i in range(m)] else _elim_fallback(Q, m)

def _elim_fallback(Q, m):
    A = list(Q); ops = []
    for c in range(m):
        if not (A[c] >> c & 1):
            r = next(r for r in range(c + 1, m) if A[r] >> c & 1); A[c] ^= A[r]; ops.append((r, c))
        for r in range(m):
            if r != c and A[r] >> c & 1: A[r] ^= A[c]; ops.append((c, r))
    return ops

def _greedy(Q, m, rng):
    """жадный спуск по суммарному весу строк; тупик -> достраиваем гауссом"""
    A = list(Q); ops = []; I = [1 << i for i in range(m)]
    pc = lambda x: bin(x).count("1")
    while A != I and len(ops) < 4 * m * m:
        bestd = 0; cands = []
        for a in range(m):
            for b in range(m):
                if a == b: continue
                d = pc(A[b]) - pc(A[b] ^ A[a])
                # учтём близость к единичной: расстояние до e_b
                d += pc(A[b] ^ I[b]) - pc(A[b] ^ A[a] ^ I[b])
                if d > bestd: bestd = d; cands = [(a, b)]
                elif d == bestd and d > 0: cands.append((a, b))
        if not cands: break
        a, b = rng.choice(cands) if rng else cands[0]
        A[b] ^= A[a]; ops.append((a, b))
    if A != I:
        # достроить: гаусс на текущей матрице
        tail = _elim_fallback(A, m)
        ops += tail
    return ops

def _pmh(Q, m, rng):
    """PMH: нижнетреугольная + верхнетреугольная часть, секции по sz столбцов."""
    sz = max(1, min(m, rng.choice([1, 2, 3]) if rng else 2))
    A = list(Q); ops = []
    def lower(A):
        for s0 in range(0, m, sz):
            cols = range(s0, min(m, s0 + sz)); seen = {}
            for r in range(s0, m):
                pat = sum(((A[r] >> c) & 1) << (c - s0) for c in cols)
                if pat == 0: continue
                if pat in seen: A[r] ^= A[seen[pat]]; ops.append((seen[pat], r))
                else: seen[pat] = r
            for c in cols:
                if not (A[c] >> c & 1):
                    cand = [r for r in range(c + 1, m) if A[r] >> c & 1]
                    if not cand: return False
                    A[c] ^= A[cand[0]]; ops.append((cand[0], c))
                for r in range(c + 1, m):
                    if A[r] >> c & 1: A[r] ^= A[c]; ops.append((c, r))
        return True
    if not lower(A): return _elim_fallback(Q, m)
    # верхняя часть: транспонирование
    def transpose(A):
        return [sum(((A[j] >> i) & 1) << j for j in range(m)) for i in range(m)]
    # после lower A верхнетреугольная с единичной диагональю; приведём построчным исключением (m^2/2 в худшем случае)
    for c in range(m - 1, -1, -1):
        for r in range(c):
            if A[r] >> c & 1: A[r] ^= A[c]; ops.append((c, r))
    return ops if A == [1 << i for i in range(m)] else _elim_fallback(Q, m)

# ---------- основной ресинтез
def _dist_table(S, m):
    T = np.full(1 << m, 99, dtype=np.int16); T[0] = 0
    ar = np.arange(1 << m)
    for w in S:
        T = np.minimum(T, T[ar ^ w] + 1)
    return T

def synth(m, terms, wires, rng, beam=4, alpha=1.0, noise=0.0, lin_tries=6, limit=10**9):
    """-> локальная схема [(имя,(локальные индексы))] или None. Целевая: terms (маска->угол), wires (маска,const)."""
    Tg = [w[0] for w in wires]; consts = [w[1] for w in wires]
    if m > 14: return None
    S0 = tuple(1 << i for i in range(m))
    def emit(S, rem, ops):
        rem = dict(rem); ops = list(ops)
        for q, w in enumerate(S):
            a = rem.pop(w, None)
            if a is not None: ops += [(nm, (q,)) for nm in SEQ[a]]
        return rem, ops
    rem0, ops0 = emit(S0, terms, [])
    states = [(0, S0, rem0, ops0)]; finished = []; best_tot = None
    seen = set()
    while states:
        nxt = []
        for cxn, S, rem, ops in states:
            if not rem:
                finished.append((cxn, S, ops)); continue
            for a in range(m):
                for b in range(m):
                    if a == b: continue
                    S2 = list(S); S2[b] ^= S[a]; S2 = tuple(S2)
                    key = (S2, frozenset(rem))
                    rem2, ops2 = emit(S2, rem, ops + [("cx", (a, b))])
                    T = _dist_table(S2, m)
                    h = sum(int(T[p]) - 1 for p in rem2)
                    sc = cxn + 1 + alpha * h + (rng.random() * noise if noise else 0) - (len(rem) - len(rem2)) * 0.0
                    nxt.append((sc, cxn + 1, S2, rem2, ops2))
        # отбор: лучшие beam различных состояний
        nxt.sort(key=lambda x: x[0]); states = []; seen = set()
        for sc, cxn, S2, rem2, ops2 in nxt:
            k = (S2, frozenset(rem2))
            if k in seen: continue
            seen.add(k); states.append((cxn, S2, rem2, ops2))
            if len(states) >= beam: break
        if best_tot is not None:
            states = [s for s in states if s[0] < best_tot]
        for st in states:
            if not st[2]:
                pass
        # завершённые считаем сразу
        fin = [s for s in states if not s[2]]
        for cxn, S, rem, ops in fin:
            lin = lin_synth(list(S), Tg, m, rng, lin_tries)
            tot = cxn + len(lin)
            if best_tot is None or tot < best_tot:
                best_tot = tot; best = ops + [("cx", p) for p in lin]
        states = [s for s in states if s[2]]
        lim = best_tot if best_tot is not None else limit
        states = [s for s in states if s[0] + 1 < lim]
    if best_tot is None and limit < 10**9:
        return None
    if best_tot is None:
        # термов нет вовсе
        lin = lin_synth(list(S0), Tg, m, rng, lin_tries); best = ops0 + [("cx", p) for p in lin]
    best = best + [("x", (q,)) for q in range(m) if consts[q]]
    return best
