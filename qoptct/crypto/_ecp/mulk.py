# SPDX-License-Identifier: Apache-2.0
"""k4/mulk.py — умножение в модели «И + снятие измерением» для схем «вычислить-применить-обратить» (мусор допускается,
его снимает обращение всей схемы):
  add_gen(C, Y, t)            t += Y (Y — список операндов W(w)/ZERO/ONE, len(Y) <= len(t)), mod 2^len(t); И = len(t)-1 - свободные
  school_mul(C, a, b, out)    out(нулевой, len a+len b) = a*b; И = 2·|a|·|b| − |a|
  school_sq(C, a, out)        out = a^2;  И ≈ m^2
  kara_mul / kara_sq          Карацуба с мусором (уровни L), база — школьное
  red_special(C, V, p)        V (2m бит) -> V mod p для p специального вида (Солинас/псевдо-Мерсенн): O(m) И, мусор
"""
from .engine import *


def add_gen(C, Y, b):
    """b += Y mod 2^len(b); Y[i] — операнд ('w',w) / ZERO / ONE; недостающие старшие — ZERO."""
    m = len(b); Y = list(Y) + [ZERO] * (m - len(Y))
    X = [W(w) for w in b]
    c = [ZERO]; R = [None]
    for i in range(m - 1):
        r = plan(X[i], Y[i], c[i]); c.append(compute(C, r, fold_ok=True)); R.append(r)
    for i in reversed(range(m)):
        if i + 1 < m and R[i + 1].kind != 'op': uncompute(C, R[i + 1])
        if Y[i] == ONE: C.x(b[i])
        elif Y[i][0] == 'w': C.cx(Y[i][1], b[i])
        folded = (i + 1 < m) and R[i + 1].kind == 'gid' and R[i + 1].fold
        if not folded:
            if c[i] == ONE: C.x(b[i])
            elif c[i][0] == 'w': C.cx(c[i][1], b[i])


def sub_gen(C, Y, b):
    """b -= Y mod 2^len(b)."""
    for w in b: C.x(w)
    add_gen(C, Y, b)
    for w in b: C.x(w)


def school_mul(C, a, b, out):
    na, nb = len(a), len(b)
    assert len(out) >= na + nb
    for i in range(nb):
        if i == 0:
            for j in range(na):
                if a[j] == b[0]: C.cx(a[j], out[j])
                else: C.AND(b[0], a[j], out[j])
            continue
        m = [C.AND(b[i], w) if w != b[i] else None for w in a]
        Y = [W(mw) if mw is not None else W(a[j]) for j, mw in enumerate(m)]
        add_gen(C, Y, out[i:i + na + 1])             # частичная сумма a*(b mod 2^(i+1)) < 2^(na+i+1)
        for j in reversed(range(na)):
            if m[j] is not None: C.erase(m[j], b[i], a[j])


def school_sq(C, a, out):
    """out = a^2: строка i = a_i·2^(2i) + Σ_{j>i} 2 a_i a_j 2^(i+j); сумма строк <= i  < 2^(m+i+2)."""
    m = len(a)
    for i in range(m):
        prods = [C.AND(a[i], a[j]) for j in range(i + 1, m)]
        Y = [W(a[i]), ZERO] + [W(t) for t in prods]
        hi = min(len(out), m + i + 2)
        add_gen(C, Y, out[2 * i:hi])
        for j, t in reversed(list(zip(range(i + 1, m), prods))): C.erase(t, a[i], a[j])


def _half_sum(C, x0, x1):
    """s = x0 + x1 в новом регистре (max+1 бит). И = len(x0) (или len(x1))."""
    L = max(len(x0), len(x1)) + 1
    s = [C.alloc() for _ in range(L)]
    lo, hi = (x0, x1) if len(x0) >= len(x1) else (x1, x0)
    for w, t in zip(lo, s): C.cx(w, t)
    add_gen(C, [W(w) for w in hi], s)
    return s


def kara_mul(C, a, b, out, L=1, base=8, gar=None):
    """out(нулевой, >= |a|+|b|) = a*b; |a| == |b| == m. Возвращает список мусорных проводов."""
    if gar is None: gar = []
    m = len(a); assert len(b) == m
    if L == 0 or m <= base:
        school_mul(C, a, b, out); return gar
    h = m // 2
    a0, a1, b0, b1 = a[:h], a[h:], b[:h], b[h:]
    kara_mul(C, a0, b0, out[:2 * h], L - 1, base, gar)
    kara_mul(C, a1, b1, out[2 * h:2 * m], L - 1, base, gar)
    sa = _half_sum(C, a0, a1); sb = _half_sum(C, b0, b1)
    T = [C.alloc() for _ in range(2 * len(sa))]
    kara_mul(C, sa, sb, T, L - 1, base, gar)
    sub_gen(C, [W(w) for w in out[:2 * h]], T)
    sub_gen(C, [W(w) for w in out[2 * h:2 * m]], T)
    # T = a0 b1 + a1 b0 < 2^(m+1): старшие биты T равны 0
    add_gen(C, [W(w) for w in T[:m + 1]], out[h:2 * m])
    gar += sa + sb + T
    return gar


def kara_sq(C, a, out, L=1, base=8, gar=None):
    if gar is None: gar = []
    m = len(a)
    if L == 0 or m <= base:
        school_sq(C, a, out); return gar
    h = m // 2
    a0, a1 = a[:h], a[h:]
    kara_sq(C, a0, out[:2 * h], L - 1, base, gar)
    kara_sq(C, a1, out[2 * h:2 * m], L - 1, base, gar)
    sa = _half_sum(C, a0, a1)
    T = [C.alloc() for _ in range(2 * len(sa))]
    kara_sq(C, sa, T, L - 1, base, gar)
    sub_gen(C, [W(w) for w in out[:2 * h]], T)
    sub_gen(C, [W(w) for w in out[2 * h:2 * m]], T)
    add_gen(C, [W(w) for w in T[:m + 1]], out[h:2 * m])
    gar += sa + T
    return gar
