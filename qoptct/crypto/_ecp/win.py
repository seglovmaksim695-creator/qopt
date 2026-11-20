# SPDX-License-Identifier: Apache-2.0
"""k4/padd_win.py — ОКОННОЕ полное сложение точек для алгоритма Шора: P <- P + T[k], k — квантовый индекс окна (w бит),
T[k] — классическая таблица точек (T[k] != O: смещение T[k] = k·G + R0, R0 — классическая точка; смещение учитывается
классически вне схемы). P — любая точка или O=(0,0).
Шаги: (1) табличный поиск (QROM) значений для k: x_Q, y_Q, −y_Q, 3x_Q, x_2Q, y_2Q, −y_2Q, x_{D+Q}, y_{D+Q};
(2) флаги краевых случаев (сравнение с квантовыми значениями); (3) подстановка безопасной точки D (классическая, безопасная
для ВСЕХ k); (4) общее сложение с квантовой Q; (5) исправление образа; (6) снятие флагов; (7) снятие поиска.
Поиск: минтермы (2^w − 2 И), снятие — повторным XOR той же таблицы (ещё 2^w − 2 И); в оценке n=256 используется
снятие измерением a4 (2^h+2^l−w−2 И), отдельно."""
from .engine import *
from .ec import ec_add, ec_neg
from .full import eq_chain, lit_and, lit_erase, xor_const
from . import padd as padd_k4
from .padd import mod_neg


def _minterms(C, g):
    """вернуть (leaves[j] = [idx==j], undo_fn). И = 2^len(g) − 2."""
    one = C.alloc(); C.x(one)
    u0 = C.alloc(); C.cx(g[0], u0); C.cx(u0, one)          # one -> ~g0
    terms = [one, u0]                                       # индекс j: бит0 = 0,1
    log = []
    for i in range(1, len(g)):
        new = [None] * (2 * len(terms))
        for j, t in enumerate(terms):
            u = C.AND(t, g[i]); C.cx(u, t)                  # t <- t & ~g_i
            new[j] = t; new[j + len(terms)] = u; log.append((t, u, g[i]))
        terms = new
    def undo():
        for t, u, gi in reversed(log):
            C.cx(u, t); C.erase(u, t, gi)
        C.cx(u0, one); C.cx(g[0], u0); C.release(u0); C.x(one); C.release(one)
    return terms, undo


FIELDS = ['xq', 'yq', 'ymq', 'x3q', 'x2q', 'y2q', 'ym2q', 'xdq', 'ydq']


def table_vals(Tk, D, cv):
    p, a = cv['p'], cv['a']
    out = []
    for Q in Tk:
        Q2 = ec_add(Q, Q, p, a); DQ = ec_add(D, Q, p, a)
        out.append(dict(xq=Q[0], yq=Q[1], ymq=(-Q[1]) % p, x3q=3 * Q[0] % p, x2q=Q2[0], y2q=Q2[1], ym2q=(-Q2[1]) % p,
                        xdq=DQ[0], ydq=DQ[1]))
    return out


def lookup(C, Kw, vals, regs):
    leaves, undo = _minterms(C, Kw)
    for j, v in enumerate(vals):
        for f in FIELDS:
            for i, w in enumerate(regs[f]):
                if (v[f] >> i) & 1: C.cx(leaves[j], w)
    undo()


def and_xor(C, f, src, dst):
    """dst ^= f·src (n И, стирание бесплатно)."""
    for s, d in zip(src, dst):
        t = C.AND(f, s); C.cx(t, d); C.erase(t, f, s)


def eq_reg_chain(C, X, Q, guard=None):
    """g = [X == Q] (Q — квантовый регистр): X ^= Q, цепочка [X == 0]; если guard — g ∧ ¬guard. возвращает (g, undo)."""
    for a, b in zip(Q, X): C.cx(a, b)
    g, und = eq_chain(C, X, 0)
    if guard is None:
        def u1():
            und()
            for a, b in zip(Q, X): C.cx(a, b)
        return g, u1
    g2 = lit_and(C, g, guard, 0)
    def u2():
        lit_erase(C, g2, g, guard, 0); und()
        for a, b in zip(Q, X): C.cx(a, b)
    return g2, u2


def parity_lit(C, g, Y0, q0, erase_target=None):
    """AND(g, [Y0 == q0]) (q0 — провод): Y0 ^= q0, литерал Y0==0."""
    C.cx(q0, Y0)
    if erase_target is None: t = lit_and(C, g, Y0, 0)
    else: lit_erase(C, erase_target, g, Y0, 0); t = None
    C.cx(q0, Y0)
    return t


def generic_q(C, X, Y, R, p, L):
    """общее смешанное сложение с КВАНТОВОЙ Q (регистры R[...]): как padd_k4.make, но константы -> регистры."""
    n = len(X)
    c2 = pow(2, 2 * n, p)
    from .padd import const_mul_mod, pmul, kaliski_forward, RED
    def prod_apply(a, b, apply):
        def fn():
            out, gar = pmul(C, a, b, p, L); return out, list(out) + gar
        compute_apply_uncompute(C, fn, apply)
    def div_apply(num, den, apply):
        def fn():
            xs, g0 = const_mul_mod(C, den, c2, p)
            A, B, K, S, g1 = kaliski_forward(C, xs, p, shift=True, trunc=True, halfB=True, red=RED['mode'])
            P, g2 = pmul(C, num, K, p, L)
            live = list(dict.fromkeys(list(P) + g2 + list(K) + g1 + list(S) + list(A) + list(B) + list(xs) + g0))
            dset = set(den) | set(num)
            return P, [w for w in live if w not in dset]
        compute_apply_uncompute(C, fn, apply)
    xq, yq, ymq, x3q = R['xq'], R['yq'], R['ymq'], R['x3q']
    mod_sub(C, xq, X, p); mod_sub(C, yq, Y, p)                         # dx, dy
    Lr = [C.alloc() for _ in range(n)]
    def xor_into(T):
        def f(P):
            for a_, b_ in zip(P, T): C.cx(a_, b_)
        return f
    div_apply(Y, X, xor_into(Lr))                                       # L = −λ
    prod_apply(Lr, X, lambda P: mod_add(C, P, Y, p))                    # Y = 0
    prod_apply(Lr, Lr, lambda P: mod_rsub_(C, P, X, p))                 # X = λ² − dx
    mod_sub(C, x3q, X, p)                                               # X = x' − xq
    prod_apply(Lr, X, lambda P: mod_add(C, P, Y, p))                    # Y = λ(xq − x') = y' + yq
    mod_neg(C, X, p)                                                    # X = xq − x'
    div_apply(Y, X, xor_into(Lr))                                       # L -> 0
    mod_neg(C, X, p); mod_add(C, xq, X, p); mod_add(C, ymq, Y, p)       # (x', y')
    for w in reversed(Lr): C.release(w)


def mod_rsub_(C, a, b, p):
    mod_neg(C, b, p); mod_add(C, a, b, p)


def point_add_window(C, X, Y, Kw, Tk, cv, D, L=1):
    p = cv['p']; n = len(X)
    vals = table_vals(Tk, D, cv)
    R = {f: [C.alloc() for _ in range(n)] for f in FIELDS}
    lookup(C, Kw, vals, R)
    # --- флаги
    fO, undO = eq_chain(C, Y, 0); undO(erase_last=False)
    g, und = eq_reg_chain(C, X, R['xq'], guard=fO)
    fQ = parity_lit(C, g, Y[0], R['yq'][0])
    fmQ = C.alloc(); C.cx(g, fmQ); C.cx(fQ, fmQ)
    und()
    g, und = eq_reg_chain(C, X, R['x2q'], guard=fO)
    fm2Q = parity_lit(C, g, Y[0], R['ym2q'][0])
    und()
    # --- подстановка D: P ^= f·(E ⊕ D)
    Dx, Dy = D
    xor_const(C, fO, X, Dx); xor_const(C, fO, Y, Dy)                     # E_O = (0,0)
    gq = C.alloc(); C.cx(fQ, gq); C.cx(fmQ, gq)                          # gq = fQ ∨ fmQ
    and_xor(C, gq, R['xq'], X); xor_const(C, gq, X, Dx)
    and_xor(C, fQ, R['yq'], Y); and_xor(C, fmQ, R['ymq'], Y); xor_const(C, gq, Y, Dy)
    and_xor(C, fm2Q, R['x2q'], X); and_xor(C, fm2Q, R['ym2q'], Y)
    xor_const(C, fm2Q, X, Dx); xor_const(C, fm2Q, Y, Dy)
    # --- общее сложение: D -> D + Q_k
    generic_q(C, X, Y, R, p, L)
    # --- исправление: R ^= F·(D+Q) ; затем образы: O->Q, Q->2Q, −Q->O, −2Q->−Q
    F = gq; C.cx(fO, F); C.cx(fm2Q, F)                                   # F = любой флаг
    and_xor(C, F, R['xdq'], X); and_xor(C, F, R['ydq'], Y)
    C.cx(fO, F); C.cx(fm2Q, F)                                           # F = gq снова
    h = C.alloc(); C.cx(fO, h); C.cx(fm2Q, h)                            # h = fO ∨ fm2Q: x образа = x_Q
    and_xor(C, h, R['xq'], X)
    and_xor(C, fO, R['yq'], Y); and_xor(C, fm2Q, R['ymq'], Y)
    and_xor(C, fQ, R['x2q'], X); and_xor(C, fQ, R['y2q'], Y)
    C.cx(fO, h); C.cx(fm2Q, h); C.release(h)
    C.cx(fQ, gq); C.cx(fmQ, gq); C.release(gq)
    # --- снятие флагов по результату (fmQ = [R==O] — охрана)
    g, und = eq_reg_chain(C, X, R['xq'], guard=fmQ)                      # R ∈ {Q, −Q}
    C.cx(g, fm2Q); C.cx(fO, fm2Q); C.release(fm2Q)
    parity_lit(C, g, Y[0], R['yq'][0], erase_target=fO)
    und()
    g, und = eq_reg_chain(C, X, R['x2q'], guard=fmQ)                     # R == 2Q
    parity_lit(C, g, Y[0], R['y2q'][0], erase_target=fQ)
    und()
    t, und = eq_chain(C, Y, 0, upto=n - 1)
    lit_erase(C, fmQ, t, Y[n - 1], 0)
    und()
    # --- снятие поиска
    lookup(C, Kw, vals, R)
    for f in reversed(FIELDS):
        for w in reversed(R[f]): C.release(w)
