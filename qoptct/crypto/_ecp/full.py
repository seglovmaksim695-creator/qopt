# SPDX-License-Identifier: Apache-2.0
"""k4/padd_full.py — ПОЛНОЕ (без исключений) смешанное сложение точек P <- P + Q, Q классическая (Q != O),
P — любая точка кривой или O (кодируется (0,0); для кривых нечётного порядка точек с y=0 нет).

Схема-обёртка над «общим» сложением (a9 point_add_mixed или улучшенные варианты k4), которое корректно только при
P not in E = {O, Q, -Q, -2Q} (знаменатели первого/второго деления не 0). Так как Q классическая, ответы на E — константы:
    O -> Q,   Q -> 2Q,   -Q -> O,   -2Q -> -Q.
1) флаги f_O=[y==0], f_Q, f_mQ (по [x==x_Q] и чётности y), f_m2Q (по [x==x_2Q] и чётности y)  — ~3n И;
2) подстановка безопасной точки D: P ^= f_j·(E_j xor D)  (только CX);
3) общее сложение: D -> D+Q;
4) исправление: R ^= f_j·((D+Q) xor (E_j+Q))  (только CX);
5) снятие флагов по результату: f_O <-> [R==Q], f_m2Q <-> [R==-Q], f_Q <-> [R==2Q], f_mQ <-> [R==O]  (~3n И).
Итого накладные ~6n И (n=256: ~1.5 тыс. против ~5 млн общей части). Требование: порядок Q не 2,3,4 (у P-256 и secp256k1
порядок группы — большое простое)."""
from .engine import *
from .ec import ec_add, ec_neg


def _lits_on(C, wires, bits):
    """X на проводах, где нужный бит = 0 (литерал «провод == бит» становится «провод == 1»)."""
    for w, b in zip(wires, bits):
        if not b: C.x(w)


def eq_chain(C, wires, val, keep_last=True, upto=None):
    """цепочка И для [wires == val]: возвращает (последний провод, undo). И = m-1 (или upto-1)."""
    m = len(wires) if upto is None else upto
    bits = [(val >> i) & 1 for i in range(len(wires))]
    _lits_on(C, wires[:m], bits[:m])
    ch = [wires[0]]
    for i in range(1, m): ch.append(C.AND(ch[-1], wires[i]))
    def undo(erase_last=True):
        top = m - 1 if erase_last else m - 2
        for i in range(top, 0, -1): C.erase(ch[i], ch[i - 1], wires[i])
        _lits_on(C, wires[:m], bits[:m])
    return ch[-1], undo


def lit_and(C, g, w, bit):
    if not bit: C.x(w)
    t = C.AND(g, w)
    if not bit: C.x(w)
    return t


def lit_erase(C, t, g, w, bit):
    if not bit: C.x(w)
    C.erase(t, g, w)
    if not bit: C.x(w)


def xor_const(C, f, wires, val):
    for i, w in enumerate(wires):
        if (val >> i) & 1: C.cx(f, w)


def point_add_full(C, X, Y, Q, cv, enc, generic, D=None):
    """X,Y — регистры в кодировке enc (enc(v) — классическое представление числа v, напр. Монтгомери v·2^n mod p).
    generic(C, X, Y, Q) — общее сложение (корректно вне E). Q=(xq,yq) в обычных координатах."""
    p, a = cv['p'], cv['a']
    n = len(X)
    Q2 = ec_add(Q, Q, p, a); mQ = ec_neg(Q, p); m2Q = ec_neg(Q2, p)
    assert Q2 is not None and m2Q not in (Q, mQ, None), "порядок Q должен быть > 4"
    E = {'O': None, 'Q': Q, 'mQ': mQ, 'm2Q': m2Q}
    img = {'O': Q, 'Q': Q2, 'mQ': None, 'm2Q': mQ}
    if D is None:
        D = ec_add(Q2, Q, p, a)                     # 3Q
    assert D is not None and D not in E.values() and ec_add(D, Q, p, a) not in (Q, mQ, None)
    def code(P):
        return (0, 0) if P is None else (enc(P[0]), enc(P[1]))
    def cat(P):
        x, y = code(P); return x | (y << n)
    XY = list(X) + list(Y)
    # --- 1) флаги (fO первым: если enc(x_Q)==0, то O=(0,0) совпадает с Q по x — тогда g <- g & ~fO, +1 И)
    fO, undO = eq_chain(C, Y, 0)
    undO(erase_last=False)
    def xcheck(xc):
        g, und = eq_chain(C, X, xc)
        if xc != 0: return g, und
        g2 = lit_and(C, g, fO, 0)
        def und2():
            lit_erase(C, g2, g, fO, 0); und()
        return g2, und2
    g, und = xcheck(enc(Q[0]))
    fQ = lit_and(C, g, Y[0], enc(Q[1]) & 1)
    fmQ = C.alloc(); C.cx(g, fmQ); C.cx(fQ, fmQ)
    und()
    g, und = xcheck(enc(Q2[0]))
    fm2Q = lit_and(C, g, Y[0], enc(m2Q[1]) & 1)
    und()
    F = {'O': fO, 'Q': fQ, 'mQ': fmQ, 'm2Q': fm2Q}
    # --- 2) подстановка D
    for k in F: xor_const(C, F[k], XY, cat(E[k]) ^ cat(D))
    if hasattr(C, 'mark'): C.mark('w_pre', fO=fO, fQ=fQ, fmQ=fmQ, fm2Q=fm2Q)
    # --- 3) общее сложение
    generic(C, X, Y, Q)
    # --- 4) исправление образов
    DQ = ec_add(D, Q, p, a)
    for k in F: xor_const(C, F[k], XY, cat(DQ) ^ cat(img[k]))
    # --- 5) снятие флагов по результату; f_mQ = [R==O] служит охраной (g <- g & ~f_mQ), если enc(x)==0
    def xcheck2(xc):
        g, und = eq_chain(C, X, xc)
        if xc != 0: return g, und
        g2 = lit_and(C, g, fmQ, 0)
        def und2():
            lit_erase(C, g2, g, fmQ, 0); und()
        return g2, und2
    g, und = xcheck2(enc(Q[0]))                    # R in {Q, -Q}: f_O (R==Q), f_m2Q (R==-Q)
    C.cx(g, fm2Q); C.cx(fO, fm2Q); C.release(fm2Q)
    lit_erase(C, fO, g, Y[0], enc(Q[1]) & 1)
    und()
    g, und = xcheck2(enc(Q2[0]))                   # R == 2Q <=> f_Q  (R == -2Q отличается чётностью y)
    lit_erase(C, fQ, g, Y[0], enc(Q2[1]) & 1)
    und()
    t, und = eq_chain(C, Y, 0, upto=n - 1)
    lit_erase(C, fmQ, t, Y[n - 1], 0)              # Y == 0  <=> R == O  <=> f_mQ
    und()
    if hasattr(C, 'mark'): C.mark('w_post')
    return dict(D=D)
