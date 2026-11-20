# SPDX-License-Identifier: Apache-2.0
"""a9/kaliski.py — обратимая инверсия по Калиски (2n шагов) с симметризованным шагом и СЛИЯНИЕМ условных обменов
соседних шагов (обмен по d = w_{k-1} xor w_k), без проверки v==0 и без условного деления пополам.
Вход: U=p (константа, загружается X-ами), V=x, R=0, S=1. Выход (стандартная система): U=1, V=0, R = -x^-1 2^(2n) mod p, S=p,
плюс мусор (w_k, cmp_k, e_k на шаг)."""
from .engine import *
from .mulk import add_gen, sub_gen


def cswap(C, d, A, B):
    for a, b in zip(A, B):
        C.cx(a, b); t = C.AND(d, b); C.cx(t, a); C.erase(t, d, b); C.cx(a, b)


def kaliski_forward(C, V, p, keep_cmp=False, final_relabel=True, shift=False, merge=True, trunc=False, halfB=False, red='div'):
    n = len(V)
    U = [C.alloc() for _ in range(n)]; R = [C.alloc() for _ in range(n)]; S = [C.alloc() for _ in range(n)]
    for i in range(n):
        if (p >> i) & 1: C.x(U[i])
    C.x(S[0])
    A, B, Rr, Ss = U, list(V), R, S
    wprev = None; garbage = []
    for k in range(2 * n):
        # ширины (усечение по доказуемым границам): u·v < 2^(2n-k) => u,v < 2^(2n-k);  max(r,s) <= 2^k
        wa = min(n, max(1, 2 * n - 1 - k)) if trunc else n   # bitlen(u)+bitlen(v) <= 2n-k, v>=1 (иначе u=1)
        wr = min(n, k + 1) if trunc else n
        wr2 = min(n, k + 2) if trunc else n
        # (a) w_k по стандартной системе; регистры в системе w_{k-1}
        U0 = C.alloc(); V0 = C.alloc(); C.cx(A[0], U0); C.cx(B[0], V0)
        if wprev is not None:
            C.cx(A[0], B[0]); t = C.AND(wprev, B[0]); C.cx(A[0], B[0])     # t = w·(A0^B0)
            C.cx(t, U0); C.cx(t, V0)                                       # U0, V0 — младшие биты u, v
            for a, b in zip(A[:wa], B[:wa]): C.cx(wprev, a); C.cx(wprev, b)
        for a in A[:wa]: C.x(a)
        X_, Y_ = [W(a) for a in A[:wa]], [W(b) for b in B[:wa]]
        ch, RR = carry_chain(C, X_, Y_, ONE)
        rf = plan(X_[wa - 1], Y_[wa - 1], ch[wa - 1])
        if keep_cmp:
            cmp = C.alloc(); compute(C, rf, t=cmp)
        else:
            cmp = compute(C, rf, t=C.alloc())[1]                           # [u <= v], сотрём ниже
        C.x(cmp); g = C.AND(V0, cmp); C.x(cmp)                             # g = V0 & ~cmp
        C.x(g); w = C.AND(U0, g); C.x(g)                                   # w = U0 & ~g
        C.x(cmp); C.erase(g, V0, cmp); C.x(cmp)
        if not keep_cmp: uncompute(C, rf, t=cmp)
        chain_undo(C, RR)
        for a in A[:wa]: C.x(a)
        if wprev is not None:
            for a, b in zip(A[:wa], B[:wa]): C.cx(wprev, a); C.cx(wprev, b)
            C.cx(t, U0); C.cx(t, V0)
            C.cx(A[0], B[0]); C.erase(t, wprev, B[0]); C.cx(A[0], B[0])
        C.cx(A[0], U0); C.cx(B[0], V0); C.release(V0); C.release(U0)
        # (b) обмен по d = w_{k-1} ^ w_k
        if wprev is None:
            d = w
        else:
            d = C.alloc(); C.cx(wprev, d); C.cx(w, d)
        cswap(C, d, A[:wa], B[:wa]); cswap(C, d, Rr[:wr], Ss[:wr])
        if wprev is not None:
            C.cx(wprev, d); C.cx(w, d); C.release(d)
        # (c) e = u'_0
        e = C.alloc(); C.cx(A[0], e)
        # (d) A -= e·B  (бит 0: e·B0 = e)
        # halfB (k4): при e=1 после обмена u>=v, значит bitlen(v) <= (2n-k)//2 — старшие частичные произведения = 0
        wb = min(wa, max(1, (2 * n - k) // 2)) if halfB else wa
        m = [e] + [C.AND(e, B[i]) for i in range(1, wb)]
        sub_gen(C, [W(t) for t in m], A[:wa])
        for i in reversed(range(1, wb)): C.erase(m[i], e, B[i])
        # (e) A /= 2: переименование
        A = A[1:] + A[:1]
        # (f) Rr += e·Ss  (при e=1 r+s <= p; при усечении r,s <= 2^k -> результат < 2^(k+2))
        m = [C.AND(e, s) for s in Ss[:wr]]
        if wr2 > wr: add(C, m, Rr[:wr], cout=Rr[wr])
        else: add(C, m, Rr[:wr])
        for i in reversed(range(wr)): C.erase(m[i], e, Ss[i])
        # (g) Ss = 2 Ss mod p   или (shift) просто сдвиг: регистр растёт до 2n бит
        if shift:
            Ss = [C.alloc()] + Ss
            if len(Ss) > 2 * n: C.release(Ss.pop())        # старший бит доказуемо 0 (r_ext < 2^(2n))
        else:
            Ss = mod_double(C, Ss, p)
        garbage += [w, e] + ([cmp] if keep_cmp else []); wprev = w
        if not merge:     # абляция: обратный обмен в конце каждого шага (как у Литински)
            cswap(C, w, A, B); cswap(C, w, Rr[:n], Ss[:n]); wprev = None
    if not merge:
        pass
    elif final_relabel:     # w_{2n} = 1 при всех x в [1,p): последний обмен — переименование
        C.x(wprev); C.release(wprev); garbage.remove(wprev)
        A, B, Rr, Ss = B, A, Ss, Rr
    else:
        cswap(C, wprev, A, B); cswap(C, wprev, Rr[:n], Ss[:n])
    if shift and red in ('sol', 'k2'):   # k4: редукция O(n) для специальных p (red_special или k2)
        if red == 'k2':
            from .padd import reduce_any
            out, g2 = reduce_any(C, Rr, p, 1 << (2 * n))
        else:
            from .redk import red_special
            out, g2 = red_special(C, Rr, p, 1 << (2 * n))
        garbage += g2
        return A, B, out, Ss, garbage
    if shift:             # редукция r_ext (< 2^(2n)) mod p делением столбиком
        R = Rr + [C.alloc() for _ in range(2 * n + 1 - len(Rr))]
        for j in range(n, -1, -1):
            win = R[j:j + n + 1]
            add_k(C, (1 << (n + 1)) - p, win)
            f = C.alloc(); C.cx(win[n], f); garbage.append(f)
            add_k(C, p, win[:n], ctrl=f)
            C.cx(f, win[n])
        for w_ in R[n:]: C.release(w_)
        Rr = R[:n]
    return A, B, Rr, Ss, garbage
