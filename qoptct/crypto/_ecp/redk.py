# SPDX-License-Identifier: Apache-2.0
"""k4/redk.py — редукция V mod p для специальных простых (2^m mod p — короткая знаковая сумма степеней двойки):
P-256 (2^256 ≡ 2^224 − 2^192 − 2^96 + 1), secp256k1 (2^256 ≡ 2^32 + 977), псевдо-Мерсенн 2^m − c для тестов.
Вне места, с мусором (для «вычислить-применить-обратить»). Раунды: V = H·2^m + L ≡ K + L + Σ s_k·(H << e_k), K — кратное p
(смещение, чтобы не уйти в минус); затем двоичный поиск частного: f=[A ≥ 2^j p], A −= f·2^j p."""
from .engine import *
from .mulk import add_gen, sub_gen


def naf(c):
    out = []; e = 0
    while c:
        if c & 1:
            d = 2 - (c & 3)          # +1 или -1
            out.append((d, e)); c -= d
        c >>= 1; e += 1
    return out


def const_ops(val, L):
    return [ONE if (val >> i) & 1 else ZERO for i in range(L)]


def cmp_ge_const(C, A, K):
    """f := [A >= K] (A — провода, K — константа) через перенос A + (2^L − K). И <= L."""
    L = len(A)
    if K >= (1 << L):
        return None
    f = C.alloc()
    X = [W(w) for w in A]; Y = const_ops((1 << L) - K, L)
    compute_cout(C, X, Y, f)
    return f


def red_special(C, V, p, Vbound, ws=None):
    """V — провода (значение < Vbound). Возвращает (out[m], garbage)."""
    m = p.bit_length()
    c = (1 << m) % p
    terms = naf(c)
    gar = []
    cur, B = list(V), Vbound
    while B >= (1 << (m + 2)) and len(cur) > m:
        Lw, Hw = cur[:m], cur[m:]
        Hmax = (B - 1) >> m
        hb = Hmax.bit_length(); Hw = Hw[:hb]
        neg = sum((Hmax << e) for s, e in terms if s < 0)
        pos = sum((Hmax << e) for s, e in terms if s > 0)
        K = -(-neg // p) * p
        Bn = K + (1 << m) + pos
        A = [C.alloc() for _ in range(Bn.bit_length())]
        for i, w in enumerate(A):
            if (K >> i) & 1: C.x(w)
        add_gen(C, [W(w) for w in Lw], A)
        for s, e in sorted(terms, key=lambda t: -t[0]):     # сначала положительные
            if e >= len(A): continue
            Y = [W(w) for w in Hw][:len(A) - e]
            (add_gen if s > 0 else sub_gen)(C, Y, A[e:])
        gar += cur
        cur, B = A, Bn
    # финальная редукция: A < B, двоичный поиск частного
    A = cur
    q = (B - 1) // p
    j = q.bit_length() - 1
    while j >= 0:
        K = p << j
        f = cmp_ge_const(C, A, K)
        if f is not None:
            add_k(C, (1 << len(A)) - K, A, ctrl=f)
            gar.append(f)
        j -= 1
    out = A[:m]
    gar += A[m:]
    return out, gar
