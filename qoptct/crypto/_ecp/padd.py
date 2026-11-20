# SPDX-License-Identifier: Apache-2.0
"""k4/padd_k4.py — общее смешанное сложение точек в ОБЫЧНОМ (не Монтгомери) представлении для специальных p
(P-256, secp256k1, псевдо-Мерсенн в тестах). Отличия от a9:
  * умножение: Карацуба с мусором (L уровней) + редукция Солинаса O(n) (вместо Монтгомери 3n^2);  λ^2 — возведение в квадрат;
  * деление: x' = den·2^(2n) mod p (умножение на разреженную константу, O(n)), Калиски k4 (halfB, Солинас вместо деления
    столбиком) -> K = −den^{-1};  λ' = num·K;
  * всё в схеме «вычислить-применить-обратить» (мусор снимается обращением).
Имена вариантов: 'k4L<L>' (например k4L3)."""
from .engine import *
from .mulk import kara_mul, kara_sq, add_gen, sub_gen
from .redk import red_special, naf
from .redk2 import reduce_k2
RED = {'mode': 'k2'}


_P256 = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff


def reduce_any(C, V, p, B):
    """V (value < B) -> V mod p: the fold reduction of qoptct.crypto.fold (mode 'k2') when it applies (V fits in 2n
    bits and c = 2^n - p satisfies c^2 + 2c + 2 < 2^n, or p = P-256), else the Solinas search red_special."""
    n = p.bit_length()
    c = (1 << n) - p
    fold_ok = (B <= (1 << (2 * n)) and c * c + 2 * c + 2 < (1 << n)) or (p == _P256 and B <= (1 << 512))
    if RED['mode'] == 'k2' and fold_ok:
        V = list(V)[:2 * n]
        pad = [C.alloc() for _ in range(2 * n - len(V))]
        out, gar = reduce_k2(C, V + pad, p)
        return out, gar + pad
    return red_special(C, V, p, B)

from .kalk import kaliski_forward


def mod_add_const(C, b, k, p):
    k %= p
    if k == 0: return
    r = [C.alloc() for _ in b]
    for i, w in enumerate(r):
        if (k >> i) & 1: C.x(w)
    mod_add(C, r, b, p)
    for i, w in enumerate(r):
        if (k >> i) & 1: C.x(w)
    for w in reversed(r): C.release(w)


def const_mul_mod(C, x, c, p):
    """out = x·c mod p (вне места, мусор). Σ s·(x<<e) по NAF(c) + смещение-кратное p, затем red_special."""
    t = naf(c % p)
    n = len(x); xm = (1 << n) - 1
    neg = sum(xm << e for s, e in t if s < 0); pos = sum(xm << e for s, e in t if s > 0)
    K = -(-neg // p) * p
    B = K + pos + 1
    A = [C.alloc() for _ in range(B.bit_length())]
    for i, w in enumerate(A):
        if (K >> i) & 1: C.x(w)
    for s, e in sorted(t, key=lambda u: -u[0]):
        (add_gen if s > 0 else sub_gen)(C, [W(w) for w in x][:len(A) - e], A[e:])
    out, gar = reduce_any(C, A, p, B)
    return out, gar + [w for w in A if w not in set(out)]


MUL = {'base': 'school'}


def pmul(C, a, b, p, L):
    n = len(a)
    sq = a is b or list(a) == list(b)
    if MUL['base'] == 'sd':
        from . import mulsd
        P, gar = mulsd.kara(C, a, None if sq else b, L, 4, sq)
    else:
        P = [C.alloc() for _ in range(2 * n)]
        gar = kara_sq(C, a, P, L) if sq else kara_mul(C, a, b, P, L)
    out, g2 = reduce_any(C, P, p, (p - 1) ** 2 + 1)
    so = set(out)
    return out, [w for w in dict.fromkeys(gar + P + g2) if w not in so]


def make(p, L):
    n = p.bit_length()
    c2 = pow(2, 2 * n, p)

    def prod_apply(C, a, b, apply):
        def fn():
            out, gar = pmul(C, a, b, p, L); return out, list(out) + gar
        compute_apply_uncompute(C, fn, apply)

    def div_apply(C, num, den, apply):
        def fn():
            xs, g0 = const_mul_mod(C, den, c2, p)            # x' = den·2^(2n)
            A, B, K, S, g1 = kaliski_forward(C, xs, p, shift=True, trunc=True, halfB=True, red=RED['mode'])   # K = −den^{-1}
            P, g2 = pmul(C, num, K, p, L)
            live = list(dict.fromkeys(list(P) + g2 + list(K) + g1 + list(S) + list(A) + list(B) + list(xs) + g0))
            dset = set(den) | set(num)
            return P, [w for w in live if w not in dset]
        compute_apply_uncompute(C, fn, apply)

    def generic(C, X, Y, Q):
        xq, yq = Q
        mod_add_const(C, X, -xq, p); mod_add_const(C, Y, -yq, p)
        Lr = [C.alloc() for _ in range(n)]
        def xor_into(T):
            def f(P):
                for a_, b_ in zip(P, T): C.cx(a_, b_)
            return f
        div_apply(C, Y, X, xor_into(Lr))                               # L = −λ
        prod_apply(C, Lr, X, lambda P: mod_add(C, P, Y, p))            # Y = dy − λ dx = 0
        prod_apply(C, Lr, Lr, lambda P: mod_rsub(C, P, X, p))          # X = λ² − dx
        mod_add_const(C, X, -3 * xq, p)                                # X = x' − xq
        prod_apply(C, Lr, X, lambda P: mod_add(C, P, Y, p))            # Y = λ(xq − x') = y' + yq
        mod_neg(C, X, p)                                               # X = xq − x'
        div_apply(C, Y, X, xor_into(Lr))                               # L ^= L -> 0
        mod_neg(C, X, p); mod_add_const(C, X, xq, p); mod_add_const(C, Y, -yq, p)   # (x', y')
        for w in reversed(Lr): C.release(w)
    return generic


def get(name, n, p, cv):
    if name.startswith('k4sL'):
        MUL['base'] = 'sd'; L = int(name[4:])
    else:
        MUL['base'] = 'school'; L = int(name[3:]) if name.startswith('k4L') else 3
    return (lambda v: v % p), make(p, L)


def mod_neg(C, b, p):
    """b <- -b mod p (0 -> 0). И = 3n-4."""
    n = len(b)
    for w in b: C.x(w)
    ch = [b[0]]; recs = []
    for i in range(1, n):            # ch[i] = AND(~b0..~bi)
        ch.append(C.AND(ch[-1], b[i]))
    z = ch[-1]                       # z = [b == 0]  (провода b инвертированы)
    for i in range(n - 2, 0, -1): C.erase(ch[i], ch[i - 1], b[i])
    for w in b: C.x(w)
    C.x(z)                           # z = [b != 0]
    for w in b: C.cx(z, w)           # ~b при z
    add_k(C, (p + 1) % (1 << n), b, ctrl=z)   # p - b
    C.x(z)                           # z = [b == 0] = [новое b == 0]
    for w in b: C.x(w)
    ch = [b[0]]
    for i in range(1, n - 1): ch.append(C.AND(ch[-1], b[i]))
    C.erase(z, ch[-1], b[n - 1])
    for i in range(n - 2, 0, -1): C.erase(ch[i], ch[i - 1], b[i])
    for w in b: C.x(w)


def mod_rsub(C, a, b, p):
    """b <- a - b mod p."""
    mod_neg(C, b, p); mod_add(C, a, b, p)
