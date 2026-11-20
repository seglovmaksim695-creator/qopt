# SPDX-License-Identifier: Apache-2.0
"""padd2: the generic mixed point addition assembled from the PROVED library blocks, with segment marks for
qoptct.verify.compose.ecadd (the default of qoptct.crypto.ecpoint.ec_point_add; legacy=True keeps padd.py).

Same structure as padd.py (lab k4: two divisions and three products in compute / apply / uncompute form, ordinary
representation), but every block is one that a prover of qoptct.verify handles at n = 256:
  * inversion: the core of qoptct.crypto.kaliski (12n^2 - 8n - 3 construction, items P1-P5; reduction='none': the raw
    2n-bit register r_ext is returned) -- proved segment by segment by compose.kal (ROBDDs + the Z3 lemmas M1..M6);
  * products / squares: qoptct.crypto.mulsq._kara (conditionally inverted additions, Karatsuba) -- proved by
    qoptct.verify.algebra (out == x*y);
  * reductions (of a product, of r_ext, of the constant multiple): qoptct.crypto.fold when it applies (p = 2^n - c
    with c^2 + 2c + 2 < 2^n, or P-256) -- proved by qoptct.verify.algebra; otherwise the signed-term folding +
    quotient search of qoptct.crypto.kaliski._red_special (any odd p) -- proved by compose (ROBDD segments);
  * modular additions, negations, constant additions, the multiplication by the sparse constant 2^(2n) mod p: the
    engine of this package -- ROBDD segments.
The marks (C.mark) are an UNTRUSTED certificate: positions and wire names only."""
from .engine import *
from .mulk import add_gen, sub_gen
from .redk import naf
from .redk2 import K2Adapter
from .padd import mod_neg
from .. import kaliski as _kal
from .. import mulsq as _mulsq
from .. import fold as _fold


def fold_params(p):
    """(n, c, method, step2) of qoptct.crypto.fold for p, or None when one fold does not apply."""
    try:
        return _fold._resolve(p, p.bit_length(), 'auto', 'auto')
    except ValueError:
        return None


def reduce2(C, V, p, B):
    """V (wires, value < B) -> (out[n] = V mod p, garbage). Mark 'fold' (library fold on 2n wires) or the marks
    sp_const / sp_addL / sp_term / sp_q of kaliski._red_special followed by 'sp_end'."""
    n = p.bit_length()
    fp = fold_params(p)
    if fp is not None and B <= (1 << (2 * n)) and len(V) <= 2 * n:
        V = list(V)
        pad = [C.alloc() for _ in range(2 * n - len(V))]
        A = K2Adapter(C)
        _, c, method, step2 = fp
        if method == 'solinas':
            out, gar = _fold._reduce_p256(A, V + pad)
        else:
            out, gar = _fold._reduce_pm(A, V + pad, n, c, step2)
        inp = set(V + pad)
        C.mark('fold', pad=list(pad), extra=[w for w in gar if w not in inp])
        return out, gar + [w for w in pad if w not in set(gar)]
    C.mark('sp_begin')
    out, gar = _kal._red_special(C, list(V), p, B)
    C.mark('sp_end')
    return out, gar


def mod_add_const(C, b, k, p):
    """b <- b + k mod p (b < p); mark 'addc'."""
    k %= p
    if k:
        r = [C.alloc() for _ in b]
        for i, w in enumerate(r):
            if (k >> i) & 1: C.x(w)
        mod_add(C, r, b, p)
        for i, w in enumerate(r):
            if (k >> i) & 1: C.x(w)
        for w in reversed(r): C.release(w)
    C.mark('addc')


def const_mul_mod(C, x, c, p):
    """out = x*c mod p (out of place, garbage): K + sum s*(x << e) over the NAF of c (K a multiple of p that keeps the
    value >= 0), then reduce2. Marks cm_const, cm_term."""
    t = naf(c % p)
    n = len(x); xm = (1 << n) - 1
    neg = sum(xm << e for s, e in t if s < 0); pos = sum(xm << e for s, e in t if s > 0)
    K = -(-neg // p) * p
    B = K + pos + 1
    A = [C.alloc() for _ in range(B.bit_length())]
    for i, w in enumerate(A):
        if (K >> i) & 1: C.x(w)
    C.mark('cm_const', A=list(A), K=K)
    for s, e in sorted(t, key=lambda u: -u[0]):
        (add_gen if s > 0 else sub_gen)(C, [W(w) for w in x][:len(A) - e], A[e:])
        C.mark('cm_term', s=s, e=e)
    out, gar = reduce2(C, A, p, B)
    so = set(out)
    return out, [w for w in dict.fromkeys(gar + A) if w not in so]


def pmul(C, a, b, p, L):
    """(a*b mod p on n wires, garbage); a is b: square. Marks 'mul' / 'sq', then those of reduce2."""
    sq = a is b or list(a) == list(b)
    gar = []
    P = _mulsq._kara(C, list(a), None if sq else list(b), L, 4, sq, gar)
    sp = set(P)
    gar = [w for w in dict.fromkeys(gar) if w not in sp]
    C.mark('sq' if sq else 'mul', out=list(P), gar=list(gar))
    out, g2 = reduce2(C, P, p, (p - 1) ** 2 + 1)
    so = set(out)
    return out, [w for w in dict.fromkeys(gar + P + g2) if w not in so]


def _cau(C, fn, apply):
    """compute / apply / uncompute with the marks 'apply' (after the apply part) and 'mirror' (after the inverse)."""
    s = len(C.ops); res, live = fn(); e = len(C.ops)
    apply(res)
    C.mark('apply')
    C.ops.append_inverse(s, e)
    C.mark('mirror')
    for w in reversed(live): C.release(w)


def make(p, L):
    n = p.bit_length()
    c2 = pow(2, 2 * n, p)

    def prod_apply(C, a, b, apply):
        def fn():
            out, gar = pmul(C, a, b, p, L); return out, list(out) + gar
        _cau(C, fn, apply)

    def div_apply(C, num, den, apply):
        def fn():
            xs, g0 = const_mul_mod(C, den, c2, p)                       # x' = den * 2^(2n) mod p
            u, v, Rr, Ss, g1 = _kal._kaliski_forward(C, xs, p, 'none')  # r_ext = -x'^{-1} 2^(2n) = -den^{-1} (mod p)
            K, g2 = reduce2(C, Rr, p, 1 << (2 * n))
            P, g3 = pmul(C, num, K, p, L)
            live = list(dict.fromkeys(list(P) + g3 + list(K) + g2 + list(Rr) + list(Ss) + g1 + list(u) + list(v)
                                      + list(xs) + g0))
            dset = set(den) | set(num)
            return P, [w for w in live if w not in dset]
        _cau(C, fn, apply)

    def generic(C, X, Y, Q):
        xq, yq = Q
        mod_add_const(C, X, -xq, p); mod_add_const(C, Y, -yq, p)
        Lr = [C.alloc() for _ in range(n)]
        C.mark('lam', Lr=list(Lr))

        def xor_into(T):
            def f(P):
                for a_, b_ in zip(P, T): C.cx(a_, b_)
            return f

        def madd(T):
            def f(P):
                mod_add(C, P, T, p)
            return f

        def rsub(T):
            def f(P):
                mod_neg(C, T, p); C.mark('neg'); mod_add(C, P, T, p)
            return f
        div_apply(C, Y, X, xor_into(Lr))                               # L = -lambda
        prod_apply(C, Lr, X, madd(Y))                                  # Y = dy - lambda dx = 0
        prod_apply(C, Lr, Lr, rsub(X))                                 # X = lambda^2 - dx
        mod_add_const(C, X, -3 * xq, p)                                # X = x' - xq
        prod_apply(C, Lr, X, madd(Y))                                  # Y = lambda (xq - x') = y' + yq
        mod_neg(C, X, p); C.mark('neg')                                # X = xq - x'
        div_apply(C, Y, X, xor_into(Lr))                               # L ^= L -> 0
        mod_neg(C, X, p); C.mark('neg')
        mod_add_const(C, X, xq, p); mod_add_const(C, Y, -yq, p)        # (x', y')
        for w in reversed(Lr): C.release(w)
    return generic
