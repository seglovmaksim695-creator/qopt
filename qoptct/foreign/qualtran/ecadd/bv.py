"""Bit vectors on ROBDD edges of qoptct (LSB first). ONE=0, ZERO=1, NOT = ^1.
This is the language of the SPECIFICATIONS: the scheme takes no part here."""
from qoptct.verify.bdd import ONE, ZERO

# Simplification under an assumption: during the proof the executor makes SIMP(f) a constant if f is constant
# under the current assumption A (else f). The values of the formulas do not change on A: the checks run under A.
SIMP = None


def _s(f):
    return SIMP(f) if SIMP is not None else f


def const(c, w):
    assert c >= 0 and c < (1 << w), (c, w)
    return [ONE if (c >> i) & 1 else ZERO for i in range(w)]


def xor(B, X, Y):
    assert len(X) == len(Y)
    return [B.XOR(x, y) for x, y in zip(X, Y)]


def mask(B, c, X):
    return [B.AND(c, x) for x in X]


def ite(B, c, X, Y):
    assert len(X) == len(Y)
    return [B.ITE(c, x, y) for x, y in zip(X, Y)]


def add(B, X, Y):
    """X+Y, width w+1 (exact sum)."""
    assert len(X) == len(Y)
    out, c = [], ZERO
    for x, y in zip(X, Y):
        out.append(B.XOR(B.XOR(x, y), c))
        c = B.MAJ(x, y, c)
    return out + [c]


def sub(B, X, Y):
    """(X-Y mod 2^w, borrow) ; borrow = [X<Y]."""
    assert len(X) == len(Y)
    out, c = [], ONE
    for x, y in zip(X, Y):
        ny = y ^ 1
        out.append(B.XOR(B.XOR(x, ny), c))
        c = B.MAJ(x, ny, c)
    return out, c ^ 1


def lt(B, X, Y):
    return sub(B, X, Y)[1]


def lt_const(B, X, c):
    return lt(B, X, const(c, len(X)))


def eq(B, X, Y):
    assert len(X) == len(Y)
    e = ONE
    for x, y in zip(X, Y):
        e = B.AND(e, B.XOR(x, y) ^ 1)
    return e


def is_zero(B, X):
    e = ONE
    for x in X:
        e = B.AND(e, x ^ 1)
    return e


def all_of(B, *fs):
    e = ONE
    for f in fs:
        e = B.AND(e, f)
    return e


def modadd(B, X, Y, p):
    """(X+Y) mod p for X,Y<p."""
    w = len(X)
    S = add(B, X, Y)
    D, br = sub(B, S, const(p, w + 1))
    return ite(B, _s(br), S[:w], D[:w])


def modsub(B, X, Y, p):
    """(X-Y) mod p for X,Y<p."""
    w = len(X)
    D, br = sub(B, X, Y)
    E = add(B, D, const(p, w))[:w]
    return ite(B, _s(br), E, D)


def modneg(B, X, p):
    w = len(X)
    N, _ = sub(B, const(p, w), X)
    return ite(B, _s(is_zero(B, X)), const(0, w), N)


def moddbl(B, X, p):
    return modadd(B, X, X, p)
