# SPDX-License-Identifier: Apache-2.0
"""compose.bv -- bit-vector arithmetic written once, over two boolean algebras.

    IntAlg(S)   a boolean is a python int used as a mask over S samples (bit-sliced execution on many inputs)
    BddAlg(B)   a boolean is an edge of the ROBDD manager B (qoptct.verify.bdd): symbolic, all inputs at once

Bit vectors are python lists, least significant bit first. The functions below are the DEFINITIONS of the word
operations used in the contracts (ripple formulas); tests/test_compose.py compares them with python integers.
"""
from ..bdd import ONE, ZERO

__all__ = ['IntAlg', 'BddAlg', 'const', 'carry', 'add', 'sub', 'le', 'lt_const', 'mux', 'allzero']


class IntAlg:
    def __init__(self, S):
        self.S = S
        self.T = (1 << S) - 1
        self.F = 0

    def and_(self, a, b):
        return a & b

    def xor(self, a, b):
        return a ^ b

    def or_(self, a, b):
        return a | b

    def not_(self, a):
        return a ^ self.T

    def ite(self, c, a, b):
        return (c & a) | ((c ^ self.T) & b)


class BddAlg:
    def __init__(self, B):
        self.B = B
        self.T = ONE
        self.F = ZERO
        self.and_ = B.AND
        self.xor = B.XOR

    def or_(self, a, b):
        return self.B.AND(a ^ 1, b ^ 1) ^ 1

    def not_(self, a):
        return a ^ 1

    def ite(self, c, a, b):
        if c == ONE:
            return a
        if c == ZERO:
            return b
        if a == b:
            return a
        A = self.B.AND
        return A(A(c, a) ^ 1, A(c ^ 1, b) ^ 1) ^ 1


def const(al, k, w):
    """the constant k as w bits (k must fit)."""
    assert 0 <= k < (1 << w), (k, w)
    return [al.T if (k >> i) & 1 else al.F for i in range(w)]


def _maj(al, a, b, c):
    return al.or_(al.and_(a, b), al.and_(c, al.xor(a, b)))


def carry(al, x, y, cin):
    """carry out of x + y + cin (equal lengths)."""
    assert len(x) == len(y)
    c = cin
    for a, b in zip(x, y):
        c = _maj(al, a, b, c)
    return c


def add(al, x, y, cin=None):
    """(x + y + cin) mod 2^w as w bits, and the carry out (equal lengths w)."""
    assert len(x) == len(y)
    c = al.F if cin is None else cin
    out = []
    for a, b in zip(x, y):
        t = al.xor(a, b)
        out.append(al.xor(t, c))
        c = al.or_(al.and_(a, b), al.and_(c, t))
    return out, c


def sub(al, x, y):
    """(x - y) mod 2^w as w bits, and the borrow [x < y]."""
    s, c = add(al, x, [al.not_(b) for b in y], al.T)
    return s, al.not_(c)


def le(al, x, y):
    """[x <= y] (equal lengths): the carry out of y + ~x + 1."""
    return carry(al, y, [al.not_(a) for a in x], al.T)


def lt_const(al, x, k):
    """[value(x) < k] for a constant k >= 0."""
    if k >= (1 << len(x)):
        return al.T
    return al.not_(carry(al, x, [al.not_(b) for b in const(al, k, len(x))], al.T))


def mux(al, c, x, y):
    """c ? x : y, bitwise."""
    return [al.ite(c, a, b) for a, b in zip(x, y)]


def allzero(al, x):
    r = al.T
    for a in x:
        r = al.and_(r, al.not_(a))
    return r
