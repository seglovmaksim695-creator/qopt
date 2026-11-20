# SPDX-License-Identifier: Apache-2.0
"""Adapter: the library reduction qoptct.crypto.fold (agent k2: p = 2^n - c by one fold, NIST P-256 by the Solinas
formula) on top of the engine of this package. The fold code is NOT copied: only the primitives (new / release / AND /
unAND / x / cx) are redirected to the engine circuit, so the reduction ANDs are counted and checked together with
the rest (2 744 ANDs for P-256, 2 175 for secp256k1)."""
from ..andcirc import Circ as _AndCirc
from .. import fold as _fold


class K2Adapter(_AndCirc):
    def __init__(self, C):
        super().__init__(); self.C = C

    def new(self, k=1):
        return [self.C.alloc() for _ in range(k)]

    def release(self, ws):
        for w in ws: self.C.release(w)

    def x(self, a): self.C.x(a)

    def cx(self, a, b): self.C.cx(a, b)

    def AND(self, a, b, t=None):
        self.nand += 1
        return self.C.AND(a, b, t)

    def unAND(self, a, b, t): self.C.erase(t, a, b, release=False)

    def z(self, a): raise NotImplementedError

    def cz(self, a, b): raise NotImplementedError

    def MZ(self, t, corr): raise NotImplementedError


def reduce_k2(C, P, p):
    """P: 2n wires (value < 2^(2n)) -> (out[n], garbage). p = P-256 or 2^n - c with c^2 + 2c + 2 < 2^n."""
    n = p.bit_length()
    assert len(P) == 2 * n, (len(P), n)
    A = K2Adapter(C)
    n, c, method, step2 = _fold._resolve(p, n, 'auto', 'auto')
    if method == 'solinas':
        return _fold._reduce_p256(A, list(P))
    return _fold._reduce_pm(A, list(P), n, c, step2)
