"""Bit-parallel classical simulation of logical programs (x, cx, and, erase).

Every wire is one Python integer whose bit L is the value of the wire on lane L.  Two lane sets are used:

* ``full``: all 2^n inputs (n free bits)                       -> a proof by exhaustion for classical circuits;
* ``w2``:   the 1 + n + n(n-1)/2 inputs of Hamming weight <= 2 -> determines every function of GF(2) degree <= 2 exactly
  (Moebius inversion; the same fact exact2 uses).

``and`` needs a clean target on EVERY lane and ``erase`` needs ``target = p AND q`` on every lane; otherwise the AND gadget
of the OpenQASM text is not the classical AND, and ``PreconditionError`` is raised (callers: verdict UNKNOWN).
The gadget itself (AND10 / erasure) is checked once by ``qoptct.verify.exact2.gadget_check``.
"""
from __future__ import annotations

__all__ = ['PreconditionError', 'masks_full', 'masks_w2', 'w2_lanes', 'simulate', 'w2_monomials']


class PreconditionError(Exception):
    def __init__(self, index, kind, wire):
        super().__init__(f'op {index}: {kind} on wire {wire}')
        self.index, self.kind, self.wire = index, kind, wire


def masks_full(n):
    """Lane masks of the n free bits over 2^n lanes (bit k of the lane number is free bit k)."""
    lanes = 1 << n
    out = []
    for k in range(n):
        unit = ((1 << (1 << k)) - 1) << (1 << k)
        m, w = unit, 1 << (k + 1)
        while w < lanes:
            m |= m << w
            w *= 2
        out.append(m)
    return out, (1 << lanes) - 1, lanes


def w2_lanes(n):
    return 1 + n + n * (n - 1) // 2


def masks_w2(n):
    """Lane masks of the n free bits over the weight <= 2 lanes: lane 0 = zero, 1+u = {u}, then the pairs {u<v}."""
    lanes = w2_lanes(n)
    masks = [0] * n
    idx = 1 + n
    for u in range(n):
        masks[u] |= 1 << (1 + u)
    for u in range(n):
        for v in range(u + 1, n):
            masks[u] |= 1 << idx
            masks[v] |= 1 << idx
            idx += 1
    return masks, (1 << lanes) - 1, lanes


def simulate(ops, nq, init, full):
    """Run ``ops`` on lanes.  ``init``: wire -> lane integer for the free wires; every other wire starts at 0."""
    st = [0] * nq
    for w, v in init.items():
        st[w] = v
    for i, op in enumerate(ops):
        k = op[0]
        if k == 'cx':
            st[op[2]] ^= st[op[1]]
        elif k == 'x':
            st[op[1]] ^= full
        elif k == 'and':
            a, b, t = op[1:]
            if st[t]:
                raise PreconditionError(i, 'AND target is not clean', t)
            st[t] = st[a] & st[b]
        elif k == 'erase':
            t = op[1]
            if len(op[2]) != 1 or op[2][0][0] != 'cz':
                raise PreconditionError(i, 'unsupported erase correction', t)
            a, b = op[2][0][1:]
            if st[t] != (st[a] & st[b]):
                raise PreconditionError(i, 'erased wire differs from p AND q', t)
            st[t] = 0
        else:
            raise PreconditionError(i, f'operation {k} is not classical', -1)
    return st


_PAIRS = {}


def _pair_index(n):
    if n not in _PAIRS:
        import numpy as np
        us, vs = np.triu_indices(n, 1)           # row-major order = the pair order of masks_w2
        _PAIRS[n] = (us.astype(np.int64), vs.astype(np.int64))
    return _PAIRS[n]


def w2_monomials(value, n):
    """Algebraic normal form (degree <= 2) of a wire given by its weight <= 2 lane integer:
    (constant, [u with a linear term], [(u, v) with a quadratic term])."""
    import numpy as np
    lanes = w2_lanes(n)
    bits = np.unpackbits(np.frombuffer(value.to_bytes((lanes + 7) // 8, 'little'), dtype=np.uint8), bitorder='little')[:lanes]
    c0 = int(bits[0])
    lin = bits[1:1 + n]
    linear = [int(u) for u in np.nonzero(lin ^ c0)[0]]
    us, vs = _pair_index(n)
    q = bits[1 + n:] ^ lin[us] ^ lin[vs] ^ c0
    idx = np.nonzero(q)[0]
    return c0, linear, [(int(us[i]), int(vs[i])) for i in idx]
