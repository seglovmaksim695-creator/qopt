# SPDX-License-Identifier: Apache-2.0
"""Specifications for qoptct.dyadic, each with an exact path sum (applied as an adjoint before the circuit).

QFTSpec(wires, swaps=True, inverse=False)
    F |x> = 2^(-m/2) sum_y exp(2 pi i X Y / 2^m) |y>,  X = sum_j 2^j x_j with x_j on wires[j] (little endian, as in
    Qiskit).  swaps=False: the "QFT with reversal" (output bit y_k on wires[m-1-k], Qiskit QFT(do_swaps=False)).
    inverse=True: the specification is F^dagger.  This is the QFT specification of Amy (arXiv:1805.06908v2,
    Example 2.2); no lemma is involved.

ArithSpec(target, terms)
    |x, s> -> |x, (s + G(x)) mod 2^m> on the m target wires (little endian), G(x) = sum_i w_i prod_{q in Q_i} x_q
    (integer weights; the wires of the terms must not be target wires).  Addition: G = sum_j 2^j a_j; a classical
    constant c: G = c; multiplication: G = sum_{i,j} 2^(i+j) a_i b_j; a control wire multiplies every term.
    Its path sum is the exact Fourier form
        A |x, s> = 2^(-m) sum_{Y, Z} exp(2 pi i (S + G(x) - Z) Y / 2^m) |x, Z>,
    and A equals the classical map by the geometric-sum lemma: sum_{Y=0}^{2^m-1} exp(2 pi i L Y / 2^m) = 2^m if
    2^m | L and 0 otherwise (L = S + G - Z is an integer; exactly one Z in [0, 2^m) has Z = S + G mod 2^m).
    The adjoint is applied first:  A^dagger |x, z> = 2^(-m) sum_{Y, S'} exp(2 pi i (Z - G(x) - S') Y / 2^m) |x, S'>.
    The frequency variables Y are frozen while the circuit is applied (strategy 'freeze'), so the rules first cancel
    the circuit's QFT against the specification and do not expand the carries of the classical map.

ClassicalSpec(n, func, zero_wires=())  /  PermutationSpec(n, func)      (small n: at most PERM_MAX_QUBITS free inputs)
    |b> -> |func(b)> for every basis state b whose zero_wires are 0 (b little endian over qubits 0..n-1); func must be
    injective there (checked).  It is checked without any lemma: the circuit's own path sum (zero wires start in 0)
    must reduce to a form without path variables whose wires are the ANF of func over the free inputs and whose phase
    is constant.
"""
from .core import ZERO, ONE, pmul

__all__ = ['QFTSpec', 'ArithSpec', 'PermutationSpec', 'ClassicalSpec', 'adder_spec', 'const_adder_spec',
           'multiplier_spec',
           'PERM_MAX_QUBITS', 'anf_over_vars']

PERM_MAX_QUBITS = 20


def _check_wires(ws, n, what):
    ws = [int(w) for w in ws]
    if len(set(ws)) != len(ws):
        raise ValueError('%s: repeated wire' % what)
    if any(w < 0 or w >= n for w in ws):
        raise ValueError('%s: wire out of range (n=%d)' % (what, n))
    return ws


def _qft_symbolic(ps, in_pos, out_pos, sign):
    """Apply 2^(-m/2) sum_z exp(sign 2 pi i X Z / 2^m) |z> where bit j of X is the wire in_pos[j]; the new path
    variable z_k goes to the wire out_pos[k].  Coefficient of [x_j] z_k: sign * 2^(K-m+j+k) (0 when j+k >= m)."""
    m = len(in_pos)
    if ps.K < m:
        raise ValueError('K=%d too small for a %d-qubit QFT' % (ps.K, m))
    W = [ps.w[q] for q in in_pos]
    zs = [ps.var('y', 'qft') for _ in range(m)]
    for j in range(m):
        for k in range(m - j):
            c = (sign << (ps.K - m + j + k)) & ps.M
            ps.add_lift(W[j], c, 1 << zs[k])
    for k in range(m):
        ps.set_wire(out_pos[k], frozenset([1 << zs[k]]))
    ps.e2 -= m
    return zs


class QFTSpec:
    kind = 'qft'
    lemma = None

    def __init__(self, wires, swaps=True, inverse=False):
        self.wires = [int(w) for w in wires]
        self.swaps = bool(swaps)
        self.inverse = bool(inverse)
        if not self.wires:
            raise ValueError('QFTSpec needs at least one wire')

    def describe(self):
        return 'QFT%s on %d wires%s' % ('^dagger' if self.inverse else '', len(self.wires),
                                        '' if self.swaps else ' (with reversal)')

    def check(self, n):
        _check_wires(self.wires, n, 'QFTSpec')

    def needed_K(self):
        return len(self.wires)

    def out_positions(self):
        m = len(self.wires)
        return list(self.wires) if self.swaps else [self.wires[m - 1 - k] for k in range(m)]

    def apply_adjoint(self, ps):
        """Apply the adjoint of the specification to the current wires; returns the variables to freeze (none)."""
        outp = self.out_positions()
        if self.inverse:        # spec = F^dagger -> adjoint = F
            _qft_symbolic(ps, self.wires, outp, +1)
        else:                   # spec = F -> adjoint = F^dagger: reads y in the output convention
            _qft_symbolic(ps, outp, self.wires, -1)
        return set()


class ArithSpec:
    kind = 'arith'
    lemma = 'fourier-shift'

    def __init__(self, target, terms, name='arith'):
        self.target = [int(w) for w in target]
        self.terms = []
        for w, qs in terms:
            if int(w) != w:
                raise ValueError('ArithSpec weights must be integers')
            self.terms.append((int(w), tuple(sorted(set(int(q) for q in qs)))))
        self.name = name
        if not self.target:
            raise ValueError('ArithSpec needs at least one target wire')

    def describe(self):
        return '%s: target %d wires += G (%d terms) mod 2^%d' % (self.name, len(self.target), len(self.terms),
                                                                 len(self.target))

    def check(self, n):
        _check_wires(self.target, n, 'ArithSpec target')
        ts = set(self.target)
        for w, qs in self.terms:
            _check_wires(qs, n, 'ArithSpec term')
            if ts & set(qs):
                raise ValueError('ArithSpec: a term uses a target wire')

    def needed_K(self):
        return len(self.target)

    def G(self, bit):
        """Integer G at a basis state given as a function q -> bit (for the independent simulator)."""
        s = 0
        for w, qs in self.terms:
            if all(bit(q) for q in qs):
                s += w
        return s

    def apply_adjoint(self, ps):
        m = len(self.target)
        K, M = ps.K, ps.M
        if K < m:
            raise ValueError('K too small')
        Zw = [ps.w[q] for q in self.target]
        Y = [ps.var('y', 'freq') for _ in range(m)]
        S = [ps.var('y', 'out') for _ in range(m)]
        for k in range(m):
            yb = 1 << Y[k]
            for j in range(m - k):
                c = 1 << (K - m + j + k)
                ps.addc((1 << S[j]) | yb, -c)            # - S' Y
                ps.add_lift(Zw[j], c, yb)                # + Z Y
            for w, qs in self.terms:                     # - G Y
                c = (-w << (K - m + k)) & M
                if not c:
                    continue
                P = ONE
                for q in qs:
                    P = pmul(P, ps.w[q])
                ps.add_lift(P, c, yb)
        for j in range(m):
            ps.set_wire(self.target[j], frozenset([1 << S[j]]))
        ps.e2 -= 2 * m
        return set(Y)


def adder_spec(a, b, carry=None):
    """b (+ carry as the top bit) += a: Qiskit DraperQFTAdder ('fixed': carry None, 'half': carry = cout wire)."""
    target = list(b) + ([carry] if carry is not None else [])
    return ArithSpec(target, [(1 << j, (q,)) for j, q in enumerate(a)], name='add')


def const_adder_spec(c, b, controls=()):
    """b += c (classical constant), optionally controlled."""
    return ArithSpec(list(b), [(int(c), tuple(controls))] if c else [], name='addconst')


def multiplier_spec(a, b, out, controls=()):
    """out += a * b (mod 2^len(out)): Qiskit RGQFTMultiplier / HRSCumulativeMultiplier on |a>|b>|out>."""
    terms = [(1 << (i + j), (qa, qb) + tuple(controls)) for i, qa in enumerate(a) for j, qb in enumerate(b)]
    return ArithSpec(list(out), terms, name='mul')


class ClassicalSpec:
    """|b> -> |func(b)> for every basis state b whose zero_wires are 0 (little endian over qubits 0..n-1), for small
    n (at most PERM_MAX_QUBITS free inputs).  func must be injective on that domain (checked), so the specification
    is an isometry.  Checked without any lemma: the circuit's own path sum (zero_wires start in 0) must reduce to a
    form without path variables whose wires are the ANF of the output bits of func as functions of the free inputs
    (fast Moebius transform of the truth table) and whose phase is constant.  Example: an out-of-place multiplier
    |a, b, 0> -> |a, b, a*b> with zero_wires = the output register."""
    kind = 'classical'
    lemma = None

    def __init__(self, n, func, zero_wires=(), name='classical'):
        self.n = int(n)
        self.zero_wires = sorted(set(int(q) for q in zero_wires))
        if any(q < 0 or q >= self.n for q in self.zero_wires):
            raise ValueError('ClassicalSpec: zero wire out of range')
        self.free = [q for q in range(self.n) if q not in self.zero_wires]
        if len(self.free) > PERM_MAX_QUBITS:
            raise ValueError('ClassicalSpec is limited to %d free inputs' % PERM_MAX_QUBITS)
        self.name = name
        table = []
        for idx in range(1 << len(self.free)):
            b = 0
            for i, q in enumerate(self.free):
                if (idx >> i) & 1:
                    b |= 1 << q
            v = int(func(b))
            if v < 0 or v >> self.n:
                raise ValueError('ClassicalSpec: func(%d) = %d is out of range' % (b, v))
            table.append(v)
        if len(set(table)) != len(table):
            raise ValueError('ClassicalSpec: func is not injective on the inputs with zero_wires = 0')
        self.table = table

    def describe(self):
        return '%s: classical map of %d qubits (%d zero wires)' % (self.name, self.n, len(self.zero_wires))

    def check(self, n):
        if n != self.n:
            raise ValueError('ClassicalSpec has %d qubits, the circuit %d' % (self.n, n))

    def needed_K(self):
        return 3

    def index(self, b):
        """Basis state with zero_wires 0 -> index into the truth table."""
        return sum(((b >> q) & 1) << i for i, q in enumerate(self.free))

    def anf(self, q):
        """ANF of output bit q over the free inputs, as masks over QUBITS (fast Moebius transform)."""
        t = [(v >> q) & 1 for v in self.table]
        h, size = 1, len(t)
        while h < size:
            for i in range(0, size, 2 * h):
                for j in range(i, i + h):
                    t[j + h] ^= t[j]
            h *= 2
        out = []
        for m in range(size):
            if t[m]:
                mq = 0
                for i, qq in enumerate(self.free):
                    if (m >> i) & 1:
                        mq |= 1 << qq
                out.append(mq)
        return out


def PermutationSpec(n, func, name='perm'):
    """A bijection of all 2^n basis states (ClassicalSpec without zero wires)."""
    return ClassicalSpec(n, func, (), name=name)


def anf_over_vars(masks, xv):
    """Masks over qubits -> ANF polynomial over the path-sum variables xv[q]."""
    out = set()
    for m in masks:
        mo, q = 0, 0
        while m:
            if m & 1:
                mo |= 1 << xv[q]
            m >>= 1
            q += 1
        out.add(mo)
    return frozenset(out) if out else ZERO
