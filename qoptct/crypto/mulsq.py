# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Maksym Shchehlov.
"""qoptct.crypto.mulsq -- out-of-place integer multiplication and squaring on conditionally inverted additions,
with optional Karatsuba levels (garbage allowed).  Pure generators: no optimisation, no search.

    multiply(n, karatsuba_levels=0, base=4, clean=False) -> MulSqCircuit    (x, y, 0) -> (x, y, x*y)
    square(n, karatsuba_levels=0, base=4, clean=False)   -> MulSqCircuit    (x, 0)    -> (x, x^2)

Cost model (Gidney 2018, arXiv:1709.06648): a logical AND into a |0> wire costs one Toffoli = 4 T; its erasure by an
X-basis measurement and a Clifford correction (CZ, Z) costs 0 T.  "AND" below = Toffoli count.

Constructions
-------------
Multiplication (base, karatsuba_levels=0).  With signed digits s_i = 2 y_i - 1 = +-1,
    2 x y = sum_i s_i 2^i x + (2^n - 1) x.
Row i adds +x or -x (two's complement: x XOR mask, carry-in = mask, mask = NOT y_i made by CX) into a window of
n+2 bits of a growing accumulator; the sign is extended by one CX per row.  No partial products are formed: a row
is a single "conditionally inverted" ripple addition with about n ANDs.  Then + 2^n x - x and the zero low bit is
dropped.  AND = n^2 + 4n exactly (n >= 1);  n = 256: 66 560.  No garbage: all ancillas return to |0>.

Squaring (base).  With L_i = x mod 2^i,
    x^2 = sum_{i>=1} s_i 2^i L_i - sum_i x_i 4^i + 2^n x,        s_i = 2 x_i - 1,
row i is a conditionally inverted addition of L_i (i bits) controlled by x_i.  AND = (n^2 + 7n - 6)/2 exactly
(n >= 1);  n = 256: 33 661.  No garbage.

Karatsuba with garbage (karatsuba_levels = L >= 1).  x = x0 + 2^h x1 (h = m // 2), z0 = x0 y0, z2 = x1 y1 written
directly into the low/high halves of the output, T = (x0 + x1)(y0 + y1) - z0 - z2 (= x0 y1 + x1 y0) added into the
middle.  Recursion stops at depth L or when the width is <= base (default 4, as used for the published numbers).
The half sums x0 + x1, y0 + y1 and the middle register T of every node are NOT uncomputed: they are the garbage
(MulSqCircuit.garbage).  It is removed by running the inverse of the whole circuit after the product has been used
(compute - use - uncompute, e.g. inside a modular multiplication / point addition that is itself reversed later);
the inverse costs one AND per measurement-based erasure of the forward circuit.  clean=True builds this directly:
forward circuit, CX copy of the product into a fresh register, inverse (Bennett); then AND = AND + ERASE of the
forward circuit and the garbage is gone.

    n = 256, base = 4                 AND        AND / n^2   qubits (incl. inputs, output, garbage)   garbage
    multiply  L=0                     66 560     1.02        1 537                                     0
    multiply  L=2                     42 727     0.65        2 707                                     1 300
    multiply  L=3                     37 880     0.58        3 911                                     2 504
    square    L=0                     33 661     0.51        1 279                                     0
    square    L=2                     23 333     0.36        2 126                                     975
    square    L=3                     22 332     0.34        3 029                                     1 878
(counted by building the circuit, no simulation; tests/test_crypto_mulsq.py.)  The base circuits contain a few
ANDs that are identically 0 (2 for multiply, 4 for square, from the first row added into the zero accumulator);
they are kept so that the counts match the published lab numbers.

Honest comparison with published circuits
-----------------------------------------
* Haner et al., arXiv:2609.05625 (2026), Sec. VII.A: integer squaring with conditionally inverted adders in
  n(n+3)/2 - 1 AND (n = 256: 33 151).  The idea of conditionally inverted additions for squaring is theirs.  Our
  base square (n^2+7n-6)/2 is 2n - 2 AND WORSE (33 661).  With Karatsuba (L=3) 22 332 AND, i.e. -33 %, but only at
  the price of many more qubits (garbage registers) and of an inverse pass to clear them; compare at equal qubits
  before claiming an improvement.  The same paper notes about n^2 for multiplication with this technique; our base
  multiplier n^2 + 4n is not better than that; Karatsuba L=3 gives 0.58 n^2 (37 880), again for qubits.
* Litinski 2023, arXiv:2306.08585, Fig. 8: MODULAR Montgomery multiplication 2.25 n^2 + 9 n (n = 256: 149 760),
  any modulus, garbage-free.  The circuits here compute the plain integer product (2n bits, no reduction mod p);
  the reduction must be added before comparing like with like.
* Karatsuba on quantum circuits itself is known (e.g. Gidney 2019, arXiv:1904.07356, garbage-free with a larger
  constant); what is new here is only the combination with conditionally inverted rows and garbage left for the
  enclosing uncomputation.

Wire layout of a MulSqCircuit (indices of the OpenQASM qubit array q[]):
    multiply:  x = 0..n-1, y = n..2n-1, out = 2n..4n-1, then garbage, then ancillas
    square:    x = 0..n-1,              out = n..3n-1,  then garbage, then ancillas
(little endian; out starts at 0; ancillas start and end at |0>; garbage starts at |0> and is left dirty.)

Checking: MulSqCircuit.simulate() runs the circuit bit-sliced on many inputs at once and verifies the
precondition of every AND (target |0>) and every erasure (target equals a*b XOR lin), inputs unchanged and
ancillas back to |0>.  prove() calls the exact prover qoptct.verify.prove (clean model, all inputs, all
measurement outcomes) against reference_qc(), an independent schoolbook Toffoli circuit.
"""
import os
import tempfile

__all__ = ['multiply', 'square', 'MulSqCircuit', 'CheckError', 'reference_qc', 'prove']

ZERO, ONE = ('c', 0), ('c', 1)


def _W(w):
    return ('w', w)


class CheckError(Exception):
    """a violated AND / erase precondition or a wrong value found by MulSqCircuit.simulate()."""


# ----------------------------------------------------------------------------------------- circuit IR

class _Circ:
    """ops: ('x', q) ('cx', c, t) ('and', a, b, t) [t must be |0>, t := a*b]
    ('erase', t, a, b, lin) [t must equal a*b XOR lin; X-basis measurement, on outcome 1: CZ a b, Z lin, X t]."""

    def __init__(self):
        self.ops = []
        self.nw = 0
        self.free = []

    def reg(self, n):
        return [self._new() for _ in range(n)]

    def _new(self):
        w = self.nw
        self.nw += 1
        return w

    def alloc(self):
        return self.free.pop() if self.free else self._new()

    def release(self, w):
        self.free.append(w)

    def x(self, q):
        self.ops.append(('x', q))

    def cx(self, c, t):
        assert c != t
        self.ops.append(('cx', c, t))

    def AND(self, a, b, t=None):
        if t is None:
            t = self.alloc()
        assert len({a, b, t}) == 3
        self.ops.append(('and', a, b, t))
        return t

    def erase(self, t, a, b, lin=(), release=True):
        assert t not in (a, b) and t not in lin and a != b
        self.ops.append(('erase', t, a, b, tuple(lin)))
        if release:
            self.release(t)


def _inverse(ops):
    out = []
    for o in reversed(ops):
        k = o[0]
        if k in ('x', 'cx'):
            out.append(o)
        elif k == 'and':
            out.append(('erase', o[3], o[1], o[2], ()))
        else:
            t, a, b, lin = o[1:]
            out.append(('and', a, b, t))
            for w in lin:
                out.append(('cx', w, t))
    return out


# ------------------------------------------------------- carry bit MAJ(x, y, c) with constant propagation

class _Rec:
    def __init__(self, kind, *args):
        self.kind = kind
        self.args = args
        self.t = None
        self.fold = False


def _plan(x, y, c):
    ops = [x, y, c]
    for i in range(3):
        for j in range(i + 1, 3):
            if ops[i] == ops[j]:
                return _Rec('op', ops[i])                  # MAJ(u, u, v) = u
    cs = [o for o in ops if o[0] == 'c']
    if len(cs) >= 2:                                       # two different constants 0, 1 -> the third operand
        return _Rec('op', [o for o in ops if o[0] == 'w'][0] if len(cs) == 2 else ZERO)
    if len(cs) == 1:
        u, v = [o[1] for o in ops if o[0] == 'w']
        return _Rec('and' if cs[0][1] == 0 else 'or', u, v)
    return _Rec('gid', x[1], y[1], c[1])                   # Gidney carry, x is the target wire of the sum


def _compute(C, r, fold_ok=False):
    k = r.kind
    if k == 'op':
        return r.args[0]
    if k == 'and':
        u, v = r.args
        r.t = C.AND(u, v)
        return _W(r.t)
    if k == 'or':
        u, v = r.args
        C.x(u); C.x(v)
        r.t = C.AND(u, v)
        C.x(u); C.x(v); C.x(r.t)
        return _W(r.t)
    x, y, c = r.args
    C.cx(c, x); C.cx(c, y)
    r.t = C.AND(x, y)
    C.cx(c, y)
    if fold_ok:
        r.fold = True                                      # x keeps x ^ c: it becomes the sum bit for free
    else:
        C.cx(c, x)
    C.cx(c, r.t)
    return _W(r.t)


def _uncompute(C, r):
    k, t = r.kind, r.t
    if k == 'and':
        u, v = r.args
        C.erase(t, u, v)
        return
    if k == 'or':
        u, v = r.args
        C.x(t); C.x(u); C.x(v)
        C.erase(t, u, v)
        C.x(u); C.x(v)
        return
    x, y, c = r.args
    if not r.fold:
        C.cx(c, x)
    C.cx(c, y)
    C.erase(t, x, y, lin=(c,))
    C.cx(c, y)
    if not r.fold:
        C.cx(c, x)


def _add(C, Y, b, cin=ZERO):
    """b += Y + cin (mod 2^len(b)); Y: operands _W(w) / ZERO / ONE (missing high operands = ZERO); cin operand.
    Gidney ripple adder with constant propagation: AND = len(b) - 1 minus the carries that are free."""
    m = len(b)
    Y = list(Y) + [ZERO] * (m - len(Y))
    X = [_W(w) for w in b]
    c = [cin]
    R = [None]
    for i in range(m - 1):
        r = _plan(X[i], Y[i], c[i])
        c.append(_compute(C, r, fold_ok=True))
        R.append(r)
    for i in reversed(range(m)):
        if i + 1 < m and R[i + 1].kind != 'op':
            _uncompute(C, R[i + 1])
        if Y[i] == ONE:
            C.x(b[i])
        elif Y[i][0] == 'w':
            C.cx(Y[i][1], b[i])
        folded = (i + 1 < m) and R[i + 1].kind == 'gid' and R[i + 1].fold
        if not folded:
            if c[i] == ONE:
                C.x(b[i])
            elif c[i][0] == 'w':
                C.cx(c[i][1], b[i])


def _sub(C, Y, b):
    """b -= Y (mod 2^len(b)) = NOT(NOT b + Y)."""
    for w in b:
        C.x(w)
    _add(C, Y, b)
    for w in b:
        C.x(w)


def _signed_row(C, A, off, src, ctrl):
    """A[off : off+len(src)+2] += (ctrl ? +src : -src): conditionally inverted addition with a 2-bit sign window."""
    C.x(ctrl)                                              # mask = NOT ctrl: subtract when ctrl = 0
    for w in src:
        C.cx(ctrl, w)
    Y = [_W(w) for w in src] + [_W(ctrl), _W(ctrl)]
    _add(C, Y, A[off:off + len(src) + 2], cin=_W(ctrl))
    for w in src:
        C.cx(ctrl, w)
    C.x(ctrl)


def _sd_mul(C, x, y):
    """-> 2n wires holding x*y; x, y disjoint registers of n wires (restored)."""
    n = len(x)
    assert len(y) == n and not (set(x) & set(y))
    A = [C.alloc() for _ in range(n + 1)]
    for i in range(n):
        top = C.alloc()
        C.cx(A[-1], top)                                   # sign extension: A now has n+i+2 bits
        A.append(top)
        _signed_row(C, A, i, x, y[i])
    # A (2n+1 bits) = sum s_i 2^i x;  + 2^n x - x = 2xy >= 0
    _add(C, [_W(w) for w in x], A[n:2 * n + 1])
    _sub(C, [_W(w) for w in x], A[:2 * n + 1])
    C.release(A[0])                                        # 2xy is even: the low bit is 0
    return A[1:2 * n + 1]


def _sd_sq(C, x):
    """-> 2n wires holding x^2."""
    n = len(x)
    A = [C.alloc() for _ in range(2)]
    for i in range(1, n):
        while len(A) < 2 * i + 2:
            t = C.alloc()
            C.cx(A[-1], t)
            A.append(t)
        _signed_row(C, A, i, x[:i], x[i])
    while len(A) < 2 * n:
        t = C.alloc()
        C.cx(A[-1], t)
        A.append(t)
    _add(C, [_W(w) for w in x], A[n:2 * n])
    D = []
    for i in range(n):
        D += [_W(x[i]), ZERO]
    _sub(C, D[:len(A)], A[:2 * n])
    return A[:2 * n]


def _half_sum(C, x0, x1):
    L = max(len(x0), len(x1)) + 1
    s = [C.alloc() for _ in range(L)]
    lo, hi = (x0, x1) if len(x0) >= len(x1) else (x1, x0)
    for w, t in zip(lo, s):
        C.cx(w, t)
    _add(C, [_W(w) for w in hi], s)
    return s


def _kara(C, a, b, L, base, sq, gar):
    m = len(a)
    if L == 0 or m <= base:
        return _sd_sq(C, a) if sq else _sd_mul(C, a, b)
    h = m // 2
    z0 = _kara(C, a[:h], None if sq else b[:h], L - 1, base, sq, gar)
    z2 = _kara(C, a[h:], None if sq else b[h:], L - 1, base, sq, gar)
    out = list(z0) + list(z2)
    sa = _half_sum(C, a[:h], a[h:])
    if sq:
        T = _kara(C, sa, None, L - 1, base, True, gar)
        gar += sa
    else:
        sb = _half_sum(C, b[:h], b[h:])
        T = _kara(C, sa, sb, L - 1, base, False, gar)
        gar += sa + sb
    _sub(C, [_W(w) for w in z0], T)
    _sub(C, [_W(w) for w in z2], T)
    _add(C, [_W(w) for w in T[:m + 1]], out[h:2 * m])     # T = a0 b1 + a1 b0 < 2^(m+1)
    gar += list(T)
    return out


# ----------------------------------------------------------------------------------------- result object

class MulSqCircuit:
    """A generated circuit.  Attributes:
        kind 'mul' | 'sq', n, karatsuba_levels, base, clean (Bennett version)
        ops        list of ('x', q) ('cx', c, t) ('and', a, b, t) ('erase', t, a, b, lin)
        num_qubits
        x, y, out  wire lists (y == [] for squaring), little endian
        garbage    wires left dirty at the end (Karatsuba half sums and middle products); [] for L = 0 or clean
        ancillas   all other wires: start and end at |0>
    """

    def __init__(self, kind, n, levels, base, clean, ops, nq, x, y, out, garbage):
        self.kind, self.n, self.karatsuba_levels, self.base, self.clean = kind, n, levels, base, clean
        self.ops, self.num_qubits = ops, nq
        self.x, self.y, self.out, self.garbage = list(x), list(y), list(out), list(garbage)
        used = set(self.x) | set(self.y) | set(self.out) | set(self.garbage)
        self.ancillas = [w for w in range(nq) if w not in used]

    def __repr__(self):
        k = self.counts()
        return 'MulSqCircuit(%s, n=%d, L=%d, base=%d, clean=%s: AND=%d, qubits=%d, garbage=%d)' % (
            self.kind, self.n, self.karatsuba_levels, self.base, self.clean, k['AND'], k['qubits'], k['garbage'])

    def counts(self):
        """Gate counts of the logical circuit (no simulation).  AND = Toffoli count; T = 4 AND; CX counts the
        explicit CX only (the emitted AND gadget adds 6 CX, each erasure one CZ)."""
        k = {'x': 0, 'cx': 0, 'and': 0, 'erase': 0}
        for o in self.ops:
            k[o[0]] += 1
        return dict(AND=k['and'], ERASE=k['erase'], T=4 * k['and'], CX=k['cx'], X=k['x'], qubits=self.num_qubits,
                    garbage=len(self.garbage))

    # ------------------------------------------------------------- checking
    def simulate(self, xs, ys=None):
        """Bit-sliced classical run on the inputs xs (and ys for 'mul'; ignored for 'sq').  Checks every AND target
        is |0>, every erasure condition, inputs restored, ancillas back to |0>.  -> list of outputs (ints).
        Raises CheckError on a violation."""
        S = len(xs)
        if self.kind == 'mul':
            if ys is None or len(ys) != S:
                raise ValueError('ys required for multiply')
        M = (1 << S) - 1
        v = [0] * self.num_qubits
        regs = [(self.x, xs)] + ([(self.y, ys)] if self.kind == 'mul' else [])
        for wires, vals in regs:
            for i, w in enumerate(wires):
                m = 0
                for s, val in enumerate(vals):
                    if val >> i & 1:
                        m |= 1 << s
                v[w] = m
        init = list(v)
        for i, o in enumerate(self.ops):
            k = o[0]
            if k == 'x':
                v[o[1]] ^= M
            elif k == 'cx':
                v[o[2]] ^= v[o[1]]
            elif k == 'and':
                if v[o[3]]:
                    raise CheckError('op %d: AND target not |0>' % i)
                v[o[3]] = v[o[1]] & v[o[2]]
            else:
                t, a, b, lin = o[1:]
                val = v[a] & v[b]
                for w in lin:
                    val ^= v[w]
                if v[t] != val:
                    raise CheckError('op %d: erase condition violated' % i)
                v[t] = 0
        for w in self.x + self.y:
            if v[w] != init[w]:
                raise CheckError('input wire %d changed' % w)
        for w in self.ancillas:
            if v[w]:
                raise CheckError('ancilla %d not returned to |0>' % w)
        out = [0] * S
        for i, w in enumerate(self.out):
            m = v[w]
            for s in range(S):
                if m >> s & 1:
                    out[s] |= 1 << i
        return out

    # ------------------------------------------------------------- output
    def to_qasm(self, num_qubits=None):
        """OpenQASM 3 text in the qoptct dialect (logical AND gadget of Gidney 2018 Fig. 3; erasure = H, measure,
        classically controlled CZ / Z / X).  num_qubits: pad the register (e.g. to match a reference)."""
        N = max(self.num_qubits, num_qubits or 0)
        L = []
        nb = 0
        for o in self.ops:
            k = o[0]
            if k == 'x':
                L.append('x q[%d];' % o[1])
            elif k == 'cx':
                L.append('cx q[%d], q[%d];' % o[1:])
            elif k == 'and':
                a, b, t = o[1:]
                L += ['h q[%d];' % t, 't q[%d];' % t, 'cx q[%d], q[%d];' % (a, t), 'cx q[%d], q[%d];' % (b, t),
                      'cx q[%d], q[%d];' % (t, a), 'cx q[%d], q[%d];' % (t, b), 'tdg q[%d];' % a, 'tdg q[%d];' % b,
                      't q[%d];' % t, 'cx q[%d], q[%d];' % (t, a), 'cx q[%d], q[%d];' % (t, b), 'h q[%d];' % t,
                      's q[%d];' % t]
            else:
                t, a, b, lin = o[1:]
                L.append('h q[%d];' % t)
                L.append('c[%d] = measure q[%d];' % (nb, t))
                body = ['cz q[%d], q[%d];' % (a, b)] + ['z q[%d];' % w for w in lin] + ['x q[%d];' % t]
                L.append('if (c[%d]) {%s}' % (nb, ' '.join(body)))
                nb += 1
        return ('OPENQASM 3.0;\ninclude "stdgates.inc";\nqubit[%d] q;\nbit[%d] c;\n' % (N, max(nb, 1))
                + '\n'.join(L) + '\n')

    def reference_qc(self):
        """independent schoolbook reference (.qc text) on this circuit's wire layout; see reference_qc()."""
        return reference_qc(self.kind, self.n, self.num_qubits)

    def reference_width(self):
        return _ref_width(self.kind, self.n, self.num_qubits)


# ----------------------------------------------------------------------------------------- generators

def _finish(kind, n, levels, base, clean, C, x, y, out, gar):
    gar = [w for w in dict.fromkeys(gar) if w not in set(out)]
    ops = C.ops
    if clean:
        # Bennett: forward, copy the product into a fresh register, inverse of the forward circuit
        out2 = [C.nw + i for i in range(2 * n)]
        C.nw += 2 * n
        ops = ops + [('cx', a, b) for a, b in zip(out, out2)] + _inverse(ops)
        out, gar = out2, []
    first = list(x) + list(y) + list(out) + list(gar)
    order = first + [w for w in range(C.nw) if w not in set(first)]
    mp = {w: i for i, w in enumerate(order)}
    new = []
    for o in ops:
        k = o[0]
        if k == 'x':
            new.append(('x', mp[o[1]]))
        elif k == 'cx':
            new.append(('cx', mp[o[1]], mp[o[2]]))
        elif k == 'and':
            new.append(('and', mp[o[1]], mp[o[2]], mp[o[3]]))
        else:
            new.append(('erase', mp[o[1]], mp[o[2]], mp[o[3]], tuple(mp[w] for w in o[4])))
    f = lambda ws: [mp[w] for w in ws]
    return MulSqCircuit(kind, n, levels, base, clean, new, C.nw, f(x), f(y), f(out), f(gar))


def _check_args(n, levels, base):
    if int(n) != n or n < 1:
        raise ValueError('n must be a positive integer')
    if int(levels) != levels or levels < 0:
        raise ValueError('karatsuba_levels must be >= 0')
    if int(base) != base or base < 1:
        raise ValueError('base must be >= 1')


def multiply(n, karatsuba_levels=0, base=4, clean=False):
    """Out-of-place multiplier (x, y, 0) -> (x, y, x*y) for n-bit x, y (2n-bit product).
    karatsuba_levels: Karatsuba recursion depth L (0 = base construction, n^2 + 4n AND, no garbage); the recursion
    stops early at width <= base.  clean=True: Bennett version without garbage (about twice the AND count of the
    forward circuit for L >= 1).  -> MulSqCircuit."""
    _check_args(n, karatsuba_levels, base)
    C = _Circ()
    x = C.reg(n)
    y = C.reg(n)
    gar = []
    out = _kara(C, x, y, karatsuba_levels, base, False, gar)
    return _finish('mul', n, karatsuba_levels, base, clean, C, x, y, out, gar)


def square(n, karatsuba_levels=0, base=4, clean=False):
    """Out-of-place square (x, 0) -> (x, x^2) for n-bit x (2n-bit result).  Base: (n^2 + 7n - 6)/2 AND, no garbage.
    Parameters as in multiply().  -> MulSqCircuit."""
    _check_args(n, karatsuba_levels, base)
    C = _Circ()
    x = C.reg(n)
    gar = []
    out = _kara(C, x, None, karatsuba_levels, base, True, gar)
    return _finish('sq', n, karatsuba_levels, base, clean, C, x, [], out, gar)


# ------------------------------------------------------------------------- independent reference + proof

def _ref_width(kind, n, num_qubits):
    io = (4 if kind == 'mul' else 3) * n
    return max(num_qubits, io + 2 * n)


def reference_qc(kind, n, num_qubits=0):
    """Independent schoolbook reference as a Feynman .qc text (Toffoli gates only), on the layout of MulSqCircuit:
    x = q0..q(n-1), [y = qn..q(2n-1)], out (2n wires), the remaining wires clean.  For every bit i of the
    multiplier: tmp_j = x_j AND y_i (for squaring y = x, and tmp_i = x_i by CX), then out[i .. i+n] += tmp with a
    ripple-carry adder (carry c_{k+1} = MAJ(tmp_k, out_k, c_k) by three Toffolis, Vedral-style), then the carries
    and tmp are uncomputed.  Its tmp and carry wires are the last 2n wires of max(num_qubits, io + 2n) wires.
    It shares no code with the generators."""
    N = _ref_width(kind, n, num_qubits)
    q = ['q%d' % i for i in range(N)]
    xs = q[:n]
    ys = q[n:2 * n] if kind == 'mul' else xs
    o0 = 2 * n if kind == 'mul' else n
    out = q[o0:o0 + 2 * n]
    tmp = q[N - 2 * n:N - n]
    car = q[N - n:N]                       # car[k] holds c_{k+1}
    ins = xs + (ys if kind == 'mul' else [])
    L = ['.v ' + ' '.join(q), '.i ' + ' '.join(ins), '', 'BEGIN']

    def partial(i):
        for j in range(n):
            if ys[i] == xs[j]:
                L.append('tof %s %s' % (xs[j], tmp[j]))
            else:
                L.append('tof %s %s %s' % (xs[j], ys[i], tmp[j]))

    def maj(k, b):
        a, c1 = tmp[k], car[k]
        L.append('tof %s %s %s' % (a, b[k], c1))
        if k >= 1:
            L.append('tof %s %s %s' % (a, car[k - 1], c1))
            L.append('tof %s %s %s' % (b[k], car[k - 1], c1))

    for i in range(n):
        partial(i)
        b = out[i:i + n + 1]
        for k in range(n):
            maj(k, b)
        L.append('tof %s %s' % (car[n - 1], b[n]))
        for k in range(n - 1, -1, -1):
            maj(k, b)
            L.append('tof %s %s' % (tmp[k], b[k]))
            if k >= 1:
                L.append('tof %s %s' % (car[k - 1], b[k]))
        partial(i)
    L.append('END')
    return '\n'.join(L) + '\n'


def prove(circ, budget_sec=600.0, max_width=22, engine='prove'):
    """Exact proof (clean model: all inputs, all measurement outcomes) of circ against the independent
    reference_qc().  Only garbage-free circuits can be proved (karatsuba_levels = 0 or clean=True).
    engine='prove': qoptct.verify.prove (path-sum normal form / trace miter; in practice multiply n <= 5,
    square n <= 6 within minutes); engine='kprove': the BDD prover qoptct.verify.kprove when the installed qoptct
    has it (ImportError otherwise; proves multiply n = 6..8 in seconds).
    -> dict with 'status' PROVED / REFUTED / UNKNOWN (and the prover's details)."""
    if circ.garbage:
        raise ValueError('circuit has garbage wires: use karatsuba_levels=0 or clean=True')
    if engine == 'kprove':
        from qoptct.verify.kprove import kprove as _kp
        verify = None
    elif engine == 'prove':
        from qoptct import verify
    else:
        raise ValueError('engine must be "prove" or "kprove"')
    N = circ.reference_width()
    d = tempfile.mkdtemp(prefix='qopt_mulsq_')
    qa, qc = os.path.join(d, 'c.qasm'), os.path.join(d, 'ref.qc')
    try:
        with open(qa, 'w') as fh:
            fh.write(circ.to_qasm(num_qubits=N))
        with open(qc, 'w') as fh:
            fh.write(circ.reference_qc())
        if verify is None:
            return _kp(qa, qc, clean=True, budget_sec=budget_sec, engine='bdd', sample=False)
        return verify.prove(qa, qc, clean=True, budget_sec=budget_sec, max_width=max_width)
    finally:
        for p in (qa, qc):
            try:
                os.unlink(p)
            except OSError:
                pass
        try:
            os.rmdir(d)
        except OSError:
            pass
