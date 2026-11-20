# SPDX-License-Identifier: Apache-2.0
"""modadd: modular adders / subtractors, doubling and comparators on logical ANDs (ported from the qoptct lab, agent a3;
checked here again exactly on every input of the domain for small n, by qoptct kprove and by canaries).

    A = mod_add(p, n=None)              (x, y) -> (x, (x + y) mod p)         in place, any p (domain x, y < p)
    A = mod_sub(p, n=None)              (x, y) -> (x, (y - x) mod p)         adjoint of mod_add (same ANDs)
    A = mod_add_out(p, n=None)          (x, y, 0) -> (x, y, (x + y) mod p)   p odd, out of place
    A = mod_add_mersenne(n)             (x, y) -> (x, (x + y) mod (2^n - 1)) in place
    A = mod_add_const(p, a, n=None)     y -> (y + a) mod p                   classical a, 0 < a < p
    A = mod_sub_const(p, a, n=None)     y -> (y - a) mod p                   (adjoint)
    A = mod_double(p, n=None)           y -> 2y mod p                        p odd; result on wires of layout['z']
    A = compare_geq(n)                  (x, y, t) -> (x, y, t ^ [x >= y])    n-bit comparator
    A = compare_const_geq(n, c)         (x, t)    -> (x, t ^ [x >= c])

Every builder returns a ModCircuit: .qasm(), .counts(), .n_and, .layout (input / output registers as wire numbers of
the emitted OpenQASM 3), .reference(point), .check(points) -> dict(ok, why): exact bit-sliced simulation of the gate list
(andcirc.simulate: the AND targets are 0, every erasure is consistent with its value for the whole batch, outputs equal
the reference, ALL other wires end at 0), .domain(): all valid inputs for small n.

COST MODEL. Gidney 2018: AND = 4 T + 3 CX; erasure by X-measurement with a CZ correction = 0 T.

AND COUNTS (generated gate lists; the tests pin them; n = bit length of the registers, p < 2^n):
    mod_add / mod_sub                4n - 2   (n = 4, 5, 6, 8: 14 / 18 / 22 / 30)
    mod_add_out (p odd)              3n - 2
    mod_add_mersenne                 3n - 1   (flag = carry of x + y + 1, increment, add, the flag is erased by the last carry)
    mod_add_const (flag-first)       about 3n - 4, depends on the bits of a and p ((n, p, a) = (5, 29, 7): 11, (6, 61, 10): 14, (8, 251, 100): 19)
    mod_double (flag-first, p odd)   2n - 2
    compare_geq                      n   (carry chain only, the sum is not formed)
    compare_const_geq                n - 1 for c = 1 and fewer for a constant with trailing zero bits, 0 for c = 2^(n-1) or c > 2^n - 1
CONSTRUCTIONS. Textbook scheme (Beauregard / Vedral-Barenco-Ekert / Haener et al. / Litinski 2023, Fig. 6): y += x (n+1 bits),
y -= p, controlled y += p, comparator restores the flag: 4 additions (4n). The lab's savings (a3), all small:
  * the carry chain is built with the majority trick (one AND per bit, carries erased by measurement on the way back);
    the comparator that restores the flag is a carry chain only (no sum);
  * "flag-first": the flag f = [y >= p - a] (or [x + y >= p]) is computed first by a carry chain, then y -= f p, then y += a
    where the carry-out erases f for free (f sits in the top bit of an (n+1)-bit register): 3n instead of 4n for a
    classical addend and 2n - 2 for the doubling (f is the low bit of 2y mod p);
  * out of place (p odd): the flag is erased by the parity z_0 ^ x_0 ^ y_0 (x + y and x + y - p have different parity).
Against the published adders of the same model (4n for the modular adder of Litinski 2023) the gain is 2 ANDs for the general
adder and subtractor (the subtractor is the adjoint: 4n - 2, against 6n for subtraction built from published pieces); the constant
adder and the doubling save about n and 2 respectively. These are small constant-factor savings, not records.
DOMAIN. Inputs must be reduced (< p); behaviour outside is not defined (and not checked). Garbage-free: every ancilla ends at 0.
"""
import itertools
import random

import numpy as np

from .andcirc import Circ as _Circ, simulate, bits_of, ints_of, text_counts

__all__ = ['ModCircuit', 'mod_add', 'mod_sub', 'mod_add_out', 'mod_add_mersenne', 'mod_add_const', 'mod_sub_const',
           'mod_double', 'compare_geq', 'compare_const_geq']


# ------------------------------------------------------------------ adaptor: literals (wire, inverted) on andcirc.Circ
def lit(w, inv=0):
    return (w, inv)


def Ls(ws, inv=0):
    return [lit(w, inv) for w in ws]


def _kbits(v, n):
    return [(v >> i) & 1 for i in range(n)]


class _A:
    """thin layer over andcirc.Circ: AND / erase take literals (wire, inverted); inversions are free X gates."""

    def __init__(s, C):
        s.C = C

    def alloc(s):
        return s.C.new()[0]

    def free(s, q):
        s.C.release([q])

    def x(s, q):
        s.C.x(q)

    def cx(s, a, b):
        assert a != b
        s.C.cx(a, b)

    def AND(s, la, lb):
        (a, ia), (b, ib) = la, lb
        if ia:
            s.C.x(a)
        if ib:
            s.C.x(b)
        t = s.C.AND(a, b)
        if ia:
            s.C.x(a)
        if ib:
            s.C.x(b)
        return t

    def erase(s, t, la, lb, extra=()):
        (a, ia), (b, ib) = la, lb
        if ia:
            s.C.x(a)
        if ib:
            s.C.x(b)
        s.C.MZ(t, [('cz', a, b)] + [('z', e) for e in extra])
        if ia:
            s.C.x(a)
        if ib:
            s.C.x(b)
        s.C.release([t])


class _Chain:
    """Ripple adder y += k (+ cin) with logical ANDs (majority trick, parallel CX 'frames'). Y: wires (LSB first), K: list of
    addends (0, 1 or a literal), cin: 0 / 1 / literal. need_cout: also compute the carry out of the top bit.
    Phase 1 (__init__): carries into fresh ancillas. Phase 2 (finish): erase the carries by measurement from the top down,
    write the sum into Y (write=True) or into sum_target, copy the carry out to cout_target."""

    def __init__(s, C, Y, K, cin=0, need_cout=False):
        s.C = C; s.Y = Y; s.K = K; n = len(Y); s.n = n; s.fr = {}
        assert len(K) == n
        s.cin = [None] * (n + 1); s.rec = [None] * n; s.need_cout = need_cout
        c = cin
        for i in range(n):
            s.cin[i] = c
            if i == n - 1 and not need_cout:
                break
            c, s.rec[i] = s._maj(i, c)
        s.cin[n] = c if need_cout else None

    def frame(s, w, cw):
        cur = s.fr.get(w)
        if cur == cw:
            return
        if cur is not None:
            s.C.cx(cur, w)
        if cw is not None:
            s.C.cx(cw, w)
        s.fr[w] = cw

    def _maj(s, i, c):
        C = s.C; y = s.Y[i]; k = s.K[i]
        ki = isinstance(k, int); ci = isinstance(c, int)
        if ki and ci:
            v = k + c
            return ({0: 0, 1: lit(y), 2: 1}[v], None)
        if ki or ci:
            v, L = (k, c) if ki else (c, k)
            s.frame(y, None); s.frame(L[0], None)
            if v == 0:
                la, lb = lit(y), lit(L[0], L[1]); t = C.AND(la, lb); return (lit(t, 0), ('and', t, la, lb))
            la, lb = lit(y, 1), lit(L[0], 1 ^ L[1]); t = C.AND(la, lb); return (lit(t, 1), ('and', t, la, lb))
        wc, ic = c; wk, ik = k
        s.frame(y, wc); s.frame(wk, wc)
        la, lb = lit(y, ic), lit(wk, ik ^ ic)
        t = C.AND(la, lb); C.cx(wc, t)
        return (lit(t, ic), ('maj', t, la, lb, wc))

    def release(s, w):
        for x in list(s.fr):
            if s.fr[x] == w and x != w:
                s.frame(x, None)

    def _erase(s, i):
        r = s.rec[i]
        if r is None:
            return
        C = s.C
        s.release(r[1])
        if r[0] == 'and':
            _, t, la, lb = r
            s.frame(la[0], None); s.frame(lb[0], None)
            C.erase(t, la, lb)
        else:
            _, t, la, lb, wc = r
            C.cx(wc, t)
            s.frame(la[0], wc); s.frame(lb[0], wc)
            C.erase(t, la, lb)

    def _copy_lit(s, l, tgt):
        C = s.C
        if isinstance(l, int):
            if l:
                C.x(tgt)
        else:
            C.cx(l[0], tgt)
            if l[1]:
                C.x(tgt)

    def finish(s, write=True, cout_target=None, sum_target=None):
        C = s.C; n = s.n
        if s.need_cout and cout_target is not None:
            s._copy_lit(s.cin[n], cout_target)
        for i in reversed(range(n)):
            s._erase(i)
            y, k, c = s.Y[i], s.K[i], s.cin[i]
            if write:
                s.release(y)
                flip = 0
                if isinstance(c, int):
                    flip ^= c
                else:
                    s.frame(y, c[0]); flip ^= c[1]
                if isinstance(k, int):
                    flip ^= k
                else:
                    s.frame(k[0], None); C.cx(k[0], y); flip ^= k[1]
                if flip:
                    C.x(y)
                s.fr[y] = None
            elif sum_target is not None:
                T = sum_target[i]
                s.frame(y, None); C.cx(y, T)
                if isinstance(k, int):
                    if k:
                        C.x(T)
                else:
                    s.frame(k[0], None); C.cx(k[0], T)
                    if k[1]:
                        C.x(T)
                s._copy_lit(c, T)
        for w in list(s.fr):
            s.frame(w, None)


def _add(C, Y, K, cin=0, cout_target=None):
    """Y += K (mod 2^len(Y)); with cout_target: it receives (xor) the carry out."""
    _Chain(C, Y, K, cin, need_cout=cout_target is not None).finish(write=True, cout_target=cout_target)


def _carry(C, Y, K, cin, target):
    """target ^= carry out of Y + K + cin; Y unchanged (the carries are erased)."""
    _Chain(C, Y, K, cin, need_cout=True).finish(write=False, cout_target=target)


# ------------------------------------------------------------------ the circuit object
class ModCircuit:
    """A generated circuit. Attributes: kind, p, n, circ (andcirc.Circ), n_and, layout = dict(inputs=..., outputs=...)
    (register name -> list of wire numbers of the emitted QASM), nref (wires that carry inputs / outputs; every other
    wire starts and ends at 0), order (internal wire ids in QASM order)."""

    def __init__(s, kind, p, n, circ, inputs, outputs, ref, extra_order=(), dom=None, consumed=()):
        s.kind, s.p, s.n, s.circ, s.n_and = kind, p, n, circ, circ.nand
        s._in, s._out, s._ref, s._dom = inputs, outputs, ref, dom
        s._consumed = tuple(consumed)
        s.a = None                                  # the classical constant of the *_const builders
        order = []
        for ws in list(inputs.values()) + list(outputs.values()) + [list(extra_order)]:
            for w in ws:
                if w not in order:
                    order.append(w)
        s.order = order; s.nref = len(order)
        pos = {w: i for i, w in enumerate(order)}
        s.layout = dict(inputs={k: [pos[w] for w in v] for k, v in inputs.items()},
                        outputs={k: [pos[w] for w in v] for k, v in outputs.items()})

    def counts(s):
        return s.circ.counts()

    def qasm(s):
        return s.circ.qasm(s.order)

    def text_counts(s):
        return text_counts(s.qasm())

    def reference(s, point):
        """independent python-int formula: dict output register -> value for the input point (dict register -> int)."""
        return s._ref(point)

    def domain(s):
        """all valid input points (small n only: product of the register ranges)."""
        return s._dom()

    def check(s, points):
        """exact simulation on `points` (list of dicts register -> int): outputs equal reference(), all other wires 0,
        every erasure consistent for the whole batch. -> dict(ok, why, bad)."""
        points = list(points)
        B = len(points)
        X = {}
        for name, ws in s._in.items():
            arr = bits_of([pt[name] for pt in points], len(ws))
            for j, w in enumerate(ws):
                X[w] = arr[j]
        v, ok, why = simulate(s.circ, X)
        if not ok:
            return dict(ok=False, why=why, bad=None)
        refs = [s._ref(pt) for pt in points]
        used = set()
        for name, ws in s._out.items():
            used |= set(ws)
            got = ints_of(v, ws)
            for pt, r, g in zip(points, refs, got):
                if g != r[name]:
                    return dict(ok=False, why='register %s, input %s: %d != %d' % (name, pt, g, r[name]), bad=pt)
        for name, ws in s._in.items():
            if name not in s._out and name not in s._consumed:      # an input that is not an output must be unchanged
                used |= set(ws)
                got = ints_of(v, ws)
                for pt, g in zip(points, got):
                    if g != pt[name]:
                        return dict(ok=False, why='input register %s changed' % name, bad=pt)
        rest = sorted(set(range(s.circ.N)) - used)
        if rest and v[rest].any():
            return dict(ok=False, why='a wire is not 0 at the end', bad=None)
        return dict(ok=True, why='', bad=None)


def _nbits(p, n):
    if p < 2:
        raise ValueError('need p >= 2')
    m = (p - 1).bit_length() if p > 2 else 1
    if n is None:
        n = max(m, 1)
    if p > (1 << n):
        raise ValueError('p must be <= 2^n')
    return n


def _dom2(p):
    return lambda: [dict(x=x, y=y) for x in range(p) for y in range(p)]


def _inverse(keep):
    """adjoint of a gate list of cx / x / AND / erase with plain CZ correction: reversed order, AND <-> MZ. `keep`: the original andcirc.Circ. -> new Circ with the same wires."""
    D = _Circ(); D.N = keep.N; D.free = []; D.live = 0; D.peak = keep.peak
    for op in reversed(keep.ops):
        k = op[0]
        if k in ('x', 'cx'):
            D.ops.append(op)
        elif k == 'AND':
            a, b, t = op[1:]
            D.ops.append(('MZ', t, (('cz', a, b),)))
        elif k == 'MZ':
            t, corr = op[1], op[2]
            if len(corr) != 1 or corr[0][0] != 'cz':
                raise ValueError('cannot invert an erasure with extra corrections')
            D.ops.append(('AND', corr[0][1], corr[0][2], t))
        else:
            raise ValueError('cannot invert op ' + k)
    D.nand = sum(1 for o in D.ops if o[0] == 'AND')
    return D


def _mk(n, names):
    C = _Circ(); A = _A(C)
    regs = {nm: C.new(1 if nm == 't' else n) for nm in names}       # 't' = a one-bit flag register
    return C, A, regs


# ------------------------------------------------------------------ builders
def _build_modadd(p, n):
    C, A, R = _mk(n, ('x', 'y')); X, Y = R['x'], R['y']
    m = A.alloc(); Yx = Y + [m]
    _add(A, Yx, Ls(X) + [0])                                   # y + x  (n + 1 bits)
    _add(A, Yx, _kbits((1 << (n + 1)) - p, n + 1))             # - p
    _add(A, Y, [lit(m) if b else 0 for b in _kbits(p, n)])     # + p if the top bit says negative
    _carry(A, Y, Ls(X, 1), 1, m)                               # m ^= [z >= x]: restores the flag to 0
    A.free(m)
    return C, X, Y


def mod_add(p, n=None):
    """(x, y) -> (x, (x + y) mod p), x, y < p. 4n - 2 ANDs. Wires: x = 0..n-1, y = n..2n-1."""
    n = _nbits(p, n)
    C, X, Y = _build_modadd(p, n)
    return ModCircuit('mod_add', p, n, C, dict(x=X, y=Y), dict(y=Y), lambda pt: dict(y=(pt['x'] + pt['y']) % p),
                      dom=_dom2(p))


def mod_sub(p, n=None):
    """(x, y) -> (x, (y - x) mod p), x, y < p: the adjoint of mod_add (ANDs and erasures exchange roles). 4n - 2 ANDs."""
    n = _nbits(p, n)
    C, X, Y = _build_modadd(p, n)
    D = _inverse(C)
    return ModCircuit('mod_sub', p, n, D, dict(x=X, y=Y), dict(y=Y), lambda pt: dict(y=(pt['y'] - pt['x']) % p),
                      dom=_dom2(p))


def mod_add_out(p, n=None):
    """(x, y, 0) -> (x, y, (x + y) mod p), p odd, x, y < p. 3n - 2 ANDs. Wires: x, y, z (n bits each)."""
    if p % 2 == 0:
        raise ValueError('out-of-place version needs an odd p')
    n = _nbits(p, n)
    C, A, R = _mk(n, ('x', 'y')); X, Y = R['x'], R['y']
    Z = C.new(n + 1)
    ch = _Chain(A, Y, Ls(X), 0, need_cout=True); ch.finish(write=False, cout_target=Z[n], sum_target=Z[:n])
    _add(A, Z, _kbits((1 << (n + 1)) - p, n + 1))
    _add(A, Z[:n], [lit(Z[n]) if b else 0 for b in _kbits(p, n)])
    s = Z[n]; C.cx(Z[0], s); C.cx(X[0], s); C.cx(Y[0], s); C.x(s)           # the flag equals z0 ^ x0 ^ y0 ^ 1
    A.free(s)
    return ModCircuit('mod_add_out', p, n, C, dict(x=X, y=Y), dict(x=X, y=Y, z=Z[:n]),
                      lambda pt: dict(x=pt['x'], y=pt['y'], z=(pt['x'] + pt['y']) % p), dom=_dom2(p))


def mod_add_mersenne(n):
    """(x, y) -> (x, (x + y) mod (2^n - 1)), x, y < 2^n - 1. 3n - 1 ANDs."""
    if n < 2:
        raise ValueError('need n >= 2')
    p = (1 << n) - 1
    C, A, R = _mk(n, ('x', 'y')); X, Y = R['x'], R['y']
    f = A.alloc()
    _carry(A, Y, Ls(X), 1, f)                                  # f = carry of x + y + 1 = [x + y >= 2^n - 1]
    _add(A, Y, [0] * n, lit(f))                                # y += f
    _add(A, Y + [f], Ls(X) + [0])                              # (y, f) += x: the carry out erases f
    A.free(f)
    return ModCircuit('mod_add_mersenne', p, n, C, dict(x=X, y=Y), dict(y=Y), lambda pt: dict(y=(pt['x'] + pt['y']) % p),
                      dom=_dom2(p))


def mod_add_const(p, a, n=None):
    """y -> (y + a) mod p for a classical a, 0 < a < p, y < p; flag-first (about 3n ANDs). Wire: y = 0..n-1."""
    n = _nbits(p, n)
    if not 0 < a < p:
        raise ValueError('need 0 < a < p')
    C, A, R = _mk(n, ('y',)); Y = R['y']
    f = A.alloc()
    _carry(A, Y, _kbits((1 << n) - (p - a), n), 0, f)          # f = [y >= p - a]
    _add(A, Y, [lit(f) if b else 0 for b in _kbits((1 << n) - p, n)])      # y -= f p
    _add(A, Y + [f], _kbits(a, n) + [0])                       # y += a, the carry out erases f
    A.free(f)
    mc = ModCircuit('mod_add_const', p, n, C, dict(y=Y), dict(y=Y), lambda pt: dict(y=(pt['y'] + a) % p),
                    dom=lambda: [dict(y=y) for y in range(p)])
    mc.a = a
    return mc


def mod_sub_const(p, a, n=None):
    """y -> (y - a) mod p (adjoint of mod_add_const)."""
    c = mod_add_const(p, a, n)
    D = _inverse(c.circ)
    mc = ModCircuit('mod_sub_const', p, c.n, D, c._in, c._out, lambda pt: dict(y=(pt['y'] - a) % p), dom=c._dom)
    mc.a = a
    return mc


def mod_double(p, n=None):
    """y -> 2y mod p for odd p, y < p. 2n - 2 ANDs. Input y on wires 0..n-1; the result is the register 'z' =
    [flag, y_0, ..., y_{n-2}] (the top input wire ends at 0)."""
    if p % 2 == 0:
        raise ValueError('doubling needs an odd p')
    n = _nbits(p, n)
    h = (p + 1) // 2
    C, A, R = _mk(n, ('y',)); Y = R['y']
    f = A.alloc()
    _carry(A, Y, _kbits((1 << n) - h, n), 0, f)                # f = [y >= h]
    _add(A, Y, [lit(f) if b else 0 for b in _kbits((1 << n) - h, n)])      # y -= f h; 2y - f p = f + 2 (y - f h)
    A.free(Y[-1])
    return ModCircuit('mod_double', p, n, C, dict(y=Y), dict(z=[f] + Y[:-1]), lambda pt: dict(z=(2 * pt['y']) % p),
                      dom=lambda: [dict(y=y) for y in range(p)], consumed=('y',))


def compare_geq(n):
    """(x, y, t) -> (x, y, t ^ [x >= y]), n-bit x, y. Wires: x, y, t."""
    if n < 1:
        raise ValueError('need n >= 1')
    C, A, R = _mk(n, ('x', 'y', 't')); X, Y, T = R['x'], R['y'], R['t']
    _carry(A, X, Ls(Y, 1), 1, T[0])
    return ModCircuit('compare_geq', None, n, C, dict(x=X, y=Y, t=T), dict(x=X, y=Y, t=T),
                      lambda pt: dict(x=pt['x'], y=pt['y'], t=pt['t'] ^ (pt['x'] >= pt['y'])),
                      dom=lambda: [dict(x=x, y=y, t=t) for x in range(1 << n) for y in range(1 << n) for t in (0, 1)])


def compare_const_geq(n, c):
    """(x, t) -> (x, t ^ [x >= c]) for a classical 0 <= c < 2^n + 1. Wires: x, t."""
    if n < 1 or not 0 <= c <= (1 << n):
        raise ValueError('need n >= 1 and 0 <= c <= 2^n')
    C, A, R = _mk(n, ('x', 't')); X, T = R['x'], R['t']
    if c == 0:
        C.x(T[0])
    elif c < (1 << n):
        _carry(A, X, _kbits((1 << n) - c, n), 0, T[0])
    return ModCircuit('compare_const_geq', None, n, C, dict(x=X, t=T), dict(x=X, t=T),
                      lambda pt: dict(x=pt['x'], t=pt['t'] ^ (pt['x'] >= c)),
                      dom=lambda: [dict(x=x, t=t) for x in range(1 << n) for t in (0, 1)])
