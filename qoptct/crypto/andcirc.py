# SPDX-License-Identifier: Apache-2.0
"""andcirc: a small engine for classical reversible circuits on logical ANDs (Gidney 2018) with measurement-based
uncomputation (ported from the qoptct lab, agent k2).

Ops (wires are ints):
    ('x', a) ('cx', a, b) ('z', a) ('cz', a, b)
    ('AND', a, b, t)        t (must be |0>) := a AND b                         -- 4 T + 3 CX (cost model)
    ('MZ', t, corr)         erase t by an X-basis measurement; on outcome 1 apply corr (a tuple of ('cz', a, b) /
                            ('z', a)) and x t                                   -- 0 T
Addends of Circ.add / Circ.sub are lists whose items are wires or the constants ZERO / ONE ('0' / '1').

simulate(C, X): bit-sliced classical simulation of a batch of inputs with the exact semantics of the gadgets:
    an AND target must be 0; at an MZ the value (correction XOR t) must be the SAME for every input of the batch
    (then the outcome-dependent phase is global, for every outcome); explicit z / cz accumulate a phase that must not
    depend on the input. On a batch containing ALL inputs of a block this is an exact check for every measurement
    outcome; on a sample it is a sampled check.
Circ.qasm(first) writes OpenQASM 3 in the dialect of qoptct.verify (AND gadget of Gidney 2018 fig. 3, `h; measure;
if (c[k]) {...}`), readable by qoptct.verify.arith.check_sim and kprove.
"""
import re

import numpy as np

ZERO, ONE = '0', '1'


def _k(x):
    if isinstance(x, str):
        return x
    if isinstance(x, (bool, np.bool_)):
        return ONE if x else ZERO
    return int(x)


def const_bits(K, n):
    """the constant K as an addend of n bits (LSB first)."""
    return [ONE if (K >> i) & 1 else ZERO for i in range(n)]


def tz(K):
    """number of trailing zero bits of K > 0."""
    return (K & -K).bit_length() - 1


class Circ:
    def __init__(s):
        s.N = 0; s.ops = []; s.free = []; s.nand = 0; s.peak = 0; s.live = 0

    # ---- wires
    def new(s, k=1):
        out = []
        for _ in range(k):
            if s.free:
                out.append(s.free.pop())
            else:
                out.append(s.N); s.N += 1
        s.live += k; s.peak = max(s.peak, s.live)
        return out

    def release(s, ws):
        """return wires that are |0> to the pool (simulate() checks that freed ancillas end at 0)."""
        for w in ws:
            s.free.append(w)
        s.live -= len(ws)

    # ---- gates
    def x(s, a): s.ops.append(('x', a))
    def cx(s, a, b): s.ops.append(('cx', a, b))
    def z(s, a): s.ops.append(('z', a))
    def cz(s, a, b): s.ops.append(('cz', a, b))

    def AND(s, a, b, t=None):
        if t is None:
            t = s.new()[0]
        s.ops.append(('AND', a, b, t)); s.nand += 1
        return t

    def MZ(s, t, corr):
        s.ops.append(('MZ', t, tuple(corr)))

    def unAND(s, a, b, t):
        s.MZ(t, [('cz', a, b)])

    # ---- Gidney ripple-carry addition in place: acc += add (mod 2^W) [carry-out into cout]
    def add(s, acc, add, cout=None):
        """acc: wires (LSB first). add: same or shorter list of ZERO / ONE / wires (wires are not changed and may
        repeat). cout: a fresh |0> wire that receives the carry. ANDs = number of positions i in [s0, W-2]
        (W-1 with cout) where s0 is the lowest non-zero addend position, minus 1 if the addend there is ONE."""
        W = len(acc); b = [_k(v) for v in add] + [ZERO] * (W - len(add))
        assert len(b) == W
        s0 = next((i for i in range(W) if b[i] != ZERO), None)
        if s0 is None:
            return
        top = W - 1 if cout is None else W
        T = {}                                     # T[i] = wire of carry c_i
        for i in range(s0, top):
            c = T.get(i); bi = b[i]; a = acc[i]
            tgt = cout if (i + 1 == W) else None
            if c is None:
                if bi == ONE:
                    t = tgt if tgt is not None else s.new()[0]; s.cx(a, t); T[i + 1] = t
                elif bi == ZERO:
                    continue
                else:
                    T[i + 1] = s.AND(a, bi, tgt)
                continue
            s.cx(c, a)                             # a := acc ^ c (kept until the backward pass)
            op = s._opnd_set(bi, c)
            t = s.AND(a, op, tgt)
            s._opnd_reset(bi, c)
            s.cx(c, t); T[i + 1] = t
        i = W - 1
        if cout is None and T.get(i) is not None:
            s.cx(T[i], acc[i])
        s._xor_b(b[i], acc[i])
        for i in range(top - 1, s0 - 1, -1):
            if i == W - 1:
                continue                           # the carry-out is kept
            if i + 1 not in T:
                continue
            c = T.get(i); bi = b[i]; a = acc[i]; t = T[i + 1]
            if c is None:
                if bi == ONE:
                    s.cx(a, t); s.release([t]); s.x(a)
                else:
                    s.unAND(a, bi, t); s.release([t]); s.cx(bi, a)
                continue
            s.cx(c, t)
            op = s._opnd_set(bi, c)
            s.unAND(a, op, t); s.release([t])
            s._opnd_reset(bi, c)
            s._xor_b(bi, a)

    def _opnd_set(s, bi, c):
        if bi == ZERO:
            return c
        if bi == ONE:
            s.x(c); return c
        s.cx(c, bi); return bi

    def _opnd_reset(s, bi, c):
        if bi == ZERO:
            return
        if bi == ONE:
            s.x(c); return
        s.cx(c, bi)

    def _xor_b(s, bi, a):
        if bi == ZERO:
            return
        if bi == ONE:
            s.x(a); return
        s.cx(bi, a)

    def sub(s, acc, add, bout=None):
        """acc -= add (mod 2^W) as ~(~acc + add). bout: borrow (acc < add)."""
        for a in acc:
            s.x(a)
        s.add(acc, add, cout=bout)
        for a in acc:
            s.x(a)

    # ---- output
    def qasm(s, first=()):
        """OpenQASM 3 text. Wires in `first` get the numbers 0..len(first)-1 (in this order), the others follow."""
        L = []; nb = 0
        mp = {w: i for i, w in enumerate(first)}
        for w in range(s.N):
            if w not in mp:
                mp[w] = len(mp)

        def g(name, *q):
            L.append('%s %s;' % (name, ', '.join('q[%d]' % mp[v] for v in q)))
        for op in s.ops:
            k = op[0]
            if k in ('x', 'z'):
                g(k, op[1])
            elif k in ('cx', 'cz'):
                g(k, op[1], op[2])
            elif k == 'AND':
                a, b, t = op[1:]
                g('h', t); g('t', t); g('cx', a, t); g('cx', b, t); g('cx', t, a); g('cx', t, b)
                g('tdg', a); g('tdg', b); g('t', t); g('cx', t, a); g('cx', t, b); g('h', t); g('s', t)
            elif k == 'MZ':
                t, corr = op[1], op[2]
                g('h', t); L.append('c[%d] = measure q[%d];' % (nb, mp[t]))
                body = ' '.join(('cz q[%d], q[%d];' % (mp[c[1]], mp[c[2]])) if c[0] == 'cz' else ('z q[%d];' % mp[c[1]])
                                for c in corr)
                L.append('if (c[%d]) {%s x q[%d];}' % (nb, body, mp[t])); nb += 1
        return ('OPENQASM 3.0;\ninclude "stdgates.inc";\nqubit[%d] q;\nbit[%d] c;\n' % (s.N, max(nb, 1))) + \
            '\n'.join(L) + '\n'

    def counts(s):
        """cost model: AND = 4 T + 3 CX; explicit cx / cz and the cz of corrections count as CX. (The emitted AND
        gadget text has 6 cx lines per AND; see text_counts.)"""
        expl = sum(1 for o in s.ops if o[0] in ('cx', 'cz'))
        corr = sum(sum(1 for c in o[2] if c[0] == 'cz') for o in s.ops if o[0] == 'MZ')
        return dict(AND=s.nand, T=4 * s.nand, CX=3 * s.nand + expl + corr, qubits=s.N, peak=s.peak)


def text_counts(qasm):
    """gate counts of an emitted text: T (t + tdg), cx lines (incl. inside if), cz, measurements, qubits."""
    t = len(re.findall(r'^(?:t|tdg) q', qasm, re.M))
    return dict(T=t, AND=t // 4, cx=len(re.findall(r'\bcx q', qasm)), cz=len(re.findall(r'\bcz q', qasm)),
                measurements=qasm.count('measure'), qubits=int(qasm.split('qubit[', 1)[1].split(']', 1)[0]))


def simulate(C, X):
    """X: {wire: np.uint8 array (B,)} (other wires start at 0). -> (bits[N, B], ok, why)."""
    B = len(next(iter(X.values())))
    v = np.zeros((C.N, B), np.uint8)
    for w, a in X.items():
        v[w] = a
    ph = np.zeros(B, np.uint8)
    for op in C.ops:
        k = op[0]
        if k == 'x':
            v[op[1]] ^= 1
        elif k == 'cx':
            v[op[2]] ^= v[op[1]]
        elif k == 'z':
            ph ^= v[op[1]]
        elif k == 'cz':
            ph ^= v[op[1]] & v[op[2]]
        elif k == 'AND':
            a, b, t = op[1:]
            if t in (a, b):
                return v, False, 'AND target = control'
            if v[t].any():
                return v, False, 'AND target %d not 0' % t
            v[t] = v[a] & v[b]
        elif k == 'MZ':
            t, corr = op[1], op[2]; cp = np.zeros(B, np.uint8)
            for c in corr:
                if t in c[1:]:
                    return v, False, 'correction touches the measured wire'
                cp ^= (v[c[1]] & v[c[2]]) if c[0] == 'cz' else v[c[1]]
            d = cp ^ v[t]
            if not (d == d[0]).all():
                return v, False, 'MZ %d: correction != value' % t
            v[t] = 0
    if not (ph == ph[0]).all():
        return v, False, 'phase depends on the input'
    return v, True, ''


def bits_of(vals, k):
    """python ints -> (k, B) uint8, LSB first."""
    vals = [int(x) for x in vals]
    out = np.zeros((k, len(vals)), np.uint8)
    for j in range(k):
        out[j] = np.fromiter(((x >> j) & 1 for x in vals), np.uint8, len(vals))
    return out


def ints_of(v, wires):
    """rows `wires` of a bit matrix -> list of python ints."""
    B = v.shape[1]; out = [0] * B
    for j, w in enumerate(wires):
        row = v[w]
        for i in np.nonzero(row)[0]:
            out[i] |= 1 << j
    return out
