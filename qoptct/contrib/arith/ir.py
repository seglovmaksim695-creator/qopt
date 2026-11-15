# SPDX-License-Identifier: Apache-2.0
"""Circuit builder on logical ANDs with exact bit-sliced checking, lowering to OpenQASM 3, and an independent
sparse simulator for the emitted text.

IR operations (wire indices: 0..n-1 are the wires of the .qc file, n.. are new qubits):
  ('x', (t,)), ('cx', (c, t)), ('cz', (a, b)), ('z'|'s'|'sdg'|'t'|'tdg', (q,))
  ('and', (a, b, t, pol))  logical AND into a |0> wire, pol in 'ab', 'a!b', '!ab', 'or'   (4 T, 3 CX)
  ('ccx', (a, b, t))       unitary Toffoli on an arbitrary target                         (7 T, 6 CX)
  ('erase', (t, corr))     X-basis measurement of t, then if outcome 1: corr gates and X(t) (0 T)
Every operation is applied to the bit-sliced values of all wires over ALL inputs of the truth table, so an
``erase`` is accepted only if its correction cancels the measurement phase exactly on every input
(4*v_t + phase(corr) = const mod 8, in units of pi/4).
"""
import cmath
import math
import random
import re

import numpy as np

PH_UNITS = {'z': 4, 's': 2, 'sdg': 6, 't': 1, 'tdg': 7}
_UNIT_GATES = {0: [], 1: ['t'], 2: ['s'], 3: ['s', 't'], 4: ['z'], 5: ['z', 't'], 6: ['sdg'], 7: ['tdg']}
# phase programs of the 4-T logical AND (H; ph; CX a; ph; CX b; ph; CX a; ph; H; S) for each polarity
_AND_PH = {'ab': (3, 5, 3, 5), '!ab': (3, 3, 5, 5), 'a!b': (1, 7, 7, 1), 'or': (3, 3, 3, 7)}


class BuildError(Exception):
    pass


class Builder:
    def __init__(self, table, reserved=()):
        """table: tt.Table. reserved: clean wires of the .qc that carry outputs (not usable as scratch)."""
        self.T = table
        self.n = table.n
        z = table.zeros()
        self.val = {w: (table.inp[w].copy() if w in table.inp else z.copy()) for w in range(self.n)}
        self.free = [w for w in table.clean if w not in set(reserved)]
        self.N = self.n
        self.ops = []
        self.gph = [z.copy(), z.copy(), z.copy()]   # global phase polynomial (units of pi/4), bit planes
        self.nand = 0

    # ------------------------------------------------------------ wires
    def alloc(self):
        if self.free:
            w = self.free.pop(0)
        else:
            w = self.N
            self.N += 1
            self.val[w] = self.T.zeros()
        if not self.T.is_zero(self.val[w]):
            raise BuildError('alloc of a dirty wire')
        return w

    def release(self, w):
        if not self.T.is_zero(self.val[w]):
            raise BuildError('release of a dirty wire %d' % w)
        if w >= self.n or w in self.T.clean:
            self.free.insert(0, w)

    def extra(self):
        return self.N - self.n

    # ------------------------------------------------------------ gates
    def x(self, t):
        self.ops.append(('x', (t,)))
        self.val[t] = self.T.mask(~self.val[t])

    def cx(self, c, t):
        if c == t:
            raise BuildError('cx on one wire')
        self.ops.append(('cx', (c, t)))
        self.val[t] = self.val[t] ^ self.val[c]

    def andv(self, a, b, pol):
        A, B = self.val[a], self.val[b]
        if pol == 'ab':
            return A & B
        if pol == 'a!b':
            return self.T.mask(A & ~B)
        if pol == '!ab':
            return self.T.mask(~A & B)
        if pol == 'or':
            return A | B
        raise BuildError(pol)

    def AND(self, a, b, t=None, pol='ab'):
        if t is None:
            t = self.alloc()
        if len({a, b, t}) != 3:
            raise BuildError('AND wires')
        if not self.T.is_zero(self.val[t]):
            raise BuildError('AND target %d not clean' % t)
        self.ops.append(('and', (a, b, t, pol)))
        self.val[t] = self.andv(a, b, pol)
        self.nand += 1
        return t

    def ccx(self, a, b, t):
        if len({a, b, t}) != 3:
            raise BuildError('ccx wires')
        self.ops.append(('ccx', (a, b, t)))
        self.val[t] = self.val[t] ^ (self.val[a] & self.val[b])

    def rtof(self, a, b, t, dagger=False):
        """relative-phase Toffoli CS^dagger(a,b)·Tof (4 T, 4 CX); dagger=True: its inverse. A pair rtof ... rtof^dagger
        with unchanged controls is an exact Toffoli pair (the phases are tracked in the global phase polynomial)."""
        if len({a, b, t}) != 3:
            raise BuildError('rtof wires')
        self.ops.append(('rtofdg' if dagger else 'rtof', (a, b, t)))
        self._addph(self.gph, self.val[a] & self.val[b], 2 if dagger else 6)
        self.val[t] = self.val[t] ^ (self.val[a] & self.val[b])

    def phase(self, name, *q):
        """diagonal phase gate outside a measurement (tracked in the global phase polynomial)."""
        self.ops.append((name, tuple(q)))
        self._addph(self.gph, self._gate_bits(name, q), PH_UNITS.get(name, 4))

    def _gate_bits(self, name, q):
        if name == 'cz':
            return self.val[q[0]] & self.val[q[1]]
        return self.val[q[0]]

    def _addph(self, P, v, u):
        a = [v if (u >> i) & 1 else None for i in range(3)]
        c = None
        for i in range(3):
            ai = a[i]
            s = P[i]
            if ai is not None and c is not None:
                ns = s ^ ai ^ c
                c = (s & ai) | (s & c) | (ai & c)
            elif ai is not None:
                ns = s ^ ai
                c = s & ai
            elif c is not None:
                ns = s ^ c
                c = s & c
            else:
                ns = s
            P[i] = ns

    def _const(self, P):
        for p in P:
            if not (self.T.is_zero(p) or self.T.eq(p, self.T.ones())):
                return False
        return True

    def erase(self, t, corr):
        """corr: list of (gate, qubits) with gates cz, z, s, sdg, t, tdg (not acting on t)."""
        z = self.T.zeros()
        P = [z.copy(), z.copy(), z.copy()]
        self._addph(P, self.val[t], 4)
        for g, q in corr:
            if t in q:
                raise BuildError('correction acts on the measured wire')
            self._addph(P, self._gate_bits(g, q), PH_UNITS.get(g, 4))
        if not self._const(P):
            raise BuildError('erase of wire %d: correction does not cancel the phase' % t)
        self.ops.append(('erase', (t, tuple(corr))))
        self.val[t] = self.T.zeros()
        self.release(t)

    def erase_and(self, t, a, b, pol='ab', lin=()):
        """erase t holding pol(a, b) ^ XOR(lin) (current values)."""
        corr = list(and_phase(a, b, pol)) + [('z', (w,)) for w in lin]
        self.erase(t, corr)

    # ------------------------------------------------------------ checks and output
    def finish(self, extra_ok=True):
        """exact check against the truth table: outputs of the .qc wires, clean wires and new qubits back to |0>,
        constant global phase."""
        T = self.T
        for w in range(self.n):
            if not T.eq(self.val[w], T.out[w]):
                raise BuildError('wire %d (%s) differs from the reference' % (w, T.names[w]))
        for w in range(self.n, self.N):
            if not T.is_zero(self.val[w]):
                raise BuildError('new qubit %d not returned to |0>' % w)
        if not self._const(self.gph):
            raise BuildError('data-dependent phase left')
        return True

    def stats(self):
        T = CX = 0
        for g, q in self.ops:
            if g == 'and':
                T += 4
                CX += 3
            elif g == 'ccx':
                T += 7
                CX += 6
            elif g in ('rtof', 'rtofdg'):
                T += 4
                CX += 4
            elif g in ('cx', 'cz'):
                CX += 1
            elif g in ('t', 'tdg'):
                T += 1
            elif g == 'erase':
                for gg, qq in q[1]:
                    if gg in ('cx', 'cz'):
                        CX += 1
                    elif gg in ('t', 'tdg'):
                        T += 1
        return T, CX


def and_phase(a, b, pol):
    """correction gates implementing (-1)^{pol(a,b)} on live wires a, b."""
    if pol == 'ab':
        return [('cz', (a, b))]
    if pol == 'a!b':
        return [('cz', (a, b)), ('z', (a,))]
    if pol == '!ab':
        return [('cz', (a, b)), ('z', (b,))]
    if pol == 'or':
        return [('cz', (a, b)), ('z', (a,)), ('z', (b,))]
    raise BuildError(pol)


# ---------------------------------------------------------------- lowering
def _ph(lines, u, w):
    for g in _UNIT_GATES[u % 8]:
        lines.append('%s q[%d];' % (g, w))


def lower_ops(ops):
    """IR -> (list of qasm statements, number of classical bits)."""
    L = []
    nb = 0
    for g, q in ops:
        if g == 'and':
            a, b, t, pol = q
            k = _AND_PH[pol]
            L.append('h q[%d];' % t)
            _ph(L, k[0], t)
            L.append('cx q[%d], q[%d];' % (a, t))
            _ph(L, k[1], t)
            L.append('cx q[%d], q[%d];' % (b, t))
            _ph(L, k[2], t)
            L.append('cx q[%d], q[%d];' % (a, t))
            _ph(L, k[3], t)
            L.append('h q[%d];' % t)
            L.append('s q[%d];' % t)
        elif g == 'ccx':
            a, b, t = q
            L += ['h q[%d];' % t, 'cx q[%d], q[%d];' % (b, t), 'tdg q[%d];' % t, 'cx q[%d], q[%d];' % (a, t),
                  't q[%d];' % t, 'cx q[%d], q[%d];' % (b, t), 'tdg q[%d];' % t, 'cx q[%d], q[%d];' % (a, t),
                  't q[%d];' % b, 't q[%d];' % t, 'h q[%d];' % t, 'cx q[%d], q[%d];' % (a, b), 't q[%d];' % a,
                  'tdg q[%d];' % b, 'cx q[%d], q[%d];' % (a, b)]
        elif g in ('rtof', 'rtofdg'):
            a, b, t = q
            seq = [('h', t), ('t', t), ('cx', b, t), ('tdg', t), ('cx', a, t), ('t', t), ('cx', b, t), ('tdg', t),
                   ('cx', a, t), ('h', t)]
            if g == 'rtofdg':
                seq = [((_INV[s[0]],) + s[1:]) for s in reversed(seq)]
            for s in seq:
                L.append(_gtxt(s[0], s[1:]))
        elif g == 'erase':
            t, corr = q
            L.append('h q[%d];' % t)
            L.append('c[%d] = measure q[%d];' % (nb, t))
            body = ' '.join(_gtxt(gg, qq) for gg, qq in corr) + (' ' if corr else '') + 'x q[%d];' % t
            L.append('if (c[%d]) { %s }' % (nb, body))
            nb += 1
        else:
            L.append(_gtxt(g, q))
    return L, nb


def _gtxt(g, q):
    return '%s %s;' % (g, ', '.join('q[%d]' % i for i in q))


_INV = {'cx': 'cx', 'cz': 'cz', 'h': 'h', 'x': 'x', 'z': 'z', 't': 'tdg', 'tdg': 't', 's': 'sdg', 'sdg': 's'}
_ZT = {'z', 's', 'sdg', 't', 'tdg', 'cz'}


def _types(it):
    if it[0] == 'ME':
        d = {it[1]: 'O'}
        for g, q in it[3]:
            for x in q:
                if x != it[1]:
                    d[x] = 'Z' if (g in _ZT and d.get(x, 'Z') == 'Z') else 'O'
        return d
    g, q = it[1], it[2]
    if g == 'cx':
        return {q[0]: 'Z', q[1]: 'X'}
    if g in _ZT:
        return {x: 'Z' for x in q}
    if g == 'x':
        return {q[0]: 'X'}
    return {x: 'O' for x in q}


def _comm(a, b):
    ta, tb = _types(a), _types(b)
    for x in ta:
        if x in tb and (ta[x] == 'O' or tb[x] == 'O' or ta[x] != tb[x]):
            return False
    return True


def _cancels(a, b):
    if a[0] != 'G' or b[0] != 'G' or _INV.get(a[1]) != b[1]:
        return False
    if a[1] == 'cz':
        return set(a[2]) == set(b[2])
    return a[2] == b[2]


def peephole(lines, maxlook=200):
    """cancel pairs of mutually inverse gates that meet after commuting through gates of the same Pauli type on
    shared qubits (Z-type: diagonal gates and CX controls; X-type: X and CX targets); t*t -> s.
    A measurement with its conditional correction is one barrier element."""
    prog = parse(lines)
    items = []
    i = 0
    while i < len(prog):
        op = prog[i]
        if op[0] == 'M':
            nxt = prog[i + 1] if i + 1 < len(prog) else None
            corr = nxt[2] if nxt is not None and nxt[0] == 'IF' and nxt[1] == op[2] else []
            items.append(('ME', op[1], op[2], corr))
            i += 2 if corr else 1
            continue
        if op[0] == 'IF':
            return list(lines)   # unexpected layout: leave untouched
        items.append(op)
        i += 1
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(items):
            a = items[i]
            j = i + 1
            done = False
            while j < len(items) and j - i < maxlook:
                b = items[j]
                if _cancels(a, b):
                    del items[j]
                    del items[i]
                    changed = done = True
                    break
                if a[0] == 'G' and b[0] == 'G' and a[1] == b[1] and a[1] in ('t', 'tdg') and a[2] == b[2]:
                    del items[j]
                    items[i] = ('G', 's' if a[1] == 't' else 'sdg', a[2])
                    changed = done = True
                    break
                if not _comm(a, b):
                    break
                j += 1
            if not done:
                i += 1
    out = []
    for it in items:
        if it[0] == 'ME':
            out.append('c[%d] = measure q[%d];' % (it[2], it[1]))
            out.append('if (c[%d]) { %s }' % (it[2], ' '.join(_gtxt(g, q) for g, q in it[3])))
        else:
            out.append(_gtxt(it[1], it[2]))
    return out


def to_qasm(body, N, nbits, plus=()):
    head = ['OPENQASM 3;', 'include "stdgates.inc";', 'qubit[%d] q;' % N]
    if nbits:
        head.append('bit[%d] c;' % nbits)
    pre = ['h q[%d];' % w for w in plus]
    return '\n'.join(head + pre + list(body)) + '\n'


def count(txt):
    T = len(re.findall(r'\b(?:t|tdg) q', txt))
    CX = len(re.findall(r'\b(?:cx|cz) q', txt))
    return T, CX


# ---------------------------------------------------------------- independent sparse simulation of the text
_W8 = cmath.exp(1j * math.pi / 4)
_PHV = {'z': -1, 's': 1j, 'sdg': -1j, 't': _W8, 'tdg': _W8.conjugate()}


def parse(lines):
    prog = []
    for s in lines:
        s = s.strip()
        if not s or s.startswith(('OPENQASM', 'include', 'qubit', 'bit')):
            continue
        m = re.match(r'c\[(\d+)\] = measure q\[(\d+)\];$', s)
        if m:
            prog.append(('M', int(m.group(2)), int(m.group(1))))
            continue
        m = re.match(r'if \(c\[(\d+)\]\) \{(.*)\}$', s)
        if m:
            body = [g for g in m.group(2).split(';') if g.strip()]
            prog.append(('IF', int(m.group(1)), [_pg(g) for g in body]))
            continue
        prog.append(('G',) + _pg(s))
    return prog


def _pg(g):
    g = g.strip().rstrip(';')
    name, args = g.split(None, 1)
    return name, tuple(int(a) for a in re.findall(r'q\[(\d+)\]', args))


def _apply(st, name, q):
    if name == 'h':
        b = 1 << q[0]
        r = 1 / math.sqrt(2)
        out = {}
        for k, a in st.items():
            k0, k1 = k & ~b, k | b
            out[k0] = out.get(k0, 0) + a * r
            out[k1] = out.get(k1, 0) + (-a * r if k & b else a * r)
        return {k: v for k, v in out.items() if abs(v) > 1e-12}
    out = {}
    for k, a in st.items():
        if name == 'x':
            k ^= 1 << q[0]
        elif name == 'cx':
            if k >> q[0] & 1:
                k ^= 1 << q[1]
        elif name == 'ccx':
            if k >> q[0] & 1 and k >> q[1] & 1:
                k ^= 1 << q[2]
        elif name == 'cz':
            if k >> q[0] & 1 and k >> q[1] & 1:
                a = -a
        elif name in _PHV:
            if k >> q[0] & 1:
                a = a * _PHV[name]
        else:
            raise ValueError(name)
        out[k] = out.get(k, 0) + a
    return out


def run_forced(prog, init, outcomes):
    st = {init: 1 + 0j}
    bits = {}
    for op in prog:
        if op[0] == 'G':
            st = _apply(st, op[1], op[2])
        elif op[0] == 'M':
            q, k = op[1], op[2]
            m = outcomes[k]
            st = {kk: a for kk, a in st.items() if bool(kk >> q & 1) == bool(m)}
            if not st:
                return None
            bits[k] = m
        else:
            if bits.get(op[1]):
                for name, q in op[2]:
                    st = _apply(st, name, q)
    return st


def check_text(body, nbits, table, samples=1500, branches=3, seed=1):
    """Independent check of the emitted statements (without the leading H of plus wires) against the truth table:
    for forced random measurement outcomes, every sampled basis input must give exactly one basis state equal to
    the reference output (all other wires |0>) with the same amplitude for all inputs. Returns (ok, message)."""
    T = table
    prog = parse(body)
    rng = random.Random(seed)
    if T.k <= 12:
        xs = list(range(1 << T.k))
    else:
        xs = [0, (1 << T.k) - 1] + [1 << j for j in range(T.k)] + [rng.getrandbits(T.k) for _ in range(samples)]
    nb1 = max(nbits, 1)
    patterns = [[0] * nb1, [1] * nb1] + [[rng.getrandbits(1) for _ in range(nb1)] for _ in range(branches)]
    for outc in patterns:
        ref = None
        for x in xs:
            init = 0
            for j, w in enumerate(T.vars):
                if x >> j & 1:
                    init |= 1 << w
            st = run_forced(prog, init, outc)
            if st is None:
                return False, 'impossible branch'
            nz = [(k, a) for k, a in st.items() if abs(a) > 1e-9]
            if len(nz) != 1:
                return False, 'input %d: %d basis states' % (x, len(nz))
            k, a = nz[0]
            want = 0
            for w in range(T.n):
                if T.at(w, x):
                    want |= 1 << w
            if k != want:
                return False, 'input %d: wrong output' % x
            if ref is None:
                ref = a
            elif abs(a - ref) > 1e-9:
                return False, 'input %d: data-dependent amplitude' % x
    return True, 'ok (%d inputs x %d branches)' % (len(xs), branches)
