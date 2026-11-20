# SPDX-License-Identifier: Apache-2.0
"""algebra -- word-level (algebraic) exact prover for large integer-arithmetic circuits on logical ANDs.

    from qoptct.verify.algebra import aprove, spec_mulsq, spec_fold
    r = aprove("mul256.qasm", spec)            # -> dict(status, method='algebra', detail, time, stats, ...)
    qoptct aprove CIRCUIT.qasm SPEC.json [--budget SEC] [--json]      (also: qoptct kprove C.qasm SPEC.json --engine algebra)

WHAT IS PROVED (status PROVED).  For EVERY computational-basis input (the declared input wires arbitrary, every other
wire |0>) and EVERY combination of measurement outcomes, the branch of the circuit maps |x, 0> to
c_m |F(x)> with an amplitude c_m = +-2^(-k/2) that does not depend on x, where F is a classical function such that
  * every claim of the spec holds: polynomial identities over the integers ("out == x*y",
    "r + p*(H + thi1 - 1 + g) == L + 2**n*H") and inequalities ("r <= p-1") between the input words and the final
    values of the output words;
  * the wires listed under "keep" end with their input value, the wires listed as garbage end with SOME classical
    function of the input (not specified further), and every other wire ends at |0>.
By linearity this is the statement for all input states.  It is the same normal form that kprove's BDD engine and
exact2 prove, with a word-level specification instead of a .qc reference.

ACCEPTED LANGUAGE (strict; anything else gives status ERROR, never a verdict):
    x q[a];   cx q[a], q[b];
    the logical-AND gadget of Gidney 2018 into a |0> target, exactly one of the two texts AND13 / AND10 below;
    AND erasure:  h q[t];  c[k] = measure q[t];  if (c[k]) {cz q[a], q[b]; ... z q[w]; ... x q[t];}
one statement per line, every classical bit measured once.

METHOD (algebraic rewriting at the word level, in the spirit of the polynomial / Groebner-basis verification of
integer multipliers; all steps are exact, nothing is sampled):
 1. Boolean layer.  Every wire holds an XOR of atoms: input bits, AND atoms (hash-consed on their operand XOR-sets),
    named carries k = AND(U, V) xor C and named sum bits s = X xor Y xor C.  An erasure is accepted only if
    (wire) xor (correction) is identically 0 (or 1): then the outcome-1 branch has no input-dependent phase; an AND is
    accepted only if its target is identically 0.  Most of these checks are syntactic (set equality after expanding
    names); the others use the facts of step 3.
 2. Full-adder cells.  For Boolean u, v, c:  AND(u, v) xor c = MAJ(u xor c, v xor c, c)  and
        2*MAJ(x, y, c) + (x xor y xor c) = x + y + c            (over the integers),
    so every carry chain found in the circuit gives the exact word identity
        sum_i 2^i s_i + 2^m K = sum_i 2^i (X_i + Y_i) + c_0 ,   K = carry out of the top cell.
    When the circuit drops the top carry (addition mod 2^m, two's complement rows) K is eliminated by RANGE analysis:
    with V = E - 2^m (a x_top + b y_top) (signed / unsigned reading of the operands), if the interval of V lies in
    [0, 2^m) the result word equals V (unsigned), if it lies in [-2^(m-1), 2^(m-1)) the two's complement value of the
    result equals V.  If no reading fits, K stays as an unknown bit (the proof then usually ends in UNKNOWN).
 3. Intervals are sound over-approximations computed from the polynomials: (A) words of earlier adders are replaced
    by interval variables (parts of a word: floor(W / 2^d)), then Horner-style interval evaluation of the multilinear
    polynomial; (B) the same after rewriting the polynomial down to the input bits with the adder identities
    (cancellations such as (x0+x1)(y0+y1) - x0 y0 - x1 y1); (C) exhaustive evaluation over a small cut of the circuit
    (the atoms of the cut are treated as free bits); case splits on top bits.  Ranges give bit FACTS (constant top
    bits, equal sign bits), which feed back into step 1 (released "zero" wires).
 4. Claims.  An identity lhs == rhs is proved by rewriting lhs - rhs with the adder identities, newest adder first
    (backward substitution), down to the zero polynomial.  An inequality is proved by the intervals of step 3.

Polynomials are multilinear over Z[1/2] in the atoms (x^2 = x).  All arithmetic is exact (Python integers).

SELF-CHECKS (never the reason for PROVED, only for stopping): the AND gadget texts are verified on the 8-dimensional
state space at each call; an independent bit-sliced simulation runs first on random and corner inputs (a violated
claim with all preconditions held is REFUTED with the input; an input-dependent erasure phase is REFUTED); every adder
identity and every fact used by the proof is re-evaluated on random inputs (a mismatch aborts with status ERROR).

LIMITS (honest).  * Not a general equivalence checker: the circuit must be built from carry chains (ripple adders of
the qoptct.crypto generators: mulsq, fold/andcirc); other structures usually end in UNKNOWN.  * UNKNOWN is not a
refutation.  * Garbage wires are only shown to be classical functions of the input; their values are not specified
unless a claim names them.  * The inequality engine is interval arithmetic: a true bound that needs a finer argument
is not found (UNKNOWN).  * Free input wires are the declared input words only; all other wires are assumed |0>.
* Explicit z / cz outside erasure corrections, H other than in the two patterns, ccx etc. are not accepted.
"""
import ast
import heapq
import json
import os
import random
import re
import sys
import threading
import time

__all__ = ['aprove', 'prove_ops', 'parse_strict', 'spec_mulsq', 'spec_fold', 'main', 'Poly']

EMPTY = frozenset()
ONE = frozenset((0,))
KI, KA, KK, KS, KV = 0, 1, 2, 3, 4          # input bit, AND atom, named carry, named sum bit, virtual carry


class Unknown(Exception):
    """the prover cannot decide (never a verdict)."""


class Internal(Exception):
    """a self-check failed: the prover stops with status ERROR."""


class Refuted(Exception):
    pass


class ParseError(ValueError):
    pass


# ============================================================================ strict parser
AND13 = (('h', 't'), ('t', 't'), ('cx', 'a', 't'), ('cx', 'b', 't'), ('cx', 't', 'a'), ('cx', 't', 'b'), ('tdg', 'a'),
         ('tdg', 'b'), ('t', 't'), ('cx', 't', 'a'), ('cx', 't', 'b'), ('h', 't'), ('s', 't'))
AND10 = (('h', 't'), ('t', 't'), ('cx', 'b', 't'), ('tdg', 't'), ('cx', 'a', 't'), ('t', 't'), ('cx', 'b', 't'),
         ('tdg', 't'), ('h', 't'), ('sdg', 't'))
_TEMPLATES = (AND13, AND10)
_GATE = re.compile(r'([a-z]+) q\[(\d+)\](?:, ?q\[(\d+)\])?;$')
_MEAS = re.compile(r'c\[(\d+)\] = measure q\[(\d+)\];$')
_IF = re.compile(r'if \(c\[(\d+)\]\) \{(.*)\}$')
_gadget_ok = {}


def gadget_check(template):
    """the AND gadget text maps |a, b, 0> to |a, b, ab> with amplitude 1 for all a, b (dense 8-dim state).
    -> max error."""
    import cmath
    import math
    s2 = 1 / math.sqrt(2)
    w = cmath.exp(1j * math.pi / 4)
    ph = {'t': w, 'tdg': w.conjugate(), 's': 1j, 'sdg': -1j}
    pos = {'a': 0, 'b': 1, 't': 2}
    err = 0.0
    for a in (0, 1):
        for b in (0, 1):
            st = [0j] * 8
            st[a | b << 1] = 1
            for g in template:
                out = [0j] * 8
                qs = [pos[r] for r in g[1:]]
                for k, amp in enumerate(st):
                    if amp == 0:
                        continue
                    if g[0] == 'cx':
                        out[k ^ (1 << qs[1]) if k >> qs[0] & 1 else k] += amp
                    elif g[0] == 'h':
                        bit = k >> qs[0] & 1
                        out[k & ~(1 << qs[0])] += amp * s2
                        out[k | (1 << qs[0])] += amp * s2 * (-1 if bit else 1)
                    else:
                        out[k] += amp * (ph[g[0]] if k >> qs[0] & 1 else 1)
                st = out
            want = a | b << 1 | (a & b) << 2
            for k in range(8):
                err = max(err, abs(st[k] - (1 if k == want else 0)))
    return err


def parse_strict(text):
    """-> (N, ops); ops: ('x', q) ('cx', c, t) ('and', a, b, t) ('erase', t, corr), corr = (('cz', a, b) | ('z', a), ...).
    ParseError on anything outside the accepted language."""
    lines = []
    for s in text.split('\n'):
        s = s.split('//')[0].strip()
        if s:
            lines.append(s)
    if len(lines) < 3 or lines[0] not in ('OPENQASM 3;', 'OPENQASM 3.0;') or lines[1] != 'include "stdgates.inc";':
        raise ParseError('header: OPENQASM 3 + stdgates.inc expected')
    mt = re.fullmatch(r'qubit\[(\d+)\] q;', lines[2])
    if not mt:
        raise ParseError('qubit[N] q; expected')
    N = int(mt.group(1))
    i = 3
    nb = 0
    if i < len(lines):
        mb = re.fullmatch(r'bit\[(\d+)\] c;', lines[i])
        if mb:
            nb = int(mb.group(1))
            i += 1
    for tpl in _TEMPLATES:
        if tpl not in _gadget_ok:
            _gadget_ok[tpl] = gadget_check(tpl) < 1e-9
        if not _gadget_ok[tpl]:
            raise ParseError('internal: AND gadget self-check failed')
    used = set()
    ops = []
    L = len(lines)
    gm = _GATE.match

    def bad(msg, k):
        raise ParseError('%s at statement %d: %s' % (msg, k, lines[k][:80] if k < L else '<end>'))

    def qq(k, *qs):
        for q in qs:
            if q >= N:
                bad('qubit index out of range', k)
        if len(set(qs)) != len(qs):
            bad('repeated qubit', k)

    while i < L:
        m = gm(lines[i])
        if not m:
            bad('unsupported statement', i)
        g, a, b = m.group(1), int(m.group(2)), m.group(3)
        if g == 'cx' and b is not None:
            b = int(b)
            qq(i, a, b)
            ops.append(('cx', a, b))
            i += 1
            continue
        if g == 'x' and b is None:
            qq(i, a)
            ops.append(('x', a))
            i += 1
            continue
        if g != 'h' or b is not None:
            bad('unsupported gate', i)
        t = a
        mm = _MEAS.match(lines[i + 1]) if i + 1 < L else None
        if mm:
            k, t1 = int(mm.group(1)), int(mm.group(2))
            mi = _IF.match(lines[i + 2]) if i + 2 < L else None
            if not mi or t1 != t or int(mi.group(1)) != k or k in used or k >= nb:
                bad('malformed erasure', i)
            used.add(k)
            body = [s.strip() for s in mi.group(2).split(';')]
            if body[-1] != '' or len(body) < 2:
                bad('malformed correction', i + 2)
            body = body[:-1]
            corr = []
            for j, s in enumerate(body):
                mg = gm(s + ';')
                if not mg:
                    bad('unsupported correction', i + 2)
                gg, p, q = mg.group(1), int(mg.group(2)), mg.group(3)
                if gg == 'x' and q is None and j == len(body) - 1 and p == t:
                    break
                if gg == 'cz' and q is not None:
                    q = int(q)
                    qq(i + 2, p, q, t)
                    corr.append(('cz', p, q))
                elif gg == 'z' and q is None:
                    qq(i + 2, p, t)
                    corr.append(('z', p))
                else:
                    bad('unsupported correction', i + 2)
            else:
                bad('correction without the final x on the measured wire', i + 2)
            qq(i, t)
            ops.append(('erase', t, tuple(corr)))
            i += 3
            continue
        done = False
        for tpl in _TEMPLATES:
            n = len(tpl)
            if i + n > L:
                continue
            role = {'t': t}
            ok = True
            for j in range(n):
                mj = gm(lines[i + j])
                if not mj or mj.group(1) != tpl[j][0]:
                    ok = False
                    break
                qs = [int(mj.group(2))] + ([int(mj.group(3))] if mj.group(3) is not None else [])
                if len(qs) != len(tpl[j]) - 1:
                    ok = False
                    break
                for r, q in zip(tpl[j][1:], qs):
                    if role.setdefault(r, q) != q:
                        ok = False
                        break
                if not ok:
                    break
            if ok and len(role) == 3:
                qq(i, role['a'], role['b'], role['t'])
                ops.append(('and', role['a'], role['b'], role['t']))
                i += n
                done = True
                break
        if not done:
            bad('unsupported construct (h that is neither an AND gadget nor an erasure)', i)
    return N, ops


# ============================================================================ polynomials
def mmul(m1, m2):
    if not m1:
        return m2
    if not m2:
        return m1
    if len(m1) == 1 and len(m2) == 1:
        a, b = m1[0], m2[0]
        if a < b:
            return (a, b)
        if a > b:
            return (b, a)
        return m1
    return tuple(sorted(set(m1) | set(m2)))


class Poly:
    """multilinear polynomial with coefficients in Z[1/2]: value = sum d[m] * prod(m) / 2^e."""
    __slots__ = ('d', 'e')

    def __init__(self, d=None, e=0):
        self.d = {} if d is None else d
        self.e = e

    @staticmethod
    def const(c):
        return Poly({(): c} if c else {})

    @staticmethod
    def var(v):
        return Poly({(v,): 1})

    def copy(self):
        return Poly(dict(self.d), self.e)

    def rescale(self, k):
        if k:
            self.d = {m: c << k for m, c in self.d.items()}
            self.e += k

    def iadd(self, o, k=1):
        if isinstance(o, int):
            o = Poly.const(o)
        if o.e > self.e:
            self.rescale(o.e - self.e)
        sh = self.e - o.e
        d = self.d
        for m, c in o.d.items():
            v = d.get(m, 0) + ((c * k) << sh)
            if v:
                d[m] = v
            elif m in d:
                del d[m]
        return self

    def __add__(self, o):
        return self.copy().iadd(o)

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -c for m, c in self.d.items()}, self.e)

    def __sub__(self, o):
        return self.copy().iadd(o, -1)

    def __rsub__(self, o):
        return (-self).iadd(o)

    def __mul__(self, o):
        if isinstance(o, int):
            return Poly({m: c * o for m, c in self.d.items()} if o else {}, self.e)
        r = {}
        for m1, c1 in self.d.items():
            for m2, c2 in o.d.items():
                m = mmul(m1, m2)
                v = r.get(m, 0) + c1 * c2
                if v:
                    r[m] = v
                elif m in r:
                    del r[m]
        return Poly(r, self.e + o.e)

    __rmul__ = __mul__

    def __pow__(self, k):
        if not isinstance(k, int) or k < 0 or k > 8:
            raise ValueError('power must be a small non-negative integer')
        r = Poly.const(1)
        for _ in range(k):
            r = r * self
        return r

    def vars(self):
        s = set()
        for m in self.d:
            s.update(m)
        return s


def _ib(D, ranges, rev):
    """interval of a multilinear polynomial (integer dict D); variables in [0, 1] unless listed in ranges.
    Horner-style: P = c0 + sum_v v * Q_v with Q_v free of the variables before v in the order."""
    lo = hi = 0
    groups = {}
    for m, c in D.items():
        if not m:
            lo += c
            hi += c
            continue
        if rev:
            v, r = m[-1], m[:-1]
        else:
            v, r = m[0], m[1:]
        g = groups.get(v)
        if g is None:
            groups[v] = g = {}
        g[r] = c
    for v, g in groups.items():
        if len(g) == 1 and () in g:
            ql = qh = g[()]
        else:
            ql, qh = _ib(g, ranges, rev)
        rg = ranges.get(v) if ranges else None
        if rg is None:
            if ql < 0:
                lo += ql
            if qh > 0:
                hi += qh
        else:
            a, b = rg
            p = (a * ql, a * qh, b * ql, b * qh)
            lo += min(p)
            hi += max(p)
    return lo, hi


def _meet(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return (max(a[0], b[0]), min(a[1], b[1]))


def _sbits(lo, hi):
    k = 1
    while lo < -(1 << (k - 1)) or hi > (1 << (k - 1)) - 1:
        k += 1
    return k


def _pow2(w):
    w = abs(w)
    return w.bit_length() - 1 if w and not (w & (w - 1)) else None


# ============================================================================ engine
class Cell:
    """full-adder cell: cout = MAJ(X, Y, C), sum = X xor Y xor C (XOR-sets of atoms)."""
    __slots__ = ('X', 'Y', 'C', 'cout', 'sumkey', 'adder', 'opi', 'prevc')

    def __init__(self, X, Y, C, cout, sumkey, opi):
        self.X, self.Y, self.C, self.cout, self.sumkey, self.opi = X, Y, C, cout, sumkey, opi
        self.adder = None
        self.prevc = None          # raw AND cell not yet used as a carry: the carry atom made just before it


class Adder:
    __slots__ = ('cells', 'order', 'E', 'bits', 'pos', 'V', 'rng', 'signed', 'stage', 'virt', 'lhs_cache', 'kindtxt',
                 'cands', 'deep', 'alt', 'vanish')

    def __init__(self, cells):
        self.cells = cells
        self.order = cells[-1].opi
        self.E = None
        self.stage = -1
        self.virt = None
        self.lhs_cache = None
        self.deep = False
        self.vanish = ()


class Engine:
    EXPAND_LIMIT = 48

    def __init__(self, N, in_wires, lanes=48, seed=7, deadline=None):
        self.N = N
        self.L = lanes
        self.mask = (1 << lanes) - 1
        self.deadline = deadline
        rng = random.Random(seed)
        self.kind = [KI]
        self.adef = [None]
        self.sim = [self.mask]
        self._fs = {0: ONE}
        self.wires = [EMPTY] * N
        self.in_atom = {}
        for w in in_wires:
            v = rng.getrandbits(lanes) & ~3 | 2          # lane 0: all zeros, lane 1: all ones
            a = self.new_atom(KI, None, v)
            self.in_atom[w] = a
            self.wires[w] = self.fs(a)
        self.and_tab = {}
        self.carry_tab = {}
        self.atom_cell = {}        # carry atom (KA / KK) -> the cell that produces it
        self.cin_cell = {}         # carry atom -> the cell that consumes it as carry-in
        self.sumkey_cell = {}      # XOR-set -> cell whose sum bit it is
        self.sum_tab = {}          # XOR-set -> named sum atom
        self.sig_cell = {}         # named sum atom -> cell
        self.virt_adder = {}
        self.facts = {}
        self.fver = 0
        self.fresh = {}
        self.wev = {}
        self.opi = 0
        self.last_gate = 0
        self.last_carry = None     # (atom, opi)
        self.pend = None
        self.aad = {}
        self.aprod = {}
        self.phantom = set()
        self.in_zero = False
        self.atom_wire = {}
        self.refining = False
        self.vint = {}             # virtual carry symbol -> its integer value on every simulation lane
        self.vrange = {}
        self.stats = dict(ands=0, erases=0, cells=0, adders=0, facts=0, stageB=0, stageC=0, splits=0, sem_zero=0,
                          virt=0, reduce_steps=0)
        self.cls = {}              # classification histogram

    # ------------------------------------------------------------------ atoms and sets
    def new_atom(self, kind, d, sim):
        self.kind.append(kind)
        self.adef.append(d)
        self.sim.append(sim)
        return len(self.kind) - 1

    def fs(self, a):
        s = self._fs.get(a)
        if s is None:
            s = self._fs[a] = frozenset((a,))
        return s

    def simset(self, S):
        x = 0
        sim = self.sim
        for a in S:
            x ^= sim[a]
        return x

    def and_set(self, U, V, prev=None):
        """XOR-set of AND(U, V)."""
        if not U or not V:
            return EMPTY
        if U == V:
            return U
        if U == ONE:
            return V
        if V == ONE:
            return U
        if U ^ V == ONE:
            return EMPTY
        key = frozenset((U, V))
        a = self.and_tab.get(key)
        if a is None:
            a = self.new_atom(KA, (U, V), self.simset(U) & self.simset(V))
            self.and_tab[key] = a
            cell = Cell(U, V, EMPTY, a, U ^ V, self.opi)
            cell.prevc = prev if prev is not None else False
            self.stats['cells'] += 1
            self.atom_cell[a] = cell
        return self.fs(a)

    def activate(self, t):
        """the raw AND atom t is used as a carry: its half-adder cell joins a chain.  MAJ(x, 0, c) = AND(x, c) with c
        the carry made just before it is read as the cell (x, 0, c) of the chain of c."""
        while True:
            cell = self.atom_cell[t]
            pc = cell.prevc
            if pc is None:
                return
            cell.prevc = None
            nxt = None
            if pc is not False and pc not in self.cin_cell:
                U, V = self.adef[t]
                ps = self.fs(pc)
                if V == ps or U == ps:
                    cell.X, cell.Y, cell.C = (U if V == ps else V), EMPTY, ps
                    self.cin_cell[pc] = cell
                    if self.kind[pc] == KA:
                        nxt = pc
            if cell.sumkey not in self.sumkey_cell:
                self.sumkey_cell[cell.sumkey] = cell
            if nxt is None:
                return
            t = nxt

    def name_carry(self, t, C):
        key = (t, C)
        k = self.carry_tab.get(key)
        if k is None:
            U, V = self.adef[t]
            k = self.new_atom(KK, key, self.sim[t] ^ self.simset(C))
            self.carry_tab[key] = k
            cell = Cell(U ^ C, V ^ C, C, k, U ^ V ^ C, self.opi)
            self.stats['cells'] += 1
            self.atom_cell[k] = cell
            if cell.sumkey not in self.sumkey_cell:
                self.sumkey_cell[cell.sumkey] = cell
            if len(C) == 1:
                (c,) = C
                if self.kind[c] in (KA, KK) and c not in self.cin_cell:
                    self.cin_cell[c] = cell
                    if self.kind[c] == KA:
                        self.activate(c)
        return k

    def name_sum(self, S, cell):
        a = self.new_atom(KS, S, self.simset(S))
        self.sum_tab[S] = a
        self.sig_cell[a] = cell
        return a

    def canon(self, S):
        if len(S) < 2:
            return S
        a = self.sum_tab.get(S)
        if a is not None:
            if self.phantom:
                self.phantom.discard(a)
            return self.fs(a)
        cell = self.sumkey_cell.get(S)
        if cell is not None:
            if len(S) > 2 or 0 not in S:
                return self.fs(self.name_sum(S, cell))
            return S
        if 0 in S and len(S) > 2:
            S2 = S ^ ONE
            a = self.sum_tab.get(S2)
            if a is None:
                cell = self.sumkey_cell.get(S2)
                if cell is not None and (len(S2) > 2 or 0 not in S2):
                    a = self.name_sum(S2, cell)
            if a is not None:
                return frozenset((0, a))
        return S

    def apply_facts(self, S):
        f = self.facts
        if not f:
            return S
        for _ in range(64):
            hit = [a for a in S if a in f]
            if not hit:
                return S
            for a in hit:
                S = S ^ self.fs(a) ^ f[a]
        raise Internal('fact substitution does not terminate')

    def read(self, w):
        S = self.wires[w]
        if len(S) > 1 or (self.facts and S and not S.isdisjoint(self.facts)):
            S2 = self.canon(self.apply_facts(S))
            if S2 is not S:
                self.wires[w] = S = S2
        if S and S is not ONE:
            x = 0
            sim = self.sim
            kind = self.kind
            for a in S:
                x ^= sim[a]
                if kind[a] == KS and not self.in_zero:
                    cell = self.sig_cell[a]
                    if cell.adder is None:
                        # first use of a sum bit of a finished adder: classify the adder now (range facts such as
                        # constant top bits must be known before the wire is reused)
                        self.in_zero = True
                        try:
                            J = self.assemble(cell, True)
                            if J is not None:
                                self.classify(J)
                        finally:
                            self.in_zero = False
                        if self.facts and not S.isdisjoint(self.facts):
                            return self.read(w)
            if (x == 0 or x == self.mask) and S != ONE and not self.in_zero and \
                    not any(kind[a] == KA or kind[a] == KK for a in S):     # (not a transient X xor carry)
                # 0 (or 1) on every simulated input: try to prove it, so that a reused wire is a clean |0>
                self.in_zero = True
                try:
                    if x == 0:
                        if self.is_zero(S, 0, False):
                            self.wires[w] = S = EMPTY
                    elif self.is_zero(S ^ ONE, 0, False):
                        self.wires[w] = S = ONE
                finally:
                    self.in_zero = False
        return S

    # ------------------------------------------------------------------ gates
    def op_cx(self, c, t):
        if c is None:
            S = ONE
        else:
            S = self.read(c)
            self.unfresh(c)
        D = self.wires[t] if t in self.fresh else self.read(t)     # a fresh AND output stays an AND atom
        self.opi += 1
        if not S:
            return
        if len(D) == 1:
            (a,) = D
            k = self.kind[a]
            if k == KA and self.fresh.get(t) == a:
                kk = self.name_carry(a, S)
                self.wires[t] = self.fs(kk)
                self.fresh.pop(t, None)
                self.last_carry = (kk, self.opi)
                self.atom_wire[kk] = t
                return
            if k == KK and self.adef[a][1] == S:
                self.wires[t] = self.fs(self.adef[a][0])
                return
        self.fresh.pop(t, None)
        new = D ^ S
        if len(D) == 1 and k == KS and S != ONE:
            cl = self.sig_cell[a]
            if S == cl.C or S == cl.Y or S == cl.X:
                new = self.adef[a] ^ S      # undo a sum bit: (X xor Y xor C) xor C = X xor Y in its original form
        lc = self.last_carry
        p = self.pend
        if lc is not None and self.last_gate <= lc[1] and len(S) == 1 and lc[0] in S:
            X, Y = D, EMPTY
            ev = self.wev.get(t)
            if ev is not None and ev[2] > lc[1] and ev[0] ^ ev[1] == D:
                X, Y = ev[0], ev[1]
            self.pend = [lc[0], t, X, Y, Y is not EMPTY]
        elif p is not None and p[1] == t:
            if p[4]:
                self.pend = None
            else:
                p[3] = S
                p[4] = True
        self.wev[t] = (D, S, self.opi)
        self.wires[t] = self.canon(new)

    def unfresh(self, w):
        """the wire w is read: an AND output that was not turned into AND xor C is a carry / data bit as it is."""
        f = self.fresh.pop(w, None)
        if f is not None:
            self.activate(f)
            S = self.wires[w]
            if len(S) == 1 and f in S:
                pass

    def finalize_pend(self):
        p = self.pend
        if p is None:
            return
        self.pend = None
        k, w, X, Y = p[0], p[1], p[2], p[3]
        if k in self.cin_cell or k in X or k in Y:
            return
        S = X ^ Y ^ self.fs(k)
        if self.wires[w] != S or S in self.sumkey_cell or S in self.sum_tab or len(S) < 2:
            return
        cell = Cell(X, Y, self.fs(k), None, S, self.opi)
        self.stats['cells'] += 1
        self.cin_cell[k] = cell
        self.sumkey_cell[S] = cell
        if self.kind[k] == KA:
            self.activate(k)
        self.wires[w] = self.fs(self.name_sum(S, cell))

    def op_and(self, a, b, t):
        U = self.read(a)
        V = self.read(b)
        self.unfresh(a)
        self.unfresh(b)
        p = self.pend
        if p is not None:
            if p[0] in U and p[0] in V:
                self.pend = None
            else:
                self.finalize_pend()
        T = self.read(t)
        if T and not self.is_zero(T):
            raise Unknown('op %d: AND target q[%d] is not proved to be |0>' % (self.opn, t))
        lc = self.last_carry
        prev = lc[0] if (lc is not None and self.last_gate <= lc[1]) else None
        self.opi += 1
        self.last_gate = self.opi
        self.stats['ands'] += 1
        R = self.and_set(U, V, prev)
        self.wires[t] = R
        self.fresh.pop(t, None)
        if len(R) == 1:
            (r,) = R
            if self.kind[r] == KA and self.adef[r] in ((U, V), (V, U)):
                self.fresh[t] = r
                self.last_carry = (r, self.opi)
                self.atom_wire[r] = t

    def op_erase(self, t, corr):
        self.finalize_pend()
        D = self.read(t)
        self.fresh.pop(t, None)
        for c in corr:
            if c[0] == 'cz':
                self.unfresh(c[1])
                self.unfresh(c[2])
                D = D ^ self.and_set(self.read(c[1]), self.read(c[2]))
            else:
                self.unfresh(c[1])
                D = D ^ self.read(c[1])
        self.opi += 1
        self.last_gate = self.opi
        self.stats['erases'] += 1
        if D and not self.is_zero(D) and not self.is_zero(D ^ ONE):
            raise Unknown('op %d: erasure of q[%d]: wire xor correction is not proved constant' % (self.opn, t))
        self.wires[t] = EMPTY

    def run(self, ops):
        self.opn = 0
        cx, x_, and_, er = self.op_cx, None, self.op_and, self.op_erase
        dl = self.deadline
        for i, o in enumerate(ops):
            self.opn = i
            k = o[0]
            if k == 'cx':
                cx(o[1], o[2])
            elif k == 'x':
                cx(None, o[1])
            elif k == 'and':
                and_(o[1], o[2], o[3])
            else:
                er(o[1], o[2])
            if dl is not None and not (i & 4095) and time.time() > dl:
                raise Unknown('time budget exceeded at op %d of %d' % (i, len(ops)))
        self.finalize_pend()

    # ------------------------------------------------------------------ Boolean zero test
    def defset(self, a):
        k = self.kind[a]
        if k == KK:
            t, C = self.adef[a]
            return self.fs(t) ^ C
        if k == KS:
            return self.adef[a]
        return None

    def is_zero(self, D, depth=0, sem=True):
        D = self.apply_facts(D)
        if not D:
            return True
        if self.simset(D):
            return False
        E = D
        kind = self.kind
        for _ in range(self.EXPAND_LIMIT):
            cand = [a for a in E if kind[a] == KK or kind[a] == KS]
            if not cand:
                break
            a = max(cand)
            E = self.apply_facts(E ^ self.fs(a) ^ self.defset(a))
            if not E:
                return True
        if len(E) == 2 and depth < 3:
            a1, a2 = E
            if kind[a1] == KA and kind[a2] == KA:
                U1, V1 = self.adef[a1]
                U2, V2 = self.adef[a2]
                if (self.is_zero(U1 ^ U2, depth + 1) and self.is_zero(V1 ^ V2, depth + 1)) or \
                        (self.is_zero(U1 ^ V2, depth + 1) and self.is_zero(V1 ^ U2, depth + 1)):
                    return True
        if depth:
            return False
        if not sem:
            return self.tt_zero(D)
        self.stats['sem_zero'] += 1
        for deep in (False, True):
            for a in D:
                J = self.var_adder(a)
                if J is not None:
                    self.ensure_facts(J, deep)
            D2 = self.apply_facts(D)
            if D2 != D:
                return self.is_zero(D2)
        return self.tt_zero(D)

    # ------------------------------------------------------------------ integer values
    def val_atom(self, a):
        if a == 0:
            return Poly.const(1)
        if self.kind[a] == KA:
            p = self.aprod.get(a)
            if p is None:
                p = False
                U, V = self.adef[a]
                kd = self.kind
                if len(U) <= 2 and len(V) <= 2 and not any(kd[b] in (KA, KK) for b in U | V):
                    try:
                        pu, pv = self.val(U), self.val(V)
                        if len(pu.d) <= 3 and len(pv.d) <= 3:
                            p = pu * pv
                    except Unknown:
                        p = False
                self.aprod[a] = p
            if p is not False:
                return p
        return Poly.var(a)

    def val(self, S):
        """integer value (0 / 1) of the XOR-set S as a multilinear polynomial in the atoms."""
        S = self.apply_facts(S)
        acc = None
        for a in sorted(S):
            p = self.val_atom(a)
            if acc is None:
                acc = p
            else:
                acc = acc + p - (acc * p) * 2
                if len(acc.d) > 64:
                    raise Unknown('XOR of %d atoms has no compact integer value' % len(S))
        return acc if acc is not None else Poly()

    def is_plain_var(self, a):
        return a != 0 and (self.kind[a] != KA or self.val_atom(a).d == {(a,): 1})

    def peval(self, P, lane):
        """value of P on one simulation lane -> (numerator, e)."""
        sim = self.sim
        vint = self.vint
        s = 0
        for m, c in P.d.items():
            for v in m:
                if v in vint:
                    c *= vint[v][lane]
                elif not (sim[v] >> lane) & 1:
                    c = 0
                    break
            s += c
        return s

    # ------------------------------------------------------------------ adders
    def sum_real(self, cell):
        S = cell.sumkey
        return S in self.sum_tab or len(S - ONE) <= 1

    def assemble(self, cell, strict=False):
        bottom = cell
        while True:
            C = bottom.C
            if len(C) != 1:
                break
            (c,) = C
            prev = self.atom_cell.get(c)
            if prev is None or self.cin_cell.get(c) is not bottom or prev.adder is not None:
                break
            if self.kind[c] == KA and not prev.C and not self.sum_real(prev):
                break
            bottom = prev
        cells = [bottom]
        cur = bottom
        while cur.cout is not None:
            nxt = self.cin_cell.get(cur.cout)
            if nxt is None or nxt.adder is not None:
                break
            if self.kind[cur.cout] == KA and not cur.C and not self.sum_real(cur):
                break
            cells.append(nxt)
            cur = nxt
        if strict:
            aw = self.atom_wire
            for c in cells[:-1]:
                w = aw.get(c.cout)
                if w is not None and self.wires[w] == self.fs(c.cout):
                    return None                     # a carry of this chain is still on its wire: adder in progress
        J = Adder(cells)
        for c in cells:
            c.adder = J
        self.stats['adders'] += 1
        return J

    def atom_adder(self, a):
        k = self.kind[a]
        if k == KS:
            cell = self.sig_cell[a]
            return cell.adder or self.assemble(cell)
        if k == KA or k == KK:
            cell = self.atom_cell.get(a)
            if cell is None or a in self.cin_cell:
                return None
            J = cell.adder or self.assemble(cell)
            return J if J.cells[-1] is cell else None
        if k == KV:
            return self.virt_adder[a]
        return None

    def var_adder(self, v):
        if v.__class__ is not int:
            return None
        r = self.aad.get(v, 0)
        if r == 0:
            r = self.atom_adder(v)
            if r is not None or self.kind[v] in (KI, KV):
                self.aad[v] = r
        return r

    def sum_set(self, cell):
        S = cell.sumkey
        a = self.sum_tab.get(S)
        if a is not None:
            return self.fs(a)
        if len(S - ONE) <= 1:
            return S
        a = self.name_sum(S, cell)                    # phantom sum bit (never on a wire so far)
        self.phantom.add(a)
        return self.fs(a)

    def classify(self, J):
        """exact word identity of the adder J: sum of weight * bit = V, plus the interval of V."""
        if J.stage >= 0:
            return
        first = True
        if first:
            cells = J.cells
            m = len(cells)
            E = Poly()
            for i, c in enumerate(cells):
                E.iadd(self.val(c.X), 1 << i)
                E.iadd(self.val(c.Y), 1 << i)
            E.iadd(self.val(cells[0].C))
            J.E = E
            J.bits = [[1 << i, self.sum_set(c)] for i, c in enumerate(cells)]
            J.pos = {}
            top = cells[-1]
            J.signed = False
            J.cands = None
            if top.cout is not None:
                J.bits.append([1 << m, self.fs(top.cout)])
                J.V = E
                J.kindtxt = 'cout'
                J.rng = _meet(self.bound(E, 0), (0, (2 << m) - 1))
            else:
                xt, yt = self.val(top.X), self.val(top.Y)
                a = self.signedlike(top.X, cells[-2].X if m > 1 else None)
                b = self.signedlike(top.Y, cells[-2].Y if m > 1 else None)
                J.cands = (a, b, E - (xt * a + yt * b) * (1 << m))
                J.alt = [(a2, b2, E - (xt * a2 + yt * b2) * (1 << m)) for a2, b2 in ((0, 0), (1, 0), (0, 1), (1, 1))
                         if (a2, b2) != (a, b) and (a2 == a or any(xt.d)) and (b2 == b or any(yt.d))]
                J.V = None
            for i, (w, S) in enumerate(J.bits):
                if len(S) == 1:
                    (a,) = S
                    J.pos.setdefault(a, i)
                    if self.kind[a] in (KA, KK):
                        self.aad[a] = J
            J.stage = 0
        if J.V is None:
            m = len(J.cells)
            a, b, V = J.cands

            def fits(iv):
                lo, hi = iv
                if lo >= 0 and hi < (1 << m):
                    J.signed = False
                    return True
                if lo >= -(1 << (m - 1)) and hi < (1 << (m - 1)):
                    J.signed = True
                    return True
                return False
            st = 0
            iv = self.bound_runs(V)
            ok = fits(iv)
            if not ok:
                st = 1
                self.stats['stageB'] += 1
                box = [iv, None]

                def hook(Q):
                    box[0] = _meet(box[0], self.bound_runs(Q))
                    if fits(box[0]):
                        box[1] = (V, box[0], 1)
                        return True
                    Qd = self.drop_mod(Q, m)
                    if Qd is not None:
                        iv3 = self.bound_runs(Qd)
                        if fits(iv3):
                            box[1] = (Qd, iv3, 3)
                            return True
                    return False
                Vr = self.reduce(V, hook=hook)
                if box[1] is None:
                    hook(Vr)
                iv = box[0]
                if box[1] is not None:
                    ok = True
                    V, iv, st = box[1]
                if not ok:
                    for a2, b2, V2 in J.alt:        # other readings: only if plainly in range
                        iv2 = self.bound_runs(V2)
                        if fits(iv2):
                            ok, a, b, V, iv, st = True, a2, b2, V2, iv2, 0
                            break
                if not ok:
                    st = 2
                    iv = _meet(iv, self.bound(V, 2))
                    ok = fits(iv)
            J.alt = None
            if ok:
                J.V, J.rng = V, iv
                J.kindtxt = ('s' if J.signed else 'u') + '%d%d/%d' % (a, b, st)
                J.bits[m - 1][0] = -(1 << (m - 1)) if J.signed else (1 << (m - 1))
            else:
                # the dropped top carry is not resolved: W = V - 2^m K with an unknown integer K (kept as a symbol)
                kv = self.new_atom(KV, None, 0)
                vals = []
                for lane in range(self.L):
                    sv = 0
                    for w, S in J.bits:
                        if (self.simset(S) >> lane) & 1:
                            sv += w
                    num = self.peval(V, lane) - (sv << V.e)
                    if num % (1 << (m + V.e)):
                        raise Internal('adder congruence fails on a simulated input')
                    vals.append(num >> (m + V.e))
                self.vint[kv] = vals
                self.vrange[kv] = (-2, 1)
                self.virt_adder[kv] = J
                self.aad[kv] = J
                J.virt = kv
                J.signed = False
                J.bits.append([1 << m, self.fs(kv)])
                J.pos[kv] = m
                J.V = V
                J.kindtxt = 'virtual%d%d' % (a, b)
                J.rng = (0, (1 << m) - 1)
                self.stats['virt'] += 1
        self.cls[J.kindtxt] = self.cls.get(J.kindtxt, 0) + 1
        self.check_identity(J)
        self.derive_facts(J)

    def signedlike(self, S, below):
        """structural guess: the top operand bit S is a sign bit (a copy of the bit below it, or the top bit of a
        two's complement word).  Only selects which reading of the adder is tried; every reading is checked."""
        if not S or S == ONE:
            return 0
        if below is not None and S == below:
            return 1
        T = S - ONE
        if len(T) == 1:
            (a,) = T
            J = self.var_adder(a)
            if J is not None:
                self.classify(J)
                if J.signed and J.bits[-1][1] == self.fs(a):
                    return 1
        return 0

    def drop_mod(self, P, m):
        """P without its monomials that contain a non-input atom and whose coefficient is a multiple of 2^m
        (a multiple of 2^m as a function): congruent to P modulo 2^m.  None if nothing is dropped."""
        mod = (1 << (m + P.e)) - 1
        kind = self.kind
        d = {}
        for mo, c in P.d.items():
            if not (c & mod) and any(kind[v] != KI for v in mo):
                continue
            d[mo] = c
        if len(d) == len(P.d):
            return None
        return Poly(d, P.e)

    def check_identity(self, J):
        """numerical self-check of the word identity on a few simulation lanes."""
        if J.virt is not None:
            return
        for lane in (0, 1, 2, 3, 4, 5):
            s = 0
            for w, S in J.bits:
                if (self.simset(S) >> lane) & 1:
                    s += w
            if (s << J.V.e) != self.peval(J.V, lane):
                raise Internal('adder identity fails on a simulated input (%s)' % J.kindtxt)

    def adder_lhs(self, J):
        """-> (lhs {atom: weight}, rhs Poly): sum lhs = rhs, after applying the known facts."""
        c = J.lhs_cache
        if c is not None and c[0] == self.fver:
            return c[1], c[2]
        lhs = {}
        rhs = self.pfacts(J.V).copy()
        for w, S in J.bits:
            S2 = self.apply_facts(S)
            if len(S2) == 1:
                (a,) = S2
                if a != 0 and self.var_adder(a) is J and self.is_plain_var(a):
                    lhs[a] = lhs.get(a, 0) + w
                    continue
            if S2:
                rhs.iadd(self.val(S2), -w)
        lhs = {a: w for a, w in lhs.items() if w}
        J.lhs_cache = (self.fver, lhs, rhs)
        return lhs, rhs

    def pfacts(self, P):
        """P with every atom that has a fact replaced by the value of its fact."""
        f = self.facts
        hit = [v for v in P.vars() if v in f]
        if not hit:
            return P
        hs = set(hit)
        out = Poly({}, P.e)
        d = out.d
        vals = {v: self.val(f[v]) for v in hit}
        for m, c in P.d.items():
            if hs.isdisjoint(m):
                v = d.get(m, 0) + c
                if v:
                    d[m] = v
                elif m in d:
                    del d[m]
                continue
            t = Poly({tuple(v for v in m if v not in hs): c})
            for v in m:
                if v in hs:
                    t = t * vals[v]
            for m2, c2 in t.d.items():
                v = d.get(m2, 0) + c2
                if v:
                    d[m2] = v
                elif m2 in d:
                    del d[m2]
        return self.pfacts(out)

    def add_fact(self, S, R):
        """the bit S (a single atom) is identically the XOR-set R."""
        if len(S) != 1:
            return False
        (a,) = S
        if a == 0 or a in self.facts or self.kind[a] == KI:
            return False
        R = self.apply_facts(R)
        if a in R:
            return False
        if self.sim[a] != self.simset(R):
            raise Internal('a derived bit fact fails on a simulated input')
        self.facts[a] = R
        self.fver += 1
        self.stats['facts'] += 1
        return True

    def derive_vanish(self, J):
        """word <= hi < 2^k - 1 (unsigned): for every 0 bit j of hi, (the bits where hi has a 1 above j) and bit j
        cannot all be 1.  These vanishing monomials make polynomials in the bits of the word canonical."""
        J.vanish = ()
        lo, hi = J.rng
        if J.virt is not None or J.signed or lo < 0 or hi <= 0:
            return
        bits = J.bits
        n = hi.bit_length()
        if n > len(bits) or hi == (1 << n) - 1:
            return
        atoms = []
        for i in range(n):
            w, S = bits[i]
            S = self.apply_facts(S)
            if w != 1 << i or len(S) != 1 or 0 in S:
                return
            atoms.append(next(iter(S)))
        out = []
        for j in range(n):
            if not (hi >> j) & 1:
                v = frozenset([atoms[j]] + [atoms[i] for i in range(j + 1, n) if (hi >> i) & 1])
                x = self.mask
                for a in v:
                    x &= self.sim[a]
                if x:
                    raise Internal('a derived vanishing monomial fails on a simulated input')
                out.append(v)
        J.vanish = tuple(out)

    def derive_facts(self, J):
        if J.virt is not None:
            return False
        self.derive_vanish(J)
        lo, hi = J.rng
        bits = J.bits
        ch = False
        nbits = len(bits)
        i = nbits - 1
        while i >= 0 and (lo >> i) == (hi >> i):
            ch |= self.add_fact(bits[i][1], ONE if (lo >> i) & 1 else EMPTY)
            i -= 1
        if J.signed and lo < 0 <= hi:
            k = _sbits(lo, hi)
            for j in range(k, nbits):
                ch |= self.add_fact(bits[j][1], bits[k - 1][1])
        return ch

    def ensure_facts(self, J, deep=False, stage=2):
        self.classify(J)
        ch = False
        if deep and not J.deep:
            J.deep = True
            if J.virt is None:
                r = J.rng
                m = len(J.cells)
                if J.cands is None:
                    r = _meet(r, self.bound(J.V, stage))
                else:
                    # W = J.V exactly; V0 = W (mod 2^m).  Every rewriting of V0, and every such rewriting without
                    # its multiples of 2^m that lies in the range of the representation, bounds W.
                    V0 = J.cands[2]
                    exact = not J.kindtxt.endswith('/3')
                    lo_r, hi_r = (-(1 << (m - 1)), (1 << (m - 1)) - 1) if J.signed else (0, (1 << m) - 1)
                    box = [r]

                    def hook(Q):
                        if exact:
                            box[0] = _meet(box[0], self.bound_runs(Q))
                        Qd = self.drop_mod(Q, m)
                        if Qd is not None:
                            iv3 = self.bound_runs(Qd)
                            if iv3[0] >= lo_r and iv3[1] <= hi_r:
                                box[0] = _meet(box[0], iv3)
                        return False
                    hook(self.reduce(V0, hook=hook))
                    r = box[0]
                    if stage >= 2 and exact:
                        r = _meet(r, self.bound(J.V, 2))
                if r[0] > r[1]:
                    raise Internal('empty range of an adder word')
                if r != J.rng:
                    J.rng = r
        ch |= self.derive_facts(J)
        return ch

    # ------------------------------------------------------------------ intervals
    def group(self, D):
        by = {}
        va = self.var_adder
        for m in D:
            for v in m:
                J = va(v)
                if J is not None:
                    s = by.get(J)
                    if s is None:
                        by[J] = s = []
                    s.append(m)
        return by

    def run_range(self, J, run, lhs):
        """interval of (sum of weight * bit over the run) / weight of its lowest bit."""
        pos = J.pos
        d = pos[run[0]]
        runset = set(run)
        w0 = lhs[run[0]]
        blo = sum(min(lhs[a], 0) for a in run)
        bhi = sum(max(lhs[a], 0) for a in run)
        lo, hi = J.rng
        top = True
        cst = 0
        for i in range(d, len(J.bits)):
            w, S = J.bits[i]
            S2 = self.apply_facts(S)
            if not S2:
                continue
            if S2 == ONE:
                cst += w
                continue
            if len(S2) == 1 and next(iter(S2)) in runset:
                continue
            top = False
            break
        if top:
            blo = max(blo, ((lo >> d) << d) - cst)
            bhi = min(bhi, ((hi >> d) << d) - cst)
        elif not J.signed and lo >= 0 and all(lhs[a] == 1 << pos[a] for a in run):
            bhi = min(bhi, min((1 << len(run)) - 1, hi >> d) << d)
        return (blo // w0, bhi // w0)

    def bound_runs(self, P):
        """stage A: words (and parts of words) of adders become interval variables; then interval evaluation."""
        if self.facts:
            P = self.pfacts(P)
        fv = self.fver
        D = dict(P.d)
        e = P.e
        ranges = dict(self.vrange)
        done = {}
        while True:
            by = self.group(D)
            todo = [J for J in by if done.get(J, 0) < 4]
            if not todo:
                break
            changed = False
            for J in sorted(todo, key=lambda j: -j.order):
                done[J] = done.get(J, 0) + 1
                self.classify(J)
                if self.fver != fv:
                    return self.bound_runs(P)       # new facts: start again
                lhs, _ = self.adder_lhs(J)
                if J.virt is not None:
                    lhs = {a: w for a, w in lhs.items() if a != J.virt}
                lin = {}
                seenm = set()
                for m in by[J]:
                    if m in seenm:
                        continue
                    seenm.add(m)
                    for a in m:
                        if a in lhs:
                            d1 = lin.get(a)
                            if d1 is None:
                                lin[a] = d1 = {}
                            d1[tuple(v for v in m if v != a)] = D[m]
                if len(lin) < 2:
                    done[J] = 4
                    continue
                pos = J.pos
                atoms = sorted(lin, key=pos.__getitem__)
                i = 0
                n = len(atoms)
                while i < n:
                    a0 = atoms[i]
                    w0 = lhs[a0]
                    c0 = lin[a0]
                    j = i + 1
                    while j < n and pos[atoms[j]] == pos[atoms[j - 1]] + 1 and _same(lin[atoms[j]], lhs[atoms[j]], c0, w0):
                        j += 1
                    run = atoms[i:j]
                    i = j
                    if len(run) < 2 or w0 <= 0:
                        continue
                    rs = set(run)
                    if any(not rs.isdisjoint(rest) for rest in c0):
                        continue                    # a square of the run: not a product with an outside factor
                    rr = self.run_range(J, run, lhs)
                    om = max(run) + 0.5
                    while om in ranges:
                        om += 0.001
                    ranges[om] = rr
                    for a in run:
                        for rest in lin[a]:
                            del D[tuple(sorted(rest + (a,)))]
                    for rest, c in c0.items():
                        key = tuple(sorted(rest + (om,)))
                        D[key] = D.get(key, 0) + c
                    changed = True
                    break
                if changed:
                    break
                done[J] = 4
        l1, h1 = _ib(D, ranges, False)
        l2, h2 = _ib(D, ranges, True)
        lo = max(l1, l2)
        hi = min(h1, h2)
        return (-((-lo) >> e), hi >> e)

    def bound(self, P, stage=0, assume=None, split=True):
        """sound interval [lo, hi] of the INTEGER-valued polynomial P over all inputs (see the module docstring)."""
        if not P.d:
            return (0, 0)
        iv = self.bound_runs(P)
        if stage >= 1:
            self.stats['stageB'] += 1
            box = [iv]

            def hook(Q):
                box[0] = _meet(box[0], self.bound_runs(Q))
                return False
            P2 = self.reduce(P, assume, hook=hook)
            iv = _meet(box[0], self.bound_runs(P2))
        if stage >= 2:
            r = self.tt_bound(P)
            if r is not None:
                self.stats['stageC'] += 1
                iv = _meet(iv, r)
            if split:
                tops = []
                for a in sorted(v for v in P.vars() if v.__class__ is int):
                    J = self.var_adder(a)
                    if J is None:
                        continue
                    lhs, _ = self.adder_lhs(J)
                    if a in lhs and self.kind[a] != KV and J.pos[a] == max(J.pos[b] for b in lhs):
                        tops.append(a)
                tops = tops[-2:]
                if tops:
                    self.stats['splits'] += 1
                    hull = None
                    for code in range(1 << len(tops)):
                        asg = dict(assume or {})
                        for i, a in enumerate(tops):
                            asg[a] = (code >> i) & 1
                        r = self.bound(self.psubst(P, asg), 2, asg, split=False)
                        if r[0] > r[1]:
                            continue
                        hull = r if hull is None else (min(hull[0], r[0]), max(hull[1], r[1]))
                    if hull is not None:
                        iv = _meet(iv, hull)
        return iv

    @staticmethod
    def psubst(P, asg):
        d = {}
        for m, c in P.d.items():
            if any(v in asg for v in m):
                if any(asg.get(v) == 0 for v in m):
                    continue
                m = tuple(v for v in m if v not in asg)
            v = d.get(m, 0) + c
            if v:
                d[m] = v
            elif m in d:
                del d[m]
        return Poly(d, P.e)

    # ------------------------------------------------------------------ backward rewriting
    def top_pattern(self, Q, monos, lhs, J):
        """P = A * [W^2] + B * W + (terms without the bits of W), W = sum of weight * bit over a top segment R of the
        word of J.  -> (R, A or None, B, monomials to delete) or None."""
        D = Q.d
        lin = {}
        quad = {}
        kill = []
        if J.virt is not None:
            lhs = {a: w for a, w in lhs.items() if a != J.virt}
            if not lhs:
                return None
        flip = False
        extra = None
        for m in monos:
            js = [v for v in m if v in lhs]
            n = len(js)
            if n == 1:
                a = js[0]
                d1 = lin.get(a)
                if d1 is None:
                    lin[a] = d1 = {}
                d1[tuple(v for v in m if v != a)] = D[m]
                kill.append(m)
            elif n == 2:
                a, b = js
                d1 = quad.get((a, b))
                if d1 is None:
                    quad[(a, b)] = d1 = {}
                d1[tuple(v for v in m if v != a and v != b)] = D[m]
                kill.append(m)
            elif n > 2:
                return None
        if not kill:
            return None
        for w in lhs.values():
            if _pow2(w) is None:
                return None
        pos = J.pos
        order = sorted(lhs, key=pos.__getitem__, reverse=True)
        if quad:
            R = set(lin)
            for p in quad:
                R.update(p)
            k = len(R)
            if set(order[:k]) != R:
                return None
            # reference pair: the two highest bits; pairs that deviate from A * 2 w_a w_b stay as they are
            Rl = order[:k]
            van2 = set(v for v in J.vanish if len(v) == 2)
            c0 = None
            for ia in range(k):
                for ib in range(ia + 1, k):
                    c0 = quad.get((Rl[ia], Rl[ib])) or quad.get((Rl[ib], Rl[ia]))
                    if c0 is not None:
                        den0 = 2 * lhs[Rl[ia]] * lhs[Rl[ib]]
                        break
                if c0 is not None:
                    break
            if c0 is None:
                return None
            good = 0
            for ia in range(k):
                for ib in range(ia + 1, k):
                    a, b = Rl[ia], Rl[ib]
                    if van2 and frozenset((a, b)) in van2:
                        good += 1
                        continue
                    q = quad.get((a, b)) or quad.get((b, a)) or {}
                    w = 2 * lhs[a] * lhs[b]
                    if len(q) == len(c0) and all(c0.get(r) is not None and c * den0 == c0[r] * w for r, c in q.items()):
                        good += 1
                        continue
                    dq = Poly(dict(q), Q.e)
                    dq.iadd(_frac(c0, den0, Q.e), -w)
                    if dq.d:
                        if extra is None:
                            extra = Poly()
                        extra.iadd(dq * Poly({mmul((a,), (b,)): 1}))
            if good < k * (k - 1) // 2 - 2 or 2 * good < k * (k - 1) // 2:
                return None
            A = _frac(c0, den0, Q.e)
            B = None
            for a in order[:k]:
                la = Poly(dict(lin.get(a, {})), Q.e)
                la.iadd(A, -lhs[a] * lhs[a])
                if B is None:
                    B = _frac(la.d, lhs[a], la.e)
                else:
                    la.iadd(B, -lhs[a])             # terms of this bit outside the pattern stay as they are
                    if la.d:
                        if extra is None:
                            extra = Poly()
                        extra.iadd(la * Poly.var(a))
        else:
            top = order[0]
            ct = lin.get(top)
            if ct is None:
                return None
            wt = lhs[top]
            if len(order) > 1:
                cn = lin.get(order[1])
                if cn is not None and not _same(cn, lhs[order[1]], ct, wt) and _same(cn, lhs[order[1]], ct, -wt):
                    flip = True             # the word is used with the other sign reading of its top bit
                    wt = -wt
            R = {top}
            for a in order[1:]:
                ca = lin.get(a)
                if ca is None or not _same(ca, lhs[a], ct, wt):
                    break
                R.add(a)
            if len(R) < len(lin):
                kill = [tuple(sorted(rest + (a,))) for a in R for rest in lin[a]]
            A = None
            B = _frac(ct, wt, Q.e)
        ph = self.phantom
        if ph:
            for a in lhs:
                if a not in R and a in ph:
                    return None
        return R, A, B, kill, flip, extra

    def reduce(self, P, assume=None, _depth=0, fallback=False, hook=None):
        """rewrite P with the adder identities, newest adder first (exact; P keeps its value on every input)."""
        Q = self.pfacts(P).copy() if self.facts else P.copy()
        fv = self.fver
        by = {}
        heap = []
        va = self.var_adder

        def reg(m):
            for v in m:
                J = va(v)
                if J is not None:
                    s = by.get(J)
                    if s is None:
                        by[J] = s = set()
                        heapq.heappush(heap, (-J.order, id(J), J))
                    s.add(m)
        for m in Q.d:
            reg(m)
        dl = self.deadline
        nck = 0
        while heap:
            J = heapq.heappop(heap)[2]
            D = Q.d
            monos = [m for m in by[J] if m in D]
            if not monos:
                continue
            self.classify(J)
            if self.fver != fv:                     # new facts: rewrite the polynomial with them
                fv = self.fver
                Q2 = self.pfacts(Q)
                if Q2 is not Q:
                    Q = Q2
                    for sset in by.values():
                        sset.clear()
                    for m in Q.d:
                        reg(m)
                    D = Q.d
                    monos = [m for m in by[J] if m in D]
                    if not monos:
                        continue
            lhs, rhs = self.adder_lhs(J)
            if J.vanish:
                keep = []
                for m in monos:
                    ms = set(m)
                    for v in J.vanish:
                        if v <= ms:
                            del D[m]                # identically 0 on the range of the word
                            break
                    else:
                        keep.append(m)
                monos = keep
                if not monos:
                    continue
            if assume:
                hit = [a for a in lhs if a in assume]
                if hit:
                    lhs = dict(lhs)
                    rhs = rhs.copy()
                    for a in hit:
                        rhs.iadd(Poly.const(assume[a]), -lhs.pop(a))
                rhs = self.psubst(rhs, assume)
            if not lhs:
                continue
            pat = self.top_pattern(Q, monos, lhs, J)
            if pat is None and not J.deep and J.virt is None and not assume and not self.refining:
                # the top bits of the word are missing in the polynomial: maybe they are constants by a sharper range
                core = [a for a in lhs]
                top = max(core, key=J.pos.__getitem__)
                pres = set(v for m in monos for v in m if v in lhs)
                if len(pres) >= 2 and top not in pres:
                    self.refining = True
                    try:
                        self.ensure_facts(J, True, 1)
                    finally:
                        self.refining = False
                    if self.fver != fv:
                        fv = self.fver
                        Q2 = self.pfacts(Q)
                        if Q2 is not Q:
                            Q = Q2
                            for sset in by.values():
                                sset.clear()
                            for m in Q.d:
                                reg(m)
                        D = Q.d
                        monos = [m for m in by[J] if m in D]
                        lhs, rhs = self.adder_lhs(J)
                        if monos and lhs:
                            pat = self.top_pattern(Q, monos, lhs, J)
            if pat is None:
                continue
            R, A, B, kill, flip, extra = pat
            if hook is not None and nck < 16 and len(lhs) - len(R) - (1 if J.virt is not None else 0) >= 2:
                # a word is about to be split (high part := (word - low part) / 2^d): bound the polynomial as it is
                nck += 1
                if hook(Q):
                    return Q
            S = rhs.copy()
            for a, w in lhs.items():
                if a not in R:
                    S.iadd(Poly.var(a), -w)
            if flip:
                top = max(R, key=J.pos.__getitem__)
                S.iadd(Poly.var(top), -2 * lhs[top])
            for m in kill:
                del D[m]
            add = B * S
            if A is not None:
                add.iadd(A * (S * S))
            if extra is not None:
                add.iadd(extra)
            if add.e and add.d:
                tz = min((c & -c).bit_length() - 1 for c in add.d.values())
                k = min(tz, add.e)
                if k:
                    add = Poly({m: c >> k for m, c in add.d.items()}, add.e - k)
            if add.e > Q.e:
                Q.rescale(max(add.e - Q.e, 64))
                D = Q.d
            sh = Q.e - add.e
            for m, c in add.d.items():
                old = D.get(m)
                if old is None:
                    D[m] = c << sh
                    reg(m)
                else:
                    v = old + (c << sh)
                    if v:
                        D[m] = v
                    else:
                        del D[m]
            self.stats['reduce_steps'] += 1
            if dl is not None and time.time() > dl:
                raise Unknown('time budget exceeded during rewriting')
        Q = self.pfacts(Q) if self.facts else Q
        if _depth < 4:
            # bit-level fallback: a left-over lowest sum bit s = x xor y of a half adder whose carry x*y is a product
            sub = {}
            for v in Q.vars():
                if v.__class__ is int and self.kind[v] == KS:
                    cell = self.sig_cell[v]
                    k = cell.cout
                    if k is not None and not cell.C and self.kind[k] == KA and self.val_atom(k).d != {(k,): 1}:
                        if not fallback and not any(len(m) == 2 and m in Q.d for m in self.val_atom(k).d):
                            continue                # (only where the product x*y is already in the polynomial)
                        try:
                            sub[v] = self.val(cell.X) + self.val(cell.Y) - self.val_atom(k) * 2
                        except Unknown:
                            pass
            if sub:
                out = Poly({}, Q.e)
                for m, c in Q.d.items():
                    t = Poly({tuple(v for v in m if v not in sub): c})
                    for v in m:
                        if v in sub:
                            t = t * sub[v]
                    out.iadd(Poly(t.d, Q.e))
                return self.reduce(out, assume, _depth + 1, fallback)
        return Q

    # ------------------------------------------------------------------ exhaustive evaluation over a small cut
    def children(self, a):
        k = self.kind[a]
        if k == KA:
            U, V = self.adef[a]
            return (U | V) - ONE
        if k == KK:
            t, C = self.adef[a]
            return (C | self.fs(t)) - ONE
        if k == KS:
            return self.adef[a] - ONE
        if k == KV:
            c = self.adef[a]
            return (c.X | c.Y | c.C) - ONE
        return EMPTY

    def tt_base(self, F):
        n = 1 << len(F)
        memo = {0: (1 << n) - 1}
        for i, a in enumerate(F):
            blk = ((1 << (1 << i)) - 1) << (1 << i)
            t = 0
            for off in range(0, n, 2 << i):
                t |= blk << off
            memo[a] = t
        return memo

    def tt_atom(self, a, memo):
        r = memo.get(a)
        if r is not None:
            return r
        k = self.kind[a]
        if k == KA:
            U, V = self.adef[a]
            r = self.tt_set(U, memo) & self.tt_set(V, memo)
        elif k == KK:
            t, C = self.adef[a]
            r = self.tt_atom(t, memo) ^ self.tt_set(C, memo)
        elif k == KS:
            r = self.tt_set(self.adef[a], memo)
        elif k == KV:
            c = self.adef[a]
            x, y, z = self.tt_set(c.X, memo), self.tt_set(c.Y, memo), self.tt_set(c.C, memo)
            r = (x & y) | (z & (x ^ y))
        else:
            raise Internal('input bit outside the cut')
        memo[a] = r
        return r

    def tt_set(self, S, memo):
        x = 0
        for a in S:
            x ^= self.tt_atom(a, memo)
        return x

    def tt_bound(self, P, maxf=10):
        vs = P.vars()
        if len(vs) > maxf or len(P.d) > 64 or any(v.__class__ is not int or self.kind[v] == KV for v in vs):
            return None
        F = set(vs) - ONE
        best = None
        kind = self.kind
        evals = 0
        for _ in range(64):
            if len(F) <= maxf and evals < 24:
                evals += 1
                Fl = sorted(F)
                memo = self.tt_base(Fl)
                vals = [0] * (1 << len(Fl))
                for m, c in P.d.items():
                    t = memo[0]
                    for v in m:
                        t &= self.tt_atom(v, memo)
                    while t:
                        b = t & -t
                        vals[b.bit_length() - 1] += c
                        t ^= b
                best = _meet(best, (-((-min(vals)) >> P.e), max(vals) >> P.e))
            cand = [a for a in F if kind[a] != KI]
            if not cand:
                break
            a = max(cand)
            F.discard(a)
            F |= self.children(a)
            if len(F) > 4 * maxf:
                break
        return best

    def tt_zero(self, D, maxf=14):
        F = set(D) - ONE
        kind = self.kind
        for _ in range(64):
            if len(F) <= maxf:
                memo = self.tt_base(sorted(F))
                if self.tt_set(D, memo) == 0:
                    self.stats['stageC'] += 1
                    return True
            cand = [a for a in F if kind[a] != KI]
            if not cand:
                break
            a = max(cand)
            F.discard(a)
            F |= self.children(a)
            if len(F) > 3 * maxf:
                break
        return False

    # ------------------------------------------------------------------ claims
    def word(self, wires):
        P = Poly()
        for i, w in enumerate(wires):
            S = self.read(w)
            if S:
                P.iadd(self.val(S), 1 << i)
        return P

    def prove_eq(self, P):
        R = self.reduce(P)
        if not R.d:
            return True, ''
        R = self.reduce(R, fallback=True)
        if not R.d:
            return True, ''
        for lane in range(min(self.L, 16)):
            if self.peval(P, lane):
                raise Refuted('an identity claim fails on a simulated input')
        return False, '%d monomials left after rewriting' % len(R.d)


def _same(ca, wa, cb, wb):
    if len(ca) != len(cb):
        return False
    for r, c in cb.items():
        x = ca.get(r)
        if x is None or x * wb != c * wa:
            return False
    return True


def _frac(cdict, den, e):
    """Poly with value (sum cdict) / (den * 2^e), den = +-2^k."""
    k = _pow2(den)
    if den < 0:
        return Poly({r: -c for r, c in cdict.items()}, e + k)
    return Poly(dict(cdict), e + k)


# ============================================================================ specification
_CMP = {ast.Eq: '==', ast.LtE: '<=', ast.Lt: '<', ast.GtE: '>=', ast.Gt: '>'}


def _parse_claim(s):
    try:
        node = ast.parse(s.strip(), mode='eval').body
    except SyntaxError as ex:
        raise ValueError('bad claim %r: %s' % (s, ex))
    if not isinstance(node, ast.Compare) or len(node.ops) != 1 or type(node.ops[0]) not in _CMP:
        raise ValueError('claim must be a single comparison (== <= < >= >): %r' % s)
    return _CMP[type(node.ops[0])], node.left, node.comparators[0]


def _ev(node, env):
    if isinstance(node, ast.Constant) and isinstance(node.value, int) and not isinstance(node.value, bool):
        return node.value
    if isinstance(node, ast.Name):
        if node.id not in env:
            raise ValueError('unknown word %r in a claim' % node.id)
        return env[node.id]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -_ev(node.operand, env)
    if isinstance(node, ast.BinOp):
        a, b = _ev(node.left, env), _ev(node.right, env)
        if isinstance(node.op, ast.Add):
            return a + b
        if isinstance(node.op, ast.Sub):
            return a - b
        if isinstance(node.op, ast.Mult):
            return a * b
        if isinstance(node.op, ast.Pow) and isinstance(b, int) and 0 <= b <= 4096 and (isinstance(a, int) or b <= 8):
            return a ** b
    raise ValueError('unsupported expression in a claim')


def _wl(x):
    if isinstance(x, dict):
        a, b = x['range']
        return list(range(int(a), int(b)))
    return [int(v) for v in x]


def normalize_spec(spec, N):
    """-> dict(inputs, outputs, claims [(text, op, lhs, rhs)], garbage set, keep set); ValueError if malformed."""
    if isinstance(spec, str):
        with open(spec) as fh:
            spec = json.load(fh)
    ins = {k: _wl(v) for k, v in spec.get('inputs', {}).items()}
    outs = {k: _wl(v) for k, v in spec.get('outputs', {}).items()}
    if set(ins) & set(outs):
        raise ValueError('a word name is used both as input and as output')
    inw = []
    for k, ws in ins.items():
        inw += ws
    if len(set(inw)) != len(inw):
        raise ValueError('input words overlap')
    allw = inw + [w for ws in outs.values() for w in ws] + _wl(spec.get('garbage', []))
    if any(w < 0 or w >= N for w in allw):
        raise ValueError('spec names a wire outside 0..%d' % (N - 1))
    for k, ws in outs.items():
        if len(set(ws)) != len(ws):
            raise ValueError('output word %s repeats a wire' % k)
    outw = set(w for ws in outs.values() for w in ws)
    gar = set(_wl(spec.get('garbage', [])))
    if 'keep' in spec:
        keep = set(_wl(spec['keep']))
        if not keep <= set(inw):
            raise ValueError('keep must list input wires')
    else:
        keep = set(inw) - outw - gar
    claims = []
    for s in spec.get('claims', []):
        op, l, r = _parse_claim(s)
        claims.append((s, op, l, r))
    if not claims:
        raise ValueError('spec without claims')
    return dict(inputs=ins, outputs=outs, claims=claims, garbage=gar, keep=keep, in_wires=inw)


# ============================================================================ independent simulation (refuter)
def _pre_sim(N, ops, sp, lanes=64, seed=11):
    """bit-sliced classical run on `lanes` inputs.  -> (status or None, detail)."""
    rng = random.Random(seed)
    M = (1 << lanes) - 1
    v = [0] * N
    words = list(sp['inputs'].items())
    for k, ws in words:
        for w in ws:
            v[w] = rng.getrandbits(lanes) & ~3 | 2
    # corner lanes: each input word all-ones with the others zero (lanes 2..), all-ones in lane 1, zeros in lane 0
    for j, (k, ws) in enumerate(words[:8]):
        bit = 1 << (2 + j)
        for k2, ws2 in words:
            for w in ws2:
                v[w] = (v[w] | bit) if k2 == k else (v[w] & ~bit)
    init = list(v)

    def inp(lane):
        return {k: sum(((init[w] >> lane) & 1) << i for i, w in enumerate(ws)) for k, ws in words}
    for i, o in enumerate(ops):
        k = o[0]
        if k == 'cx':
            v[o[2]] ^= v[o[1]]
        elif k == 'x':
            v[o[1]] ^= M
        elif k == 'and':
            if v[o[3]]:
                lane = (v[o[3]] & -v[o[3]]).bit_length() - 1
                return 'UNKNOWN', 'op %d: AND target q[%d] is not |0> on the input %s (precondition of the gadget ' \
                    'violated; the gate model does not apply)' % (i, o[3], inp(lane))
            v[o[3]] = v[o[1]] & v[o[2]]
        else:
            t = o[1]
            d = v[t]
            for c in o[2]:
                d ^= (v[c[1]] & v[c[2]]) if c[0] == 'cz' else v[c[1]]
            if d != 0 and d != M:
                l1 = (d & -d).bit_length() - 1
                l0 = ((~d & M) & -(~d & M)).bit_length() - 1
                return 'REFUTED', 'op %d: erasure of q[%d]: the phase of the outcome-1 branch depends on the input ' \
                    '(inputs %s and %s differ)' % (i, t, inp(l1), inp(l0))
            v[t] = 0
    outw = set(w for ws in sp['outputs'].values() for w in ws)
    for w in range(N):
        if w in sp['keep']:
            bad = v[w] ^ init[w]
        elif w in outw or w in sp['garbage']:
            continue
        else:
            bad = v[w]
        if bad:
            lane = (bad & -bad).bit_length() - 1
            return 'REFUTED', 'wire q[%d] ends %s on the input %s' % (
                w, 'changed' if w in sp['keep'] else 'not at |0>', inp(lane))
    for lane in range(lanes):
        env = inp(lane)
        for k, ws in sp['outputs'].items():
            env[k] = sum(((v[w] >> lane) & 1) << i for i, w in enumerate(ws))
        for text, op, l, r in sp['claims']:
            a, b = _ev(l, env), _ev(r, env)
            if not {'==': a == b, '<=': a <= b, '<': a < b, '>=': a >= b, '>': a > b}[op]:
                return 'REFUTED', 'claim "%s" fails on the input %s' % (text[:60], inp(lane))
    return None, ''


# ============================================================================ driver
def prove_ops(N, ops, spec, budget_sec=600.0, lanes=48):
    """exact proof of the op list (see parse_strict) against a word-level spec.  -> result dict."""
    t0 = time.time()
    res = dict(status='UNKNOWN', method='algebra', detail='')
    try:
        sp = normalize_spec(spec, N)
    except (ValueError, KeyError, OSError, TypeError) as ex:
        res.update(status='ERROR', detail='spec: %s' % ex)
        return res
    st, det = _pre_sim(N, ops, sp)
    res['counts'] = dict(qubits=N, ops=len(ops), AND=sum(1 for o in ops if o[0] == 'and'),
                         erasures=sum(1 for o in ops if o[0] == 'erase'))
    if st is not None:
        res.update(status=st, detail=det + ' [independent bit-sliced simulation]', time=round(time.time() - t0, 2))
        return res
    eng = Engine(N, sp['in_wires'], lanes=lanes, deadline=(t0 + budget_sec) if budget_sec else None)
    try:
        eng.run(ops)
        outw = set(w for ws in sp['outputs'].values() for w in ws)
        for w in range(N):
            S = eng.read(w)
            if w in sp['keep']:
                if not eng.is_zero(S ^ eng.fs(eng.in_atom[w])):
                    raise Unknown('input wire q[%d] is not proved unchanged' % w)
            elif w in outw or w in sp['garbage']:
                continue
            elif S and not eng.is_zero(S):
                raise Unknown('wire q[%d] is not proved to end at |0>' % w)
        env = {}
        for k, ws in sp['inputs'].items():
            env[k] = Poly({(eng.in_atom[w],): 1 << i for i, w in enumerate(ws)})
        for k, ws in sp['outputs'].items():
            env[k] = eng.word(ws)
        proved = []
        for text, op, l, r in sp['claims']:
            a, b = _ev(l, env), _ev(r, env)
            a = Poly.const(a) if isinstance(a, int) else a
            if op == '==':
                ok, why = eng.prove_eq(a - b)
                if not ok:
                    raise Unknown('claim "%s": %s' % (text[:60], why))
            else:
                if op in ('<=', '<'):
                    P = (b - a) if not isinstance(b, int) else (-a + b)
                else:
                    P = a - b
                if op in ('<', '>'):
                    P = P - 1
                lo, hi = eng.bound(P, 0)
                if lo < 0:
                    lo, hi = eng.bound(P, 2)
                if lo < 0:
                    raise Unknown('claim "%s": interval analysis gives only >= %d' % (text[:60], lo))
            proved.append(text if len(text) <= 80 else text[:77] + '...')
        res.update(status='PROVED', detail='all inputs, all measurement outcomes: %d AND gadgets and %d erasures with '
                   'their preconditions, %d word identities (carry chains), claims: %s' % (
                       eng.stats['ands'], eng.stats['erases'], eng.stats['adders'], '; '.join(proved)))
    except Unknown as ex:
        res.update(status='UNKNOWN', detail=str(ex))
    except Refuted as ex:
        res.update(status='REFUTED', detail=str(ex))
    except Internal as ex:
        res.update(status='ERROR', detail='internal self-check: %s' % ex)
    except ValueError as ex:
        res.update(status='ERROR', detail='spec: %s' % ex)
    res['stats'] = dict(eng.stats, atoms=len(eng.kind), classes=dict(eng.cls))
    res['time'] = round(time.time() - t0, 2)
    return res


def _peak_mb():
    try:
        import resource
        p = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round((p if sys.platform == 'darwin' else p * 1024) / 1e6, 1)
    except Exception:
        return None


def aprove(qasm, spec, budget_sec=600.0):
    """Word-level exact proof.  qasm: path or OpenQASM text; spec: dict or path of a JSON file:
        {"inputs":  {"x": [wires], "y": [wires]},          free input words (LSB first); all other wires start |0>
         "outputs": {"out": [wires]},                      words read at the end
         "claims":  ["out == x*y"],                        integer identities / inequalities over these words
         "garbage": [wires],                               wires that may end with any classical function
         "keep":    [wires]}                               input wires that must be unchanged (default: all input
                                                           wires that are neither output nor garbage)
    A wire list may be given as {"range": [a, b]} (a..b-1).  -> dict(status PROVED | REFUTED | UNKNOWN | ERROR,
    method='algebra', detail, time, counts, stats, peak_mb)."""
    t0 = time.time()
    try:
        text = qasm if ('\n' in qasm) else open(qasm).read()
        N, ops = parse_strict(text)
    except (ParseError, OSError) as ex:
        return dict(status='ERROR', method='algebra', detail='algebra engine does not apply: %s' % ex,
                    time=round(time.time() - t0, 2))
    del text
    out = {}

    def run():
        try:
            out['r'] = prove_ops(N, ops, spec, budget_sec=budget_sec)
        except RecursionError:
            out['r'] = dict(status='UNKNOWN', method='algebra', detail='recursion limit')
        except MemoryError:
            out['r'] = dict(status='UNKNOWN', method='algebra', detail='out of memory')
    old = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old, 200000))
    try:
        threading.stack_size(256 * 1024 * 1024)
        th = threading.Thread(target=run)
        th.start()
        th.join()
    finally:
        sys.setrecursionlimit(old)
        threading.stack_size(0)
    r = out.get('r') or dict(status='ERROR', method='algebra', detail='prover thread died')
    r['total_time'] = round(time.time() - t0, 2)
    r['peak_mb'] = _peak_mb()
    return r


# ============================================================================ specs of the qoptct.crypto generators
def spec_mulsq(circ):
    """spec of a qoptct.crypto.mulsq.MulSqCircuit: out == x*y (x*x), inputs unchanged, ancillas |0>, garbage free."""
    if circ.kind == 'mul':
        return dict(inputs=dict(x=list(circ.x), y=list(circ.y)), outputs=dict(out=list(circ.out)),
                    claims=['out == x*y'], garbage=list(circ.garbage))
    return dict(inputs=dict(x=list(circ.x)), outputs=dict(out=list(circ.out)), claims=['out == x*x'],
                garbage=list(circ.garbage))


def spec_fold(R):
    """spec of a qoptct.crypto.fold.FoldReduction: r == P mod p with 0 <= r < p (stated with the quotient made of the
    garbage registers), H unchanged, ancillas |0>."""
    n, p, c = R.n, R.p, R.c
    lay = R.layout
    outs = dict(r=list(lay['r']), thi1=list(lay['thi1']), g=list(lay['g']))
    gar = list(lay['thi1']) + list(lay['g'])
    claims = []
    if 'M' in lay:
        outs['M'] = list(lay['M'])
        gar += list(lay['M'])
        claims.append('M == %d*thi1' % c)
    if R.method == 'solinas':
        from qoptct.crypto.fold import p256_formula
        ins = {'A%d' % i: list(range(32 * i, 32 * i + 32)) for i in range(16)}
        qs = []
        for i in range(16):
            d = (1 << (32 * i)) - p256_formula(1 << (32 * i))
            if d % p:
                raise ValueError('internal: NIST formula coefficient')
            qs.append('%d*A%d' % (d // p, i))
        claims.append('r + %d*(thi1 - 1 + g - 5 + %s) == %s' % (
            p, ' + '.join(qs), ' + '.join('%d*A%d' % (1 << (32 * i), i) for i in range(16))))
    else:
        ins = dict(L=list(range(n)), H=list(range(n, 2 * n)))
        claims.append('r + %d*(H + thi1 - 1 + g) == L + %d*H' % (p, 1 << n))
    claims.append('r <= %d' % (p - 1))
    return dict(inputs=ins, outputs=outs, claims=claims, garbage=gar, keep=list(range(n, 2 * n)))


# ============================================================================ command line
def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(prog='qoptct aprove', description='Word-level algebraic exact proof of an integer-'
                                 'arithmetic circuit on logical ANDs against a JSON specification (see '
                                 'qoptct.verify.algebra).',
                                 epilog='Exit code: 0 PROVED, 1 REFUTED, 2 UNKNOWN (not proved), 3 usage or input '
                                        'error.')
    ap.add_argument('qasm')
    ap.add_argument('spec', help='JSON file: inputs, outputs, claims, garbage, keep')
    ap.add_argument('--budget', type=float, default=3600.0, help='time budget in seconds (default 3600; 0 = none)')
    ap.add_argument('--json', action='store_true')
    try:
        a = ap.parse_args(argv)
    except SystemExit as e:
        return 3 if e.code else 0
    for p in (a.qasm, a.spec):
        if not os.path.isfile(p):
            print('error: no such file: %s' % p, file=sys.stderr)
            return 3
    r = aprove(a.qasm, a.spec, budget_sec=a.budget or None)
    if a.json:
        print(json.dumps(r, default=str))
    else:
        print('%s  (algebra, %.1f s, peak %s MB)  %s' % (r['status'], r.get('total_time', 0.0), r.get('peak_mb'),
                                                          r.get('counts', '')))
        print('  ' + str(r.get('detail')))
    return {'PROVED': 0, 'REFUTED': 1, 'ERROR': 3}.get(r['status'], 2)


if __name__ == '__main__':
    sys.exit(main())
