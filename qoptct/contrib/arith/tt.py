# SPDX-License-Identifier: Apache-2.0
"""Truth tables of classical reversible .qc circuits (bit-sliced, numpy uint64 words).

A .qc circuit of the benchmark suites is written with H-frames (``H t; Z a b t; H t`` = Toffoli). ``lower`` turns it
into a classical gate list (x, cx, ccx); circuits with genuine phases or entangling H-structure are rejected.

Model ``clean``: wires not listed on the ``.i`` line start in |0>.  A clean wire whose very first gate is ``H`` and
that cannot be lowered otherwise is treated as a "plus" wire: the circuit is ``H(w)`` followed by a classical circuit
in which ``w`` is an extra input variable (this is how ``qcla_com_7`` leaves a carry wire in |+>).
Model ``strict``: every wire is an input.
"""
import itertools

import numpy as np

import qoptct

ONES = np.uint64(0xFFFFFFFFFFFFFFFF)
_LOW = [0xAAAAAAAAAAAAAAAA, 0xCCCCCCCCCCCCCCCC, 0xF0F0F0F0F0F0F0F0, 0xFF00FF00FF00FF00,
        0xFFFF0000FFFF0000, 0xFFFFFFFF00000000]


class NotClassical(Exception):
    pass


def lower(hl, n, skip_first_h=()):
    """H-frame lowering of qoptct.load_qc gates -> list of (name, qubits) with names x, cx, ccx.
    ``skip_first_h``: wires whose first H is dropped (plus wires). Raises NotClassical."""
    fr = [0] * n
    out = []
    pending = set(skip_first_h)
    for g, q in hl:
        q = tuple(q)
        if g == 'h':
            if q[0] in pending:
                pending.discard(q[0])
                continue
            fr[q[0]] ^= 1
            continue
        pending.difference_update(q)
        if g in ('ccz', 'cz', 'z'):
            f = [x for x in q if fr[x]]
            if len(f) != 1:
                raise NotClassical('phase gate %s' % g)
            t = f[0]
            c = tuple(x for x in q if x != t)
            out.append(({'ccz': 'ccx', 'cz': 'cx', 'z': 'x'}[g], c + (t,)))
        elif g == 'x':
            if fr[q[0]]:
                raise NotClassical('Z phase')
            out.append(('x', q))
        elif g in ('cx', 'ccx'):
            if any(fr[x] for x in q):
                raise NotClassical('gate inside a frame')
            out.append((g, q))
        else:
            raise NotClassical(g)
    if any(fr):
        raise NotClassical('open frame')
    return out


def input_names(path):
    for line in open(path):
        if line.startswith('.i'):
            return line.split()[1:]
    return None


class Table:
    """Bit-sliced truth table. ``vars``: wire indices that are input variables (variable j = wire vars[j]);
    ``out[w]``: words of the output value of wire w; ``plus``: plus wires (subset of vars);
    ``clean``: wires that start in |0> (not variables)."""

    def __init__(self, names, gates, vars_, plus, clean_wires):
        self.names = names
        self.n = len(names)
        self.gates = gates
        self.vars = list(vars_)
        self.plus = list(plus)
        self.clean = list(clean_wires)
        self.k = len(self.vars)
        self.W = max(1, (1 << self.k) >> 6)
        self.valid = ONES if self.k >= 6 else np.uint64((1 << (1 << self.k)) - 1)
        self.inp = {w: self.var(j) for j, w in enumerate(self.vars)}
        self.pos = {w: j for j, w in enumerate(self.vars)}
        zero = self.zeros()
        st = [self.inp.get(w, zero).copy() for w in range(self.n)]
        self.out = self.run(gates, st)

    def run(self, gates, st):
        st = list(st)
        for g, q in gates:
            if g == 'x':
                st[q[0]] = ~st[q[0]]
            elif g == 'cx':
                st[q[1]] = st[q[1]] ^ st[q[0]]
            else:
                st[q[2]] = st[q[2]] ^ (st[q[0]] & st[q[1]])
        return [self.mask(s) for s in st]

    def zeros(self):
        return np.zeros(self.W, dtype=np.uint64)

    def ones(self):
        return self.mask(np.full(self.W, ONES, dtype=np.uint64))

    def mask(self, a):
        if self.k >= 6:
            return a
        return a & self.valid

    def var(self, j):
        if j < 6:
            return self.mask(np.full(self.W, np.uint64(_LOW[j]), dtype=np.uint64))
        idx = np.arange(self.W, dtype=np.uint64)
        return np.where((idx >> np.uint64(j - 6)) & np.uint64(1), ONES, np.uint64(0)).astype(np.uint64)

    def at(self, w, x):
        """value (0/1) of output wire w on input assignment x (int over variables)."""
        return int((int(self.out[w][x >> 6]) >> (x & 63)) & 1)

    def assign(self, wires_on):
        """input assignment (int) with the given variable wires set to 1."""
        x = 0
        for w in wires_on:
            x |= 1 << self.pos[w]
        return x

    def outset(self, x):
        """set of wires that are 1 at the output on input assignment x."""
        return {w for w in range(self.n) if self.at(w, x)}

    def eq(self, a, b):
        return not np.any(self.mask(a ^ b))

    def is_zero(self, a):
        return not np.any(self.mask(a))

    def changed(self):
        """wires whose output differs from their input value (clean wires: from 0)."""
        res = []
        for w in range(self.n):
            iv = self.inp.get(w)
            if iv is None:
                if not self.is_zero(self.out[w]):
                    res.append(w)
            elif not self.eq(iv, self.out[w]):
                res.append(w)
        return res


def load(path, clean, max_vars=24):
    """-> Table or None (not classical / too many variables)."""
    try:
        names, hl = qoptct.load_qc(path)
    except Exception:
        return None
    n = len(names)
    ins = input_names(path) if clean else None
    if ins is None:
        var_w = list(range(n))
    else:
        s = set(ins)
        var_w = [i for i, nm in enumerate(names) if nm in s]
    clean_w = [i for i in range(n) if i not in set(var_w)]
    plus = []
    try:
        gates = lower(hl, n)
    except NotClassical:
        if not clean_w:
            return None
        first = {}
        for g, q in hl:
            for x in q:
                first.setdefault(x, g)
        cand = [w for w in clean_w if first.get(w) == 'h']
        gates = None
        for size in (1, 2):
            for sub in itertools.combinations(cand, size):
                try:
                    gates = lower(hl, n, skip_first_h=sub)
                except NotClassical:
                    continue
                plus = list(sub)
                break
            if gates is not None:
                break
        if gates is None:
            return None
    vars_ = sorted(var_w + plus)
    if len(vars_) > max_vars:
        return None
    clean_w = [w for w in clean_w if w not in plus]
    return Table(names, gates, vars_, plus, clean_w)
