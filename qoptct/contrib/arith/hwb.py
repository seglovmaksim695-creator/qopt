# SPDX-License-Identifier: Apache-2.0
"""Hidden weighted bit: y[(i + w) mod n] = x[i], w = popcount(x)  (cyclic rotation of the data by its weight).

Recognition: one-bit probes give the successor of every data wire (a single bit moves by one position); the
successors must form one n-cycle, and the whole truth table must equal the definition.
Synthesis (lab agent b4): weight counter (full/half adders, one AND each, garbage erased at once),
staged barrel shifter (controlled cyclic shifts by 2^k made of controlled swaps, one AND per swap), and erasure of
the weight bits WITHOUT a second counter: (-1)^{bit k of w} is a phase polynomial of the (rotated) data wires,
applied as T/S/Z corrections after the X-basis measurement.  'classes': counters of even and odd positions (both
invariant under even shifts), so that only weight bits <= 2 are erased by phase gadgets (needed when n > 8).
"""
import math

import numpy as np

from .ir import Builder, BuildError
from .fam import flips, confirm

_GN = {1: 't', 2: 's', 4: 'z', 6: 'sdg', 7: 'tdg'}


def _words_from_bool(T, arr):
    b = np.packbits(arr.astype(np.uint8), bitorder='little')
    if b.size < 8:
        b = np.concatenate([b, np.zeros(8 - b.size, dtype=np.uint8)])
    return T.mask(b.view('<u8').astype(np.uint64))


def hwb_expect(T, ring):
    n = len(ring)
    k = T.k
    xs = np.arange(1 << k, dtype=np.uint64)
    bits = [((xs >> np.uint64(T.pos[w])) & np.uint64(1)).astype(np.int64) for w in ring]
    wt = np.zeros(1 << k, dtype=np.int64)
    for b in bits:
        wt += b
    exp = {}
    for p in range(n):
        y = np.zeros(1 << k, dtype=bool)
        for i in range(n):
            y |= (bits[i] == 1) & (((i + wt) % n) == p)
        exp[ring[p]] = _words_from_bool(T, y)
    return exp


def rec_hwb(T):
    if T.k > 20:
        return []
    if flips(T, []):
        return []
    succ = {}
    for u in T.vars:
        d = flips(T, [u]) ^ {u}
        if len(d) != 1:
            return []
        v = d.pop()
        if v not in T.inp or v == u:
            return []
        succ[u] = v
    if len(set(succ.values())) != len(succ):
        return []
    ring = [T.vars[0]]
    while True:
        nx = succ[ring[-1]]
        if nx == ring[0]:
            break
        ring.append(nx)
        if len(ring) > len(succ):
            return []
    n = len(ring)
    if n != len(succ) or n < 3:
        return []
    if not confirm(T, hwb_expect(T, ring)):
        return []
    return [dict(kind='hwb', ring=ring)]


# ---------------------------------------------------------------- builder helpers
class HB:
    def __init__(self, T):
        self.T = T
        self.bd = Builder(T)
        self.alive = set()

    def _int(self, w):
        return int.from_bytes(self.bd.val[w].astype('<u8').tobytes(), 'little')

    def alloc(self):
        w = self.bd.alloc()
        self.alive.add(w)
        return w

    def is_alive(self, w):
        return w in self.alive

    def xor_into_new(self, ws):
        t = self.alloc()
        for w in ws:
            self.bd.cx(w, t)
        return t

    def AND(self, a, b):
        t = self.alloc()
        self.bd.AND(a, b, t)
        return t

    def maj_and(self, a, b, c):
        bd = self.bd
        bd.cx(c, a)
        bd.cx(c, b)
        t = self.AND(a, b)
        bd.cx(c, a)
        bd.cx(c, b)
        bd.cx(c, t)
        return t

    def erase(self, w, corr):
        self.bd.erase(w, corr)
        self.alive.discard(w)

    def erase_auto(self, w):
        """find a correction (linear Z's, then CZ pairs) for the value of w over the live wires, exactly."""
        T = self.T
        L = 1 << T.k
        full = (1 << L) - 1
        target = self._int(w)
        live = [q for q in T.vars] + sorted(q for q in self.alive if q != w)
        cols = [(('z', q), self._int(q)) for q in live] + [(('c',), full)]
        r = _solve(target, cols)
        if r is None:
            vals = {q: self._int(q) for q in live}
            pcols = [(('cz', a, b), vals[a] & vals[b]) for i, a in enumerate(live) for b in live[i + 1:]]
            pcols = [pc for pc in pcols if pc[1]]
            for pc in pcols:                       # fewest CZ first: one pair, then two pairs
                r = _solve(target, cols + [pc])
                if r is not None:
                    break
            if r is None and len(pcols) <= 400:
                for i, p1 in enumerate(pcols):
                    for p2 in pcols[i + 1:]:
                        r = _solve(target, cols + [p1, p2])
                        if r is not None:
                            break
                    if r is not None:
                        break
            if r is None:
                r = _solve(target, cols + pcols)
        if r is None:
            raise BuildError('wire %d not erasable' % w)
        corr = []
        for key in r:
            if key[0] == 'z':
                corr.append(('z', (key[1],)))
            elif key[0] == 'cz':
                corr.append(('cz', (key[1], key[2])))
        self.erase(w, corr)

    def cswap(self, c, i, j):
        bd = self.bd
        bd.cx(i, j)
        t = self.AND(c, j)
        bd.cx(t, i)
        self.erase(t, [('cz', (c, j))])
        bd.cx(i, j)

    def ctrl_cycle(self, c, pos):
        for k in range(1, len(pos)):
            self.cswap(c, pos[0], pos[k])

    def reduce(self, levels, nbits, garbage):
        levels = [list(l) for l in levels]
        out = []
        j = 0
        while j < nbits:
            while len(levels) <= j + 1:
                levels.append([])
            L = levels[j]
            if j == nbits - 1:
                if len(L) == 1:
                    out.append(L[0])
                else:
                    w = self.xor_into_new(L)
                    garbage.append(w)
                    out.append(w)
                break
            queue = list(L)
            while len(queue) >= 3:
                a, b, c = queue[:3]
                queue = queue[3:]
                sm = self.xor_into_new([a, b, c])
                k = self.maj_and(a, b, c)
                garbage += [sm, k]
                queue.append(sm)
                levels[j + 1].append(k)
            if len(queue) == 2:
                sm = self.xor_into_new(queue)
                k = self.AND(queue[0], queue[1])
                garbage += [sm, k]
                queue = [sm]
                levels[j + 1].append(k)
            if not queue:
                z = self.alloc()
                garbage.append(z)
                out.append(z)
            else:
                out.append(queue[0])
            j += 1
        return out

    def counter(self, wires, nbits):
        G = []
        out = self.reduce([wires], nbits, G)
        for w in reversed(G):
            if w in out or not self.is_alive(w):
                continue
            self.erase_auto(w)
        return out

    def erase_weight_bit(self, w, k, wires, low):
        if k > 2:
            raise BuildError('weight bit > 2')
        units = [(q, 4 >> k) for q in wires] + [(low[j], (-(1 << j) * (4 >> k)) % 8) for j in range(k)]
        corr = [(_GN[u], (q,)) for q, u in units if u % 8]
        self.erase(w, corr)


def _solve(target, cols):
    """GF(2): subset of cols (key, int) whose XOR is target -> list of keys, or None."""
    basis = {}   # pivot bit -> (vector, set of keys)
    for key, v in cols:
        comb = {key}
        while v:
            p = v.bit_length() - 1
            if p in basis:
                bv, bc = basis[p]
                v ^= bv
                comb ^= bc
            else:
                basis[p] = (v, comb)
                break
    v = target
    comb = set()
    while v:
        p = v.bit_length() - 1
        if p not in basis:
            return None
        bv, bc = basis[p]
        v ^= bv
        comb ^= bc
    return sorted(comb, key=str)


def _cycles(n, shift):
    g = math.gcd(n, shift)
    return [[(r + k * shift) % n for k in range(n // g)] for r in range(g)]


def gen_generic(T, h):
    ring = h['ring']
    n = len(ring)
    nb = (n - 1).bit_length()
    if nb > 3:
        raise BuildError('generic hwb needs n <= 8')
    H = HB(T)
    W = H.counter(list(ring), nb)
    for k in range(1, nb):
        sh = (1 << k) % n
        if sh == 0:
            continue
        for cyc in _cycles(n, sh):
            H.ctrl_cycle(W[k], [ring[i] for i in cyc])
    for cyc in _cycles(n, 1):
        H.ctrl_cycle(W[0], [ring[i] for i in cyc])
    for k in reversed(range(nb)):
        H.erase_weight_bit(W[k], k, list(ring), W[:k])
    return H.bd


def gen_classes(T, h):
    ring = h['ring']
    n = len(ring)
    if n % 2:
        raise BuildError('classes need even n')
    H = HB(T)
    ev = ring[0::2]
    od = ring[1::2]
    bl = (n // 2).bit_length()
    if bl > 3:
        raise BuildError('class counters too wide')
    E = H.counter(ev, bl)
    O = H.counter(od, bl)
    nb = (n - 1).bit_length()
    G3 = []
    levels = [[E[j], O[j]] if j < bl else [] for j in range(nb)]
    W = H.reduce(levels, nb, G3)
    for k in range(1, nb):
        sh = (1 << k) % n
        if sh == 0:
            continue
        for cyc in _cycles(n, sh):
            H.ctrl_cycle(W[k], [ring[i] for i in cyc])
    for w in reversed(G3):
        if w == W[0] or not H.is_alive(w):
            continue
        H.erase_auto(w)
    for cls, Z in ((ev, E), (od, O)):
        for k in reversed(range(bl)):
            H.erase_weight_bit(Z[k], k, cls, Z[:k])
    for cyc in _cycles(n, 1):
        H.ctrl_cycle(W[0], [ring[i] for i in cyc])
    if H.is_alive(W[0]):
        H.erase_auto(W[0])
    return H.bd


def variants_hwb(T, h):
    return [('weight-barrel', lambda: gen_generic(T, h)), ('class-barrel', lambda: gen_classes(T, h))]
