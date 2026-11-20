# SPDX-License-Identifier: Apache-2.0
"""Narrow table lookups  z ^= T[a]  (address a: c bits, word: w bits) on logical ANDs, checked exactly at build time.

Cost model: logical AND of Gidney (arXiv:1709.06648) = 4 T, 3 CX into a clean qubit; erasure by X-measurement and a
Clifford (CZ) correction = 0 T.  "ANDs" below = Toffoli-equivalents in this model.

Constructions
-------------
narrow_lookup(table, c, w)      BILINEAR split of the address into hi (h bits) and lo (l = c-h bits):
    T(a) = sum_{S,U} C[S,U] M_S(hi) N_U(lo)  (ANF, M_S / N_U = monomials of the two halves).
    1. all monomials of degree >= 2 of each half:            (2^h-h-1) + (2^l-l-1) ANDs  ("base", 8 for c=6, h=3)
    2. terms with S=0 or U=0 are XORs of monomial wires: free (CX fan-out);
    3. the mixed part (core matrices C_k, one per output bit, (2^h-1) x (2^l-1)) is covered by R products
       f*g, f = alpha.M + beta.N, g = alpha'.M + beta'.N (ONE AND each); f*g has the core alpha beta'^T + alpha' beta^T
       (rank <= 2) and pure hi/lo parts (free).  For fixed (alpha, alpha') the condition "core in V = span{C_k}" is
       linear in (beta, beta'), so candidates are kernels of small GF(2) systems; a basis of V among them gives
       R = dim V (= w when w <= rank).  Helpers outside V are added if no basis is found (R > dim V).
    ANDs = base + R  >=  base + dim V.   c=6, w=32: 8 + 32 = 40 (ANF/Davio and Gidney 2025's count 2^c-c-1 = 57);
    c=6, w=21: 33 with h=2 (12 + 21); c=5, w=16: 21 instead of 26.  Splits h = 1..min(3, c//2) are tried, the
    cheapest wins; if no split beats Davio (2^c-c-1), the Davio lookup is returned (method 'davio').
    Measured with the defaults on random tables (Random(10000+t), lab script, qoptct 0.5):
        c=6 w=32: 40 ANDs (= 8 + w) in 38/40 tables, 41 in 2/40 (h=3)          Davio 57
        c=6 w=21: 33 (= 12 + w) in 15/20, 34 in 4/20, 35 in 1/20 (h=2)         Davio 57
        c=6 w=16: 29..31 in 19/20, 34 in 1/20 (lower bound 28)                  Davio 57
        c=5 w=16: 21 in 37/40, 20 in 2/40, 22 in 1/40;  c=5 w=8: 13..15         Davio 26
        c=4 w=4: 6 in 35/40, 7 in 5/40;  c=4 w=8: 8..10, Davio (11) in 1/40      Davio 11
        c=7 w=32: 64..65 (5 tables)                                              Davio 120
    "Success" (R = dim V, i.e. base + w for w <= rank) is thus ~95 % at c=6 w=32; otherwise 1-2 ANDs more; the
    Davio fallback never triggered for c >= 5 in these samples.  The classical search is pure Python: c=6 w=32
    0.25 s per table (max 0.4 s), c=6 w=16..21 ~2-3 s, c=7 w=32 ~25 s, c <= 5 < 0.01 s: fine for a few tables;
    for the ~10^7 tables of a full RSA run a C implementation would be needed (not provided).  Ancillas: base
    monomial wires + 1 product wire (c=6: h=3 9 extra qubits, h=2 13; Davio: c-1 = 5).
davio_lookup(table, c, w)       ANF subset tree: one AND per needed monomial of degree >= 2 (<= 2^c-c-1 ANDs).
qroam1_lookup(table, c, w)      QROAM with one select bit (lambda = 2): ANF lookup of [T[2h] | T[2h]^T[2h+1]] over the
    c-1 high bits into [z | R], then z_i ^= a0 R_i (w ANDs), then R is cleared by measurement-based unlookup:
    (2^(c-1)-c) + w + U(c-1) ANDs  (c=6, w=21: 26+21+5 = 52), +w qubits while it runs.
unlookup_meas / lookup_pair     measurement-based uncomputation (X-measure z, phase fix-up (-1)^{m.T[a]} by products
    of hi x lo monomials): U(c) = (2^ceil(c/2)-ceil(c/2)-1) + (2^floor(c/2)-floor(c/2)-1) ANDs (c=6: 8).
    lookup_pair(table, c, w, mode): y ^= T[a] through a temporary z: lookup, CX, unlookup ('qroam1' clears [z | R]
    with ONE unlookup of the combined table: c=6, w=21: 26+21+8 = 55 instead of 57+8 = 65; 'narrow': 33+8 = 41).
fanout='gray' (davio_lookup, lookup_add; emit_davio_gray)   CX fan-out by a Gray-code accumulator ("Russian fours",
    idea of lab agent a4, measured by k1): the products M_S N_T (hi / lo monomials, 1 AND each) are made in batches of
    b wires that are alive together and fanned into the z register through ONE accumulator wire walked along the
    reflected Gray code: 2^b + w CX per batch instead of ~ b w / 2.  ANDs are those of the Davio lookup (never more; a
    sparse table can save one or two), T is unchanged, +1..+30 qubits (all lo monomials stay alive, + the accumulator).
    For wide words only (c = 10, w = 1024, dense random table: 528 201 -> 167 231 CX = -68 %, 1013 ANDs either way;
    c = 8, w = 256: 33 520 -> 16 192; c = 6, w = 64: 2 243 -> 1 444); batch by batch the direct fan-out is kept when it is
    cheaper, so narrow words (w <= 16) are not worse.  The default fanout='plain' reproduces the previous circuits
    byte for byte.  The -78 % of the lab claim (186 514 vs 1 077 908 CX per lookup-add at n = 2048, c = 9..10) is for the
    whole lookup-addition of the lab's runway-chain structure, which this module does not contain; here only the
    lookup itself is built.
lookup_add(table, c, w)         acc(w) += T[a] mod 2^w: narrow lookup + Gidney adder (w-1 ANDs) + unlookup
    (c=6, w=32: 40 + 31 + 8 = 79; one step of loop1 in Gidney 2025).

Priority and scope (honest attribution)
---------------------------------------
* Using a narrow (sliced) lookup inside a windowed lookup-addition (Gidney-Ekera framework) is prior work:
  Luongo, Narasimhachar, Sireesh 2025 (arXiv:2502.17325, App. B, "sliced windowing").
* QROAM / select-swap with lambda = 2 is Berry et al. 2019 (arXiv:1902.02134) and Low, Kliuchnikov, Schaeffer 2018
  (arXiv:1812.00954); only its use for these narrow 6-bit lookups is new here.  Measurement-based unlookup:
  Berry et al. 2019, Gidney 2019 (arXiv:1905.07682).
* The bilinear decomposition of one lookup (mixed products, R = dim V) is the lab's own construction (k1); we did not
  find it published, but no exhaustive literature search was made.
* No gain in the structure of Sharma-Lau 2026 (arXiv:2611.08907): a shared lookup for 3 primes (~3 x 32 output bits)
  has a full-rank core, so the bilinear split costs as much as ANF there (57/3 ~ 19 ANDs per prime < 40).
* Gidney 2025 (arXiv:2505.15917, Sec. 3.1) counts a lookup as 2^c-c-1; the savings quoted above are against that count
  for one lookup, not end-to-end figures.

Checks
------
Every builder runs the circuit on a batch of classical inputs (ALL 2^c addresses x z_samples values of z/acc) for
several random measurement-outcome patterns, tracking the phase exactly (units of pi/4): outputs, clean ancillas and a
constant global phase are checked; every AND erasure must cancel its phase on the whole batch.  The emitted OpenQASM 3
text (qoptct subset, parity-conditioned blocks ``if (c[i] ^ c[j]) {...}``) can be checked independently against
``reference_qc`` with qoptct.verify.prove (exact, small sizes) or qoptct.verify.arith.check_sim.

All public functions are pure: no global state, no environment variables, deterministic for a given seed; the input
table is not modified.
"""
import random
import time
from dataclasses import dataclass, field

import numpy as np

from qoptct.contrib.arith import ir as _ir

__all__ = ['LookupBuildError', 'Circ', 'Plan', 'LookupCircuit', 'anf_table', 'plan_narrow', 'emit_narrow',
           'emit_davio', 'emit_davio_gray', 'gray_fanout_cx', 'emit_qroam1', 'emit_unlookup_meas', 'emit_add_gidney', 'narrow_lookup', 'davio_lookup',
           'qroam1_lookup', 'lookup_pair', 'lookup_add', 'reference_qc', 'davio_ands', 'unlookup_ands',
           'qroam1_ands', 'base_ands']

_PHU = {'z': 4, 's': 2, 'sdg': 6, 't': 1, 'tdg': 7}


class LookupBuildError(Exception):
    """A build-time exact check failed (an internal error: no circuit is returned)."""


# ================================================================================================ circuit builder
class Circ:
    """Circuit on logical ANDs with an exact batch simulation of B classical inputs (bit arrays) and of the phase
    (units of pi/4, per input).  X-measurement outcomes are drawn from Random(seed); the emitted gate list does not
    depend on them (corrections are classically conditioned), so several seeds check several outcome patterns."""

    def __init__(self, B, seed=1):
        self.B = B
        self.val = {}
        self.N = 0
        self.free = []
        self.ops = []
        self.ph = np.zeros(B, np.int64)
        self.nand = 0
        self.nbits = 0
        self.m = {}
        self.rng = random.Random(seed)
        self.live = 0
        self.peak = 0
        self.inputs = []

    # ---------------- wires
    def inp(self, arr):
        w = self.N
        self.N += 1
        self.val[w] = np.asarray(arr, np.uint8).copy()
        self.inputs.append(w)
        self._grow()
        return w

    def new(self):
        if self.free:
            w = self.free.pop()
        else:
            w = self.N
            self.N += 1
        self.val[w] = np.zeros(self.B, np.uint8)
        self._grow()
        return w

    def _grow(self):
        self.live += 1
        self.peak = max(self.peak, self.live)

    def rel(self, w):
        if self.val[w].any():
            raise LookupBuildError('release of dirty wire %d' % w)
        self.free.append(w)
        self.live -= 1

    # ---------------- gates
    def x(self, t):
        self.ops.append(('x', (t,)))
        self.val[t] ^= 1

    def cx(self, c, t):
        if c == t:
            raise LookupBuildError('cx on one wire')
        self.ops.append(('cx', (c, t)))
        self.val[t] ^= self.val[c]

    def AND(self, a, b, t=None):
        """t = a AND b into a clean wire (4 T)."""
        if t is None:
            t = self.new()
        if len({a, b, t}) != 3:
            raise LookupBuildError('AND wires')
        if self.val[t].any():
            raise LookupBuildError('AND target not clean')
        self.ops.append(('and', (a, b, t, 'ab')))
        self.val[t] = (self.val[a] & self.val[b]).astype(np.uint8)
        self.nand += 1
        return t

    def _corrph(self, corr):
        p = np.zeros(self.B, np.int64)
        for g, q in corr:
            v = self.val[q[0]] & (self.val[q[1]] if g == 'cz' else 1)
            p += _PHU.get(g, 4) * v.astype(np.int64)
        return p

    def erase(self, t, corr):
        """X-measure t; if the outcome is 1 apply corr (diagonal Cliffords) and X(t).  Accepted only if
        4*v(t) + phase(corr) is constant mod 8 on the whole batch (both outcomes are then valid)."""
        if any(t in q for _, q in corr):
            raise LookupBuildError('correction acts on the measured wire')
        p = (4 * self.val[t].astype(np.int64) + self._corrph(corr)) % 8
        if (p != p[0]).any():
            raise LookupBuildError('erase of wire %d: correction does not cancel the phase' % t)
        self.ops.append(('erase', (t, tuple(corr))))
        self.nbits += 1
        self.val[t] = np.zeros(self.B, np.uint8)
        self.rel(t)

    def erase_and(self, t, a, b, lin=()):
        """erase t holding a AND b XOR (XOR of lin), current values."""
        self.erase(t, [('cz', (a, b))] + [('z', (w,)) for w in lin])

    def mx(self, t):
        """X-basis measurement of a data wire t (reset by a conditional X): phase (-1)^{m v(t)} is left for the
        caller to fix (ifpar).  Returns the classical bit index."""
        k = self.nbits
        self.nbits += 1
        m = self.rng.getrandbits(1)
        self.m[k] = m
        self.ph += 4 * m * self.val[t].astype(np.int64)
        self.val[t] = np.zeros(self.B, np.uint8)
        self.ops.append(('mx', (t, k)))
        self.rel(t)
        return k

    def ifpar(self, ks, gates):
        """apply diagonal Clifford gates if the parity of the classical bits ks is 1."""
        ks = tuple(sorted(ks))
        if not ks or not gates:
            return
        self.ops.append(('ifp', (ks, tuple(gates))))
        if sum(self.m[k] for k in ks) & 1:
            self.ph += self._corrph(gates)

    # ---------------- checks / output
    def finish(self, expect):
        """expect: wire -> bit array.  Other non-input wires must be 0, inputs not listed unchanged (they are never
        written by the builders unless listed); the global phase must be constant over the batch."""
        for w in range(self.N):
            if w in expect:
                if (self.val[w] != expect[w]).any():
                    raise LookupBuildError('wire %d wrong on %d inputs' % (w, int((self.val[w] != expect[w]).sum())))
            elif w not in self.inputs and self.val[w].any():
                raise LookupBuildError('wire %d not returned to 0' % w)
        p = self.ph % 8
        if (p != p[0]).any():
            raise LookupBuildError('data-dependent phase')
        return True

    def counts(self):
        T = CX = 0
        for g, q in self.ops:
            if g == 'and':
                T += 4
                CX += 3
            elif g == 'cx':
                CX += 1
            elif g == 'erase':
                CX += sum(1 for gg, _ in q[1] if gg == 'cz')
            elif g == 'ifp':
                CX += sum(1 for gg, _ in q[1] if gg == 'cz')
        return dict(AND=self.nand, T=T, CX=CX, qubits=self.N, peak=self.peak)

    def qasm(self, n_qubits=None, peephole=True):
        """OpenQASM 3 text (qoptct subset).  n_qubits: declare at least this many qubits (idle extra wires, e.g. to
        match a reference with more clean wires)."""
        L = []
        segs = []
        cur = []
        nb = 0
        for g, q in self.ops:
            if g in ('and', 'erase'):
                body, _ = _ir.lower_ops([(g, q)])
                if g == 'erase':
                    body = [s.replace('c[0]', 'c[%d]' % nb) for s in body]
                    nb += 1
                cur += body
            elif g == 'mx':
                t, k = q
                if k != nb:
                    raise LookupBuildError('classical bit order')
                nb += 1
                cur += ['h q[%d];' % t, 'c[%d] = measure q[%d];' % (k, t), 'if (c[%d]) { x q[%d]; }' % (k, t)]
            elif g == 'ifp':
                ks, gates = q
                segs.append((True, cur))
                cur = []
                cond = ' ^ '.join('c[%d]' % k for k in ks)
                segs.append((False, ['if (%s) { %s }' % (cond, ' '.join(_ir._gtxt(gg, qq) for gg, qq in gates))]))
            else:
                cur.append(_ir._gtxt(g, q))
        segs.append((True, cur))
        for ok, s in segs:
            if peephole and ok and s:
                s = _ir.peephole(s)
            L += s
        N = max(self.N, n_qubits or 0)
        head = ['OPENQASM 3;', 'include "stdgates.inc";', 'qubit[%d] q;' % N]
        if nb:
            head.append('bit[%d] c;' % nb)
        return '\n'.join(head + L) + '\n'


# ================================================================================================ cost formulas
def davio_ands(c):
    """ANDs of an ANF/Davio lookup with a full core (= Gidney 2025's lookup count 2^c - c - 1)."""
    return max(0, (1 << c) - c - 1)


def base_ands(h, l):
    """ANDs for all monomials of degree >= 2 of two address halves of h and l bits."""
    return davio_ands(h) + davio_ands(l)


def unlookup_ands(c):
    """U(c): ANDs of the measurement-based unlookup of a c-bit-address table (generic table)."""
    h = (c + 1) // 2
    return base_ands(h, c - h)


def qroam1_ands(c, w):
    """ANDs of qroam1_lookup for a generic table: (2^(c-1) - c) + w + U(c-1)."""
    return ((1 << (c - 1)) - c) + w + unlookup_ands(c - 1)


# ================================================================================================ ANF and GF(2) helpers
def anf_table(T, c):
    """Moebius transform: coefficient (as a w-bit int) of each monomial (bit mask of address bits)."""
    a = list(T)
    for i in range(c):
        bit = 1 << i
        for x in range(1 << c):
            if x & bit:
                a[x] ^= a[x ^ bit]
    return a


class _Basis:
    def __init__(self):
        self.piv = {}

    def red(self, v):
        while v:
            p = v.bit_length() - 1
            if p in self.piv:
                v ^= self.piv[p]
            else:
                return v
        return 0

    def add(self, v):
        v = self.red(v)
        if v:
            self.piv[v.bit_length() - 1] = v
            return True
        return False

    def dim(self):
        return len(self.piv)


def _cores(T, h, l, w):
    """core matrix of each output bit k: bit (S-1)*B + (U-1) = ANF coefficient of M_S N_U (S, U != 0)."""
    anf = anf_table(T, h + l)
    B = (1 << l) - 1
    cores = []
    for k in range(w):
        v = 0
        for S in range(1, 1 << h):
            for U in range(1, 1 << l):
                if (anf[(S << l) | U] >> k) & 1:
                    v |= 1 << ((S - 1) * B + (U - 1))
        cores.append(v)
    return cores


def _in_V(V, h, l, limit, rng):
    """mixed elements (X, alpha, beta, alpha', beta') with core X in V (kernels of linear systems, see module doc)."""
    A = (1 << h) - 1
    B = (1 << l) - 1
    pairs = [(a, a2) for a in range(1, 1 << A) for a2 in range(0, 1 << A)]
    rng.shuffle(pairs)
    out = []
    for a, a2 in pairs:
        cols = []
        for src in (a2, a):            # beta_j -> entries (i, j), i in alpha';  beta'_j -> entries (i, j), i in alpha
            for j in range(B):
                v = 0
                for i in range(A):
                    if (src >> i) & 1:
                        v |= 1 << (i * B + j)
                cols.append(v)
        piv = {}
        kern = []
        for idx, cv in enumerate(cols):
            r = V.red(cv)
            comb = 1 << idx
            while r:
                p = r.bit_length() - 1
                if p in piv:
                    r ^= piv[p][0]
                    comb ^= piv[p][1]
                else:
                    piv[p] = (r, comb)
                    break
            if r == 0:
                kern.append(comb)
        for comb in kern:
            beta = comb & ((1 << B) - 1)
            beta2 = comb >> B
            X = 0
            for idx in range(2 * B):
                if (comb >> idx) & 1:
                    X ^= cols[idx]
            if X == 0 or (a2 == 0 and beta == 0):
                continue
            out.append((X, a, beta, a2, beta2))
        if len(out) >= limit:
            break
    return out


def _random_mixed(h, l, rng):
    A = (1 << h) - 1
    B = (1 << l) - 1
    while True:
        a, a2 = rng.randrange(1, 1 << A), rng.randrange(0, 1 << A)
        b, b2 = rng.randrange(0, 1 << B), rng.randrange(1, 1 << B)
        X = 0
        for i in range(A):
            for j in range(B):
                if (((a >> i) & (b2 >> j)) ^ ((a2 >> i) & (b >> j))) & 1:
                    X |= 1 << (i * B + j)
        if X:
            return (X, a, b, a2, b2)


def _prune(cores, sel, rng):
    order = list(range(len(sel)))
    rng.shuffle(order)
    keep = set(order)
    for i in order:
        rest = keep - {i}
        Bs = _Basis()
        for j in rest:
            Bs.add(sel[j][0])
        if all(Bs.red(c) == 0 for c in cores):
            keep = rest
    return [sel[j] for j in sorted(keep)]


def _decompose(cores, h, l, rng, limit=10 ** 9, max_helpers=12):
    """-> list of mixed elements whose cores span all cores (helpers outside V if needed), or None."""
    U = _Basis()
    for c in cores:
        U.add(c)
    helpers = []
    for _ in range(max_helpers + 1):
        els = _in_V(U, h, l, limit, rng) + helpers
        S = _Basis()
        sel = []
        for e in els:
            if S.add(e[0]):
                sel.append(e)
        if all(S.red(c) == 0 for c in cores):
            return _prune(cores, sel, rng)
        cand = [e for e in (_random_mixed(h, l, rng) for _ in range(64)) if U.red(e[0])]
        if not cand:
            return None
        helpers.append(cand[0])
        U.add(cand[0][0])
    return None


def _coeffs(cores, vecs):
    """lam[i] (w-bit int): bit k set iff vecs[i] is used in the expansion of core k."""
    piv = {}
    for idx, v in enumerate(vecs):
        comb = 1 << idx
        while v:
            p = v.bit_length() - 1
            if p in piv:
                v ^= piv[p][0]
                comb ^= piv[p][1]
            else:
                piv[p] = (v, comb)
                break
    lam = [0] * len(vecs)
    for k, cm in enumerate(cores):
        v = cm
        comb = 0
        while v:
            p = v.bit_length() - 1
            if p not in piv:
                raise LookupBuildError('core outside the span of the decomposition')
            v ^= piv[p][0]
            comb ^= piv[p][1]
        for idx in range(len(vecs)):
            if (comb >> idx) & 1:
                lam[idx] |= 1 << k
    return lam


@dataclass
class Plan:
    """Bilinear plan of a narrow lookup: address split (hi = top h bits, lo = low l bits) and the mixed products
    (alpha, beta, alpha', beta', lam): f = alpha.M + beta.N, g = alpha'.M + beta'.N, lam = output bits fed by f*g."""
    c: int
    w: int
    h: int
    l: int
    products: list
    dim_core: int
    seed: int
    seconds: float

    @property
    def base(self):
        return base_ands(self.h, self.l)

    @property
    def R(self):
        return len(self.products)

    @property
    def n_and(self):
        return self.base + self.R


def _check_table(table, c, w):
    if not isinstance(c, int) or c < 1:
        raise ValueError('address width c must be a positive int')
    if not isinstance(w, int) or w < 1:
        raise ValueError('word width w must be a positive int')
    T = [int(v) for v in table]
    if len(T) != 1 << c:
        raise ValueError('table must have 2^c = %d entries, got %d' % (1 << c, len(T)))
    if any(v < 0 or v >> w for v in T):
        raise ValueError('table entries must be in [0, 2^w)')
    return T


def plan_narrow(table, c, w, seed=0, splits=None, tries=None):
    """Classical search for the bilinear plan (pure Python, single thread).
    splits: candidate h (hi bits); default 1..min(3, c//2) (h > l is the transposed problem, h > 3 is too slow).
    tries: randomised searches per split, best kept, early stop at R = dim V; default 3 for h <= 2 (fast) and 1 for
    h = 3 (each try ~0.3-2 s at c=6).  A split whose lower bound base + dim V cannot beat the best plan is skipped.
    Measured search time per table (1 core): c <= 5: < 0.01 s; c=6: w=32 0.3 s (max 1.4 s), w=16..21 ~2-3 s (the
    h=3 split is tried and loses); c=7, w=32: ~25 s.  Returns the cheapest Plan, or None (c < 2 / nothing found)."""
    T = _check_table(table, c, w)
    if c < 2:
        return None
    if splits is None:
        splits = [h for h in (1, 2, 3) if h <= c - h]
    t_all = time.time()
    best = None
    for h in splits:
        l = c - h
        if h < 1 or l < 1:
            raise ValueError('bad split h=%d for c=%d' % (h, c))
        t0 = time.time()
        cores = _cores(T, h, l, w)
        V = _Basis()
        for cm in cores:
            V.add(cm)
        if best is not None and base_ands(h, l) + V.dim() >= best.n_and:
            continue                        # this split cannot beat the best plan found so far
        if V.dim() == 0:
            p = Plan(c, w, h, l, [], 0, seed, time.time() - t0)
        else:
            sel = None
            nt = tries if tries is not None else (3 if h <= 2 else 1)
            for t in range(nt):             # best of nt randomised searches; R = dim V is optimal: stop
                s = _decompose(cores, h, l, random.Random(seed * 1000 + t))
                if s is not None and (sel is None or len(s) < len(sel)):
                    sel = s
                if sel is not None and len(sel) == V.dim():
                    break
            if sel is None:
                continue
            lam = _coeffs(cores, [e[0] for e in sel])
            prods = [(e[1], e[2], e[3], e[4], lv) for e, lv in zip(sel, lam) if lv]
            p = Plan(c, w, h, l, prods, V.dim(), seed, time.time() - t0)
        if best is None or p.n_and < best.n_and:
            best = p
    if best is not None:
        best.seconds = time.time() - t_all          # total search time over all splits
    return best


# ================================================================================================ emitters (on a Circ)
def _fan(k, src, vec, z):
    """z[i] ^= src for each set bit i of vec (src None: X)."""
    i = 0
    while vec:
        if vec & 1:
            if src is None:
                k.x(z[i])
            else:
                k.cx(src, z[i])
        vec >>= 1
        i += 1


def _build_mono(k, bits, W):
    """monomial wires of `bits`: W[mask] (degree 1 = the bit itself, degree >= 2 by ANDs). Returns the made list."""
    made = []
    for m in range(1, 1 << len(bits)):
        if m & (m - 1) == 0:
            W[m] = bits[m.bit_length() - 1]
    for m in range(1, 1 << len(bits)):
        if m & (m - 1):
            top = m.bit_length() - 1
            par = m ^ (1 << top)
            W[m] = k.AND(W[par], bits[top])
            made.append((m, par, top))
    return made


def _erase_mono(k, bits, W, made):
    for m, par, top in reversed(made):
        k.erase_and(W[m], W[par], bits[top])


def _need(anf, c):
    need = [False] * (1 << c)
    for S in range((1 << c) - 1, -1, -1):
        if anf[S]:
            need[S] = True
        else:
            for j in range(S.bit_length(), c):
                if need[S | (1 << j)]:
                    need[S] = True
                    break
    return need


def _davio_count(T, c):
    need = _need(anf_table(T, c), c)
    return sum(1 for S in range(1 << c) if need[S] and bin(S).count('1') >= 2)


def emit_davio(k, addr, table, z):
    """z ^= table[addr] by the ANF subset tree (one AND per needed monomial of degree >= 2)."""
    c = len(addr)
    anf = anf_table(table, c)
    need = _need(anf, c)

    def visit(S, wv):
        _fan(k, wv, anf[S], z)
        for j in range(S.bit_length(), c):
            ch = S | (1 << j)
            if not need[ch]:
                continue
            if S == 0:
                visit(ch, addr[j])
            else:
                t = k.AND(wv, addr[j])
                visit(ch, t)
                k.erase_and(t, wv, addr[j])
    visit(0, None)


def gray_fanout_cx(naive, nz, b):
    """CX of a Gray-code fan-out batch: 2^b - 1 (accumulator walk) + (b - 1 at most for the undo, here exactly the
    weight of the last Gray word, 1) + nz (one CX per output row with a nonzero pattern); `naive` = direct CX count.
    Returned: the cheaper of the two, which is what _gray_fanout emits."""
    return min(naive, (1 << b) - 1 + 1 + nz)


def _gray_fanout(k, wires, vecs, z, force=False):
    """z[i] ^= XOR_{j: bit i of vecs[j]} wires[j] for a batch of b = len(wires) simultaneously live wires.
    Direct CX per (row, wire), or ("Russian fours", Four Russians) ONE extra accumulator wire walked along the reflected
    Gray code over all 2^b - 1 non-empty subsets of the batch: row i with pattern p (the set of wires it needs) takes one CX
    from the accumulator at the moment it holds exactly p.  CX = (2^b - 1) + 1 + #rows with p != 0 instead of
    sum_i |p_i| (about b w / 2 for a dense random table); the cheaper is chosen unless force.  No AND, T = 0; +1 qubit
    (released at the end; the accumulator returns to 0)."""
    b = len(wires)
    w = len(z)
    pat = [0] * w
    for j, v in enumerate(vecs):
        i = 0
        while v:
            if v & 1:
                pat[i] |= 1 << j
            v >>= 1
            i += 1
    naive = sum(bin(p).count('1') for p in pat)
    nz = sum(1 for p in pat if p)
    if not force and naive <= (1 << b) + nz:
        for i, p in enumerate(pat):
            for j in range(b):
                if (p >> j) & 1:
                    k.cx(wires[j], z[i])
        return
    bucket = {}
    for i, p in enumerate(pat):
        if p:
            bucket.setdefault(p, []).append(i)
    acc = k.new()
    cur = 0
    for t in range(1, 1 << b):
        g = t ^ (t >> 1)
        j = (g ^ cur).bit_length() - 1
        k.cx(wires[j], acc)
        cur = g
        for i in bucket.get(g, ()):
            k.cx(acc, z[i])
    for j in range(b):
        if (cur >> j) & 1:
            k.cx(wires[j], acc)
    k.rel(acc)


def emit_davio_gray(k, addr, table, z, l=None, B=None, force=False):
    """z ^= table[addr] like emit_davio (same AND count for a dense table: one AND per needed monomial of degree
    >= 2) but with the CX fan-out done in batches by the Gray-code accumulator (_gray_fanout).  The address is split
    into lo (l bits, default c // 2: ALL lo monomials stay alive, 2^l - l - 1 ANDs) and hi (h = c - l bits, monomials by
    a DFS); for each hi monomial M_S the products M_S N_T (one AND each) are made in batches of B wires that are alive
    together and fanned out; batch size B (default: from the word width, 2..8) = 2^B + w CX per batch. Pays off for wide words
    (w >> 2^B / B); for narrow words the direct fan-out is chosen automatically batch by batch, so the CX count is never
    larger than that of the same bilinear order with direct CX.  force=True always walks the Gray code (tests)."""
    c = len(addr)
    if l is None:
        l = c // 2
    if B is None:
        B = max(2, min(8, len(z).bit_length() - 2))          # w = 32: 4, 64: 5, 256: 7, >= 1024: 8
    if not 0 <= l <= c or B < 1:
        raise LookupBuildError('bad split / batch size')
    h = c - l
    anf = anf_table(table, c)
    lo, hi = addr[:l], addr[l:]
    Nw = {}
    madeN = _build_mono(k, lo, Nw)
    hneed = [False] * (1 << h)
    for S in range((1 << h) - 1, -1, -1):
        if any(anf[(S << l) | Tm] for Tm in range(1 << l)):
            hneed[S] = True
        else:
            for j in range(S.bit_length(), h):
                if hneed[S | (1 << j)]:
                    hneed[S] = True
                    break

    def node(S, mS):
        items = [(Tm, anf[(S << l) | Tm]) for Tm in range(1 << l) if anf[(S << l) | Tm]]
        if S == 0:
            for Tm, v in items:
                if Tm == 0:
                    _fan(k, None, v, z)
            items = [(Tm, v) for Tm, v in items if Tm]
            for st in range(0, len(items), B):
                bt = items[st:st + B]
                _gray_fanout(k, [Nw[Tm] for Tm, _ in bt], [v for _, v in bt], z, force)
            return
        for st in range(0, len(items), B):
            bt = items[st:st + B]
            ws = []
            for Tm, v in bt:
                if Tm == 0:
                    ws.append(mS)
                else:
                    ws.append(k.AND(mS, Nw[Tm]))
            _gray_fanout(k, ws, [v for _, v in bt], z, force)
            for (Tm, v), wv in zip(bt, ws):
                if Tm:
                    k.erase_and(wv, mS, Nw[Tm])

    def visit(S, wv):
        node(S, wv)
        for j in range(S.bit_length(), h):
            ch = S | (1 << j)
            if not hneed[ch]:
                continue
            if S == 0:
                visit(ch, hi[j])
            else:
                t = k.AND(wv, hi[j])
                visit(ch, t)
                k.erase_and(t, wv, hi[j])
    visit(0, None)
    _erase_mono(k, lo, Nw, madeN)


def emit_narrow(k, addr, table, z, plan):
    """z ^= table[addr] by the bilinear plan (plan_narrow).  Returns the number of mixed products (ANDs beyond base).
    The free part (pure hi / pure lo monomials) is computed from truth tables and the absence of a mixed residue is
    asserted (LookupBuildError otherwise)."""
    c = len(addr)
    h, l = plan.h, plan.l
    if h + l != c:
        raise LookupBuildError('plan for another address width')
    lo, hi = addr[:l], addr[l:]
    Mw, Nw = {}, {}
    madeM = _build_mono(k, hi, Mw)
    madeN = _build_mono(k, lo, Nw)
    A = (1 << h) - 1
    B = (1 << l) - 1
    ttM = {}
    ttN = {}
    for S in range(1, 1 << h):
        ttM[S] = sum(1 << a for a in range(1 << c) if ((a >> l) & S) == S)
    for U in range(1, 1 << l):
        ttN[U] = sum(1 << a for a in range(1 << c) if (a & U) == U)
    L = list(table)
    prods = []
    for a, b, a2, b2, lv in plan.products:
        sf = [('M', s + 1) for s in range(A) if (a >> s) & 1] + [('N', t + 1) for t in range(B) if (b >> t) & 1]
        sg = [('M', s + 1) for s in range(A) if (a2 >> s) & 1] + [('N', t + 1) for t in range(B) if (b2 >> t) & 1]
        if set(sf) <= set(sg):
            sg = sorted(set(sf) ^ set(sg))          # f*g = f*(f+g) + f   (f is free)
        elif set(sg) <= set(sf):
            sf = sorted(set(sf) ^ set(sg))
        if not sf or not sg:
            raise LookupBuildError('degenerate product')
        tf = 0
        for kind, idx in sf:
            tf ^= ttM[idx] if kind == 'M' else ttN[idx]
        tg = 0
        for kind, idx in sg:
            tg ^= ttM[idx] if kind == 'M' else ttN[idx]
        P = tf & tg
        for x in range(1 << c):
            if (P >> x) & 1:
                L[x] ^= lv
        prods.append((sf, sg, lv))
    anf = anf_table(L, c)
    for S in range(1, 1 << h):
        for U in range(1, 1 << l):
            if anf[(S << l) | U]:
                raise LookupBuildError('mixed-monomial residue: decomposition wrong')

    def wire(kind, idx):
        return Mw[idx] if kind == 'M' else Nw[idx]
    for U in range(1 << l):
        _fan(k, Nw[U] if U else None, anf[U], z)
    for S in range(1, 1 << h):
        _fan(k, Mw[S], anf[S << l], z)
    for sf, sg, lv in prods:
        pf = next(x for x in sf if x not in sg)
        pg = next(x for x in sg if x not in sf)
        fw, gw = wire(*pf), wire(*pg)
        of = [wire(*x) for x in sf if x != pf]
        og = [wire(*x) for x in sg if x != pg]
        if fw in og or gw in of:
            raise LookupBuildError('pivot clash')
        for wv in of:
            k.cx(wv, fw)
        for wv in og:
            k.cx(wv, gw)
        p = k.AND(fw, gw)
        _fan(k, p, lv, z)
        k.erase_and(p, fw, gw)
        for wv in reversed(og):
            k.cx(wv, gw)
        for wv in reversed(of):
            k.cx(wv, fw)
    _erase_mono(k, lo, Nw, madeN)
    _erase_mono(k, hi, Mw, madeM)
    return len(prods)


def emit_unlookup_meas(k, addr, table, z, h=None):
    """z holds table[addr] (plus nothing else): X-measure every z wire (outcomes m), then fix the phase
    (-1)^{m.table[addr]} with conditional Z / CZ on products of hi x lo monomials (U(c) ANDs).  z ends clean."""
    c = len(addr)
    if h is None:
        h = (c + 1) // 2
    l = c - h
    anf = anf_table(table, c)
    ks = [k.mx(wv) for wv in z]
    lo, hi = addr[:l], addr[l:]
    Nw = {0: None}
    made = _build_mono(k, lo, Nw)

    def cond(mono):
        v = anf[mono]
        return [ks[i] for i in range(len(z)) if (v >> i) & 1]

    def apply(S, ws):
        for U in range(1 << l):
            kk = cond((S << l) | U)
            if not kk:
                continue
            wt = Nw[U]
            if ws is None and wt is None:
                continue                      # outcome-dependent global phase
            if ws is None:
                k.ifpar(kk, [('z', (wt,))])
            elif wt is None:
                k.ifpar(kk, [('z', (ws,))])
            else:
                k.ifpar(kk, [('cz', (ws, wt))])
    hneed = [False] * (1 << h)
    for S in range((1 << h) - 1, -1, -1):
        if any(anf[(S << l) | U] for U in range(1 << l)):
            hneed[S] = True
        else:
            for j in range(S.bit_length(), h):
                if hneed[S | (1 << j)]:
                    hneed[S] = True
                    break

    def visit(S, wv):
        apply(S, wv)
        for j in range(S.bit_length(), h):
            ch = S | (1 << j)
            if not hneed[ch]:
                continue
            if S == 0:
                visit(ch, hi[j])
            else:
                t = k.AND(wv, hi[j])
                visit(ch, t)
                k.erase_and(t, wv, hi[j])
    visit(0, None)
    _erase_mono(k, lo, Nw, made)


def emit_qroam1(k, addr, table, z):
    """z ^= table[addr]: QROAM with one select bit (Berry et al. 2019 / Low et al. 2018, lambda = 2)."""
    c = len(addr)
    w = len(z)
    if c < 2:
        raise ValueError('qroam1 needs c >= 2')
    a0, hi = addr[0], addr[1:]
    R = [k.new() for _ in range(w)]
    H = 1 << (c - 1)
    T2 = [table[2 * j] | ((table[2 * j] ^ table[2 * j + 1]) << w) for j in range(H)]
    emit_davio(k, hi, T2, list(z) + R)
    for i in range(w):
        t = k.AND(a0, R[i])
        k.cx(t, z[i])
        k.erase_and(t, a0, R[i])
    emit_unlookup_meas(k, hi, [table[2 * j] ^ table[2 * j + 1] for j in range(H)], R)


def emit_add_gidney(k, Y, Z):
    """Y += Z mod 2^len(Y) (Gidney 1709.06648 adder; Z may be shorter: zero high bits); len(Y)-1 ANDs."""
    m = len(Y)
    s = len(Z)
    if m == 0 or s == 0:
        return
    if m == 1:
        k.cx(Z[0], Y[0])
        return
    cr = {1: k.AND(Y[0], Z[0])}
    for i in range(1, m - 1):
        if i < s:
            k.cx(cr[i], Y[i])
            k.cx(cr[i], Z[i])
            t = k.AND(Y[i], Z[i])
            k.cx(cr[i], t)
        else:
            t = k.AND(Y[i], cr[i])
        cr[i + 1] = t
    i = m - 1
    if i < s:
        k.cx(Z[i], Y[i])
    k.cx(cr[i], Y[i])
    for i in range(m - 2, 0, -1):
        if i < s:
            k.erase_and(cr[i + 1], Y[i], Z[i], [cr[i]])
            k.cx(cr[i], Z[i])
            k.cx(Z[i], Y[i])
        else:
            k.erase_and(cr[i + 1], Y[i], cr[i])
            k.cx(cr[i], Y[i])
    k.erase_and(cr[1], Y[0], Z[0])
    k.cx(Z[0], Y[0])


# ================================================================================================ references
def _ref_lookup(g, addr, T, z, chain):
    """reference Toffoli network: z ^= T[addr] by enumerating addresses (AND chain of literals into `chain`)."""
    c = len(addr)
    for a in range(1 << c):
        bits = [i for i in range(len(z)) if (T[a] >> i) & 1]
        if not bits:
            continue
        neg = [addr[i] for i in range(c) if not (a >> i) & 1]
        for q in neg:
            g.append('tof ' + q)
        if c == 1:
            for i in bits:
                g.append('tof %s %s' % (addr[0], z[i]))
        else:
            g.append('tof %s %s %s' % (addr[0], addr[1], chain[0]))
            for i in range(2, c):
                g.append('tof %s %s %s' % (chain[i - 2], addr[i], chain[i - 1]))
            for i in bits:
                g.append('tof %s %s' % (chain[c - 2], z[i]))
            for i in range(c - 1, 1, -1):
                g.append('tof %s %s %s' % (chain[i - 2], addr[i], chain[i - 1]))
            g.append('tof %s %s %s' % (addr[0], addr[1], chain[0]))
        for q in neg:
            g.append('tof ' + q)


def _ref_cuccaro(g, Y, Z, anc):
    """reference adder Y += Z mod 2^m (Cuccaro et al. 2004 with Toffolis), Z restored, anc clean."""
    m = len(Y)
    if m == 1:
        g.append('tof %s %s' % (Z[0], Y[0]))
        return
    xs = [anc] + list(Z[:m - 1])
    for i in range(m - 1):
        g += ['tof %s %s' % (Z[i], Y[i]), 'tof %s %s' % (Z[i], xs[i]), 'tof %s %s %s' % (xs[i], Y[i], Z[i])]
    g += ['tof %s %s' % (Z[m - 1], Y[m - 1]), 'tof %s %s' % (xs[m - 1], Y[m - 1])]
    for i in range(m - 2, -1, -1):
        g += ['tof %s %s %s' % (xs[i], Y[i], Z[i]), 'tof %s %s' % (Z[i], xs[i]), 'tof %s %s' % (xs[i], Y[i])]


def reference_qc(table, c, w, kind='xor', n_wires=None):
    """Independent reference (.qc text, Toffoli network) for the circuits of this module.
    kind 'xor': wires a0..a(c-1), z0..z(w-1) (both inputs), z ^= T[a]   (lookups; lookup_pair: y = the z wires)
    kind 'add': wires a.., y0..y(w-1) inputs, y += T[a] mod 2^w           (lookup_add)
    then clean wires e0.. (n_wires in total, at least what the reference needs)."""
    T = _check_table(table, c, w)
    names = ['a%d' % i for i in range(c)] + ['z%d' % i for i in range(w)]
    ins = list(names)
    need = c + w + max(c - 1, 1) + (w + 1 if kind == 'add' else 0)
    n = max(need, n_wires or 0)
    names += ['e%d' % i for i in range(n - c - w)]
    addr, z, ex = names[:c], names[c:c + w], names[c + w:]
    g = []
    if kind == 'xor':
        _ref_lookup(g, addr, T, z, ex[:max(c - 1, 1)])
    elif kind == 'add':
        t = ex[:w]
        chain = ex[w:w + max(c - 1, 1)]
        anc = ex[w + max(c - 1, 1)]
        _ref_lookup(g, addr, T, t, chain)
        _ref_cuccaro(g, z, t, anc)
        _ref_lookup(g, addr, T, t, chain)
    else:
        raise ValueError(kind)
    return '.v %s\n.i %s\n\nBEGIN\n%s\nEND\n' % (' '.join(names), ' '.join(ins), '\n'.join(g))


# ================================================================================================ public builders
@dataclass
class LookupCircuit:
    """A checked circuit.  Wires 0..c-1: address (unchanged), c..c+w-1: z (lookups: z ^= T[a]) or y (lookup_add:
    y += T[a] mod 2^w); all other qubits start and end in |0>.  qasm declares max(qubits, reference wires) qubits."""
    kind: str
    method: str
    c: int
    w: int
    table: tuple
    qasm: str
    n_and: int
    t_count: int
    cx_count: int
    qubits: int
    ancillas_peak: int
    plan: object = None
    plan_seconds: float = 0.0
    checked: str = ''
    ref_kind: str = 'xor'
    extra: dict = field(default_factory=dict)

    def reference_qc(self):
        n = int(self.qasm.split('qubit[', 1)[1].split(']', 1)[0])
        return reference_qc(self.table, self.c, self.w, self.ref_kind, n_wires=n)

    def write(self, qasm_path, qc_path=None):
        with open(qasm_path, 'w') as fh:
            fh.write(self.qasm)
        if qc_path:
            with open(qc_path, 'w') as fh:
                fh.write(self.reference_qc())


def _z_samples(w, n, rng):
    full = (1 << w) - 1
    base = [0, full, 0x5555555555555555 & full, 1 << (w - 1)]
    out = base[:n]
    while len(out) < n:
        out.append(rng.getrandbits(w))
    return out


def _assemble(kind, method, T, c, w, emit, ref_kind, z_samples, outcome_patterns, seed, plan=None, pt=0.0,
              mode='xor'):
    """build `emit(k, addr, zw)` on all 2^c addresses x z_samples values (x outcome_patterns seeds), check, export."""
    rng = random.Random(1000003 * seed + 17)
    zs = _z_samples(w, max(1, z_samples), rng)
    A = np.array([a for a in range(1 << c) for _ in zs], np.int64)
    Z = [zv for _ in range(1 << c) for zv in zs]
    Zbits = np.array([[(zv >> i) & 1 for i in range(w)] for zv in Z], np.uint8)
    if mode == 'xor':
        want = [zv ^ T[a] for a, zv in zip(A.tolist(), Z)]
    else:
        want = [(zv + T[a]) % (1 << w) for a, zv in zip(A.tolist(), Z)]
    Wb = np.array([[(v >> i) & 1 for i in range(w)] for v in want], np.uint8)
    k = None
    for ms in range(max(1, outcome_patterns)):
        k = Circ(len(A), seed=ms + 1)
        aw = [k.inp(((A >> i) & 1).astype(np.uint8)) for i in range(c)]
        zw = [k.inp(Zbits[:, i]) for i in range(w)]
        emit(k, aw, zw)
        k.finish({zw[i]: Wb[:, i] for i in range(w)})
    need = c + w + max(c - 1, 1) + (w + 1 if ref_kind == 'add' else 0)
    txt = k.qasm(n_qubits=need)
    T_, CX_ = _ir.count(txt)
    cn = k.counts()
    return LookupCircuit(kind=kind, method=method, c=c, w=w, table=tuple(T), qasm=txt, n_and=cn['AND'], t_count=T_,
                         cx_count=CX_, qubits=cn['qubits'], ancillas_peak=cn['peak'] - (c + w), plan=plan,
                         plan_seconds=pt, ref_kind=ref_kind,
                         checked='exact on all %d addresses x %d values x %d outcome patterns'
                                 % (1 << c, len(zs), max(1, outcome_patterns)))


def _fanout_mode(fanout):
    if fanout not in ('plain', 'gray'):
        raise ValueError("fanout must be 'plain' or 'gray'")
    return fanout == 'gray'


def davio_lookup(table, c, w, z_samples=4, outcome_patterns=3, seed=0, fanout='plain', gray_b=None, gray_l=None,
                 gray_force=False):
    """z ^= T[a] by the ANF/Davio subset tree (baseline). fanout='gray': the same ANDs (for a dense table), CX fan-out
    by the Gray-code accumulator in batches of gray_b wires (None: automatic, 2..8) (emit_davio_gray; method 'davio-gray'); default 'plain'
    reproduces the previous circuits exactly."""
    T = _check_table(table, c, w)
    if _fanout_mode(fanout):
        return _assemble('lookup', 'davio-gray', T, c, w,
                         lambda k, a, z: emit_davio_gray(k, a, T, z, l=gray_l, B=gray_b, force=gray_force), 'xor',
                         z_samples, outcome_patterns, seed)
    return _assemble('lookup', 'davio', T, c, w, lambda k, a, z: emit_davio(k, a, T, z), 'xor', z_samples,
                     outcome_patterns, seed)


def narrow_lookup(table, c, w, seed=0, splits=None, tries=None, z_samples=4, outcome_patterns=3):
    """z ^= T[a] with the bilinear split (see the module docstring); falls back to Davio when that is not cheaper.
    Returns a LookupCircuit (method 'bilinear' or 'davio'; plan, plan_seconds = classical search time)."""
    T = _check_table(table, c, w)
    t0 = time.time()
    plan = plan_narrow(T, c, w, seed=seed, splits=splits, tries=tries)
    pt = time.time() - t0
    dv = _davio_count(T, c)
    if plan is None or plan.n_and >= dv:
        r = _assemble('lookup', 'davio', T, c, w, lambda k, a, z: emit_davio(k, a, T, z), 'xor', z_samples,
                      outcome_patterns, seed, plan=plan, pt=pt)
    else:
        r = _assemble('lookup', 'bilinear', T, c, w, lambda k, a, z: emit_narrow(k, a, T, z, plan), 'xor',
                      z_samples, outcome_patterns, seed, plan=plan, pt=pt)
        if r.n_and != plan.n_and:
            raise LookupBuildError('AND count %d != plan %d' % (r.n_and, plan.n_and))
    r.extra['davio_ands'] = dv
    return r


def qroam1_lookup(table, c, w, z_samples=4, outcome_patterns=3, seed=0):
    """z ^= T[a] by QROAM with one select bit, R cleared by measurement-based unlookup."""
    T = _check_table(table, c, w)
    if c < 2:
        raise ValueError('qroam1 needs c >= 2')
    return _assemble('lookup', 'qroam1', T, c, w, lambda k, a, z: emit_qroam1(k, a, T, z), 'xor', z_samples,
                     outcome_patterns, seed)


def lookup_pair(table, c, w, mode='narrow', seed=0, splits=None, z_samples=4, outcome_patterns=3):
    """y ^= T[a] through a temporary register z: lookup (mode 'narrow' | 'davio' | 'qroam1'), y ^= z, and the
    measurement-based unlookup (for 'qroam1': one unlookup of the combined table [T[a] | T[a&~1] ^ T[a|1]] clears
    z and R together).  Wires: a, y (= the 'z' wires of reference_qc kind 'xor')."""
    T = _check_table(table, c, w)
    plan = None
    pt = 0.0
    method = mode
    if mode == 'narrow':
        t0 = time.time()
        plan = plan_narrow(T, c, w, seed=seed, splits=splits)
        pt = time.time() - t0
        if plan is None or plan.n_and >= _davio_count(T, c):
            method = 'davio'
        else:
            method = 'bilinear'
    elif mode not in ('davio', 'qroam1'):
        raise ValueError(mode)
    if mode == 'qroam1' and c < 2:
        raise ValueError('qroam1 needs c >= 2')

    def emit(k, a, y):
        z = [k.new() for _ in range(w)]
        if method == 'qroam1':
            a0, hi = a[0], a[1:]
            R = [k.new() for _ in range(w)]
            H = 1 << (c - 1)
            emit_davio(k, hi, [T[2 * j] | ((T[2 * j] ^ T[2 * j + 1]) << w) for j in range(H)], list(z) + R)
            for i in range(w):
                t = k.AND(a0, R[i])
                k.cx(t, z[i])
                k.erase_and(t, a0, R[i])
            for i in range(w):
                k.cx(z[i], y[i])
            Tc = [T[x] | ((T[x & ~1] ^ T[x | 1]) << w) for x in range(1 << c)]
            emit_unlookup_meas(k, a, Tc, list(z) + R)
            return
        if method == 'bilinear':
            emit_narrow(k, a, T, z, plan)
        else:
            emit_davio(k, a, T, z)
        for i in range(w):
            k.cx(z[i], y[i])
        emit_unlookup_meas(k, a, T, z)
    return _assemble('pair', method, T, c, w, emit, 'xor', z_samples, outcome_patterns, seed, plan=plan, pt=pt)


def lookup_add(table, c, w, seed=0, splits=None, z_samples=8, outcome_patterns=3, fanout='plain', gray_b=None,
               gray_l=None, gray_force=False):
    """acc(w) += T[a] mod 2^w (one step of Gidney 2025's loop1): narrow lookup into a temporary z, Gidney adder
    (w-1 ANDs), measurement-based unlookup.  Wires: a, acc (= 'z' wires of reference_qc kind 'add').
    fanout='gray' (wide words, e.g. the lookup-additions of a windowed RSA multiplication, w ~ 1000-2000): the lookup is
    the Davio-count lookup with the Gray-code fan-out (emit_davio_gray, method 'davio-gray'; no bilinear plan is
    searched: for w >> 2^c the plan has no advantage, the ANF core is full); default 'plain' is unchanged."""
    T = _check_table(table, c, w)
    if _fanout_mode(fanout):
        def emit_g(k, a, y):
            z = [k.new() for _ in range(w)]
            emit_davio_gray(k, a, T, z, l=gray_l, B=gray_b, force=gray_force)
            emit_add_gidney(k, y, z)
            emit_unlookup_meas(k, a, T, z)
        return _assemble('add', 'davio-gray', T, c, w, emit_g, 'add', z_samples, outcome_patterns, seed, mode='add')
    t0 = time.time()
    plan = plan_narrow(T, c, w, seed=seed, splits=splits)
    pt = time.time() - t0
    method = 'bilinear' if plan is not None and plan.n_and < _davio_count(T, c) else 'davio'

    def emit(k, a, y):
        z = [k.new() for _ in range(w)]
        if method == 'bilinear':
            emit_narrow(k, a, T, z, plan)
        else:
            emit_davio(k, a, T, z)
        emit_add_gidney(k, y, z)
        emit_unlookup_meas(k, a, T, z)
    return _assemble('add', method, T, c, w, emit, 'add', z_samples, outcome_patterns, seed, plan=plan, pt=pt,
                     mode='add')
