# SPDX-License-Identifier: Apache-2.0
"""Exact sum-over-paths with dyadic phases (coefficients in Z_{2^K}): the engine of qoptct.dyadic.

Representation (M. Amy, "Towards Large-scale Functional Verification of Universal Quantum Circuits", QPL 2018,
EPTCS 287 (2019), doi:10.4204/EPTCS.287.1, arXiv:1805.06908v2, Definition 2.1).  An operator acts on basis states as

    |x>  ->  sqrt(2)^e * sum_{y in {0,1}^m} exp(2 pi i F(x, y) / 2^K) |W(x, y)>

* x are free (input) variables, y are path variables (summed over {0, 1});
* F is a multilinear polynomial over the Boolean variables with coefficients in Z_{2^K}.  Amy takes coefficients in
  the dyadic fractions modulo 1 (D/Z); every phase of a circuit with angles pi*p/2^k lies in the subgroup
  2^-K Z / Z, which is Z_{2^K}, so all arithmetic here is exact integer arithmetic (Python integers: no overflow,
  unlike the 31-bit limit reported for the QFT in Amy, Sec. 5.2);
* W holds one Boolean (F2) polynomial per wire in algebraic normal form (ANF);
* e is an integer: the amplitude factor sqrt(2)^e.
F and W are canonical: a function {0,1}^v -> Z_{2^K} has exactly one multilinear representation (Moebius inversion
is an integer formula, so it works over any commutative ring) and a Boolean function has exactly one ANF.
A Boolean polynomial P enters F through its integer lift [P] (Amy, Sec. 2.1): [P xor Q] = [P] + [Q] - 2 [P][Q],
which is exact on every Boolean point; terms whose coefficient becomes 0 mod 2^K are dropped.

Rewrite rules (Amy, Fig. 2 in Sec. 3.2; soundness: Proposition 3.1), generalised from Z_8 to Z_{2^K}.  In every rule
y is a path variable that occurs on no wire and is not frozen; the side condition is checked exactly, and a rule whose
condition does not hold is not applied:

  [Elim]   y does not occur in F:   sum_y g = 2 g.                                                     e += 2
  [HH]     every term of F containing y has coefficient 2^(K-1), so exp(2 pi i F / 2^K) = (-1)^(y Q) exp(..R..)
           where Q is the F2 sum of the quotients (only Q mod 2 matters at coefficient 2^(K-1)).  If Q = z xor Q'
           for a path variable z that occurs in Q only as the linear monomial z (so z is not in Q'), then
               sum_{y, z} (-1)^(y (z xor Q')) g(z) = 2 g(Q'):
           the y terms are dropped and z := Q' is substituted in F (through the lift) and in W.          e += 2
           Q = 1 (the constant) makes the whole sum 0 (flag `zero`; impossible for a unitary circuit).
  [omega]  the linear term y has coefficient c = 2^(K-2) or 3*2^(K-2) and every other term with y has
           coefficient 2^(K-1):  sum_y exp(2 pi i y (c + 2^(K-1) Q) / 2^K) = 1 +- i (-1)^Q
                                 = sqrt(2) exp(+- i pi / 4) exp(-+ i pi Q / 2):
           the y terms are dropped and F += +-2^(K-3) -+ 2^(K-2) [Q].                                     e += 1

Why the HH condition is also necessary for a delta: sum_y exp(2 pi i y G / 2^K) = 1 + exp(2 pi i G / 2^K) is 0 or 2
for every point iff 2 G = 0 mod 2^K pointwise, iff (canonical representation) every coefficient of G is 0 or 2^(K-1).
Each rule removes at least one path variable, so the reduction terminates; nothing here ever answers PROVED by
itself - the callers (qoptct.dyadic.api) only accept a form with no path variables left (canonical) or an
exact enumeration.
"""
import heapq
import time

ONE = frozenset([0])
ZERO = frozenset()

__all__ = ['DPathSum', 'Blowup', 'Timeout', 'ONE', 'ZERO', 'bits', 'pmul', 'peval', 'v2']


class Blowup(Exception):
    """The representation exceeded a size limit (terms of F, size of a lift): the caller answers UNKNOWN."""


class Timeout(Exception):
    """The wall-clock budget was exceeded: the caller answers UNKNOWN."""


def bits(m):
    """Indices of the set bits of m (low to high)."""
    r = []
    while m:
        lo = m & -m
        r.append(lo.bit_length() - 1)
        m ^= lo
    return r


PMUL_CAP = 4000000      # largest number of monomial products in one ANF multiplication (else Blowup)


def pmul(a, b):
    """Product of two ANF polynomials (frozensets of monomial bit masks; 0 is the constant monomial 1)."""
    if len(a) * len(b) > PMUL_CAP:
        raise Blowup('ANF product of %d x %d monomials' % (len(a), len(b)))
    r = set()
    for i in a:
        for j in b:
            k = i | j
            if k in r:
                r.remove(k)
            else:
                r.add(k)
    return frozenset(r)


def peval(P, val):
    """Value (0/1) of the ANF polynomial P at the assignment val (bit v of val = value of variable v)."""
    s = 0
    for mo in P:
        if mo & val == mo:
            s ^= 1
    return s


def v2(c):
    """2-adic valuation of a positive integer."""
    return (c & -c).bit_length() - 1


class DPathSum:
    """Path sum with phases in Z_{2^K}.  Primitive operations (all exact):
         h(q)                 Hadamard: new path variable y, F += 2^(K-1) [W_q] y, W_q := y, e -= 1
         xor(t, cs)           W_t ^= prod_{c in cs} W_c          (X, CX, Toffoli, multi-controlled X)
         swap(a, b), cswap(c, a, b)
         phase(qs, c)         F += c [prod_{q in qs} W_q]        (Z, S, T, P, CP, MCP, CZ, ...; qs = () global)
         phase_xor(qs, c)     F += c [xor_{q in qs} W_q]         (RZZ)
       Variables: kind 'x' (free input) or 'y' (path).  Only 'y' variables are ever summed or substituted."""

    def __init__(self, K, max_terms=400000, deadline=None, max_wire=100000):
        if K < 3:
            raise ValueError('K must be >= 3 (the omega rule needs 2^(K-3))')
        self.K = K
        self.N = 1 << K
        self.M = self.N - 1
        self.nv = 0
        self.kind = []
        self.label = []
        self.ymask = 0          # bit mask of all path variables ever created
        self.ph = {}            # monomial -> coefficient in [1, 2^K)
        self.w = []             # wires: frozenset of monomials (ANF)
        self.e2 = 0
        self.ys = set()         # live path variables
        self.frozen = set()     # path variables excluded from the rules (until unfrozen)
        self.occ = {}           # path variable -> set of monomials of F that contain it
        self.wcnt = {}          # path variable -> number of (wire, monomial) occurrences
        self.dirty = set()      # path variables whose terms or wire occurrences changed
        self.zero = False
        self.max_terms = max_terms
        self.max_wire = max_wire
        self.wtotal = 0         # monomials on all wires together (bounded by max_terms)
        self.deadline = deadline
        self.order = 'min'      # order in which candidate variables are examined: 'min' or 'max' index first
        self.stats = dict(h=0, elim=0, hh=0, omega=0, peak_terms=0)
        self._ticks = 0

    # ------------------------------------------------------------------ variables, limits
    def var(self, kind, label=None):
        v = self.nv
        self.nv += 1
        self.kind.append(kind)
        self.label.append(label)
        if kind == 'y':
            self.ymask |= 1 << v
            self.ys.add(v)
        elif kind != 'x':
            raise ValueError('variable kind must be x or y')
        return v

    def check_limits(self):
        if self.deadline is not None and time.time() > self.deadline:
            raise Timeout('time budget exceeded')
        n = len(self.ph)
        if n > self.stats['peak_terms']:
            self.stats['peak_terms'] = n
        if n > self.max_terms:
            raise Blowup('phase polynomial exceeds %d terms' % self.max_terms)

    def _tick(self):
        self._ticks += 1
        if not self._ticks & 1023:
            self.check_limits()

    # ------------------------------------------------------------------ wires
    def add_wire(self, P):
        self.w.append(ZERO)
        self.set_wire(len(self.w) - 1, P)
        return len(self.w) - 1

    def set_wire(self, q, P):
        old = self.w[q]
        if old == P:
            return
        if len(P) > self.max_wire:
            raise Blowup('wire polynomial exceeds %d monomials' % self.max_wire)
        self.wtotal += len(P) - len(old)
        if self.wtotal > self.max_terms:
            raise Blowup('wire polynomials exceed %d monomials in total' % self.max_terms)
        ym = self.ymask
        for mo in old:
            for v in bits(mo & ym):
                c = self.wcnt[v] - 1
                if c:
                    self.wcnt[v] = c
                else:
                    del self.wcnt[v]
                    self.dirty.add(v)
        for mo in P:
            for v in bits(mo & ym):
                self.wcnt[v] = self.wcnt.get(v, 0) + 1
        self.w[q] = P

    def on_wire(self, v):
        return v in self.wcnt

    # ------------------------------------------------------------------ phase polynomial
    def addc(self, mo, c):
        """F += c * mo  (mod 2^K)."""
        c &= self.M
        if not c:
            return
        ph = self.ph
        old = ph.get(mo)
        ym = mo & self.ymask
        if old is None:
            ph[mo] = c
            if ym:
                for v in bits(ym):
                    s = self.occ.get(v)
                    if s is None:
                        self.occ[v] = {mo}
                    else:
                        s.add(mo)
                    self.dirty.add(v)
            self._tick()
            return
        nc = (old + c) & self.M
        if nc:
            ph[mo] = nc
            if ym:
                self.dirty.update(bits(ym))
        else:
            del ph[mo]
            if ym:
                for v in bits(ym):
                    s = self.occ[v]
                    s.discard(mo)
                    if not s:
                        del self.occ[v]
                    self.dirty.add(v)

    def lift(self, P, kb):
        """Integer multilinear lift [P] of the F2 polynomial P modulo 2^kb -> dict monomial -> coefficient.
        [m1 xor ... xor mL] is built by [A xor m] = [A] + m - 2 [A] m (Amy, Sec. 2.1), exact on Boolean points."""
        if kb <= 0 or not P:
            return {}
        if kb == 1 or len(P) == 1:
            return {mo: 1 for mo in P}
        mm = (1 << kb) - 1
        L = {}
        for m in P:
            upd = {}
            for u, cu in L.items():
                k = u | m
                upd[k] = upd.get(k, 0) - 2 * cu
            upd[m] = upd.get(m, 0) + 1
            for k, d in upd.items():
                nv = (L.get(k, 0) + d) & mm
                if nv:
                    L[k] = nv
                else:
                    L.pop(k, None)
            if len(L) > self.max_terms:
                raise Blowup('lift exceeds %d terms' % self.max_terms)
            self._tick()
        return L

    def add_lift(self, P, c, extra=0):
        """F += c * [P] * extra  (P an F2 polynomial, extra a monomial)."""
        c &= self.M
        if not c or not P:
            return
        if len(P) == 1:
            for mo in P:
                self.addc(mo | extra, c)
            return
        kb = self.K - v2(c)
        for u, cu in self.lift(P, kb).items():
            self.addc(u | extra, c * cu)

    # ------------------------------------------------------------------ primitive gates
    def h(self, q):
        y = self.var('y', 'h')
        yb = 1 << y
        half = self.N >> 1
        for mo in self.w[q]:
            self.addc(mo | yb, half)
        self.set_wire(q, frozenset([yb]))
        self.e2 -= 1
        self.stats['h'] += 1
        return y

    def xor(self, t, cs):
        P = frozenset([0])
        for c in cs:
            P = pmul(P, self.w[c])
        if P:
            self.set_wire(t, self.w[t] ^ P)

    def swap(self, a, b):
        self.w[a], self.w[b] = self.w[b], self.w[a]

    def cswap(self, c, a, b):
        d = pmul(self.w[c], self.w[a] ^ self.w[b])
        if d:
            self.set_wire(a, self.w[a] ^ d)
            self.set_wire(b, self.w[b] ^ d)

    def phase(self, qs, c):
        P = frozenset([0])
        for q in qs:
            P = pmul(P, self.w[q])
        self.add_lift(P, c)

    def phase_xor(self, qs, c):
        P = ZERO
        for q in qs:
            P = P ^ self.w[q]
        self.add_lift(P, c)

    # ------------------------------------------------------------------ substitution z := R
    def subst(self, z, R):
        zb = 1 << z
        items = [(mo, self.ph[mo]) for mo in list(self.occ.get(z, ()))]
        for mo, c in items:
            self.addc(mo, -c)
        for mo, c in items:
            self.add_lift(R, c, mo ^ zb)
        for q, P in enumerate(self.w):
            if any(mo & zb for mo in P):
                r = set()
                for mo in P:
                    if len(r) > self.max_wire:
                        raise Blowup('wire polynomial exceeds %d monomials' % self.max_wire)
                    if mo & zb:
                        for k in pmul(frozenset([mo ^ zb]), R):
                            if k in r:
                                r.remove(k)
                            else:
                                r.add(k)
                    elif mo in r:
                        r.remove(mo)
                    else:
                        r.add(mo)
                self.set_wire(q, frozenset(r))
        if z in self.occ or z in self.wcnt:
            raise AssertionError('substitution left the variable in place')

    # ------------------------------------------------------------------ rules
    def _hh_candidates(self, Q):
        cnt = {}
        for mo in Q:
            for v in bits(mo & self.ymask):
                cnt[v] = cnt.get(v, 0) + 1
        return [v for v, k in cnt.items()
                if k == 1 and (1 << v) in Q and v in self.ys and v not in self.frozen]

    def unfreeze(self, vs=None):
        vs = set(self.frozen) if vs is None else set(vs)
        self.frozen -= vs
        self.dirty |= vs

    def reduce(self):
        """Apply [Elim], [HH], [omega] until none applies to a dirty variable (worklist)."""
        N = self.N
        half, q1 = N >> 1, N >> 2
        q3 = 3 * q1
        sign = 1 if self.order == 'min' else -1
        heap = []
        queued = set()

        def pull():
            d = self.dirty
            self.dirty = set()
            for v in d:
                if v in self.ys and v not in queued:
                    if v in self.frozen:
                        self.dirty.add(v)       # stays dirty for unfreeze()
                        continue
                    heapq.heappush(heap, sign * v)
                    queued.add(v)
        pull()
        while heap:
            y = sign * heapq.heappop(heap)
            queued.discard(y)
            self._tick()
            if y not in self.ys or y in self.frozen or y in self.wcnt:
                continue
            yb = 1 << y
            terms = self.occ.get(y)
            if not terms:
                self.ys.discard(y)
                self.e2 += 2
                self.stats['elim'] += 1
                pull()
                continue
            cs = [(mo, self.ph[mo]) for mo in terms]
            if all(c == half for _, c in cs):
                Q = frozenset(mo ^ yb for mo, _ in cs)
                if Q == frozenset([0]):
                    self.zero = True
                    return
                cand = self._hh_candidates(Q)
                if not cand:
                    continue
                z = min(cand, key=lambda v: (len(self.occ.get(v, ())), v))
                R = Q ^ frozenset([1 << z])
                for mo, c in cs:
                    self.addc(mo, -c)
                self.ys.discard(y)
                self.subst(z, R)
                self.ys.discard(z)
                self.e2 += 2
                self.stats['hh'] += 1
                if y in self.occ or y in self.wcnt:
                    raise AssertionError('HH left the summed variable in place')
                pull()
                continue
            lin = self.ph.get(yb)
            if lin in (q1, q3) and all(c == half for mo, c in cs if mo != yb):
                Q = frozenset(mo ^ yb for mo, _ in cs if mo != yb)
                for mo, c in cs:
                    self.addc(mo, -c)
                self.ys.discard(y)
                if lin == q1:
                    self.addc(0, N >> 3)
                    self.add_lift(Q, q3)
                else:
                    self.addc(0, 7 * (N >> 3))
                    self.add_lift(Q, q1)
                self.e2 += 1
                self.stats['omega'] += 1
                if y in self.occ:
                    raise AssertionError('omega left the summed variable in place')
                pull()
                continue
        self.check_limits()

    # ------------------------------------------------------------------ inspection
    def size(self):
        return dict(ys=len(self.ys), terms=len(self.ph), wmax=max((len(P) for P in self.w), default=0))

    def evaluate(self, val):
        """(output basis state as int over the wires, F mod 2^K) at a full assignment val of all live variables."""
        out = 0
        for q, P in enumerate(self.w):
            if peval(P, val):
                out |= 1 << q
        f = 0
        for mo, c in self.ph.items():
            if mo & val == mo:
                f += c
        return out, f & self.M
