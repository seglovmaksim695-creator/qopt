"""Scalable "logical AND + measurement-based uncomputation" compiler for reversible (Toffoli-level) circuits.

Every non-linear product L_a*L_b is computed with a logical AND into a fresh ancilla (4 T, 3 CX) and erased
by an X-basis measurement with a Clifford correction (CZ/Z on the CURRENT wire values, 0 T).
All wire values are tracked exactly as polynomials over GF(2) in the inputs, so reuse of products, moving the
"carrier" of a logical wire into an ancilla, linear expressions over live wires and quadratic corrections are
found automatically.

Strict model: all .qc wires carry data in any state; extra ancillas start in |0>.
Clean model (``inputs`` given): .qc wires not listed in ``.i`` start in |0> as well.

Scalability changes w.r.t. the original c2 compiler (lab agent b1):
interned monomials (a polynomial is an int bitset over a table of seen monomials, not over 2^n masks),
optional cut K (a wire value with more than K monomials is renamed to a fresh variable; exact), incremental
quad_rep/lin_rep around the support of the target, versioned state, fixes of the lazy mode, and ``lctl``
(lazy controls: cx/ccx take the delta from LOGICAL values and realize controls only when needed).

Note: the monomial table is module-global state (``reset_monos()`` before each compilation); run one
compilation per process at a time (the bandc worker does).
"""
import random
import re

import numpy as _np

from qoptct.portfolio import evo_lib
from qoptct.qopt_core import load_qc
from . import peep


class Fail(Exception): pass


def logical(hl, n):
    """H-frames -> logical classical circuit over x, z, cx, cz, ccx, ccz."""
    fr = [0] * n; out = []
    for g, q in hl:
        q = tuple(q)
        if g == 'h': fr[q[0]] ^= 1; continue
        if g in ('ccz', 'cz', 'z'):
            F = [x for x in q if fr[x]]
            if not F: out.append((g, q))
            elif len(F) == 1:
                t = F[0]; c = tuple(x for x in q if x != t)
                out.append(({'ccz': 'ccx', 'cz': 'cx', 'z': 'x'}[g], c + (t,)))
            else: raise Fail('2 frames')
        elif g == 'x':
            out.append(('z' if fr[q[0]] else 'x', q))
        elif g in ('cx', 'ccx'):
            if any(fr[x] for x in q[:-1]): raise Fail('ctrl in frame')
            out.append((({'cx': 'cz', 'ccx': 'ccz'}[g]) if fr[q[-1]] else g, q))
        else: raise Fail(g)
    if any(fr): raise Fail('open frame')
    return out


# ---------- polynomials: int bitset over interned monomials (monomial = variable mask) ----------
# bit j of a polynomial = monomial MONO[j] (variable mask); ONE = monomial 0 = bit 0.
# The int length is the number of DISTINCT monomials seen in a run (not 2^n), so 30-1500 variables work.
ONE = 1
MAXMONO = 2500000
COMPACT_AT = 150000
MONO = [0]; MID = {0: 0}
def mid(m):
    j = MID.get(m)
    if j is None:
        j = len(MONO)
        if j > MAXMONO: raise Fail('monos')   # memory guard
        MONO.append(m); MID[m] = j
    return j
def var(i): return 1 << mid(1 << i)
def reset_monos():
    MONO[:] = [0]; MID.clear(); MID[0] = 0; _mc.clear(); _sc.clear()
def bits(x):
    if x.bit_length() > 2048:   # long ints: via numpy instead of x & -x (O(size) per bit)
        b = x.to_bytes((x.bit_length() + 7) // 8, 'little')
        return _np.flatnonzero(_np.unpackbits(_np.frombuffer(b, dtype=_np.uint8), bitorder='little')).tolist()
    r = []
    while x:
        lo = x & -x; r.append(lo.bit_length() - 1); x ^= lo
    return r
def mask_key(p):
    """polynomial key independent of the monomial numbering"""
    return frozenset(MONO[j] for j in bits(p))
def from_masks(ms):
    r = 0
    for m in ms: r |= 1 << mid(m)
    return r
def nterms(p): return bin(p).count('1')
_sc = {}
def supp(p):
    """variable mask of a polynomial (cached)"""
    r = _sc.get(p)
    if r is None:
        r = 0
        for j in bits(p): r |= MONO[j]
        if len(_sc) > 200000: _sc.clear()
        _sc[p] = r
    return r
def degree(p): return max((bin(MONO[j]).count('1') for j in bits(p)), default=0)
_DENSE_MAXV = 20
def _supp_list(p):
    m = 0
    for j in bits(p): m |= MONO[j]
    return m
def _dense_ok(p, q):
    """dense multiplication (Moebius transform) pays off when |A||B| is large and there are few variables"""
    a = nterms(p); b = nterms(q)
    if a * b < 600: return False
    s = bin(supp(p) | supp(q)).count('1')
    return s <= _DENSE_MAXV and a * b > 3.5 * (1 << s)
def _moebius(arr, v):
    for i in range(v):
        a = arr.reshape(-1, 2, 1 << i); a[:, 1, :] ^= a[:, 0, :]
    return arr
def _mul_dense(p, q):
    S = supp(p) | supp(q); vs = []; x = S
    while x:
        lo = x & -x; vs.append(lo); x ^= lo
    v = len(vs)
    def loc(m):
        r = 0
        for i, b in enumerate(vs):
            if m & b: r |= 1 << i
        return r
    A = _np.zeros(1 << v, dtype=_np.uint8); B = _np.zeros(1 << v, dtype=_np.uint8)
    for j in bits(p): A[loc(MONO[j])] = 1
    for j in bits(q): B[loc(MONO[j])] = 1
    C = _moebius(_moebius(A, v) & _moebius(B, v), v)   # ANF -> table, AND, table -> ANF (Moebius is an involution)
    r = 0
    for idx in _np.nonzero(C)[0].tolist():
        m = 0; i = 0
        while idx:
            if idx & 1: m |= vs[i]
            idx >>= 1; i += 1
        r |= 1 << mid(m)
    return r
_mc = {}
def mul(p, q):
    if p > q: p, q = q, p
    k = (p, q)
    r = _mc.get(k)
    if r is not None: return r
    if p == ONE: r = q
    elif p == 0: r = 0
    elif _dense_ok(p, q): r = _mul_dense(p, q)
    else:
        A = [MONO[j] for j in bits(p)]; B = [MONO[j] for j in bits(q)]
        s = set()
        for a in A:
            for b in B:
                m = a | b
                if m in s: s.discard(m)
                else: s.add(m)
        r = 0
        for m in s: r |= 1 << mid(m)
    if len(_mc) > 400000: _mc.clear()
    _mc[k] = r
    return r


def solve(target, cols):
    """target as XOR of a subset of cols (list of (key, vector)); cols order = priority. -> list of keys or None."""
    basis = {}  # pivot -> (vec, set_of_keys)
    for key, v in cols:
        s = {key}
        while v:
            p = v.bit_length() - 1
            if p in basis:
                bv, bs = basis[p]; v ^= bv; s = s ^ bs
            else:
                basis[p] = (v, s); break
    v = target; s = set()
    while v:
        p = v.bit_length() - 1
        if p not in basis: return None
        bv, bs = basis[p]; v ^= bv; s ^= bs
    return list(s)


def solve_min(target, cols, wt):
    """like solve, then greedily reduces the weight wt(key) of the solution using kernel vectors."""
    basis = {}; kernel = []
    for key, v in cols:
        s = {key}
        while v:
            p = v.bit_length() - 1
            if p in basis:
                bv, bs = basis[p]; v ^= bv; s = s ^ bs
            else:
                basis[p] = (v, s); break
        if not v: kernel.append(frozenset(s))
    v = target; s = set()
    while v:
        p = v.bit_length() - 1
        if p not in basis: return None
        bv, bs = basis[p]; v ^= bv; s ^= bs
    W = lambda z: sum(wt(k) for k in z)
    cur = W(s); imp = True
    while imp:
        imp = False
        for k in kernel:
            z = s ^ k; w = W(z)
            if w < cur: s, cur, imp = z, w, True
    return list(s)


class _Ver:
    __slots__ = ('v',)
    def __init__(self): self.v = 0


class VList(list):
    """list of wire values with a version number (O(1) state key for caches) and an inverted index
    variable -> live wires (owner.inv), updated on every write"""
    def __init__(self, it, ver, owner): super().__init__(it); self.ver = ver; self.owner = owner
    def __setitem__(self, i, x):
        old = list.__getitem__(self, i)
        if old != x:
            self.ver.v += 1
            list.__setitem__(self, i, x)
            self.owner._reidx(i, old, x)


class VSet(set):
    def __init__(self, it, ver, owner): super().__init__(it); self.ver = ver; self.owner = owner
    def add(self, x):
        if x not in self: self.ver.v += 1; set.add(self, x); self.owner._reidx(x, 0, self.owner.V[x], live_change=True)
    def discard(self, x):
        if x in self: self.ver.v += 1; set.discard(self, x); self.owner._reidx(x, self.owner.V[x], 0, live_change=True)


class Compiler:
    def __init__(self, n, lg, max_anc, rng, eps=0.0, keep_prob=1.0, lazy=True, linmax=3, inputs=None, K=None, qfull=False, qmax=48, lctl=False, dtries=1, lin_force=True):
        self.n, self.lg, self.rng, self.eps, self.keep_prob = n, lg, rng, eps, keep_prob
        self.lazy, self.linmax, self.pend, self.guard, self.k = lazy, linmax, {}, False, 0
        self.depth, self.MAXD, self.busy, self.locked = 0, 40, set(), set()
        self.qfull, self.QFULL_MAX = qfull, qmax
        self.lctl = lctl
        self.dtries = dtries
        self.lin_force = lin_force
        self.prot = []   # (polynomial, excl): representations nested erasures must not break
        self._vi, self._vi_ver = None, -1
        self._pc = {}
        self._ck, self._cv, self._memo, self._msk = None, None, {}, None
        self.N = n + max_anc
        # clean model: wires not in inputs start in |0> (value 0); variables are inputs only
        if inputs is None: inputs = list(range(n))
        vi = {q: k for k, q in enumerate(inputs)}
        init = [var(vi[i]) if i in vi else 0 for i in range(n)]
        self._ver = _Ver()
        self.inv = {}
        self.live = set()
        self.V = VList(init + [0] * max_anc, self._ver, self)
        self.L = list(init)
        self.car = list(range(n))
        self.live = VSet((), self._ver, self)
        for _w in range(n): self.live.add(_w)   # physical data wires (data is always live)
        self.free = list(range(n, self.N))
        self.ops = []; self.nand = 0; self.used = set()
        # cut: if the new logical value of a wire has more than K monomials it is replaced by a fresh
        # variable y (y := that value). Identities in the free algebra (x, y) stay valid after substitution,
        # so everything remains exact (only some reuse is lost). Cut points (step, wire) are fixed
        # deterministically here; run() realizes and renames the wire at the same point.
        self.K = K; self.cuts = {}; nv = len(inputs)
        # future needs
        L = list(self.L); seq = []
        for k, (g, q) in enumerate(lg):
            need = set()
            if g in ('ccx', 'ccz'):
                a, b, c = q
                need |= {mul(L[a], L[b]), mul(L[a], L[c]), mul(L[b], L[c])}
                if g == 'ccx': L[c] ^= mul(L[a], L[b])
            elif g == 'cx': L[q[1]] ^= L[q[0]]
            elif g == 'x': L[q[0]] ^= ONE
            if K and g in ('cx', 'ccx'):
                t = q[-1]
                if nterms(L[t]) > K:
                    self.cuts[k] = (t, nv); L[t] = var(nv); nv += 1
            need |= set(L)
            seq.append(need)
        self.final = list(L)
        # last step at which a value is needed (instead of a set per step: O(G * |needs|) memory);
        # keys are monomial sets (not ints) so they survive compaction of the monomial table
        self.last = {}
        for k, s in enumerate(seq):
            for v in s: self.last[mask_key(v)] = k
        for v in L: self.last[mask_key(v)] = len(lg)
        del seq
        self.G = len(lg); self._mkc = {}; self.compact_at = COMPACT_AT; self.ncompact = 0

    def _mk(self, v):
        r = self._mkc.get(v)
        if r is None:
            r = mask_key(v)
            if len(self._mkc) > 100000: self._mkc.clear()
            self._mkc[v] = r
        return r

    def needed(self, v, k):
        return self.last.get(self._mk(v), -1) >= min(k, self.G)

    def _compact(self):
        """compact the monomial table: keep only monomials of live polynomials (V, L, pending, final),
        renumber, reset caches. Keeps int polynomials short on long circuits."""
        old = list(MONO)
        def keys(p): return [old[j] for j in bits(p)]
        V = [keys(p) for p in self.V]; L = [keys(p) for p in self.L]; F = [keys(p) for p in self.final]
        P = {t: keys(D) for t, D in self.pend.items()}
        reset_monos()
        lv = list(self.live); self.inv = {}; self.live = VSet((), self._ver, self)
        self.V = VList([from_masks(x) for x in V], self._ver, self); self._ver.v += 1
        for _w in lv: self.live.add(_w)
        self.L = [from_masks(x) for x in L]; self.final = [from_masks(x) for x in F]
        self.pend = {t: from_masks(x) for t, x in P.items()}
        self._mkc.clear(); self._memo = {}; self._msk = None; self._ck = None; self._cv = None; self._pc = {}
        self.ncompact += 1

    def carriers(self): return set(self.car)

    def emit(self, g, q): self.ops.append((g, q))

    def pmul(self, a, b):
        """product of wire values a, b, cached per wire pair (checked by int identity, no hashing of long ints)"""
        e = self._pc.get((a, b))
        va, vb = self.V[a], self.V[b]
        if e is not None and e[0] is va and e[1] is vb: return e[2]
        r = mul(va, vb)
        if len(self._pc) > 300000: self._pc.clear()
        self._pc[(a, b)] = (va, vb, r)
        return r

    def _reidx(self, w, old, new, live_change=False):
        if not live_change and w not in self.live: return
        so, sn = supp(old) if old else 0, supp(new) if new else 0
        x = so & ~sn
        while x:
            lo = x & -x; x ^= lo; st = self.inv.get(lo)
            if st is not None: st.discard(w)
        x = sn & ~so
        while x:
            lo = x & -x; x ^= lo; self.inv.setdefault(lo, set()).add(w)

    def near(self, sm, excl=()):
        """live wires whose variables intersect the mask sm (via the inverted index)"""
        r = set(); x = sm
        while x:
            lo = x & -x; x ^= lo; st = self.inv.get(lo)
            if st: r |= st
        return sorted(r - set(excl))

    def _vidx(self):
        if self._vi_ver != self._ver.v:
            d = {}
            for w in self.live: d.setdefault(self.V[w], []).append(w)
            self._vi, self._vi_ver = d, self._ver.v
        return self._vi

    def lin_rep(self, target, excl=(), maxw=None):
        if target == 0: return [], False
        if target == ONE: return [], True
        d = self._vidx()
        for t2, cst in ((target, False), (target ^ ONE, True)):
            for w in sorted(d.get(t2, ())):
                if w not in excl: return [w], cst
        sm = supp(target); w1 = self.near(sm, excl)
        for t2, cst in ((target, False), (target ^ ONE, True)):
            for w in w1:   # at least one wire of the pair intersects the support of target
                for u in d.get(t2 ^ self.V[w], ()):
                    if u != w and u not in excl: return [w, u], cst
        if maxw is not None and maxw <= 2: return None
        best = None
        # first wires sharing variables (levels 1, 2); the full set only for small circuits
        sm2 = sm
        for w in w1: sm2 |= supp(self.V[w])
        w2 = self.near(sm2, excl)
        small = self.qfull or len(self.live) <= self.QFULL_MAX
        ws = sorted(w for w in self.live if w not in excl) if small else []
        levels = [w1] + ([w2] if len(w2) > len(w1) else []) + ([ws] if len(ws) > len(w2) else [])
        for wl in levels:
            for t2, cst in ((target, False), (target ^ ONE, True)):
                r = solve(t2, [(w, self.V[w]) for w in wl])
                if r is not None and (best is None or len(r) < len(best[0])): best = (r, cst)
            if best: break
        if best and maxw is not None and len(best[0]) > maxw: return None
        return best

    def _skey(self): return self._ver.v

    def _cols(self):
        """columns (linear, constant, pairwise products) for ALL live wires, cached per state."""
        key = self._skey()
        if self._ck != key:
            ws = sorted(self.live)
            cols = [(('l', w), self.V[w]) for w in ws] + [(('c',), ONE)]
            for i, a in enumerate(ws):
                for b in ws[i + 1:]:
                    cols.append((('q', a, b), self.pmul(a, b)))
            self._ck, self._cv = key, cols
        return self._cv

    def quad_rep(self, target, excl, need=False):
        """target = const + sum lin + sum quad over live wires except excl -> (lin, pairs) or None.
        need=False: feasibility only (no CZ minimization over the kernel), memoized within a state."""
        r = self.lin_rep(target, excl)
        if r is not None: return r[0], []
        mk = (target, frozenset(excl), need)
        sk = self._skey()
        if sk != self._msk: self._msk, self._memo = sk, {}
        if mk in self._memo: return self._memo[mk]
        ex = set(excl)
        # first only wires whose variables intersect the variables of target (sound; usually enough)
        sm = supp(target)
        ws = self.near(sm, ex)
        if len(ws) + len(ex) < len(self.live):
            cs = [(('l', w), self.V[w]) for w in ws] + [(('c',), ONE)]
            for i, a in enumerate(ws):
                for b in ws[i + 1:]:
                    cs.append((('q', a, b), self.pmul(a, b)))
            s = solve_min(target, cs, lambda k: 1 if k[0] == 'q' else 0) if need else solve(target, cs)
            if s is not None:
                res = ([k[1] for k in s if k[0] == 'l'], [(k[1], k[2]) for k in s if k[0] == 'q'])
                self._memo[mk] = res
                return res
        if not self.qfull and len(self.live) > self.QFULL_MAX:
            # second level: wires intersecting the joint support of the first level (at most QFULL_MAX)
            sm2 = sm
            for w in ws: sm2 |= supp(self.V[w])
            ws2 = self.near(sm2, ex)
            if len(ws2) > len(ws) and len(ws2) <= self.QFULL_MAX:
                cs = [(('l', w), self.V[w]) for w in ws2] + [(('c',), ONE)]
                for i, a in enumerate(ws2):
                    for b in ws2[i + 1:]:
                        cs.append((('q', a, b), self.pmul(a, b)))
                s = solve_min(target, cs, lambda k: 1 if k[0] == 'q' else 0) if need else solve(target, cs)
                if s is not None:
                    res = ([k[1] for k in s if k[0] == 'l'], [(k[1], k[2]) for k in s if k[0] == 'q'])
                    self._memo[mk] = res
                    return res
            # on large circuits a full search over all pairs of live wires is too costly (O(w^2) columns);
            # answering "no" here is sound (the compiler just acts more conservatively)
            self._memo[mk] = None
            return None
        cols = self._cols()
        cs = [c for c in cols if not (c[0][0] == 'l' and c[0][1] in ex or c[0][0] == 'q' and (c[0][1] in ex or c[0][2] in ex))]
        s = solve_min(target, cs, lambda k: 1 if k[0] == 'q' else 0) if need else solve(target, cs)
        res = None if s is None else ([k[1] for k in s if k[0] == 'l'], [(k[1], k[2]) for k in s if k[0] == 'q'])
        self._memo[mk] = res
        return res

    def mz(self, a):
        r = self.quad_rep(self.V[a], {a}, need=True)
        if r is None: return False
        lin, pairs = r
        corr = [('cz', p) for p in pairs] + [('z', (w,)) for w in lin] + [('x', (a,))]
        self.emit('MZ', (a, tuple(corr)))
        self.V[a] = 0; self.live.discard(a)
        if a >= self.n: self.free.append(a)
        return True

    def garbage(self):
        # an "orphaned" data wire (its logical value lives in an ancilla) is garbage too: its old
        # value is erased by measurement once unneeded, otherwise it would get stuck after the swap
        cs = self.carriers()
        return [w for w in self.live if w not in cs and (w >= self.n or self.car[w] != w)]

    def _rel(self, target, chg):
        """can changing wires with variables chg affect the representability of target (on large circuits
        quad_rep/lin_rep only look at level 1-2 wires around the support of target)"""
        if chg is None or self.qfull or len(self.live) <= self.QFULL_MAX: return True
        sm = supp(target)
        if chg & sm: return True
        sm2 = sm
        for w in self.near(sm): sm2 |= supp(self.V[w])
        return bool(chg & sm2)

    def ok_state(self, skip=None, chg=None):
        """all garbage (except skip) is quadratically expressible and all pending deltas are representable."""
        for a in self.garbage():
            if a == skip: continue
            if self._rel(self.V[a], chg) and self.quad_rep(self.V[a], {a}) is None: return False
        for t, D in self.pend.items():
            if D and self.car[t] != skip and self._rel(D, chg) and self.quad_rep(D, {self.car[t]}) is None: return False
        return True

    def _prot_ok(self, a):
        """erasing wire a does not break protected representations (the delta being realized etc.)"""
        if not self.prot: return True
        self.live.discard(a)
        try:
            return all(self.quad_rep(p, ex) is not None for p, ex in self.prot)
        finally:
            self.live.add(a)

    def safe_mz(self, a):
        if a in self.locked: return False
        if not self._prot_ok(a): return False
        if self.quad_rep(self.V[a], {a}) is None: return False
        chg = supp(self.V[a])
        self.live.discard(a); ok = self.ok_state(chg=chg); self.live.add(a)
        if not ok: return False
        return self.mz(a)

    def mz_force(self, a):
        """erase a; first realize pending deltas and erase garbage that cannot be expressed without a."""
        # recursive: erasing dependent garbage b may itself break other pending deltas/garbage
        self._mzf(a, set(), 0)

    def _mzf(self, a, removing, d):
        if d > 30: raise Fail('mzf depth')
        removing.add(a)
        for it in range(6):   # until stable: realizing deltas changes wires -> recompute dependencies
            back = [x for x in removing if x in self.live]
            for x in back: self.live.discard(x)
            dep_p = [t for t, D in self.pend.items() if D and t not in self.busy and self.quad_rep(D, {self.car[t]}) is None]
            dep_g = [b for b in self.garbage() if b not in self.locked and b not in removing
                     and self.quad_rep(self.V[b], {b}) is None]
            for x in back: self.live.add(x)
            if not dep_p and not dep_g: break
            for t in dep_p: self.realize(t, self.k)
            if dep_p: continue
            for b in dep_g:
                if b in self.live and b not in self.carriers() and b not in removing and self._prot_ok(b): self._mzf(b, removing, d + 1)
        if a in self.live and self._prot_ok(a): self.mz(a)

    def before_modify(self, w, newval):
        # nested realizations are protected too (depth <= MAXD, targets in progress are skipped)
        # instead of being disabled by a guard flag (in c2 nested erasures broke other pending deltas).
        if self.guard: return
        if self.depth >= self.MAXD: raise Fail("depth")
        # repeat until stable: fixing one item (realize/erase) may break the representation
        # of another one that relied on the fixed wire (seen on qcla_mod_7)
        for it in range(6):
            old = self.V[w]
            bad_g, bad_p = [], []
            chg = supp(old) | supp(newval)
            pre_g = {a: self.quad_rep(self.V[a], {a}) is not None for a in self.garbage()
                     if a != w and self._rel(self.V[a], chg)}
            pre_p = {t: self.quad_rep(D, {self.car[t]}) is not None for t, D in self.pend.items()
                     if D and self.car[t] != w and t not in self.busy and self._rel(D, chg)}
            self.V[w] = newval
            for a, ok in pre_g.items():
                if ok and self.quad_rep(self.V[a], {a}) is None: bad_g.append(a)
            for t, ok in pre_p.items():
                if ok and self.quad_rep(self.pend[t], {self.car[t]}) is None: bad_p.append(t)
            self.V[w] = old
            if not bad_g and not bad_p: return
            self.depth += 1
            try:
                for t in bad_p:
                    if self.pend.get(t): self.realize(t, self.k)
                for a in bad_g:
                    if a in self.live and a not in self.carriers() and a not in self.locked: self.mz_force(a)
            finally:
                self.depth -= 1
            if self.V[w] != old: return   # the wire was changed by nested actions; the caller re-checks

    def cx(self, c, t):
        self.before_modify(t, self.V[t] ^ self.V[c])
        self.emit('cx', (c, t)); self.V[t] ^= self.V[c]

    def rawcx(self, c, t):
        self.emit('cx', (c, t)); self.V[t] ^= self.V[c]

    def helper(self, target, excl):
        """add one AND(W_i, W_j) into a fresh ancilla so that target becomes quadratically expressible."""
        ws = sorted(w for w in self.live if w not in excl)
        if not self.free: return False
        a = self.free[-1]
        for i, x in enumerate(ws):
            for y in ws[i + 1:]:
                self.V[a] = mul(self.V[x], self.V[y]); self.live.add(a)
                ok = self.quad_rep(target, set(excl)) is not None
                self.live.discard(a); self.V[a] = 0
                if ok:
                    b = self.fresh(); assert b == a
                    self.emit('AND', (x, y, a)); self.nand += 1; self.V[a] = mul(self.V[x], self.V[y])
                    return True
        return False

    def fresh(self):
        if not self.free:
            for a in sorted(self.garbage()):
                if self.safe_mz(a): break
        if not self.free: raise Fail('no ancilla')
        a = self.free.pop(); self.live.add(a); return a

    def try_AND(self, a, b):
        if not self.free:
            for x in sorted(self.garbage()):
                if self.safe_mz(x): break
        if not self.free: return None
        return self.AND(a, b)

    def AND(self, a, b):
        t = self.fresh(); self.emit('AND', (a, b, t)); self.nand += 1
        self.V[t] = self.pmul(a, b); return t

    def cleanup(self, k):
        for a in sorted(self.garbage()):
            if self.needed(self.V[a], k) and self.rng.random() < self.keep_prob: continue
            self.safe_mz(a)

    def decompose(self, pairs):
        """several choices of the pivot pair (random tie breaks); keep the decomposition
        with the minimum (number of products, number of CX for linear combinations)"""
        best = None
        for it in range(self.dtries):
            r = self._decompose(pairs, None if it == 0 else random.Random(it))
            cost = (len(r[0]), sum(len(a) + len(b) - 2 for a, b in r[0]) * 2 + len(r[1]))
            if best is None or cost < best[0]: best = (cost, r)
        return best[1]

    @staticmethod
    def _decompose(pairs, rng=None):
        Q = set(pairs); prods = []; lin = set()
        while Q:
            cnt = {}
            for i, j in Q: cnt[i] = cnt.get(i, 0) + 1; cnt[j] = cnt.get(j, 0) + 1
            if rng is None: i, j = min(Q, key=lambda p: (cnt[p[0]] + cnt[p[1]], p))
            else: i, j = min(Q, key=lambda p: (cnt[p[0]] + cnt[p[1]], rng.random()))
            A = {k for p in Q for k in p if i in p and j not in p and k != i}
            B = {k for p in Q for k in p if j in p and i not in p and k != j}
            for p in list(Q):
                if i in p or j in p: Q.discard(p)
            for a in A:
                for b in B:
                    if a == b: lin ^= {a}
                    else: Q ^= {(min(a, b), max(a, b))}
            prods.append(({i} | B, {j} | A))
        return prods, lin

    def realize(self, t, k):
        if not self.pend.get(t, 0) or t in self.busy: return
        self.busy.add(t); item = (self.pend[t], frozenset({self.car[t]}))
        self.prot.append(item)
        try: self._realize(t, k)
        finally:
            self.busy.discard(t)
            if item in self.prot: self.prot.remove(item)

    def _realize(self, t, k):
        D = self.pend.get(t, 0)
        if not D: return
        newL = self.L[t]; ct = self.car[t]; old = self.V[ct]
        cs = self.carriers()
        for w in sorted(self.live):
            if w in cs and w != ct: continue
            if w < self.n and w != t: continue
            if self.V[w] in (newL, newL ^ ONE):
                if w != ct and ct >= self.n:
                    self.car[t] = w
                    ok = self.quad_rep(self.V[ct], {ct}) is not None
                    self.car[t] = ct
                    if not ok: continue
                if self.V[w] != newL: self.emit('x', (w,)); self.V[w] ^= ONE
                self.car[t] = w; self.pend[t] = 0; return
        lr = self.lin_rep(D, {ct}, maxw=self.linmax)
        if newL in (0, ONE) and old not in (0, ONE):
            # erase by measurement: the new value is a constant, the old one is quadratic in the others
            r = self.quad_rep(old, {ct}, need=True)
            if r is not None and (lr is None or len(lr[0]) > len(r[1])):
                self.before_modify(ct, 0)
                if self.V[ct] != old: raise Fail('perturbed-e1')
                r = self.quad_rep(old, {ct}, need=True)   # recompute after nested actions
                if r is not None:
                    corr = [('cz', p) for p in r[1]] + [('z', (w,)) for w in r[0]] + [('x', (ct,))]
                    self.emit('MZ', (ct, tuple(corr))); self.V[ct] = 0
                    if newL == ONE: self.emit('x', (ct,)); self.V[ct] = ONE
                    self.pend[t] = 0; return
                lr = self.lin_rep(D, {ct}, maxw=self.linmax)   # erasure became impossible: regular path
        for _ in range(3 if self.lin_force else 0):   # lin_force=False: compiler version b1/v7
            # first realize pending deltas whose carriers occur in lr (otherwise nested protection changes them)
            if lr is None: break
            cs_ = {self.car[s_]: s_ for s_, D_ in self.pend.items() if D_ and s_ not in self.busy}
            dep = [cs_[w] for w in lr[0] if w in cs_]
            if not dep: break
            for s_ in dep: self.realize(s_, k)
            if self.V[ct] != old or self.pend.get(t) != D: raise Fail('perturbed-lin0')
            lr = self.lin_rep(D, {ct}, maxw=self.linmax)
        if lr is not None:
            # wires of lr are locked (nested erasures do not touch them); afterwards check they did not change
            lk = set(lr[0]) - self.locked; self.locked |= lk
            snap = {w: self.V[w] for w in lr[0]}
            try:
                for w in lr[0]:
                    if self.V[w] != snap[w]: raise Fail('perturbed-lin')
                    self.cx(w, ct)
            finally:
                self.locked -= lk
            if lr[1]: self.emit('x', (ct,)); self.V[ct] ^= ONE
            if self.V[ct] != newL: raise Fail('perturbed-lin2')
            self.pend[t] = 0; return
        move = self.needed(old, k + 1)
        if not move and len(self.prot) > 1 and self.free:
            # if changing the carrier would break a protected representation (a delta being realized above relies
            # on the old value), the new value goes into a fresh ancilla (move) and the old one stays
            me = (D, frozenset({ct})); self.V[ct] = newL
            try: move = any(self.quad_rep(p, ex) is None for p, ex in self.prot if (p, ex) != me)
            finally: self.V[ct] = old
        if self.rng.random() < self.eps: move = not move
        if old in (0, ONE): move = False
        # protect BEFORE computing the representation (c2 did it after: nested realizations changed wires
        # already chosen for the decomposition -> wrong result / assert)
        if not move: self.before_modify(ct, newL)
        if self.V[ct] != old or self.pend.get(t) != D: raise Fail('perturbed-q %s %s' % (self.V[ct] != old, self.pend.get(t) != D))
        qr = self.quad_rep(D, {ct}, need=True)
        if qr is None and self.helper(D, {ct}): qr = self.quad_rep(D, {ct}, need=True)
        if qr is None: raise Fail('unrep')
        lin, pairs = qr
        prods, lx = self.decompose([tuple(sorted(p)) for p in pairs])
        lin = set(lin) ^ lx
        lk = set(lin) | {w for p in pairs for w in p}
        self.locked |= lk
        try:
            self._apply(t, ct, newL, prods, lin, move, lk)
        finally:
            self.locked -= lk

    def _apply(self, t, ct, newL, prods, lin, move, lk):
        tgt = ct
        for idx, (Ls, Ms) in enumerate(prods):
            u = min(Ls, key=lambda w: (w in Ms, w)); v = min(Ms, key=lambda w: (w in Ls, w))
            assert u not in Ms and v not in Ls
            mods = [(w, u) for w in Ls if w != u] + [(w, v) for w in Ms if w != v]
            for c_, t_ in mods: self.rawcx(c_, t_)
            if self.V[tgt] in (0, ONE) and not (idx == 0 and move):
                # the target is clean now: AND directly into it (no ancilla, no uncomputation)
                if self.V[tgt] == ONE: self.emit('x', (tgt,)); self.V[tgt] = 0
                self.emit('AND', (u, v, tgt)); self.nand += 1; self.V[tgt] = mul(self.V[u], self.V[v])
                for c_, t_ in reversed(mods): self.rawcx(c_, t_)
                continue
            g = self.try_AND(u, v)
            if g is None:
                # protect here too (in c2 with move=True and no ancillas CCX changed a wire unchecked)
                pr = mul(self.V[u], self.V[v]); vt = self.V[tgt]
                if move or idx > 0:
                    sv = {w: self.V[w] for w in lk}
                    self.before_modify(tgt, self.V[tgt] ^ pr)
                    if self.V[tgt] != vt or any(self.V[w] != sv[w] for w in lk): raise Fail('perturbed-ccx')
                self.emit('CCX', (u, v, tgt)); self.V[tgt] ^= mul(self.V[u], self.V[v])
            elif idx == 0 and move:
                self.rawcx(tgt, g); self.car[t] = g; tgt = g
            else:
                self.rawcx(g, tgt)
                if mods: self.mz(g)
            for c_, t_ in reversed(mods): self.rawcx(c_, t_)
        for w in sorted(lin): self.rawcx(w, tgt)
        if self.V[tgt] != newL:
            assert self.V[tgt] ^ ONE == newL
            self.emit('x', (tgt,)); self.V[tgt] ^= ONE
        self.pend[t] = 0

    def _post(self, k, g, q):
        if k in self.cuts:
            t, y = self.cuts[k]
            self.realize(t, k)
            ct = self.car[t]
            if self.V[ct] != self.L[t]: raise Fail('cut mismatch')
            Y = var(y)
            self.before_modify(ct, Y)
            self.emit('CUT', (ct, y, k, t))
            self.V[ct] = Y; self.L[t] = Y
        self.cleanup(k + 1)

    def add(self, t, delta, k):
        lz = self.lazy if isinstance(self.lazy, bool) else (self.rng.random() < self.lazy)
        newp = self.pend.get(t, 0) ^ delta
        if not lz or (newp and self.quad_rep(newp, {self.car[t]}) is None):
            self.realize(t, k); newp = delta
        self.pend[t] = newp; self.L[t] ^= delta
        if not lz: self.realize(t, k)

    def run(self):
        for k, (g, q) in enumerate(self.lg):
            self.k = k
            if len(MONO) > self.compact_at:
                self._compact(); self.compact_at = max(COMPACT_AT, 3 * len(MONO))
            if g in ('z', 'cz', 'ccz'): rd = q
            elif g in ('cx', 'ccx'): rd = q[:-1]
            else: rd = ()
            # lctl: in lazy mode cx/ccx do NOT realize their control wires first; the delta is taken from
            # the logical values L; realization happens only if the target delta is otherwise unrepresentable
            if self.lctl and self.lazy is not False and g in ('cx', 'ccx'):
                t = q[-1]
                dl = self.L[q[0]] if g == 'cx' else mul(self.L[q[0]], self.L[q[1]])
                newp = self.pend.get(t, 0) ^ dl
                if newp == 0 or self.quad_rep(newp, {self.car[t]}) is not None:
                    self.pend[t] = newp; self.L[t] ^= dl
                    self._post(k, g, q)
                    continue
            for w in rd: self.realize(w, k)
            c = self.car
            if g == 'x':
                self.emit('x', (c[q[0]],)); self.V[c[q[0]]] ^= ONE; self.L[q[0]] ^= ONE
            elif g == 'z': self.emit('z', (c[q[0]],))
            elif g == 'cz': self.emit('cz', (c[q[0]], c[q[1]]))
            elif g == 'cx': self.add(q[1], self.V[c[q[0]]], k)
            elif g == 'ccx': self.add(q[2], mul(self.V[c[q[0]]], self.V[c[q[1]]]), k)
            elif g == 'ccz':
                opts = []
                for (x, y, z) in ((q[0], q[1], q[2]), (q[0], q[2], q[1]), (q[1], q[2], q[0])):
                    r = self.lin_rep(mul(self.L[x], self.L[y]), {c[z]}, maxw=3)
                    opts.append((len(r[0]) if r else 99, self.rng.random(), x, y, z, r))
                opts.sort(key=lambda o: o[:2])
                w, _, x, y, z, r = opts[0]
                if r is not None and w <= 3:
                    for u in r[0]: self.emit('cz', (u, c[z]))
                    if r[1]: self.emit('z', (c[z],))
                else:
                    a = self.try_AND(c[x], c[y])
                    if a is None: self.emit('CCZ', (c[x], c[y], c[z]))
                    else: self.emit('cz', (a, c[z]))
            else: raise Fail(g)
            self._post(k, g, q)
        self.k = len(self.lg)
        for t in range(self.n): self.realize(t, self.k)
        for w in range(self.n):
            if w not in self.live: self.live.add(w)   # erased orphaned data wire (= 0)
        for w in range(self.n):
            if self.car[w] != w:
                d = self.V[w] ^ self.L[w]
                r = self.lin_rep(d, {w})
                if r is not None:
                    for u in r[0]: self.cx(u, w)
                    if r[1]: self.emit('x', (w,)); self.V[w] ^= ONE
                else:
                    a = self.car[w]
                    for (x, y) in ((a, w), (w, a), (a, w)): self.rawcx(x, y)
                self.car[w] = w
        for _ in range(6):
            prog = True
            while prog and self.garbage():
                prog = False
                for a in sorted(self.garbage()):
                    if self.safe_mz(a): prog = True
            gs = [a for a in sorted(self.garbage()) if self.quad_rep(self.V[a], {a}) is None]
            if not gs: continue
            if not self.helper(self.V[gs[0]], {gs[0]}): break
        if self.garbage(): raise Fail('stuck garbage')
        for w in range(self.n):
            if self.V[w] != self.final[w]: raise Fail('final mismatch')
        return self.ops



def to_items(ops):
    it = []
    for g, q in ops:
        if g == 'AND':
            a, b, t = q
            for gg, qq in evo_lib.rccx(a, b, t) + [('sdg', (t,))]: it.append(('G', gg, tuple(qq)))
        elif g == 'MZ':
            # H before the measurement is part of the M element (the peephole must not cancel it against another H:
            # that would turn the X-measurement into a Z-measurement, equivalent but breaking sampled-check normalization)
            it.append(('M', q[0], q[1]))
        elif g in ('CCX', 'CCZ'):
            for gg, qq in evo_lib.lower_gate(g.lower(), q): it.append(('G', gg, tuple(qq)))
        elif g == 'CUT': continue   # cut marker (for lcheck), no gates
        else: it.append(('G', g, tuple(q)))
    return it


FOLDABLE = {'h', 'x', 'cx', 't', 'tdg', 's', 'sdg', 'z'}


def fold_items(items, N):
    """phase folding (evo_lib.fold_cancel) on unitary segments between measurements/other elements."""
    out, seg = [], []
    def flush():
        if seg:
            out.extend(('G', g, tuple(q)) for g, q in evo_lib.fold_cancel([(g, q) for _, g, q in seg], N))
            seg.clear()
    for it in items:
        if it[0] == 'G' and it[1] in FOLDABLE: seg.append(it)
        else: flush(); out.append(it)
    flush()
    return out


def to_qasm(ops, N, path=None, peephole=True, fold=False, ndata=None):
    """ndata: number of data wires; after the X-measurement of a data wire (erasing an "orphaned" wire) a reset
    is appended (the wire is already |0>; lets the path-sum checker treat the data wire as live; T/CX unchanged)."""
    items = to_items(ops)
    if fold: items = fold_items(items, N)
    if peephole: items = peep.run(items)
    lines = []; nb = 0
    for it in items:
        if it[0] == 'M':
            a, corr = it[1], it[2]
            lines.append(f"h q[{a}];"); lines.append(f"c[{nb}] = measure q[{a}];")
            body = " ".join(f"{gg} " + ", ".join(f"q[{x}]" for x in qq) + ";" for gg, qq in corr)
            lines.append(f"if (c[{nb}]) {{{body}}}"); nb += 1
            if ndata is not None and a < ndata: lines.append(f"reset q[{a}];")
        else:
            lines.append(f"{it[1]} " + ", ".join(f"q[{x}]" for x in it[2]) + ";")
    txt = "OPENQASM 3;\ninclude \"stdgates.inc\";\n" + f"qubit[{N}] q;\nbit[{max(nb,1)}] c;\n" + "\n".join(lines) + "\n"
    if path: open(path, 'w').write(txt)
    return txt


def count(txt):
    T = len(re.findall(r"\b(t|tdg) q", txt)); CX = len(re.findall(r"\b(cx|cz) q", txt))
    return T, CX


def compact(ops, n):
    """renumber used ancillas consecutively -> (ops, N)."""
    used = sorted({x for g, q in ops for x in (q[:3] if g == 'AND' else ((q[0],) + tuple(y for _, qq in q[1] for y in qq)) if g == 'MZ' else q[:1] if g == 'CUT' else q) if x >= n})
    mp = {u: n + i for i, u in enumerate(used)}
    f = lambda x: mp.get(x, x)
    out = []
    for g, q in ops:
        if g == 'MZ': out.append((g, (f(q[0]), tuple((gg, tuple(map(f, qq))) for gg, qq in q[1]))))
        elif g == 'CUT': out.append((g, (f(q[0]),) + tuple(q[1:])))
        else: out.append((g, tuple(map(f, q))))
    return out, n + len(used)


def inputs_idx(path, names):
    ins = None
    for line in open(path):
        if line.startswith('.i'): ins = line.split()[1:]
    return None if ins is None else [i for i, x in enumerate(names) if x in ins]


def compile_path(path, seed=0, max_anc=None, eps=0.0, keep_prob=1.0, lazy=True, linmax=3, clean=False):
    names, hl = load_qc(path); n = len(names)
    lg = logical(hl, n)
    if max_anc is None: max_anc = 22 - n
    C = Compiler(n, lg, max_anc, random.Random(seed), eps=eps, keep_prob=keep_prob, lazy=lazy, linmax=linmax,
                 inputs=inputs_idx(path, names) if clean else None)
    ops = C.run()
    ops, N = compact(ops, n)
    return ops, N, n
