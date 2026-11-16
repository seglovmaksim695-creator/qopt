# SPDX-License-Identifier: Apache-2.0
"""compose.kal -- the contract chain of the Kaliski inversion circuit (qoptct.crypto.kaliski, legacy=False).

One class, Chain, walks the algorithm P_std (Kaliski's almost-inverse in the "frame of w", exact integer semantics)
segment by segment. Every segment has a CONTRACT written once over a boolean algebra (compose.bv):

    spec(al, val) -> (posts, asserts, prefacts, postfacts)
        val(wire)   value of a wire at the entry of the segment
        posts       {wire: value at the exit}; every other wire keeps its value
        asserts     [(name, bit)]: MATHEMATICAL ASSUMPTIONS on the entry state (properties of the algorithm that are
                    proved on paper, compose/PROOF.md: M1..M5), used as the precondition
        prefacts    facts (see fact_bit) that must have been ESTABLISHED by earlier segments (machine-proved
                    postconditions carried over wires that were not touched in between), used as the precondition
        postfacts   facts proved at the exit of this segment

Two modes of the same walk:
  * prove (lifted circuit given): for every segment the slice of the REAL lifted circuit between two certificate
    marks is executed symbolically on ROBDDs over all entry values at full width, and it is proved that under the
    precondition  (i) every AND target is 0,  (ii) every erase condition holds,  (iii) every wire ends with the value
    of the contract (ancillas: 0),  (iv) the postfacts hold;  in clean mode also  (v) the mirror slice of the inverse
    half, started from the exit state, satisfies (i), (ii) and returns every wire to its entry value.
    The certificate (marks: positions and wire names, from KaliskiCircuit.hints()) is NOT trusted: a wrong mark
    makes a proof fail, it cannot make a wrong circuit pass.
  * run (inputs xs given): P_std itself is executed bit-sliced on the inputs; all asserts, facts and the final
    result are checked. This is the TEST of the mathematical assumptions (exhaustive for small p, sampled at n=256).
"""
import time

from ..bdd import BDD, ONE, ZERO
from . import bv
from .bv import IntAlg, BddAlg
from .lift import K_X, K_CX, K_AND

__all__ = ['Chain', 'ProofFailure', 'wspec', 'fact_bit', 'simulate_flat']


class ProofFailure(Exception):
    def __init__(self, msg, **detail):
        super().__init__(msg)
        self.detail = detail


def _inter(*lists):
    out = []
    for i in range(max(len(x) for x in lists)):
        for x in lists:
            if i < len(x):
                out.append(x[i])
    return out


def wspec(al, a, b, wp):
    """w_k of Kaliski's step in the frame of wp = w_{k-1}: (u, v) = (b, a) if wp else (a, b);
    w = [the step operates on v] = [u odd] and ([v even] or [u <= v])."""
    u0 = al.ite(wp, b[0], a[0])
    v0 = al.ite(wp, a[0], b[0])
    lee = al.ite(wp, bv.le(al, b, a), bv.le(al, a, b))
    return al.and_(u0, al.or_(al.not_(v0), lee))


def fact_wires(key):
    t = key[0]
    if t in ('lt', 'nand', 'imp'):
        return key[1]
    if t == 'geodd':
        return key[1] + key[2]
    if t == 'wdef':
        return (key[1],) + ((key[2],) if isinstance(key[2], int) else ()) + key[3] + key[4]
    raise ValueError(t)


def fact_bit(al, val, key):
    """the meaning of a fact, as a boolean over the wire values val:
    ('lt', wires, c)          value(wires) < c
    ('nand', (a, b))          not (a and b)
    ('imp', (a, b))           a implies b
    ('geodd', A, B)           A_0 implies value(B) <= value(A)
    ('wdef', w, wprev, A, B)  w = wspec(A, B, wprev)"""
    t = key[0]
    if t == 'lt':
        return bv.lt_const(al, [val(w) for w in key[1]], key[2])
    if t == 'nand':
        return al.not_(al.and_(val(key[1][0]), val(key[1][1])))
    if t == 'imp':
        return al.or_(al.not_(val(key[1][0])), val(key[1][1]))
    if t == 'geodd':
        a = [val(w) for w in key[1]]
        b = [val(w) for w in key[2]]
        return al.or_(al.not_(a[0]), bv.le(al, b, a))
    if t == 'wdef':
        a = [val(w) for w in key[3]]
        b = [val(w) for w in key[4]]
        return al.not_(al.xor(val(key[1]), wspec(al, a, b, val(key[2]))))
    raise ValueError(t)


class Chain:
    def __init__(self, p, reduction, hints, lifted=None, xs=None, clean=False, progress=None):
        self.p = p
        self.n = n = p.bit_length()
        self.red = reduction
        self.clean = bool(clean)
        self.marks = hints['marks']
        self.mi = 0
        self.prove = lifted is not None
        self.progress = progress
        self.owned = set(range(n))
        self.reserved = set(range(n, 2 * n)) if clean else set()
        self.owned |= self.reserved
        self.facts = set()
        self.assumed = {}
        self.admitted = []
        self.nseg = 0
        self.nodes = 0
        self.maxnodes = 0
        if self.prove:
            self.ops, self.nops, self.nq = lifted.ops, lifted.n, lifted.nq
            self.pos = 0
            self.cinv = 0
            self.Z = set(range(n, lifted.nq))
            self.F = int(hints.get('fwd_len', lifted.n)) if clean else lifted.n
            if not (0 <= self.F <= self.nops):
                raise ProofFailure('certificate: bad forward length')
        else:
            self.nq = hints['nqubits']
            self.al = IntAlg(len(xs))
            self.xs = list(xs)
            self.store = {}
            for i in range(n):
                m = 0
                for j, x in enumerate(self.xs):
                    if (x >> i) & 1:
                        m |= 1 << j
                self.store[i] = m
            self.viol = []

    # ------------------------------------------------------------------ certificate and bookkeeping
    def mark(self, label):
        if self.mi >= len(self.marks):
            raise ProofFailure('certificate: marks exhausted, expected %r' % label)
        pos, lab, info = self.marks[self.mi]
        if lab != label:
            raise ProofFailure('certificate: expected mark %r, found %r' % (label, lab))
        self.mi += 1
        if self.prove and not (isinstance(pos, int) and self.pos <= pos <= self.F):
            raise ProofFailure('certificate: bad position at mark %r' % label)
        return pos, info

    def peek(self):
        return self.marks[self.mi][1] if self.mi < len(self.marks) else None

    def _violation(self, label, k, name, bad):
        j = (bad & -bad).bit_length() - 1
        self.viol.append(dict(segment=label, k=k, what=name, x=self.xs[j]))

    def fresh(self, w):
        """a wire named by the certificate as freshly allocated: must be a valid, unowned wire that is known to be 0."""
        if not isinstance(w, int) or isinstance(w, bool) or not (0 <= w < self.nq) or w in self.owned:
            raise ProofFailure('certificate: wire %r is not a free wire' % (w,))
        if self.prove:
            if w not in self.Z:
                raise ProofFailure('certificate: wire %d is not known to be zero' % w)
        elif self.store.get(w, 0):
            self._violation('fresh', None, 'wire %d not zero' % w, self.store[w])
        self.owned.add(w)
        return w

    def fresh_list(self, ws):
        if not isinstance(ws, list):
            raise ProofFailure('certificate: wire list expected')
        return [self.fresh(w) for w in ws]

    def assume_zero(self, wires, name):
        """mathematical assumption: these register bits are 0 here."""
        if self.prove:
            new = [w for w in wires if w not in self.Z]
            if new:
                self.assumed[name] = self.assumed.get(name, 0) + 1
                self.Z.update(new)
        else:
            for w in wires:
                if self.store.get(w, 0):
                    self._violation('assume_zero', None, name, self.store[w])

    def require_zero(self, w, why):
        """a derived fact: this wire is 0 (must already be known from the proofs)."""
        if self.prove:
            if w not in self.Z:
                raise ProofFailure('wire %d is not proved to be zero (%s)' % (w, why))
        elif self.store.get(w, 0):
            self._violation('require_zero', None, why, self.store[w])

    def need(self, cond, msg):
        if not cond:
            raise ProofFailure(msg)

    def establish(self, key):
        """a fact that is a tautology given the wires known to be zero."""
        if self.prove:
            B = BDD()
            al = BddAlg(B)
            var = {}
            for w in fact_wires(key):
                if w not in self.Z and w not in var:
                    var[w] = B.new_var()[1]
            bit = fact_bit(al, lambda w: var.get(w, ZERO), key)
            self.need(bit == ONE, 'fact %s is not a tautology' % (key[0],))
        else:
            bit = fact_bit(self.al, self._val_run, key)
            if bit != self.al.T:
                self._violation('establish', None, key[0], bit ^ self.al.T)
        self.facts.add(key)

    def admit(self, key, why):
        """a fact that follows from a lemma proved on paper whose numeric side conditions were just checked."""
        if self.prove:
            self.admitted.append(why)
        else:
            bit = fact_bit(self.al, self._val_run, key)
            if bit != self.al.T:
                self._violation('admit', None, why, bit ^ self.al.T)
        self.facts.add(key)

    def _val_run(self, w):
        if w is None:
            return 0
        if w == 'ONE':
            return self.al.T
        return self.store.get(w, 0)

    # ------------------------------------------------------------------ one segment
    def seg(self, label, k, pos, order, spec):
        self.nseg += 1
        if self.prove:
            self._seg_prove(label, k, pos, order, spec)
            self.pos = pos
        else:
            al, st = self.al, self.store
            posts, asserts, prefacts, postfacts = spec(al, self._val_run)
            for name, bit in asserts:
                if bit != al.T:
                    self._violation(label, k, name, bit ^ al.T)
            for key in prefacts:
                bit = fact_bit(al, self._val_run, key)
                if bit != al.T:
                    self._violation(label, k, 'prefact ' + key[0], bit ^ al.T)
            st.update(posts)
            for key in postfacts:
                bit = fact_bit(al, self._val_run, key)
                if bit != al.T:
                    self._violation(label, k, 'postfact ' + key[0], bit ^ al.T)
                self.facts.add(key)

    def _support(self, lo, hi):
        ops = self.ops
        sup = set()
        nlin = 0
        for i in range(5 * lo, 5 * hi, 5):
            kd = ops[i]
            sup.add(ops[i + 1])
            if kd != K_X:
                sup.add(ops[i + 2])
                if kd >= K_AND:
                    sup.add(ops[i + 3])
                    if ops[i + 4] >= 0:
                        sup.add(ops[i + 4])
                        nlin += 1
        return sup, nlin

    def _sim(self, B, pre, v, base, lo, hi):
        """symbolic run of ops[lo:hi] on the state v (dict wire -> edge, default base(wire)).
        -> None, or (op index, message, witness BDD) at the first violated obligation."""
        ops = self.ops
        AND, XOR = B.AND, B.XOR
        get = v.get
        for i in range(5 * lo, 5 * hi, 5):
            kd = ops[i]
            a = ops[i + 1]
            if kd == K_CX:
                b = ops[i + 2]
                va = get(a)
                if va is None:
                    va = base(a)
                vb = get(b)
                if vb is None:
                    vb = base(b)
                v[b] = XOR(va, vb)
            elif kd == K_X:
                va = get(a)
                if va is None:
                    va = base(a)
                v[a] = va ^ 1
            else:
                b, t, l = ops[i + 2], ops[i + 3], ops[i + 4]
                va = get(a)
                if va is None:
                    va = base(a)
                vb = get(b)
                if vb is None:
                    vb = base(b)
                vt = get(t)
                if vt is None:
                    vt = base(t)
                if kd == K_AND:
                    if vt != ZERO and AND(pre, vt) != ZERO:
                        return i // 5, 'AND target is not zero', AND(pre, vt)
                    v[t] = AND(va, vb)
                else:
                    c = AND(va, vb)
                    if l >= 0:
                        vl = get(l)
                        if vl is None:
                            vl = base(l)
                        c = XOR(c, vl)
                    c = XOR(c, vt)
                    if c != ZERO and AND(pre, c) != ZERO:
                        return i // 5, 'erase condition violated', AND(pre, c)
                    v[t] = ZERO
        return None

    def _seg_prove(self, label, k, hi, order, spec):
        lo = self.pos
        Z = self.Z
        B = BDD()
        al = BddAlg(B)
        sup, nlin = self._support(lo, hi)
        if sup & self.reserved:
            raise ProofFailure('segment touches the output register', segment=label, k=k)
        var = {}
        lev = {}
        oset = set(order)

        def newvar(w):
            lv, e = B.new_var()
            var[w] = e
            lev[lv] = w
        for w in sorted(sup):
            if w not in Z and w not in oset:
                newvar(w)
        for w in order:
            if isinstance(w, int) and w not in Z and w not in var:
                newvar(w)

        def val(w):
            if w is None:
                return ZERO
            if w == 'ONE':
                return ONE
            e = var.get(w)
            if e is None:
                if w in Z:
                    return ZERO
                raise ProofFailure('internal: contract of %r reads undeclared wire %r' % (label, w))
            return e

        def fail(msg, f=None, **kw):
            cex = None
            if f is not None:
                asg = B.pick_sat(f) or {}
                cex = {lev[lv]: bit for lv, bit in asg.items()}
            raise ProofFailure('segment %s (k=%s, ops %d..%d): %s' % (label, k, lo, hi, msg),
                               segment=label, k=k, lo=lo, hi=hi, entry_counterexample=cex, **kw)
        posts, asserts, prefacts, postfacts = spec(al, val)
        AND, XOR = B.AND, B.XOR
        pre = ONE
        for name, bit in asserts:
            if bit != ONE:
                self.assumed[name] = self.assumed.get(name, 0) + 1
                pre = AND(pre, bit)
        for key in prefacts:
            if key not in self.facts:
                fail('fact %s was not established' % key[0])
            pre = AND(pre, fact_bit(al, val, key))
        if pre == ZERO:
            fail('the precondition is unsatisfiable (vacuous contract)')
        v = {}
        r = self._sim(B, pre, v, val, lo, hi)
        if r is not None:
            fail('op %d: %s' % (r[0], r[1]), r[2], op=r[0])
        touched = sup | set(posts)
        exit_ = {}
        for w in touched:
            ent = val(w)
            exp = posts.get(w, ent)
            got = v.get(w, ent)
            if got != exp and AND(pre, XOR(got, exp)) != ZERO:
                fail('wire %d does not end with the value of the contract' % w, AND(pre, XOR(got, exp)), wire=w)
            exit_[w] = exp
        newz = [w for w in touched if exit_[w] == ZERO or AND(pre, exit_[w]) == ZERO]
        if self.facts:
            self.facts = {key for key in self.facts if touched.isdisjoint(fact_wires(key))}
        if postfacts:
            def valx(w):
                e = exit_.get(w)
                return val(w) if e is None else e
            for key in postfacts:
                bit = fact_bit(al, valx, key)
                if bit != ONE and AND(pre, bit ^ 1) != ZERO:
                    fail('postcondition fact %s does not hold' % key[0], AND(pre, bit ^ 1))
                self.facts.add(key)
        if self.clean:
            cnt = (hi - lo) + nlin
            mhi = self.nops - self.cinv
            mlo = mhi - cnt
            if mlo < self.F + self.n:
                fail('the inverse half is too short for the mirror of this segment')
            sup2, _ = self._support(mlo, mhi)
            if sup2 & self.reserved:
                fail('the mirror segment touches the output register')
            zafter = set(newz)

            def base2(w):
                e = var.get(w)
                if e is not None:
                    return e
                if (w in Z and w not in touched) or w in zafter:
                    return ZERO
                raise ProofFailure('segment %s (k=%s): the mirror ops %d..%d touch a foreign wire %d'
                                   % (label, k, mlo, mhi, w), segment=label, k=k)
            v2 = dict(exit_)
            r = self._sim(B, pre, v2, base2, mlo, mhi)
            if r is not None:
                fail('mirror op %d: %s' % (r[0], r[1]), r[2], op=r[0])
            for w in touched | sup2:
                ent = var.get(w, ZERO)
                got = v2.get(w)
                if got is None:
                    got = base2(w)
                if got != ent and AND(pre, XOR(got, ent)) != ZERO:
                    fail('mirror: wire %d does not return to its entry value' % w, AND(pre, XOR(got, ent)), wire=w)
            self.cinv += cnt
        for w in touched:
            Z.discard(w)
        Z.update(newz)
        nn = B.nodes()
        self.nodes += nn
        if nn > self.maxnodes:
            self.maxnodes = nn

    # ------------------------------------------------------------------ the chain: P_std
    def run(self):
        n, p = self.n, self.p
        need = self.need
        X = list(range(n))
        pos, info = self.mark('init')
        U, R, S = self.fresh_list(info['U']), self.fresh_list(info['R']), self.fresh_list(info['S'])
        need(len(U) == n and len(R) == n and len(S) == n, 'certificate: bad register lengths')

        def spec(al, val):
            posts = {U[i]: (al.T if (p >> i) & 1 else al.F) for i in range(n)}
            posts[S[0]] = al.T
            return posts, [], [], []
        self.seg('init', None, pos, U + [S[0]], spec)                  # u = p, v = x, r = 0, s = 1
        A, B, wprev = U, X, None
        for k in range(2 * n):
            if self.progress and k % 16 == 0:
                self.progress('step %d/%d' % (k, 2 * n))
            wa = min(n, max(1, 2 * n - 1 - k))
            wb = min(wa, max(1, (2 * n - k) // 2))
            self.assume_zero(A[wa:] + B[wa:], 'M1-width')              # u, v < 2^wa
            if k == 2 * n - 1:
                # M5: the state is (u, v) = (1, 0) with w_{2n-2} = 1: A = v is even, w = 1, d = 0, e = 0: the step
                # is A /= 2 (relabelling) and s *= 2 (relabelling) only
                pos, info = self.mark('last')
                self.assume_zero([A[0]], 'M5-last-even')

                def spec(al, val, A=A, B=B, wprev=wprev):
                    wv = wspec(al, [val(x) for x in A], [val(x) for x in B], val(wprev))
                    return {}, [('M5-last-w', wv)], [], []
                self.seg('last', k, pos, [wprev] + A + B, spec)
                A = A[1:] + A[:1]
                S = self._double(S, info.get('snew'))
                continue
            if k == 2 * n - 2:
                w, wkey = None, None
            else:
                pos, info = self.mark('cmp')
                w = self.fresh(info.get('w'))
                wkey = ('wdef', w, wprev, tuple(A), tuple(B))

                def spec(al, val, A=A, B=B, wprev=wprev, w=w, wkey=wkey):
                    wv = wspec(al, [val(x) for x in A], [val(x) for x in B], val(wprev))
                    return {w: wv}, [], [], [wkey]
                self.seg('cmp', k, pos, [wprev] + _inter(A[:wa], B[:wa]) + [w], spec)
            pos, info = self.mark('swap')

            def spec(al, val, A=A, B=B, R=R, S=S, wprev=wprev, w=w, wkey=wkey, k=k):
                a = [val(x) for x in A]
                b = [val(x) for x in B]
                r = [val(x) for x in R]
                s = [val(x) for x in S]
                wp = val(wprev)
                asserts, prefacts = [], []
                if w is None:
                    wv = al.T
                    asserts.append(('M5-w', wspec(al, a, b, wp)))      # w_{2n-2} = 1
                else:
                    wv = val(w)
                    prefacts.append(wkey)
                d = al.xor(wp, wv)
                nd = al.not_(d)
                if k >= 1:
                    asserts.append(('M2-P3', al.or_(nd, al.not_(al.xor(a[0], b[0])))))
                asserts.append(('M4-rs', al.or_(nd, bv.allzero(al, s[n:]))))
                posts = {}
                for i in range(n):
                    posts[A[i]] = al.ite(d, b[i], a[i])
                    posts[B[i]] = al.ite(d, a[i], b[i])
                    posts[R[i]] = al.ite(d, s[i], r[i])
                    posts[S[i]] = al.ite(d, r[i], s[i])
                return posts, asserts, prefacts, [('imp', (A[0], B[0])), ('geodd', tuple(A), tuple(B)),
                                                  ('nand', (R[0], S[0]))]
            self.seg('swap', k, pos, [wprev, w] + _inter(A, B, R, S[:n]) + S[n:], spec)
            pos, info = self.mark('sub')
            e = self.fresh(info.get('e'))

            def spec(al, val, A=A, B=B, e=e, wb=wb):
                a = [val(x) for x in A]
                b = [val(x) for x in B]
                ev = a[0]
                asserts = [('M3-vshort', al.or_(al.not_(ev), bv.allzero(al, b[wb:])))]
                diff, _ = bv.sub(al, a, [al.and_(ev, x) for x in b])
                posts = {A[i]: diff[i] for i in range(n)}
                posts[e] = ev
                return posts, asserts, [('imp', (A[0], B[0])), ('geodd', tuple(A), tuple(B))], []
            self.seg('sub', k, pos, _inter(A, B) + [e], spec)          # e = u_0; u -= e*v (exact: u >= v if e)
            self.require_zero(A[0], 'u is even after the subtraction')
            A = A[1:] + A[:1]                                          # u /= 2
            pos, info = self.mark('add')

            def spec(al, val, R=R, S=S, e=e):
                ev = val(e)
                r = [val(x) for x in R]
                s = [val(x) for x in S]
                asserts = [('M4-sfit', al.or_(al.not_(ev), bv.allzero(al, s[n:])))]
                sm, c = bv.add(al, r, [al.and_(ev, x) for x in s[:n]])
                asserts.append(('M4-rfit', al.not_(c)))
                return {R[i]: sm[i] for i in range(n)}, asserts, [('nand', (R[0], S[0]))], []
            self.seg('add', k, pos, [e] + _inter(R, S[:n]) + S[n:], spec)   # r += e*s (exact, fits n bits)
            S = self._double(S, info.get('snew'))                      # s *= 2
            wprev = w if w is not None else 'ONE'
        need(len(S) == 2 * n, 'internal: s register length')
        # frame of w_{2n-1} = 1: u = B, v = A, r = S, s = R
        self.u, self.v, self.sreg = B, A, R
        if self.red == 'generic':
            out = self._reduce_generic(S)
        else:
            out = self._reduce_special(S)
        self.out = out
        need(self.mi == len(self.marks), 'certificate: unused marks')
        if self.prove:
            need(self.pos == self.F, 'the forward part is not covered by the segments (%d of %d ops)'
                 % (self.pos, self.F))
            if self.clean:
                O = list(range(n, 2 * n))
                ops = self.ops
                need(self.F + n <= self.nops, 'no room for the copy')
                for i in range(n):
                    j = 5 * (self.F + i)
                    need(ops[j] == K_CX and ops[j + 1] == out[i] and ops[j + 2] == O[i] and O[i] in self.Z,
                         'copy of the result: op %d is not cx out[%d] -> q[%d]' % (self.F + i, i, O[i]))
                need(self.nops - self.cinv == self.F + n, 'the inverse half is not covered by the mirror segments')
        return self

    def _double(self, S, snew):
        n = self.n
        S = [self.fresh(snew)] + S
        if len(S) > 2 * n:
            top = S.pop()
            self.assume_zero([top], 'M6-rtop')                         # r < 2^k: bit 2n of the doubled register is 0
            self.owned.discard(top)
        return S

    def _reduce_generic(self, S):
        """long division of r (2n bits) by p: n+1 steps  W = 2*rem + next bit < 2p,  rem' = W mod p."""
        n, p = self.n, self.p
        pos, info = self.mark('red_init')
        ext = self.fresh_list(info.get('ext'))
        R2 = S + ext
        self.need(len(R2) == 2 * n + 1, 'certificate: bad extension of r')
        self.need((1 << (n - 1)) < p < (1 << n), 'p must have exactly n bits')
        self.establish(('lt', tuple(R2[n + 1:2 * n + 1]), p))
        for j in range(n, -1, -1):
            win = R2[j:j + n + 1]
            pos, info = self.mark('div')
            f = self.fresh(info.get('f'))

            def spec(al, val, win=win, f=f):
                wv = [val(x) for x in win]
                d, bor = bv.sub(al, wv, bv.const(al, p, n + 1))
                posts = {win[i]: al.ite(bor, wv[i], d[i]) for i in range(n)}
                posts[win[n]] = al.F
                posts[f] = bor
                return posts, [], [('lt', tuple(win[1:]), p)], [('lt', tuple(win[:n]), p)]
            self.seg('div', j, pos, win + [f], spec)
        return R2[:n]

    def _reduce_special(self, S):
        """folding rounds (literal modular additions; exactness by the fold lemma, PROOF.md, with its numeric side
        conditions checked here), then the binary search of the quotient."""
        n, p = self.n, self.p
        need = self.need
        m = n
        cur, Bd = S, 1 << (2 * n)
        self.establish(('lt', tuple(cur), Bd))
        while self.peek() == 'sp_const':
            pos, info = self.mark('sp_const')
            K = info.get('K')
            Aw = self.fresh_list(info.get('A'))
            L = len(Aw)
            need(isinstance(K, int) and 0 <= K < (1 << L) and L >= m, 'certificate: bad fold constant')
            need(('lt', tuple(cur), Bd) in self.facts, 'internal: bound fact lost')
            Hmax = (Bd - 1) >> m
            Lw, Hw = cur[:m], cur[m:][:Hmax.bit_length()]

            def spec(al, val, Aw=Aw, K=K):
                return {Aw[i]: (al.T if (K >> i) & 1 else al.F) for i in range(len(Aw))}, [], [], []
            self.seg('sp_const', None, pos, Aw, spec)
            pos, info = self.mark('sp_addL')

            def spec(al, val, Aw=Aw, Lw=Lw):
                a = [val(x) for x in Aw]
                l = [val(x) for x in Lw] + [al.F] * (len(Aw) - len(Lw))
                sm, _ = bv.add(al, a, l)
                return {Aw[i]: sm[i] for i in range(len(Aw))}, [], [], []
            self.seg('sp_addL', None, pos, _inter(Aw, Lw), spec)
            terms = []
            while self.peek() == 'sp_term':
                pos, info = self.mark('sp_term')
                s_, e_ = info.get('s'), info.get('e')
                need(s_ in (1, -1) and isinstance(e_, int) and 0 <= e_ < L, 'certificate: bad fold term')
                sub_, h = Aw[e_:], Hw[:L - e_]

                def spec(al, val, sub_=sub_, h=h, s_=s_):
                    a = [val(x) for x in sub_]
                    hv = [val(x) for x in h] + [al.F] * (len(sub_) - len(h))
                    res = bv.add(al, a, hv)[0] if s_ > 0 else bv.sub(al, a, hv)[0]
                    return {sub_[i]: res[i] for i in range(len(sub_))}, [], [], []
                self.seg('sp_term', e_, pos, _inter(sub_, h), spec)
                terms.append((s_, e_))
            neg = sum(Hmax << e_ for s_, e_ in terms if s_ < 0)
            pos_ = sum(Hmax << e_ for s_, e_ in terms if s_ > 0)
            Bn = K + (1 << m) + pos_
            need((sum(s_ << e_ for s_, e_ in terms) - (1 << m)) % p == 0, 'fold lemma: terms are not 2^m mod p')
            need(K % p == 0 and K >= neg and Bn <= (1 << L), 'fold lemma: side conditions fail')
            self.admit(('lt', tuple(Aw), Bn), 'fold lemma')
            cur, Bd = Aw, Bn
        while self.peek() == 'sp_q':
            pos, info = self.mark('sp_q')
            f = self.fresh(info.get('f'))
            Kq = info.get('K')
            L = len(cur)
            need(isinstance(Kq, int) and 0 < Kq < (1 << L) and Kq % p == 0, 'certificate: bad quotient constant')
            Bn = max(Kq, Bd - Kq) if Bd > Kq else Bd

            def spec(al, val, cur=cur, f=f, Kq=Kq, Bd=Bd, Bn=Bn):
                a = [val(x) for x in cur]
                ge = al.not_(bv.lt_const(al, a, Kq))
                d, _ = bv.sub(al, a, bv.const(al, Kq, len(a)))
                posts = {cur[i]: al.ite(ge, d[i], a[i]) for i in range(len(a))}
                posts[f] = ge
                return posts, [], [('lt', tuple(cur), Bd)], [('lt', tuple(cur), Bn)]
            self.seg('sp_q', None, pos, cur + [f], spec)
            Bd = Bn
        need(Bd <= p, 'special reduction: the final bound %d bits is not <= p' % Bd.bit_length())
        return cur[:m]


def simulate_flat(lifted, init, S):
    """bit-sliced classical run of a lifted program on S samples; init {wire: mask}.
    -> (values list, None) or (values, (op index, message, mask of violating samples))."""
    M = (1 << S) - 1
    v = [0] * lifted.nq
    for w, mk in init.items():
        v[w] = mk
    ops = lifted.ops
    for i in range(0, len(ops), 5):
        kd = ops[i]
        a = ops[i + 1]
        if kd == K_CX:
            v[ops[i + 2]] ^= v[a]
        elif kd == K_X:
            v[a] ^= M
        elif kd == K_AND:
            t = ops[i + 3]
            if v[t]:
                return v, (i // 5, 'AND target not zero', v[t])
            v[t] = v[a] & v[ops[i + 2]]
        else:
            t, l = ops[i + 3], ops[i + 4]
            c = (v[a] & v[ops[i + 2]]) ^ v[t]
            if l >= 0:
                c ^= v[l]
            if c:
                return v, (i // 5, 'erase condition violated', c)
            v[t] = 0
    return v, None
