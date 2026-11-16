# SPDX-License-Identifier: Apache-2.0
"""compose.walker -- the generic engine of the compositional proofs: segments, contracts, facts, zero wires.

A client (compose.kal for the Kaliski inversion) subclasses Walker and walks its algorithm segment by segment:

    pos, info = self.mark(label)            next mark of the (untrusted) certificate: end position, wire names
    w = self.fresh(info['w'])               a wire the certificate calls freshly allocated (checked: free and zero)
    self.assume_zero(wires, name)           a mathematical assumption "these register bits are 0 here"
    self.seg(label, k, pos, order, spec, subst)     one segment with its contract (see compose.kal)
    self.establish(key) / self.admit(key, why)      facts: tautologies / consequences of paper lemmas
    self.finish(out)                        coverage of all ops; clean mode: the copy and the mirror tiling

Two modes: prove (ROBDD proof of every segment on the real lifted ops) and run (bit-sliced execution of the contracts
on concrete inputs: the test of the mathematical assumptions). See compose/PROOF.md, sections 2 and 3.

Facts are predicates over wires, registered by name with register_fact(name, wires_of_key, bit_of_key); the built-in
ones are ('lt', wires, c): value < c;  ('nand', (a, b));  ('imp', (a, b)).
"""
import random

from ..bdd import BDD, ONE, ZERO
from . import bv
from .bv import IntAlg, BddAlg
from .lift import K_X, K_CX, K_AND

__all__ = ['Walker', 'ProofFailure', 'register_fact', 'fact_bit', 'fact_wires', 'simulate_flat', 'inter']


class ProofFailure(Exception):
    def __init__(self, msg, **detail):
        super().__init__(msg)
        self.detail = detail


def inter(*lists):
    """interleave lists by index (the ROBDD variable order of registers that are combined bitwise)."""
    out = []
    for i in range(max(len(x) for x in lists)):
        for x in lists:
            if i < len(x):
                out.append(x[i])
    return out


FACTS = {}


def register_fact(name, wires, bit):
    """a kind of fact: wires(key) -> the wires it speaks about; bit(al, val, key) -> its truth value."""
    FACTS[name] = (wires, bit)


def fact_wires(key):
    return FACTS[key[0]][0](key)


def fact_bit(al, val, key):
    return FACTS[key[0]][1](al, val, key)


register_fact('lt', lambda key: key[1], lambda al, val, key: bv.lt_const(al, [val(w) for w in key[1]], key[2]))
register_fact('nand', lambda key: key[1], lambda al, val, key: al.not_(al.and_(val(key[1][0]), val(key[1][1]))))
register_fact('imp', lambda key: key[1], lambda al, val, key: al.or_(al.not_(val(key[1][0])), val(key[1][1])))


def _models(B, f, rnd, count):
    """`count` random total assignments {level: bit} satisfying f (f != ZERO)."""
    var, hi, lo = B.var, B.hi, B.lo
    out = []
    nv = B.nvars
    for _ in range(count):
        bits = rnd.getrandbits(nv + 1)
        asg = {lv: (bits >> lv) & 1 for lv in range(nv)}
        e = f
        while e != ONE:
            n, c = e >> 1, e & 1
            h, l = hi[n] ^ c, lo[n] ^ c
            if h == ZERO:
                b = 0
            elif l == ZERO:
                b = 1
            else:
                b = asg[var[n]]
            asg[var[n]] = b
            e = h if b else l
        out.append(asg)
    return out


class Walker:
    def __init__(self, hints, inputs, lifted=None, init=None, labels=None, clean=False, reserved=(), progress=None):
        """inputs: the wires that carry the input (every other wire starts in |0>). prove mode: lifted (compose.lift).
        run mode: init {input wire: mask over the samples}, labels = the samples (for violation reports).
        clean: the circuit is forward + copy to `reserved` + inverse; hints['fwd_len'] = length of the forward part."""
        self.clean = bool(clean)
        self.marks = hints['marks']
        self.mi = 0
        self.prove = lifted is not None
        self.progress = progress
        self.owned = set(inputs)
        self.reserved = set(reserved)
        self.owned |= self.reserved
        self.facts = set()
        self.assumed = {}
        self.admitted = []
        self.nseg = 0
        self.nodes = 0
        self.maxnodes = 0
        self.rnd = random.Random(20261116)
        if self.prove:
            self.ops, self.nops, self.nq = lifted.ops, lifted.n, lifted.nq
            self.pos = 0
            self.cinv = 0
            self.Z = set(range(lifted.nq)) - set(inputs)
            self.F = hints.get('fwd_len', lifted.n) if clean else lifted.n
            if not (isinstance(self.F, int) and 0 <= self.F <= self.nops):
                raise ProofFailure('certificate: bad forward length')
        else:
            self.nq = hints['nqubits']
            self.xs = list(labels)
            self.al = IntAlg(len(self.xs))
            self.store = dict(init)
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
    def seg(self, label, k, pos, order, spec, subst=()):
        """one segment. subst: preconditions given as PARAMETRISATIONS: tuples (wire, expr, name, key, pred) meaning
        "the precondition pred (a mathematical assumption `name`, or the established fact `key`) holds, and under it
        the entry value of `wire` equals expr(al, val)". In prove mode the wire's entry value is replaced by expr
        (after checking that this parametrises exactly the states satisfying the preconditions), which makes the
        circuit and the contract equal as ROBDDs instead of equal under a side condition."""
        self.nseg += 1
        if self.prove:
            self._seg_prove(label, k, pos, order, spec, subst)
            self.pos = pos
        else:
            al, st = self.al, self.store
            for w, expr, name, key, pred in subst:
                bit = pred(al, self._val_run)
                if bit != al.T:
                    self._violation(label, k, name or ('prefact ' + key[0]), bit ^ al.T)
            posts, asserts, prefacts, postfacts, proved = spec(al, self._val_run)
            for name, bit in asserts + proved:
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

    def _seg_prove(self, label, k, hi, order, spec, subst=()):
        lo = self.pos
        Z = self.Z
        B = BDD()
        al = BddAlg(B)
        sup, nlin = self._support(lo, hi)
        if sup & self.reserved:
            raise ProofFailure('segment touches the output register', segment=label, k=k)
        var = {}
        lev = {}
        lev_of_lit = set()
        oset = set(order)

        def newvar(w):
            lv, e = B.new_var()
            var[w] = e
            lev[lv] = w
            lev_of_lit.add(e >> 1)
        for w in sorted(sup):
            if w not in Z and w not in oset:
                newvar(w)
        for w in order:
            if isinstance(w, int) and w not in Z and w not in var:
                newvar(w)

        def raw(w):
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
        val = raw
        AND, XOR = B.AND, B.XOR
        if subst:
            # entry values under the parametrised preconditions; checked: (1) every parametrised precondition holds
            # identically on the substituted values, (2) every raw state that satisfies all of them is a fixed point
            # of the substitution. Hence "for all variables, on the substituted values" = "for all entry states
            # satisfying these preconditions".
            sub = {}
            for w, expr, name, key, pred in subst:
                if w in sub or not isinstance(w, int):
                    raise ProofFailure('internal: bad substitution in %r' % label)
                if w not in Z:
                    sub[w] = expr(al, raw)

            def val(w):
                e = sub.get(w)
                return raw(w) if e is None else e
            phi = ONE
            used = set()
            for w, expr, name, key, pred in subst:
                if key is not None and key not in self.facts:
                    raise ProofFailure('segment %s (k=%s): fact %s was not established' % (label, k, key[0]),
                                       segment=label, k=k)
                if pred(al, val) != ONE:
                    raise ProofFailure('internal: substitution of %r does not imply its precondition' % label)
                pr = pred(al, raw)
                if name is not None and pr != ONE and name not in used:
                    used.add(name)
                    self.assumed[name] = self.assumed.get(name, 0) + 1
                phi = AND(phi, pr)
            for w, e in sub.items():
                if AND(phi, XOR(e, raw(w))) != ZERO:
                    raise ProofFailure('internal: substitution of %r is not a parametrisation' % label)

        def fail(msg, f=None, **kw):
            cex = None
            if f is not None:
                asg = B.pick_sat(f) or {}
                cex = {lev[lv]: bit for lv, bit in asg.items()}
            raise ProofFailure('segment %s (k=%s, ops %d..%d): %s' % (label, k, lo, hi, msg),
                               segment=label, k=k, lo=lo, hi=hi, entry_counterexample=cex, **kw)
        posts, asserts, prefacts, postfacts, proved = spec(al, val)
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
        for name, bit in proved:
            if bit != ONE and AND(pre, bit ^ 1) != ZERO:
                fail('the contract property "%s" does not hold' % name, AND(pre, bit ^ 1))
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
        # wires that are 0 at the exit under the precondition. Cheap pre-test: if the exit value is 1 on one of a few
        # random models of pre the wire is not zero; only otherwise the (exact) conjunction is built
        models = _models(B, pre, self.rnd, 6) if pre != ONE else None
        ev = B.eval
        newz = []
        for w in touched:
            f = exit_[w]
            if f == ZERO:
                newz.append(w)
            elif f == ONE or f >> 1 in lev_of_lit:
                continue
            elif models is None:
                continue
            elif not any(ev(f, m) for m in models) and AND(pre, f) == ZERO:
                newz.append(w)
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
                if w in var:
                    return val(w)
                if (w in Z and w not in touched) or w in zafter:
                    return ZERO
                raise ProofFailure('segment %s (k=%s): the mirror ops %d..%d touch a foreign wire %d'
                                   % (label, k, mlo, mhi, w), segment=label, k=k)
            v2 = dict(exit_)
            r = self._sim(B, pre, v2, base2, mlo, mhi)
            if r is not None:
                fail('mirror op %d: %s' % (r[0], r[1]), r[2], op=r[0])
            for w in touched | sup2:
                ent = val(w) if (w in var or w in Z) else ZERO
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

    def finish(self, out=()):
        """all marks used; prove mode: every forward op is covered by a segment; clean mode: the next len(out) ops
        are cx out[i] -> reserved[i] and the mirror segments cover exactly the rest."""
        need = self.need
        need(self.mi == len(self.marks), 'certificate: unused marks')
        if self.prove:
            need(self.pos == self.F, 'the forward part is not covered by the segments (%d of %d ops)'
                 % (self.pos, self.F))
            if self.clean:
                O = sorted(self.reserved)
                m = len(O)
                ops = self.ops
                need(m == len(out) and self.F + m <= self.nops, 'no room for the copy of the result')
                for i in range(m):
                    j = 5 * (self.F + i)
                    need(ops[j] == K_CX and ops[j + 1] == out[i] and ops[j + 2] == O[i] and O[i] in self.Z,
                         'copy of the result: op %d is not cx out[%d] -> q[%d]' % (self.F + i, i, O[i]))
                need(self.nops - self.cinv == self.F + m, 'the inverse half is not covered by the mirror segments')
        return self


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
