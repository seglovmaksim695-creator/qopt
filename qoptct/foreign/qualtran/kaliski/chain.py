"""chain -- a Walker client (modelled on qoptct compose.kal) for the Qualtran circuit KaliskiModInverse(N, p).

Segments: 'init' (loaders f=1, u=p, s=1); for each of the 2N iterations the six real Qualtran steps S1..S6;
'fin' (NOT r, AddK(p+1), XorK, three controlled XorK).  Every segment has a CONTRACT spec(al, val):
    posts      values of the output wires of the step as functions of the input wires (literally from the code of the
               Qualtran steps, arithmetic mod 2^N); wires the step releases -> 0
    asserts    MATHEMATICAL assumptions on the input state (proved separately by the Z3 lemmas for all N, p, x:
               lemmas_z3.py): M-carry, M-dbl, M-a, M-s, M-end; plus assume_zero 'M-end-vf'
    prefacts   facts PROVED (ROBDD) at the output of the previous segments (R1..R3, Q1, P1..P4, nand)
    postfacts  facts proved at the output
    proved     properties proved under the precondition (no borrow; a and b cleared)
Mode prove: every segment is proved by the unchanged Walker on the real lifted operations (ROBDD, all inputs at once).
Mode run: the contracts are executed bit-parallel on a sample of x; all asserts/facts and the result are checked.
The certificate (segment positions, wire names) is taken from the expansion and is NOT trusted: a wrong one fails the proof.
"""
from qoptct.verify.compose import bv
from qoptct.verify.compose.walker import Walker, ProofFailure, register_fact, inter

__all__ = ['QChain', 'ProofFailure', 'build_marks', 'RELS', 'FlowChecker', 'check_flow', 'expected_flow_counts', 'FLOW_VERSION']


def _halves(x):
    h = len(x) // 2
    return x[:h], x[h:]


def _gt(al, u, v):
    return al.not_(bv.le(al, u, v))


def _q1(al, x):          # a = f & (~u0 | (u0 & v0 & [u > v]))
    a, f = x[0], x[1]
    u, v = _halves(x[2:])
    rhs = al.and_(f, al.or_(al.not_(u[0]), al.and_(al.and_(u[0], v[0]), _gt(al, u, v))))
    return al.not_(al.xor(a, rhs))


def _p1(al, x):          # f & ~b -> u <= v
    b, f = x[0], x[1]
    u, v = _halves(x[2:])
    return al.or_(al.not_(al.and_(f, al.not_(b))), bv.le(al, u, v))


RELS = {
    'R1': lambda al, x: al.not_(al.xor(x[0], al.xor(x[1], x[2]))),                         # (b, a, m): b = a ^ m
    'R3': lambda al, x: al.not_(al.xor(x[0], al.and_(x[1], al.not_(x[2])))),               # (a, f, u0): a = f & ~u0
    'R2': lambda al, x: al.or_(al.not_(x[1]), al.not_(al.xor(x[0], al.not_(al.and_(x[2], x[3]))))),
    #                                                                  (b, f, u0, v0): f -> (b = ~(u0 & v0))
    'Q1': _q1,                                                                             # (a, f, U, V)
    'P1': _p1,                                                                             # (b, f, U, V)
    'P2': lambda al, x: al.or_(al.not_(al.and_(x[1], al.not_(x[0]))), al.and_(x[2], x[3])),  # (b, f, u0, v0)
    'P3': lambda al, x: al.or_(al.not_(al.and_(x[1], x[0])), al.not_(x[2])),               # (b, f, v0): f&b -> ~v0
    'P4': lambda al, x: al.or_(al.not_(x[0]), al.not_(x[1])),                              # (f, v0): f -> ~v0
}
register_fact('qrel', lambda key: key[2], lambda al, val, key: RELS[key[1]](al, [val(w) for w in key[2]]))


def rel(name, *wires):
    flat = []
    for w in wires:
        flat.extend(w if isinstance(w, (list, tuple)) else [w])
    return ('qrel', name, tuple(flat))


def _lsb(msb):
    return list(reversed([int(q) for q in msb]))


def build_marks(C):
    """the certificate: [(end position, label, info)], info = wire lists (LSB first) at entry/exit of the segment."""
    P, n = C.piece, C.n
    marks = []
    g0 = C.gmaps[0]
    u0 = [int(g0[w]) for w in _lsb(P.ins['u'])]
    s0 = [int(g0[w]) for w in _lsb(P.ins['s'])]
    f0 = int(g0[P.ins['f'][0]])
    marks.append((C.iter_start[0], 'init', dict(u=u0, s=s0, f=f0)))
    for k in range(2 * n):
        g = C.gmaps[k]
        base = C.iter_start[k]
        ab = [int(g[w]) for w in P.allocs]
        for st in P.steps:
            ins = {r: [int(g[w]) for w in _lsb(ws)] for r, ws in st['ins'].items()}
            outs = {r: [int(g[w]) for w in _lsb(ws)] for r, ws in st['outs'].items()}
            marks.append((base + st['hi'], st['name'][-5:].replace('Step', 'S'), dict(ins=ins, outs=outs, ab=ab)))
    gl = C.gmaps[-1]
    last = {r: [int(gl[w]) for w in _lsb(ws)] for r, ws in P.outs.items()}
    t0 = int(C.gmaps[0][P.outs['is_terminal'][0]])
    marks.append((C.nops, 'fin', dict(u=last['u'], v=last['v'], r=last['r'], s=last['s'], f=last['f'][0], t0=t0)))
    return marks


FLOW_VERSION = 'flow-1'
_STEP_REGS = {'S1': ('v', 'm', 'f', 'is_terminal'), 'S2': ('u', 'v', 'a', 'b', 'm', 'f'), 'S3': ('u', 'v', 'a', 'b', 'm', 'f'),
              'S4': ('u', 'v', 'r', 's', 'a'), 'S5': ('u', 'v', 'r', 's', 'b', 'f'),
              'S6': ('u', 'v', 'r', 's', 'a', 'b', 'm', 'f')}
_STEP_NEW = {'S1': ('m', 'is_terminal'), 'S2': ('a', 'b'), 'S3': (), 'S4': (), 'S5': (), 'S6': ()}   # r: new in iteration 0, S4


class FlowChecker:
    """WIRE FLOW (closes the hole AUD-D-01).  The certificate (marks) gives the wires of the registers of every step anew; this class keeps
    the CURRENT wires of every register and checks:
      * the input of every step = the current wires of the register (u, v, r, s, f, a, b, m, t) -- between steps and between iterations;
      * m, t (S1), a, b (S2) and r (S4 of iteration 0) -- only NEW wires (not occupied by another register; m, t, r and the initial
        u, s, f never used before; that they are 0 is checked by Walker.fresh); a, b cannot be junk wires;
      * the wires of one step and of its outputs are pairwise distinct, the outputs do not touch foreign registers;
      * the end: u, v, r, s, f = the current ones, t0 = the t wire of iteration 0;
      * outputs['junk'] = the t/m wires in the documented order (junk = M*2^(2N) + T; most significant bit first:
        m_0..m_(2N-1), t_0..t_(2N-1); the wire list is least significant bit first).
    Any mismatch -> ProofFailure('wire flow: ...').  The same class works inside the Walker client (QChain._begin/_fin)
    and in the preliminary static pass over the whole certificate (check_flow)."""

    def __init__(self, n, xwires):
        self.n = n
        self.cur = {'v': list(xwires)}
        self.used = set(xwires)
        self.tw, self.mw = [], []
        self.k = -1
        self.eq_checks = 0
        self.new_checks = 0
        self.disjoint_checks = 0
        self.junk_ok = False
        self.fin_ok = False

    def _fail(self, msg):
        raise ProofFailure('wire flow: ' + msg)

    def _width(self, reg):
        return self.n if reg in ('u', 'v', 'r', 's') else 1

    def _live(self, except_regs=()):
        return {w for r, ws in self.cur.items() if r not in except_regs for w in ws}

    def _wires(self, reg, ws, where):
        if not isinstance(ws, (list, tuple)) or len(ws) != self._width(reg) or \
                not all(isinstance(w, int) and not isinstance(w, bool) and w >= 0 for w in ws) or len(set(ws)) != len(ws):
            self._fail('%s: register %s: wrong wire list' % (where, reg))
        return [int(w) for w in ws]

    def _new(self, reg, ws, where, strict):
        ws = self._wires(reg, ws, where)
        if reg in self.cur:
            self._fail('%s: register %s already exists but must be new' % (where, reg))
        live = self._live()
        if any(w in live for w in ws):
            self._fail('%s: new register %s occupies a wire of another register' % (where, reg))
        if any(w in self.tw or w in self.mw for w in ws):
            self._fail('%s: new register %s occupies a junk wire' % (where, reg))
        if strict and any(w in self.used for w in ws):
            self._fail('%s: new register %s uses an "old" wire' % (where, reg))
        self.cur[reg] = ws
        self.used.update(ws)
        self.new_checks += 1

    def _eq(self, reg, ws, where):
        ws = self._wires(reg, ws, where)
        if reg not in self.cur:
            self._fail('%s: register %s does not exist yet' % (where, reg))
        if self.cur[reg] != ws:
            self._fail('%s: input of register %s differs from its current wires (the input of the step is not the output of the previous one)' % (where, reg))
        self.eq_checks += 1

    def feed(self, label, info):
        k = self.k
        if label == 'init':
            where = 'init'
            for reg in ('u', 's'):
                self._new(reg, info[reg], where, True)
            self._new('f', [info['f']], where, True)
            return
        if label == 'fin':
            return self._fin(info)
        if label not in _STEP_REGS:
            self._fail('unknown label %r' % (label,))
        if label == 'S1':
            self.k = k = k + 1
            for r in ('m', 'is_terminal', 'a', 'b'):
                self.cur.pop(r, None)
        where = '%s iteration %d' % (label, k)
        ins, outs = info['ins'], info['outs']
        exp = set(_STEP_REGS[label])
        if set(ins) != exp or set(outs) != exp:
            self._fail('%s: registers of the step %r/%r, expected %r' % (where, sorted(ins), sorted(outs), sorted(exp)))
        allin = [w for ws in ins.values() for w in ws]
        if len(set(allin)) != len(allin):
            self._fail('%s: input wires of the step repeat' % where)
        new = set(_STEP_NEW[label])
        if label == 'S4':
            if (k == 0) != ('r' not in self.cur):
                self._fail('%s: r must appear exactly in iteration 0' % where)
            if k == 0:
                self._new('r', ins['r'], where, True)
        for reg in sorted(ins):
            if reg == 'r' and label == 'S4' and k == 0:
                continue
            if reg in new:
                self._new(reg, ins[reg], where, reg in ('m', 'is_terminal'))
            else:
                self._eq(reg, ins[reg], where)
        # outputs: width, pairwise distinct, do not touch registers that are not in the step
        for reg in outs:
            self._wires(reg, outs[reg], where)
        allout = [w for ws in outs.values() for w in ws]
        if len(set(allout)) != len(allout):
            self._fail('%s: output wires of the step repeat' % where)
        if set(allout) & self._live(except_regs=tuple(exp)):
            self._fail('%s: an output of the step touches a wire of a register that is not in the step' % where)
        self.disjoint_checks += 1
        for reg in outs:
            self.cur[reg] = [int(w) for w in outs[reg]]
        self.used.update(allout)
        if label == 'S1':
            self.tw.append(self.cur['is_terminal'][0])
        if label == 'S6':
            self.mw.append(self.cur['m'][0])
            for r in ('a', 'b'):
                self.cur.pop(r)

    def _fin(self, info):
        n = self.n
        where = 'end'
        if self.k != 2 * n - 1 or len(self.tw) != 2 * n or len(self.mw) != 2 * n:
            self._fail('%s: number of iterations %d, expected %d' % (where, self.k + 1, 2 * n))
        for reg in ('u', 'v', 'r', 's'):
            self._eq(reg, info[reg], where)
        self._eq('f', [info['f']], where)
        if info['t0'] != self.tw[0] or not isinstance(info['t0'], int):
            self._fail('%s: t0 differs from the t wire of iteration 0' % where)
        self.eq_checks += 1
        self.fin_ok = True

    def expected_junk(self):
        n = self.n
        ts, ms = self.tw, self.mw
        return [ts[2 * n - 1 - i] if i < 2 * n else ms[4 * n - 1 - i] for i in range(4 * n)]

    def check_outputs(self, outputs):
        if not self.fin_ok:
            self._fail('the end was not checked')
        exp = self.expected_junk()
        got = [int(w) for w in outputs['junk']]
        if got != exp:
            bad = [i for i in range(min(len(got), len(exp))) if got[i] != exp[i]][:5]
            self._fail('outputs[junk] differs from the t/m wires in the documented order (lengths %d/%d, first differences %r)'
                       % (len(got), len(exp), bad))
        if len(set(exp)) != len(exp):
            self._fail('junk wires repeat')
        xo = [int(w) for w in outputs['x']]
        if xo != self.cur['r'] or set(xo) & set(exp):
            self._fail('outputs[x] differs from the r wires or meets the junk')
        self.junk_ok = True

    def summary(self):
        return dict(version=FLOW_VERSION, ok=bool(self.fin_ok and self.junk_ok), iterations=self.k + 1,
                    eq_checks=self.eq_checks, new_checks=self.new_checks, disjoint_checks=self.disjoint_checks,
                    junk_ok=self.junk_ok)


def expected_flow_counts(n):
    """how many flow checks a full pass must make (compared by the verdict): 31 comparisons per iteration (minus one: r is new
    in iteration 0), 4 new registers (m, t, a, b); r and the initial u, s, f -- 4 more new ones; the end: 5 comparisons + t0."""
    return dict(iterations=2 * n, eq_checks=31 * 2 * n - 1 + 6, new_checks=4 * 2 * n + 1 + 3, disjoint_checks=6 * 2 * n)


def check_flow(C, marks):
    """static pass over the whole certificate (fast, before the proof): wire flow + outputs."""
    fc = FlowChecker(C.n, list(C.inputs['x']))
    for pos, lab, info in marks:
        fc.feed(lab, info)
    fc.check_outputs(C.outputs)
    return fc


class QChain(Walker):
    def __init__(self, C, lifted=None, xs=None, progress=None, marks=None):
        self.C, self.n, self.p = C, C.n, C.p
        X = list(C.inputs['x'])
        self.flow = FlowChecker(C.n, X)
        init = None
        if lifted is None:
            xs = list(xs)
            init = {}
            for i, w in enumerate(X):
                mk = 0
                for j, x in enumerate(xs):
                    if (x >> i) & 1:
                        mk |= 1 << j
                init[w] = mk
        hints = dict(marks=marks if marks is not None else build_marks(C), nqubits=C.nq)
        Walker.__init__(self, hints, X, lifted=lifted, init=init, labels=xs, progress=progress)
        self.trace = []

    # ------------------------------------------------------------------ helpers
    def _own_flow(self, ins, outs):
        """wires entering a contract that nobody owns must be fresh zero wires; wires that appear only at the exit
        are fresh; wires that disappear are released (they must end at 0: their post is 0)."""
        iw = {w for ws in ins.values() for w in ws}
        ow = {w for ws in outs.values() for w in ws}
        for w in sorted(iw | ow):
            if w not in self.owned:
                self.fresh(w)
        return sorted(iw - ow)

    def _posts(self, al, outs, vals, released):
        posts = {}
        for r, ws in outs.items():
            vs = vals[r]
            self.need(len(vs) == len(ws), 'internal: width of %s' % r)
            for w, x in zip(ws, vs):
                posts[w] = x
        for w in released:
            self.need(w not in posts, 'internal: released wire is an output')
            posts[w] = al.F
        return posts

    # ------------------------------------------------------------------ the chain
    def run(self):
        n, p = self.n, self.p
        check_flow(self.C, self.marks)           # static pass over the whole certificate: fails fast, before the proof
        pos, info = self.mark('init')
        self.flow.feed('init', info)
        U, S, F = info['u'], info['s'], info['f']
        self.need(len(U) == n and len(S) == n, 'certificate: init widths')
        for w in U + S + [F]:
            self.fresh(w)

        def spec(al, val):
            posts = {U[i]: al.xor(val(U[i]), al.T if (p >> i) & 1 else al.F) for i in range(n)}
            posts.update({S[i]: al.xor(val(S[i]), al.T if i == 0 else al.F) for i in range(n)})
            posts[F] = al.not_(val(F))
            return posts, [], [], [], []
        self.seg('init', None, pos, [F] + inter(U, S), spec)
        for k in range(2 * n):
            if self.progress and k % 16 == 0:
                self.progress('iteration %d/%d' % (k, 2 * n))
            self._iteration(k)
        self._fin()
        self.finish()
        self._final_zero_check()
        return self

    def _begin(self, label):
        pos, info = self.mark(label)
        self.flow.feed(label, info)              # input of the step = current wires of the registers; new wires are new
        ins, outs = info['ins'], info['outs']
        rel_ = self._own_flow(ins, outs)
        one = lambda d: {r: (ws[0] if len(ws) == 1 else ws) for r, ws in d.items()}
        return pos, info, one(ins), one(outs), outs, rel_

    def _end(self, rel_):
        for w in rel_:
            self.require_zero(w, 'wire released by a step')
            self.owned.discard(w)

    def _iteration(self, k):
        n, p = self.n, self.p
        # ---- S1 (Step1): q = f & [v = 0]; m ^= q; f ^= m; t ^= m
        pos, info, i, o, outs, rel_ = self._begin('S1')
        V, M, Fw, T = i['v'], i['m'], i['f'], i['is_terminal']

        def spec(al, val):
            v = [val(w) for w in V]
            f, m, t = val(Fw), val(M), val(T)
            q = al.and_(f, bv.allzero(al, v))
            m1 = al.xor(m, q)
            f1 = al.xor(f, m1)
            vals = dict(v=v, m=[m1], f=[f1], is_terminal=[al.xor(t, m1)])
            return self._posts(al, outs, vals, rel_), [], [], [('nand', (o['f'], o['m']))], []
        self.seg('S1', k, pos, [Fw, M, T] + V, spec)
        self._end(rel_)
        T1 = o['is_terminal']
        # ---- S2 (Step2): a ^= f & ~u0; m ^= f & ~v0 & ~a; b ^= a; b ^= m
        pos, info, i, o, outs, rel_ = self._begin('S2')
        U, V, A, B, M, Fw = i['u'], i['v'], i['a'], i['b'], i['m'], i['f']
        self.need(sorted(info['ab']) == sorted([A, B]), 'certificate: a, b')

        def spec(al, val):
            u = [val(w) for w in U]
            v = [val(w) for w in V]
            a, b, m, f = val(A), val(B), val(M), val(Fw)
            a1 = al.xor(a, al.and_(f, al.not_(u[0])))
            m1 = al.xor(m, al.and_(al.and_(f, al.not_(v[0])), al.not_(a1)))
            b1 = al.xor(al.xor(b, a1), m1)
            vals = dict(u=u, v=v, a=[a1], b=[b1], m=[m1], f=[f])
            post = [rel('R1', o['b'], o['a'], o['m']), rel('R3', o['a'], o['f'], o['u'][0]),
                    rel('R2', o['b'], o['f'], o['u'][0], o['v'][0])]
            return self._posts(al, outs, vals, rel_), [], [('nand', (Fw, M))], post, []
        self.seg('S2', k, pos, [Fw, M, A, B] + inter(U, V), spec)
        self._end(rel_)
        # ---- S3 (Step3): c = [u > v] & f & ~b; a ^= c; m ^= c
        pos, info, i, o, outs, rel_ = self._begin('S3')
        U, V, A, B, M, Fw = i['u'], i['v'], i['a'], i['b'], i['m'], i['f']

        def spec(al, val):
            u = [val(w) for w in U]
            v = [val(w) for w in V]
            a, b, m, f = val(A), val(B), val(M), val(Fw)
            c = al.and_(al.and_(_gt(al, u, v), f), al.not_(b))
            vals = dict(u=u, v=v, a=[al.xor(a, c)], b=[b], m=[al.xor(m, c)], f=[f])
            pre = [rel('R1', B, A, M), rel('R3', A, Fw, U[0]), rel('R2', B, Fw, U[0], V[0])]
            post = [rel('R1', o['b'], o['a'], o['m']), rel('Q1', o['a'], o['f'], o['u'], o['v']),
                    rel('R2', o['b'], o['f'], o['u'][0], o['v'][0])]
            return self._posts(al, outs, vals, rel_), [], pre, post, []
        self.seg('S3', k, pos, [Fw, A, B, M] + inter(U, V), spec)
        self._end(rel_)
        Bk, Mk, Fk = o['b'], o['m'], o['f']                      # not touched by S4
        # ---- S4 (Step4): if a: swap(u, v), swap(r, s)
        pos, info, i, o, outs, rel_ = self._begin('S4')
        U, V, R, S, A = i['u'], i['v'], i['r'], i['s'], i['a']

        def spec(al, val):
            u, v, r, s = ([val(w) for w in X] for X in (U, V, R, S))
            a = val(A)
            vals = dict(u=bv.mux(al, a, v, u), v=bv.mux(al, a, u, v), r=bv.mux(al, a, s, r), s=bv.mux(al, a, r, s), a=[a])
            pre = [rel('R1', Bk, A, Mk), rel('Q1', A, Fk, U, V), rel('R2', Bk, Fk, U[0], V[0])]
            post = [rel('R1', Bk, o['a'], Mk), rel('P1', Bk, Fk, o['u'], o['v']),
                    rel('P2', Bk, Fk, o['u'][0], o['v'][0]), rel('P3', Bk, Fk, o['v'][0])]
            return self._posts(al, outs, vals, rel_), [], pre, post, []
        self.seg('S4', k, pos, [A, Bk, Fk, Mk] + inter(U, V, R, S), spec)
        self._end(rel_)
        Ak = o['a']                                               # not touched by S5
        # ---- S5 (Step5): c = f & ~b; v = v - c*u (NOT, CAdd, NOT); s = s + c*r (CAdd)
        pos, info, i, o, outs, rel_ = self._begin('S5')
        U, V, R, S, B, Fw = i['u'], i['v'], i['r'], i['s'], i['b'], i['f']

        def spec(al, val):
            u, v, r, s = ([val(w) for w in X] for X in (U, V, R, S))
            b, f = val(B), val(Fw)
            c = al.and_(f, al.not_(b))
            v1, bor = bv.sub(al, v, [al.and_(c, x) for x in u])
            s1, car = bv.add(al, s, [al.and_(c, x) for x in r])
            vals = dict(u=u, v=v1, r=r, s=s1, b=[b], f=[f])
            pre = [rel('R1', B, Ak, Mk), rel('P1', B, Fw, U, V), rel('P2', B, Fw, U[0], V[0]), rel('P3', B, Fw, V[0])]
            post = [rel('R1', o['b'], Ak, Mk), rel('P4', o['f'], o['v'][0])]
            return (self._posts(al, outs, vals, rel_), [('M-carry', al.not_(car))], pre, post,
                    [('no borrow in v - c*u', al.not_(bor))])
        self.seg('S5', k, pos, [Fw, B, Ak, Mk] + inter(U, V, R, S), spec)
        self._end(rel_)
        # ---- S6 (Step6): b ^= m ^ a; if f: v = rotr(v); r = ModDbl(r); if a: swap(u,v), swap(r,s); a ^= [s even]
        pos, info, i, o, outs, rel_ = self._begin('S6')
        U, V, R, S, A, B, M, Fw = i['u'], i['v'], i['r'], i['s'], i['a'], i['b'], i['m'], i['f']

        def spec(al, val):
            u, v, r, s = ([val(w) for w in X] for X in (U, V, R, S))
            a, b, m, f = val(A), val(B), val(M), val(Fw)
            b1 = al.xor(b, al.xor(m, a))
            vr = [al.ite(f, v[(j + 1) % n], v[j]) for j in range(n)]
            t = [al.F] + r                                    # 2r on n+1 bits
            d, bor = bv.sub(al, t, bv.const(al, p, n + 1))     # bor = [2r < p]
            r2 = [al.ite(bor, t[j], d[j]) for j in range(n)]
            u2, v2 = bv.mux(al, a, vr, u), bv.mux(al, a, u, vr)
            r3, s3 = bv.mux(al, a, s, r2), bv.mux(al, a, r2, s)
            a1 = al.xor(a, al.not_(s3[0]))
            vals = dict(u=u2, v=v2, r=r3, s=s3, a=[al.F], b=[al.F], m=[m], f=[f])
            asserts = [('M-dbl', bv.lt_const(al, r, p)), ('M-a', al.or_(al.not_(a), bor)),
                       ('M-s', al.or_(a, s[0]))]
            proved = [('b cleared', al.not_(b1)), ('a cleared', al.not_(a1))]
            return (self._posts(al, outs, vals, rel_), asserts, [rel('R1', B, A, M), rel('P4', Fw, V[0])], [], proved)
        self.seg('S6', k, pos, [A, B, M, Fw] + inter(U, V, R, S), spec)
        self._end(rel_)
        for w in (A, B):
            self.require_zero(w, 'a/b freed at the end of the iteration')
            self.owned.discard(w)
        if not self.prove:                       # run mode: the state after the iteration (for the P_Q comparison)
            st = self.store
            self.trace.append(dict((r, [st.get(q, 0) for q in outs[r]]) for r in ('u', 'v', 'r', 's', 'f', 'm'))
                              | {'t': [st.get(T1, 0)]})

    def _fin(self):
        n, p = self.n, self.p
        pos, info = self.mark('fin')
        self.flow.feed('fin', info)              # u, v, r, s, f = the current ones; t0 = t of iteration 0
        U, V, R, S, Fw, T0 = info['u'], info['v'], info['r'], info['s'], info['f'], info['t0']
        X = list(self.C.outputs['x'])
        self.need(R == X, 'certificate: the output register is not r')
        self.flow.check_outputs(self.C.outputs)  # outputs[junk] = the t/m wires in the documented order
        # mathematical assumption (Z3 FINAL): v = 0 and f = 0 after the 2N iterations
        self.assume_zero(V + [Fw], 'M-end-vf')
        pm1 = p - 1
        bit = lambda al, k, i: al.T if (k >> i) & 1 else al.F

        def m_end(al, val):
            # M-end: (t0 -> u = p, s = 1, r = 0) and (not t0 -> u = 1, s = p)
            u, r, s = ([val(w) for w in W] for W in (U, R, S))
            t0 = val(T0)
            eq = lambda x, k: bv.allzero(al, [al.xor(x[i], bit(al, k, i)) for i in range(n)])
            return al.and_(al.or_(al.not_(t0), al.and_(al.and_(eq(u, p), eq(s, 1)), eq(r, 0))),
                           al.or_(t0, al.and_(eq(u, 1), eq(s, p))))
        # M-end as a PARAMETRISATION (Walker.seg subst): the entry values of u, s, r are replaced by their values under
        # M-end; the walker checks that this parametrises exactly the states satisfying M-end
        sb = []
        for i, w in enumerate(U):
            sb.append((w, lambda al, val, i=i: al.ite(val(T0), bit(al, p, i), bit(al, 1, i)), 'M-end', None, m_end))
        for i, w in enumerate(S):
            sb.append((w, lambda al, val, i=i: al.ite(val(T0), bit(al, 1, i), bit(al, p, i)), 'M-end', None, m_end))
        for i, w in enumerate(R):
            sb.append((w, lambda al, val, w=w: al.and_(al.not_(val(T0)), val(w)), 'M-end', None, m_end))

        def spec(al, val):
            u, r, s = ([val(w) for w in W] for W in (U, R, S))
            t0 = val(T0)
            rr, _ = bv.add(al, [al.not_(x) for x in r], bv.const(al, (p + 1) % (1 << n), n))
            r1 = [al.xor(rr[i], al.and_(t0, bit(al, p, i))) for i in range(n)]
            u1 = [al.xor(al.xor(u[i], bit(al, 1, i)), al.and_(t0, bit(al, pm1, i))) for i in range(n)]
            s1 = [al.xor(al.xor(s[i], bit(al, p, i)), al.and_(t0, bit(al, pm1, i))) for i in range(n)]
            posts = {}
            for W, x in ((R, r1), (U, u1), (S, s1)):
                posts.update(zip(W, x))
            return posts, [], [], [], []
        self.seg('fin', None, pos, [T0] + inter(U, R, S), spec, sb)
        self.out = R

    def _final_zero_check(self):
        if not self.prove:
            return
        keep = set(self.C.outputs['x']) | set(self.C.outputs['junk'])
        bad = [w for w in range(self.nq) if w not in keep and w not in self.Z]
        self.need(not bad, 'at the end %d wires outside x and junk are not proved zero (first %r)' % (len(bad), bad[:5]))
