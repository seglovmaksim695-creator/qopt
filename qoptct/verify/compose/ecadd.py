# SPDX-License-Identifier: Apache-2.0
"""compose.ecadd -- the contract chain of the exact elliptic-curve point addition P -> P + Q (Q classical) of
qoptct.crypto.ecpoint.ec_point_add (default construction, legacy=False).

The circuit is walked segment by segment (compose.walker); there are four kinds of segments:

  ROBDD segments (Walker.seg, as in compose.kal): modular additions / negations / constant additions, the shift-and-add
      multiplication by the constant 2^(2n) mod p, the flags and the substitution of the exceptional-case wrapper, the
      steps of the Kaliski inversion (compose.kal.Chain.kaliski_core on the wires of the denominator) and, for primes
      without a one-fold reduction, the signed-term reduction.  Proved for ALL entry values that satisfy the
      precondition, at the real width.
  ALGEBRA segments (seg_algebra): integer multiplication / squaring (out == x*y) and the fold reduction
      (r + p*q == V, r <= p-1).  The real ops of the segment are handed to qoptct.verify.algebra.prove_ops with a
      word-level spec built HERE; every wire of the segment that is not an input word must be known to be |0>.
      Afterwards the output words are arbitrary (unknown) values for the following segments; only the proved bound
      r < p is carried as a fact.
  MIRROR segments (_mirror): the uncompute half of every "compute / apply / uncompute" block must be, op by op, the
      syntactic inverse of the compute half (AND <-> erase), and the apply part must not change any wire that the
      compute half mentions in ANY role.  Then (mirror lemma, PROOF.md section 9) every mirrored gate meets the state
      of its original: all its preconditions hold and every wire of the compute half returns to its value before it.
  ZERO-LENGTH segments: named mathematical assumptions on the state (e.g. "the denominator is not 0").

Together the segments prove that the circuit REFINES the word-level program W_ec (the same chain executed on concrete
inputs in run mode: modular additions, products mod p, the Kaliski contract, flags), under the named assumptions:
    D-input   the input (x, y) is a point of the curve or O = (0, 0)             [the domain of the theorem]
    G1        the denominators x - xq and xq - x' are not 0 (the Kaliski inversion needs x in [1, p))
    G2        the second division returns lambda again (the lambda register returns to 0)
    G3, G4    after the generic addition: flagged inputs carry D + Q, the others a curve point not in {Q, -Q, 2Q, O}
    M1..M6    the lemmas of the Kaliski algorithm (compose/PROOF.md, proved in Z3 for all n, p, x in [1, p))
D-input is the hypothesis of the theorem; G1..G4 are consequences of elementary field / curve algebra (Theorem T,
PROOF.md section 10: written proof, polynomial identities checked by exact arithmetic in compose.ecmath, the whole of
W_ec executed exhaustively on small curves and sampled at n = 256 by wec_test()).
"""
import random
import sys
import threading
from array import array

from . import bv
from .bv import BddAlg
from .kal import Chain
from .lift import Lifted, K_X, K_CX, K_AND, K_ERASE
from .walker import Walker, ProofFailure, register_fact, inter as _inter

__all__ = ['ECChain', 'ECWinChain', 'ec_add', 'ec_neg', 'ec_mul', 'on_curve', 'lift_circuit', 'curve_points',
           'verify_ecadd', 'verify_ecadd_window', 'wec_test', 'wec_test_window', 'ecadd_hints', 'check_parameters']


# ------------------------------------------------------------------ classical reference (independent of the generator)
def ec_neg(P, p):
    return None if P is None else (P[0], (-P[1]) % p)


def ec_add(P, Q, p, a):
    """the textbook affine addition on y^2 = x^3 + a x + b over F_p (None = O). This function IS the specification."""
    if P is None:
        return Q
    if Q is None:
        return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        l = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        l = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (l * l - x1 - x2) % p
    return x3, (l * (x1 - x3) - y1) % p


def ec_mul(k, P, p, a):
    R = None
    while k:
        if k & 1:
            R = ec_add(R, P, p, a)
        P = ec_add(P, P, p, a)
        k >>= 1
    return R


def on_curve(P, p, a, b):
    return P is None or (0 <= P[0] < p and 0 <= P[1] < p and (P[1] * P[1] - P[0] ** 3 - a * P[0] - b) % p == 0)


def curve_points(p, a, b):
    """all affine points, by brute force (small p)."""
    sq = {}
    for y in range(p):
        sq.setdefault(y * y % p, []).append(y)
    return [(x, y) for x in range(p) for y in sq.get((x ** 3 + a * x + b) % p, [])]


def _naf(c):
    out = []
    e = 0
    while c:
        if c & 1:
            d = 2 - (c & 3)
            out.append((d, e))
            c -= d
        c >>= 1
        e += 1
    return out


# ------------------------------------------------------------------ facts
def _eqc(al, bits, c):
    """[value(bits) == c] for a constant c."""
    if c >> len(bits):
        return al.F
    r = al.T
    for i, b in enumerate(bits):
        r = al.and_(r, b if (c >> i) & 1 else al.not_(b))
    return r


def _ecdom_bit(al, x, y, p, consts):
    """the consequences of "(x, y) is O = (0, 0) or a point of the curve" that the wrapper uses:
    x, y < p;  y = 0 -> x = 0 (no curve point has y = 0);  x = xc and y != 0 -> y = yc or y = p - yc, for the two
    constants (xc, yc) = Q, 2Q (the points of the curve with abscissa xc are (xc, +-yc))."""
    y0 = bv.allzero(al, y)
    r = al.and_(bv.lt_const(al, x, p), bv.lt_const(al, y, p))
    r = al.and_(r, al.or_(al.not_(y0), bv.allzero(al, x)))
    for xc, yc in consts:
        hyp = al.and_(_eqc(al, x, xc), al.not_(y0))
        r = al.and_(r, al.or_(al.not_(hyp), al.or_(_eqc(al, y, yc), _eqc(al, y, p - yc))))
    return r


def _ecdom(al, val, key):
    return _ecdom_bit(al, [val(w) for w in key[1]], [val(w) for w in key[2]], key[3], key[4])


def _onehot0(al, val, key):
    f = [val(w) for w in key[1]]
    r = al.T
    for i in range(len(f)):
        for j in range(i + 1, len(f)):
            r = al.and_(r, al.not_(al.and_(f[i], f[j])))
    return r


_MUXC = [None]


def _muxbits(al, kb, consts, nbits):
    """(bits, eqk): bit i of consts[k] selected by the index bits kb, and the indicators [k == j]. One-entry cache
    that holds the algebra object itself (never its id: ids are reused)."""
    c = _MUXC[0]
    if c is not None and c[0] is al and c[1] == kb and c[2] is consts and c[3] == nbits:
        return c[4], c[5]
    eqk = [_eqc(al, kb, j) for j in range(len(consts))]
    bits = []
    for i in range(nbits):
        r = al.F
        for j, cst in enumerate(consts):
            if (cst >> i) & 1:
                r = al.or_(r, eqk[j])
        bits.append(r)
    _MUXC[0] = (al, kb, consts, nbits, bits, eqk)
    return bits, eqk


def _tab(al, val, key):
    kb = tuple(val(w) for w in key[1])
    bits, _ = _muxbits(al, kb, key[3], len(key[2]))
    r = al.T
    for w, m in zip(key[2], bits):
        r = al.and_(r, al.not_(al.xor(val(w), m)))
    return r


def _krestrict(al, kb):
    """ROBDD algebra only: restr(j, bits) = the functions `bits` with the index k fixed to j (Shannon cofactors).
    Used as  [k = j] and P(x)  ==  [k = j] and P(x | k = j): the per-k predicates are then evaluated on functions
    without k (16 times smaller). Requires the index bits to be the top variables of the order (checked).
    Bit-sliced algebra: None (concrete values need no restriction)."""
    if not isinstance(al, BddAlg):
        return None
    B = al.B
    vs = []
    for i, e in enumerate(kb):
        if e == al.T or e == al.F:
            continue
        v = B.var[e >> 1]
        if B.cof(e, v) != (al.T, al.F):
            raise ProofFailure('internal: an index bit is not a variable')
        vs.append((v, i))
    vs.sort(reverse=True)
    if [v for v, _ in vs] != list(range(B.nvars - 1, B.nvars - 1 - len(vs), -1)):
        raise ProofFailure('internal: the index bits are not the top ROBDD variables')
    cof = B.cof

    def restr(j, bits):
        out = []
        for f in bits:
            for v, i in vs:
                f1, f0 = cof(f, v)
                f = f1 if (j >> i) & 1 else f0
            out.append(f)
        return out
    return restr


def _fcode(al, val, key, iff):
    """flags against k-dependent points: for every flag j,  f_j -> (x, y) = codes[j][k]   (iff: <->)."""
    f = [val(w) for w in key[1]]
    x = [val(w) for w in key[2]]
    y = [val(w) for w in key[3]]
    kb = [val(w) for w in key[4]]
    eqk = [_eqc(al, kb, j) for j in range(1 << len(kb))]
    restr = _krestrict(al, kb)
    xy = {}
    r = al.T
    for j, row in enumerate(key[5]):
        hit = al.F
        for k, (cx_, cy_) in enumerate(row):
            if eqk[k] == al.F:
                continue
            if k not in xy:
                xy[k] = (restr(k, x), restr(k, y)) if restr else (x, y)
            xr, yr = xy[k]
            hit = al.or_(hit, al.and_(eqk[k], al.and_(_eqc(al, xr, cx_), _eqc(al, yr, cy_))))
        r = al.and_(r, al.not_(al.xor(f[j], hit)) if iff else al.or_(al.not_(f[j]), hit))
    return r


def _ecdomw(al, val, key):
    x = [val(w) for w in key[1]]
    y = [val(w) for w in key[2]]
    kb = [val(w) for w in key[3]]
    restr = _krestrict(al, kb)
    r = al.T
    for j, consts in enumerate(key[5]):
        ej = _eqc(al, kb, j)
        if ej == al.F:
            continue
        xr, yr = (restr(j, x), restr(j, y)) if restr else (x, y)
        r = al.and_(r, al.or_(al.not_(ej), _ecdom_bit(al, xr, yr, key[4], consts)))
    return r


# ('ecdom', X, Y, p, ((xq, yq), (x2q, y2q)));   ('onehot0', flags): at most one flag is set
# window: ('tab', K, R, consts): the wires R hold consts[k];  ('ecdomw', X, Y, K, p, per-k consts): ecdom for T[k];
#         ('xor3', (g, f1, f2)): g = f1 xor f2
register_fact('ecdom', lambda key: key[1] + key[2], _ecdom)
register_fact('onehot0', lambda key: key[1], _onehot0)
register_fact('tab', lambda key: key[1] + key[2], _tab)
register_fact('ecdomw', lambda key: key[1] + key[2] + key[3], _ecdomw)
# ('fimp' / 'fiff', flags, X, Y, K, codes): flag j implies / is equivalent to (x, y) = codes[j][k]
register_fact('fimp', lambda key: key[1] + key[2] + key[3] + key[4], lambda al, val, key: _fcode(al, val, key, False))
register_fact('fiff', lambda key: key[1] + key[2] + key[3] + key[4], lambda al, val, key: _fcode(al, val, key, True))
register_fact('xor3', lambda key: key[1],
              lambda al, val, key: al.not_(al.xor(val(key[1][0]), al.xor(val(key[1][1]), val(key[1][2])))))


def lift_circuit(C):
    """the lifted form (compose.lift.Lifted) of an engine circuit of qoptct.crypto._ecp (ECPointAddCircuit.circ)
    WITHOUT going through the OpenQASM text: for tests and mutants. The proofs of record lift the text."""
    a = C.ops.a
    out = array('i')
    ext = out.extend
    nm = 0
    for i in range(0, len(a), 5):
        k = a[i]
        if k == 1:
            ext((K_CX, a[i + 1], a[i + 2], 0, -1))
        elif k == 0:
            ext((K_X, a[i + 1], 0, 0, -1))
        elif k == 3:
            ext((K_AND, a[i + 1], a[i + 2], a[i + 3], -1))
        elif k == 4:
            ext((K_ERASE, a[i + 2], a[i + 3], a[i + 1], a[i + 4]))
            nm += 1
        else:
            raise ValueError('ccx is not a gadget of the lifted language')
    return Lifted(C.nw, out, nm)


def _run_deep(fn):
    """run fn() in a thread with a deep stack (the algebra prover recurses deeply)."""
    out = {}

    def run():
        try:
            out['r'] = fn()
        except RecursionError:
            out['r'] = dict(status='UNKNOWN', detail='recursion limit')
        except MemoryError:
            out['r'] = dict(status='UNKNOWN', detail='out of memory')
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
    return out.get('r') or dict(status='ERROR', detail='prover thread died')


def _fold_spec(p, n, P, extra):
    """the word-level spec of the library fold reduction on the 2n wires P (in place: r on P[:n]) with the extra
    wires thi1, g, M named by the certificate. The claims have the form  r + p * (integer expression) == value(P)
    and r <= p - 1: for ANY integer expression they say r = value(P) mod p."""
    from ...crypto import fold as _fold
    try:
        _, c, method, step2 = _fold._resolve(p, n, 'auto', 'auto')
    except ValueError as ex:
        raise ProofFailure('the fold reduction does not apply to p: %s' % ex)
    if method == 'solinas':
        e1, wM = 5, 0
    else:
        e1 = _fold._fold_e(n, c) + 1
        wM = _fold._fold_wM(n, c) if step2 == 'M' else 0
    if len(extra) != e1 + 1 + wM or len(P) != 2 * n:
        raise ProofFailure('certificate: bad layout of the fold reduction')
    thi1, g, M = list(extra[:e1]), [extra[e1]], list(extra[e1 + 1:])
    outs = dict(r=list(P[:n]), thi1=thi1, g=g)
    gar = thi1 + g
    claims = []
    if wM:
        outs['M'] = M
        gar = gar + M
        claims.append('M == %d*thi1' % c)
    if method == 'solinas':
        ins = {'A%d' % i: list(P[32 * i:32 * i + 32]) for i in range(16)}
        qs = []
        for i in range(16):
            d = (1 << (32 * i)) - _fold.p256_formula(1 << (32 * i))
            if d % p:
                raise ProofFailure('internal: NIST formula coefficient')
            qs.append('%d*A%d' % (d // p, i))
        claims.append('r + %d*(thi1 - 1 + g - 5 + %s) == %s' % (
            p, ' + '.join(qs), ' + '.join('%d*A%d' % (1 << (32 * i), i) for i in range(16))))
    else:
        ins = dict(L=list(P[:n]), H=list(P[n:]))
        claims.append('r + %d*(H + thi1 - 1 + g) == L + %d*H' % (p, 1 << n))
    claims.append('r <= %d' % (p - 1))
    return dict(inputs=ins, outputs=outs, claims=claims, garbage=gar, keep=list(P[n:]))


class ECChain(Chain):
    """the walk of the point addition; see the module docstring. prove mode: lifted given; run mode: points given
    (the word-level program W_ec executed bit-sliced on these inputs, every assumption and fact checked)."""

    def __init__(self, p, a, b, Q, D, hints, lifted=None, points=None, progress=None, budget=900.0):
        self.p, self.a, self.b = p, a, b
        self.n = n = p.bit_length()
        self.red = 'none'
        self.Q, self.D = (int(Q[0]), int(Q[1])), (int(D[0]), int(D[1]))
        self.budget = budget
        self._written = None
        self.alg = []
        self.nmirror = 0
        X, Y = list(range(n)), list(range(n, 2 * n))
        init = None
        if lifted is None:
            points = list(points)
            init = {}
            for j, P in enumerate(points):
                x, y = (0, 0) if P is None else P
                for i in range(n):
                    if (x >> i) & 1:
                        init[i] = init.get(i, 0) | (1 << j)
                    if (y >> i) & 1:
                        init[n + i] = init.get(n + i, 0) | (1 << j)
        Walker.__init__(self, hints, X + Y, lifted=lifted, init=init, labels=points, progress=progress)
        self.X, self.Y = X, Y
        Q = self.Q
        Q2 = ec_add(Q, Q, p, a)
        if Q2 is None:
            raise ProofFailure('2Q = O')
        self.Q2, self.mQ, self.m2Q = Q2, ec_neg(Q, p), ec_neg(Q2, p)
        self.DQ = ec_add(self.D, Q, p, a)
        self.E = [None, Q, self.mQ, self.m2Q]                 # the flags fO, fQ, fmQ, fm2Q
        self.img = [Q, Q2, None, self.mQ]                     # their images E_j + Q
        self.consts = (Q, Q2)

    # ------------------------------------------------------------------ small tools
    def _here(self):
        return self.pos if self.prove else None

    def seg(self, label, k, pos, order, spec, subst=()):
        """Walker.seg, plus: (i) the wires with a post are recorded while an apply part is walked (frame condition
        of the mirror lemma); (ii) a fact about wires WITHOUT a post survives the segment: the segment proof shows
        that every touched wire without a post ends with its entry value (under the precondition, which holds on
        the real runs), so a predicate over such wires that was true at the entry is true at the exit."""
        from .walker import fact_wires
        wr = self._written
        posted = [None]

        def spec2(al, val, spec=spec):
            r = spec(al, val)
            posted[0] = set(r[0])
            if wr is not None:
                wr.update(r[0])
            return r
        before = set(self.facts) if self.prove else None
        Walker.seg(self, label, k, pos, order, spec2, subst)
        if before is not None and posted[0] is not None:
            for key in before - self.facts:
                if posted[0].isdisjoint(fact_wires(key)):
                    self.facts.add(key)

    def assume_fact(self, key, name):
        """a named mathematical assumption given as a fact."""
        from .walker import fact_bit
        if self.prove:
            self.assumed[name] = self.assumed.get(name, 0) + 1
        else:
            bit = fact_bit(self.al, self._val_run, key)
            if bit != self.al.T:
                self._violation('assume', None, name, bit ^ self.al.T)
        self.facts.add(key)

    def disown(self, wires):
        """wires released by the circuit: they must be known to be |0>."""
        for w in wires:
            self.require_zero(w, 'released wire')
            self.owned.discard(w)

    def _get(self, wires):
        out = [0] * len(self.xs)
        st = self.store
        for i, w in enumerate(wires):
            m = st.get(w, 0)
            while m:
                j = (m & -m).bit_length() - 1
                out[j] |= 1 << i
                m &= m - 1
        return out

    def _put(self, wires, vals):
        for i, w in enumerate(wires):
            m = 0
            for j, v in enumerate(vals):
                if (v >> i) & 1:
                    m |= 1 << j
            self.store[w] = m

    # ------------------------------------------------------------------ algebra segments
    def seg_algebra(self, label, pos, spec, runf, newfacts=()):
        """one segment proved by qoptct.verify.algebra against the word-level `spec` (inputs: free words; outputs /
        garbage: wires that are |0> at the entry; claims). runf(get, put): the same contract on concrete values."""
        from .walker import fact_wires, fact_bit
        self.nseg += 1
        ins = [w for ws in spec['inputs'].values() for w in ws]
        outw = set(w for ws in spec['outputs'].values() for w in ws) | set(spec.get('garbage', []))
        inset = set(ins)
        self.need(len(inset) == len(ins), 'internal: input words of %r overlap' % label)
        if not self.prove:
            runf(self._get, self._put)
            self.facts = {k for k in self.facts if outw.isdisjoint(fact_wires(k))}
            for key in newfacts:
                bit = fact_bit(self.al, self._val_run, key)
                if bit != self.al.T:
                    self._violation(label, None, 'fact ' + key[0], bit ^ self.al.T)
                self.facts.add(key)
            return
        lo, hi = self.pos, pos
        sup, _ = self._support(lo, hi)
        Z = self.Z
        for w in (sup | outw) - inset:
            if w not in Z:
                raise ProofFailure('segment %s (ops %d..%d): wire %d is neither an input word nor known to be |0>'
                                   % (label, lo, hi, w), segment=label, lo=lo, hi=hi)
        ops = self.ops
        alg = []
        app = alg.append
        nand = 0
        for i in range(5 * lo, 5 * hi, 5):
            kd = ops[i]
            if kd == K_CX:
                app(('cx', ops[i + 1], ops[i + 2]))
            elif kd == K_X:
                app(('x', ops[i + 1]))
            elif kd == K_AND:
                app(('and', ops[i + 1], ops[i + 2], ops[i + 3]))
                nand += 1
            else:
                l = ops[i + 4]
                app(('erase', ops[i + 3], (('cz', ops[i + 1], ops[i + 2]),) + ((('z', l),) if l >= 0 else ())))
        from .. import algebra
        nq, budget = self.nq, self.budget
        r = _run_deep(lambda: algebra.prove_ops(nq, alg, spec, budget_sec=budget))
        del alg
        self.alg.append(dict(label=label, lo=lo, hi=hi, AND=nand, status=r.get('status'), time=r.get('time')))
        if self.progress:
            self.progress('algebra %s ops %d..%d: %s (%s s)' % (label, lo, hi, r.get('status'), r.get('time')))
        if r.get('status') != 'PROVED':
            raise ProofFailure('segment %s (ops %d..%d): the algebra sub-goal is %s: %s'
                               % (label, lo, hi, r.get('status'), str(r.get('detail'))[:300]),
                               segment=label, lo=lo, hi=hi, algebra=r.get('status'))
        Z -= outw
        self.facts = {k for k in self.facts if outw.isdisjoint(fact_wires(k))}
        self.facts.update(newfacts)
        self.pos = pos

    # ------------------------------------------------------------------ compute / apply / uncompute
    def _cau(self, compute, apply):
        if self.prove:
            snap = dict(pos=self.pos, Z=set(self.Z), owned=set(self.owned), facts=set(self.facts))
        else:
            snap = dict(store=dict(self.store), owned=set(self.owned), facts=set(self.facts))
        res = compute()
        mid = self.pos if self.prove else dict(self.store)
        self._written = set()
        try:
            apply(res)
        finally:
            written, self._written = self._written, None
        pos, _ = self.mark('mirror')
        self._mirror(snap, mid, pos, written)

    def _mirror(self, snap, mid, pos, written):
        from .walker import fact_wires
        self.nseg += 1
        self.nmirror += 1
        if not self.prove:
            st, st0 = self.store, snap['store']
            changed = [w for w in mid if mid[w] != st0.get(w, 0)]
            if written & set(changed):
                self.viol.append(dict(segment='mirror', what='the apply part writes a wire of the compute part',
                                      x=None))
            for w in changed:
                st[w] = st0.get(w, 0)
            keep = {k for k in snap['facts'] if written.isdisjoint(fact_wires(k))}
            ch = set(changed)
            self.facts = keep | {k for k in self.facts if ch.isdisjoint(fact_wires(k))}
            self.owned = set(snap['owned'])
            return
        s, e, lo = snap['pos'], mid, self.pos
        ops = self.ops
        j, end = 5 * lo, 5 * pos

        def fail(i):
            raise ProofFailure('mirror segment (ops %d..%d): op %d is not the inverse of op %d of the compute part '
                               '%d..%d' % (lo, pos, j // 5, i // 5, s, e), segment='mirror', lo=lo, hi=pos, op=j // 5)
        for i in range(5 * (e - 1), 5 * s - 1, -5):
            if j >= end:
                fail(i)
            kd = ops[i]
            if kd == K_CX or kd == K_X:
                if ops[j] != kd or ops[j + 1] != ops[i + 1] or ops[j + 2] != ops[i + 2]:
                    fail(i)
                j += 5
            elif kd == K_AND:
                if (ops[j] != K_ERASE or ops[j + 1] != ops[i + 1] or ops[j + 2] != ops[i + 2]
                        or ops[j + 3] != ops[i + 3] or ops[j + 4] != -1):
                    fail(i)
                j += 5
            else:
                if (ops[j] != K_AND or ops[j + 1] != ops[i + 1] or ops[j + 2] != ops[i + 2]
                        or ops[j + 3] != ops[i + 3]):
                    fail(i)
                j += 5
                l = ops[i + 4]
                if l >= 0:
                    if j >= end or ops[j] != K_CX or ops[j + 1] != l or ops[j + 2] != ops[i + 3]:
                        fail(i)
                    j += 5
        if j != end:
            raise ProofFailure('mirror segment (ops %d..%d): %d extra ops after the inverse of the compute part'
                               % (lo, pos, (end - j) // 5), segment='mirror', lo=lo, hi=pos)
        supF, _ = self._support(s, e)
        bad = written & supF
        if bad:
            raise ProofFailure('mirror segment (ops %d..%d): the apply part changes wire %d, which the compute part '
                               'mentions' % (lo, pos, min(bad)), segment='mirror', lo=lo, hi=pos)
        self.Z = (self.Z - supF) | (snap['Z'] & supF)
        keep = {k for k in snap['facts'] if written.isdisjoint(fact_wires(k))}
        self.facts = keep | {k for k in self.facts if supF.isdisjoint(fact_wires(k))}
        self.owned = set(snap['owned'])
        self.pos = pos
        if self.progress:
            self.progress('mirror ops %d..%d = inverse of %d..%d' % (lo, pos, s, e))

    # ------------------------------------------------------------------ ROBDD contracts of the modular layer
    def _modsum(self, al, s, c):
        """(s + 2^n c) mod p for a value < 2p given as n bits s and the carry c."""
        n, p = self.n, self.p
        ge = al.or_(c, al.not_(bv.lt_const(al, s, p)))
        d, _ = bv.sub(al, s, bv.const(al, p, n))
        return [al.ite(ge, d[i], s[i]) for i in range(n)]

    def _addc(self, R, k):
        """R <- R + k mod p (mark 'addc')."""
        n, p = self.n, self.p
        k %= p
        pos, _ = self.mark('addc')
        lt = ('lt', tuple(R), p)

        def spec(al, val):
            if k == 0:
                return {}, [], [lt], [lt], []
            r = [val(w) for w in R]
            s, c = bv.add(al, r, bv.const(al, k, n))
            res = self._modsum(al, s, c)
            return {R[i]: res[i] for i in range(n)}, [], [lt], [lt], []
        self.seg('addc', None, pos, list(R), spec)

    def _msub(self, A, R):
        """R <- R - A mod p (mark 'msub'; window variant: A is a looked-up register)."""
        n, p = self.n, self.p
        pos, _ = self.mark('msub')
        ltA, ltR = ('lt', tuple(A), p), ('lt', tuple(R), p)

        def spec(al, val):
            a = [val(w) for w in A]
            r = [val(w) for w in R]
            d, bor = bv.sub(al, r, a)
            dp, _ = bv.add(al, d, bv.const(al, p, n))
            return {R[i]: al.ite(bor, dp[i], d[i]) for i in range(n)}, [], [ltA, ltR], [ltA, ltR], []
        self.seg('msub', None, pos, _inter(A, R), spec)

    def _madd(self, P, R, label='apply'):
        """R <- R + P mod p (the apply part ends with the mark 'apply'; window variant: mark 'madd')."""
        n, p = self.n, self.p
        pos, _ = self.mark(label)
        ltP, ltR = ('lt', tuple(P), p), ('lt', tuple(R), p)

        def spec(al, val):
            a = [val(w) for w in P]
            r = [val(w) for w in R]
            s, c = bv.add(al, a, r)
            res = self._modsum(al, s, c)
            return {R[i]: res[i] for i in range(n)}, [], [ltP, ltR], [ltP, ltR], []
        self.seg('madd', None, pos, _inter(P, R), spec)

    def _neg(self, R):
        """R <- -R mod p (0 -> 0) (mark 'neg')."""
        n, p = self.n, self.p
        pos, _ = self.mark('neg')
        lt = ('lt', tuple(R), p)

        def spec(al, val):
            r = [val(w) for w in R]
            nz = al.not_(bv.allzero(al, r))
            d, _ = bv.sub(al, bv.const(al, p, n), r)
            return {R[i]: al.and_(nz, d[i]) for i in range(n)}, [], [lt], [lt], []
        self.seg('neg', None, pos, list(R), spec)

    def _xor(self, P, Lr, first):
        """Lr ^= P (mark 'apply'). first: Lr is |0> and receives P. second: G2 (P = Lr), Lr returns to |0>."""
        n, p = self.n, self.p
        pos, _ = self.mark('apply')
        ltP, ltL = ('lt', tuple(P), p), ('lt', tuple(Lr), p)
        if first:
            for w in Lr:
                self.require_zero(w, 'the lambda register before the first division')

            def spec(al, val):
                return {Lr[i]: al.xor(val(Lr[i]), val(P[i])) for i in range(n)}, [], [ltP], [ltP, ltL], []
        else:
            def spec(al, val):
                eq = al.T
                for i in range(n):
                    eq = al.and_(eq, al.not_(al.xor(val(Lr[i]), val(P[i]))))
                return {Lr[i]: al.F for i in range(n)}, [('G2-lambda-recomputed', eq)], [], [], []
        self.seg('xor', None, pos, _inter(P, Lr), spec)

    def _nonzero(self, R, name):
        def spec(al, val):
            return {}, [(name, al.not_(bv.allzero(al, [val(w) for w in R])))], [], [], []
        self.seg('assume', None, self._here(), list(R), spec)

    # ------------------------------------------------------------------ the arithmetic blocks
    def _reduce(self, V, Bd):
        """V (wires, value < Bd) -> the n wires of V mod p; fact ('lt', out, p)."""
        n, p = self.n, self.p
        need = self.need
        V = list(V)
        if self.peek() == 'fold':
            pos, info = self.mark('fold')
            pad = self.fresh_list(info.get('pad'))
            extra = self.fresh_list(info.get('extra'))
            P = V + pad
            need(len(P) == 2 * n and Bd <= (1 << (2 * n)), 'certificate: bad input of the fold reduction')
            spec = _fold_spec(p, n, P, extra)
            out = P[:n]

            def runf(get, put):
                put(out, [v % p for v in get(P)])
            self.seg_algebra('fold', pos, spec, runf, newfacts=[('lt', tuple(out), p)])
            return out
        self.mark('sp_begin')
        if Bd == (1 << len(V)):
            self.establish(('lt', tuple(V), Bd))
        out = self._reduce_special(V, Bd)
        cur, Bl = self._sp_last
        need(Bl <= p and list(cur[:n]) == list(out), 'internal: special reduction')
        key0, key1 = ('lt', tuple(cur), Bl), ('lt', tuple(out), p)

        def spec(al, val):
            return {}, [], [key0], [key0, key1], []
        self.seg('sp_out', None, self._here(), list(cur), spec)
        self.mark('sp_end')
        return out

    def _constmul(self, Xr):
        """the n wires of 2^(2n) * value(Xr) mod p. Shift-and-add rounds: literal additions mod 2^(L-e) into the
        suffix A[e:], i.e. A <- A + s 2^e x (mod 2^L); the exact integer K + sum s 2^e x lies in [0, Bn) with
        Bn <= 2^L (numeric side conditions checked here), so A holds it exactly (fact lt, admitted)."""
        n, p = self.n, self.p
        need = self.need
        c2 = pow(2, 2 * n, p)
        pos, info = self.mark('cm_const')
        K = info.get('K')
        A = self.fresh_list(info.get('A'))
        L = len(A)
        need(isinstance(K, int) and not isinstance(K, bool) and L >= 1 and 0 <= K < (1 << L) and K % p == 0,
             'certificate: bad constant of the constant multiplication')

        def spec(al, val):
            return {A[i]: (al.T if (K >> i) & 1 else al.F) for i in range(L)}, [], [], [], []
        self.seg('cm_const', None, pos, A, spec)
        ltX = ('lt', tuple(Xr), p)
        terms = []
        while self.peek() == 'cm_term':
            pos, info = self.mark('cm_term')
            s_, e_ = info.get('s'), info.get('e')
            need(s_ in (1, -1) and isinstance(e_, int) and not isinstance(e_, bool) and 0 <= e_ < L,
                 'certificate: bad term of the constant multiplication')
            sub_, h = A[e_:], Xr[:L - e_]

            def spec(al, val, sub_=sub_, h=h, s_=s_):
                a = [val(x) for x in sub_]
                hv = [val(x) for x in h] + [al.F] * (len(sub_) - len(h))
                res = bv.add(al, a, hv)[0] if s_ > 0 else bv.sub(al, a, hv)[0]
                return {sub_[i]: res[i] for i in range(len(sub_))}, [], [ltX], [ltX], []
            self.seg('cm_term', e_, pos, _inter(sub_, h), spec)
            terms.append((int(s_), e_))
        xm = (1 << n) - 1                                    # x < p <= xm: the (looser) bound the generator uses
        neg = sum(xm << e_ for s_, e_ in terms if s_ < 0)
        pos_ = sum(xm << e_ for s_, e_ in terms if s_ > 0)
        Bn = K + pos_ + 1
        need((sum(s_ << e_ for s_, e_ in terms) - c2) % p == 0, 'constant multiplication: the terms are not 2^(2n) mod p')
        need(K >= neg and Bn <= (1 << L), 'constant multiplication: side conditions fail')
        need(ltX in self.facts, 'internal: bound of the multiplicand lost')
        self.admit(('lt', tuple(A), Bn), 'shift-and-add lemma')
        return self._reduce(A, Bn)

    def _pmul(self, a, b):
        """the n wires of value(a) * value(b) mod p (a is b: square)."""
        n, p = self.n, self.p
        sq = a is b
        pos, info = self.mark('sq' if sq else 'mul')
        out = self.fresh_list(info.get('out'))
        gar = self.fresh_list(info.get('gar'))
        self.need(len(out) == 2 * n, 'certificate: bad product register')
        if sq:
            spec = dict(inputs=dict(x=list(a)), outputs=dict(out=out), claims=['out == x*x'], garbage=gar)

            def runf(get, put):
                put(out, [v * v for v in get(a)])
        else:
            spec = dict(inputs=dict(x=list(a), y=list(b)), outputs=dict(out=out), claims=['out == x*y'], garbage=gar)

            def runf(get, put):
                put(out, [u * v for u, v in zip(get(a), get(b))])
        self.seg_algebra('sq' if sq else 'mul', pos, spec, runf)
        Bd = (p - 1) ** 2 + 1
        if self.peek() != 'fold':
            self.need(('lt', tuple(a), p) in self.facts and ('lt', tuple(b), p) in self.facts,
                      'internal: bounds of the factors lost')
            self.admit(('lt', tuple(out), Bd), 'product bound: out == x*y with x, y <= p - 1')
        return self._reduce(out, Bd)

    def _div(self, num, den):
        """the n wires of -num / den mod p:  x' = 2^(2n) den;  r_ext = Kaliski(x') = -x'^{-1} 2^(2n) = -den^{-1};
        K = r_ext mod p;  num * K mod p."""
        n = self.n
        xs = self._constmul(den)
        self._nonzero(xs, 'G1-denominator-nonzero')
        if self.progress:
            self.progress('kaliski inversion')
        S = self.kaliski_core(xs)
        K = self._reduce(S, 1 << (2 * n))
        return self._pmul(num, K)

    # ------------------------------------------------------------------ the chain: W_ec
    def run(self):
        n, p = self.n, self.p
        need = self.need
        X, Y = self.X, self.Y
        xq, yq = self.Q
        E, D, img, consts = self.E, self.D, self.img, self.consts
        code = lambda P: (0, 0) if P is None else P
        bit = lambda al, v, i: al.T if (v >> i) & 1 else al.F
        dom = ('ecdom', tuple(X), tuple(Y), p, consts)
        self.assume_fact(dom, 'D-input')
        ltX, ltY = ('lt', tuple(X), p), ('lt', tuple(Y), p)
        # ---- wrapper, part 1: flags f_j = [P == E_j], substitution of the safe point D
        pos, info = self.mark('w_pre')
        F = [self.fresh(info.get(k)) for k in ('fO', 'fQ', 'fmQ', 'fm2Q')]
        oh = ('onehot0', tuple(F))

        def spec(al, val):
            x = [val(w) for w in X]
            y = [val(w) for w in Y]
            ind = [al.and_(_eqc(al, x, code(P)[0]), _eqc(al, y, code(P)[1])) for P in E]
            anyf = al.F
            for f in ind:
                anyf = al.or_(anyf, f)
            posts = {F[j]: ind[j] for j in range(4)}
            for i in range(n):
                posts[X[i]] = al.ite(anyf, bit(al, D[0], i), x[i])
                posts[Y[i]] = al.ite(anyf, bit(al, D[1], i), y[i])
            return posts, [], [dom], [ltX, ltY, oh], []
        self.seg('w_pre', None, pos, _inter(X, Y) + F, spec)
        # ---- the generic addition (chord formula) on P' = D or P
        self._addc(X, -xq)                                   # X = dx = x - xq
        self._addc(Y, -yq)                                   # Y = dy = y - yq
        pos, info = self.mark('lam')
        Lr = self.fresh_list(info.get('Lr'))
        need(len(Lr) == n, 'certificate: bad lambda register')
        self._cau(lambda: self._div(Y, X), lambda P: self._xor(P, Lr, True))        # Lr = -dy/dx = -lambda
        self._cau(lambda: self._pmul(Lr, X), lambda P: self._madd(P, Y))            # Y = dy - lambda dx (= 0)
        self._cau(lambda: self._pmul(Lr, Lr), lambda P: (self._neg(X), self._madd(P, X)))   # X = lambda^2 - dx
        self._addc(X, -3 * xq)                               # X = x' - xq
        self._cau(lambda: self._pmul(Lr, X), lambda P: self._madd(P, Y))            # Y = lambda (xq - x') = y' + yq
        self._neg(X)                                         # X = xq - x'
        self._cau(lambda: self._div(Y, X), lambda P: self._xor(P, Lr, False))       # Lr ^= -(y' + yq)/(xq - x') -> 0
        self._neg(X)
        self._addc(X, xq)                                    # X = x'
        self._addc(Y, -yq)                                   # Y = y'
        self.disown(Lr)
        # ---- wrapper, part 2: images of the exceptional inputs, erasure of the flags from the result
        pos, _ = self.mark('w_post')
        DQ, Q, mQ, Q2 = self.DQ, self.Q, self.mQ, self.Q2

        def spec(al, val):
            x = [val(w) for w in X]
            y = [val(w) for w in Y]
            f = [val(w) for w in F]
            eqp = lambda P: al.and_(_eqc(al, x, code(P)[0]), _eqc(al, y, code(P)[1]))
            anyf = al.F
            for t in f:
                anyf = al.or_(anyf, t)
            a1 = al.or_(al.not_(anyf), eqp(DQ))
            bad = al.or_(al.or_(eqp(Q), eqp(mQ)), al.or_(eqp(Q2), eqp(None)))
            a2 = al.or_(anyf, al.and_(_ecdom_bit(al, x, y, p, consts), al.not_(bad)))
            posts = {w: al.F for w in F}
            for i in range(n):
                vx, vy = x[i], y[i]
                for j in range(4):
                    cx_, cy_ = code(img[j])
                    vx = al.ite(f[j], bit(al, cx_, i), vx)
                    vy = al.ite(f[j], bit(al, cy_, i), vy)
                posts[X[i]] = vx
                posts[Y[i]] = vy
            return posts, [('G3-flagged-result', a1), ('G4-generic-result', a2)], [oh], [ltX, ltY], []
        self.seg('w_post', None, pos, F + _inter(X, Y), spec)
        self.disown(F)
        self.finish(())
        if self.prove:
            dirty = set(range(self.nq)) - self.Z - set(X) - set(Y)
            need(not dirty, 'ancilla wires are not proved to end at |0> (%d wires, e.g. %s)'
                 % (len(dirty), sorted(dirty)[:4]))
        return self


# ====================================================================== window variant: (P, k) -> (P + T[k], k)
WFIELDS = ('xq', 'yq', 'ymq', 'x3q', 'x2q', 'y2q', 'ym2q', 'xdq', 'ydq')


class ECWinChain(ECChain):
    """the walk of the windowed addition (qoptct.crypto.ecpoint.ec_point_add_window, default construction): wires
    0..n-1 x, n..2n-1 y, 2n..2n+w-1 the index k; Q = table[k] is looked up into nine registers (xq, yq, -yq, 3xq,
    x(2Q), y(2Q), -y(2Q), x(D+Q), y(D+Q)). Same chain as ECChain with: the lookup / unlookup segments (ROBDD over
    k), the wrapper contracts with k-dependent constants (the looked-up registers are substituted by their table
    values: fact 'tab'), register subtractions / additions instead of constant additions. Assumptions as ECChain,
    for every k."""

    def __init__(self, p, a, b, table, w, D, hints, lifted=None, inputs=None, progress=None, budget=900.0):
        self.p, self.a, self.b = p, a, b
        self.n = n = p.bit_length()
        self.red = 'none'
        self.w = w
        self.table = [(int(T[0]), int(T[1])) for T in table]
        self.D = (int(D[0]), int(D[1]))
        self.budget = budget
        self._written = None
        self.alg = []
        self.nmirror = 0
        if len(self.table) != (1 << w) or w < 1:
            raise ProofFailure('the table must have 2^w entries')
        X, Y, K = list(range(n)), list(range(n, 2 * n)), list(range(2 * n, 2 * n + w))
        init = None
        if lifted is None:
            inputs = list(inputs)
            init = {}
            for j, (P, k) in enumerate(inputs):
                x, y = (0, 0) if P is None else P
                v = x | (y << n) | (k << (2 * n))
                for i in range(2 * n + w):
                    if (v >> i) & 1:
                        init[i] = init.get(i, 0) | (1 << j)
        Walker.__init__(self, hints, X + Y + K, lifted=lifted, init=init, labels=inputs, progress=progress)
        self.X, self.Y, self.K = X, Y, K
        self.tabs = []
        for Q in self.table:
            Q2 = ec_add(Q, Q, p, a)
            if Q2 is None:
                raise ProofFailure('2Q = O for a table point')
            DQ = ec_add(self.D, Q, p, a)
            if DQ is None:
                raise ProofFailure('D + Q = O for a table point')
            mQ, m2Q = ec_neg(Q, p), ec_neg(Q2, p)
            self.tabs.append(dict(Q=Q, Q2=Q2, mQ=mQ, m2Q=m2Q, DQ=DQ, E=[None, Q, mQ, m2Q], img=[Q, Q2, None, mQ],
                                  vals=dict(xq=Q[0], yq=Q[1], ymq=mQ[1], x3q=3 * Q[0] % p, x2q=Q2[0], y2q=Q2[1],
                                            ym2q=m2Q[1], xdq=DQ[0], ydq=DQ[1])))

    def run(self):
        n, p, w = self.n, self.p, self.w
        need = self.need
        X, Y, K, D, tabs = self.X, self.Y, self.K, self.D, self.tabs
        NK = 1 << w
        code = lambda P: (0, 0) if P is None else P
        bit = lambda al, v, i: al.T if (v >> i) & 1 else al.F
        dom = ('ecdomw', tuple(X), tuple(Y), tuple(K), p, tuple((t['Q'], t['Q2']) for t in tabs))
        self.assume_fact(dom, 'D-input')
        ltX, ltY = ('lt', tuple(X), p), ('lt', tuple(Y), p)
        # ---- lookup of the nine fields of T[k]
        pos, info = self.mark('lk')
        R = {f: self.fresh_list(info.get(f)) for f in WFIELDS}
        need(all(len(R[f]) == n for f in WFIELDS), 'certificate: bad lookup registers')
        Rall = [x for f in WFIELDS for x in R[f]]
        consts = tuple(sum(t['vals'][f] << (n * i) for i, f in enumerate(WFIELDS)) for t in tabs)
        tab = ('tab', tuple(K), tuple(Rall), consts)
        NR = len(Rall)

        def mux(al, val):
            return _muxbits(al, tuple(val(k) for k in K), consts, NR)[0]

        ltR = [('lt', tuple(R[f]), p) for f in WFIELDS]

        def spec(al, val):
            m = mux(al, val)
            return {Rall[i]: m[i] for i in range(NR)}, [], [], [tab] + ltR, []
        self.seg('lk', None, pos, Rall + K, spec)
        # the looked-up registers as functions of k (parametrised precondition: fact tab). ROBDD order of the
        # wrapper segments: Rall (bottom), X / Y interleaved, flags, K (top); the list follows the wire order, so
        # every new precondition sits above the conjunction of the earlier ones
        idx = {x: i for i, x in enumerate(Rall)}

        def sbf(fields):
            ws = [x for f in fields for x in R[f]]
            return [(x, (lambda al, val, i=idx[x]: mux(al, val)[i]), None, tab,
                     (lambda al, val, x=x, i=idx[x]: al.not_(al.xor(val(x), mux(al, val)[i]))))
                    for x in sorted(ws, key=lambda x: idx[x])]

        def reg(al, val, f):
            """the n bits of the looked-up field f as functions of k."""
            m = mux(al, val)
            o = n * WFIELDS.index(f)
            return m[o:o + n]
        codes = lambda what: tuple(tuple(code(what(t, j)) for t in tabs) for j in range(4))
        # ---- wrapper, part 1a: the flags f_j = [P == E_j(k)]
        pos, info = self.mark('w_flags')
        F = [self.fresh(info.get(k)) for k in ('fO', 'fQ', 'fmQ', 'fm2Q')]
        oh = ('onehot0', tuple(F))
        fimp = ('fimp', tuple(F), tuple(X), tuple(Y), tuple(K), codes(lambda t, j: t['E'][j]))

        def spec(al, val):
            x = [val(v) for v in X]
            y = [val(v) for v in Y]
            eqk = _muxbits(al, tuple(val(k) for k in K), consts, NR)[1]
            posts = {}
            for j in range(4):
                hit = al.F
                for k in range(NK):
                    cx_, cy_ = code(tabs[k]['E'][j])
                    hit = al.or_(hit, al.and_(eqk[k], al.and_(_eqc(al, x, cx_), _eqc(al, y, cy_))))
                posts[F[j]] = hit
            return posts, [], [dom], [oh, fimp, ltX, ltY], []
        self.seg('w_flags', None, pos, Rall + _inter(X, Y) + F + K, spec, sbf(('xq', 'yq', 'x2q', 'ym2q')))
        # ---- wrapper, part 1b: substitution of the safe point D (literal contract; flags are free variables that
        #      satisfy fimp and onehot0)
        pos, info = self.mark('w_pre')
        gq = self.fresh(info.get('gq'))
        gqf = ('xor3', (gq, F[1], F[2]))

        def spec(al, val):
            x = [val(v) for v in X]
            y = [val(v) for v in Y]
            fO, fQ, fmQ, fm2Q = (val(v) for v in F)
            g = al.xor(fQ, fmQ)
            xq, yq, ymq, x2q, ym2q = (reg(al, val, f) for f in ('xq', 'yq', 'ymq', 'x2q', 'ym2q'))
            A, XR = al.and_, al.xor
            posts = {gq: XR(val(gq), g)}
            for i in range(n):
                dx_, dy_ = bit(al, D[0], i), bit(al, D[1], i)
                vx = XR(x[i], A(fO, dx_))
                vx = XR(vx, A(g, XR(xq[i], dx_)))
                vx = XR(vx, A(fm2Q, XR(x2q[i], dx_)))
                vy = XR(y[i], A(fO, dy_))
                vy = XR(vy, XR(A(fQ, yq[i]), A(fmQ, ymq[i])))
                vy = XR(vy, A(g, dy_))
                vy = XR(vy, A(fm2Q, XR(ym2q[i], dy_)))
                posts[X[i]] = vx
                posts[Y[i]] = vy
            return posts, [], [ltX, ltY, oh, fimp], [ltX, ltY, oh, gqf], []
        self.seg('w_pre', None, pos, Rall + _inter(X, Y) + F + [gq] + K, spec,
                 sbf(('xq', 'yq', 'ymq', 'x2q', 'ym2q')))
        # ---- the generic addition with the looked-up Q
        self._msub(R['xq'], X)
        self._msub(R['yq'], Y)
        pos, info = self.mark('lam')
        Lr = self.fresh_list(info.get('Lr'))
        need(len(Lr) == n, 'certificate: bad lambda register')
        self._cau(lambda: self._div(Y, X), lambda P: self._xor(P, Lr, True))
        self._cau(lambda: self._pmul(Lr, X), lambda P: self._madd(P, Y))
        self._cau(lambda: self._pmul(Lr, Lr), lambda P: (self._neg(X), self._madd(P, X)))
        self._msub(R['x3q'], X)
        self._cau(lambda: self._pmul(Lr, X), lambda P: self._madd(P, Y))
        self._neg(X)
        self._cau(lambda: self._div(Y, X), lambda P: self._xor(P, Lr, False))
        self._neg(X)
        self._madd(R['xq'], X, 'madd')
        self._madd(R['ymq'], Y, 'madd')
        self.disown(Lr)
        # ---- wrapper, part 2a: images of the flagged inputs (literal contract); G3, G4
        pos, _ = self.mark('w_fix')
        cpk = dom[5]
        fiff = ('fiff', tuple(F), tuple(X), tuple(Y), tuple(K), codes(lambda t, j: t['img'][j]))

        def spec(al, val):
            x = [val(v) for v in X]
            y = [val(v) for v in Y]
            f = [val(v) for v in F]
            fO, fQ, fmQ, fm2Q = f
            eqk = _muxbits(al, tuple(val(k) for k in K), consts, NR)[1]
            A, XR = al.and_, al.xor
            anyf = al.F
            for t in f:
                anyf = al.or_(anyf, t)
            hitD = al.F
            good = al.T
            eqO = al.and_(_eqc(al, x, 0), _eqc(al, y, 0))
            for k in range(NK):
                t = tabs[k]
                hitD = al.or_(hitD, al.and_(eqk[k], al.and_(_eqc(al, x, t['DQ'][0]), _eqc(al, y, t['DQ'][1]))))
                bad = eqO
                for P in (t['Q'], t['mQ'], t['Q2']):
                    bad = al.or_(bad, al.and_(_eqc(al, x, P[0]), _eqc(al, y, P[1])))
                ok = al.and_(_ecdom_bit(al, x, y, p, cpk[k]), al.not_(bad))
                good = al.and_(good, al.or_(al.not_(eqk[k]), ok))
            a1 = al.or_(al.not_(anyf), hitD)
            a2 = al.or_(anyf, good)
            xdq, ydq, xq, yq, ymq, x2q, y2q = (reg(al, val, g) for g in ('xdq', 'ydq', 'xq', 'yq', 'ymq', 'x2q', 'y2q'))
            fa = XR(val(gq), XR(fO, fm2Q))                   # the circuit's "any flag": gq xor fO xor fm2Q
            h = XR(fO, fm2Q)
            posts = {gq: al.F}
            for i in range(n):
                vx = XR(x[i], A(fa, xdq[i]))
                vx = XR(vx, A(h, xq[i]))
                vx = XR(vx, A(fQ, x2q[i]))
                vy = XR(y[i], A(fa, ydq[i]))
                vy = XR(vy, XR(A(fO, yq[i]), A(fm2Q, ymq[i])))
                vy = XR(vy, A(fQ, y2q[i]))
                posts[X[i]] = vx
                posts[Y[i]] = vy
            return (posts, [('G3-flagged-result', a1), ('G4-generic-result', a2)], [oh, gqf],
                    [ltX, ltY, dom, fiff], [])
        self.seg('w_fix', None, pos, Rall + _inter(X, Y) + F + [gq] + K, spec,
                 sbf(('xq', 'yq', 'ymq', 'x2q', 'y2q', 'xdq', 'ydq')))
        self.disown([gq])
        # ---- wrapper, part 2b: the flags are erased from the result (f_j <-> R == img_j(k), fact fiff)
        pos, _ = self.mark('w_post')

        def spec(al, val):
            return {v: al.F for v in F}, [], [dom, fiff], [ltX, ltY], []
        self.seg('w_post', None, pos, Rall + _inter(X, Y) + F + K, spec, sbf(('xq', 'yq', 'x2q', 'y2q')))
        self.disown(F)
        # ---- unlookup
        pos, _ = self.mark('unlk')

        def spec(al, val):
            return {Rall[i]: al.F for i in range(NR)}, [], [], [], []
        self.seg('unlk', None, pos, Rall + K, spec, sbf(WFIELDS))
        self.disown(Rall)
        self.finish(())
        if self.prove:
            dirty = set(range(self.nq)) - self.Z - set(X) - set(Y) - set(K)
            need(not dirty, 'ancilla wires are not proved to end at |0> (%d wires, e.g. %s)'
                 % (len(dirty), sorted(dirty)[:4]))
        return self


# ====================================================================== drivers
def _two_torsion(p, a, b):
    """True iff x^3 + a x + b has a root in F_p (a point with y = 0):  gcd(x^p - x, f) != 1.  Own implementation
    (polynomials of degree < 3 modulo the monic cubic f), independent of the circuit generator."""
    if p < 2000:
        return any((x * x * x + a * x + b) % p == 0 for x in range(p))
    f0, f1 = b % p, a % p                                   # x^3 = -(f0 + f1 x)

    def mul(u, v):
        r = [0] * 5
        for i in range(3):
            if u[i]:
                for j in range(3):
                    r[i + j] = (r[i + j] + u[i] * v[j]) % p
        for d in (4, 3):
            c = r[d]
            if c:
                r[d - 3] = (r[d - 3] - c * f0) % p
                r[d - 2] = (r[d - 2] - c * f1) % p
                r[d] = 0
        return r[:3]
    res, base, e = [1, 0, 0], [0, 1, 0], p
    while e:
        if e & 1:
            res = mul(res, base)
        base = mul(base, base)
        e >>= 1
    g = [res[0], (res[1] - 1) % p, res[2]]                  # x^p - x mod f
    u, v = [f0, f1, 0, 1], g

    def trim(w):
        while w and w[-1] == 0:
            w.pop()
        return w
    u, v = trim(list(u)), trim(list(v))
    while v:
        inv = pow(v[-1], -1, p)
        w = list(u)
        while len(w) >= len(v):
            c = w[-1] * inv % p
            sft = len(w) - len(v)
            for i, t in enumerate(v):
                w[sft + i] = (w[sft + i] - c * t) % p
            trim(w)
        u, v = v, w
    return len(u) > 1


def check_parameters(p, a, b, Q, D):
    """the numeric side conditions of Theorem T (PROOF.md section 10). -> None, or the text of the violated one."""
    from . import primality
    if not all(isinstance(v, int) and not isinstance(v, bool) for v in (p, a, b)) or p < 5 or p % 2 == 0:
        return 'p must be an odd prime >= 5'
    if primality(p)[0] == 'composite':
        return 'p is not prime'
    if not (0 <= a < p and 0 <= b < p):
        return 'a, b must be reduced mod p'
    if (4 * a ** 3 + 27 * b * b) % p == 0:
        return 'singular curve'
    if _two_torsion(p, a, b):
        return 'the curve has a point with y = 0 (even group order): O = (0, 0) is not representable'
    for nm, P in (('Q', Q), ('D', D)):
        if P is None or len(P) != 2 or not on_curve(tuple(P), p, a, b):
            return '%s is not a point of the curve' % nm
    Q, D = tuple(Q), tuple(D)
    Q2 = ec_add(Q, Q, p, a)
    if Q2 is None or Q2[0] == Q[0]:
        return 'the order of Q must be greater than 4 (2Q = +-Q or O)'
    if D[0] == Q[0] or D == ec_neg(Q2, p):
        return 'D is one of the exceptional points Q, -Q, -2Q'
    DQ = ec_add(D, Q, p, a)
    if DQ is None or DQ[0] == Q[0]:
        return 'D + Q is one of Q, -Q, O'
    return None


def sample_points(p, a, b, Q, D, count, seed=1):
    """points for the sampled tests: all points and O if p is small, else O, the exceptional points, D, D + Q and
    their negatives, small multiples of Q and random multiples of Q."""
    if p <= 4096:
        return [None] + curve_points(p, a, b)
    rnd = random.Random(seed)
    Q, D = tuple(Q), tuple(D)
    Q2 = ec_add(Q, Q, p, a)
    DQ = ec_add(D, Q, p, a)
    pts = [None, Q, ec_neg(Q, p), Q2, ec_neg(Q2, p), D, ec_neg(D, p), DQ, ec_neg(DQ, p)]
    R = Q2
    for _ in range(6):
        R = ec_add(R, Q, p, a)
        pts += [R, ec_neg(R, p)]
    pts = list(dict.fromkeys(pts))
    while len(pts) < count:
        pts.append(ec_mul(rnd.getrandbits(p.bit_length() + 8) + 5, Q, p, a))
    return pts[:max(count, 9)]


def ecadd_hints(p, a, b, Q, karatsuba_levels=3, D=None):
    """the (untrusted) certificate of the library circuit for these parameters."""
    from ...crypto import ecpoint
    return ecpoint.ec_point_add(p, a, b, tuple(Q), karatsuba_levels=karatsuba_levels, D=D).hints()


def _circuit_run(L, n, pts, p, a, Q):
    """classical run of the lifted circuit on the points. -> None or dict(point, reason)."""
    from .walker import simulate_flat
    S = len(pts)
    init = {}
    for j, P in enumerate(pts):
        x, y = (0, 0) if P is None else P
        for i in range(n):
            if (x >> i) & 1:
                init[i] = init.get(i, 0) | (1 << j)
            if (y >> i) & 1:
                init[n + i] = init.get(n + i, 0) | (1 << j)
    v, bad = simulate_flat(L, init, S)
    if bad is not None:
        j = (bad[2] & -bad[2]).bit_length() - 1
        return dict(point=pts[j], reason='op %d: %s' % (bad[0], bad[1]))
    for w in range(2 * n, L.nq):
        if v[w]:
            j = (v[w] & -v[w]).bit_length() - 1
            return dict(point=pts[j], reason='wire %d is not returned to zero' % w)
    for j, P in enumerate(pts):
        gx = sum(((v[i] >> j) & 1) << i for i in range(n))
        gy = sum(((v[n + i] >> j) & 1) << i for i in range(n))
        R = ec_add(P, Q, p, a)
        if (gx, gy) != ((0, 0) if R is None else R):
            return dict(point=P, reason='wrong result on the output wires')
    return None


def wec_test(p, a, b, Q, hints=None, points=None, count=64, lifted=None, karatsuba_levels=3):
    """execute the word-level program W_ec (the contract chain itself, run mode) on the points and check every
    mathematical assumption (D-input, G1..G4, M1..M6), every carried fact and the result (= P + Q, every other wire
    0). With `lifted`: also run the circuit itself classically on the same points (result and clean ancillas).
    -> dict(ok, inputs, violations[, circuit_agrees])."""
    Q = tuple(Q)
    if hints is None:
        hints = ecadd_hints(p, a, b, Q, karatsuba_levels)
    D = tuple(hints['D'])
    if points is None:
        points = sample_points(p, a, b, Q, D, count)
    points = list(points)
    n = p.bit_length()
    viol = []
    try:
        ch = ECChain(p, a, b, Q, D, hints, points=points).run()
        viol = list(ch.viol)
        X, Y = ch.X, ch.Y
        gx, gy = ch._get(X), ch._get(Y)
        io = set(X) | set(Y)
        for j, P in enumerate(points):
            R = ec_add(P, Q, p, a)
            if (gx[j], gy[j]) != ((0, 0) if R is None else R) and len(viol) < 20:
                viol.append(dict(segment='final', what='result (%d, %d)' % (gx[j], gy[j]), x=P))
        for w, m in ch.store.items():
            if m and w not in io and len(viol) < 20:
                viol.append(dict(segment='final', what='wire %d is not 0' % w, x=points[(m & -m).bit_length() - 1]))
    except ProofFailure as e:
        viol.append(dict(segment='chain', what=str(e), x=None))
    res = dict(ok=not viol, inputs=len(points), violations=viol[:20])
    if lifted is not None:
        cex = _circuit_run(lifted, n, points, p, a, Q)
        res['circuit_agrees'] = cex is None
        if cex is not None:
            res['circuit_counterexample'] = cex
        res['ok'] = res['ok'] and cex is None
    return res


def verify_ecadd(qasm, p, a, b, Q, hints=None, lifted=None, refute_samples=64, test_samples=64, progress=None,
                 karatsuba_levels=3, budget=900.0):
    """Compositional proof that the OpenQASM circuit `qasm` (path or text) is the exact point addition P -> P + Q on
    y^2 = x^3 + a x + b over F_p: the wires 0..n-1 hold x and n..2n-1 hold y (n = p.bit_length(); O is (0, 0)), all
    other wires |0>;  for EVERY point P of the curve and for O the circuit ends with P + Q on the same wires and
    every other wire back at |0>, every AND target was |0> and every measurement-based erasure is exact (no phase,
    for every outcome string).  P + Q is the textbook affine addition ec_add (the specification).
    hints: the certificate (ECPointAddCircuit.hints(), untrusted; also gives the safe point D, which is checked);
    default: regenerated from the library (karatsuba_levels).
    -> dict: status 'PROVED' | 'REFUTED' (with a point) | 'UNKNOWN' | 'ERROR' (bad parameters), stage, detail,
       assumptions (name -> number of segments using it), segments, algebra (sub-goals), mirrors, time, ...
    PROVED means: machine-proved refinement of the word-level program W_ec by the circuit (module docstring), under
    D-input and the assumptions G1..G4, M1..M6 that Theorem T / the Kaliski lemmas discharge (PROOF.md 9-10)."""
    import time
    from . import primality, check_gadget_lemmas, lift, LiftError
    t0 = time.time()
    n = p.bit_length() if isinstance(p, int) else 0
    res = dict(status='UNKNOWN', p=p, a=a, b=b, Q=list(Q) if Q else Q, n=n, stage='input')

    def done(**kw):
        res.update(kw)
        res['time'] = round(time.time() - t0, 2)
        return res
    if hints is None:
        try:
            hints = ecadd_hints(p, a, b, Q, karatsuba_levels)
        except (ValueError, TypeError) as e:
            return done(status='ERROR', detail='parameters: %s' % e)
    D = hints.get('D')
    bad = check_parameters(p, a, b, Q, D)
    if bad:
        return done(status='ERROR', detail='parameters: ' + bad)
    Q, D = tuple(Q), tuple(D)
    pv, phow = primality(p)
    res['p_primality'] = '%s: %s' % (pv, phow)
    res['D'] = list(D)
    res['stage'] = 'gadgets'
    try:
        check_gadget_lemmas()
    except AssertionError as e:
        return done(detail='gadget lemma failed: %r' % (e,))
    res['stage'] = 'lift'
    try:
        L = lifted if lifted is not None else lift(qasm)
    except LiftError as e:
        return done(detail='not a composition of the AND / ERASE gadgets: %s' % e)
    res.update(ops=L.n, qubits=L.nq, measurements=L.nmeas, t_lift=round(time.time() - t0, 2))
    res['stage'] = 'chain'
    try:
        if not isinstance(hints.get('marks'), list):
            raise ProofFailure('certificate: no marks')
        if L.nq < 2 * n:
            raise ProofFailure('interface: fewer than 2n wires')
        ch = ECChain(p, a, b, Q, D, hints, lifted=L, progress=progress, budget=budget).run()
    except ProofFailure as e:
        det = {k: v for k, v in e.detail.items() if k != 'entry_counterexample'}
        cex = None
        if refute_samples:
            cex = _circuit_run(L, n, sample_points(p, a, b, Q, D, refute_samples, seed=7), p, a, Q)
        if cex is not None:
            return done(status='REFUTED', detail=str(e), where=det, counterexample=cex)
        return done(detail=str(e), where=det)
    except (KeyError, TypeError, IndexError, AttributeError) as e:      # malformed certificate
        return done(detail='certificate: malformed (%s: %s)' % (type(e).__name__, e))
    res.update(segments=ch.nseg, mirrors=ch.nmirror, assumptions=dict(sorted(ch.assumed.items())),
               lemma_facts=sorted(set(ch.admitted)), algebra=ch.alg, bdd_nodes_total=ch.nodes,
               bdd_nodes_max=ch.maxnodes)
    if test_samples:
        res['stage'] = 'wec_test'
        t = wec_test(p, a, b, Q, hints=hints, count=test_samples, lifted=L)
        res['wec_test'] = t
        if not t['ok']:
            if not t.get('circuit_agrees', True):
                return done(status='REFUTED', detail='the circuit is wrong on a concrete point',
                            counterexample=t.get('circuit_counterexample'))
            return done(detail='an assumption of W_ec is violated on a concrete point: %r' % (t['violations'][:1],))
    note = '' if pv == 'prime' else ' PROVED under the condition that p is prime (only a probabilistic test was run).'
    return done(status='PROVED', stage='done',
                proved_modulo='D-input (domain), G1..G4 (Theorem T, PROOF.md section 10: elementary curve algebra, '
                              'tested by wec_test) and the Kaliski lemmas M1..M6 (Z3, PROOF.md sections 4-5)',
                detail='every segment of the circuit is machine-proved against its contract at full width: ROBDD '
                       'segments, algebra sub-goals (products, fold reductions), mirror segments.' + note)


# ---------------------------------------------------------------------- window variant: drivers
def check_window_parameters(p, a, b, table, w, D):
    """the side conditions of Theorem T for every table entry. -> None or the text of the violated one."""
    if not isinstance(w, int) or isinstance(w, bool) or not 1 <= w <= 8 or len(table) != (1 << w):
        return 'the table must have 2^w entries, 1 <= w <= 8'
    for k, T in enumerate(table):
        bad = check_parameters(p, a, b, T, D)
        if bad:
            return 'table entry %d: %s' % (k, bad)
    return None


def sample_window_inputs(p, a, b, table, D, count, seed=1):
    """inputs (P, k): all pairs if p is small; else for every k the exceptional points of T[k], D, D + T[k] and
    their negatives, and random pairs."""
    NK = len(table)
    if p <= 600:
        return [(P, k) for k in range(NK) for P in [None] + curve_points(p, a, b)]
    rnd = random.Random(seed)
    out = []
    for k, Q in enumerate(table):
        Q2 = ec_add(Q, Q, p, a)
        DQ = ec_add(D, Q, p, a)
        out += [(P, k) for P in (None, Q, ec_neg(Q, p), Q2, ec_neg(Q2, p), D, ec_neg(DQ, p))]
    while len(out) < count:
        k = rnd.randrange(NK)
        out.append((ec_mul(rnd.getrandbits(p.bit_length() + 8) + 5, table[rnd.randrange(NK)], p, a), k))
    return out


def ecadd_window_hints(p, a, b, table, w, karatsuba_levels=3, D=None):
    from ...crypto import ecpoint
    return ecpoint.ec_point_add_window(p, a, b, [tuple(T) for T in table], w, karatsuba_levels=karatsuba_levels,
                                       D=D).hints()


def _circuit_run_window(L, n, w, inputs, p, a, table):
    from .walker import simulate_flat
    init = {}
    for j, (P, k) in enumerate(inputs):
        x, y = (0, 0) if P is None else P
        v = x | (y << n) | (k << (2 * n))
        for i in range(2 * n + w):
            if (v >> i) & 1:
                init[i] = init.get(i, 0) | (1 << j)
    v, bad = simulate_flat(L, init, len(inputs))
    if bad is not None:
        j = (bad[2] & -bad[2]).bit_length() - 1
        return dict(input=inputs[j], reason='op %d: %s' % (bad[0], bad[1]))
    for q in range(2 * n + w, L.nq):
        if v[q]:
            j = (v[q] & -v[q]).bit_length() - 1
            return dict(input=inputs[j], reason='wire %d is not returned to zero' % q)
    for j, (P, k) in enumerate(inputs):
        got = sum(((v[i] >> j) & 1) << i for i in range(2 * n + w))
        R = ec_add(P, table[k], p, a)
        x, y = (0, 0) if R is None else R
        if got != (x | (y << n) | (k << (2 * n))):
            return dict(input=(P, k), reason='wrong result on the output wires')
    return None


def wec_test_window(p, a, b, table, w, hints=None, inputs=None, count=160, lifted=None, karatsuba_levels=3):
    """wec_test for the windowed addition: the word-level program on the inputs (P, k)."""
    table = [tuple(T) for T in table]
    if hints is None:
        hints = ecadd_window_hints(p, a, b, table, w, karatsuba_levels)
    D = tuple(hints['D'])
    if inputs is None:
        inputs = sample_window_inputs(p, a, b, table, D, count)
    inputs = list(inputs)
    n = p.bit_length()
    viol = []
    try:
        ch = ECWinChain(p, a, b, table, w, D, hints, inputs=inputs).run()
        viol = list(ch.viol)
        gx, gy, gk = ch._get(ch.X), ch._get(ch.Y), ch._get(ch.K)
        io = set(ch.X) | set(ch.Y) | set(ch.K)
        for j, (P, k) in enumerate(inputs):
            R = ec_add(P, table[k], p, a)
            if (gx[j], gy[j], gk[j]) != ((0, 0) if R is None else R) + (k,) and len(viol) < 20:
                viol.append(dict(segment='final', what='result (%d, %d, k=%d)' % (gx[j], gy[j], gk[j]), x=(P, k)))
        for q, m in ch.store.items():
            if m and q not in io and len(viol) < 20:
                viol.append(dict(segment='final', what='wire %d is not 0' % q, x=inputs[(m & -m).bit_length() - 1]))
    except ProofFailure as e:
        viol.append(dict(segment='chain', what=str(e), x=None))
    res = dict(ok=not viol, inputs=len(inputs), violations=viol[:20])
    if lifted is not None:
        cex = _circuit_run_window(lifted, n, w, inputs, p, a, table)
        res['circuit_agrees'] = cex is None
        if cex is not None:
            res['circuit_counterexample'] = cex
        res['ok'] = res['ok'] and cex is None
    return res


def verify_ecadd_window(qasm, p, a, b, table, w, hints=None, lifted=None, refute_samples=160, test_samples=160,
                        progress=None, karatsuba_levels=3, budget=900.0):
    """As verify_ecadd, for the windowed addition (P, k) -> (P + table[k], k) of ec_point_add_window (default
    construction): wires 0..n-1 x, n..2n-1 y, 2n..2n+w-1 the index k; for EVERY point P of the curve and O and EVERY
    k in [0, 2^w). Every table entry must satisfy the conditions of Theorem T with the certificate's D (order > 4,
    D safe for every entry): checked."""
    import time
    from . import primality, check_gadget_lemmas, lift, LiftError
    t0 = time.time()
    n = p.bit_length() if isinstance(p, int) else 0
    res = dict(status='UNKNOWN', p=p, a=a, b=b, w=w, n=n, stage='input')

    def done(**kw):
        res.update(kw)
        res['time'] = round(time.time() - t0, 2)
        return res
    try:
        table = [(int(T[0]), int(T[1])) for T in table]
    except (TypeError, IndexError, ValueError):
        return done(status='ERROR', detail='parameters: bad table')
    if hints is None:
        try:
            hints = ecadd_window_hints(p, a, b, table, w, karatsuba_levels)
        except (ValueError, TypeError) as e:
            return done(status='ERROR', detail='parameters: %s' % e)
    D = hints.get('D')
    bad = check_window_parameters(p, a, b, table, w, D)
    if bad:
        return done(status='ERROR', detail='parameters: ' + bad)
    D = tuple(D)
    pv, phow = primality(p)
    res['p_primality'] = '%s: %s' % (pv, phow)
    res['D'] = list(D)
    res['stage'] = 'gadgets'
    try:
        check_gadget_lemmas()
    except AssertionError as e:
        return done(detail='gadget lemma failed: %r' % (e,))
    res['stage'] = 'lift'
    try:
        L = lifted if lifted is not None else lift(qasm)
    except LiftError as e:
        return done(detail='not a composition of the AND / ERASE gadgets: %s' % e)
    res.update(ops=L.n, qubits=L.nq, measurements=L.nmeas, t_lift=round(time.time() - t0, 2))
    res['stage'] = 'chain'
    try:
        if not isinstance(hints.get('marks'), list):
            raise ProofFailure('certificate: no marks')
        if L.nq < 2 * n + w:
            raise ProofFailure('interface: fewer than 2n + w wires')
        ch = ECWinChain(p, a, b, table, w, D, hints, lifted=L, progress=progress, budget=budget).run()
    except ProofFailure as e:
        det = {k: v for k, v in e.detail.items() if k != 'entry_counterexample'}
        cex = None
        if refute_samples:
            cex = _circuit_run_window(L, n, w, sample_window_inputs(p, a, b, table, D, refute_samples, seed=7),
                                      p, a, table)
        if cex is not None:
            return done(status='REFUTED', detail=str(e), where=det, counterexample=cex)
        return done(detail=str(e), where=det)
    except (KeyError, TypeError, IndexError, AttributeError) as e:
        return done(detail='certificate: malformed (%s: %s)' % (type(e).__name__, e))
    res.update(segments=ch.nseg, mirrors=ch.nmirror, assumptions=dict(sorted(ch.assumed.items())),
               lemma_facts=sorted(set(ch.admitted)), algebra=ch.alg, bdd_nodes_total=ch.nodes,
               bdd_nodes_max=ch.maxnodes)
    if test_samples:
        res['stage'] = 'wec_test'
        t = wec_test_window(p, a, b, table, w, hints=hints, count=test_samples, lifted=L)
        res['wec_test'] = t
        if not t['ok']:
            if not t.get('circuit_agrees', True):
                return done(status='REFUTED', detail='the circuit is wrong on a concrete input',
                            counterexample=t.get('circuit_counterexample'))
            return done(detail='an assumption of W_ec is violated on a concrete input: %r' % (t['violations'][:1],))
    note = '' if pv == 'prime' else ' PROVED under the condition that p is prime (only a probabilistic test was run).'
    return done(status='PROVED', stage='done',
                proved_modulo='D-input (domain), G1..G4 (Theorem T for every table entry, PROOF.md section 10) and '
                              'the Kaliski lemmas M1..M6 (Z3, PROOF.md sections 4-5)',
                detail='every segment of the circuit is machine-proved against its contract at full width: ROBDD '
                       'segments (incl. the table lookup), algebra sub-goals, mirror segments.' + note)
