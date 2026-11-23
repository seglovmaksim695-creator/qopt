"""lemmas_z3 -- the mathematical assumptions of the ``chain`` client as Z3 conditions, for ALL N, p, x.

Claimed (every implication is an unsat of the negation): for all integers N, all odd p with 2^(nb-1) < p < 2^nb, p < 2^N,
p prime (ONE instance of primality is used), all 0 <= x < p, by induction on the iteration number k:
    INIT    Inv(0) for the initial state (u, v, r, s, f) = (p, x, 0, 1, 1)
    STEP    Inv(k), k <= 2N-1  ==>  all assumptions of step k (M-dbl, M-a, M-s, M-carry), cross checks
            (no borrow, parity before the shift, clearing of a and b, all words < 2^N) and Inv(k+1)
    FINAL   Inv(2N) ==> M-end-vf (v = 0, f = 0), M-end, and the result: x*(p - r) = 2^(2N) (mod p), 1 <= p - r <= p-1
The iteration step is written exactly as the contracts of ``chain`` give it (the functions below; the same text is
``pq.step``, which is compared with the execution of the contracts).  Axioms: pow2(0) = 1, pow2(i+1) = 2 pow2(i), pow2
monotone.  Auxiliary instances (identities and inequalities of integers) are proved by separate LEMMA queries.
    python -m qoptct.foreign.qualtran.kaliski.lemmas_z3 [--timeout S] [--dump DIR] [--out FILE]
"""
import hashlib, itertools, json, os, shutil, subprocess, sys, tempfile, time

BASE = """(set-logic ALL)
(declare-fun pow2 (Int) Int)
(assert (= (pow2 0) 1))
(assert (forall ((i Int)) (! (=> (>= i 0) (= (pow2 (+ i 1)) (* 2 (pow2 i)))) :pattern ((pow2 i)))))
(assert (forall ((i Int) (j Int)) (! (=> (and (<= 0 i) (<= i j)) (<= (pow2 i) (pow2 j))) :pattern ((pow2 i) (pow2 j)))))
(declare-const N Int) (declare-const nb Int) (declare-const p Int) (declare-const x Int)
(assert (>= nb 2)) (assert (>= N 2))
(assert (= (mod p 2) 1))
(assert (< (pow2 (- nb 1)) p)) (assert (< p (pow2 nb))) (assert (< p (pow2 N)))
(assert (<= 0 x)) (assert (< x p))
(define-fun odd ((t Int)) Bool (= (mod t 2) 1))
"""
NB_LEMMA = '(assert (<= nb N))\n'          # proved separately (query NB)

# the invariant, written once as templates (so that its conjuncts can also be stated one by one)
INV_COMMON = ['(<= 0 {k})', '(<= {k} (* 2 N))', '(= {T} (pow2 {k}))']
INV_ACT = ['{f}', '(>= x 1)', '(>= {u} 1)', '(>= {v} 1)', '(>= {r} 0)', '(>= {s} 1)', '(<= {u} p)', '(< {v} p)',
           '(= (+ (* {u} {s}) (* {v} {r})) p)',
           '(=> (= {r} 0) (and (= {u} p) (= {s} 1)))',
           '(= (+ (* x {r}) (* {u} {T})) (* p {cR}))',
           '(= (- (* x {s}) (* {v} {T})) (* p {cS}))',
           '(<= 1 {ga})', '(<= {ga} nb)', '(<= 1 {gb})', '(<= {gb} nb)', '(= (+ {ga} {gb}) (- (* 2 nb) {k}))',
           '(< {u} (pow2 {ga}))', '(< {v} (pow2 {gb}))']
INV_DET = ['{f}', '(>= x 1)', '(= {v} 0)', '(= {u} 1)', '(= {s} p)', '(<= 1 {r})', '(< {r} p)', '(<= {k} (- (* 2 nb) 1))',
           '(= (+ (* x {r}) {T}) (* p {cR}))']
INV_POST = ['(not {f})', '(>= x 1)', '(= {v} 0)', '(= {u} 1)', '(= {s} p)', '(<= 1 {r})', '(< {r} p)', '(>= {k} 1)',
            '(= (+ (* x {r}) {T}) (* p {cR}))']
INV_ZERO = ['(= x 0)', '(= {v} 0)', '(= {u} p)', '(= {s} 1)', '(= {r} 0)', '(= {f} (= {k} 0))']
PARAMS = ('k', 'T', 'f', 'u', 'v', 'r', 's', 'ga', 'gb', 'cR', 'cS')


def inst(templates, args):
    d = dict(zip(PARAMS, args))
    return [t.format(**d) for t in templates]


def _conj(xs):
    return '(and %s)' % ' '.join(xs)


INV = ('(define-fun inv ((k Int) (T Int) (f Bool) (u Int) (v Int) (r Int) (s Int) (ga Int) (gb Int) (cR Int) (cS Int)) Bool\n'
       '  (and %s\n   (or %s)))\n' % (' '.join(inst(INV_COMMON, PARAMS)),
                                     '\n       '.join(_conj(inst(t, PARAMS)) for t in (INV_ACT, INV_DET, INV_POST, INV_ZERO))))

# the iteration exactly as the contracts of qkal_chain (S1..S6) define it; arguments = the state (u v r s f)
_DEFS = [('z1', 'Bool', '(= v 0)'),                                   # S1
         ('m1', 'Bool', '(and f (z1 A))'),
         ('f1', 'Bool', '(xor f (m1 A))'),
         ('t1', 'Bool', '(m1 A)'),
         ('a2', 'Bool', '(and (f1 A) (not (odd u)))'),                   # S2
         ('m2', 'Bool', '(xor (m1 A) (and (f1 A) (not (odd v)) (not (a2 A))))'),
         ('b2', 'Bool', '(xor (a2 A) (m2 A))'),
         ('c3', 'Bool', '(and (> u v) (f1 A) (not (b2 A)))'),            # S3
         ('a3', 'Bool', '(xor (a2 A) (c3 A))'),
         ('m3', 'Bool', '(xor (m2 A) (c3 A))'),
         ('u4', 'Int', '(ite (a3 A) v u)'), ('v4', 'Int', '(ite (a3 A) u v)'),          # S4
         ('r4', 'Int', '(ite (a3 A) s r)'), ('s4', 'Int', '(ite (a3 A) r s)'),
         ('c5', 'Bool', '(and (f1 A) (not (b2 A)))'),                    # S5
         ('v5', 'Int', '(ite (c5 A) (- (v4 A) (u4 A)) (v4 A))'),
         ('s5', 'Int', '(ite (c5 A) (+ (s4 A) (r4 A)) (s4 A))'),
         ('b6', 'Bool', '(xor (b2 A) (xor (m3 A) (a3 A)))'),             # S6
         ('v6', 'Int', '(ite (f1 A) (div (v5 A) 2) (v5 A))'),
         ('r6', 'Int', '(ite (< (* 2 (r4 A)) p) (* 2 (r4 A)) (- (* 2 (r4 A)) p))'),
         ('u7', 'Int', '(ite (a3 A) (v6 A) (u4 A))'), ('v7', 'Int', '(ite (a3 A) (u4 A) (v6 A))'),
         ('r7', 'Int', '(ite (a3 A) (s5 A) (r6 A))'), ('s7', 'Int', '(ite (a3 A) (r6 A) (s5 A))'),
         ('a7', 'Bool', '(xor (a3 A) (not (odd (s7 A))))')]
ARGS = '(u Int) (v Int) (r Int) (s Int) (f Bool)'
DEFS = ''.join('(define-fun %s (%s) %s %s)\n' % (nm, ARGS, ty, body.replace(' A)', ' u v r s f)'))
               for nm, ty, body in _DEFS)


def S(name):
    return '(%s u v r s f)' % name


STATE = """(declare-const k Int) (declare-const T Int) (declare-const f Bool)
(declare-const u Int) (declare-const v Int) (declare-const r Int) (declare-const s Int)
(declare-const ga Int) (declare-const gb Int) (declare-const cR Int) (declare-const cS Int)
(assert (inv k T f u v r s ga gb cR cS))
"""
STEP_HYP = '(assert (<= k (- (* 2 N) 1)))\n'
PRIME = '(assert (=> (and (= (+ (* u s) (* u r)) p) (<= 1 u) (< u p)) (= u 1)))\n'    # ONE instance of primality

ACT = '(and f (>= x 1) (>= v 1))'

# instances of the lemma schemas proved below (LEMMA queries); each is a valid formula, so adding it is sound
def eq2(X, h, Z):        # L-half: 2h = X -> X*Z = 2*(h*Z)
    return '(assert (=> (= (* 2 %s) %s) (= (* %s %s) (* 2 (* %s %s)))))\n' % (h, X, X, Z, h, Z)


def dist(a, b, c, op='+'):  # L-dist: a*(b op c) = a*b op a*c ; (b op c)*a = b*a op c*a
    return ('(assert (= (* %s (%s %s %s)) (%s (* %s %s) (* %s %s))))\n(assert (= (* (%s %s %s) %s) (%s (* %s %s) (* %s %s))))\n'
            % (a, op, b, c, op, a, b, a, c, op, b, c, a, op, b, a, c, a))


def mono(a, c, b):       # L-mono: a >= c, b >= 0 -> a*b >= c*b
    return '(assert (=> (and (>= %s %d) (>= %s 0)) (>= (* %s %s) (* %d %s))))\n' % (a, c, b, a, b, c, b)


def one1(a, b):          # L-one: a >= 1, b >= 1 -> a*b >= 1
    return '(assert (=> (and (>= %s 1) (>= %s 1)) (>= (* %s %s) 1)))\n' % (a, b, a, b)


def unit(u, Z):          # L-unit: u = 1 -> u*Z = Z
    return '(assert (=> (= %s 1) (= (* %s %s) %s)))\n' % (u, u, Z, Z)


H = '(declare-const h Int)\n'
Q = '(declare-const q Int)\n'
AUX_RED = dist('x', '(* 2 r)', 'p', '-') + dist('p', '(* 2 cR)', 'x', '-')
AUX = {
    'B1': H + '(assert (= (* 2 h) u))\n' + eq2('u', 'h', 's') + eq2('u', 'h', 'T') + mono('u', 2, 's') + mono('v', 0, 'r'),
    'B2': (H + '(assert (= (* 2 h) v))\n' + eq2('v', 'h', 'r') + eq2('v', 'h', 'T') + mono('v', 2, 'r') + one1('u', 's')
           + Q + '(assert (= q (div s 2)))\n' + eq2('s', 'q', 'u')),
    'B3': (H + '(assert (= (* 2 h) (- u v)))\n' + eq2('(- u v)', 'h', 's') + dist('s', 'u', 'v', '-') + eq2('(- u v)', 'h', 'T')
           + dist('T', 'u', 'v', '-') + dist('v', 'r', 's') + dist('x', 'r', 's') + dist('p', 'cR', 'cS')
           + mono('u', 3, 's') + mono('v', 0, 'r') + mono('u', 1, 's') + mono('v', 1, 'r')),
    'B4a': (H + '(assert (= (* 2 h) (- v u)))\n' + eq2('(- v u)', 'h', 'r') + dist('r', 'v', 'u', '-') + dist('u', 's', 'r')
            + eq2('(- v u)', 'h', 'T') + dist('T', 'v', 'u', '-') + dist('x', 's', 'r') + dist('p', 'cS', 'cR')
            + mono('v', 3, 'r') + one1('u', 's') + mono('u', 1, 's') + mono('v', 1, 'r')
            + Q + '(assert (= q (div (+ s r) 2)))\n' + eq2('(+ s r)', 'q', 'u')),
    'B4b': unit('u', 's') + unit('u', 'r') + unit('u', 'T') + unit('v', 'r') + AUX_RED,
    'RED': AUX_RED,
}
CASES = [  # name, hypothesis, ghosts (ga', gb', cR', cS'), auxiliary instances
    ('ACT-B1 u even', '(and %s (not (odd u)))' % ACT, ('(- ga 1)', 'gb', 'cR', '(* 2 cS)'), AUX['B1']),
    ('ACT-B2 u odd, v even', '(and %s (odd u) (not (odd v)))' % ACT, ('ga', '(- gb 1)', '(* 2 cR)', 'cS'), AUX['B2']),
    ('ACT-B3 both odd, u > v', '(and %s (odd u) (odd v) (> u v))' % ACT, ('(- ga 1)', 'gb', '(+ cR cS)', '(* 2 cS)'), AUX['B3']),
    ('ACT-B4a both odd, u < v', '(and %s (odd u) (odd v) (< u v))' % ACT, ('ga', '(- gb 1)', '(* 2 cR)', '(+ cS cR)'), AUX['B4a']),
    ('ACT-B4b u = v, 2r < p', '(and %s (odd u) (odd v) (= u v) (< (* 2 r) p))' % ACT, ('ga', 'gb', '(* 2 cR)', 'cS'), AUX['B4b']),
    ('ACT-B4b u = v, 2r >= p', '(and %s (odd u) (odd v) (= u v) (>= (* 2 r) p))' % ACT, ('ga', 'gb', '(- (* 2 cR) x)', 'cS'), AUX['B4b']),
    ('DET 2r < p', '(and f (>= x 1) (= v 0) (< (* 2 r) p))', ('ga', 'gb', '(* 2 cR)', 'cS'), AUX['RED']),
    ('DET 2r >= p', '(and f (>= x 1) (= v 0) (>= (* 2 r) p))', ('ga', 'gb', '(- (* 2 cR) x)', 'cS'), AUX['RED']),
    ('POST 2r < p', '(and (not f) (>= x 1) (< (* 2 r) p))', ('ga', 'gb', '(* 2 cR)', 'cS'), AUX['RED']),
    ('POST 2r >= p', '(and (not f) (>= x 1) (>= (* 2 r) p))', ('ga', 'gb', '(- (* 2 cR) x)', 'cS'), AUX['RED']),
    ('ZERO x = 0', '(= x 0)', ('ga', 'gb', '0', '0'), ''),
]
EXPL = {'B1': ('h', 'v', 'r', '(* 2 s)'), 'B2': ('u', 'h', '(* 2 r)', 's'),        # explicit next states (u, v, r, s)
        'B3': ('h', 'v', '(+ r s)', '(* 2 s)'), 'B4a': ('u', 'h', '(* 2 r)', '(+ s r)')}
# the lemma schemas, proved for all integers (fresh variables)
_V = '(declare-const a Int) (declare-const b Int) (declare-const c Int) (declare-const X Int) (declare-const h Int) (declare-const Z Int)\n'
LEMMAS = [
    ('L-half: 2h = X -> X*Z = 2*(h*Z)', '(set-logic ALL)\n' + _V + '(assert (not (=> (= (* 2 h) X) (= (* X Z) (* 2 (* h Z))))))\n(check-sat)\n'),
    ('L-dist: a*(b+c) = a*b + a*c, (b+c)*a = b*a + c*a, the same with -', '(set-logic ALL)\n' + _V +
     '(assert (not (and (= (* a (+ b c)) (+ (* a b) (* a c))) (= (* (+ b c) a) (+ (* b a) (* c a))) '
     '(= (* a (- b c)) (- (* a b) (* a c))) (= (* (- b c) a) (- (* b a) (* c a))))))\n(check-sat)\n'),
    ('L-one: a >= 1, b >= 1 -> a*b >= 1', '(set-logic ALL)\n' + _V + '(assert (not (=> (and (>= a 1) (>= b 1)) (>= (* a b) 1))))\n(check-sat)\n'),
    ('L-unit: a = 1 -> a*Z = Z', '(set-logic ALL)\n' + _V + '(assert (not (=> (= a 1) (= (* a Z) Z))))\n(check-sat)\n'),
] + [('L-mono c=%d: a >= c, b >= 0 -> a*b >= c*b' % c, '(set-logic ALL)\n' + _V +
      '(assert (not (=> (and (>= a %d) (>= b 0)) (>= (* a b) (* %d b)))))\n(check-sat)\n' % (c, c)) for c in (0, 1, 2, 3)]

GOALS = [
    ('M-dbl: 0 <= r4 < p', '(and (<= 0 %s) (< %s p))' % (S('r4'), S('r4'))),
    ('M-a: a3 -> 2 r4 < p', '(=> %s (< (* 2 %s) p))' % (S('a3'), S('r4'))),
    ('M-s: not a3 -> s5 odd', '(=> (not %s) (odd %s))' % (S('a3'), S('s5'))),
    ('M-carry: c5 -> s4 + r4 < 2^N', '(=> %s (< (+ %s %s) (pow2 N)))' % (S('c5'), S('s4'), S('r4'))),
    ('X no borrow: c5 -> u4 <= v4', '(=> %s (<= %s %s))' % (S('c5'), S('u4'), S('v4'))),
    ('X v5 even when f1', '(=> %s (not (odd %s)))' % (S('f1'), S('v5'))),
    ('X b cleared', '(not %s)' % S('b6')),
    ('X a cleared', '(not %s)' % S('a7')),
    ('X all words in [0, 2^N)', '(and %s)' % ' '.join('(<= 0 %s) (< %s (pow2 N))' % (S(w), S(w))
                                                        for w in ('u7', 'v7', 'r7', 's7', 'v5', 's5'))),
]



def inv_next(gh):
    return '(inv (+ k 1) (* 2 T) %s %s %s %s %s %s)' % (S('f1'), S('u7'), S('v7'), S('r7'), S('s7'), ' '.join(gh))


def head(extra=''):
    return BASE + NB_LEMMA + INV + DEFS + STATE + extra


def queries():
    out = []
    out.append(('NB: nb <= N', BASE + '(assert (not (<= nb N)))\n(check-sat)\n'))
    for nm, txt in LEMMAS:
        out.append(('LEMMA ' + nm, txt))
    out.append(('INIT: Inv(0) for (p, x, 0, 1, f=1)', BASE + NB_LEMMA + INV +
                '(assert (not (inv 0 1 true p x 0 1 nb nb (ite (= x 0) 0 1) 0)))\n(check-sat)\n'))
    out.append(('T0: the t bit of iteration 0 is [x = 0]', BASE + NB_LEMMA + INV + DEFS +
                '(assert (not (= (t1 p x 0 1 true) (= x 0))))\n(check-sat)\n'))
    fin = BASE + NB_LEMMA + INV + DEFS + STATE.replace('(assert (inv k T', '(assert (= k (* 2 N)))\n(assert (inv k T')
    out.append(('FINAL: M-end-vf (v = 0, f = 0)', fin + '(assert (not (and (= v 0) (not f))))\n(check-sat)\n'))
    out.append(('FINAL: M-end (x = 0 -> u = p, s = 1, r = 0; x > 0 -> u = 1, s = p)', fin +
                '(assert (not (and (=> (= x 0) (and (= u p) (= s 1) (= r 0))) (=> (>= x 1) (and (= u 1) (= s p))))))\n'
                '(check-sat)\n'))
    out.append(('FINAL: x > 0 -> out = p - r in [1, p-1] and x*out = 2^(2N) (mod p)', fin +
                '(assert (not (=> (>= x 1) (and (<= 1 (- p r)) (<= (- p r) (- p 1)) (= T (pow2 (* 2 N))) '
                '(= (- (* x (- p r)) T) (* p (- x cR)))))))\n(check-sat)\n'))
    for cname, hyp, gh, aux in CASES:
        h = head(STEP_HYP + PRIME + aux + '(assert %s)\n' % hyp)
        for gname, g in GOALS:
            out.append(('STEP [%s] %s' % (cname, gname), h + '(assert (not %s))\n(check-sat)\n' % g))
        key = cname.split()[0][4:]
        if key not in EXPL:
            out.append(('STEP [%s] Inv(k+1)' % cname, h + '(assert (not %s))\n(check-sat)\n' % inv_next(gh)))
            continue
        # Inv(k+1) in three stages: (1) the iteration gives the explicit next state; (2) every conjunct of the active
        # part of Inv(k+1) for the explicit state; (3) assembly: (1) and (2) imply Inv(k+1) of the generic step
        U, V, R, Sx = EXPL[key]
        # names of the next state of the generic step (definitions): U7 = (u7 u v r s f), ...
        h = h + ''.join('(declare-const %s %s) (assert (= %s %s))\n' % (nm, ty, nm, S(df)) for nm, ty, df in
                        (('U7', 'Int', 'u7'), ('V7', 'Int', 'v7'), ('R7', 'Int', 'r7'), ('S7', 'Int', 's7'), ('F7', 'Bool', 'f1')))
        cv = '(and (= U7 %s) (= V7 %s) (= R7 %s) (= S7 %s) F7)' % (U, V, R, Sx)
        out.append(('STEP [%s] CASE-VALUE: next state = (%s, %s, %s, %s), f = 1' % (cname, U, V, R, Sx),
                    h + '(assert (not %s))\n(check-sat)\n' % cv))
        args = ('(+ k 1)', '(* 2 T)', 'true', U, V, R, Sx) + tuple(gh)
        parts = inst(INV_COMMON, args) + inst(INV_ACT, args)
        for j, c in enumerate(parts):
            out.append(('STEP [%s] Inv(k+1) part %02d: %s' % (cname, j, c), h + '(assert (not %s))\n(check-sat)\n' % c))
        # assembly as a standalone lemma (no hypothesis about the current state is needed): the explicit next state
        # satisfies every conjunct, the names U7.. equal it  ==>  Inv(k+1, U7, ..)
        dec = ''.join('(declare-const %s Int)\n' % v_ for v_ in ('k', 'T', 'u', 'v', 'r', 's', 'ga', 'gb', 'cR', 'cS', 'h',
                                                                'U7', 'V7', 'R7', 'S7')) + '(declare-const F7 Bool)\n'
        out.append(('STEP [%s] Inv(k+1) of the next state (U7, V7, R7, S7, F7): assembly from CASE-VALUE and the parts' % cname,
                    BASE + NB_LEMMA + INV + dec + '(assert %s)\n' % cv + ''.join('(assert %s)\n' % c for c in parts) +
                    '(assert (not (inv (+ k 1) (* 2 T) F7 U7 V7 R7 S7 %s)))\n(check-sat)\n' % ' '.join(gh)))
    out.append(('COVER: the cases cover every state', head(STEP_HYP) +
                '(assert (not (or %s)))\n(check-sat)\n' % ' '.join(c[1] for c in CASES)))
    return out


def controls():
    """must NOT be unsat: consistency of the case hypotheses, slightly too strong claims, no primality."""
    out = []
    # every case head (hypotheses, Inv, primality, auxiliary instances) is SATISFIED by a real state of a real run
    # (ground query, pow2 as a table): the STEP proofs are not vacuous. Expected answer: sat.
    from . import pq as qkal_pq
    pyc = {'ACT-B1': lambda f, x, u, v, r, p: f and x and v >= 1 and u % 2 == 0,
           'ACT-B2': lambda f, x, u, v, r, p: f and x and v >= 1 and u % 2 == 1 and v % 2 == 0,
           'ACT-B3': lambda f, x, u, v, r, p: f and x and v >= 1 and u % 2 == 1 and v % 2 == 1 and u > v,
           'ACT-B4a': lambda f, x, u, v, r, p: f and x and v >= 1 and u % 2 == 1 and v % 2 == 1 and u < v,
           'ACT-B4b lo': lambda f, x, u, v, r, p: f and x and v >= 1 and u == v and 2 * r < p,
           'ACT-B4b hi': lambda f, x, u, v, r, p: f and x and v >= 1 and u == v and 2 * r >= p,
           'DET lo': lambda f, x, u, v, r, p: f and x and v == 0 and 2 * r < p,
           'DET hi': lambda f, x, u, v, r, p: f and x and v == 0 and 2 * r >= p,
           'POST lo': lambda f, x, u, v, r, p: (not f) and x and 2 * r < p,
           'POST hi': lambda f, x, u, v, r, p: (not f) and x and 2 * r >= p,
           'ZERO': lambda f, x, u, v, r, p: x == 0}
    names = list(pyc)
    for (cname, hyp, gh, aux), key in zip(CASES, names):
        found = None
        for p in (7, 11, 13, 17, 19, 23, 29, 31, 61, 127, 251):
            N = p.bit_length() + 1
            for x in range(p):
                for (k, f, u, v, r, s, ga, gb, cR, cS) in qkal_pq.ghosts(N, p, x)[:-1]:
                    if found is None and pyc[key](f, x, u, v, r, p):
                        found = (p, N, x, k, f, u, v, r, s, ga, gb, cR, cS)
        if found is None:
            out.append(('CONTROL realizable [%s]: NO REAL STATE FOUND' % cname, '(set-logic ALL)\n(assert false)\n(check-sat)\n'))
            continue
        p, N, x, k, f, u, v, r, s, ga, gb, cR, cS = found
        g = _ground_base(p, N) + '(assert (= x %d))\n' % x
        st = ('(assert (and (= k %d) (= T %d) (= f %s) (= u %d) (= v %d) (= r %d) (= s %d) (= ga %d) (= gb %d) (= cR %d) (= cS %d)))\n'
              % (k, 1 << k, 'true' if f else 'false', u, v, r, s, ga, gb, cR, cS))
        out.append(('CONTROL realizable [%s] by p=%d N=%d x=%d k=%d' % (cname, p, N, x, k),
                    g + NB_LEMMA + INV + DEFS + STATE + STEP_HYP + PRIME + aux + '(assert %s)\n' % hyp + st + '(check-sat)\n'))
    strong = [('ACT-B1 u even', 'M-a with 4 r4 < p', '(=> %s (< (* 4 %s) p))' % (S('a3'), S('r4'))),
              ('ACT-B3 both odd, u > v', 'M-carry with 2^(N-1)', '(=> %s (< (+ %s %s) (pow2 (- N 1))))' % (S('c5'), S('s4'), S('r4'))),
              ('POST 2r >= p', 'M-dbl with r4 < p - 1', '(< %s (- p 1))' % S('r4')),
              ('ACT-B2 u odd, v even', 'Inv(k+1) with the wrong ghost cR', inv_next(('ga', '(- gb 1)', 'cR', 'cS'))),
              ('ACT-B4a both odd, u < v', 'v5 odd', '(odd %s)' % S('v5'))]
    cmap = {c[0]: c for c in CASES}
    for cname, nm, g in strong:
        hyp, aux = cmap[cname][1], cmap[cname][3]
        out.append(('CONTROL too strong [%s] %s' % (cname, nm),
                    head(STEP_HYP + PRIME + aux + '(assert %s)\n' % hyp) + '(assert (not %s))\n(check-sat)\n' % g))
    # the same too-strong claims REFUTED on real states of real runs (ground, pow2 as a table): expected answer sat
    # st = (k, f, u, v, r, s, ga, gb, cR, cS)
    tests = [('ACT-B1 u even', 'ACT-B1', lambda st, p, N: 4 * st[5] >= p, strong[0][2]),          # 4 s >= p
             ('ACT-B3 both odd, u > v', 'ACT-B3', lambda st, p, N: st[4] + st[5] >= 2 ** (N - 1), strong[1][2]),
             ('POST 2r >= p', 'POST hi', lambda st, p, N: st[4] == p - 1, strong[2][2]),
             ('ACT-B2 u odd, v even', 'ACT-B2', lambda st, p, N: st[8] != 0, strong[3][2]),
             ('ACT-B4a both odd, u < v', 'ACT-B4a', lambda st, p, N: True, strong[4][2])]
    for cname, key, extra_pred, g in tests:
        hyp, aux = cmap[cname][1], cmap[cname][3]
        found = None
        for p, N in itertools.product((7, 11, 13, 17, 19, 23, 29, 31, 61, 127, 251), (0, 1)):
            N += p.bit_length()
            for x in range(p):
                for st in qkal_pq.ghosts(N, p, x)[:-1]:
                    (k, f, u, v, r, s) = st[:6]
                    if found is None and pyc[key](f, x, u, v, r, p) and extra_pred(st, p, N):
                        found = (p, N, x) + st
        if found is None:
            out.append(('CONTROL counterexample [%s]: NO REAL STATE FOUND' % cname, '(set-logic ALL)\n(assert false)\n(check-sat)\n'))
            continue
        p, N, x, k, f, u, v, r, s, ga, gb, cR, cS = found
        st = ('(assert (and (= k %d) (= T %d) (= f %s) (= u %d) (= v %d) (= r %d) (= s %d) (= ga %d) (= gb %d) (= cR %d) (= cS %d)))\n'
              % (k, 1 << k, 'true' if f else 'false', u, v, r, s, ga, gb, cR, cS))
        out.append(('CONTROL realizable counterexample to the too strong claim [%s] p=%d N=%d x=%d k=%d' % (cname, p, N, x, k),
                    _ground_base(p, N) + '(assert (= x %d))\n' % x + NB_LEMMA + INV + DEFS + STATE + STEP_HYP + PRIME + aux +
                    '(assert %s)\n' % hyp + st + '(assert (not %s))\n(check-sat)\n' % g))
    c = cmap['ACT-B4b u = v, 2r < p']
    out.append(('CONTROL [%s] Inv(k+1) without the primality instance' % c[0],
                head(STEP_HYP + c[3] + '(assert %s)\n' % c[1]) + '(assert (not %s))\n(check-sat)\n' % inv_next(c[2])))
    fin = BASE + NB_LEMMA + INV + DEFS + STATE.replace('(assert (inv k T', '(assert (= k (* 2 N)))\n(assert (inv k T')
    out.append(('CONTROL FINAL with 2^(2N-1)', fin + '(assert (not (=> (>= x 1) (= (- (* x (- p r)) (pow2 (- (* 2 N) 1))) '
                '(* p (- x cR))))))\n(check-sat)\n'))
    # the ground header of the CONCRETE queries is consistent, and a wrong state is detected (expected: sat)
    out.append(('CONTROL CONCRETE header p=13 N=5 is consistent', _ground_base(13, 5) + '(check-sat)\n'))
    gs = qkal_pq.ghosts(5, 13, 7)
    (k, f, u, v, r, s, ga, gb, cR, cS) = gs[3]
    inv = INV.replace('(define-fun inv ((k Int)', '(define-fun invx ((x Int) (k Int)')
    body = _ground_base(13, 5).replace('(declare-const x Int)', '').replace('(assert (<= 0 x)) (assert (< x p))', '')
    out.append(('CONTROL CONCRETE a wrong state (u+2) is detected', body + inv +
                '(assert (not (invx 7 %d %d %s %d %d %d %d %d %d %d %d)))\n(check-sat)\n'
                % (k, 1 << k, 'true' if f else 'false', u + 2, v, r, s, ga, gb, cR, cS)))
    return out


def _ground_base(p, N):
    """BASE with pow2 as a table (no quantifiers) and N, nb, p fixed."""
    tab = '0'
    for i in range(2 * N + 2, -1, -1):
        tab = '(ite (= i %d) %d %s)' % (i, 1 << i, tab)
    base = BASE.replace('(declare-fun pow2 (Int) Int)', '(define-fun pow2 ((i Int)) Int %s)' % tab)
    base = '\n'.join(ln for ln in base.split('\n') if 'forall' not in ln and '(assert (= (pow2 0) 1))' not in ln) + '\n'
    return base + '(assert (= N %d)) (assert (= nb %d)) (assert (= p %d))\n' % (N, p.bit_length(), p)


def concrete(cases):
    """Inv and the step functions evaluated on REAL runs (qkal_pq), pow2 as a table: unsat iff all hold."""
    from . import pq as qkal_pq
    out = []
    for p, N in cases:
        nb = p.bit_length()
        tab = '0'
        for i in range(2 * N + 2, -1, -1):
            tab = '(ite (= i %d) %d %s)' % (i, 1 << i, tab)
        base = BASE.replace('(declare-fun pow2 (Int) Int)', '(define-fun pow2 ((i Int)) Int %s)' % tab)
        base = '\n'.join(ln for ln in base.split('\n') if 'forall' not in ln and '(assert (= (pow2 0) 1))' not in ln) + '\n'
        body = base.replace('(declare-const x Int)', '').replace('(assert (<= 0 x)) (assert (< x p))', '')
        inv = INV.replace('(define-fun inv ((k Int)', '(define-fun invx ((x Int) (k Int)')
        calls, steps = [], []
        for x in range(p):
            gs = qkal_pq.ghosts(N, p, x)
            for (k, f, u, v, r, s, ga, gb, cR, cS) in gs:
                b = 'true' if f else 'false'
                calls.append('(invx %d %d %d %s %d %d %d %d %d %d %d %d)' % (x, k, 1 << k, b, u, v, r, s, ga, gb, cR, cS))
                if k < 2 * N:
                    u7, v7, r7, s7, f1, m3, t1, _ = qkal_pq.step(N, p, u, v, r, s, f)
                    a = '%d %d %d %d %s' % (u, v, r, s, b)
                    steps.append('(= (u7 %s) %d) (= (v7 %s) %d) (= (r7 %s) %d) (= (s7 %s) %d) (= (f1 %s) %s) (= (m3 %s) %s) '
                                 '(= (t1 %s) %s)' % (a, u7, a, v7, a, r7, a, s7, a, 'true' if f1 else 'false',
                                                    a, 'true' if m3 else 'false', a, 'true' if t1 else 'false'))
        hdr = body + '(assert (= N %d)) (assert (= nb %d)) (assert (= p %d))\n' % (N, nb, p)
        out.append(('CONCRETE Inv on all runs p=%d N=%d' % (p, N), hdr + inv + '(assert (not (and %s)))\n(check-sat)\n' % '\n'.join(calls)))
        out.append(('CONCRETE step functions = qkal_pq.step, p=%d N=%d' % (p, N),
                    hdr + DEFS + '(assert (not (and %s)))\n(check-sat)\n' % '\n'.join(steps)))
    return out


def z3_version():
    """Version string of the Z3 that is used: the ``z3`` binary if there is one, else the z3-solver package."""
    try:
        return subprocess.run(['z3', '--version'], capture_output=True, text=True, timeout=30).stdout.strip() or 'z3 (no version)'
    except (OSError, subprocess.SubprocessError):
        try:
            import z3 as _z3
            return 'Z3 version %s (z3-solver)' % _z3.get_version_string()
        except ImportError:
            return None


def z3_run(path, to):
    """Run one SMT-LIB file: ``unsat`` / ``sat`` / ``unknown`` / ``timeout``; ``error`` if no Z3 is available."""
    try:
        r = subprocess.run(['z3', '-T:%d' % to, path], capture_output=True, text=True, timeout=to + 30)
        return (r.stdout.strip().splitlines() or ['?'])[0]
    except subprocess.TimeoutExpired:
        return 'timeout'
    except OSError:                       # no z3 binary: use the z3-solver package on the same SMT-LIB text
        try:
            import z3 as _z3
        except ImportError:
            return 'error'
        s = _z3.Solver()
        s.set('timeout', to * 1000)
        s.from_file(path)
        return str(s.check())


def run(timeout=120, dump=None, only=None, ctl_timeout=20, conc=((3, 2), (5, 3), (7, 3), (7, 4), (11, 4), (13, 4), (13, 5),
                                                                  (17, 5), (19, 6), (23, 5), (29, 5), (31, 5), (31, 7))):
    qs = queries() + concrete(conc)
    nq = len(qs)
    qs += controls()
    if only:
        qs = [q for q in qs if only in q[0]]
    d = dump or tempfile.mkdtemp(prefix='qkal_lemmas_')
    os.makedirs(d, exist_ok=True)
    rows = []
    for i, (name, txt) in enumerate(qs):
        path = os.path.join(d, 'q%03d.smt2' % i)
        with open(path, 'w') as f:
            f.write('; %s\n' % name + txt)
        ctl = name.startswith('CONTROL')
        t = time.time()
        ans = z3_run(path, ctl_timeout if ctl else timeout)
        exact_sat = name.startswith('CONTROL realizable') or name.startswith('CONTROL CONCRETE')
        ok = (ans == 'sat') if exact_sat else (ans != 'unsat') if ctl else (ans == 'unsat')
        rows.append(dict(i=i, name=name, answer=ans, ok=ok, sec=round(time.time() - t, 2),
                         sha256=hashlib.sha256(txt.encode()).hexdigest()))
        print('%-8s %-5s %6.2fs %s' % (ans, 'ok' if ok else 'FAIL', rows[-1]['sec'], name), flush=True)
    if dump is None:
        shutil.rmtree(d, ignore_errors=True)
    impl = [r for r in rows if not r['name'].startswith('CONTROL')]
    ctls = [r for r in rows if r['name'].startswith('CONTROL')]
    return dict(ok=all(r['ok'] for r in rows), implications=len(impl), proved=sum(r['ok'] for r in impl),
                controls=len(ctls), controls_ok=sum(r['ok'] for r in ctls), rows=rows)


if __name__ == '__main__':
    a = sys.argv[1:]
    to = int(a[a.index('--timeout') + 1]) if '--timeout' in a else 120
    dmp = a[a.index('--dump') + 1] if '--dump' in a else None
    onl = a[a.index('--only') + 1] if '--only' in a else None
    outp = a[a.index('--out') + 1] if '--out' in a else None
    t0 = time.time()
    res = run(to, dmp, onl)
    zv = z3_version()
    src = open(__file__, 'rb').read()
    res.update(z3=zv, source_sha256=hashlib.sha256(src).hexdigest(), sec=round(time.time() - t0, 1),
               finished=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    print('%d of %d implications unsat; controls not unsat: %d of %d; ok=%s; %.0fs' % (
        res['proved'], res['implications'], res['controls_ok'], res['controls'], res['ok'], res['sec']))
    if outp:
        with open(outp, 'w') as f:
            json.dump(res, f, indent=1)
    sys.exit(0 if res['ok'] else 1)


def ensure_lemmas(path=None, timeout=120, refresh=False, verbose=False):
    """The result of ``run`` for THIS source file and THIS Z3: read from the cache ($QOPT_CACHE/kaliski) if a file with the
    same source hash and Z3 version is there, else computed (a few minutes) and stored.  ``path`` reads/writes that file."""
    from .._cache import cache_dir
    zv = z3_version()
    sha = hashlib.sha256(open(__file__, 'rb').read()).hexdigest()
    path = path or os.path.join(str(cache_dir('kaliski')), 'lemmas_%s.json' % sha[:16])
    if not refresh and os.path.exists(path):
        try:
            res = json.load(open(path))
            if res.get('source_sha256') == sha and res.get('z3') == zv and zv:
                return res
        except (OSError, ValueError):
            pass
    if zv is None:
        raise RuntimeError('no Z3 (neither the z3 binary nor the z3-solver package): the lemmas cannot be checked')
    t0 = time.time()
    if verbose:
        res = run(timeout)
    else:
        import contextlib, io
        with contextlib.redirect_stdout(io.StringIO()):
            res = run(timeout)
    res.update(z3=zv, source_sha256=sha, sec=round(time.time() - t0, 1), finished=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    tmp = path + '.tmp%d' % os.getpid()
    with open(tmp, 'w') as f:
        json.dump(res, f, indent=1)
    os.replace(tmp, path)
    return res
