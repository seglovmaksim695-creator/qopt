# SPDX-License-Identifier: Apache-2.0
"""compose.lemmas_z3 -- the mathematical assumptions M1..M6 of the Kaliski chain as verification conditions for Z3.

The contracts of compose.kal use properties of Kaliski's algorithm on INTEGERS as preconditions (M1..M6, PROOF.md
section 4). Here they are proved for EVERY n >= 2, every odd prime p with 2^(n-1) < p < 2^n and every x in [1, p), by
induction over the step number k with an explicit inductive invariant Inv; every implication is a quantifier-free
(up to the three axioms of 2^x) formula over the integers that the SMT solver Z3 refutes the negation of:

    INIT   Inv holds for the initial state (A, B, R, S) = (p, x, 0, 1), k = 0
    STEP   Inv(k) and k <= 2n-1  ==>  every assumption used by the contracts of step k (M1..M6), the facts that
           make the step exact (operated register even after the subtraction, no borrow), and Inv(k+1)
    FINAL  Inv(2n)  ==>  u = 1, v = 0, 0 <= r < 2^(2n), x*r + 2^(2n) = 0 (mod p)   (r = -x^{-1} 2^{2n} mod p)

The step is the frame step of P_std exactly as in compose.kal (w = wspec, d = w_prev xor w, swap, e = A_0,
A = (A - e*B)/2, R = R + e*S, S = 2*S); the 2 x 4 x 2 cases (frame bit, branch of Kaliski's step, v >= 1 or v = 0)
are separate queries.

What is assumed (axioms, all true statements):
  * pow2(0) = 1, pow2(i+1) = 2 pow2(i) for i >= 0, pow2 monotone on i >= 0   (pow2(i) stands for 2^i);
  * p odd, pow2(n-1) < p < pow2(n), and ONE instance of the primality of p: if B*(R+S) = p with 1 <= B < p then B = 1.
A register "fits w bits" is written value < pow2(w), with the widths wa, wb computed exactly as in the walker.

    python -m qoptct.verify.compose.lemmas_z3 [--timeout SEC] [--dump DIR]      (needs the z3 binary on PATH)
Exit code 0 iff every query is unsat (= every implication is proved).
"""
import itertools
import os
import shutil
import subprocess
import sys
import tempfile

__all__ = ['queries', 'run', 'z3_available']

PRELUDE = """(set-logic ALL)
(declare-fun pow2 (Int) Int)
(assert (= (pow2 0) 1))
(assert (forall ((i Int)) (! (=> (>= i 0) (= (pow2 (+ i 1)) (* 2 (pow2 i)))) :pattern ((pow2 i)))))
(assert (forall ((i Int) (j Int)) (! (=> (and (<= 0 i) (<= i j)) (<= (pow2 i) (pow2 j))) :pattern ((pow2 i) (pow2 j)))))
(declare-const n Int) (declare-const p Int) (declare-const x Int) (declare-const k Int)
(declare-const wp Bool)
(declare-const A Int) (declare-const B Int) (declare-const R Int) (declare-const S Int)
(declare-const a Int) (declare-const b Int) (declare-const cR Int) (declare-const cS Int)
(assert (>= n 2))
(assert (= (mod p 2) 1))
(assert (< (pow2 (- n 1)) p)) (assert (< p (pow2 n)))
(assert (<= 1 x)) (assert (< x p))
(define-fun odd ((t Int)) Bool (= (mod t 2) 1))
(define-fun min2 ((s Int) (t Int)) Int (ite (<= s t) s t))
(define-fun max2 ((s Int) (t Int)) Int (ite (<= s t) t s))
(define-fun inv ((k Int) (wp Bool) (A Int) (B Int) (R Int) (S Int) (a Int) (b Int) (cR Int) (cS Int)) Bool
  (let ((u (ite wp B A)) (v (ite wp A B)) (r (ite wp S R)) (s (ite wp R S)) (sg (ite wp (- 1) 1)))
    (and (<= 0 k) (<= k (* 2 n)) (>= A 0) (>= B 0) (>= R 0) (>= S 0)
         (= (+ (* A S) (* B R)) p)
         (>= u 1) (< r (pow2 k)) (<= s (pow2 k))
         (= (+ (* x R) (* sg A (pow2 k))) (* p cR))
         (= (- (* x S) (* sg B (pow2 k))) (* p cS))
         (=> (>= k 1) (odd B))
         (=> (= k 0) (not wp))
         (or (and (>= v 1) (< v p) (= (+ a b) (- (* 2 n) k)) (<= 0 a) (<= a n) (<= 0 b) (<= b n)
                  (< A (pow2 a)) (< B (pow2 b)))
             (and (= v 0) (= u 1) wp (>= k 1))))))
"""

# the frame step of P_std (compose.kal) for one case. The case fixes the frame bit wp, the branch of Kaliski's step
# and whether v >= 1; then w, d = wp xor w and e are constants, which the query CASE proves from the generic
# definitions; the other queries of the case use the constants (so that every formula is free of if-then-else).
STEP = """(assert (inv k wp A B R S a b cR cS))
(assert (<= k (- (* 2 n) 1)))
(define-fun u () Int (ite wp B A))
(define-fun v () Int (ite wp A B))
(define-fun wa () Int (min2 n (max2 1 (- (* 2 n) 1 k))))
(define-fun wb () Int (min2 wa (max2 1 (div (- (* 2 n) k) 2))))
; one instance of the primality of p (used when u = v)
(assert (=> (and (= (* B (+ R S)) p) (<= 1 B) (< B p)) (= B 1)))
"""
W_GEN = '(and (odd u) (or (not (odd v)) (<= u v)))'
CASES_BR = [('u even', '(assert (not (odd u)))', 0, 0),                 # name, hypothesis, w, e
            ('u odd, v even', '(assert (and (odd u) (not (odd v))))', 1, 0),
            ('both odd, u<=v', '(assert (and (odd u) (odd v) (<= u v)))', 1, 1),
            ('both odd, u>v', '(assert (and (odd u) (odd v) (> u v)))', 0, 1)]
CASES_PH = [('v>=1', '(assert (>= v 1))'), ('v=0', '(assert (= v 0))')]

# standalone lemmas of integer arithmetic (congruence of multiplication), proved by Z3 without any hypothesis and
# then used as instances (Z = S1 and Z = pow2(k)):
LEMMAS = [('LEMMA: 2*A2 = Y -> 2*A2*Z = Y*Z',
           '(declare-const A2 Int) (declare-const Y Int) (declare-const Z Int)\n(assert (= (* 2 A2) Y))\n'
           '(assert (not (= (* 2 A2 Z) (* Y Z))))\n(check-sat)\n'),
          ('LEMMA: 2*A2 = A1 - B1 -> 2*A2*Z = A1*Z - B1*Z',
           '(declare-const A2 Int) (declare-const A1 Int) (declare-const B1 Int) (declare-const Z Int)\n'
           '(assert (= (* 2 A2) (- A1 B1)))\n(assert (not (= (* 2 A2 Z) (- (* A1 Z) (* B1 Z)))))\n(check-sat)\n'),
          ('LEMMA: A = B -> A*Z = B*Z',
           '(declare-const A Int) (declare-const B Int) (declare-const Z Int)\n(assert (= A B))\n'
           '(assert (not (= (* A Z) (* B Z))))\n(check-sat)\n')]


def _case(wp, br, ph):
    """-> (case name, head text, list of (goal name, extra text, goal formula))."""
    nb, cb, w, e = CASES_BR[br]
    nph, cph = CASES_PH[ph]
    d = wp ^ w
    A1, B1, R1, S1, a1, b1, cR1, cS1 = ('B', 'A', 'S', 'R', 'b', 'a', 'cS', 'cR') if d else \
        ('A', 'B', 'R', 'S', 'a', 'b', 'cR', 'cS')
    T, F = 'true', 'false'
    dd, ww = (T if d else F), (T if w else F)
    sg = '(- 1)' if w else '1'
    name = 'wp=%d; %s; %s' % (wp, nb, nph)
    head = PRELUDE + STEP + ('(assert wp)\n' if wp else '(assert (not wp))\n') + cb + '\n' + cph + '\n'
    sub = '(- %s %s)' % (A1, B1) if e else A1
    R2 = '(+ %s %s)' % (R1, S1) if e else R1
    S2 = '(* 2 %s)' % S1
    cR2 = '(+ %s %s)' % (cR1, cS1) if e else cR1
    cS2 = '(* 2 %s)' % cS1
    T1, T2 = '(pow2 k)', '(* 2 (pow2 k))'           # pow2(k+1) = 2 pow2(k): the axiom, k >= 0
    post = '(declare-const A2 Int)\n(assert (= (* 2 A2) %s))\n' % sub
    lem = ''
    for Z in (S1, T1):
        lem += '(assert (= (* 2 A2 %s) %s))\n' % (
            Z, '(- (* %s %s) (* %s %s))' % (A1, Z, B1, Z) if e else '(* %s %s)' % (A1, Z))
    tie = '(assert (=> (= A B) (and (= (* A S) (* B S)) (= (* A R) (* B R)))))\n'      # third lemma, Z = S, R
    a2 = '(- %s 1)' % a1
    u2, v2 = (B1, 'A2') if w else ('A2', B1)
    r2, s2 = (S2, R2) if w else (R2, S2)
    goals = [
        ('CASE: w, d, e are the constants of the case', '',
         '(and (= %s %s) (= (xor wp %s) %s) (= (mod %s 2) %d))' % (W_GEN, ww, W_GEN, dd, A1, e)),
        ('M1 width: A, B < 2^wa', '', '(and (< A (pow2 wa)) (< B (pow2 wa)))'),
        ('M2 (P3): k>=1 and d -> A, B odd', '', '(=> (and (>= k 1) %s) (and (odd A) (odd B)))' % dd),
        ('M4-rs: d -> R, S < 2^n', '', '(=> %s (and (< S (pow2 n)) (< R (pow2 n))))' % dd),
        ('M3 short v: e -> B1 < 2^wb', '', '(< %s (pow2 wb))' % B1 if e else 'true'),
        ('exact step: e -> B1 odd and B1 <= A1; A1 - e*B1 even and >= 0', '',
         '(and %s (= (mod %s 2) 0) (>= %s 0))' % ('(odd %s) (<= %s %s)' % (B1, B1, A1) if e else 'true', sub, sub)),
        ('M4-sfit: e -> S1 < 2^n', '', '(< %s (pow2 n))' % S1 if e else 'true'),
        ('M5: k = 2n-2 -> w', '', '(=> (= k (- (* 2 n) 2)) %s)' % ww),
        ('M5: k = 2n-1 -> w, A even, d = 0, e = 0', '',
         '(=> (= k (- (* 2 n) 1)) (and %s (not (odd A)) (not %s) %s))' % (ww, dd, F if e else T)),
        ('M4-rfit: R + e*S < 2^n', post + lem, '(< %s (pow2 n))' % R2),
        ('M6: 2*S < 2^(2n)', post, '(< %s (pow2 (* 2 n)))' % S2),
        ('Inv(k+1): range, signs', post,
         '(and (<= (+ k 1) (* 2 n)) (>= A2 0) (>= %s 0) (>= %s 0) (>= %s 0))' % (B1, R2, S2)),
        ('Inv(k+1): A*S + B*R = p', post + lem, '(= (+ (* A2 %s) (* %s %s)) p)' % (S2, B1, R2)),
        ('Inv(k+1): u >= 1', post, '(>= %s 1)' % u2),
        ('Inv(k+1): r < 2^(k+1), s <= 2^(k+1)', post,
         '(and (< %s (pow2 (+ k 1))) (<= %s (pow2 (+ k 1))))' % (r2, s2)),
        ('Inv(k+1): x*R + sg*A*2^(k+1) = p*cR', post + lem,
         '(= (+ (* x %s) (* %s A2 %s)) (* p %s))' % (R2, sg, T2, cR2)),
        ('Inv(k+1): x*S - sg*B*2^(k+1) = p*cS', post,
         '(= (- (* x %s) (* %s %s %s)) (* p %s))' % (S2, sg, B1, T2, cS2)),
        ('Inv(k+1): B odd', post, '(odd %s)' % B1),
        ('Inv(k+1): bounds while v >= 1, or (u, v) = (1, 0) in the frame of v', post + lem + tie,
         '(or (and (>= %s 1) (< %s p) (= (+ %s %s) (- (* 2 n) (+ k 1))) (<= 0 %s) (<= %s n) (<= 0 %s) (<= %s n) '
         '(< A2 (pow2 %s)) (< %s (pow2 %s))) (and (= %s 0) (= %s 1) %s))'
         % (v2, v2, a2, b1, a2, a2, b1, b1, a2, B1, b1, v2, u2, ww)),
    ]
    return name, head, goals


def queries():
    """-> list of (name, smt2 text); every text must be unsat."""
    out = list(LEMMAS)
    out.append(('INIT: Inv(0) for (A, B, R, S) = (p, x, 0, 1)',
                PRELUDE + '(assert (not (inv 0 false p x 0 1 n n 1 0)))\n(check-sat)\n'))
    out.append(('FINAL: Inv(2n) -> u = 1, v = 0, 0 <= r < 2^(2n), x*r + 2^(2n) = p*c',
                PRELUDE + '(assert (inv (* 2 n) wp A B R S a b cR cS))\n'
                '(assert (not (and wp (= A 0) (= B 1) (>= S 0) (< S (pow2 (* 2 n))) '
                '(= (+ (* x S) (pow2 (* 2 n))) (* p cS)))))\n(check-sat)\n'))
    for wp, br, ph in itertools.product((0, 1), range(4), range(2)):
        name, head, goals = _case(wp, br, ph)
        for gname, extra, g in goals:
            if g == 'true':
                continue
            out.append(('STEP [%s] %s' % (name, gname), head + extra + '(assert (not %s))\n(check-sat)\n' % g))
    return out


def controls():
    """negative controls: queries that must NOT be unsat (the hypotheses of a case alone, and goals that are
    slightly too strong). They guard against a vacuous proof (contradictory hypotheses)."""
    out = []
    for wp, br, ph in ((0, 3, 0), (1, 2, 0), (1, 1, 1), (0, 0, 0)):
        name, head, goals = _case(wp, br, ph)
        out.append(('CONTROL [%s] the hypotheses of the case are consistent' % name, head + '(check-sat)\n'))
    name, head, goals = _case(0, 3, 0)
    g = dict((a, (b, c)) for a, b, c in goals)
    strong = [('M1 with wa - 1', '', '(and (< A (pow2 (- wa 1))) (< B (pow2 (- wa 1))))'),
              ('M3 with wb - 1', '', '(< B (pow2 (- wb 1)))'),
              ('M4-sfit with n - 2', '', '(< S (pow2 (- n 2)))'),
              ('r < 2^k instead of 2^(k+1)',) + (g['Inv(k+1): r < 2^(k+1), s <= 2^(k+1)'][0], '(< (+ R S) (pow2 k))'),
              ('wrong sign in the congruence',) + (g['Inv(k+1): x*R + sg*A*2^(k+1) = p*cR'][0],
                                                    '(= (- (* x (+ R S)) (* A2 (* 2 (pow2 k)))) (* p (+ cR cS)))')]
    for nm, extra, goal in strong:
        out.append(('CONTROL [%s] too strong: %s' % (name, nm), head + extra + '(assert (not %s))\n(check-sat)\n' % goal))
    name, head, goals = _case(0, 2, 0)                    # the tie u = v without the primality instance
    g = dict((a, (b, c)) for a, b, c in goals)
    extra, goal = g['Inv(k+1): bounds while v >= 1, or (u, v) = (1, 0) in the frame of v']
    h2 = head.replace('(assert (=> (and (= (* B (+ R S)) p) (<= 1 B) (< B p)) (= B 1)))\n', '')
    assert h2 != head
    out.append(('CONTROL [%s] without the primality of p' % name, h2 + extra + '(assert (not %s))\n(check-sat)\n' % goal))
    return out


def kaliski_trace(p, x):
    """the states (k, wp, A, B, R, S, a, b, cR, cS) of P_std before every step k = 0..2n, with the ghost values
    of the invariant (a, b: exponents of the bounds, frozen once v = 0; cR, cS: quotients of the congruences)."""
    n = p.bit_length()
    k, wp, A, B, R, S, a, b, cR, cS = 0, 0, p, x, 0, 1, n, n, 1, 0
    out = [(k, wp, A, B, R, S, a, b, cR, cS)]
    for k in range(2 * n):
        u, v = (B, A) if wp else (A, B)
        w = int(u % 2 == 1 and (v % 2 == 0 or u <= v))
        if wp ^ w:
            A, B, R, S, a, b, cR, cS = B, A, S, R, b, a, cS, cR
        e = A % 2
        A, R, S, a, cR, cS = (A - e * B) // 2, R + e * S, 2 * S, a - 1, cR + e * cS, 2 * cS
        wp = w
        out.append((k + 1, wp, A, B, R, S, a, b, cR, cS))
    return out


def concrete_queries(primes):
    """the invariant Inv of the SMT text evaluated by Z3 on real runs (all x in [1, p) for the given small
    primes), with pow2 given as a table instead of the axioms: one ground query per prime, unsat iff Inv holds in
    every state of every run. Ties the SMT text of Inv to the algorithm (and shows it is not contradictory)."""
    out = []
    for p in primes:
        n = p.bit_length()
        tab = '0'
        for i in range(2 * n, -1, -1):
            tab = '(ite (= i %d) %d %s)' % (i, 1 << i, tab)
        pre = PRELUDE.replace('(declare-fun pow2 (Int) Int)', '(define-fun pow2 ((i Int)) Int %s)' % tab)
        pre = '\n'.join(ln for ln in pre.split('\n') if 'forall' not in ln) + '\n'
        # inv refers to the global x: make it a parameter by a local definition
        body = pre.replace('(declare-const x Int)', '').replace('(assert (<= 1 x)) (assert (< x p))', '')
        body = body.replace('(define-fun inv ((k Int)', '(define-fun invx ((x Int) (k Int)')
        calls = []
        for x in range(1, p):
            for st in kaliski_trace(p, x):
                vals = [('true' if t else 'false') if i == 1 else ('(- %d)' % -t if t < 0 else str(t))
                        for i, t in enumerate(st)]
                calls.append('(invx %d %s)' % (x, ' '.join(vals)))
        out.append(('CONCRETE: Inv holds along all runs for p = %d' % p,
                    body + '(assert (= n %d)) (assert (= p %d))\n(assert (not (and %s)))\n(check-sat)\n'
                    % (n, p, '\n'.join(calls))))
    return out


def z3_available():
    return shutil.which('z3') is not None


def run(timeout=60, dump=None, verbose=False, only=None, control_timeout=10):
    """run every query with the z3 binary. -> dict(ok, total, proved, failed=[(name, answer)], controls_ok).
    ok: every implication is unsat AND no negative control is unsat."""
    qs = queries() + concrete_queries((5, 7, 11, 13, 17, 19, 23, 29, 31, 61))
    nq = len(qs)
    qs += controls()
    if only:
        qs = [q for q in qs if only in q[0]]
        nq = len([q for q in qs if not q[0].startswith('CONTROL')])
    failed = []
    bad_controls = []
    d = dump or tempfile.mkdtemp(prefix='kaliski_lemmas_')
    os.makedirs(d, exist_ok=True)
    for i, (name, txt) in enumerate(qs):
        path = os.path.join(d, 'q%03d.smt2' % i)
        with open(path, 'w') as f:
            f.write('; %s\n' % name + txt)
        ctl = name.startswith('CONTROL')
        to = control_timeout if ctl else timeout
        try:
            r = subprocess.run(['z3', '-T:%d' % to, path], capture_output=True, text=True, timeout=to + 30)
            ans = (r.stdout.strip().splitlines() or ['?'])[0]
        except subprocess.TimeoutExpired:
            ans = 'timeout'
        if ctl:
            if ans == 'unsat':
                bad_controls.append((name, ans))
        elif ans != 'unsat':
            failed.append((name, ans))
        if verbose:
            print('%-8s %s' % (ans, name), flush=True)
    if dump is None:
        shutil.rmtree(d, ignore_errors=True)
    return dict(ok=not failed and not bad_controls, total=nq, proved=nq - len(failed), failed=failed,
                controls=len(qs) - nq, bad_controls=bad_controls)


if __name__ == '__main__':
    a = sys.argv[1:]
    to = int(a[a.index('--timeout') + 1]) if '--timeout' in a else 60
    dmp = a[a.index('--dump') + 1] if '--dump' in a else None
    onl = a[a.index('--only') + 1] if '--only' in a else None
    if not z3_available():
        print('z3 binary not found')
        sys.exit(3)
    res = run(to, dmp, verbose='--quiet' not in a, only=onl)
    print('%d of %d implications proved by Z3; %d negative controls, %d of them wrongly unsat'
          % (res['proved'], res['total'], res['controls'], len(res['bad_controls'])))
    for name, ans in res['failed']:
        print('  NOT PROVED (%s): %s' % (ans, name))
    sys.exit(0 if res['ok'] else 1)
