# SPDX-License-Identifier: Apache-2.0
"""Public API of qoptct.dyadic: exact equivalence of circuits with dyadic angles pi*p/2^k.

    prove_equal(circuit, reference, phase='global', ancillas=None, budget_sec=120.0)
    prove_spec(circuit, spec, phase='global', ancillas=None, budget_sec=120.0)
    prove_qasm(qasm, reference, ...)          # OpenQASM 2/3 text or path; reference: QASM, .qc path or DCircuit
        -> dict(verdict='PROVED' | 'REFUTED' | 'UNKNOWN' | 'ERROR' | 'INTERNAL', method, global_phase, ...)

Meaning of PROVED.  prove_equal: circuit |x, 0> = exp(i g) (reference (x) I) |x, 0> for every basis input x of the
non-ancilla wires (ancillas start in 0; extra wires beyond the reference must end in 0); g is reported exactly as
global_phase (a multiple of pi); phase='exact' also demands g = 0.  prove_spec: the same with the specification in
place of the reference (QFTSpec, ArithSpec, PermutationSpec; see qoptct.dyadic.spec).
How.  The path sum of the composition with the adjoint of the reference/specification is reduced with the exact
rules of qoptct.dyadic.core.  PROVED only from (a) a form with no path variables left whose wires are the
expected polynomials and whose phase is a constant (both canonical), with the amplitude factor exactly 1, or (b) an
exact enumeration of the remaining small sum in Z[zeta] (zeta = exp(2 pi i / 2^K), canonical power basis
1, zeta, ..., zeta^(2^(K-1)-1)).  REFUTED from a canonical violation of that form (with a basis-state witness) or from
the exact enumeration.  Anything else (rules stuck on a large sum, time or size limit, unsupported gate, angle not
dyadic) is UNKNOWN.  Never raises: bad input is ERROR, a self-inconsistency INTERNAL (reported, never PROVED).
ArithSpec answers carry the assumption 'fourier-shift lemma' (proved in docs/DYADIC.md).  With selfcheck=True the
independent dense simulator re-checks every answer up to MAX_SELFCHECK qubits; a disagreement gives INTERNAL.
"""
import cmath
import math
import time
from fractions import Fraction

from .core import DPathSum, Blowup, Timeout, ZERO, peval
from .circuit import DCircuit, NotDyadic, Unsupported, invert_prims, needed_K, parse_qasm, from_qc, _mod2
from .spec import QFTSpec, ArithSpec, ClassicalSpec, anf_over_vars

__all__ = ['prove_dyadic', 'prove_equal', 'prove_spec', 'prove_qasm', 'VERDICTS', 'ENUM_MAX_VARS', 'MAX_SELFCHECK',
           'LEMMA']

VERDICTS = ('PROVED', 'REFUTED', 'UNKNOWN', 'ERROR', 'INTERNAL')
ENUM_MAX_VARS = 16       # exact enumeration of a stuck sum: inputs + path variables
MAX_SELFCHECK = 10       # dense re-check of every answer up to this many qubits
LEMMA = ('fourier-shift lemma: the exact Fourier path sum of an ArithSpec equals the classical map '
         '|x,s> -> |x,(s+G(x)) mod 2^m> (geometric sum; docs/DYADIC.md, Lemma 2)')


def _coef(ps, r):
    v = Fraction(ps.N, 2) * _mod2(r)
    if v.denominator != 1:
        raise AssertionError('K=%d too small for the phase pi*%s' % (ps.K, r))
    return int(v) & ps.M


def _run(ps, prims, step):
    for p in prims:
        op = p[0]
        if op == 'h':
            ps.h(p[1])
        elif op == 'xor':
            ps.xor(p[1], p[2])
        elif op == 'swap':
            ps.swap(p[1], p[2])
        elif op == 'cswap':
            ps.cswap(p[1], p[2], p[3])
        elif op == 'ph':
            ps.phase(p[1], _coef(ps, p[2]))
        elif op == 'phx':
            ps.phase_xor(p[1], _coef(ps, p[2]))
        else:
            raise AssertionError('unknown primitive %r' % (op,))
        if step:
            ps.reduce()


def _new(K, n, anc, deadline, max_terms, order):
    ps = DPathSum(K, max_terms=max_terms, deadline=deadline)
    ps.order = order
    xv = {}
    for q in range(n):
        if q in anc:
            ps.add_wire(ZERO)
        else:
            v = ps.var('x', 'q%d' % q)
            xv[q] = v
            ps.add_wire(frozenset([1 << v]))
    return ps, xv


def _violations(ps, want):
    viol = []
    for q, W in enumerate(want):
        if ps.w[q] != W:
            viol.append(('wire', q, ps.w[q] ^ W))
    nonconst = {mo: c for mo, c in ps.ph.items() if mo}
    if nonconst:
        viol.append(('phase', None, nonconst))
    return viol


def _witness(v, xv):
    """Basis input (dict qubit -> bit) exhibiting a canonical violation (only input variables are left).
    Wire: a minimal-degree monomial m of the nonzero polynomial P; at the indicator of m, P = 1.
    Phase: a minimal-degree nonconstant monomial m with c_m != 0; F(ind m) - F(0) = c_m != 0."""
    kind, q, P = v
    monos = list(P.keys()) if isinstance(P, dict) else list(P)
    m = min(monos, key=lambda mo: (bin(mo).count('1'), mo))
    return {qq: (m >> var) & 1 for qq, var in xv.items()}


def _gphase(ps):
    return _mod2(Fraction(2 * ps.ph.get(0, 0), ps.N))


def _fmt(r):
    if r is None:
        return None
    r = _mod2(r)
    return '0' if not r else 'pi*%d/%d' % (r.numerator, r.denominator)


def _enumerate(ps, xv, want_out, limit):
    """Exact operator of a small path sum.  want_out(xbits) -> expected output basis state (int over the wires).
    -> ('PROVED', phase or None) | ('REFUTED', witness) | None (too large)."""
    xs = sorted(xv.items(), key=lambda t: t[1])
    ys = sorted(ps.ys)
    if len(xs) + len(ys) > limit:
        return None
    half = ps.N >> 1
    diag0 = None
    for xa in range(1 << len(xs)):
        ps.check_limits()
        xval = 0
        xbits = {}
        for i, (q, v) in enumerate(xs):
            b = (xa >> i) & 1
            xbits[q] = b
            if b:
                xval |= 1 << v
        acc = {}
        for ya in range(1 << len(ys)):
            val = xval
            for i, v in enumerate(ys):
                if (ya >> i) & 1:
                    val |= 1 << v
            out, f = ps.evaluate(val)
            d = acc.setdefault(out, {})
            sgn = 1
            if f >= half:          # zeta^(k + 2^(K-1)) = -zeta^k: canonical exponents in [0, 2^(K-1))
                f -= half
                sgn = -1
            d[f] = d.get(f, 0) + sgn
        exp_out = want_out(xbits)
        diag = None
        for out, d in acc.items():
            d = {k: c for k, c in d.items() if c}
            if out != exp_out:
                if d:
                    return 'REFUTED', xbits
                continue
            diag = d
        if not diag:
            return 'REFUTED', xbits
        if diag0 is None:
            diag0 = diag
        elif diag != diag0:
            return 'REFUTED', xbits
    # (circuit o adjoint) = lambda I with lambda = sqrt(2)^e2 * diag0; read lambda as a root of unity when it is one
    ph = None
    if len(diag0) == 1:
        (k, c), = diag0.items()
        if Fraction(c) ** 2 * Fraction(2) ** ps.e2 == 1:
            ph = _mod2(Fraction(2 * k, ps.N) + (1 if c < 0 else 0))
    return 'PROVED', ph


def _canon(d, N):
    """Canonical element of Z[zeta] (zeta = exp(2 pi i / N)): exponents in [0, N/2), zero coefficients dropped."""
    half = N >> 1
    out = {}
    for k, c in d.items():
        k %= N
        if k >= half:
            k -= half
            c = -c
        out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c}


def _times_sqrt2(d, N):
    """d * sqrt(2), with sqrt(2) = zeta^(N/8) - zeta^(3N/8)."""
    e1, e3 = N >> 3, 3 * (N >> 3)
    out = {}
    for k, c in d.items():
        out[k + e1] = out.get(k + e1, 0) + c
        out[k + e3] = out.get(k + e3, 0) - c
    return _canon(out, N)


def _lam_equal(a, b, N):
    """Exact equality of sqrt(2)^e1 * d1 and sqrt(2)^e2 * d2."""
    (d1, e1), (d2, e2) = a, b
    if e1 < e2:
        (d1, e1), (d2, e2) = (d2, e2), (d1, e1)
    while e1 - e2 >= 2:
        d1 = {k: 2 * c for k, c in d1.items()}
        e1 -= 2
    if e1 - e2 == 1:
        d1 = _times_sqrt2(d1, N)
    return _canon(d1, N) == _canon(d2, N)


WITNESS_TRIES = 24
WITNESS_ENUM = 18


def _inputs(n, anc, tries, seed):
    """Candidate basis inputs (dict qubit -> bit; ancillas 0): all-0, all-1, seeded random ones interleaved with the
    alternating patterns and single bits."""
    import random
    free = [q for q in range(n) if q not in anc]
    rng = random.Random(seed)
    rnd = [{q: rng.randint(0, 1) for q in free} for _ in range(max(tries, 4))]
    singles = [{q: int(j == i) for j, q in enumerate(free)} for i in rng.sample(range(len(free)), min(4, len(free)))]
    cand = [{q: 0 for q in free}, {q: 1 for q in free}, rnd[0], {q: i & 1 for i, q in enumerate(free)}, rnd[1],
            {q: 1 - (i & 1) for i, q in enumerate(free)}, rnd[2]] + singles + rnd[3:]
    seen, out = set(), []
    for x in cand:
        key = tuple(x[q] for q in free)
        if key not in seen:
            seen.add(key)
            out.append(x)
    return out[:tries]


def _conj_mul(d, N):
    """d * conj(d) in canonical form (conj(zeta^k) = zeta^-k)."""
    out = {}
    for k1, c1 in d.items():
        for k2, c2 in d.items():
            k = k1 - k2
            out[k] = out.get(k, 0) + c1 * c2
    return _canon(out, N)


def _unit_modulus(d, e2, N):
    """|sqrt(2)^e2 * d|^2 == 1, exactly."""
    nn = _conj_mul(d, N)
    if e2 <= 0:
        return nn == {0: 1 << (-e2)}
    return {k: c << e2 for k, c in nn.items()} == {0: 1}


def _scalar_sum(ps):
    """Exact sum over the (few) live path variables of exp(2 pi i F / N) -> canonical element of Z[zeta]."""
    ys = sorted(ps.ys)
    acc = {}
    for ya in range(1 << len(ys)):
        val = 0
        for i, v in enumerate(ys):
            if (ya >> i) & 1:
                val |= 1 << v
        f = 0
        for mo, c in ps.ph.items():
            if mo & val == mo:
                f += c
        f &= ps.M
        acc[f] = acc.get(f, 0) + 1
    return _canon(acc, ps.N)


def _fix_outputs(ps, exp_out):
    """Restrict the sum to the paths whose output is the basis state exp_out (then <exp_out|psi> is a scalar sum).
    A wire W_q = y xor R with y a path variable occurring only there and only linearly: y := R xor t_q (one term of
    the sum over y survives, e unchanged).  -> 'ok', 'zero' (a constant wire differs: amplitude 0) or None."""
    from .core import bits
    one = frozenset([0])
    for q in range(len(ps.w)):
        t = one if (exp_out >> q) & 1 else ZERO
        W = ps.w[q]
        if W == t:
            continue
        if all(mo == 0 for mo in W):
            return 'zero'
        cnt = {}
        for mo in W:
            for v in bits(mo & ps.ymask):
                cnt[v] = cnt.get(v, 0) + 1
        cand = [v for v, k in cnt.items() if k == 1 and (1 << v) in W and ps.wcnt.get(v, 0) == 1 and v in ps.ys]
        if not cand:
            return None
        y = min(cand)
        ps.subst(y, W ^ frozenset([1 << y]) ^ t)
        ps.ys.discard(y)
        if ps.w[q] != t:
            raise AssertionError('output fixing failed')
    return 'ok'


def _witness_search(circ, pre, post, n, anc, K, deadline, max_terms, expected, tries=WITNESS_TRIES, seed=1):
    """Exact per-input refutation (never used to prove anything).  For a concrete basis input x the composition is a
    unit vector psi(x) given by a path sum over path variables only.  It must equal lambda |expected(x)> with one
    lambda (|lambda| = 1) for all x.  The amplitude A = <expected(x)|psi(x)> is computed exactly: the output wires
    are fixed to expected(x), the rules reduce the (now internal) variables, a remaining sum of at most WITNESS_ENUM
    variables is enumerated in Z[zeta].  |A| != 1, or A(x) != A(x') -> witness.  -> (input, reason) or None."""
    prims = circ.prims()
    lam0 = None
    for x in _inputs(n, anc, tries, seed):
        if deadline is not None and time.time() > deadline:
            return None
        try:
            ps = DPathSum(K, max_terms=max_terms, deadline=deadline)
            for q in range(n):
                ps.add_wire(frozenset([0]) if x.get(q) else ZERO)
            if pre is not None:
                pre(ps)
            ps.reduce()
            _run(ps, prims, True)
            if post:
                _run(ps, post, True)
            ps.reduce()
            if ps.zero:
                continue
            exp_out = expected(x)
            fx = _fix_outputs(ps, exp_out)
            if fx == 'zero':
                return x, 'the expected output basis state has amplitude 0'
            if fx is None:
                continue
            ps.unfreeze()
            ps.reduce()
            if ps.zero:
                return x, 'the expected output basis state has amplitude 0'
            if any(mo for mo in ps.ph) and not ps.ys:
                raise AssertionError('phase depends on a variable that is not live')
            if not ps.ys:
                lam = (_canon({ps.ph.get(0, 0): 1}, ps.N), ps.e2)
            elif len(ps.ys) <= WITNESS_ENUM:
                lam = (_scalar_sum(ps), ps.e2)
            else:
                continue
        except (Timeout, Blowup):
            continue
        if not lam[0] or not _unit_modulus(lam[0], lam[1], ps.N):
            return x, '|<expected|output>| != 1: the output is not the expected basis state'
        if lam0 is None:
            lam0 = (lam, x)
        elif not _lam_equal(lam, lam0[0], ps.N):
            return x, 'phase differs from the phase at input %s' % ([lam0[1].get(q, 0) for q in range(n)],)
    return None


def _widen(c, n):
    w = DCircuit(n)
    w.ops = list(c.ops)
    w.gphase = c.gphase
    return w


def _dense_check(res, circ, ref=None, spec=None, anc=()):
    """Independent dense re-check (small n); a disagreement makes the answer INTERNAL."""
    try:
        from . import sim
        import numpy as np
        n = circ.n
        if n > MAX_SELFCHECK:
            res['dense_check'] = 'skipped (n > %d)' % MAX_SELFCHECK
            return
        U = sim.unitary(circ)
        V = sim.unitary(ref if ref.n == n else _widen(ref, n)) if ref is not None else sim.spec_unitary(spec, n)
        keep = [b for b in range(1 << n) if not any((b >> a) & 1 for a in anc)]
        Uk, Vk = U[:, keep], V[:, keep]
        lam = complex(np.vdot(Vk[:, 0], Uk[:, 0]))
        ok = abs(abs(lam) - 1) < 1e-7 and float(np.max(np.abs(Uk - lam * Vk))) < 1e-7
        v = res['verdict']
        if v == 'PROVED':
            good = ok
            if good and res.get('global_phase_pi') is not None:
                good = abs(cmath.exp(1j * math.pi * float(res['global_phase_pi'])) - lam) < 1e-7
        elif v == 'REFUTED':
            good = (not ok) or (res.get('phase_mode') == 'exact' and abs(lam - 1) > 1e-7)
        else:
            res['dense_check'] = 'equal' if ok else 'different'
            return
        res['dense_check'] = 'agrees' if good else 'DISAGREES'
        if not good:
            res['internal_verdict'] = v
            res['verdict'] = 'INTERNAL'
            res['reason'] = 'the dense simulator disagrees with the path-sum verdict'
    except Exception as e:      # the re-check is auxiliary
        res['dense_check'] = 'error: %s' % (repr(e)[:160],)


def _finish(res, t0):
    res['time'] = round(time.time() - t0, 3)
    return res


def _prove(circ, pre, post, want_fn, K, n, anc, phase, budget_sec, max_terms, strategies, lemma, extra_check,
           expected=None):
    """pre(ps): applies what comes before the circuit (the adjoint of a reference or a specification) and returns the
    path variables to freeze in strategy 'freeze'; post: primitive list applied after the circuit (or None).
    expected(x) -> output basis state for a basis input x (dict qubit -> bit), used by the witness search.
    Order: strategies until one leaves no path variable; after the first failure an exact per-input witness search
    (refutation only); then an exact enumeration of a small stuck sum; otherwise UNKNOWN."""
    t0 = time.time()
    deadline = None if budget_sec is None else t0 + budget_sec
    res = dict(verdict='UNKNOWN', method=None, K=K, n=n, phase_mode=phase, strategies_tried=[])
    prims = circ.prims()
    if expected is None:
        def expected(x):
            return sum(b << q for q, b in x.items())
    state = dict(searched=False)

    def witness():
        state['searched'] = True
        left = None if deadline is None else deadline - time.time()
        if left is not None and left <= 0:
            return False
        wit = _witness_search(circ, pre, post, n, anc, K, None if left is None else time.time() + 0.5 * left,
                              max_terms, expected)
        res['strategies_tried'].append(('witness', 'found' if wit else 'none'))
        if wit is None:
            return False
        res.update(verdict='REFUTED', method='witness', reason=wit[1],
                   counterexample=dict(input=[wit[0].get(q, 0) for q in range(n)], kind='witness'))
        return True
    best = None
    # quick exact refutation first (4 structured/random inputs, at most 10% of the budget): a wrong circuit is
    # usually caught here in a fraction of a second, before a symbolic strategy can diverge
    qdl = None if deadline is None else t0 + 0.1 * (deadline - t0)
    wit = _witness_search(circ, pre, post, n, anc, K, qdl, max_terms, expected, tries=7)
    if wit is not None:
        res.update(verdict='REFUTED', method='witness', reason=wit[1],
                   counterexample=dict(input=[wit[0].get(q, 0) for q in range(n)], kind='witness'))
        res['strategies_tried'].append(('quick-witness', 'found'))
        strategies = []
    for i, strat in enumerate(strategies):
        if i > 0 and not state['searched'] and witness():
            break
        left = None if deadline is None else deadline - time.time()
        if left is not None and left <= 0:
            break
        last = i == len(strategies) - 1 and state['searched']
        share = None if left is None else time.time() + (left if last else 0.6 * left)
        step = strat in ('step', 'freeze', 'step-max')
        try:
            ps, xv = _new(K, n, anc, share, max_terms, 'max' if strat == 'step-max' else 'min')
            fz = pre(ps) if pre is not None else set()
            if strat == 'freeze':
                ps.frozen |= set(fz)
            if step:
                ps.reduce()
            _run(ps, prims, step)
            if post:
                _run(ps, post, step)
            ps.unfreeze()
            ps.reduce()
        except (Timeout, Blowup) as e:
            res['strategies_tried'].append((strat, 'limit: %s' % e))
            continue
        res['strategies_tried'].append((strat, 'ys=%d' % len(ps.ys)))
        if ps.zero:
            res.update(verdict='INTERNAL', reason='the path sum is identically zero (impossible for a unitary)')
            return _finish(res, t0)
        if best is None or len(ps.ys) < len(best[0].ys):
            best = (ps, xv, strat)
        if not ps.ys:
            break
    if res['verdict'] == 'REFUTED':
        pass
    elif best is not None and not best[0].ys:
        ps, xv, strat = best
        res.update(strategy=strat, stats=dict(ps.stats), size=ps.size(), e2=ps.e2, method='normal-form')
        if ps.e2 != 0:
            res.update(verdict='INTERNAL',
                       reason='no path variables left but the amplitude factor is sqrt(2)^%d' % ps.e2)
            return _finish(res, t0)
        viol = _violations(ps, want_fn(ps, xv))
        if viol:
            w = _witness(viol[0], xv)
            res.update(verdict='REFUTED', reason=[(v[0], v[1]) for v in viol[:6]],
                       counterexample=dict(input=[w.get(q, 0) for q in range(n)], kind=viol[0][0],
                                           wire=viol[0][1]))
            if viol[0][0] == 'phase':
                res['counterexample']['compare_with'] = [0] * n
        else:
            g = _gphase(ps)
            res.update(global_phase_pi=g, global_phase=_fmt(g))
            if phase == 'exact' and g:
                res.update(verdict='REFUTED', reason='equal only up to the global phase exp(i %s)' % _fmt(g))
            else:
                res['verdict'] = 'PROVED'
    else:
        r = None
        if best is not None:
            ps, xv, strat = best
            res.update(strategy=strat, stats=dict(ps.stats), size=ps.size(), e2=ps.e2)
            want = want_fn(ps, xv)

            def want_out(xbits):
                val = 0
                for q, b in xbits.items():
                    if b:
                        val |= 1 << xv[q]
                return sum(1 << q for q, W in enumerate(want) if peval(W, val))
            try:
                ps.deadline = deadline
                r = _enumerate(ps, xv, want_out, ENUM_MAX_VARS)
            except (Timeout, Blowup) as e:
                res['enum_error'] = str(e)
        if r is None and not state['searched'] and witness():
            pass
        elif r is None:
            res.update(verdict='UNKNOWN', method='rules-stuck' if best is not None else 'limit',
                       ys_left=None if best is None else len(best[0].ys),
                       reason='all strategies hit a limit' if best is None else
                       '%d path variables left after the rules; no exact enumeration or witness' % len(best[0].ys))
        elif r[0] == 'REFUTED':
            res.update(verdict='REFUTED', method='enumeration',
                       counterexample=dict(input=[r[1].get(q, 0) for q in range(n)], kind='amplitude'))
        else:
            res.update(method='enumeration', global_phase_pi=r[1], global_phase=_fmt(r[1]))
            if phase == 'exact' and r[1] != 0:
                if r[1] is None:
                    res.update(verdict='UNKNOWN', reason='global phase is not read as a root of unity')
                else:
                    res.update(verdict='REFUTED', reason='equal only up to the global phase exp(i %s)' % _fmt(r[1]))
            else:
                res['verdict'] = 'PROVED'
    if lemma and res['verdict'] == 'PROVED':
        res['assumptions'] = [lemma]
    if circ.angle_rec.get('n_float'):
        res['float_angles'] = dict(circ.angle_rec)
    if extra_check is not None:
        extra_check(res)
    return _finish(res, t0)


def _safe(fn):
    def wrap(*a, **k):
        t0 = time.time()
        try:
            return fn(*a, **k)
        except (NotDyadic, Unsupported) as e:
            return dict(verdict='UNKNOWN', method='input', reason=str(e), time=round(time.time() - t0, 3))
        except (ValueError, TypeError, KeyError, IndexError, OSError) as e:
            return dict(verdict='ERROR', method='input', reason='%s: %s' % (type(e).__name__, e),
                        time=round(time.time() - t0, 3))
        except (Timeout, Blowup) as e:
            return dict(verdict='UNKNOWN', method='limit', reason=str(e), time=round(time.time() - t0, 3))
        except Exception as e:      # never PROVED on an unexpected error
            return dict(verdict='INTERNAL', method='error', reason='%s: %s' % (type(e).__name__, e),
                        time=round(time.time() - t0, 3))
    wrap.__name__ = fn.__name__
    wrap.__doc__ = fn.__doc__
    return wrap


def _identity_want(n):
    def want_fn(ps, xv):
        return [frozenset([1 << xv[q]]) if q in xv else ZERO for q in range(n)]
    return want_fn


@_safe
def prove_equal(circuit, reference, phase='global', ancillas=None, budget_sec=120.0, max_terms=400000,
                strategies=('step', 'end', 'step-max'), selfcheck=True):
    """circuit |x, 0> = exp(i g) (reference (x) I) |x, 0> for all x?  The reference acts on wires
    0..reference.n-1; ancillas (start in 0) default to the extra wires reference.n..circuit.n-1, which must also end
    in 0.  Ancillas inside the reference (the qoptct clean model) only start in 0; the reference decides their end.
    If no ancilla is inside the reference, its adjoint is applied before the circuit (circuit o ref^dagger = I),
    otherwise after it (ref^dagger o circuit = I on the inputs with ancillas 0)."""
    if phase not in ('global', 'exact'):
        raise ValueError("phase must be 'global' or 'exact'")
    if not isinstance(circuit, DCircuit) or not isinstance(reference, DCircuit):
        raise TypeError('circuit and reference must be DCircuit objects (see parse_qasm, from_qc)')
    n = circuit.n
    if reference.n > n:
        raise ValueError('the reference has more qubits than the circuit')
    anc = set(range(reference.n, n)) if ancillas is None else {int(a) for a in ancillas}
    if any(a < 0 or a >= n for a in anc):
        raise ValueError('ancilla wire out of range')
    if not set(range(reference.n, n)) <= anc:
        raise ValueError('every wire beyond the reference must be an ancilla')
    rp = invert_prims(reference.prims())
    K = max(needed_K(circuit.prims()), needed_K(rp))
    inside = any(a < reference.n for a in anc)
    if inside:
        pre, post = None, rp
    else:
        def pre(ps):
            _run(ps, rp, False)
            return set()
        post = None
    ext = (lambda res: _dense_check(res, circuit, ref=reference, anc=anc)) if selfcheck else None
    return _prove(circuit, pre, post, _identity_want(n), K, n, anc, phase, budget_sec, max_terms,
                  list(strategies), None, ext)


@_safe
def prove_spec(circuit, spec, phase='global', ancillas=None, budget_sec=120.0, max_terms=400000,
               strategies=None, selfcheck=True):
    """circuit |x, 0> = exp(i g) (spec (x) I) |x, 0>?  spec: QFTSpec, ArithSpec or PermutationSpec; ancillas start in 0
    and must end in 0 (the specification must not use them).  ArithSpec answers carry the 'fourier-shift lemma'."""
    if phase not in ('global', 'exact'):
        raise ValueError("phase must be 'global' or 'exact'")
    if not isinstance(circuit, DCircuit):
        raise TypeError('circuit must be a DCircuit')
    n = circuit.n
    anc = set() if ancillas is None else {int(a) for a in ancillas}
    if any(a < 0 or a >= n for a in anc):
        raise ValueError('ancilla wire out of range')
    spec.check(n)          # parameter domain before any computation (AUD7)
    if isinstance(spec, ClassicalSpec):
        if anc:
            raise ValueError('ClassicalSpec: give the zero wires in the specification, not as ancillas')
        K = max(needed_K(circuit.prims()), 3)
        zero = set(spec.zero_wires)

        def want_fn(ps, xv):
            return [anf_over_vars(spec.anf(q), xv) for q in range(n)]
        ext = (lambda res: _dense_check(res, circuit, spec=spec, anc=zero)) if selfcheck else None

        def expected(x):
            return spec.table[spec.index(sum(b << q for q, b in x.items()))]
        return _prove(circuit, None, None, want_fn, K, n, zero, phase, budget_sec, max_terms,
                      list(strategies or ('step', 'end', 'step-max')), None, ext, expected)
    if isinstance(spec, QFTSpec):
        touched = set(spec.wires)
    elif isinstance(spec, ArithSpec):
        touched = set(spec.target) | {q for _, qs in spec.terms for q in qs}
    else:
        raise TypeError('unknown specification type')
    if touched & anc:
        raise ValueError('the specification uses an ancilla wire')
    K = max(needed_K(circuit.prims()), spec.needed_K(), 3)
    lemma = LEMMA if isinstance(spec, ArithSpec) else None
    ext = (lambda res: _dense_check(res, circuit, spec=spec, anc=anc)) if selfcheck else None
    if strategies is None:
        # without freezing, the rules expand the carries of the classical map (exponential in the target width)
        strategies = (('freeze', 'step') if len(spec.target) <= 12 else ('freeze',)) if isinstance(spec, ArithSpec) \
            else ('step', 'end')
    strategies = list(strategies)
    return _prove(circuit, spec.apply_adjoint, None, _identity_want(n), K, n, anc, phase, budget_sec, max_terms,
                  strategies, lemma, ext)


def prove_qasm(qasm, reference, phase='global', clean=False, budget_sec=120.0, **kw):
    """OpenQASM 2/3 circuit (text or path) against a reference: DCircuit, OpenQASM text/path, or a qoptct .qc path.
    clean=True with a .qc reference: wires not on its .i line start in 0 (qoptct clean model); the extra wires of the
    circuit are ancillas in both models."""
    t0 = time.time()
    try:
        c = parse_qasm(qasm) if isinstance(qasm, str) else qasm
        anc = None
        if isinstance(reference, DCircuit):
            ref = reference
        elif isinstance(reference, str) and reference.endswith('.qc') and '\n' not in reference:
            ref, names, ins = from_qc(reference)
            if clean and ins is not None:
                anc = {i for i, nm in enumerate(names) if nm not in ins} | set(range(len(names), c.n))
        else:
            ref = parse_qasm(reference)
    except (NotDyadic, Unsupported) as e:
        return dict(verdict='UNKNOWN', method='input', reason=str(e), time=round(time.time() - t0, 3))
    except (ValueError, TypeError, KeyError, IndexError, OSError) as e:
        return dict(verdict='ERROR', method='input', reason='%s: %s' % (type(e).__name__, e),
                    time=round(time.time() - t0, 3))
    return prove_equal(c, ref, phase=phase, ancillas=anc, budget_sec=budget_sec, **kw)


def prove_dyadic(circuit, target, **kw):
    """One entry point for the bridges: target is a reference DCircuit (-> prove_equal) or a specification
    (QFTSpec, ArithSpec / adder_spec / const_adder_spec / multiplier_spec, ClassicalSpec / PermutationSpec;
    -> prove_spec).  Keyword arguments as there (phase, ancillas, budget_sec, max_terms, strategies, selfcheck).
    -> dict with verdict PROVED | REFUTED | UNKNOWN | ERROR | INTERNAL; never raises."""
    if isinstance(target, DCircuit):
        return prove_equal(circuit, target, **kw)
    return prove_spec(circuit, target, **kw)
