# SPDX-License-Identifier: Apache-2.0
"""qoptct.verify.compose -- compositional verification of large reversible arithmetic circuits on logical ANDs
(Gidney 2018) with measurement-based erasure; first client: the Kaliski modular inversion at full size (n = 256).

    r = verify_kaliski(qasm, p, reduction='generic', clean=False)      # qasm: path or text
    r['status']   'PROVED' | 'REFUTED' | 'UNKNOWN'
    python -m qoptct.verify.compose gen    --prime secp256k1 --reduction special --clean -o DIR
    python -m qoptct.verify.compose verify --prime secp256k1 --reduction special --clean DIR/....qasm

WHAT 'PROVED' MEANS (details and the paper proofs: compose/PROOF.md).
Machine-checked, for the given OpenQASM text, at full width, for all inputs at once:
  L0  the text is a sequence of x / cx / AND gadget / ERASE gadget (strict parser, compose.lift), and the two gadget
      lemmas hold (exact arithmetic);
  L1  the lifted circuit is cut into segments (certificate marks, untrusted); for EVERY segment and ALL values of its
      entry wires that satisfy the precondition, the real ops of the segment keep every AND target at 0, satisfy
      every erase condition, and implement the word-level contract of the corresponding step of P_std (ROBDD proof
      at the real width, up to 256 bits; every ancilla returns to 0); in clean mode the mirror segment of the inverse
      half is proved to undo it;
  L2  the bookkeeping of the chain (registers, wires known to be zero, facts carried between segments, the copy of
      the result, coverage of all ops) is checked by the walker.
Proved on paper (PROOF.md) and TESTED by executing P_std (exhaustively for all primes p < 2^11 and all x; sampled
at n = 256): the mathematical assumptions M1..M6 on Kaliski's algorithm (register bounds, parities) that appear as
preconditions, the fold lemma of the special reduction, and that P_std returns -x^{-1} 2^{2n} mod p with u = 1, v = 0.
Together: for every x in [1, p) the circuit maps |x>|0...> to the specified state for every measurement outcome,
with no residual phase.
"""
import random
import time

from .kal import Chain
from .walker import Walker, ProofFailure, simulate_flat, register_fact
from .lift import lift, Lifted, LiftError, check_gadget_lemmas

__all__ = ['verify_kaliski', 'pstd_test', 'kaliski_hints', 'lift', 'LiftError', 'check_gadget_lemmas', 'Chain',
           'ProofFailure', 'PRIMES', 'named_prime', 'Walker', 'register_fact']

PRIMES = {
    'secp256k1': 2 ** 256 - 2 ** 32 - 977,
    'p256': 2 ** 256 - 2 ** 224 + 2 ** 192 + 2 ** 96 - 1,
    'curve25519': 2 ** 255 - 19,
    'p384': 2 ** 384 - 2 ** 128 - 2 ** 96 + 2 ** 32 - 1,
    'p521': 2 ** 521 - 1,
}


def _is_probable_prime(m, rounds=64, seed=12345):
    if m < 2 or m % 2 == 0:
        return m == 2
    d, s = m - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    rnd = random.Random(seed)
    for _ in range(rounds):
        a = rnd.randrange(2, m - 1)
        x = pow(a, d, m)
        if x in (1, m - 1):
            continue
        for _ in range(s - 1):
            x = x * x % m
            if x == m - 1:
                break
        else:
            return False
    return True


_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)      # deterministic for every m < 3.18e23 > 2^64


def _mr_witness(m, a, d, s):
    """True if the base a proves m composite."""
    x = pow(a, d, m)
    if x in (1, m - 1):
        return False
    for _ in range(s - 1):
        x = x * x % m
        if x == m - 1:
            return False
    return True


def primality(m, rounds=40):
    """-> (verdict, how): verdict 'prime' (certain: a named prime, or m < 2^64 by deterministic Miller-Rabin),
    'probable' (m >= 2^64: the 12 fixed bases and `rounds` further pseudo-random bases found no witness; a composite
    passes with probability < 4^-rounds on average, never claimed as proof) or 'composite'."""
    if not isinstance(m, int) or m < 2:
        return 'composite', 'm < 2'
    for q in _MR_BASES:
        if m == q:
            return 'prime', 'small prime'
        if m % q == 0:
            return 'composite', 'divisible by %d' % q
    for name, v in PRIMES.items():
        if m == v:
            return 'prime', 'known prime (%s)' % name
    if m < 2 ** 64:
        d, s = m - 1, 0
        while d % 2 == 0:
            d //= 2
            s += 1
        if any(_mr_witness(m, a, d, s) for a in _MR_BASES):
            return 'composite', 'Miller-Rabin witness'
        return 'prime', 'deterministic Miller-Rabin (bases 2..37, exact for m < 2^64)'
    d, s = m - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    rnd = random.Random(m & 0xffffffff)
    bases = list(_MR_BASES) + [rnd.randrange(2, m - 1) for _ in range(rounds)]
    if any(_mr_witness(m, a, d, s) for a in bases):
        return 'composite', 'Miller-Rabin witness'
    return 'probable', 'probable prime: Miller-Rabin, %d bases (a probabilistic test, not a proof of primality)' % len(bases)


def named_prime(name):
    """a prime by name: secp256k1, p256, curve25519, p384, p521, or rand256 = a 256-bit prime without special form
    (the first probable prime, Miller-Rabin with 64 bases, after floor(2^254 * pi); it shows that nothing depends on
    the shape of p)."""
    if name == 'rand256':
        if 'rand256' not in _CACHE:
            pi80 = 31415926535897932384626433832795028841971693993751058209749445923078164062862089
            m = ((1 << 254) * pi80 // 10 ** 79) | 1
            while not _is_probable_prime(m):
                m += 2
            _CACHE['rand256'] = m
        return _CACHE['rand256']
    return PRIMES[name]


_CACHE = {}


def kaliski_hints(p, reduction='generic', clean=False):
    """the (untrusted) certificate of the library circuit for these parameters."""
    from ...crypto.kaliski import kaliski_inverse
    return kaliski_inverse(p, reduction=reduction, clean=clean).hints()


def _expected(p, n, x):
    try:
        return (-pow(x, -1, p) * pow(2, 2 * n, p)) % p
    except ValueError:                                   # x not invertible mod p: outside the domain
        return None


def _sample_xs(p, count, seed=1):
    rnd = random.Random(seed)
    xs = [1, 2, 3, p - 1, p - 2, (p - 1) // 2, (p + 1) // 2, 1 << (p.bit_length() - 1), (1 << (p.bit_length() - 1)) - 1]
    if p <= count:
        return list(range(1, p))
    if count >= 1024 and pow(2, p - 1, p) == 1:
        # structured inputs, where the bounds M1, M3, M4 are tight: powers of two and their neighbours, their
        # negatives and modular inverses, small numbers and their inverses
        n = p.bit_length()
        try:
            for k in range(n):
                t = 1 << k
                xs += [t, t - 1, t + 1, p - t, p - t - 1, p - t + 1, pow(t, -1, p), p - pow(t, -1, p), (t * 3) % p]
            for c in range(2, 40):
                xs += [c, p - c, pow(c, -1, p), p - pow(c, -1, p)]
        except ValueError:                               # p is not prime
            pass
    xs = [x for x in dict.fromkeys(xs) if 1 <= x < p]
    while len(xs) < count:
        xs.append(rnd.randrange(1, p))
    return xs[:max(count, 1)] if p > count else list(range(1, p))


def _value(store, wires, j):
    return sum(((store.get(w, 0) >> j) & 1) << i for i, w in enumerate(wires))


def pstd_test(p, reduction='generic', xs=None, hints=None, count=256, lifted=None, clean=False):
    """execute P_std (the contract chain itself) on the inputs xs and check every mathematical assumption, every
    carried fact and the final result (out = -x^{-1} 2^{2n} mod p, u = 1, v = 0). -> dict(ok, inputs, violations).
    With `lifted` (the lifted circuit the hints belong to): also run the circuit itself classically on the same
    inputs and compare EVERY wire with the state of P_std (forward circuit), resp. with |x>|K>|0> (clean circuit):
    an end-to-end consistency check of the chain against the real circuit (key 'circuit_agrees')."""
    if hints is None:
        hints = kaliski_hints(p, reduction)
    if xs is None:
        xs = _sample_xs(p, count)
    xs = list(xs)
    ch = Chain(p, reduction, hints, xs=xs).run()
    viol = list(ch.viol)
    n = p.bit_length()
    for j, x in enumerate(xs):
        got = (_value(ch.store, ch.out, j), _value(ch.store, ch.u, j), _value(ch.store, ch.v, j))
        if got != (_expected(p, n, x), 1, 0) and len(viol) < 20:
            viol.append(dict(segment='final', what='result (out, u, v) = %r' % (got,), x=x))
    res = dict(ok=not viol, inputs=len(xs), violations=viol[:20])
    if lifted is not None:
        S = len(xs)
        init = {i: sum(((x >> i) & 1) << j for j, x in enumerate(xs)) for i in range(n)}
        v, bad = simulate_flat(lifted, init, S)
        agree = bad is None
        if agree:
            if clean:
                want = dict(init)
                for i in range(n):
                    want[n + i] = sum(((_expected(p, n, x) or 0) >> i & 1) << j for j, x in enumerate(xs))
            else:
                want = ch.store
            agree = all(v[w] == want.get(w, 0) for w in range(lifted.nq))
        res['circuit_agrees'] = agree
        res['ok'] = res['ok'] and agree
    return res


def _refute(L, p, n, clean, out, count, seed=7):
    """search for a concrete counterexample by a classical run of the lifted circuit."""
    xs = _sample_xs(p, count, seed)
    S = len(xs)
    init = {}
    for i in range(n):
        m = 0
        for j, x in enumerate(xs):
            if (x >> i) & 1:
                m |= 1 << j
        init[i] = m
    v, bad = simulate_flat(L, init, S)
    if bad is not None:
        j = (bad[2] & -bad[2]).bit_length() - 1
        return dict(x=xs[j], reason='op %d: %s' % (bad[0], bad[1]))
    if out is None or any((not isinstance(w, int)) or w < 0 or w >= L.nq for w in out):
        return None
    full = (1 << S) - 1
    for j, x in enumerate(xs):
        got = sum(((v[w] >> j) & 1) << i for i, w in enumerate(out))
        if got != _expected(p, n, x):
            return dict(x=x, reason='wrong result on the output wires')
    if clean:
        keep = set(range(n)) | set(out)
        for w in range(L.nq):
            if w not in keep and v[w]:
                j = (v[w] & -v[w]).bit_length() - 1
                return dict(x=xs[j], reason='wire %d is not returned to zero' % w)
        for i in range(n):
            if v[i] != init[i]:
                j = ((v[i] ^ init[i]) & -(v[i] ^ init[i])).bit_length() - 1
                return dict(x=xs[j], reason='the input register is changed')
    del full
    return None


def verify_kaliski(qasm, p, reduction='generic', clean=False, hints=None, lifted=None, refute_samples=64,
                   test_samples=128, progress=None):
    """Compositional proof that the OpenQASM circuit `qasm` (path or text) is the Kaliski inversion for the odd
    prime p (n = p.bit_length(); wires 0..n-1 hold x in [1, p), all other wires |0>):
      clean=False: the wires r['out'] end with -x^{-1} 2^{2n} mod p, r['u'] with 1, r['v'] with 0 (other wires:
                   garbage), every AND target was |0>, every measurement-based erasure is exact (no phase);
      clean=True:  |x>|0> -> |x>|-x^{-1} 2^{2n} mod p> on wires 0..n-1, n..2n-1, every other wire back to |0>.
    hints: the certificate (KaliskiCircuit.hints(), untrusted); default: regenerated from the library.
    PRIMALITY of p is checked first (key 'p_primality'): a composite p -> status 'ERROR' (never PROVED); a named
    prime or p < 2^64 is certain; a larger unnamed p gets PROVED only 'under the condition that p is prime'
    (Miller-Rabin is probabilistic) and the detail says so.
    -> dict: status ('PROVED' | 'REFUTED' | 'UNKNOWN' | 'ERROR'), stage, detail, counterexample, assumptions (name ->
       number of segments where the mathematical assumption is used), segments, time, ...
    See the module docstring for the exact meaning of PROVED."""
    t0 = time.time()
    n = p.bit_length()
    res = dict(status='UNKNOWN', p=p, n=n, reduction=reduction, clean=bool(clean), stage='gadgets')

    def done(**kw):
        res.update(kw)
        res['time'] = round(time.time() - t0, 2)
        return res
    if not isinstance(p, int) or isinstance(p, bool) or p < 5 or p % 2 == 0:
        return done(stage='input', detail='p must be an odd prime >= 5')
    pv, phow = primality(p)
    res['p_primality'] = '%s: %s' % (pv, phow)
    if pv == 'composite':
        return done(status='ERROR', stage='input', detail='p is not prime (%s): the theorem is for a prime p; the '
                    'inverse does not exist for x sharing a factor with p' % phow)
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
    if hints is None:
        hints = kaliski_hints(p, reduction, clean)
    out_claim = hints.get('fwd_out' if clean else 'out')
    res['stage'] = 'chain'
    try:
        if not isinstance(hints.get('marks'), list):
            raise ProofFailure('certificate: no marks')
        ch = Chain(p, reduction, hints, lifted=L, clean=clean, progress=progress).run()
        if not clean and list(ch.out) != list(hints.get('out', [])):
            raise ProofFailure('interface: the result is on wires other than the declared output wires')
    except ProofFailure as e:
        det = {k: v for k, v in e.detail.items() if k != 'entry_counterexample'}
        cex = _refute(L, p, n, clean, list(range(n, 2 * n)) if clean else out_claim, refute_samples) \
            if refute_samples else None
        if cex is not None:
            return done(status='REFUTED', detail=str(e), where=det, counterexample=cex)
        return done(detail=str(e), where=det, local_counterexample=e.detail.get('entry_counterexample'))
    res.update(segments=ch.nseg, assumptions=dict(sorted(ch.assumed.items())), lemma_facts=len(ch.admitted),
               bdd_nodes_total=ch.nodes, bdd_nodes_max=ch.maxnodes,
               out=list(range(n, 2 * n)) if clean else list(ch.out), u=list(ch.u), v=list(ch.v))
    if test_samples:
        res['stage'] = 'pstd_test'
        t = pstd_test(p, reduction, hints=hints, count=test_samples, lifted=L, clean=clean)
        res['pstd_test'] = t
        if not t['ok']:
            return done(detail='P_std disagrees with the circuit on a concrete input' if not t.get('circuit_agrees', 1)
                        else 'a mathematical assumption of P_std is violated on a concrete input: %r'
                        % (t['violations'][:1],))
    note = '' if pv == 'prime' else ' PROVED under the condition that p is prime (only a probabilistic primality test was run).'
    return done(status='PROVED', stage='done', proved_modulo='paper lemmas M1..M6 on the integer algorithm (and the fold '
                'lemma for reduction=special), compose/PROOF.md sections 4-5; tested by pstd_test',
                detail='every segment of the circuit is machine-proved against its contract at full width.' + note)


from .ecadd import verify_ecadd, wec_test, ecadd_hints, ECChain, lift_circuit  # noqa: E402

__all__ += ['verify_ecadd', 'wec_test', 'ecadd_hints', 'ECChain', 'lift_circuit']
