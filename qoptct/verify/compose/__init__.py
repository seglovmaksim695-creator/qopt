# SPDX-License-Identifier: Apache-2.0
"""qoptct.verify.compose -- compositional verification of large reversible arithmetic circuits on logical ANDs
(Gidney 2018) with measurement-based erasure; first client: the Kaliski modular inversion at full size (n = 256).

    r = verify_kaliski(qasm, p, reduction='generic', clean=False)      # qasm: path or text
    r['status']   'PROVED' | 'REFUTED' | 'NOT_PROVED'
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
at n = 256): the mathematical assumptions M1..M5 on Kaliski's algorithm (register bounds, parities) that appear as
preconditions, the fold lemma of the special reduction, and that P_std returns -x^{-1} 2^{2n} mod p with u = 1, v = 0.
Together: for every x in [1, p) the circuit maps |x>|0...> to the specified state for every measurement outcome,
with no residual phase.
"""
import random
import time

from .kal import Chain, ProofFailure, simulate_flat
from .lift import lift, Lifted, LiftError, check_gadget_lemmas

__all__ = ['verify_kaliski', 'pstd_test', 'kaliski_hints', 'lift', 'LiftError', 'check_gadget_lemmas', 'Chain',
           'ProofFailure', 'PRIMES']

PRIMES = {
    'secp256k1': 2 ** 256 - 2 ** 32 - 977,
    'p256': 2 ** 256 - 2 ** 224 + 2 ** 192 + 2 ** 96 - 1,
}


def kaliski_hints(p, reduction='generic', clean=False):
    """the (untrusted) certificate of the library circuit for these parameters."""
    from ...crypto.kaliski import kaliski_inverse
    return kaliski_inverse(p, reduction=reduction, clean=clean).hints()


def _expected(p, n, x):
    return (-pow(x, -1, p) * pow(2, 2 * n, p)) % p


def _sample_xs(p, count, seed=1):
    rnd = random.Random(seed)
    xs = [1, 2, 3, p - 1, p - 2, (p - 1) // 2, (p + 1) // 2, 1 << (p.bit_length() - 1), (1 << (p.bit_length() - 1)) - 1]
    xs = [x for x in dict.fromkeys(xs) if 1 <= x < p]
    while len(xs) < count:
        xs.append(rnd.randrange(1, p))
    return xs[:max(count, 1)] if p > count else list(range(1, p))


def _value(store, wires, j):
    return sum(((store.get(w, 0) >> j) & 1) << i for i, w in enumerate(wires))


def pstd_test(p, reduction='generic', xs=None, hints=None, count=256):
    """execute P_std (the contract chain itself) on the inputs xs and check every mathematical assumption, every
    carried fact and the final result (out = -x^{-1} 2^{2n} mod p, u = 1, v = 0). -> dict(ok, inputs, violations)."""
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
    return dict(ok=not viol, inputs=len(xs), violations=viol[:20])


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
    -> dict: status ('PROVED' | 'REFUTED' | 'NOT_PROVED'), stage, detail, counterexample, assumptions (name ->
       number of segments where the mathematical assumption is used), segments, time, ...
    See the module docstring for the exact meaning of PROVED."""
    t0 = time.time()
    n = p.bit_length()
    res = dict(status='NOT_PROVED', p=p, n=n, reduction=reduction, clean=bool(clean), stage='gadgets')

    def done(**kw):
        res.update(kw)
        res['time'] = round(time.time() - t0, 2)
        return res
    if not isinstance(p, int) or p < 5 or p % 2 == 0:
        return done(stage='input', detail='p must be an odd prime >= 5')
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
        fh = hints if not clean else dict(hints, marks=hints['marks'])
        t = pstd_test(p, reduction, hints=fh, count=test_samples)
        res['pstd_test'] = t
        if not t['ok']:
            return done(detail='a mathematical assumption of P_std is violated on a concrete input: %r'
                        % (t['violations'][:1],))
    return done(status='PROVED', stage='done', detail='machine-checked modulo the lemmas M1..M5 of compose/PROOF.md')
