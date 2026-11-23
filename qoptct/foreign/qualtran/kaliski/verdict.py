"""The verdict for one pair (n, p): PROVED only if EVERY item of the criterion holds.

``prove_kaliski(n, p)`` runs the whole pipeline (Walker chain proof, execution mode, independent simulation, Z3 lemmas, primality
of p, canaries) and returns the verdict dictionary; ``judge(evidence, ...)`` makes the verdict from ready evidence (this is what
the negative tests call with damaged evidence).

What PROVED means (see docs/QUALTRAN_KALISKI.md): the Qualtran circuit ``KaliskiModInverse(bitsize=n, mod=p)`` lowered gate by
gate computes, for every input 0 <= x < p, x^-1 * 2^(2n) mod p (0 for x = 0) with the junk register of the documented closed form,
every other wire clean, under the assumptions that the Z3 lemmas close for all n, p, x (p prime).

Items (all must hold; the names are the keys of ``checks``):
  chain, assumptions_covered   the Walker chain proof ended CHAIN_PROVED with only the assumptions the lemmas cover
  run                          the contracts executed on the sample inputs agree with Qualtran's classical action and with P_Q
  flow                         the wire-flow check ran in both modes with the number of checks the formula demands
  client_hash                  the evidence was made by exactly this client and this lemma file
  sim                          independent bit-parallel simulation of the lowered circuit on the sample inputs
  sources_unchanged            no source file changed during the run
  lemmas                       the Z3 file is for THIS lemma source, THIS list of queries, all proved, controls as expected
  prime                        a certificate for p (or a named standard prime)
  canaries                     (n = 256 or when given) damaged circuits are all NOT PROVED, enough of them, no false PROVED
  samples                      enough sample inputs
"""
import hashlib
import inspect
import json
import os
import time

FLOW_VERSION = 'flow-1'
ASSUMPTIONS = {'M-dbl', 'M-a', 'M-s', 'M-carry', 'M-end', 'M-end-vf'}
HERE = os.path.dirname(os.path.abspath(__file__))


def _sha_file(name):
    with open(os.path.join(HERE, name), 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def client_hashes():
    return dict(chain_sha256=_sha_file('chain.py'), lemmas_sha256=_sha_file('lemmas_z3.py'), flow_version=FLOW_VERSION)


def flow_expected(n):
    """The number of wire-flow checks (written here separately, so that the verdict does not trust the client's own number)."""
    return dict(iterations=2 * n, eq_checks=31 * 2 * n - 1 + 6, new_checks=4 * 2 * n + 1 + 3, disjoint_checks=6 * 2 * n)


def flow_ok(fl, n):
    return bool(fl) and fl.get('version') == FLOW_VERSION and fl.get('ok') is True and fl.get('junk_ok') is True \
        and all(fl.get(k) == v for k, v in flow_expected(n).items())


_QCACHE = []


def expected_queries():
    """The full list of lemma queries of the CURRENT lemma source: (name, text) of the 243 implications, then the controls."""
    if not _QCACHE:
        from . import lemmas_z3 as L
        conc = inspect.signature(L.run).parameters['conc'].default
        qs = L.queries() + L.concrete(conc)
        nimpl = len(qs)
        _QCACHE.append((qs + L.controls(), nimpl))
    return _QCACHE[0]


def lemmas_ok(lem):
    """The lemma file must be for this source and the complete query list (an empty or foreign file is NOT PROVED)."""
    try:
        qs, nimpl = expected_queries()
        src = _sha_file('lemmas_z3.py')
        good = bool(lem) and lem.get('ok') is True and lem['proved'] == lem['implications'] and lem['controls_ok'] == lem['controls']
        hashes = (lem.get('source_sha256') == src and len(lem['rows']) == len(qs) and lem['implications'] == nimpl
                  and lem['controls'] == len(qs) - nimpl and lem['implications'] > 0
                  and all(r['name'] == q[0] and r['sha256'] == hashlib.sha256(q[1].encode()).hexdigest() and r['i'] == i
                          and r['ok'] is True for i, (r, q) in enumerate(zip(lem['rows'], qs))))
        return bool(good and hashes and lem.get('z3'))
    except (KeyError, TypeError, ValueError, AttributeError):
        return False


def canaries_summary(can):
    """-> (ok, summary).  Counted canary = damaged circuit that is NOT PROVED and has a witness (a wrong behaviour found)."""
    if not can or not isinstance(can.get('canaries'), list):
        return False, dict(total=0, counted=0, targeted=0, false_proved=[], reason='no canary evidence')
    rows = can['canaries']
    cnt = sum(bool(c.get('counted')) for c in rows)
    tcnt = sum(bool(c.get('counted')) for c in rows if c['name'].startswith('T'))
    false_proved = [c['name'] for c in rows if c.get('prove') == 'CHAIN_PROVED' and c.get('witness')]
    for c in can.get('extra', []):
        if c.get('prove') == 'CHAIN_PROVED' and c.get('witness'):
            false_proved.append(c['name'] + ' (extra search)')
    from ._client import same_client as _same
    same_client = _same(can.get('manifest'))      # the 6 client files, the whole bridge and qoptct/verify: all equal to the code judging now
    ok = cnt >= 20 and 2 * tcnt >= cnt and tcnt >= 10 and not false_proved and same_client and can.get('n') == 256
    return bool(ok), dict(total=len(rows), counted=cnt, targeted=tcnt, false_proved=false_proved, made_by_this_client=same_client)


def judge(rep, *, lemmas, prime, canaries=None, canaries_required=None, min_run=None, min_sim=None):
    """Verdict for one evidence dictionary ``rep`` (the result of ``run.run_pair``)."""
    n, p = rep['n'], rep['p']
    pr, ru, si = rep.get('prove') or {}, rep.get('run') or {}, rep.get('sim') or {}
    if 'error' in ru:                 # the execution mode failed (e.g. the wire flow): never PROVED
        ru = dict(violations=1, n_wrong_output=1, n_differs=1, n_pq_mismatch=1, pq_assumption_failures=[ru['error']],
                  iterations_traced=0, inputs=0)
    from ._client import same_client as _same_client
    ch = client_hashes()
    each = (rep.get('manifest_before') or {}).get('src', {}).get('each', {})
    need_can = (n == 256) if canaries_required is None else canaries_required
    can_ok, can_sum = canaries_summary(canaries) if canaries is not None else (False, dict(total=0, counted=0, targeted=0, false_proved=[]))
    min_run = min(p, 256) if min_run is None else min_run
    min_sim = min(p, 127) if min_sim is None else min_sim
    try:
        checks = dict(
            chain=pr.get('status') == 'CHAIN_PROVED' and not pr['admitted'] and pr['pos'] == pr['F'],
            assumptions_covered=bool(pr) and set(pr['assumed']) <= ASSUMPTIONS,
            run=ru['violations'] == 0 and ru['n_wrong_output'] == 0 and ru['n_differs'] == 0 and ru['n_pq_mismatch'] == 0
            and not ru['pq_assumption_failures'] and ru['iterations_traced'] == 2 * n,
            flow=flow_ok(pr.get('flow'), n) and flow_ok(ru.get('flow'), n),
            client_hash=each.get('kaliski/chain.py') == ch['chain_sha256'] and each.get('kaliski/lemmas_z3.py') == ch['lemmas_sha256']
            and _same_client(rep.get('manifest_before')),
            sim=si.get('gadget_error') is None and si.get('dirty_wires') == 0 and si.get('n_mismatch') == 0 and si.get('inputs', 0) >= min_sim > 0,
            sources_unchanged=rep.get('sources_unchanged_during_run') is True,
            lemmas=lemmas_ok(lemmas),
            prime=prime_ok(prime),
            samples=ru.get('inputs', 0) >= min_run > 0,
            canaries=can_ok if (need_can or canaries is not None) else True)
    except (KeyError, TypeError):
        checks = dict(evidence=False)
    status = 'PROVED' if all(checks.values()) else 'NOT PROVED (%s)' % ', '.join(k for k, v in checks.items() if not v)
    comp = rep.get('compose', {})
    return dict(
        status=status, n=n, p=p, checks=checks,
        program=dict(top_ops_sha256=comp.get('top_ops_sha256'), piece_qasm_sha256=comp.get('piece_qasm_sha256'),
                     ops=comp.get('nops'), wires=comp.get('nq'), measurements=comp.get('nmeas')),
        prove=dict(sec=pr.get('sec'), segments=pr.get('segments'), maxnodes=pr.get('maxnodes')),
        prime=prime,
        lemmas=dict(implications=(lemmas or {}).get('implications'), proved=(lemmas or {}).get('proved'),
                    controls=(lemmas or {}).get('controls'), controls_ok=(lemmas or {}).get('controls_ok'),
                    z3=(lemmas or {}).get('z3'), source_sha256=(lemmas or {}).get('source_sha256')),
        canaries=can_sum, client=ch, qualtran=(rep.get('manifest_before') or {}).get('qualtran_version'),
        manifest_sha256={g: (rep.get('manifest_before') or {}).get(g, {}).get('sha256') for g in ('src', 'bridge', 'qoptct', 'qualtran')},
        peak_rss_mb=rep.get('peak_rss_mb'), evidence_sec=rep.get('sec'))


def prime_ok(row):
    from .primes import prime_accepted
    return prime_accepted(row)


def default_canaries_path(n):
    from .._cache import cache_dir
    return os.path.join(str(cache_dir('kaliski')), 'canaries_n%d.json' % n)


def load_canaries(canaries, n):
    """``canaries``: None / 'auto' (the cache file written by ``canary.py`` for n = 256, else the evidence shipped with the
    package) / a path / a dict.  Whatever the source, the verdict accepts it only if this very client made it."""
    if isinstance(canaries, dict):
        return canaries
    path = default_canaries_path(n) if canaries in (None, 'auto') else canaries
    if canaries in (None, 'auto') and n == 256 and not os.path.exists(path):
        from ..registry import packaged_canaries      # evidence shipped with the package; judged by the same client-hash check
        path = packaged_canaries()
    if path and os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return None


def prove_kaliski(n, p, *, run_inputs=256, sim_inputs=127, lemmas=None, canaries='auto', prime=None, out=None, verbose=False):
    """Prove that ``KaliskiModInverse(bitsize=n, mod=p)`` of Qualtran 0.7.0 is correct, gate by gate.

    n        bit width of the register (p < 2^n - 1, p odd prime; x lives on n bits)
    p        the modulus
    lemmas   Z3 lemma result (dict or path); default: computed once and cached in $QOPT_CACHE/kaliski (a few minutes)
    canaries canary evidence (dict, path or 'auto' = cache file); REQUIRED at n = 256; never run here (about 45 minutes), see
             ``python -m qoptct.foreign.qualtran.kaliski canaries``
    prime    primality evidence (dict from ``primes.certify_prime``); default: computed
    out      write the verdict (with the evidence) to this JSON file
    Returns the verdict dictionary: ``status`` is ``'PROVED'`` or ``'NOT PROVED (<failed checks>)'``.  Without Z3 the lemma check
    cannot be made and the status is never PROVED.
    """
    from . import lemmas_z3
    from .primes import certify_prime
    from .run import run_pair
    t0 = time.time()
    if prime is None:
        prime = certify_prime(p)
    if lemmas is None:
        try:
            lemmas = lemmas_z3.ensure_lemmas(verbose=verbose)
        except RuntimeError as e:
            lemmas = dict(ok=False, error=str(e))
    elif not isinstance(lemmas, dict):
        with open(lemmas) as f:
            lemmas = json.load(f)
    can = load_canaries(canaries, n) if (n == 256 or canaries not in (None, 'auto')) else None
    rep = run_pair(n, p, prove=True, run='sample', sim=sim_inputs, verbose=verbose, run_inputs=run_inputs)
    v = judge(rep, lemmas=lemmas, prime=prime, canaries=can)
    v['sec'] = round(time.time() - t0, 1)
    v['finished'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    if out:
        with open(out, 'w') as f:
            json.dump(dict(verdict=v, evidence=rep), f, indent=1, default=str)
    return v
