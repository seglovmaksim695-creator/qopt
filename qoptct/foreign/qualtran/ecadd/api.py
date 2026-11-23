"""``prove_ecadd`` and ``prove_ecwindowaddr``: orchestration of the child-process proofs (see the package docstring)."""
from __future__ import annotations

import json
import os
import tempfile
import time

from .._variant import expected_kind_name, run_json
from . import _z3
from .contracts import contracts_status

_PKG = 'qoptct.foreign.qualtran.ecadd'
STEPS = (1, 2, 3, 4, 5, 6)


def _part(res, keys=('status', 'detail', 'wall_seconds', 'peak_rss_mb', 'kind', 'logical_ops', 'stats')):
    return {k: res[k] for k in keys if k in res}


def _check_variant(res, kind):
    """A part only counts when the child ran on exactly the Qualtran variant that was asked for."""
    want = expected_kind_name(kind)
    if res.get('kind') != want:
        return dict(res, status='NOT_PROVED', detail=f"child ran on Qualtran variant {res.get('kind')!r}, expected {want!r}")
    return res


def _assumed_parts(cs):
    """Names of the contract parts of ``contracts_status`` that are only ASSUMED for this instance."""
    return [f'{b} {part}' for b, st in cs.items() for part in ('forward', 'pair')
            if st.get(part) and st[part]['status'] != 'PROVED_FOR_THIS_INSTANCE']


def prove_missing_boxes(n, p, cs, *, log=None, timeout=None):
    """Prove, here and now, every box contract of ``cs`` (``contracts_status``) that no registry result covers for this (n, p).

    Runs ``prove_kaliski`` / ``prove_montgomery`` (w = 1) / ``check_pair_kaliski`` / ``check_pair_montgomery`` as needed and marks
    a part ``PROVED_FOR_THIS_INSTANCE`` (evidence ``proved in this run``) only when its status is exactly PROVED.  Returns the
    updated copy of ``cs``; a part that fails keeps ``ASSUMED_FOR_THIS_INSTANCE`` and gets ``attempt`` with the verdict status."""
    import copy
    say = log or (lambda m: None)
    cs = copy.deepcopy(cs)

    def attempt(box, part, fn):
        st = cs[box][part]
        if st is None or st['status'] == 'PROVED_FOR_THIS_INSTANCE':
            return
        say(f'proving {box} {part} for n={n} p={p}')
        try:
            r = fn()
            status = r.get('verdict') or r.get('status')
        except Exception as e:                       # noqa: BLE001  (a crash is never a proof)
            r, status = {}, f'ERROR {type(e).__name__}: {e}'
        st['attempt'] = status
        if status == 'PROVED':
            st['status'] = 'PROVED_FOR_THIS_INSTANCE'
            st['evidence'] = 'proved in this run'
    from ..kaliski import prove_kaliski
    from ..montgomery import prove_montgomery
    from ..pairs import check_pair_kaliski, check_pair_montgomery
    if 'MM' in cs:
        attempt('MM', 'forward', lambda: prove_montgomery(n, 1, p))
        attempt('MM', 'pair', lambda: check_pair_montgomery(n, 1, p, timeout=timeout))
    if 'KINV' in cs:
        attempt('KINV', 'forward', lambda: prove_kaliski(n, p))
        attempt('KINV', 'pair', lambda: check_pair_kaliski(n, p, timeout=timeout))
    return cs


def _domain_refusal(out, kind, n, p):
    from ..domain import check_domain
    d = check_domain(kind, n, p)
    out['domain'] = d
    if not d['ok']:
        out.update(status='NOT_PROVED', detail='outside the domain of the proof: ' + '; '.join(d['problems']))
        return True
    return False


def prove_ecadd(n, p, *, patched=True, steps=STEPS, timeout=None, log=None, prove_boxes=False):
    """Prove ``ECAdd(n, p, window_size=1)``.

    ``patched=True``: fixed scheme; the verdict ``PROVED`` means that the six steps and their gluing prove
    "ECAdd = point addition in Montgomery form (all cases), a, b, lam_r unchanged, everything else clean", for all
    inputs with every word < p, under the contracts of the multiplier and the inverse (``contracts_status``).
    ``patched=False``: original 0.7.0; ``PROVED_OUTSIDE_G1_G3`` means the same outside the groups G1-G3 of
    finding 05 (on the groups the scheme leaves a dirty flag; that is not part of this verdict).
    ``steps``: subset of 1..6 (a partial run is never PROVED).  ``log(msg)``: optional progress callback.

    The verdict ``PROVED`` additionally needs every box contract (multiplier and inverse, forward and pair) to be PROVED for
    exactly this (n, p) by a registry result or, with ``prove_boxes=True``, by a proof made in this call.  Otherwise a
    successful run is ``PROVED_CONDITIONAL`` (the field ``conditional_on`` lists the assumed parts).  The domain of p is checked
    first (``qoptct.foreign.qualtran.domain``: odd prime, 2p <= 2^n, p < 2^n - 1); outside it the status is ``NOT_PROVED``.
    """
    say = log or (lambda m: None)
    kind = 'e' if patched else None
    t0 = time.time()
    out = dict(api='prove_ecadd', n=n, p=str(p), patched=bool(patched), qualtran_kind=expected_kind_name(kind),
               z3=_z3.version(), contracts_status=contracts_status(n, p, ('MM', 'KINV')), steps={})
    if out['z3'] is None:
        out.update(status='NOT_PROVED', detail='Z3 not found (verdicts without Z3 are never PROVED)')
        return out
    if _domain_refusal(out, 'ecadd', n, p):
        out['claim'] = None
        return out
    if prove_boxes and _assumed_parts(out['contracts_status']):
        out['contracts_status'] = prove_missing_boxes(n, p, out['contracts_status'], log=say, timeout=timeout)
    ps = str(p)
    w = _check_variant(run_json(kind, f'{_PKG}.wiring', ['--n', n, '--p', ps], timeout=timeout), kind)
    out['wiring'] = dict(status=w.get('status'), kind=w.get('kind'), connections=w.get('connections'),
                         order=w.get('order'), detail=w.get('detail'))
    say(f'wiring: {w.get("status")}')
    ok = w.get('status') == 'OK'
    for k in steps:
        r = _check_variant(run_json(kind, f'{_PKG}.run_step',
                                    ['--step', k, '--n', n, '--p', ps], timeout=timeout), kind)
        out['steps'][str(k)] = _part(r)
        say(f'step {k}: {r.get("status")} ({r.get("wall_seconds")} s)')
        ok &= r.get('status') == 'PROVED'
    glue = None
    if ok or set(steps) != set(STEPS):
        with tempfile.TemporaryDirectory(prefix='qoptct_') as td:
            wf = os.path.join(td, 'wiring.json')
            json.dump(w, open(wf, 'w'))
            args = ['--n', n, '--p', ps, '--wiring', wf] + (['--mode', 'E2'] if patched else
                                                          ['--mode', 'E3', '--orig', '1,2,6'])
            glue = run_json(None, f'{_PKG}.zglue', args, timeout=timeout)
    else:
        glue = dict(status='NOT_RUN', detail='a step failed')
    out['glue'] = {k: glue.get(k) for k in ('status', 'mode', 'checks', 'answers', 'failed', 'z3_seconds', 'z3_version',
                                           'wiring_source', 'wiring_equals_reference', 'lemmas', 'detail')}
    say(f'glue: {glue.get("status")}')
    complete = set(steps) == set(STEPS)
    ok &= glue.get('status') == 'PROVED' and complete
    assumed = _assumed_parts(out['contracts_status'])
    out['conditional_on'] = assumed
    if ok:
        full = 'PROVED' if patched else 'PROVED_OUTSIDE_G1_G3'
        out['status'] = full if not assumed else 'PROVED_CONDITIONAL'
        out['claim'] = ('ECAdd = point addition in Montgomery form' if patched else
                        'ECAdd (original 0.7.0) = point addition outside the groups G1-G3 of finding 05 '
                        '(that the scheme is dirty on G1-G3 was shown in round 3, not by this package)')
        if assumed:
            out['detail'] = ('the six steps and their gluing are proved, but the box contracts are not proved for exactly this (n, p, w=1): '
                             + ', '.join(assumed) + f' (status {full} needs them; use --prove-boxes or the registry)')
    else:
        out['status'] = 'NOT_PROVED'
        if not complete:
            out['detail'] = 'partial run (not all six steps)'
    out['wall_seconds'] = round(time.time() - t0, 1)
    return out


def prove_ecwindowaddr(n, p, window=1, *, timeout=None, log=None, prove_boxes=False):
    """Prove ``ECWindowAddR`` (fixed scheme, ``add_window_size=window``, ``mul_window_size=1``): (c, Q) -> (c, Q + c*R).

    The curve is y^2 = x^3 + 2x + 3 and R a fixed point (``e4.base_point``), as in the proof.  The lookup of the
    table row c and its removal (rule C1) are proved by gates; the inner ``ECAdd`` is a box with the contract of
    ``prove_ecadd`` (``contracts_status(...)`` does not cover it: call ``prove_ecadd`` for the same n, p and combine).
    """
    say = log or (lambda m: None)
    t0 = time.time()
    out = dict(api='prove_ecwindowaddr', n=n, p=str(p), window=window, qualtran_kind='e4', z3=_z3.version(),
               depends_on=dict(box='ECADD', proved_by='ecadd', note='call prove_ecadd(n, p, patched=True) for this (n, p)'))
    if _domain_refusal(out, 'ecwindowaddr', n, p):
        out['claim'] = None
        return out
    from .contracts import evidence_id
    cs = contracts_status(n, p, ('ECADD',))
    out['contracts_status'] = cs
    if prove_boxes and _assumed_parts(cs):
        say(f'proving ECAdd (with its boxes) for n={n} p={p}')
        ea = prove_ecadd(n, p, patched=True, timeout=timeout, log=log, prove_boxes=True)
        out['ecadd'] = {k: ea.get(k) for k in ('status', 'detail', 'conditional_on', 'wall_seconds')}
        if ea.get('status') == 'PROVED':
            cs['ECADD']['forward'].update(status='PROVED_FOR_THIS_INSTANCE', evidence='proved in this run')
        else:
            cs['ECADD']['forward']['attempt'] = ea.get('status')
    r = _check_variant(run_json('e4', f'{_PKG}.e4', ['prove', '--n', n, '--p', str(p), '--w', window], timeout=timeout), 'e4')
    out['e4'] = {k: r.get(k) for k in ('status', 'detail', 'wall_seconds', 'peak_rss_mb', 'logical_ops', 'regions', 'R',
                                      'curve', 'stats', 'kind')}
    say(f'ecwindowaddr: {r.get("status")}')
    assumed = _assumed_parts(cs)
    out['conditional_on'] = assumed
    if r.get('status') != 'PROVED':
        out['status'] = 'NOT_PROVED'
    else:
        out['status'] = 'PROVED' if not assumed else 'PROVED_CONDITIONAL'
        if assumed:
            out['detail'] = ('the table lookup and its removal are proved, but ECAdd for exactly this (n, p) is not proved: '
                             + ', '.join(assumed) + ' (use --prove-boxes, or run prove_ecadd for this (n, p))')
    out['claim'] = '(c, Q) -> (c, Q + c*R) for all c < 2^window and all x, y < p, given the ECAdd contract'
    out['wall_seconds'] = round(time.time() - t0, 1)
    return out
