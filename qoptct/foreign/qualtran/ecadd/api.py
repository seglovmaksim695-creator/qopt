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


def prove_ecadd(n, p, *, patched=True, steps=STEPS, timeout=None, log=None):
    """Prove ``ECAdd(n, p, window_size=1)``.

    ``patched=True``: fixed scheme; the verdict ``PROVED`` means that the six steps and their gluing prove
    "ECAdd = point addition in Montgomery form (all cases), a, b, lam_r unchanged, everything else clean", for all
    inputs with every word < p, under the contracts of the multiplier and the inverse (``contracts_status``).
    ``patched=False``: original 0.7.0; ``PROVED_OUTSIDE_G1_G3`` means the same exactly outside the groups G1-G3 of
    finding 05 (on the groups the scheme leaves a dirty flag; that is not part of this verdict).
    ``steps``: subset of 1..6 (a partial run is never PROVED).  ``log(msg)``: optional progress callback.
    """
    say = log or (lambda m: None)
    kind = 'e' if patched else None
    t0 = time.time()
    out = dict(api='prove_ecadd', n=n, p=str(p), patched=bool(patched), qualtran_kind=expected_kind_name(kind),
               z3=_z3.version(), contracts_status=contracts_status(n, p, ('MM', 'KINV')), steps={})
    if out['z3'] is None:
        out.update(status='NOT_PROVED', detail='Z3 not found (verdicts without Z3 are never PROVED)')
        return out
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
    if ok:
        out['status'] = 'PROVED' if patched else 'PROVED_OUTSIDE_G1_G3'
        out['claim'] = ('ECAdd = point addition in Montgomery form' if patched else
                        'ECAdd (original 0.7.0) = point addition outside the groups G1-G3 of finding 05')
    else:
        out['status'] = 'NOT_PROVED'
        if not complete:
            out['detail'] = 'partial run (not all six steps)'
    out['wall_seconds'] = round(time.time() - t0, 1)
    return out


def prove_ecwindowaddr(n, p, window=1, *, timeout=None, log=None):
    """Prove ``ECWindowAddR`` (fixed scheme, ``add_window_size=window``, ``mul_window_size=1``): (c, Q) -> (c, Q + c*R).

    The curve is y^2 = x^3 + 2x + 3 and R a fixed point (``e4.base_point``), as in the proof.  The lookup of the
    table row c and its removal (rule C1) are proved by gates; the inner ``ECAdd`` is a box with the contract of
    ``prove_ecadd`` (``contracts_status(...)`` does not cover it: call ``prove_ecadd`` for the same n, p and combine).
    """
    say = log or (lambda m: None)
    t0 = time.time()
    out = dict(api='prove_ecwindowaddr', n=n, p=str(p), window=window, qualtran_kind='e4', z3=_z3.version(),
               depends_on=dict(box='ECADD', proved_by='ecadd', note='call prove_ecadd(n, p, patched=True) for this (n, p)'))
    r = _check_variant(run_json('e4', f'{_PKG}.e4', ['prove', '--n', n, '--p', str(p), '--w', window], timeout=timeout), 'e4')
    out['e4'] = {k: r.get(k) for k in ('status', 'detail', 'wall_seconds', 'peak_rss_mb', 'logical_ops', 'regions', 'R',
                                      'curve', 'stats', 'kind')}
    say(f'ecwindowaddr: {r.get("status")}')
    out['status'] = 'PROVED' if r.get('status') == 'PROVED' else 'NOT_PROVED'
    out['claim'] = '(c, Q) -> (c, Q + c*R) for all c < 2^window and all x, y < p, given the ECAdd contract'
    out['wall_seconds'] = round(time.time() - t0, 1)
    return out
