# SPDX-License-Identifier: Apache-2.0
"""Cross-checks with the results of 0.7.  Nothing here proves anything: these functions only establish that the
generic mechanism and the block-specific clients of 0.7 talk about the same program.

* ``compare_ecadd_step(n, p, k)`` lowers step ``k`` of ``ECAdd(n, p, 1)`` twice, with the 0.7 ECAdd client
  (``qualtran.ecadd.core.lower_step``) and with ``qualtran_adapter.lower``, and compares the operation lists one by
  one (needs Qualtran kind 'e').
* ``hashes_in_07(result_07, digests)`` says which of the given SHA-256 values occur in a 0.7 result file (the program
  hash and the QASM hashes of the pieces that 0.7 checked).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

__all__ = ['compare_ecadd_step', 'hashes_in_07', 'registry_07_entry']


def compare_ecadd_step(n, p, k, *, w=1):
    """dict(equal=..., ops_07, ops_auto, digest_07, digest_auto, same_inputs, opaque_in_07)"""
    from qualtran.bloqs.cryptography.ecc.ec_add import ECAdd
    from ..qualtran.ecadd.core import lower_step
    from .qualtran_adapter import lower, ops_digest
    steps = [bi.bloq for bi in sorted(ECAdd(n, p, w).decompose_bloq().bloq_instances, key=lambda b: b.i)
             if type(bi.bloq).__name__.startswith('_ECAddStep')]
    if not 1 <= k <= len(steps):
        raise ValueError(f'ECAdd has {len(steps)} steps, not {k}')
    step = steps[k - 1]
    low7, tr7 = lower_step(step)
    ops7 = tr7['ops']
    P = lower(step)
    inputs7 = {kk: [int(q) for q in v] for kk, v in low7['inputs'].items()}
    return dict(kind='compare_ecadd_step', n=n, p=p, w=w, step=k, step_bloq=type(step).__name__,
                ops_07=len(ops7), ops_auto=len(P.ops), digest_07=ops_digest(ops7), digest_auto=P.meta['ops_sha256'],
                equal=(list(ops7) == list(P.ops)), same_inputs=(inputs7 == P.meta['program_inputs']),
                opaque_in_07=sum(1 for o in ops7 if o[0] == 'opaque'))


_HEX = re.compile(r'\b[0-9a-f]{64}\b')


def hashes_in_07(result_07, digests):
    """the subset of ``digests`` (SHA-256 hex strings) that occur anywhere in the 0.7 result file ``result_07``"""
    text = Path(result_07).read_text()
    found = set(_HEX.findall(text))
    return sorted(d for d in set(digests) if d in found)


def registry_07_entry(name, n, p, w=None):
    """the entry of the 0.7 registry (``qoptct.foreign.qualtran.registry``) for one block and instance, or None"""
    from ..qualtran import registry as R07
    for e in R07.entries(name, n):
        if str(e['p']) == str(p) and (w is None or e.get('w') == w):
            return e
    return None


def main_compare(target, n, p, out=None):
    """``python -m qoptct.foreign.auto compare ecadd_stepK --n N --inprocess``"""
    if not target.startswith('ecadd_step'):
        raise ValueError('compare is defined for ecadd_step1 .. ecadd_step6')
    res = compare_ecadd_step(n, p, int(target[len('ecadd_step'):]))
    res['p'] = str(res['p'])
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text(json.dumps(res, indent=1))
    return res
