"""Clifford+T programs: the program that ProjectQ decomposed to H, S, T, CNOT is compared with the same program at the NCT
level (multi-controlled X, no decomposition) by the exact path-sum prover ``qoptct.verify.prove`` (clean-ancilla model).

Wires of the two recordings differ (decompositions allocate their own ancillas), so both are renumbered the same way:
input registers (sorted by name), then output-only registers, then the remaining wires.  The remaining wires of both
recordings are clean ancillas: they start in |0> and must end in |0> in the circuit and in the reference.
"""
from __future__ import annotations

import os
import tempfile

from ._logical import emit_qasm

__all__ = ['canonical_perm', 'check_clifford_t']


def canonical_perm(inputs, outputs, n):
    order = [w for name in sorted(inputs) for w in inputs[name]]
    seen = set(order)
    for name in sorted(outputs):
        for w in outputs[name]:
            if w not in seen:
                order.append(w)
                seen.add(w)
    order += [w for w in range(n) if w not in seen]
    return {old: new for new, old in enumerate(order)}


def _same_layout(a, b):
    return ({k: len(v) for k, v in a.items()} == {k: len(v) for k, v in b.items()})


def check_clifford_t(lg_ct, rec_ref, budget_sec=300.0, max_width=22):
    """``lg_ct``: Logical in Clifford+T mode; ``rec_ref``: Recording at the NCT level with unbounded controls."""
    from ...verify import prove
    out = dict(method='pathsum-vs-nct-reference')
    if not (_same_layout(lg_ct.inputs, rec_ref.inputs) and _same_layout(lg_ct.outputs, rec_ref.outputs)):
        return dict(out, status='UNKNOWN', detail='the two recordings have different register layouts')
    if lg_ct.garbage or rec_ref.garbage:
        return dict(out, status='UNKNOWN', detail='declared garbage is not supported in Clifford+T mode')
    p_ct = canonical_perm(lg_ct.inputs, lg_ct.outputs, lg_ct.nq)
    p_ref = canonical_perm(rec_ref.inputs, rec_ref.outputs, rec_ref.n_wires)
    N = max(lg_ct.nq, rec_ref.n_wires)
    ops = []
    for op in lg_ct.ops:
        k = op[0]
        ops.append((k,) + tuple(p_ct[w] for w in op[1:]))
    free = [p_ref[w] for name in sorted(rec_ref.inputs) for w in rec_ref.inputs[name]]
    lines = ['.v ' + ' '.join(f'q{i}' for i in range(N)), '.i ' + ' '.join(f'q{i}' for i in sorted(free)), 'BEGIN']
    for op in rec_ref.ops:
        if op[0] == 'x':
            lines.append('tof ' + ' '.join(f'q{p_ref[w]}' for w in tuple(op[1]) + (op[2],)))
        elif op[0] not in ('alloc', 'free'):
            return dict(out, status='UNKNOWN', detail=f'the NCT reference contains {op[0]}')
    lines.append('END')
    fd, path = tempfile.mkstemp(suffix='.qc')
    try:
        with os.fdopen(fd, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')
        r = prove(emit_qasm(N, ops), path, clean=True, budget_sec=budget_sec, max_width=max_width)
    finally:
        os.unlink(path)
    return dict(out, status=r['status'], detail=str(r.get('detail'))[:300], prover=r.get('method'), time=r.get('time'),
                witness=r.get('counterexample'), nq=N, n_ops=len(ops))
