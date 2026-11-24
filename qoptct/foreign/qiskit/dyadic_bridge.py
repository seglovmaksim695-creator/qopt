# SPDX-License-Identifier: Apache-2.0
"""Qiskit circuits with rotations pi*p/2^k: hand them to the dyadic prover (``qoptct.dyadic``).

``check_qiskit`` calls :func:`via_dyadic` when the translation meets a rotation that is not a multiple of pi/4
(``NeedsFractionalPhases``).  The circuit is turned into a ``DCircuit`` (``qoptct.dyadic.qiskit_io.from_qiskit``;
library blocks are unrolled to their Qiskit definition) and proved against

* a reference ``QuantumCircuit`` (``check_qiskit(circuit, reference=ref)``), or
* a dyadic specification (``QFTSpec``, ``ArithSpec`` / ``adder_spec`` / ``multiplier_spec``, ``ClassicalSpec``;
  ``check_qiskit(circuit, spec)``).

Without such a target nothing is decided: the result stays UNKNOWN with ``needs == 'fractional_phases'``.
Float angles are read as exact multiples of pi/2^k by ``qoptct.dyadic`` (a tolerance of a few units in the last place;
an angle that is not close to such a multiple makes the answer UNKNOWN, never PROVED).  ``clean`` registers start in
|0> and must end in |0>.  The verdict holds up to a global phase (the contract of the Qiskit bridge).

Assumption of the answer for an ``ArithSpec`` (addition or multiplication in the Fourier basis): the fourier-shift lemma
(geometric sum, ``docs/DYADIC.md``, Lemma 2); it is listed in ``result['assumptions']``.
"""
from __future__ import annotations

__all__ = ['via_dyadic', 'dyadic_target']

_STATUS = {'PROVED': 'PROVED', 'REFUTED': 'REFUTED'}


def _wires(circuit, names):
    names = set(names)
    have = {r.name for r in circuit.qregs}
    missing = names - have
    if missing:
        raise ValueError('unknown clean registers: %s' % sorted(missing))
    return sorted(circuit.find_bit(q).index for r in circuit.qregs if r.name in names for q in r)


def is_dyadic_spec(spec):
    from qoptct.dyadic import ArithSpec, ClassicalSpec, QFTSpec
    return isinstance(spec, (QFTSpec, ArithSpec, ClassicalSpec))


def dyadic_target(spec, reference):
    """The target for ``prove_dyadic`` or None when the call carries nothing the dyadic prover can use."""
    from qoptct.dyadic import ArithSpec, ClassicalSpec, QFTSpec
    from qoptct.dyadic.qiskit_io import from_qiskit
    if reference is not None:
        return from_qiskit(reference)
    if isinstance(spec, (QFTSpec, ArithSpec, ClassicalSpec)):
        return spec
    return None


def via_dyadic(out, circuit, spec, reference, clean=(), ref_clean=(), budget_sec=300, phase='global'):
    """Update the result dict ``out`` of ``check_qiskit`` with the verdict of the dyadic prover."""
    try:
        target = dyadic_target(spec, reference)
    except Exception as e:     # Unsupported gates of the reference etc.: UNKNOWN, never PROVED
        out.update(detail='%s; dyadic: reference not translated (%s: %s)' % (out.get('detail', ''), type(e).__name__, e))
        return out
    if target is None:
        out['detail'] = (out.get('detail', '') + '; give a QuantumCircuit reference or a dyadic spec '
                         '(qoptct.dyadic.QFTSpec / ArithSpec / ClassicalSpec) to prove it').lstrip('; ')
        return out
    try:
        from qoptct.dyadic import DCircuit, prove_dyadic
        from qoptct.dyadic.qiskit_io import from_qiskit
        c = from_qiskit(circuit)
        anc = set(_wires(circuit, clean))
        if isinstance(target, DCircuit):
            anc |= set(range(target.n, c.n))
        r = prove_dyadic(c, target, phase=phase, ancillas=sorted(anc), budget_sec=budget_sec)
    except Exception as e:
        out.update(detail='%s; dyadic: not run (%s: %s)' % (out.get('detail', ''), type(e).__name__, e))
        return out
    v = r.get('verdict')
    out['dyadic'] = {k: r[k] for k in r if k not in ('rules_log',)}
    out['method'] = 'dyadic:' + str(r.get('method'))
    out['assumptions'] = list(r.get('assumptions') or [])
    if v in _STATUS:
        out['status'] = _STATUS[v]
        out.pop('needs', None)
        out['detail'] = 'dyadic prover: %s' % v
        for k in ('counterexample', 'witness'):
            if r.get(k) is not None:
                out['witness'] = r[k]
                break
    else:
        out['status'] = 'UNKNOWN'          # UNKNOWN, ERROR and INTERNAL never become a proof
        out['detail'] = 'dyadic prover: %s %s' % (v, str(r.get('reason', ''))[:300])
        if v == 'INTERNAL':
            out['internal_error'] = True
    return out
