# SPDX-License-Identifier: Apache-2.0
"""Qiskit circuits with rotations pi*p/2^k: hand them to the dyadic prover (``qoptct.dyadic``).

``check_qiskit`` calls :func:`via_dyadic` when the translation meets a rotation that is not a multiple of pi/4
(``NeedsFractionalPhases``).  The circuit is turned into a ``DCircuit`` (``dyadic_io.from_qiskit_strict``: a gate is
taken by name only when it is an instance of a standard Qiskit gate class with the matching matrix; everything else is
unrolled to its Qiskit definition) and proved against

* a reference ``QuantumCircuit`` (``check_qiskit(circuit, reference=ref)``); its registers are matched to those of
  the circuit BY NAME (same name and size = same wires; every other reference wire is a fresh clean wire; ``ref_clean``
  as in the Clifford+T path); qubits without a unique register -> UNKNOWN, or
* a dyadic specification (``QFTSpec``, ``ArithSpec`` / ``adder_spec`` / ``multiplier_spec``, ``ClassicalSpec``;
  ``check_qiskit(circuit, spec)``).

Without such a target nothing is decided: the result stays UNKNOWN with ``needs == 'fractional_phases'``.
Float angles are read as multiples of pi/2^k by ``qoptct.dyadic`` (a tolerance of 32 units in the last place of the
angle reduced exactly modulo 4 pi; an angle that is not within it of such a multiple makes the answer UNKNOWN).  ``clean`` registers start in
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


def _register_layout(circuit, what):
    """name -> list of wire indices (``find_bit`` order); every qubit must belong to exactly one register."""
    lay = {}
    for q in circuit.qubits:
        regs = circuit.find_bit(q).registers
        if len(regs) != 1:
            raise ValueError('%s: every qubit must belong to exactly one named register (found %d)' % (what, len(regs)))
    for r in circuit.qregs:
        lay[r.name] = [circuit.find_bit(q).index for q in r]
    return lay


def _reference_wires(circuit, reference, clean, ref_clean):
    """Place the reference on the wires of the circuit BY REGISTER NAME (as ``check.py::_map_reference`` does for the
    Clifford+T path): registers with the same name and size are shared (a register that is in ``ref_clean`` but not in
    ``clean`` is not shared), every other reference wire is a fresh clean wire after the circuit.
    -> (wire map of the reference, number of wires, number of fresh wires)."""
    lay = _register_layout(circuit, 'circuit')
    rlay = _register_layout(reference, 'reference')
    for nm in ref_clean:
        if nm not in rlay:
            raise ValueError('unknown ref_clean registers: %s' % [nm])
    wmap = [None] * reference.num_qubits
    for name, ws in rlay.items():
        if name in lay and len(lay[name]) == len(ws) and not (name in set(ref_clean) and name not in set(clean)):
            for a, b in zip(ws, lay[name]):
                wmap[a] = b
    nxt = circuit.num_qubits
    for i in range(len(wmap)):
        if wmap[i] is None:
            wmap[i] = nxt
            nxt += 1
    return wmap, nxt, nxt - circuit.num_qubits


def dyadic_target(spec, reference, circuit=None, clean=(), ref_clean=()):
    """The target for ``prove_dyadic`` or None when the call carries nothing the dyadic prover can use.  A reference
    circuit is read together with ``circuit`` (register names decide the wires): -> (target, fresh wires)."""
    from qoptct.dyadic import ArithSpec, ClassicalSpec, QFTSpec
    from .dyadic_io import from_qiskit_strict
    if reference is not None:
        if circuit is None:
            raise ValueError('a reference needs the circuit to match registers by name')
        wmap, nq, fresh = _reference_wires(circuit, reference, clean, ref_clean)
        return from_qiskit_strict(reference, wire_map=wmap, n=nq), fresh
    if isinstance(spec, (QFTSpec, ArithSpec, ClassicalSpec)):
        return spec, 0
    return None


def via_dyadic(out, circuit, spec, reference, clean=(), ref_clean=(), budget_sec=300, phase='global'):
    """Update the result dict ``out`` of ``check_qiskit`` with the verdict of the dyadic prover."""
    try:
        got = dyadic_target(spec, reference, circuit, clean, ref_clean)
        target, fresh = got if got is not None else (None, 0)
    except Exception as e:     # Unsupported gates, registers that cannot be matched by name etc.: UNKNOWN, never PROVED
        out.update(detail='%s; dyadic: reference not translated (%s: %s)' % (out.get('detail', ''), type(e).__name__, e))
        return out
    if target is None:
        out['detail'] = (out.get('detail', '') + '; give a QuantumCircuit reference or a dyadic spec '
                         '(qoptct.dyadic.QFTSpec / ArithSpec / ClassicalSpec) to prove it').lstrip('; ')
        return out
    try:
        from qoptct.dyadic import DCircuit, prove_dyadic
        from .dyadic_io import from_qiskit_strict
        anc = set(_wires(circuit, clean))
        if isinstance(target, DCircuit):
            c = from_qiskit_strict(circuit, n=target.n)           # the fresh reference wires are idle wires of the circuit
            anc |= set(range(circuit.num_qubits, target.n))
        else:
            c = from_qiskit_strict(circuit)
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
