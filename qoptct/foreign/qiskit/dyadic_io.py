# SPDX-License-Identifier: Apache-2.0
"""Strict Qiskit -> DCircuit reader for the dyadic path of ``check_qiskit``.

Differences from ``qoptct.dyadic.qiskit_io.from_qiskit`` (which stays as it is):

* an instruction is taken as a dyadic gate only when it is an instance of EXACTLY one of the standard gate classes of
  Qiskit (or of the immutable singleton variant Qiskit builds from it) (``qiskit.circuit.library.standard_gates``) and its matrix equals the matrix of the dyadic gate of that name;
  everything else (a user gate that merely has the name ``h``, a subclass of a standard gate, a library block) is
  unrolled by its definition; no definition -> ``Unsupported`` (UNKNOWN);
* the wires of the circuit are ``find_bit`` indices; a reference circuit is placed on the wires of the circuit by
  REGISTER NAME (``wire_map``), as the Clifford+T path of ``check_qiskit`` does.
"""
from __future__ import annotations

import inspect

import numpy as np

from qoptct.dyadic.circuit import DCircuit, GATES, Unsupported
from qoptct.dyadic.qiskit_io import _ALIAS, _CONTROLLED

__all__ = ['from_qiskit_strict', 'standard_gate_classes']

_STD = None
_MATRIX_OK = {}
_MAX_MATRIX_QUBITS = 3      # larger gates (mcx, ...) are checked by class, name, arity and control state only


def standard_gate_classes():
    """The set of the classes defined in ``qiskit.circuit.library.standard_gates`` (exact classes, no subclasses)."""
    global _STD
    if _STD is None:
        import qiskit.circuit.library.standard_gates as sg
        from qiskit.circuit import Gate
        found = set()
        for _, mod in inspect.getmembers(sg, inspect.ismodule):
            for _, cls in inspect.getmembers(mod, inspect.isclass):
                if issubclass(cls, Gate) and cls.__module__.startswith('qiskit.circuit.library.standard_gates'):
                    found.add(cls)
        _STD = frozenset(found)
    return _STD


def _std_class(op):
    """The standard Qiskit class of ``op`` (an exact instance, or the immutable ``_Singleton<Name>`` variant that Qiskit
    builds from it), else None.  A user subclass of a standard gate gives None."""
    t = type(op)
    std = standard_gate_classes()
    if t in std:
        return t
    bases = t.__bases__
    if (len(bases) == 2 and bases[1] in std and t.__name__ == '_Singleton' + bases[1].__name__
            and bases[0].__module__ == 'qiskit.circuit.singleton'):
        return bases[1]
    return None


def _direct(op):
    if _std_class(op) is None:
        return None
    nm = _ALIAS.get(op.name, op.name)
    if nm not in GATES:
        return None
    npar, arity = GATES[nm]
    if len(op.params) != npar:
        return None
    if arity is not None and op.num_qubits != arity:
        return None
    if nm in _CONTROLLED:
        k = getattr(op, 'num_ctrl_qubits', None)
        cs = getattr(op, 'ctrl_state', None)
        if k is None or cs is None or cs != (1 << k) - 1:
            return None
        if nm == 'cswap':
            if op.num_qubits != k + 2:
                return None
        elif op.num_qubits != k + 1:
            return None
    return nm


def _matrix_matches(op, nm):
    """The matrix Qiskit gives for ``op`` equals the dense matrix of the dyadic gate ``nm`` (same angles)."""
    if nm == 'gphase' or op.num_qubits > _MAX_MATRIX_QUBITS:
        return True
    from qoptct.dyadic import sim
    try:
        params = tuple(float(p) for p in op.params)
        key = (_std_class(op), nm, params, op.num_qubits)
    except (TypeError, ValueError):
        return False
    if key in _MATRIX_OK:
        return _MATRIX_OK[key]
    ok = False
    try:
        c = DCircuit(op.num_qubits)
        c.add(nm, range(op.num_qubits), *params)
        U = sim.unitary(c)
        V = np.asarray(op.to_matrix(), dtype=complex)
        ok = U.shape == V.shape and bool(np.max(np.abs(U - V)) < 1e-9)
    except Exception:
        ok = False
    _MATRIX_OK[key] = ok
    return ok


def from_qiskit_strict(qc, wire_map=None, n=None, max_depth=64):
    """QuantumCircuit -> DCircuit with ``n`` wires.  Qubit i of ``qc`` (``find_bit`` index) is wire ``wire_map[i]``
    (default: i); ``n`` defaults to the number of qubits of ``qc``."""
    wire_map = list(range(qc.num_qubits)) if wire_map is None else list(wire_map)
    c = DCircuit(qc.num_qubits if n is None else n, name=getattr(qc, 'name', None))

    def walk(circ, qmap, depth):
        if depth > max_depth:
            raise Unsupported('definitions nested too deep')
        gp = circ.global_phase
        if gp:
            c.add_global_phase(float(gp))
        for inst in circ.data:
            op = inst.operation
            if op.name == 'barrier':
                continue
            if inst.clbits or op.name in ('measure', 'reset', 'delay') or getattr(op, 'condition', None):
                raise Unsupported('non-unitary instruction: %s' % op.name)
            qs = [qmap[circ.find_bit(q).index] for q in inst.qubits]
            nm = _direct(op)
            if nm is not None and _matrix_matches(op, nm):
                c.add(nm, qs, *[float(p) for p in op.params])
                continue
            d = op.definition
            if d is None:
                raise Unsupported('gate %s has no definition' % op.name)
            walk(d, qs, depth + 1)
    walk(qc, wire_map, 0)
    return c
