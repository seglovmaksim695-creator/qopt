# SPDX-License-Identifier: Apache-2.0
"""Optional helper: Qiskit QuantumCircuit -> DCircuit (needs qiskit; used by the dyadic runs and tests).

A gate is taken as a dyadic gate only when it is an instance of a standard Qiskit gate class and (up to 3 qubits) its
matrix equals the matrix of the dyadic gate of that name (``qoptct.foreign.qiskit.dyadic_io.from_qiskit_strict``); every
other instruction - a user gate that merely carries a standard name, a library block (QFTGate, adders, multipliers,
controlled gates with open controls, ...) - is replaced by its Qiskit definition (recursively, with the definition's
global phase) down to the supported set.  A controlled gate is
taken directly only when all its controls are closed (ctrl_state all ones) and it acts on exactly controls + 1 qubits;
otherwise its definition is used.  Measurements, resets, classical conditions and gates without a definition raise
Unsupported.  Float angles go through the dyadic recognition of qoptct.dyadic.circuit (recorded in
DCircuit.angle_rec).
"""
from .circuit import DCircuit, GATES, Unsupported

__all__ = ['from_qiskit']

_ALIAS = {'cu1': 'cp', 'mcu1': 'mcp', 'mcphase': 'mcp', 'cphase': 'cp', 'c3x': 'mcx', 'c4x': 'mcx',
          'toffoli': 'ccx', 'cnot': 'cx', 'global_phase': 'gphase', 'u1': 'p', 'fredkin': 'cswap'}
_CONTROLLED = {'cx', 'cy', 'cz', 'cs', 'csdg', 'csx', 'cp', 'crz', 'ccx', 'ccz', 'cswap', 'mcx', 'mcz', 'mcp'}


def _direct(op):
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


def from_qiskit(qc, max_depth=64):
    """QuantumCircuit -> DCircuit with the same qubit order (qubit i of qc = wire i).

    Gate names are not trusted: this is ``qoptct.foreign.qiskit.dyadic_io.from_qiskit_strict`` (see the module text).
    """
    from qoptct.foreign.qiskit.dyadic_io import from_qiskit_strict
    return from_qiskit_strict(qc, max_depth=max_depth)
