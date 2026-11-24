# SPDX-License-Identifier: Apache-2.0
"""Optional helper: Qiskit QuantumCircuit -> DCircuit (needs qiskit; used by the dyadic runs and tests).

Instructions whose name is in the supported gate set are taken as they are; every other instruction is replaced by
its Qiskit definition (recursively, with the definition's global phase), so library blocks (QFTGate, adders,
multipliers, controlled gates with open controls, ...) are unrolled down to the supported set.  A controlled gate is
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
    """QuantumCircuit -> DCircuit with the same qubit order (qubit i of qc = wire i)."""
    c = DCircuit(qc.num_qubits, name=getattr(qc, 'name', None))

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
            if nm is not None:
                c.add(nm, qs, *[float(p) for p in op.params])
                continue
            d = op.definition
            if d is None:
                raise Unsupported('gate %s has no definition' % op.name)
            walk(d, qs, depth + 1)
    walk(qc, list(range(qc.num_qubits)), 0)
    return c
