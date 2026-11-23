"""Strict gate-by-gate export of Qiskit circuits, without optimization or resynthesis.

No qoptct gadgets are built here: the result is the four allowed kinds of
leaf gates and the exact origin of every leaf. Composite operations are
expanded only from their own ``definition``.
"""

from __future__ import annotations

import math
import numbers
from typing import Any

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import ControlledGate
from qiskit.circuit.library import CCXGate, CXGate, SwapGate, XGate


class TranslationError(ValueError):
    """The source gate is not in the strict supported language."""


_LEAF_CLASSES = {
    XGate: ("x", 1),
    CXGate: ("cx", 2),
    CCXGate: ("ccx", 3),
    SwapGate: ("swap", 2),
}
# Qiskit uses its own immutable singleton subclasses. We accept exactly the
# types created by the library, not an arbitrary isinstance subclass.
_LEAF_TYPES = {
    typ: (cls, name, width)
    for cls, (name, width) in _LEAF_CLASSES.items()
    for typ in (cls, type(cls()))
}

_FORBIDDEN_NAMES = frozenset(
    {
        "id", "i", "barrier", "delay", "measure", "reset", "store",
        "if_else", "while_loop", "for_loop", "switch_case", "break_loop", "continue_loop",
        "h", "y", "z", "s", "sdg", "t", "tdg", "sx", "sxdg",
        "rx", "ry", "rz", "p", "u", "u1", "u2", "u3", "r",
        "cz", "cy", "cp", "crx", "cry", "crz", "rxx", "ryy", "rzz", "rzx",
        "global_phase",
    }
)


def _error(path: str, op: Any, message: str) -> TranslationError:
    name = getattr(op, "name", type(op).__name__)
    return TranslationError(f"{path}: gate {name!r} ({type(op).__name__}): {message}")


def _check_circuit(circuit: QuantumCircuit, path: str) -> None:
    phase = circuit.global_phase
    if not isinstance(phase, numbers.Real) or not math.isfinite(float(phase)) or phase != 0:
        raise TranslationError(f"{path}: global phase {phase!r}; only an exact zero is allowed")
    if circuit.num_clbits:
        raise TranslationError(f"{path}: classical bits are not supported by this adapter")
    if getattr(circuit, "num_vars", 0) or getattr(circuit, "num_stretches", 0):
        raise TranslationError(f"{path}: real-time/duration variables are not supported")
    if circuit.parameters:
        raise TranslationError(f"{path}: unbound circuit parameters: {tuple(circuit.parameters)!r}")


def _expected_matrix(name: str, width: int) -> np.ndarray:
    """Matrix of the rule in Qiskit numbering: the first argument is the least significant bit."""
    result = np.zeros((1 << width, 1 << width), dtype=complex)
    for column in range(1 << width):
        row = column
        if name == "x":
            row ^= 1
        elif name == "cx" and column & 1:
            row ^= 2
        elif name == "ccx" and column & 3 == 3:
            row ^= 4
        elif name == "swap" and ((column & 1) != ((column >> 1) & 1)):
            row ^= 3
        result[row, column] = 1
    return result


_MATRICES = {name: _expected_matrix(name, width) for name, width in _LEAF_CLASSES.values()}


def lower(circuit: QuantumCircuit, *, max_depth: int = 128, max_operations: int = 200_000) -> dict:
    """Expand the source circuit to X/CX/CCX/SWAP or report the exact place of refusal.

    Returns ``nq``, ``leaves``, ``registers``, ``source`` and ``flat_source``.
    Indices in ``registers`` follow the Qiskit order (LSB first).
    Clean inputs are not assumed here: a separate specification states them.
    No operation disappears through optimization or recognition of arithmetic.
    """
    if not isinstance(circuit, QuantumCircuit):
        raise TranslationError("root: a QuantumCircuit is expected")
    if type(max_depth) is not int or max_depth < 1 or max_depth > 128:
        raise TranslationError("max_depth must be an integer from 1 to 128")
    if type(max_operations) is not int or not 1 <= max_operations <= 200_000:
        raise TranslationError("max_operations must be an integer from 1 to 200000")

    _check_circuit(circuit, "root")
    leaves: list[dict] = []
    flat = QuantumCircuit(circuit.num_qubits, name="strict_source_leaves")
    visited_operations = 0

    def walk(block: QuantumCircuit, mapping: list[int], path: str, depth: int) -> None:
        nonlocal visited_operations
        if depth > max_depth:
            raise TranslationError(f"{path}: decomposition depth exceeded {max_depth}")
        _check_circuit(block, path)
        if len(mapping) != block.num_qubits or len(set(mapping)) != len(mapping):
            raise TranslationError(f"{path}: wrong or non-injective wire layout")

        for index, instruction in enumerate(block.data):
            op = instruction.operation
            current = f"{path}/{index}:{getattr(op, 'name', type(op).__name__)}"
            visited_operations += 1
            if visited_operations > max_operations:
                raise _error(current, op, f"more than {max_operations} operations to expand")
            if instruction.clbits or getattr(op, "num_clbits", 0):
                raise _error(current, op, "classical arguments are not supported")
            if getattr(op, "condition", None) is not None or getattr(op, "condition_bits", ()):
                raise _error(current, op, "classical control is not supported")
            if getattr(op, "params", ()):
                raise _error(current, op, "gate parameters are not supported")
            if getattr(op, "name", "") in _FORBIDDEN_NAMES:
                raise _error(current, op, "not in the table of allowed rules")

            try:
                qargs = [mapping[block.find_bit(q).index] for q in instruction.qubits]
            except (KeyError, IndexError, ValueError) as exc:
                raise _error(current, op, f"could not match the wires: {exc}") from exc
            if len(qargs) != op.num_qubits or len(set(qargs)) != len(qargs):
                raise _error(current, op, "number or repetition of quantum arguments")
            if any(q < 0 or q >= circuit.num_qubits for q in qargs):
                raise _error(current, op, "wire index outside the source circuit")

            leaf_info = _LEAF_TYPES.get(type(op))
            open_control = isinstance(op, ControlledGate) and op.ctrl_state != (1 << op.num_ctrl_qubits) - 1
            if leaf_info is not None and not open_control:
                cls, name, width = leaf_info
                if op.num_qubits != width or op.name != name or op.base_class is not cls:
                    raise _error(current, op, "wrong attributes of a known primitive")
                if name in ("cx", "ccx") and op.num_ctrl_qubits != width - 1:
                    raise _error(current, op, "wrong number of control wires")
                try:
                    native = np.asarray(op.to_matrix(), dtype=complex)
                except Exception as exc:
                    raise _error(current, op, f"primitive matrix is unavailable: {exc}") from exc
                if not np.array_equal(native, _MATRICES[name]):
                    raise _error(current, op, "the matrix of the source primitive does not match the rule")
                leaves.append({"name": name, "qubits": qargs, "path": current})
                # Keep the source primitive itself; never rebuild it from the name.
                flat.append(op, qargs)
                continue

            # Also for an open control we take exactly the own
            # definition of Qiskit. There is no home-made MCX decomposition here.
            try:
                definition = op.definition
            except Exception as exc:
                raise _error(current, op, f"own definition is unavailable: {exc}") from exc
            if not isinstance(definition, QuantumCircuit):
                raise _error(current, op, "no own decomposition as a QuantumCircuit")
            if definition.num_qubits != len(qargs):
                raise _error(current, op, "definition width does not match the gate")
            if not definition.data:
                raise _error(current, op, "an empty definition of a non-primitive gate is not accepted")
            walk(definition, qargs, current, depth + 1)

    walk(circuit, list(range(circuit.num_qubits)), "root", 0)
    registers = {
        reg.name: [circuit.find_bit(qubit).index for qubit in reg]
        for reg in circuit.qregs
    }
    if len(registers) != len(circuit.qregs):
        raise TranslationError("root: ambiguous names of quantum registers")
    return {
        "nq": circuit.num_qubits,
        "leaves": leaves,
        "registers": registers,
        "source": circuit,
        "flat_source": flat,
    }
