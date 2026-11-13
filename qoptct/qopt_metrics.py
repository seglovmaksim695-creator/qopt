"""Resource counts for flat, unitary Clifford+T Qiskit circuits.

The two depths are weighted longest paths through the supplied circuit's wire
dependencies. Every instruction propagates dependencies, including Clifford
gates with zero T weight. Barriers synchronize their participating wires.
``t_depth`` is therefore an upper bound for this ordering, not the minimum
T-depth obtainable by algebraic rewrites, commutation, or ancillas.

Supported gates: h, x, y, z, s, sdg, t, tdg, cx, cz, swap, id, and barrier.
SWAP is counted as one supplied two-qubit gate; no decomposition is performed.
Global phase is ignored. Measurements, classical bits, control flow, arbitrary
rotations, and composite/custom gates raise ValueError. This module neither
optimizes a circuit nor verifies its equivalence to another circuit.
"""

from qiskit import QuantumCircuit
from qiskit.circuit import Barrier
from qiskit.circuit.library import (
    CXGate,
    CZGate,
    HGate,
    IGate,
    SGate,
    SdgGate,
    SwapGate,
    TGate,
    TdgGate,
    XGate,
    YGate,
    ZGate,
)


_GATES = {
    "h": HGate,
    "x": XGate,
    "y": YGate,
    "z": ZGate,
    "s": SGate,
    "sdg": SdgGate,
    "t": TGate,
    "tdg": TdgGate,
    "cx": CXGate,
    "cz": CZGate,
    "swap": SwapGate,
    "id": IGate,
    "barrier": Barrier,
}


def circuit_metrics(qc: QuantumCircuit) -> dict:
    """Return counts and weighted dependency depths without changing ``qc``.

    Returns integer ``n_qubits``, ``cx_count``, ``cz_count``,
    ``two_qubit_count``, ``t_count`` (T plus T-dagger), ``t_depth``, and
    ``two_qubit_depth``, plus a ``depth_semantics`` description identifying
    these as dependency-based upper bounds for the supplied ordering.
    ``two_qubit_depth`` counts layers of all supplied two-qubit gates.

    Raises ValueError if the input or any instruction is unsupported. Classical
    registers are rejected even when unused, to keep this API strictly unitary.
    """
    if not isinstance(qc, QuantumCircuit):
        raise ValueError("Expected a Qiskit QuantumCircuit")
    if qc.num_clbits or getattr(qc, "num_vars", 0):
        raise ValueError("Classical bits and variables are unsupported")

    result = {
        "n_qubits": qc.num_qubits,
        "cx_count": 0,
        "cz_count": 0,
        "two_qubit_count": 0,
        "t_count": 0,
        "t_depth": 0,
        "two_qubit_depth": 0,
        "depth_semantics": "weighted_dependency_upper_bounds_for_supplied_ordering",
    }
    t_levels = [0] * qc.num_qubits
    two_qubit_levels = [0] * qc.num_qubits
    indices = {qubit: index for index, qubit in enumerate(qc.qubits)}

    for position, instruction in enumerate(qc.data):
        operation = instruction.operation
        name = operation.name
        expected = _GATES.get(name)
        if expected is None or operation.base_class is not expected:
            raise ValueError(f"Unsupported instruction {name!r} at index {position}")
        if operation.params or instruction.clbits or getattr(operation, "condition", None) is not None:
            raise ValueError(f"Parameterized or classical instruction {name!r} at index {position}")

        qubits = [indices[qubit] for qubit in instruction.qubits]
        is_t = int(name in ("t", "tdg"))
        is_two_qubit = int(name != "barrier" and len(qubits) == 2)
        result["t_count"] += is_t
        result["two_qubit_count"] += is_two_qubit
        result["cx_count"] += int(name == "cx")
        result["cz_count"] += int(name == "cz")

        t_level = max((t_levels[index] for index in qubits), default=0) + is_t
        two_qubit_level = max((two_qubit_levels[index] for index in qubits), default=0) + is_two_qubit
        for index in qubits:
            t_levels[index] = t_level
            two_qubit_levels[index] = two_qubit_level

    result["t_depth"] = max(t_levels, default=0)
    result["two_qubit_depth"] = max(two_qubit_levels, default=0)
    return result
