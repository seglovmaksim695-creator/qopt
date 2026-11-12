"""Resource accounting checks, including dependencies carried by zero-weight gates."""

import pathlib
import sys

import pytest
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit import Gate


from qopt.qopt_metrics import circuit_metrics


def test_empty_circuit():
    metrics = circuit_metrics(QuantumCircuit(3))
    assert metrics == {
        "n_qubits": 3,
        "cx_count": 0,
        "cz_count": 0,
        "two_qubit_count": 0,
        "t_count": 0,
        "t_depth": 0,
        "two_qubit_depth": 0,
        "depth_semantics": "weighted_dependency_upper_bounds_for_supplied_ordering",
    }
    assert circuit_metrics(QuantumCircuit())["t_depth"] == 0


def test_parallel_t_and_tdag_are_one_layer():
    qc = QuantumCircuit(3)
    qc.t(0)
    qc.tdg(1)
    qc.t(2)
    metrics = circuit_metrics(qc)
    assert metrics["t_count"] == 3
    assert metrics["t_depth"] == 1
    assert metrics["two_qubit_depth"] == 0


def test_sequential_t_layers_count_through_clifford_gates():
    qc = QuantumCircuit(2)
    qc.t(0)
    qc.h(0)
    qc.tdg(0)
    qc.t(1)
    metrics = circuit_metrics(qc)
    assert metrics["t_count"] == 3
    assert metrics["t_depth"] == 2


def test_cx_transfers_t_dependency_between_registers():
    left = QuantumRegister(1, "left")
    right = QuantumRegister(1, "right")
    qc = QuantumCircuit(left, right)
    qc.t(left[0])
    qc.cx(left[0], right[0])
    qc.t(right[0])
    metrics = circuit_metrics(qc)
    assert metrics["t_depth"] == 2
    assert metrics["two_qubit_depth"] == 1


def test_two_qubit_layers_and_swap_without_decomposition():
    qc = QuantumCircuit(4)
    qc.cx(0, 1)
    qc.cz(2, 3)
    qc.h(1)
    qc.cx(1, 2)
    qc.swap(0, 2)
    metrics = circuit_metrics(qc)
    assert metrics["cx_count"] == 2
    assert metrics["cz_count"] == 1
    assert metrics["two_qubit_count"] == 4
    assert metrics["two_qubit_depth"] == 3
    assert metrics["t_depth"] == 0


def test_barriers_synchronize_only_participating_qubits():
    qc = QuantumCircuit(3)
    qc.t(0)
    qc.barrier(1, 2)
    qc.t(1)
    assert circuit_metrics(qc)["t_depth"] == 1
    qc.barrier(0, 1)
    qc.tdg(0)
    metrics = circuit_metrics(qc)
    assert metrics["t_depth"] == 2
    assert metrics["two_qubit_count"] == 0


def test_barrier_transfers_two_qubit_dependency():
    qc = QuantumCircuit(4)
    qc.cx(0, 1)
    qc.barrier(1, 2)
    qc.cz(2, 3)
    assert circuit_metrics(qc)["two_qubit_depth"] == 2


def test_metrics_does_not_optimize_or_mutate_circuit():
    qc = QuantumCircuit(1, global_phase=0.3)
    qc.id(0)
    qc.x(0)
    qc.y(0)
    qc.z(0)
    qc.s(0)
    qc.sdg(0)
    qc.t(0)
    qc.tdg(0)
    original = qc.copy()
    metrics = circuit_metrics(qc)
    assert metrics["t_count"] == 2  # Even adjacent inverse gates are retained.
    assert metrics["t_depth"] == 2
    assert qc == original


@pytest.mark.parametrize("angle", [0.2, 0.0])
def test_arbitrary_rotations_are_explicitly_unsupported(angle):
    qc = QuantumCircuit(1)
    qc.rx(angle, 0)
    with pytest.raises(ValueError, match="Unsupported instruction 'rx'"):
        circuit_metrics(qc)


def test_composite_cannot_impersonate_a_supported_gate():
    qc = QuantumCircuit(1)
    qc.append(Gate("t", 1, []), [0])
    with pytest.raises(ValueError, match="Unsupported instruction 't'"):
        circuit_metrics(qc)


def test_measurement_and_classical_control_are_unsupported():
    qc = QuantumCircuit(1, 1)
    qc.measure(0, 0)
    with qc.if_test((qc.clbits[0], 1)):
        qc.x(0)
    with pytest.raises(ValueError, match="Classical"):
        circuit_metrics(qc)


def test_control_flow_without_classical_bits_is_unsupported():
    qc = QuantumCircuit(1)
    with qc.for_loop(range(2)):
        qc.t(0)
    with pytest.raises(ValueError, match="Unsupported instruction 'for_loop'"):
        circuit_metrics(qc)
