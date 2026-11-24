import pytest


def test_from_qiskit_does_not_trust_gate_names():
    pytest.importorskip('qiskit')
    from qiskit import QuantumCircuit
    from qiskit.circuit import Gate

    from qoptct.dyadic import DCircuit, prove_dyadic
    from qoptct.dyadic.qiskit_io import from_qiskit

    ref = DCircuit(1)
    ref.add('h', [0])
    good = QuantumCircuit(1)
    good.h(0)
    assert prove_dyadic(from_qiskit(good), ref)['verdict'] == 'PROVED'          # control: the standard H
    body = QuantumCircuit(1)
    body.s(0)
    fake = Gate('h', 1, [])
    fake.definition = body                                                    # named "h", acts as S
    bad = QuantumCircuit(1)
    bad.append(fake, [0])
    assert prove_dyadic(from_qiskit(bad), ref)['verdict'] != 'PROVED'
