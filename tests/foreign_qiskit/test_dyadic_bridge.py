"""Qiskit -> qoptct.dyadic: rotations that are not multiples of pi/4 are proved (or refuted) by the dyadic prover."""
import math
import warnings

import pytest

pytest.importorskip("qiskit")
warnings.filterwarnings("ignore")

from qiskit import QuantumCircuit
from qiskit.circuit.library import DraperQFTAdder, QFTGate

from qoptct.dyadic import DCircuit, QFTSpec, adder_spec
from qoptct.foreign.qiskit import check_qiskit


def qft(n):
    qc = QuantumCircuit(n)
    qc.append(QFTGate(n), range(n))
    return qc


@pytest.mark.parametrize('n', [2, 3, 5, 8])
def test_qft_against_spec_proved(n):
    r = check_qiskit(qft(n), QFTSpec(range(n)))
    assert r['status'] == 'PROVED', r
    assert r['method'].startswith('dyadic') and 'needs' not in r


def test_qft_against_reference_circuit_proved():
    from qiskit.synthesis.qft import synth_qft_full
    n = 4
    r = check_qiskit(qft(n), reference=synth_qft_full(n, do_swaps=True))
    assert r['status'] == 'PROVED', r


def test_without_target_stays_unknown_with_needs():
    from qoptct.foreign.qiskit import Spec
    r = check_qiskit(qft(3), Spec('nothing', lambda v: {}))
    assert r['status'] == 'UNKNOWN' and r.get('needs') == 'fractional_phases', r
    assert 'dyadic' not in r


def test_wrong_spec_is_refuted():
    r = check_qiskit(qft(4), QFTSpec(range(4), inverse=True))
    assert r['status'] == 'REFUTED', r
    assert r.get('witness') is not None


def test_draper_adder_proved_with_lemma_assumption():
    n = 4
    r = check_qiskit(DraperQFTAdder(n, kind='fixed').decompose(), adder_spec(range(n), range(n, 2 * n)))
    assert r['status'] == 'PROVED', r
    assert any('fourier-shift' in a for a in r['assumptions']), r['assumptions']


def test_draper_adder_wrong_spec_refuted():
    n = 4
    r = check_qiskit(DraperQFTAdder(n, kind='fixed').decompose(), adder_spec(range(n - 1), range(n, 2 * n)))
    assert r['status'] == 'REFUTED', r


def _bent(n, bend_at=None, angle=None):
    """Hand-built QFT (H, controlled phases pi/2^d, swaps); the angle of the controlled phase number ``bend_at`` is replaced."""
    qc = QuantumCircuit(n)
    idx = 0
    for j in reversed(range(n)):
        qc.h(j)
        for k in reversed(range(j)):
            ang = math.pi / 2 ** (j - k)
            if idx == bend_at:
                ang = angle
            qc.cp(ang, k, j)
            idx += 1
    for i in range(n // 2):
        qc.swap(i, n - 1 - i)
    return qc


def test_hand_built_qft_is_proved():
    r = check_qiskit(_bent(5), QFTSpec(range(5)))
    assert r['status'] == 'PROVED', r


def test_corrupted_angle_to_another_dyadic_angle_is_refuted():
    n = 4
    bad = _bent(n, 3, math.pi / 2 ** 7)
    r = check_qiskit(bad, QFTSpec(range(n)))
    assert r['status'] == 'REFUTED', r


@pytest.mark.parametrize('delta', [1e-3, 1e-6, 1e-9])
def test_corrupted_angle_off_the_dyadic_grid_is_unknown_never_proved(delta):
    n = 4
    bad = _bent(n, 3, math.pi / 2 ** 3 + delta)
    r = check_qiskit(bad, QFTSpec(range(n)))
    assert r['status'] == 'UNKNOWN', r
    assert 'dyadic' in r['detail'] or r['status'] != 'PROVED'


def test_clean_register_is_an_ancilla():
    n = 3
    from qiskit import QuantumRegister
    data, anc = QuantumRegister(n, 'd'), QuantumRegister(1, 'anc')
    c2 = QuantumCircuit(data, anc)
    c2.append(QFTGate(n), range(n))
    c2.x(3); c2.x(3)
    r = check_qiskit(c2, QFTSpec(range(n)), clean=['anc'])
    assert r['status'] == 'PROVED', r


def test_dcircuit_reference_not_accepted_as_qiskit_reference():
    r = check_qiskit(qft(3), reference=DCircuit(3))
    assert r['status'] == 'UNKNOWN', r
