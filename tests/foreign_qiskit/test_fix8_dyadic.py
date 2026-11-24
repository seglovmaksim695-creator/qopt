"""Holes H1-H3 of the 0.8 audit: register matching by name, standard gate classes only, large float angles."""
import math
import warnings

import pytest

pytest.importorskip("qiskit")
warnings.filterwarnings("ignore")

from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library import QFTGate
from qiskit.quantum_info import Operator

from qoptct.dyadic import DCircuit, Pi, prove_equal
from qoptct.dyadic.circuit import NotDyadic, recognise_float
from qoptct.foreign.qiskit import check_qiskit


def test_h1_auditor_example_not_proved():
    # R: cp on (b, a) and H on register b: a different circuit by register names
    a, b = QuantumRegister(1, 'a'), QuantumRegister(1, 'b')
    C = QuantumCircuit(a, b)
    C.cp(math.pi / 8, a[0], b[0])
    C.h(a[0])
    a2, b2 = QuantumRegister(1, 'a'), QuantumRegister(1, 'b')
    R = QuantumCircuit(b2, a2)
    R.cp(math.pi / 8, b2[0], a2[0])
    R.h(b2[0])
    r = check_qiskit(C, reference=R)
    assert r['status'] in ('REFUTED', 'UNKNOWN'), r
    assert r['status'] != 'PROVED'


def test_h1_same_names_other_order_is_proved_when_equal():
    # the same circuit, registers declared in the other order: equal by name
    a, b = QuantumRegister(2, 'a'), QuantumRegister(2, 'b')
    C = QuantumCircuit(a, b)
    C.cp(math.pi / 8, a[0], b[1])
    C.h(a[1])
    a2, b2 = QuantumRegister(2, 'a'), QuantumRegister(2, 'b')
    R = QuantumCircuit(b2, a2)
    R.cp(math.pi / 8, a2[0], b2[1])
    R.h(a2[1])
    r = check_qiskit(C, reference=R)
    assert r['status'] == 'PROVED', r
    # and a one-gate change is refuted
    R2 = QuantumCircuit(b2, a2)
    R2.cp(math.pi / 8, a2[0], b2[1])
    R2.h(b2[1])
    assert check_qiskit(C, reference=R2)['status'] == 'REFUTED'


def test_h1_register_sizes_differ_is_not_proved():
    a = QuantumRegister(2, 'a')
    C = QuantumCircuit(a)
    C.cp(math.pi / 8, 0, 1)
    b = QuantumRegister(2, 'b')                     # other name: the reference acts on fresh wires
    R = QuantumCircuit(b)
    R.cp(math.pi / 8, 0, 1)
    assert check_qiskit(C, reference=R)['status'] != 'PROVED'


def test_h1_qubit_in_two_registers_is_unknown():
    a = QuantumRegister(2, 'a')
    C = QuantumCircuit(a)
    C.cp(math.pi / 8, 0, 1)
    R = QuantumCircuit(a)
    R.cp(math.pi / 8, 0, 1)
    R.add_register(QuantumRegister(1, 'extra'))
    # a plain circuit still works; a bit shared by two registers is the case that must not be guessed
    from qiskit.circuit import Qubit
    q0, q1 = Qubit(), Qubit()
    R2 = QuantumCircuit([q0, q1], QuantumRegister(name='r1', bits=[q0, q1]), QuantumRegister(name='r2', bits=[q0]))
    R2.cp(math.pi / 8, q0, q1)
    r = check_qiskit(C, reference=R2)
    assert r['status'] == 'UNKNOWN', r


def test_h1_ref_clean_register_is_a_fresh_wire():
    # the reference's own 'anc' register (ref_clean) is not shared with the circuit's 'anc'
    d1, an1 = QuantumRegister(2, 'd'), QuantumRegister(1, 'anc')
    C = QuantumCircuit(d1, an1)
    C.cp(math.pi / 8, 0, 1)
    C.x(2)
    d2, an2 = QuantumRegister(2, 'd'), QuantumRegister(1, 'anc')
    R = QuantumCircuit(d2, an2)
    R.cp(math.pi / 8, 0, 1)
    r = check_qiskit(C, reference=R, clean=['anc'], ref_clean=['anc'])
    assert r['status'] != 'PROVED', r             # the circuit leaves anc in |1>
    r2 = check_qiskit(C, reference=R, clean=['anc'], ref_clean=['nope'])
    assert r2['status'] == 'UNKNOWN', r2


def test_h2_custom_gate_named_h_is_not_h():
    sub = QuantumCircuit(1, name='h')
    sub.s(0)
    C = QuantumCircuit(2)
    C.append(sub.to_gate(), [0])
    C.cp(math.pi / 8, 0, 1)
    R = QuantumCircuit(2)
    R.h(0)
    R.cp(math.pi / 8, 0, 1)
    assert not Operator(C).equiv(Operator(R))
    assert check_qiskit(C, reference=R)['status'] != 'PROVED'
    # the custom gate is proved equal to its own definition, i.e. it is unrolled
    R2 = QuantumCircuit(2)
    R2.s(0)
    R2.cp(math.pi / 8, 0, 1)
    assert check_qiskit(C, reference=R2)['status'] == 'PROVED'


_STD = {'t': ('t', 1, ()), 'x': ('x', 1, ()), 's': ('s', 1, ()), 'p': ('p', 1, (0.25,)), 'rz': ('rz', 1, (0.25,)),
        'swap': ('swap', 2, ()), 'rzz': ('rzz', 2, (0.25,)), 'cp': ('cp', 2, (0.25,)), 'u': ('u', 1, (0.25, 0.5, 0.75))}


@pytest.mark.parametrize('name', sorted(_STD))
def test_h2_custom_gate_with_a_standard_name_is_unrolled(name):
    from qiskit.circuit import Gate
    _, nq, params = _STD[name]

    class Fake(Gate):
        def __init__(self):
            super().__init__(name, nq, list(params))

        def _define(self):
            q = QuantumCircuit(nq)
            q.h(0)                              # not the gate its name says
            self.definition = q
    C = QuantumCircuit(2)
    C.append(Fake(), list(range(nq)))
    C.cp(math.pi / 8, 0, 1)
    R = QuantumCircuit(2)
    getattr(R, name)(*params, *range(nq))
    R.cp(math.pi / 8, 0, 1)
    assert not Operator(C).equiv(Operator(R))
    assert check_qiskit(C, reference=R)['status'] != 'PROVED'
    # the same circuit with the real definition spelled out is what the prover sees
    R2 = QuantumCircuit(2)
    R2.h(0)
    R2.cp(math.pi / 8, 0, 1)
    assert check_qiskit(C, reference=R2)['status'] == 'PROVED'


def test_h2_subclass_of_standard_gate_is_unrolled():
    from qiskit.circuit.library import HGate

    class MyH(HGate):
        def to_matrix(self):                       # lies about itself
            import numpy as np
            return np.eye(2, dtype=complex)

        def _define(self):
            q = QuantumCircuit(1)
            q.s(0)
            self.definition = q
    C = QuantumCircuit(2)
    C.append(MyH(), [0])
    C.cp(math.pi / 8, 0, 1)
    R = QuantumCircuit(2)
    R.h(0)
    R.cp(math.pi / 8, 0, 1)
    assert check_qiskit(C, reference=R)['status'] != 'PROVED'


def test_h2_library_blocks_still_unrolled():
    n = 3
    C = QuantumCircuit(n)
    C.append(QFTGate(n), range(n))
    from qiskit.synthesis.qft import synth_qft_full
    assert check_qiskit(C, reference=synth_qft_full(n, do_swaps=True))['status'] == 'PROVED'


def test_h3_large_float_angle_is_not_read_as_dyadic():
    th = 2 * math.pi * 2 ** 26 + math.pi / 8 + 1e-6
    with pytest.raises(NotDyadic):
        recognise_float(th)
    c = DCircuit(1)
    with pytest.raises(NotDyadic):
        c.add('rz', (0,), th)
    # the angle that is really dyadic and large is refused as well (the float cannot carry it), not read approximately
    with pytest.raises(NotDyadic):
        recognise_float(2 * math.pi * 2 ** 26 + math.pi / 8)


def test_h3_angles_up_to_4pi_unchanged_and_small_reduce_exactly():
    from fractions import Fraction
    assert recognise_float(math.pi / 8) == Fraction(1, 8)
    c = DCircuit(1)
    c.add('rz', (0,), 3 * math.pi + math.pi / 8)
    ref = DCircuit(1)
    ref.add('rz', (0,), Pi(25, 8))
    assert prove_equal(c, ref, phase='exact')['verdict'] == 'PROVED'
