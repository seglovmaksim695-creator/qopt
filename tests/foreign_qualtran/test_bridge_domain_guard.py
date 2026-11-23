"""The exact partial-gate domain of AND/ERASE cannot be replaced by 'equal up to a global phase'."""
import itertools

import pytest

pytest.importorskip("qualtran")
pytest.importorskip("cirq")

from qoptct.foreign.qualtran.bridge.controlled_and import expand
from qoptct.foreign.qualtran.bridge.domain_guard import check_domains
from qoptct.foreign.qualtran.bridge.emit import emit_qasm
from qoptct.verify.bdd import ZERO

LAYOUT = dict(nq=5, inputs={'c': [0], 'a': [1], 'b': [2]}, outputs={'c': [0], 'a': [1], 'b': [2]},
              initial_zero=[3, 4], garbage=[], partial_gate_domains_required=True,
              domain_policy='exact_logical_AND_ERASE_v1')
CONTRACT = {'meaning': 'controlled And followed by its partial adjoint'}


@pytest.mark.parametrize('cvs', list(itertools.product((0, 1), repeat=3)))
def test_correct_pair_passes_and_constant_phase_wrong_domain_is_refuted(cvs):
    leaf = dict(name='controlled_and', qubits=[0, 1, 2, 3], controls=list(cvs))
    forward = expand(leaf, ancilla=4, nq=5, known_zero=[3, 4])['ops']
    backward = expand(dict(leaf, name='controlled_erase'), ancilla=4, nq=5, known_zero=[4])['ops']
    good = check_domains(emit_qasm(5, forward + backward), LAYOUT, contract=CONTRACT)
    assert good['status'] == 'PASS', good
    # A target complement leaves only a constant measurement phase but is OUTSIDE the declared partial domain.
    bad = check_domains(emit_qasm(5, forward + [('x', 3)] + backward), LAYOUT, contract=CONTRACT)
    assert bad['status'] == 'REFUTED', bad


def test_empty_domain_is_an_error():
    leaf = dict(name='controlled_and', qubits=[0, 1, 2, 3], controls=[1, 1, 1])
    forward = expand(leaf, ancilla=4, nq=5, known_zero=[3, 4])['ops']
    empty = check_domains(emit_qasm(5, forward), LAYOUT, contract=CONTRACT, pre=lambda S: ZERO)
    assert empty['status'] == 'ERROR', empty
