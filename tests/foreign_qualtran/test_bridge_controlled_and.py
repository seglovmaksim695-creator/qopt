"""A1: raw QASM, source tensors, all domains/branches and coherent control.
"""
import pytest
pytest.importorskip("qualtran")
pytest.importorskip("cirq")

from hashlib import sha256
from importlib.metadata import version
from datetime import datetime, timezone
import inspect
import itertools
import json
import math
from pathlib import Path
import random
import resource
import sys
import time
import traceback

from qoptct.foreign.qualtran.bridge.controlled_and import ControlledAndError, expand, match, match_gate, require_domain
from qoptct.foreign.qualtran.bridge.emit import emit_qasm
from qoptct.foreign.qualtran.bridge.sim import run
from test_bridge_rules import basis, compare, dense_branch, random_state


ROOT = Path(__file__).resolve().parents[1]
TOLERANCE = 1e-10


class BitMasks:
    def __init__(self, count):
        self.F, self.T = 0, (1 << count) - 1

    def and_(self, a, b):
        return a & b

    def not_(self, a):
        return a ^ self.T


def make_leaf(cvs, uncompute=False):
    return dict(name="controlled_erase" if uncompute else "controlled_and",
                qubits=[0, 1, 2, 3], controls=list(cvs), path="test:A1")


def embed(vector):
    return list(vector) + [0j] * (32 - len(vector))


def source_partial(vector, tensor, cv, uncompute=False):
    """Controlled partial isometry from the REAL And.tensor_contract.

    Our index: c,a,b,t in increasing bit order. Qualtran tensor: MSB a,b,t.
    The inactive branch keeps the data and requires t=0. There is no full
    controlled-And-dagger Cirq instrument; the partial adjoint is used explicitly.
    """
    result = [0j] * 32
    for index, amplitude in enumerate(vector):
        if not amplitude:
            continue
        c, a, b, t = ((index >> i) & 1 for i in range(4))
        if c != cv:
            if t:
                raise ValueError("Inactive partial-isometry input has target=1")
            result[index] += amplitude
        elif uncompute:
            column = (a << 2) | (b << 1) | t
            for ao, bo in itertools.product((0, 1), repeat=2):
                result[c | (ao << 1) | (bo << 2)] += tensor[(ao << 1) | bo, column] * amplitude
        else:
            if t:
                raise ValueError("Compute input has target=1")
            column = (a << 1) | b
            for ao, bo, to in itertools.product((0, 1), repeat=3):
                result[c | (ao << 1) | (bo << 2) | (to << 3)] += tensor[(ao << 2) | (bo << 1) | to, column] * amplitude
    return result


def source_erase_branch(gate, vector, outcome):
    """Full branch of the uncontrolled source And-dagger, without replacing its gates.

    The matrices of the real unitary operations come from Cirq. The only source
    measurement is projected on the given outcome; the reset applies the matching
    real Kraus operator. This does not simulate a non-existent controlled reset.
    """
    import cirq
    import numpy as np
    a, b, t = map(cirq.NamedQubit, ("ctrl_0", "ctrl_1", "target"))
    order = [t, b, a]
    state = np.asarray(vector, dtype=np.complex128)
    measured = False
    for op in gate.to_clifford_t_circuit().all_operations():
        if type(op) is cirq.ClassicallyControlledOperation:
            conditions = tuple(op.classical_controls)
            if len(conditions) != 1 or type(conditions[0]) is not cirq.KeyCondition or str(conditions[0].key) != "target":
                raise ValueError("Неожиданное source классическое условие")
            if not measured:
                raise ValueError("Source читает измерение до его выполнения")
            if not outcome:
                continue
            op = op.without_classical_controls()
        if type(op.gate) is cirq.MeasurementGate:
            if measured or op.qubits != (t,) or cirq.measurement_key_name(op) != "target" or any(op.gate.invert_mask):
                raise ValueError("Неожиданное source измерение")
            state = np.asarray([amp if ((i >> 2) & 1) == outcome else 0 for i, amp in enumerate(state)])
            measured = True
        elif type(op.gate) is cirq.ResetChannel:
            if not measured or op.qubits != (t,):
                raise ValueError("Неожиданный source reset")
            matrix = cirq.kraus(op.gate)[outcome]
            # t — старший разряд в order=[t,b,a]. Kraus-оператор не унитарен:
            # нельзя оборачивать его в MatrixGate или нормировать эту ветвь.
            state = (matrix @ state.reshape(2, 4)).reshape(8)
        else:
            state = cirq.Circuit(op).unitary(qubit_order=order) @ state
    if not measured:
        raise ValueError("Source And† не содержит измерения")
    return state


def test_domains_and_recognition():
    import cirq
    from qualtran import CtrlSpec
    from qualtran._infra.controlled import Controlled
    from qualtran.bloqs.mcmt import And
    q = cirq.LineQubit.range(4)
    checks = 0
    for cvs in itertools.product((0, 1), repeat=3):
        cv, cv1, cv2 = cvs
        for inverse in (False, True):
            gate = And(cv1, cv2, uncompute=inverse).controlled(control_values=[cv])
            leaf = match(gate.on(*q), path="source")
            if leaf is None or leaf["controls"] != list(cvs) or leaf["name"] != ("controlled_erase" if inverse else "controlled_and"):
                raise AssertionError("Официальный controlled-And не распознан точно")
            leaf["qubits"] = [0, 1, 2, 3]
            for physical in range(32):
                values = {w: (physical >> w) & 1 for w in range(5)}
                predicate = int(tuple(values[w] for w in (0, 1, 2)) == cvs)
                valid = values[4] == 0 and values[3] == (predicate if inverse else 0)
                try:
                    require_domain(leaf, values, BitMasks(1), ancilla=4)
                except ControlledAndError:
                    if valid:
                        raise
                else:
                    if not valid:
                        raise AssertionError("Принята неверная область физического входа")
                checks += 1
            values = {w: sum(((i >> w) & 1) << i for i in range(8)) for w in range(3)}
            values[3] = sum(int(tuple((i >> w) & 1 for w in range(3)) == cvs) << i for i in range(8)) if inverse else 0
            values[4] = 0
            require_domain(leaf, values, BitMasks(8), ancilla=4)
    for bad_gate in (Controlled(And(), CtrlSpec(cvs=(1, 1))),
                     cirq.ControlledGate(And(), num_controls=1, control_values=[(0, 1)])):
        try:
            match_gate(bad_gate, list(range(cirq.num_qubits(bad_gate))), "reject")
        except ControlledAndError:
            pass
        else:
            raise AssertionError("Не отклонена неподдержанная форма управления")
    if match_gate(cirq.X, [0], "not-And") is not None:
        raise AssertionError("Принят другой вентиль")
    for leaf, kwargs in (
        (make_leaf((1, 1, 1)), dict(ancilla=4, nq=4, known_zero=[])),
        (make_leaf((1, 1, 1)), dict(ancilla=4, nq=5, known_zero=[3])),
        (make_leaf((1, 1, 1)), dict(ancilla=3, nq=4, known_zero=[3])),
    ):
        try:
            expand(leaf, **kwargs)
        except ControlledAndError:
            pass
        else:
            raise AssertionError("expand принял неподтверждённую чистоту/совпадение проводов")
    return dict(physical_inputs_checked=checks, exact_truth_mask_domains=16,
                symbolic_or_multi_control_rejected=True)


def test_source_and_all_branches():
    import cirq
    import numpy as np
    from qualtran.bloqs.mcmt import And
    from qoptct.verify.compose.lift import check_gadget_lemmas
    if not __debug__:
        raise ValueError("Python -O запрещён для точных qoptct gadget lemmas")
    if check_gadget_lemmas() is not True:
        raise AssertionError("Точные базовые фазовые леммы не прошли")
    rng = random.Random(20261121)
    errors, native_limits = [], []
    branches = 0
    for cvs in itertools.product((0, 1), repeat=3):
        cv, cv1, cv2 = cvs
        and_gate = And(cv1, cv2)
        tensor = np.asarray(and_gate.tensor_contract())
        dagger = np.asarray(and_gate.adjoint().tensor_contract())
        if tensor.shape != (8, 4) or dagger.shape != (4, 8):
            raise AssertionError("Неожиданная форма source And tensor")
        errors.append(compare(dagger.ravel(), tensor.conj().T.ravel()))
        compute = expand(make_leaf(cvs), ancilla=4, nq=4, known_zero=[3])
        erase = expand(make_leaf(cvs, True), ancilla=4, nq=5, known_zero=[4])
        source_controlled = and_gate.controlled(control_values=[cv])
        q = cirq.LineQubit.range(4)
        source_matrix = cirq.Circuit(source_controlled.on(*q)).unitary(qubit_order=list(reversed(q)))
        # qoptct.bdd повышает общий recursionlimit до 1<<20. Native Cirq
        # protocol может рекурсивно возвращаться к Controlled(And†); эта
        # дополнительная проба доступности не должна расходовать гигабайты.
        previous_recursion_limit = sys.getrecursionlimit()
        native_recursion_limit = min(previous_recursion_limit, 1500)
        try:
            sys.setrecursionlimit(native_recursion_limit)
            try:
                native_dagger = cirq.unitary(source_controlled.adjoint(), default=None)
                native_limits.append(dict(controls=list(cvs), native_controlled_adjoint_unitary=native_dagger is not None,
                                          recursion_limit=native_recursion_limit))
            except (ValueError, TypeError, RecursionError) as exc:
                native_limits.append(dict(controls=list(cvs), native_controlled_adjoint_unitary=False,
                                          recursion_limit=native_recursion_limit,
                                          detail=f"{type(exc).__name__}: {exc}"))
        finally:
            sys.setrecursionlimit(previous_recursion_limit)
        coherent = [0j] * 8
        index = (cv1 << 1) | (cv2 << 2)
        coherent[index], coherent[index | 1] = 1 / math.sqrt(2), 1j / math.sqrt(2)
        vectors = [basis(8, i) for i in range(8)] + [random_state(8, rng) for _ in range(12)] + [coherent]
        text_compute = emit_qasm(5, compute["ops"])
        text_erase = emit_qasm(5, erase["ops"])
        text_pair = emit_qasm(5, compute["ops"] + erase["ops"])
        for vector in vectors:
            initial = embed(vector)
            produced = source_partial(initial, tensor, cv)
            errors.append(compare(list(source_matrix @ np.asarray(initial[:16])) + [0j] * 16, produced))
            restored = source_partial(produced, dagger, cv, uncompute=True)
            errors.append(compare(restored, initial))
            for outcomes in itertools.product((0, 1), repeat=1):
                errors.append(compare(dense_branch(text_compute, initial, outcomes), [v / math.sqrt(2) for v in produced]))
                branches += 1
            for outcomes in itertools.product((0, 1), repeat=2):
                errors.append(compare(dense_branch(text_erase, produced, outcomes), [v / 2 for v in initial]))
                branches += 1
            for outcomes in itertools.product((0, 1), repeat=3):
                errors.append(compare(dense_branch(text_pair, initial, outcomes), [v / math.sqrt(8) for v in initial]))
                branches += 1
        for index in range(8):
            expected = index | (int(tuple((index >> w) & 1 for w in range(3)) == cvs) << 3)
            if run(5, compute["ops"], index) != (expected, [0]):
                raise AssertionError("Логическое compute/phase certificate нарушено")
            if run(5, erase["ops"], expected) != (index, [0, 0]):
                raise AssertionError("Логическое erase/phase certificate нарушено")
        # Реальный неконтролируемый And† source: отдельно все его ветви.
        for vector in [basis(4, i) for i in range(4)] + [random_state(4, rng) for _ in range(3)]:
            source_input = [0j] * 8
            for a, b in itertools.product((0, 1), repeat=2):
                amplitude = vector[a | (b << 1)]
                for ao, bo, t in itertools.product((0, 1), repeat=3):
                    source_input[ao | (bo << 1) | (t << 2)] += tensor[(ao << 2) | (bo << 1) | t, (a << 1) | b] * amplitude
            for outcome in (0, 1):
                errors.append(compare(source_erase_branch(and_gate.adjoint(), source_input, outcome),
                                      [v / math.sqrt(2) for v in vector + [0j] * 4]))
    return dict(control_patterns=8, valid_basis_inputs_per_rule=8,
                complex_vectors_per_pattern=13, raw_qasm_branches_checked=branches,
                exact_base_gadget_lemmas=True, exact_boolean_phase_checks=128,
                external_control_superposition=True, global_phase_compared=True,
                max_error=max(errors), source_controlled_adjoint_limits=native_limits,
                uncompute_reference="adjoint of source And partial tensor, with explicit coherent control")


def test_twenty_phase_canaries():
    rng = random.Random(21)
    vector = embed(random_state(8, rng))
    witnesses = []
    seen = set()
    for cvs in itertools.product((0, 1), repeat=3):
        compute = expand(make_leaf(cvs), ancilla=4, nq=4, known_zero=[3])
        erase = expand(make_leaf(cvs, True), ancilla=4, nq=5, known_zero=[4])
        original = compute["ops"] + erase["ops"]
        for pos, op in enumerate(original):
            if op[0] != "erase":
                continue
            mutated = list(original)
            mutated[pos] = ("erase", op[1], ())
            text = emit_qasm(5, mutated)
            fingerprint = sha256(text.encode()).hexdigest()
            if fingerprint in seen:
                continue
            seen.add(fingerprint)
            phase_vectors = {tuple(run(5, mutated, index)[1]) for index in range(8)}
            if len(phase_vectors) < 2:
                raise AssertionError("Мутант не показал входозависимую фазу")
            got = dense_branch(text, vector, (1, 1, 1))
            error = max(abs(a - b / math.sqrt(8)) for a, b in zip(got, vector))
            if error <= TOLERANCE:
                raise AssertionError("Удаление фазовой поправки не поймано")
            witnesses.append(dict(controls=list(cvs), operation=pos, sha256=fingerprint,
                                  defect="removed_CZ_in_measurement_correction", max_vector_error=float(error),
                                  phase_patterns=[list(x) for x in sorted(phase_vectors)]))
            if len(witnesses) == 20:
                return dict(count=20, targeted_measurement_canaries=20, witnesses=witnesses)
    raise AssertionError(f"Получено только {len(witnesses)} различных канареек")
