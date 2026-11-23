"""Numerical check of the translation rules, unnormalized branches and phase included.

A small own vector simulator parses the emitted text; the qoptct provers are not imported.
"""
import pytest
pytest.importorskip("qualtran")
pytest.importorskip("cirq")
import cmath
import itertools
import json
import math
from pathlib import Path
import random
import re
import time
import traceback

from qoptct.foreign.qualtran.bridge.emit import emit_qasm
from qoptct.foreign.qualtran.bridge.sim import SimulationError, run

TOLERANCE = 1e-10
SEED = 20261120
ROOT = Path(__file__).resolve().parents[1]


def _gate(statement):
    match = re.fullmatch(r'([a-z]+)\s+(q\[\d+\](?:\s*,\s*q\[\d+\])*)', statement.strip())
    if match is None:
        raise ValueError('Неподдержанное векторным симулятором выражение: ' + statement)
    return match.group(1), tuple(int(q) for q in re.findall(r'q\[(\d+)\]', match.group(2)))


def parse_qasm(text):
    """Independent strict parse of a small finite subset of QASM."""
    nq = None
    nb = 0
    program = []
    for line in text.splitlines():
        line = line.split('//', 1)[0].strip()
        if not line or line in ('OPENQASM 3;', 'OPENQASM 3.0;', 'include "stdgates.inc";'):
            continue
        if match := re.fullmatch(r'qubit\[(\d+)\] q;', line):
            if nq is not None:
                raise ValueError('Повторное объявление квантового регистра')
            nq = int(match.group(1))
            continue
        if match := re.fullmatch(r'bit\[(\d+)\] c;', line):
            nb = int(match.group(1))
            continue
        if match := re.fullmatch(r'c\[(\d+)\] = measure q\[(\d+)\];', line):
            k, q = map(int, match.groups())
            if not 0 <= k < nb:
                raise ValueError('Измерение вне классического регистра')
            program.append(('measure', q, k))
            continue
        if match := re.fullmatch(r'if \(c\[(\d+)\]\) \{(.*)\}', line):
            k, body = int(match.group(1)), match.group(2)
            statements = body.split(';')
            if statements[-1].strip():
                raise ValueError('Нет ; в условном блоке')
            program.append(('if', k, [_gate(s) for s in statements[:-1] if s.strip()]))
            continue
        if not line.endswith(';'):
            raise ValueError('Нет ;: ' + line)
        program.append(('gate', *_gate(line[:-1])))
    if nq is None:
        raise ValueError('Нет объявления квантовых проводов')
    return nq, program


def apply_gate(state, name, qs):
    """Small dense simulator; q[0] is the least significant bit of the index."""
    size = len(state)
    nq = size.bit_length() - 1
    arities = {'x': 1, 'z': 1, 'h': 1, 's': 1, 'sdg': 1, 't': 1, 'tdg': 1, 'cx': 2, 'cz': 2}
    if name not in arities or len(qs) != arities[name] or len(set(qs)) != len(qs):
        raise ValueError('Неподдержанный вентиль: %r' % ((name, qs),))
    if any(not 0 <= q < nq for q in qs):
        raise ValueError('Индекс провода вне диапазона')
    result = [0j] * size
    phases = {'z': -1, 's': 1j, 'sdg': -1j,
              't': cmath.exp(1j * math.pi / 4), 'tdg': cmath.exp(-1j * math.pi / 4)}
    for index, amplitude in enumerate(state):
        if name == 'x':
            result[index ^ (1 << qs[0])] += amplitude
        elif name == 'cx':
            result[index ^ (((index >> qs[0]) & 1) << qs[1])] += amplitude
        elif name == 'cz':
            result[index] += amplitude * (-1 if all((index >> q) & 1 for q in qs) else 1)
        elif name == 'h':
            q = qs[0]
            result[index & ~(1 << q)] += amplitude / math.sqrt(2)
            result[index | (1 << q)] += amplitude * (-1 if (index >> q) & 1 else 1) / math.sqrt(2)
        else:
            result[index] += amplitude * (phases[name] if (index >> qs[0]) & 1 else 1)
    return result


def dense_branch(text, state, outcomes=()):
    """Unnormalized vector of the given branch; keeps the global phase."""
    nq, program = parse_qasm(text)
    if len(state) != 1 << nq:
        raise ValueError('Неверный размер вектора')
    state = list(state)
    cbits = {}
    measurements = 0
    for op in program:
        if op[0] == 'gate':
            state = apply_gate(state, op[1], op[2])
        elif op[0] == 'measure':
            q, k = op[1:]
            if measurements >= len(outcomes) or outcomes[measurements] not in (0, 1):
                raise ValueError('Нет допустимого исхода измерения')
            outcome = outcomes[measurements]
            measurements += 1
            cbits[k] = outcome
            state = [a if (i >> q) & 1 == outcome else 0j for i, a in enumerate(state)]
        else:
            if op[1] not in cbits:
                raise ValueError('Условие читает ещё не измеренный бит')
            if cbits[op[1]]:
                for name, qs in op[2]:
                    state = apply_gate(state, name, qs)
    if measurements != len(outcomes):
        raise ValueError('Лишние исходы измерения')
    return state


def basis(size, index):
    return [complex(i == index) for i in range(size)]


def random_state(size, rng):
    vector = [complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(size)]
    norm = math.sqrt(sum(abs(x) ** 2 for x in vector))
    return [x / norm for x in vector]


def compare(got, expected):
    error = max((abs(a - b) for a, b in zip(got, expected)), default=0.0)
    if len(got) != len(expected) or error > TOLERANCE:
        raise AssertionError('Вектор/фаза не совпали: максимальная ошибка %.16g' % error)
    return error


def permuted(vector, permutation, output_size=None, factor=1):
    result = [0j] * (output_size or len(vector))
    for index, amplitude in enumerate(vector):
        result[permutation(index)] += factor * amplitude
    return result


def test_x_cx_swap():
    cases = [('X', 1, [('x', 0)], lambda i: i ^ 1),
             ('CX', 2, [('cx', 0, 1)], lambda i: i ^ ((i & 1) << 1)),
             ('SWAP', 2, [('cx', 0, 1), ('cx', 1, 0), ('cx', 0, 1)],
              lambda i: ((i & 1) << 1) | ((i >> 1) & 1))]
    rng = random.Random(SEED)
    errors = []
    checked = 0
    for name, nq, ops, permutation in cases:
        text = emit_qasm(nq, ops)
        for vector in [basis(1 << nq, i) for i in range(1 << nq)] + [random_state(1 << nq, rng) for _ in range(10)]:
            errors.append(compare(dense_branch(text, vector), permuted(vector, permutation)))
            checked += 1
        for i in range(1 << nq):
            assert run(nq, ops, i) == (permutation(i), [])
    return {'rules': ['X', 'CX', 'SWAP'], 'vectors': checked, 'max_error': max(errors)}


def controlled_and_ops(cv1, cv2, erase=False):
    wrappers = [('x', q) for q, cv in enumerate((cv1, cv2)) if cv == 0]
    operation = ('erase', 2, (('cz', 0, 1),)) if erase else ('and', 0, 1, 2)
    return wrappers + [operation] + wrappers


def test_and_and_erase():
    rng = random.Random(SEED + 1)
    errors = []
    checked = 0
    for cv1, cv2 in itertools.product((0, 1), repeat=2):
        compute = controlled_and_ops(cv1, cv2)
        erase = controlled_and_ops(cv1, cv2, erase=True)
        text_c = emit_qasm(3, compute)
        text_e = emit_qasm(3, erase)
        def image(i):
            return i | (int((i & 1) == cv1 and ((i >> 1) & 1) == cv2) << 2)
        for vector in [basis(4, i) for i in range(4)] + [random_state(4, rng) for _ in range(10)]:
            extended = vector + [0j] * 4
            computed = permuted(vector, image, output_size=8)
            errors.append(compare(dense_branch(text_c, extended), computed))
            for outcome in (0, 1):
                errors.append(compare(dense_branch(text_e, computed, (outcome,)),
                                      [x / math.sqrt(2) for x in extended]))
            checked += 3
        for i in range(4):
            assert run(3, compute, i) == (image(i), [])
            assert run(3, erase, image(i)) == (i, [0])
    try:
        run(3, [('and', 0, 1, 2)], 4)
    except SimulationError:
        pass
    else:
        raise AssertionError('Классический симулятор принял грязную цель И')
    return {'control_patterns': 4, 'vectors_and_branches': checked,
            'clean_target_only': True, 'max_error': max(errors)}


def test_ccx_arbitrary_target():
    ops = [('and', 0, 1, 3), ('cx', 3, 2), ('erase', 3, (('cz', 0, 1),))]
    text = emit_qasm(4, ops)
    rng = random.Random(SEED + 2)
    errors = []
    def image(i):
        return i ^ (((i & 1) & ((i >> 1) & 1)) << 2)
    vectors = [basis(8, i) for i in range(8)] + [random_state(8, rng) for _ in range(12)]
    for vector in vectors:
        expected = permuted(vector, image, output_size=16, factor=1 / math.sqrt(2))
        for outcome in (0, 1):
            errors.append(compare(dense_branch(text, vector + [0j] * 8, (outcome,)), expected))
    for i in range(8):
        assert run(4, ops, i) == (image(i), [0])
    return {'basis_inputs': 8, 'complex_superpositions': 12, 'branches_per_input': 2,
            'global_phase_compared': True, 'max_error': max(errors)}


def test_linear_correction_and_all_branches():
    # Два стирания: первое t=ab XOR l, второе u=ab. Перебираются все 4 ветви.
    ops = [('and', 0, 1, 3), ('cx', 2, 3),
           ('erase', 3, (('cz', 0, 1), ('z', 2))),
           ('and', 0, 1, 4), ('erase', 4, (('cz', 0, 1),))]
    text = emit_qasm(5, ops)
    rng = random.Random(SEED + 3)
    errors = []
    vectors = [basis(8, i) for i in range(8)] + [random_state(8, rng) for _ in range(10)]
    for vector in vectors:
        extended = vector + [0j] * 24
        for outcomes in itertools.product((0, 1), repeat=2):
            errors.append(compare(dense_branch(text, extended, outcomes), [x / 2 for x in extended]))
    for i in range(8):
        assert run(5, ops, i) == (i, [0, 0])
    # Удалённая CZ должна давать входозависимый коэффициент фазы, а не исчезнуть.
    broken = [('and', 0, 1, 3), ('cx', 3, 2), ('erase', 3, ())]
    assert run(4, broken, 0)[1] == [0]
    assert run(4, broken, 3)[1] == [1]
    vector = [1 / math.sqrt(2), 0j, 0j, 1j / math.sqrt(2)] + [0j] * 12
    good_text = emit_qasm(4, [('and', 0, 1, 3), ('cx', 3, 2), ('erase', 3, (('cz', 0, 1),))])
    bad_text = emit_qasm(4, broken)
    good = dense_branch(good_text, vector, (1,))
    bad = dense_branch(bad_text, vector, (1,))
    assert max(abs(a - b) for a, b in zip(good, bad)) > 0.9
    return {'vectors': len(vectors), 'outcome_strings': 4, 'missing_CZ_detected': True,
            'max_error': max(errors)}


def test_source_qualtran_and():
    """Compare the actual Qualtran And13 with the emitted And10 on the clean domain."""
    import cirq
    from qualtran.bloqs.mcmt import And
    qubit_order = [cirq.NamedQubit('target'), cirq.NamedQubit('ctrl_1'), cirq.NamedQubit('ctrl_0')]
    rng = random.Random(SEED + 4)
    errors = []
    checked = 0
    for cv1, cv2 in itertools.product((0, 1), repeat=2):
        circuit = And(cv1=cv1, cv2=cv2).to_clifford_t_circuit()
        matrix = circuit.unitary(qubit_order=qubit_order)
        text = emit_qasm(3, controlled_and_ops(cv1, cv2))
        for vector in [basis(4, i) for i in range(4)] + [random_state(4, rng) for _ in range(10)]:
            source_result = [sum(matrix[row, col] * vector[col] for col in range(4)) for row in range(8)]
            errors.append(compare(dense_branch(text, vector + [0j] * 4), source_result))
            checked += 1
    return {'source': 'qualtran.bloqs.mcmt.And.to_clifford_t_circuit',
            'control_patterns': 4, 'vectors': checked, 'matrix_columns_checked_per_pattern': 4,
            'domain': 'clean target |0> only; the full operator on a dirty target is not claimed',
            'global_phase_compared': True, 'max_error': max(errors)}
