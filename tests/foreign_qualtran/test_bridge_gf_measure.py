"""A2: local tomography, all branches and the native cross-check of GF2MulMBUC.

No test imports the qoptct provers or their parser/simulator.
"""
from __future__ import annotations
import pytest
pytest.importorskip("qualtran")
pytest.importorskip("cirq")

import argparse
import itertools
import json
from pathlib import Path
import time
import traceback

import numpy as np

from qoptct.foreign.qualtran.bridge.gf_measure import (
    GFMeasureError, basis_branch, emit_qasm, lower_mbuc, multiply_poly,
    parse_qasm, validate_ir, vector_branch,
)

ROOT = Path(__file__).resolve().parents[1]
TOL = 2e-10


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def close(a, b, message):
    error = float(np.max(np.abs(np.asarray(a) - np.asarray(b))))
    require(error < TOL, f'{message}: error={error}')
    return error


def tomography_states(dimension):
    """All bases and two relative phases for every pair of bases."""
    eye = np.eye(dimension, dtype=np.complex128)
    yield from eye
    for a in range(dimension):
        for b in range(a + 1, dimension):
            yield (eye[a] + eye[b]) / np.sqrt(2)
            yield (eye[a] + 1j * eye[b]) / np.sqrt(2)


def local_instruments():
    from qualtran import CBit, CtrlSpec
    from qualtran.bloqs.basic_gates import CZ, MeasureX

    # The native MeasureX gives the projectors as tensor[q_ket,q_bra,c].
    native_tensor = next(t for t in MeasureX().my_tensors(
        incoming={'q': 'input_q'}, outgoing={'c': 'output_c'}) if hasattr(t, 'data'))
    data = np.asarray(native_tensor.data)
    require(data.shape == (2, 2, 2), 'The shape of the source tensor of MeasureX changed')
    states, branches, max_error = 0, 0, 0.0
    for cv in (None, 0, 1):
        ops = [('measure_x', 2, 0)]
        if cv is not None:
            ops.append(('classical_cz', 0, 0, 1, cv))
            native_cz = CZ().controlled(CtrlSpec(qdtypes=CBit(), cvs=cv))
        ops.append(('discard_classical', 0))
        program = {'nq': 3, 'ops': ops}
        qasm = emit_qasm(program)
        for state in tomography_states(8):
            states += 1
            psi = state.reshape(2, 4)
            total_probability = 0.0
            reduced_channel = np.zeros((4, 4), dtype=np.complex128)
            native_channel = np.zeros_like(reduced_channel)
            for m in (0, 1):
                diagonal = np.ones(4, dtype=np.complex128)
                if cv is not None:
                    diagonal = np.asarray([
                        native_cz.basis_state_phase(ctrl=m, q1=a & 1, q2=(a >> 1) & 1) or 1
                        for a in range(4)], dtype=np.complex128)
                # Exact convention K_m=|0><m|H, followed by source diagonal.
                expected = diagonal * (psi[0] + (-1)**m * psi[1]) / np.sqrt(2)
                actual_full = vector_branch(qasm, state, (m,))
                max_error = max(max_error, close(actual_full[:4], expected, 'MeasureX/reset relative phase'))
                close(actual_full[4:], np.zeros(4), 'The measured wire must become |0>')
                source_rho = np.einsum('ta,ub,tu->ab', psi, psi.conj(), data[:, :, m])
                source_rho = diagonal[:, None] * source_rho * diagonal.conj()[None, :]
                target_rho = np.outer(actual_full[:4], actual_full[:4].conj())
                close(source_rho, target_rho, 'The native instrument MeasureX and the classical CZ')
                probability = float(np.vdot(actual_full, actual_full).real)
                close(probability, np.trace(source_rho).real, 'Probability of the branch')
                total_probability += probability
                reduced_channel += target_rho
                native_channel += source_rho
                branches += 1
            close(total_probability, 1.0, 'Normalisation after Discard classical outcome')
            close(reduced_channel, native_channel, 'Discard sums densities, not amplitudes')
    return {'status': 'PASS', 'states': states, 'branches': branches,
            'max_vector_error': max_error, 'source_tensor_checked': True,
            'scope': 'MeasureX/reset; classically controlled CZ with cv=0,1; classical Discard'}


def parity_rules():
    from qualtran.bloqs.gf_arithmetic.gf2_multiplication import Parity

    tables, branches = 0, 0
    for width in range(5):
        for values in itertools.product((0, 1), repeat=width):
            actual = Parity(width).on_classical_vals(x=np.asarray(values, dtype=np.uint8))
            require(np.array_equal(actual['x'], values), 'Parity changed the preserved classical inputs')
            parity = sum(values) % 2
            require(int(actual['parity']) == parity, 'The native Parity differs from XOR')
            tables += 1
            for cv in (0, 1):
                ops = [('measure_x', 2 + j, j) for j in range(width)]
                ops += [('parity', tuple(range(width)), width), ('classical_cz', width, 0, 1, cv),
                        ('discard_classical', width)]
                ops.extend(('discard_classical', j) for j in range(width))
                qasm = emit_qasm({'nq': width + 2, 'ops': ops})
                for state in tomography_states(4):
                    big = np.zeros(1 << (width + 2), dtype=np.complex128)
                    big[:4] = state
                    expected = state.copy() / np.sqrt(1 << width)
                    if parity == cv:
                        expected[3] *= -1
                    output = vector_branch(qasm, big, values)
                    close(output[:4], expected, 'Parity expression conditional phase')
                    close(output[4:], np.zeros(len(output) - 4), 'Parity reset work wires') if len(output) > 4 else None
                    branches += 1
    return {'status': 'PASS', 'truth_table_rows': tables, 'branches': branches,
            'widths': list(range(5)), 'superpositions': 'all real/imaginary pairs on both quantum controls'}


def malformed_rejections():
    rejected = []
    cases = {
        'unknown_quantum_gate': {'nq': 3, 'ops': [('t', 0)]},
        'read_before_definition': {'nq': 3, 'ops': [('classical_cz', 0, 0, 1, 1)]},
        'read_after_discard': {'nq': 3, 'ops': [('measure_x', 2, 0), ('discard_classical', 0), ('classical_cz', 0, 0, 1, 1)]},
        'source_qubit_reuse': {'nq': 3, 'ops': [('measure_x', 2, 0), ('classical_cz', 0, 2, 1, 1), ('discard_classical', 0)]},
        'classical_ssa_overwrite': {'nq': 3, 'ops': [('measure_x', 2, 0), ('parity', (0,), 0)]},
        'repeated_quantum_control': {'nq': 3, 'ops': [('measure_x', 2, 0), ('classical_cz', 0, 1, 1, 1), ('discard_classical', 0)]},
    }
    for label, program in cases.items():
        try:
            emit_qasm(program)
        except GFMeasureError:
            rejected.append(label)
        else:
            raise AssertionError('A wrong layout was accepted: ' + label)
    header = 'OPENQASM 3;\ninclude "stdgates.inc";\nqubit[3] q;\nbit[1] c;\n'
    for text in ('if (c[0]) {cz q[0], q[1];}', 'h q[-1];',
                 'c[0] = measure q[2];\nif (c[0] || 1) {x q[1];}', 'reset q[0];'):
        try:
            parse_qasm(header + text)
        except GFMeasureError:
            rejected.append(text)
        else:
            raise AssertionError('The parser accepted text outside the A2 language: ' + text)
    return {'status': 'PASS', 'count': len(rejected), 'cases': rejected}


def _encode(values, registers):
    return sum(((values[key] >> i) & 1) << q
               for key, wires in registers.items() for i, q in enumerate(wires))


def source_mbuc(n):
    from qualtran.bloqs.gf_arithmetic.gf2_multiplication import GF2MulMBUC
    from qualtran.simulation.classical_sim import PhasedClassicalSimState

    source = GF2MulMBUC(n)
    program = lower_mbuc(source)
    qasm = emit_qasm(program)
    ins, outs = program['layout']['inputs'], program['layout']['outputs']
    measure_entries = [entry for entry in program['trace'] if entry['rule'] == 'measure_x_reset_consumed_wire']
    require(len(measure_entries) == n, 'Exactly n source X-measurements')
    require([entry['binst_i'] for entry in program['trace']] == sorted(entry['binst_i'] for entry in program['trace']), 'The order of bb.add changed')
    outcomes_list = list(itertools.product((0, 1), repeat=n))
    native_phases = np.empty((len(outcomes_list), 1 << program['nq']), dtype=np.int8)
    native_output = np.empty(1 << program['nq'], dtype=np.int64)
    cases, raw_basis = 0, 0
    for x, y, product_register in itertools.product(range(1 << n), repeat=3):
        values = {'x': x, 'y': y, 'result': product_register}
        initial = _encode(values, ins)
        field_values = {key: source.qgf.gf_type(value) for key, value in values.items()}
        independent_product = multiply_poly(x, y, program['field_poly'], n)
        for branch, outcomes in enumerate(outcomes_list):
            fixed = {entry['binst_i']: value for entry, value in zip(measure_entries, outcomes)}
            native = PhasedClassicalSimState.from_cbloq(program['_source'], field_values, fixed_random_vals=fixed)
            result = native.simulate()
            require(set(result) == {'x', 'y'}, 'The measured result must not be a source output')
            actual_values = {key: int(value) for key, value in result.items()}
            require(actual_values == {'x': x, 'y': y}, 'The source block changed the operands')
            expected_bits = _encode(actual_values, outs)
            bits, sign, count = basis_branch(program, initial, outcomes)
            require(bits == expected_bits and count == n, 'The IR changed the basis action of the source')
            close(sign, native.phase, 'IR and the native phase of the source branch')
            if product_register == independent_product:
                close(sign, 1, 'A correct product must be removed without a dependent phase')
            native_phases[branch, initial] = sign
            native_output[initial] = expected_bits
            if n <= 3:
                state = np.zeros(1 << program['nq'], dtype=np.complex128)
                state[initial] = 1
                expected = np.zeros_like(state)
                expected[expected_bits] = sign / np.sqrt(1 << n)
                close(vector_branch(qasm, state, outcomes), expected, 'The actual QASM and the source basis branch')
                raw_basis += 1
            cases += 1

    rng = np.random.default_rng(20261121 + n)
    vector_cases = []
    for supported in (False, True):
        # Two independent complex superpositions: an arbitrary result and
        # the promised domain result=x*y. In the first case the probabilities of the branches
        # need not be equal at all.
        state = rng.normal(size=1 << program['nq']) + 1j * rng.normal(size=1 << program['nq'])
        if supported:
            valid = np.zeros_like(state)
            for x, y in itertools.product(range(1 << n), repeat=2):
                product = multiply_poly(x, y, program['field_poly'], n)
                initial = _encode({'x': x, 'y': y, 'result': product}, ins)
                valid[initial] = state[initial]
            state = valid
        state /= np.linalg.norm(state)
        total = 0.0
        for branch, outcomes in enumerate(outcomes_list):
            expected = np.zeros_like(state)
            np.add.at(expected, native_output, state * native_phases[branch] / np.sqrt(1 << n))
            actual = vector_branch(qasm, state, outcomes)
            close(actual, expected, 'The native linear superposition and QASM, including the phase')
            probability = float(np.vdot(actual, actual).real)
            if supported:
                close(probability, 1 / (1 << n), 'Probability on the promised domain result=x*y')
            total += probability
            vector_cases.append({'result_equals_product': supported, 'outcomes': outcomes, 'probability': probability})
        close(total, 1.0, 'Sum of the probabilities of all branches')

    return {'status': 'PASS', 'n': n, 'field_poly': program['field_poly'],
            'source_basis_branches': cases, 'actual_qasm_basis_branches': raw_basis,
            'superposition_branches': vector_cases, 'source_bloqs': len(program['trace'])}


def test_local_instruments():
    local_instruments()


def test_parity_rules():
    parity_rules()


def test_malformed_rejections():
    malformed_rejections()


@pytest.mark.parametrize('n', [1, 2])
def test_gf2mulmbuc(n):
    assert source_mbuc(n)['status'] == 'PASS'
