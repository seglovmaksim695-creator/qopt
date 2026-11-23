"""Opt-in Fredkin quantum rule tests. Execute only through lab/tools/one.

No whole CModMulK or arithmetic PROVED verdict is produced here. Numerical
quantum checks retain branch normalization and absolute phase. Independent
raw-QASM simulation and the existing exact domain guard are both required.
"""
import pytest
pytest.importorskip("qualtran")
pytest.importorskip("cirq")
from copy import deepcopy
from datetime import datetime, timezone
from hashlib import sha256
import argparse
import itertools
import json
import math
from pathlib import Path
import resource
import sys
import time
import traceback


ROOT = Path(__file__).resolve().parents[1]
TOLERANCE = 1e-10
GROUPS = ('native_quantum_protocols', 'atom_all_basis_and_tomography',
          'external_entanglement_and_reuse', 'strict_rejections_and_phase_canaries')


def need(condition, message):
    if not condition:
        raise AssertionError(message)


def _compile(source):
    from qoptct.qualtran.bridge.fredkin_rule import FredkinRules
    from qoptct.qualtran.bridge.qualtran_layout import lower_bloqgraph
    from qoptct.qualtran.bridge.rules import translate
    book = FredkinRules()
    low = lower_bloqgraph(source, extra_atom_rule=book.chain_with_basis_and_a1)
    prog = translate(low, low['inputs'], low['outputs'])
    return low, prog, book, book.bind(low, prog)


def _wire_order(low, side):
    return [q for reg in low['registers'][side] for element in reg['elements']
            for q in element['wires_msb']]


def _index(wires_msb, logical):
    return sum(((logical >> (len(wires_msb) - 1 - i)) & 1) << q for i, q in enumerate(wires_msb))


def _tomography(dimension):
    import numpy as np
    for i in range(dimension):
        v = np.zeros(dimension, complex)
        v[i] = 1
        yield f'basis:{i}', v
    # Both real and imaginary coherences of every pair, including controls.
    for i in range(dimension):
        for j in range(i + 1, dimension):
            for phase, label in ((1, 'real'), (1j, 'imaginary')):
                v = np.zeros(dimension, complex)
                v[i], v[j] = 1 / math.sqrt(2), phase / math.sqrt(2)
                yield f'pair:{i}:{j}:{label}', v


def _error(got, expected):
    import numpy as np
    got, expected = np.asarray(got), np.asarray(expected)
    need(got.shape == expected.shape, 'Quantum vector/matrix shape differs')
    value = float(np.max(np.abs(got - expected)))
    need(math.isfinite(value) and value <= TOLERANCE, f'Absolute complex amplitude error {value}')
    return value


def test_native_quantum_protocols():
    import cirq
    import numpy as np
    from qualtran import DecomposeTypeError
    from qualtran.bloqs.basic_gates import TwoBitCSwap
    source = TwoBitCSwap()
    need(tuple(r.name for r in source.signature) == ('ctrl', 'x', 'y'), 'Native register order differs')
    tensor = np.asarray(source.tensor_contract(), dtype=complex)
    need(tensor.shape == (8, 8), 'Native tensor is not an 8x8 quantum operator')
    native_circuit = source.to_clifford_t_circuit()
    unitary = native_circuit.unitary(qubit_order=[cirq.NamedQubit(n) for n in ('ctrl', 'x', 'y')])
    expected = np.zeros((8, 8), dtype=complex)
    for i in range(8):
        c, x, y = (i >> 2) & 1, (i >> 1) & 1, i & 1
        j = (c << 2) | ((y if c else x) << 1) | (x if c else y)
        expected[j, i] = 1
        need(source.call_classically(ctrl=c, x=x, y=y) == (c, y if c else x, x if c else y),
             'Native classical corroboration differs')
    errors = [_error(tensor, expected), _error(unitary, tensor),
              _error(tensor.conj().T @ tensor, np.eye(8))]
    try:
        source.decompose_bloq()
    except DecomposeTypeError as exc:
        atomic = dict(type=type(exc).__name__, message=str(exc))
    else:
        raise AssertionError('Native atomic-decomposition assumption changed')
    return dict(basis_columns=8, native_tensor=True, native_clifford_t_unitary=True,
                actual_signature_order=['ctrl', 'x', 'y'], native_atomic_diagnostic=atomic,
                global_phase_quotient_used=False, max_error=max(errors))


def _quantum_check(source, low, prog, binding, *, expected_atoms):
    import numpy as np
    from qoptct.qualtran.bridge.domain_guard import check_domains
    from qoptct.qualtran.bridge.sim import run
    from qoptct.qualtran.bridge.sim_parse import parse_qasm
    from test_bridge_rules import dense_branch
    ins, outs = _wire_order(low, 'inputs'), _wire_order(low, 'outputs')
    nq, ops = parse_qasm(prog['qasm'])
    need(nq == prog['nq'] and ops == prog['ops'], 'Independent whole raw parse differs')
    need(nq <= 8 and len(ins) <= 5 and len(outs) == len(ins), 'Dense tensor allocation cap exceeded')
    need(len(binding['instances']) == expected_atoms and binding['semantic_status'] == 'UNPROVED',
         'Rule binding count or structural-only qualification differs')
    measurements = sum(op[0] == 'erase' for op in ops)
    need(measurements == expected_atoms <= 2, 'Unexpected measurement count')
    tensor = np.asarray(source.tensor_contract(), dtype=complex).reshape(1 << len(outs), 1 << len(ins))
    contract = dict(rule=binding['rule'], binding=binding, requires_exact_partial_domains=True,
                    inputs='all physical source input states; scratch initially zero',
                    output='native source matrix; all physical nonoutputs zero')
    domains = check_domains(prog['qasm'], prog['layout'], contract=contract,
                            budget=30, max_nodes=100000, max_mb=300)
    need(domains['status'] == 'PASS', 'Exact partial AND/ERASE domain guard failed: ' + str(domains))
    outputs_mask = sum(1 << q for q in outs)
    size = 1 << nq
    max_error, tested, branches = 0.0, 0, 0
    for label, logical in _tomography(1 << len(ins)):
        physical, wanted = np.zeros(size, complex), np.zeros(size, complex)
        for i, value in enumerate(logical):
            physical[_index(ins, i)] = value
        for i, value in enumerate(tensor @ logical):
            wanted[_index(outs, i)] = value
        probabilities = []
        for outcomes in itertools.product((0, 1), repeat=measurements):
            actual = np.asarray(dense_branch(prog['qasm'], physical, outcomes), complex)
            desired = wanted / math.sqrt(1 << measurements)
            max_error = max(max_error, _error(actual, desired))
            need(not any(abs(a) > TOLERANCE for i, a in enumerate(actual) if i & ~outputs_mask),
                 'A physical working wire is not finally zero')
            probabilities.append(float(np.vdot(actual, actual).real))
            branches += 1
        need(abs(sum(probabilities) - 1) <= TOLERANCE, 'Measurement branches do not preserve total norm')
        if label.startswith('basis:'):
            index = int(label.split(':')[1])
            final, phase = run(nq, ops, _index(ins, index))
            native_column = tensor[:, index]
            nonzero = np.flatnonzero(np.abs(native_column) > TOLERANCE)
            need(len(nonzero) == 1 and abs(native_column[int(nonzero[0])] - 1) <= TOLERANCE,
                 'Native source is not the claimed phase-free basis permutation')
            need(final == _index(outs, int(nonzero[0])) and not any(phase),
                 'Independent raw all-wire result or phase differs')
        tested += 1
    return dict(vectors=tested, basis_vectors=1 << len(ins), pair_phases=['1', 'i'],
                all_pairs=True, branch_vectors=branches, outcomes_per_vector=1 << measurements,
                expected_branch_amplitude=1 / math.sqrt(1 << measurements),
                all_working_wires_initial_final_zero=True, domain_guard=domains,
                absolute_phase_compared=True, max_error=max_error)


def test_atom_all_basis_and_tomography():
    from qoptct.qualtran.bridge.rules import translate
    from qualtran.bloqs.basic_gates import TwoBitCSwap
    source = TwoBitCSwap()
    low, prog, _, bound = _compile(source)
    ordinary = _quantum_check(source, low, prog, bound, expected_atoms=1)
    # A second nontrivial physical layout, with two unused clean wires.
    # Invoke the exact source atom callback on that layout; no source gates
    # are altered and ctrl/x/y keep the native register ordering.
    import numpy as np
    from qoptct.qualtran.bridge.fredkin_rule import FredkinRules
    from qoptct.qualtran.bridge.qualtran_layout import flatten_registers
    book = FredkinRules()
    inputs = {name: np.asarray([q], dtype=object) for name, q in zip(('ctrl', 'x', 'y'), (4, 0, 2))}
    def forbidden_allocate(*args):
        raise AssertionError('Fredkin must not allocate source wires')
    expansion = book(source, inputs, forbidden_allocate, 'TwoBitCSwap')
    words, desc = flatten_registers(source.signature, inputs)
    perm = deepcopy(low)
    perm.update(nq=5, inputs=words, outputs=deepcopy(words), leaves=expansion.leaves,
                wire_names=['test physical q%d' % i for i in range(5)],
                registers=dict(inputs=desc, outputs=deepcopy(desc)), initial_zero=[1, 3], final_zero=[1, 3],
                graph_events=[dict(path='TwoBitCSwap', bloq_type='TwoBitCSwap',
                                   inputs_msb=deepcopy(words), outputs_msb=deepcopy(words))])
    pp = translate(perm, perm['inputs'], perm['outputs'])
    pb = book.bind(perm, pp)
    alternate = _quantum_check(source, perm, pp, pb, expected_atoms=1)
    return dict(layouts=[ordinary, alternate], permuted_wires={'ctrl': 4, 'x': 0, 'y': 2},
                max_error=max(ordinary['max_error'], alternate['max_error']))


def _external_source(repeat):
    import attrs
    from qualtran import Bloq, Register, QBit, Signature
    from qualtran.bloqs.basic_gates import TwoBitCSwap
    @attrs.frozen
    class WithExternal(Bloq):
        repetitions: int

        @property
        def signature(self):
            return Signature([Register(name, QBit()) for name in ('external', 'ctrl', 'x', 'y')])

        def build_composite_bloq(self, bb, *, external, ctrl, x, y):
            for _ in range(self.repetitions):
                ctrl, x, y = bb.add(TwoBitCSwap(), ctrl=ctrl, x=x, y=y)
            return dict(external=external, ctrl=ctrl, x=x, y=y)
    return WithExternal(repeat)


def test_external_entanglement_and_reuse():
    rows = []
    for repeat in (1, 2):
        source = _external_source(repeat)
        low, prog, _, bound = _compile(source)
        result = _quantum_check(source, low, prog, bound, expected_atoms=repeat)
        need(len(prog['layout']['bridge_scratch']) == 1, 'Existing translator scratch policy changed')
        result.update(repeated_atoms=repeat, external_qubit=True,
                      all_pairs_include_control_external_entanglement=True,
                      scratch_reused=repeat == 2)
        rows.append(result)
    return dict(cases=rows, max_error=max(row['max_error'] for row in rows))


def test_strict_rejections_and_phase_canaries():
    import numpy as np
    from qualtran.bloqs.basic_gates import TwoBitCSwap, XGate
    from qoptct.qualtran.bridge.fredkin_rule import FredkinRules, FredkinError, match
    from qoptct.qualtran.bridge.emit import emit_qasm
    from qoptct.qualtran.bridge.domain_guard import check_domains
    from test_bridge_rules import dense_branch
    source = TwoBitCSwap()
    inputs = {name: np.asarray([q], object) for q, name in enumerate(('ctrl', 'x', 'y'))}
    class Forged:
        name = 'TwoBitCSwap'
    class Subclass(TwoBitCSwap):
        pass
    need(match(Forged(), inputs, 'forged') is None and match(Subclass(), inputs, 'subclass') is None
         and match(XGate(), inputs, 'other') is None, 'Matcher trusts a name, subclass or unknown atom')
    rejected = []
    for label, changed in (('alias', dict(inputs, y=np.asarray([0], object))),
                           ('bool', dict(inputs, ctrl=np.asarray([True], object))),
                           ('negative', dict(inputs, x=np.asarray([-1], object))),
                           ('shape', dict(inputs, y=np.asarray([[2]], object))),
                           ('missing', {k: v for k, v in inputs.items() if k != 'y'})):
        try:
            match(source, changed, label)
        except FredkinError:
            rejected.append(label)
        else:
            raise AssertionError('Malformed physical mapping accepted: ' + label)
    low, prog, book, bound = _compile(source)
    try:
        book.bind(low, prog)
    except FredkinError:
        rejected.append('book_reuse')
    else:
        raise AssertionError('Rule book binding was reusable')
    # Both a gate mutation and a genuine phase-only error must be detected.
    # The latter keeps every classical output but loses the actual CZ.
    missing_cz = deepcopy(prog['ops'])
    position = next(i for i, op in enumerate(missing_cz) if op[0] == 'erase')
    missing_cz[position] = ('erase', missing_cz[position][1], ())
    wrong_direction = deepcopy(prog['ops'])
    wrong_direction[0] = ('cx', wrong_direction[0][2], wrong_direction[0][1])
    witnesses = []
    for label, changed in (('missing_actual_CZ', missing_cz), ('reversed_first_CX', wrong_direction)):
        text = emit_qasm(prog['nq'], changed)
        found = None
        for vector_label, v in _tomography(8):
            if not vector_label.endswith('imaginary'):
                continue
            physical = np.zeros(1 << prog['nq'], complex)
            for i, amplitude in enumerate(v):
                physical[_index(_wire_order(low, 'inputs'), i)] = amplitude
            actual = np.asarray(dense_branch(text, physical, (1,)), complex)
            reference = np.asarray(dense_branch(prog['qasm'], physical, (1,)), complex)
            if float(np.max(np.abs(actual - reference))) > 0.1:
                found = dict(mutation=label, vector=vector_label, outcome=[1],
                             max_error=float(np.max(np.abs(actual - reference))))
                break
        need(found is not None, 'Complex branch canary was not detected: ' + label)
        witnesses.append(found)
    guard = check_domains(emit_qasm(prog['nq'], missing_cz), prog['layout'],
                          contract=dict(rule='missing correction negative control'),
                          budget=30, max_nodes=100000, max_mb=300)
    need(guard['status'] == 'REFUTED', 'Exact partial-domain guard accepted missing CZ')
    # Malformed same-sized raw/trace/source map cannot pass structural bind.
    tamper_checks = []
    for label in ('source_map', 'trace', 'leaf', 'raw', 'scratch_clean'):
        original_low, original_prog, fresh, _ = _compile(source)
        # _compile bound that book already; use a new independently populated
        # callback book, while retaining the same source/lowered program.
        fresh = FredkinRules()
        fresh(source, inputs, lambda *args: {}, 'TwoBitCSwap')
        if label == 'source_map':
            original_low['graph_events'][0]['inputs_msb']['x'] = [2]
        elif label == 'trace':
            original_prog['trace'][1]['logical_start'] += 1
        elif label == 'leaf':
            original_low['leaves'][0]['qubits'].reverse()
        elif label == 'raw':
            original_prog['qasm'] = emit_qasm(original_prog['nq'], wrong_direction)
        else:
            original_prog['layout']['initial_zero'].remove(original_prog['layout']['bridge_scratch'][0])
        try:
            fresh.bind(original_low, original_prog)
        except (ValueError, KeyError, IndexError):
            tamper_checks.append(label)
        else:
            raise AssertionError('Changed structural source binding accepted: ' + label)
    return dict(rejected_inputs=rejected, tamper_rejected=tamper_checks,
                quantum_canary_witnesses=witnesses, missing_correction_domain_guard=guard,
                arithmetic_prover_invoked=False)
