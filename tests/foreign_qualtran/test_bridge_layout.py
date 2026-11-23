"""A3 regressions. Run ONLY in an allocated slot through lab/tools/one.

No arithmetic prover is imported. Source tensor contraction, explicit integer
contracts, and an independent emitted-QASM branch simulator check the new wiring.
"""
import pytest
pytest.importorskip("qualtran")
pytest.importorskip("cirq")
import itertools
import json
import math
from pathlib import Path
import random
import resource
import time
import unittest

import attrs
import numpy as np
from qualtran import Bloq, BloqBuilder, QAny, QBit, QUInt, Register, Side, Signature
from qualtran.bloqs.basic_gates import CNOT, XGate
from qualtran.bloqs.bookkeeping import Cast, Partition
from qualtran.bloqs.mcmt import And

from qoptct.foreign.qualtran.bridge.qualtran_layout import (ClassicalContract, LayoutError, expr, flatten_registers,
                              lower_bloqgraph, montgomery_junk_contract, ref,
                              require_cirq_thru_identity, controlled_and_atom_rule)
from qoptct.foreign.qualtran.bridge.rules import translate
from qoptct.foreign.qualtran.bridge.sim import run
from test_bridge_rules import dense_branch, random_state

ROOT = Path(__file__).resolve().parents[1]


@attrs.frozen
class ShapedPermutation(Bloq):
    @property
    def signature(self):
        return Signature([Register('a', QUInt(2), shape=(2, 2))])

    def build_composite_bloq(self, bb, *, a):
        a = a.copy()
        bits = bb.split(a[0, 0])
        bits[0] = bb.add(XGate(), q=bits[0])
        bits[0], bits[1] = bb.add(CNOT(), ctrl=bits[0], target=bits[1])
        a[0, 0] = bb.join(bits, dtype=QUInt(2))
        # Source explicitly returns reordered soquets; no SWAP gate is inserted.
        return {'a': a[[1, 0], :]}


@attrs.frozen
class PartitionCast(Bloq):
    @property
    def signature(self):
        return Signature([Register('x', QUInt(4), side=Side.LEFT),
                          Register('a', QAny(2), shape=(2,), side=Side.RIGHT)])

    def build_composite_bloq(self, bb, *, x):
        x = bb.add(Cast(QUInt(4), QAny(4)), reg=x)
        a = bb.add(Partition(4, (Register('a', QAny(2), shape=(2,)),)), x=x)
        return {'a': a}


@attrs.frozen
class FreeAncilla(Bloq):
    clean: bool = True

    @property
    def signature(self):
        return Signature([Register('x', QBit())])

    def build_composite_bloq(self, bb, *, x):
        anc = bb.allocate(1)
        x, anc = bb.add(CNOT(), ctrl=x, target=anc)
        if self.clean:
            x, anc = bb.add(CNOT(), ctrl=x, target=anc)
        bb.free(anc)
        return {'x': x}


@attrs.frozen
class CopyJunk(Bloq):
    @property
    def signature(self):
        return Signature([Register('x', QUInt(2)),
                          Register('junk', QBit(), shape=(2,), side=Side.RIGHT)])

    def build_composite_bloq(self, bb, *, x):
        bits = bb.split(x)
        junk = bb.split(bb.allocate(2))
        for i in range(2):
            bits[i], junk[i] = bb.add(CNOT(), ctrl=bits[i], target=junk[i])
        return {'x': bb.join(bits, dtype=QUInt(2)), 'junk': junk}


@attrs.frozen
class AndRoundTrip(Bloq):
    cv1: int
    cv2: int

    @property
    def signature(self):
        return Signature([Register('ctrl', QBit(), shape=(2,))])

    def build_composite_bloq(self, bb, *, ctrl):
        ctrl, target = bb.add(And(self.cv1, self.cv2), ctrl=ctrl)
        ctrl = bb.add(And(self.cv1, self.cv2, uncompute=True), ctrl=ctrl, target=target)
        return {'ctrl': ctrl}


def _index(wires_msb, value):
    return sum(((value >> (len(wires_msb) - 1 - i)) & 1) << q
               for i, q in enumerate(wires_msb))


def _native_vectors(case, lowered):
    """Compare exact source tensor to every emitted-QASM measurement branch.

    No global phase quotient. Output tensor coordinates are embedded in the
    actual graph's returned physical registers; all nonoutputs are set to zero.
    """
    compiled = translate(lowered, lowered['inputs'], lowered['outputs'])
    in_msb = [q for reg in lowered['registers']['inputs'] for el in reg['elements'] for q in el['wires_msb']]
    out_msb = [q for reg in lowered['registers']['outputs'] for el in reg['elements'] for q in el['wires_msb']]
    tensor = np.asarray(case.tensor_contract(), dtype=np.complex128).reshape(1 << len(out_msb), 1 << len(in_msb))
    measurements = sum(op[0] == 'erase' for op in compiled['ops'])
    rng = random.Random(20261121)
    identity = np.eye(1 << len(in_msb), dtype=complex)
    vectors = [identity[:, i] for i in range(1 << len(in_msb))]
    vectors.extend(np.asarray(random_state(1 << len(in_msb), rng)) for _ in range(4))
    max_error, checks = 0.0, 0
    for logical in vectors:
        physical = np.zeros(1 << compiled['nq'], dtype=complex)
        for i, amplitude in enumerate(logical):
            physical[_index(in_msb, i)] = amplitude
        expected = np.zeros_like(physical)
        for i, amplitude in enumerate(tensor @ logical):
            expected[_index(out_msb, i)] = amplitude / math.sqrt(2 ** measurements)
        for outcomes in itertools.product((0, 1), repeat=measurements):
            actual = np.asarray(dense_branch(compiled['qasm'], physical, outcomes))
            error = float(np.max(np.abs(actual - expected)))
            if error > 1e-9:
                raise AssertionError(f'source tensor / physical QASM mismatch including phase: {error}')
            max_error = max(max_error, error)
            checks += 1
    return {'vectors_and_branches': checks, 'max_error_including_phase': max_error}


def _check_contract_basis(lowered, contract):
    contract.bind(lowered)
    compiled = translate(lowered, lowered['inputs'], lowered['outputs'])
    names = list(contract.input_widths)
    checks = 0
    for values in itertools.product(*(range(1 << contract.input_widths[k]) for k in names)):
        values = dict(zip(names, values))
        initial = sum(((values[k] >> i) & 1) << q for k, word in lowered['inputs'].items()
                      for i, q in enumerate(word))
        actual, phase = run(compiled['nq'], compiled['ops'], initial)
        expected_bits = contract.expected_physical_bits(lowered, values)
        expected = sum(bit << i for i, bit in enumerate(expected_bits))
        if actual != expected or any(phase):
            raise AssertionError(f'full physical contract failed: {values}, {actual} != {expected}, phases={phase}')
        checks += 1
    return checks


class LayoutTests(unittest.TestCase):
    evidence = {}

    def test_shape_axis_and_name_collision(self):
        reg = Register('a', QUInt(2), shape=(2, 2))
        words, metadata = flatten_registers([reg], {'a': np.arange(8).reshape(2, 2, 2)})
        self.assertEqual(words, {'a[0,0]': [1, 0], 'a[0,1]': [3, 2],
                                 'a[1,0]': [5, 4], 'a[1,1]': [7, 6]})
        self.assertEqual(metadata[0]['shape'], [2, 2])
        with self.assertRaises(LayoutError):
            flatten_registers([Register('a[0]', QBit())], {'a[0]': np.array([0])})
        with self.assertRaises(LayoutError):
            flatten_registers([reg], {'a': np.arange(8)})
        with self.assertRaises(LayoutError):
            flatten_registers([reg], {'a': np.zeros((2, 2, 2), dtype=int)})

    def test_shaped_source_order_and_permutation(self):
        case = ShapedPermutation()
        lowered = lower_bloqgraph(case)
        self.assertEqual([leaf['name'] for leaf in lowered['leaves']], ['x', 'cx'])
        self.assertEqual([leaf['qubits'] for leaf in lowered['leaves']], [[0], [0, 1]])
        self.assertEqual(lowered['outputs']['a[0,0]'], [5, 4])
        self.assertEqual(lowered['outputs']['a[1,0]'], [1, 0])
        changed = expr('xor', ref('a[0,0]'), 2)
        changed = expr('xor', changed, expr('shr', changed, 1))
        widths = {name: 2 for name in lowered['inputs']}
        contract = ClassicalContract(widths, widths.copy(),
            {'a[0,0]': ref('a[1,0]'), 'a[0,1]': ref('a[1,1]'),
             'a[1,0]': changed, 'a[1,1]': ref('a[0,1]')})
        self.evidence['shaped'] = {'basis_contracts': _check_contract_basis(lowered, contract),
                                  **_native_vectors(case, lowered)}
        with self.assertRaisesRegex(LayoutError, 'THRU permutation'):
            require_cirq_thru_identity(lowered)

    def test_split_join_cast_partition(self):
        case = PartitionCast()
        lowered = lower_bloqgraph(case)
        self.assertEqual(lowered['leaves'], [])
        contract = ClassicalContract({'x': 4}, {'a[0]': 2, 'a[1]': 2},
                                     {'a[0]': expr('shr', ref('x'), 2),
                                      'a[1]': expr('and', ref('x'), 3)})
        self.evidence['partition_cast'] = {'basis_contracts': _check_contract_basis(lowered, contract),
                                          **_native_vectors(case, lowered)}
        inverse = case.adjoint()
        reversed_layout = lower_bloqgraph(inverse)
        reversed_contract = ClassicalContract({'a[0]': 2, 'a[1]': 2}, {'x': 4},
            {'x': expr('add', expr('shl', ref('a[0]'), 2), ref('a[1]'))})
        self.evidence['partition_cast_inverse'] = {
            'basis_contracts': _check_contract_basis(reversed_layout, reversed_contract),
            **_native_vectors(inverse, reversed_layout)}

    def test_free_is_an_obligation_and_never_a_reset(self):
        contract = ClassicalContract({'x': 1}, {'x': 1}, {'x': ref('x')})
        good = lower_bloqgraph(FreeAncilla())
        self.assertEqual(good['retired_wires'], [1])
        self.assertEqual(good['final_zero'], [1])
        self.evidence['clean_free'] = {'basis_contracts': _check_contract_basis(good, contract),
                                       **_native_vectors(FreeAncilla(), good)}
        bad = lower_bloqgraph(FreeAncilla(clean=False))
        with self.assertRaisesRegex(AssertionError, 'full physical contract failed'):
            _check_contract_basis(bad, contract)
        missing = dict(good, final_zero=[])
        with self.assertRaisesRegex(LayoutError, 'clean-output obligation'):
            contract.bind(missing)

    def test_explicit_junk_total_contract(self):
        case = CopyJunk()
        lowered = lower_bloqgraph(case)
        contract = ClassicalContract({'x': 2}, {'x': 2, 'junk[0]': 1, 'junk[1]': 1},
            {'x': ref('x'), 'junk[0]': expr('shr', ref('x'), 1),
             'junk[1]': expr('and', ref('x'), 1)})
        self.evidence['explicit_junk'] = {'basis_contracts': _check_contract_basis(lowered, contract),
                                         **_native_vectors(case, lowered)}
        with self.assertRaisesRegex(LayoutError, 'arbitrary garbage'):
            ClassicalContract(contract.input_widths, contract.output_widths, {'x': ref('x')}).validate()
        with self.assertRaisesRegex(LayoutError, 'arbitrary garbage'):
            contract.bind(dict(lowered, garbage=lowered['outputs']['junk[0]']))
        with self.assertRaisesRegex(LayoutError, 'arbitrary garbage'):
            contract.bind(dict(lowered, free_wires=[0]))
        with self.assertRaises(LayoutError):
            ClassicalContract({'x': 2}, {'junk': 1}, {'junk': None}).validate()
        with self.assertRaises(LayoutError):
            ClassicalContract({'x': 2}, {'junk': 1}, {'junk': ref('unbound')}).validate()
        with self.assertRaisesRegex(LayoutError, 'overflows'):
            ClassicalContract({'x': 2}, {'junk': 1}, {'junk': ref('x')}).evaluate({'x': 3})
        wrong = ClassicalContract(contract.input_widths, contract.output_widths,
                                  {'x': ref('x'), 'junk[0]': 0, 'junk[1]': 0})
        with self.assertRaises(AssertionError):
            _check_contract_basis(lowered, wrong)

    def test_and_roundtrip_shaped_controls_all_branches(self):
        evidence = []
        for cv1, cv2 in itertools.product((0, 1), repeat=2):
            case = AndRoundTrip(cv1, cv2)
            lowered = lower_bloqgraph(case)
            self.assertEqual(lowered['inputs'], {'ctrl[0]': [0], 'ctrl[1]': [1]})
            contract = ClassicalContract({'ctrl[0]': 1, 'ctrl[1]': 1},
                {'ctrl[0]': 1, 'ctrl[1]': 1}, {'ctrl[0]': ref('ctrl[0]'), 'ctrl[1]': ref('ctrl[1]')})
            evidence.append({'cvs': [cv1, cv2], 'basis_contracts': _check_contract_basis(lowered, contract),
                             **_native_vectors(case, lowered)})
        self.evidence['and_shaped_all_branches'] = evidence

    def test_finding04_rejected_structurally_before_arithmetic(self):
        from qualtran.bloqs.gf_arithmetic import GF2Square
        case = GF2Square(3, k=1)
        lowered = lower_bloqgraph(case)
        synth = [p for p in lowered['thru_permutations'] if p['bloq_type'] == 'SynthesizeLRCircuit']
        self.assertTrue(synth, 'finding04 source did not expose the expected THRU permutation')
        with self.assertRaisesRegex(LayoutError, 'THRU permutation'):
            require_cirq_thru_identity(lowered)
        self.evidence['finding04'] = {'legacy_cirq_route': 'REJECTED', 'arithmetic_prover_used': False,
                                      'structural_evidence': synth,
                                      'separate_graph_backend': _native_vectors(case, lowered)}

    def test_controlled_and_callback_wiring(self):
        from qualtran import CtrlSpec
        from qualtran._infra.controlled import Controlled
        from qualtran.bloqs.mcmt.controlled_via_and import ControlledViaAnd
        checks = 0
        for cls in (Controlled, ControlledViaAnd):
            for cv, cv1, cv2 in itertools.product((0, 1), repeat=3):
                for uncompute in (False, True):
                    source = cls(And(cv1, cv2, uncompute=uncompute), CtrlSpec(cvs=cv))
                    lowered = lower_bloqgraph(source, extra_atom_rule=controlled_and_atom_rule)
                    self.assertEqual(lowered['nq'], 4)
                    self.assertEqual(len(lowered['leaves']), 1)
                    leaf = lowered['leaves'][0]
                    self.assertEqual(leaf['qubits'], [0, 1, 2, 3])
                    self.assertEqual(leaf['controls'], [cv, cv1, cv2])
                    self.assertEqual(leaf['name'], 'controlled_erase' if uncompute else 'controlled_and')
                    self.assertEqual(lowered['initial_zero'], [] if uncompute else [3])
                    self.assertEqual(lowered['final_zero'], [3] if uncompute else [])
                    self.assertEqual(lowered['retired_wires'], [3] if uncompute else [])
                    checks += 1
        self.evidence['a1_callback_layout_only'] = {
            'cases': checks, 'semantics_checked_elsewhere': 'A1 tests; this regression checks physical routing only'}

    def test_montgomery_junk_function_is_independent_and_complete(self):
        checks = 0
        for n, w, p in ((4, 1, 3), (4, 2, 7), (6, 2, 13)):
            contract = montgomery_junk_contract(n, w, p)
            for x, y in itertools.product(range(p), repeat=2):
                got = contract.evaluate({'x': x, 'y': y})
                self.assertEqual(got['x'], x)
                self.assertEqual(got['y'], y)
                self.assertEqual(got['target'], (x * y * pow(1 << n, -1, p)) % p)
                self.assertIn(got['reduced'], (0, 1))
                # A second recurrence derives each correction by divisibility,
                # without using the contract's modular inverse expression.
                t, h = 0, 0
                for shift in range(0, n, w):
                    u = t + ((x >> shift) & ((1 << w) - 1)) * y
                    m = u % (1 << w)
                    candidates = [k for k in range(1 << w) if (u + k * p) % (1 << w) == 0]
                    self.assertEqual(len(candidates), 1)
                    t = (u + candidates[0] * p) // (1 << w)
                    h = (h << w) | m
                self.assertEqual(got['qrom_indices'], h)
                self.assertEqual(got['reduced'], int(t >= p))
                checks += 1
            with self.assertRaisesRegex(LayoutError, 'precondition'):
                contract.evaluate({'x': p, 'y': 0})
            json.dumps(contract.to_dict())
        self.evidence['montgomery_spec_only'] = {'inputs': checks,
            'foreign_circuit_claim': 'NONE; this test checks the independent specification only'}

    def test_strict_rejections(self):
        from qualtran.bloqs.bookkeeping import Allocate, Free
        from qualtran.bloqs.basic_gates import ZGate
        with self.assertRaisesRegex(LayoutError, 'dirty Allocate'):
            lower_bloqgraph(Allocate(QBit(), dirty=True))
        with self.assertRaisesRegex(LayoutError, 'dirty Free'):
            lower_bloqgraph(Free(QBit(), dirty=True))
        with self.assertRaises(LayoutError):
            lower_bloqgraph(ZGate())
        with self.assertRaisesRegex(LayoutError, 'BloqGraph'):
            require_cirq_thru_identity({'thru_permutations': []})
        identity = lower_bloqgraph(XGate())
        self.assertEqual(require_cirq_thru_identity(identity)['status'], 'PASS')
        with self.assertRaisesRegex(LayoutError, 'too many leaves'):
            lower_bloqgraph(ShapedPermutation(), max_leaves=1)
