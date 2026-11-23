"""Public Qualtran adapters: mandatory old-path guard and explicit graph path.

Run only through lab/tools/one in a delegated compute slot. No arithmetic
prover is imported. This module does not modify common rules.py.
"""
import pytest
pytest.importorskip("qualtran")
pytest.importorskip("cirq")
import itertools
import json
from pathlib import Path
import resource
import time
import unittest
from unittest.mock import patch

import attrs
from qualtran import Bloq, CompositeBloq, CtrlSpec, QBit, QUInt, Register, Signature
from qualtran._infra.controlled import Controlled
from qualtran.bloqs.arithmetic import Add
from qualtran.bloqs.gf_arithmetic import GF2Square
from qualtran.bloqs.mcmt import And
from qualtran.bloqs.mcmt.controlled_via_and import ControlledViaAnd

from qoptct.foreign.qualtran.bridge.from_qualtran import lower, lower_graph
from qoptct.foreign.qualtran.bridge.qualtran_layout import ClassicalContract, LayoutError, expr, ref
from qoptct.foreign.qualtran.bridge.rules import translate
from qoptct.foreign.qualtran.bridge.sim import run
from test_bridge_layout import _check_contract_basis, _native_vectors

ROOT = Path(__file__).resolve().parents[1]


@attrs.frozen
class ExportOnly(Bloq):
    """Cirq-only source must not evade the newly mandatory layout audit."""
    @property
    def signature(self):
        return Signature([Register('x', QBit())])

    def as_cirq_op(self, qubit_manager, **quregs):
        import cirq
        return cirq.X(quregs['x'][0]), quregs


class IntegrationTests(unittest.TestCase):
    evidence = {}

    def test_public_lower_rejects_finding04_before_export(self):
        with patch.object(CompositeBloq, 'to_cirq_circuit_and_quregs',
                          side_effect=AssertionError('export must not occur before the guard')) as export:
            with self.assertRaisesRegex(LayoutError, 'THRU permutation'):
                lower(GF2Square(3, k=1))
            export.assert_not_called()
        self.evidence['finding04_public_lower'] = {
            'status': 'REJECTED', 'export_called': False, 'arithmetic_prover_used': False}

    def test_explicit_graph_square_matches_tensor(self):
        source = GF2Square(3, k=1)
        lowered = lower_graph(source)
        self.assertEqual(lowered['source_backend'], 'qualtran_bloqgraph')
        self.assertNotIn('source', lowered, 'do not falsely label a graph result as a Cirq circuit')
        self.assertTrue(lowered['thru_permutations'])
        self.evidence['square_graph_native_tensor'] = _native_vectors(source, lowered)

    def test_add_public_graph_and_cirq_have_distinct_sources(self):
        evidence = []
        for n in (2, 3):
            source = Add(QUInt(n))
            graph = lower_graph(source)
            self.assertEqual(graph['source_backend'], 'qualtran_bloqgraph')
            contract = ClassicalContract({'a': n, 'b': n}, {'a': n, 'b': n},
                {'a': ref('a'), 'b': expr('mod', expr('add', ref('a'), ref('b')), 1 << n)})
            checked_graph = _check_contract_basis(graph, contract)
            old = lower(source)
            self.assertEqual(old['source_backend'], 'qualtran_cirq_export')
            self.assertEqual(old['layout_guard']['status'], 'PASS')
            self.assertIn('source', old)
            compiled = translate(old, old['inputs'], old['outputs'])
            checked_old = 0
            for a, b in itertools.product(range(1 << n), repeat=2):
                values = {'a': a, 'b': b}
                initial = sum(((values[k] >> i) & 1) << q for k, word in old['inputs'].items()
                              for i, q in enumerate(word))
                actual, phases = run(compiled['nq'], compiled['ops'], initial)
                expected_values = {'a': a, 'b': (a + b) % (1 << n)}
                expected = sum(((expected_values[k] >> i) & 1) << q
                               for k, word in old['outputs'].items() for i, q in enumerate(word))
                self.assertEqual(actual, expected)
                self.assertFalse(any(phases))
                checked_old += 1
            evidence.append({'n': n, 'graph_full_basis': checked_graph, 'cirq_full_basis': checked_old,
                             'graph_native_tensor': _native_vectors(source, graph)})
        self.evidence['add'] = evidence

    def test_public_graph_controlled_and_interface(self):
        checks = 0
        for cls in (Controlled, ControlledViaAnd):
            for cv, cv1, cv2 in itertools.product((0, 1), repeat=3):
                for uncompute in (False, True):
                    source = cls(And(cv1, cv2, uncompute=uncompute), CtrlSpec(cvs=cv))
                    graph = lower_graph(source)
                    self.assertEqual(graph['source_backend'], 'qualtran_bloqgraph')
                    self.assertEqual(graph['nq'], 4)
                    self.assertEqual(len(graph['leaves']), 1)
                    leaf = graph['leaves'][0]
                    self.assertEqual(leaf['name'], 'controlled_erase' if uncompute else 'controlled_and')
                    self.assertEqual(leaf['controls'], [cv, cv1, cv2])
                    self.assertEqual(leaf['qubits'], [0, 1, 2, 3])
                    self.assertEqual(graph['initial_zero'], [] if uncompute else [3])
                    self.assertEqual(graph['final_zero'], [3] if uncompute else [])
                    checks += 1
        self.evidence['a1_public_graph_layout'] = {'cases': checks,
            'scope': 'API routing; atom phase/domain semantics are covered by A1 tests'}

    def test_cirq_only_source_fails_closed_before_export(self):
        with patch.object(ExportOnly, 'as_cirq_op',
                          side_effect=AssertionError('unaudited export must not run')) as export:
            with self.assertRaisesRegex(LayoutError, 'no supported atom or source BloqGraph'):
                lower(ExportOnly())
            export.assert_not_called()
        self.evidence['unknown_graph_public_lower'] = 'REJECTED before exporter'
