"""Local source-vector rules and CModAddK context; run only via one.

Positive checks are source-tensor/raw-QASM rule tests, not new arithmetic
PROVED claims. A deliberately wrong IntEffect reaches the unchanged prover
as a negative control to check the complete physical zero-output obligation.
"""
import pytest
pytest.importorskip("qualtran")
pytest.importorskip("cirq")

from copy import deepcopy
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import version
import itertools
import json
from pathlib import Path
import random
import resource
import sys
import time
import traceback

import attrs
import numpy as np
from qualtran import Bloq, QBit, QUInt, Register, Side, Signature
from qualtran.bloqs.basic_gates import IntState, IntEffect, ZeroState, OneState, ZeroEffect, OneEffect

from qoptct.qualtran.bridge.basis_vectors import BasisVectorError, BasisVectorRules
from qoptct.qualtran.bridge.domain_guard import check_domains
from qoptct.qualtran.bridge.emit import emit_qasm
from qoptct.qualtran.bridge.qualtran_layout import lower_bloqgraph
from qoptct.qualtran.bridge.rules import translate
from qoptct.qualtran.bridge.sim import run
from qoptct.qualtran.bridge.sim_parse import parse_qasm
from qoptct.qualtran.bridge.strict_proof import digest_json, prove_program
from test_bridge_rules import compare, dense_branch, random_state




@attrs.frozen
class WithExternal(Bloq):
    bit: int
    state: bool

    @property
    def signature(self):
        return Signature([Register("e", QUInt(2)),
                          Register("q", QBit(), side=Side.RIGHT if self.state else Side.LEFT)])

    def build_composite_bloq(self, bb, *, e, **kwargs):
        cls = (OneState if self.bit else ZeroState) if self.state else (OneEffect if self.bit else ZeroEffect)
        result = bb.add_d(cls(), **({} if self.state else {"q": kwargs["q"]}))
        return dict(e=e, **result)


@attrs.frozen
class IntRoundTrip(Bloq):
    n: int
    state_value: int
    effect_value: int
    repeat: int = 1

    @property
    def signature(self):
        return Signature([Register("e", QUInt(2))])

    def build_composite_bloq(self, bb, *, e):
        for _ in range(self.repeat):
            value = bb.add(IntState(self.state_value, self.n))
            bb.add_d(IntEffect(self.effect_value, self.n), val=value)
        return dict(e=e)


def _compile(source, *, bind=True):
    book = BasisVectorRules()
    lowered = lower_bloqgraph(source, extra_atom_rule=book.chain_with_a1)
    compiled = translate(lowered, lowered["inputs"], lowered["outputs"])
    binding = book.bind(lowered, compiled) if bind else None
    return lowered, compiled, book, binding


def _msb_wires(lowered, side):
    return [q for reg in lowered["registers"][side] for element in reg["elements"]
            for q in element["wires_msb"]]


def _embed_index(wires_msb, logical):
    return sum(((logical >> (len(wires_msb) - 1 - j)) & 1) << q for j, q in enumerate(wires_msb))


def _encode(words, values):
    return sum(((values[name] >> j) & 1) << q for name, wires in words.items() for j, q in enumerate(wires))


def _tensor_vectors(source, lowered, compiled, allowed):
    ins, outs = _msb_wires(lowered, "inputs"), _msb_wires(lowered, "outputs")
    if compiled["nq"] > 10 or len(ins) > 8 or len(outs) > 8:
        raise ValueError("Small native tensor/dense allocation limit exceeded")
    tensor = np.asarray(source.tensor_contract(), dtype=complex).reshape(1 << len(outs), 1 << len(ins))
    weights = [[complex(i == j) for i in range(len(allowed))] for j in range(len(allowed))]
    rng = random.Random(20261121)
    weights.extend(random_state(len(allowed), rng) for _ in range(4))
    # A Bell state with a non-real relative phase in the external two-qubit
    # register. The promised effect qubit itself is fixed, as its domain requires.
    bell = [0j] * len(allowed)
    bell[0], bell[-1] = 1 / np.sqrt(2), 1j / np.sqrt(2)
    weights.append(bell)
    size = 1 << compiled["nq"]
    errors = []
    for coefficients in weights:
        logical = np.zeros(1 << len(ins), dtype=complex)
        physical = np.zeros(size, dtype=complex)
        for index, coefficient in zip(allowed, coefficients):
            logical[index] = coefficient
            physical[_embed_index(ins, index)] = coefficient
        desired = np.zeros(size, dtype=complex)
        for index, amplitude in enumerate(tensor @ logical):
            desired[_embed_index(outs, index)] = amplitude
        errors.append(compare(dense_branch(compiled["qasm"], physical, ()), desired))
    return dict(vectors=len(weights), basis_vectors=len(allowed), random_complex_vectors=4,
                external_entangled_vector=True, global_phase_compared=True, max_error=max(errors)), tensor


def test_four_atoms():
    results = []
    for state, bit in itertools.product((False, True), (0, 1)):
        source = WithExternal(bit, state)
        low, prog, book, binding = _compile(source)
        allowed = list(range(4)) if state else [index for index in range(8) if index & 1 == bit]
        evidence, tensor = _tensor_vectors(source, low, prog, allowed)
        nq, ops = parse_qasm(prog["qasm"])
        contract = dict(rule="basis vectors", bit=bit, state=state, binding=binding,
                        post="external e unchanged; source outputs exact; all other wires zero")
        domain_only = check_domains(prog["qasm"], prog["layout"], contract=contract)
        if domain_only["status"] != "PASS" or binding["semantic_status"] != "UNPROVED":
            raise AssertionError("Unexpected domain-only or structural qualification")
        invalid = 0
        if not state:
            output_mask = sum(1 << q for wires in prog["layout"]["outputs"].values() for q in wires)
            in_msb = _msb_wires(low, "inputs")
            for index in range(8):
                physical = _embed_index(in_msb, index)
                final, phase = run(nq, ops, physical)
                admissible = index & 1 == bit
                if bool(final & ~output_mask) == admissible or phase:
                    raise AssertionError("Full-zero output does not exactly separate the effect domain")
                if not admissible:
                    if np.any(tensor[:, index]):
                        raise AssertionError("Native effect does not annihilate the excluded basis input")
                    invalid += 1
            # Explicit entanglement between target and external register leaves
            # the domain. Native projection loses norm; physical X/free model
            # retains a nonzero consumed wire. It is therefore NOT admitted.
            mixed = np.zeros(8, dtype=complex)
            mixed[0], mixed[7] = 1 / np.sqrt(2), 1j / np.sqrt(2)
            projected_norm = float(np.linalg.norm(tensor @ mixed))
            if not projected_norm < 1 - 1e-10:
                raise AssertionError("Out-of-domain coherent input was not projected by native effect")
            evidence["out_of_domain_control_entangled_input_rejected"] = True
        results.append(dict(state=state, bit=bit, status="PASS", tensors=evidence,
                            invalid_basis_inputs_rejected=invalid, binding=binding,
                            and_erase_domain_guard=domain_only,
                            domain_guard_does_not_certify_effects=True))
    return dict(cases=results, status="PASS")


def test_int_roundtrips_and_structure():
    rows = []
    for n in range(1, 5):
        for value in range(1 << n):
            source = IntRoundTrip(n, value, value)
            low, prog, _, binding = _compile(source)
            tensors, _ = _tensor_vectors(source, low, prog, list(range(4)))
            kinds = [row["kind"] for row in binding["obligations"]]
            if kinds.count("state") != n or kinds.count("effect") != n:
                raise AssertionError("IntState/Effect did not use their actual bit decomposition")
            rows.append(dict(n=n, value=value, status="PASS", tensors=tensors))
    source = IntRoundTrip(2, 3, 3, repeat=2)
    low, prog, _, binding = _compile(source)
    allocated = [row["wire"] for row in binding["obligations"] if row["kind"] == "state"]
    if len(allocated) != 4 or len(set(allocated)) != 4:
        raise AssertionError("Consumed source wires were recycled")
    for defect in ("remove_final_zero", "remove_initial_zero", "later_physical_use"):
        low, prog, book, _ = _compile(source, bind=False)
        q = book.obligations[0]["wire"]
        if defect == "remove_final_zero":
            low["final_zero"].remove(q)
        elif defect == "remove_initial_zero":
            prog["layout"]["initial_zero"].remove(q)
        else:
            prog["ops"].append(("x", q))
            prog["qasm"] = emit_qasm(prog["nq"], prog["ops"])
        try:
            book.bind(low, prog)
        except BasisVectorError:
            pass
        else:
            raise AssertionError(f"Structural binding accepted {defect}")
    return dict(status="PASS", roundtrips=rows, distinct_reallocated_wires=allocated,
                structural_negative_cases=3)


def test_wrong_effect_full_contract():
    correct = IntRoundTrip(3, 5, 5)
    wrong = IntRoundTrip(3, 5, 4)
    _, baseline, _, _ = _compile(correct)
    low, prog, _, binding = _compile(wrong)
    raw_nq, raw_ops = parse_qasm(prog["qasm"])
    final, phases = run(raw_nq, raw_ops, 0)
    allowed = sum(1 << q for wires in prog["layout"]["outputs"].values() for q in wires)
    if not final & ~allowed or any(phases):
        raise AssertionError("Wrong constant effect did not leave a concrete dirty wire")
    # Native tensor of mismatched state/effect is the zero operator, not identity.
    native = np.asarray(wrong.tensor_contract(), dtype=complex)
    if np.any(native):
        raise AssertionError("Native mismatched IntEffect tensor is not zero")
    contract = dict(rule="mismatched source IntEffect negative control", expected="e'=e; all other wires zero",
                    obligations=binding, requires_exact_partial_domains=True)
    domain_only = check_domains(prog["qasm"], prog["layout"], contract=contract)
    if domain_only["status"] != "PASS":
        raise AssertionError("Expected AND/ERASE guard to expose its narrower scope")
    witness = dict(reason="wrong source effect leaves consumed physical wire nonzero",
                   inputs={"e": 0}, dirty_bits=final & ~allowed, checked_inputs=1)
    receipt = dict(status="WITNESS", independent_inputs=1, all_wires_checked=True, phase_checked=True,
                   qasm_sha256=sha256(prog["qasm"].encode()).hexdigest(),
                   baseline_qasm_sha256=sha256(baseline["qasm"].encode()).hexdigest(),
                   layout_sha256=digest_json(prog["layout"]), contract_sha256=digest_json(contract), witness=witness)
    symbolic = lambda S: S.out(prog["layout"]["outputs"]["e"], S.reg(prog["layout"]["inputs"]["e"]))
    verdict = prove_program(prog["qasm"], prog["layout"], symbolic, validation=receipt,
                            contract=contract, negative_control=True, budget=30, max_mb=500)
    if verdict["status"] == "PROVED":
        raise AssertionError("HOLE: false PROVED for a wrong constant effect")
    if verdict["status"] != "REFUTED" or verdict.get("prover_invoked") is not True:
        raise AssertionError(f"Wrong effect did not receive actual full-spec REFUTED: {verdict}")
    return dict(status="PASS", wrong_effect_binding=binding, and_erase_guard=domain_only,
                raw_witness=witness, full_physical_spec=verdict)
