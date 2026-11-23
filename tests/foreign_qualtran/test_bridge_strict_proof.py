"""Strict proof: the full output contract, the promise and the canary mutants go through the real prover.

Ported from the 0.6 bridge round (test_strict_proof): a Qualtran Add(QUInt(4)) is lowered, proved, and
every targeted mutant of the program must be NOT PROVED; incomplete or tampered contracts are refused.
"""
import copy
import itertools
from hashlib import sha256

import pytest

pytest.importorskip("qualtran")
pytest.importorskip("cirq")

from qualtran import QUInt
from qualtran.bloqs.arithmetic import Add

from qoptct.foreign.qualtran.bridge.domain_guard import check_domains
from qoptct.foreign.qualtran.bridge.emit import emit_qasm
from qoptct.foreign.qualtran.bridge.from_qualtran import lower_graph
from qoptct.foreign.qualtran.bridge.rules import translate
from qoptct.foreign.qualtran.bridge.strict_proof import digest_json, prove_program
from qoptct.foreign.qualtran.bridge.targeted_canaries import generate
from qoptct.verify.bdd import ZERO

N = 4


@pytest.fixture(scope="module")
def setup():
    low = lower_graph(Add(QUInt(N)))
    prog = translate(low, low["inputs"], low["outputs"])
    layout = prog["layout"]
    samples = [dict(a=a, b=b) for a, b in itertools.product(range(1 << N), repeat=2)]
    contract = dict(source="Qualtran Add", n=N, formula="a_out=a; b_out=(a+b) mod 2^n; all other wires=0",
                    requires_exact_partial_domains=True)
    receipt = dict(status="PASS", qasm_sha256=sha256(prog["qasm"].encode()).hexdigest(),
                   contract_sha256=digest_json(contract), independent_inputs=len(samples),
                   layout_sha256=digest_json(layout), all_wires_checked=True, phase_checked=True)

    def spec(S):
        a = S.reg(layout["inputs"]["a"])
        b = S.reg(layout["inputs"]["b"])
        out = dict(S.out(layout["outputs"]["a"], a))
        out.update(S.out(layout["outputs"]["b"], a + b))
        return out

    return prog, layout, samples, contract, receipt, spec


def test_baseline_is_proved_and_mutants_are_not(setup):
    prog, layout, samples, contract, receipt, spec = setup
    domains = check_domains(prog["qasm"], layout, contract=contract)
    assert domains["status"] == "PASS", domains
    baseline = prove_program(prog["qasm"], layout, spec, validation=receipt, contract=contract, domains=domains)
    assert baseline["status"] == "PROVED", baseline
    assert "/" not in baseline["prover_file"][:1], "no absolute paths in results"

    expected = lambda v: dict(a=v["a"], b=(v["a"] + v["b"]) % (1 << N))
    rare = [dict(label="zero", match={"a": 0}, witness={"a": 0, "b": 0}),
            dict(label="carry", match={"a": 15, "b": 1}, witness={"a": 15, "b": 1}),
            dict(label="maximum", match={"a": 15}, witness={"a": 15, "b": 0})]
    mutants = generate(prog, expected, lambda v: True, samples, rare_cases=rare)
    assert len(mutants) >= 10
    for m in mutants:
        qasm = emit_qasm(m["nq"], m["ops"])
        changed = copy.deepcopy(layout)
        changed["nq"] = m["nq"]
        changed["initial_zero"] += list(range(prog["nq"], m["nq"]))
        v = dict(receipt, status="WITNESS", qasm_sha256=m["qasm_sha256"], layout_sha256=digest_json(changed),
                 baseline_qasm_sha256=m["baseline_qasm_sha256"], witness=m["witness"])
        verdict = prove_program(qasm, changed, spec, validation=v, contract=contract, negative_control=True)
        assert verdict["prover_invoked"], verdict
        assert verdict["status"] != "PROVED", ("false PROVED", m["kind"], m["witness"])


@pytest.mark.parametrize("label", ["empty_pre", "missing_outputs", "wrong_receipt", "arbitrary_garbage",
                                   "bad_clean", "missing_domain_policy", "stale_partial_domains"])
def test_fail_closed(setup, label):
    prog, layout, samples, contract, receipt, spec = setup
    domains = check_domains(prog["qasm"], layout, contract=contract)
    change = {
        "empty_pre": dict(pre=lambda S: ZERO),
        "missing_outputs": dict(symbolic=lambda S: {}),
        "wrong_receipt": dict(validation=dict(receipt, qasm_sha256="0" * 64)),
        "arbitrary_garbage": dict(layout=dict(layout, garbage=[0])),
        "bad_clean": dict(layout=dict(layout, initial_zero=[])),
        "missing_domain_policy": dict(layout={k: v for k, v in layout.items() if k != "domain_policy"}),
        "stale_partial_domains": dict(domains=dict(domains, qasm_sha256="0" * 64)),
    }[label]
    kwargs = dict(qasm=prog["qasm"], layout=layout, symbolic=spec, validation=receipt, contract=contract,
                  domains=domains)
    kwargs.update(change)
    verdict = prove_program(**kwargs)
    assert verdict["status"] != "PROVED", (label, verdict)
