"""Qiskit bridge: translation, check_qiskit, library proofs, mutants.  Fast by default; large n marked slow."""
import math
import warnings

import numpy as np
import pytest

pytest.importorskip("qiskit")
warnings.filterwarnings("ignore")

from qiskit import QuantumCircuit
from qiskit.circuit.library import (CDKMRippleCarryAdder, DraperQFTAdder, PhaseGate, VBERippleCarryAdder)
from qiskit.quantum_info import Statevector

from qoptct.foreign.qiskit import (NeedsFractionalPhases, QiskitTranslationError, check_qiskit, run_classical,
                                   translate)
from qoptct.foreign.qiskit.library import FAMILIES, controlled, prove_case, cdkm
from qoptct.foreign.qiskit.mutants import make_mutants

FAST = ["cdkm-half", "cdkm-full", "cdkm-fixed", "vbe-half", "vbe-full", "half-adder-gate", "full-adder-gate",
        "modular-adder-gate", "comparator-geq", "comparator-lt", "comparator-gate", "and", "or", "xor-const",
        "c-cdkm-half", "c-cdkm-full", "c-cdkm-fixed", "c-vbe-half", "c-vbe-full"]


def test_translate_leaves_and_phase():
    qc = QuantumCircuit(3)
    qc.x(0); qc.cx(0, 1); qc.ccx(0, 1, 2); qc.swap(0, 2); qc.h(1); qc.t(1); qc.s(2); qc.p(math.pi / 2, 0)
    qc.p(7 * math.pi / 4, 1); qc.cz(0, 1)
    tr = translate(qc)
    assert tr.counts()["ccx"] == 1 and tr.counts()["s"] == 2 and tr.counts()["tdg"] == 1
    assert "ccx q[0], q[1], q[2];" in tr.qasm()


def test_fractional_phase_is_unknown_not_proved():
    qc = QuantumCircuit(1)
    qc.p(0.3, 0)
    with pytest.raises(NeedsFractionalPhases):
        translate(qc)
    r = check_qiskit(DraperQFTAdder(3, kind="half"), lambda v: {})
    assert r["status"] == "UNKNOWN" and r["needs"] == "fractional_phases"
    qc = QuantumCircuit(1)
    qc.p(math.pi / 8, 0)
    # 0.8: a reference circuit is a target for the dyadic prover, which decides (a rotation by pi/8 is not the identity)
    r = check_qiskit(qc, reference=QuantumCircuit(1))
    assert r["status"] == "REFUTED" and r["method"].startswith("dyadic") and r.get("witness") is not None, r
    same = QuantumCircuit(1)
    same.p(math.pi / 8, 0)
    assert check_qiskit(qc, reference=same)["status"] == "PROVED"
    off = QuantumCircuit(1)
    off.p(math.pi / 8 + 1e-6, 0)           # not a multiple of pi/2^k: never PROVED, never REFUTED by guessing
    assert check_qiskit(qc, reference=off)["status"] == "UNKNOWN"


def test_rejects_unsupported():
    qc = QuantumCircuit(1)
    qc.reset(0)
    with pytest.raises(QiskitTranslationError):
        translate(qc)
    assert check_qiskit("not a circuit", lambda v: {})["status"] == "UNKNOWN"
    assert check_qiskit(QuantumCircuit(1))["status"] == "UNKNOWN"


def test_measure_and_if_translate():
    qc = QuantumCircuit(2, 1)
    qc.h(0)
    qc.measure(0, 0)
    with qc.if_test((qc.clbits[0], 1)):
        qc.x(1)
    tr = translate(qc)
    assert [o[0] for o in tr.ops] == ["h", "measure", "if"]
    assert "if (c[0]) {x q[1];}" in tr.qasm()


@pytest.mark.parametrize("fam", FAST)
@pytest.mark.parametrize("n", [5, 12])
def test_library_proved(fam, n):
    r = prove_case(FAMILIES[fam](n), budget_sec=120, exhaustive_bits=0)
    assert r["status"] == "PROVED" and r["method"] == "kprove-bdd", r


def test_exhaustive_path_small():
    r = prove_case(FAMILIES["cdkm-half"](3))
    assert r["status"] == "PROVED" and r["method"] == "exhaustive"


@pytest.mark.parametrize("n", [2, 3, 5])
def test_vbe_fixed_is_refuted_with_qiskit_witness(n):
    r = prove_case(FAMILIES["vbe-fixed"](n), exhaustive_bits=0)
    assert r["status"] == "REFUTED" and r["witness"]
    # independent: Qiskit's own state-vector simulation of the ORIGINAL circuit shows the same defect
    a, b = r["witness"]["a"], r["witness"]["b"]
    qc = VBERippleCarryAdder(n, kind="fixed")
    sv = Statevector.from_int(a | (b << n), 2 ** qc.num_qubits).evolve(qc)
    k = int(np.argmax(np.abs(sv.data)))
    assert (k >> n) & ((1 << n) - 1) != (a + b) % (1 << n)


@pytest.mark.parametrize("fam", ["cdkm-half", "vbe-full", "c-cdkm-fixed", "comparator-geq"])
def test_translation_matches_qiskit_simulation(fam):
    """The classical replay of the translated leaves equals Qiskit's own state-vector evolution (n=3)."""
    case = FAMILIES[fam](3)
    qc = case.circuit
    tr = translate(qc)
    rng = np.random.default_rng(3)
    for x in rng.integers(0, 1 << qc.num_qubits, 24):
        sv = Statevector.from_int(int(x), 2 ** qc.num_qubits).evolve(qc)
        k = int(np.argmax(np.abs(sv.data)))
        assert abs(abs(sv.data[k]) - 1) < 1e-9
        assert k == run_classical(tr.ops, int(x))


@pytest.mark.parametrize("fam,n", [("cdkm-half", 10), ("vbe-half", 10), ("c-cdkm-half", 8), ("comparator-geq", 10),
                                   ("modular-adder-gate", 6)])
def test_mutants_never_proved(fam, n):
    case = FAMILIES[fam](n)
    clean = tuple(r for r in case.clean if r in {q.name for q in case.circuit.qregs})
    reg = next(q for q in case.circuit.qregs if q.name not in clean and q.name != "ctrl" and len(q) >= 3)
    wires = [case.circuit.find_bit(q).index for q in reg]
    mutants = make_mutants(case.circuit, rare_wires=wires)
    assert len(mutants) >= 20
    refuted = 0
    for label, mut in mutants:
        r = check_qiskit(mut, case.spec, clean=clean, budget_sec=60, exhaustive_bits=0)
        if r["status"] == "PROVED":       # allowed only for a genuinely equivalent mutant
            q = check_qiskit(mut, reference=case.circuit, clean=clean, ref_clean=clean, budget_sec=60)
            assert q["status"] == "PROVED", (label, r, q)
        else:
            assert r["status"] == "REFUTED", (label, r)
            refuted += 1
    assert refuted >= 20


def test_reference_mode():
    a = cdkm(6, "half").circuit
    b = VBERippleCarryAdder(6, kind="half")
    r = check_qiskit(a, reference=b, clean=("help",), ref_clean=("helper",), exhaustive_bits=0)
    assert r["status"] == "PROVED", r
    assert check_qiskit(a, reference=a, clean=("help",), ref_clean=("help",))["status"] == "PROVED"
    bad = QuantumCircuit(*b.qregs)
    bad.compose(b, inplace=True)
    bad.x(b.qubits[0])
    r = check_qiskit(a, reference=bad, clean=("help",), ref_clean=("helper",), exhaustive_bits=0)
    assert r["status"] == "REFUTED"


def test_spec_function_wrong_is_refuted():
    r = check_qiskit(CDKMRippleCarryAdder(4, kind="fixed"),
                     lambda v: {"b": (v["a"] + v["b"] + 1) % 16}, clean=("help",))
    assert r["status"] == "REFUTED" and "witness" in r


def test_controlled_is_identity_when_off():
    case = controlled(cdkm(4, "half"))
    tr = translate(case.circuit)
    assert tr.primitives.get("mcx3", 0) > 0
    for x in range(0, 1 << 8, 7):
        s = x << 1                      # ctrl = 0
        assert run_classical(tr.ops, s) == s


@pytest.mark.slow
@pytest.mark.parametrize("fam", ["cdkm-half", "cdkm-full", "cdkm-fixed", "vbe-half", "vbe-full", "half-adder-gate",
                                 "full-adder-gate", "modular-adder-gate", "comparator-geq", "comparator-lt",
                                 "and", "or", "xor-const", "c-cdkm-half", "c-cdkm-full", "c-vbe-half"])
def test_library_proved_n256(fam):
    r = prove_case(FAMILIES[fam](256), budget_sec=600)
    assert r["status"] == "PROVED", r


def test_registry_matches_families():
    import json, pathlib, qoptct.foreign.qiskit as pkg
    d = json.loads((pathlib.Path(pkg.__file__).parent / "data" / "registry.json").read_text())
    assert d["rows"]
    for r in d["rows"]:
        assert r["family"] in FAMILIES
        if r["status"] == "PROVED":
            assert r["max_n_proved"] >= 16
            # every PROVED row has canaries at its largest proved size (audit 0.8, X4)
            if not r["canaries"]:                     # allowed only when the row says so, in words (and so do the documents)
                assert r.get("canaries_note", "").startswith("none:"), r["family"]
                continue
            if r["canaries"]["n"] != r["max_n_proved"]:    # only with an explicit note in the row
                assert r.get("canaries_note", "").startswith("canaries at n = %d" % r["canaries"]["n"]), r["family"]
            assert r["canaries"]["unknown"] <= 2, r["family"]
            assert r["canaries"]["proved_false"] == 0
            assert r["canaries"]["refuted"] >= 20, r["family"]
