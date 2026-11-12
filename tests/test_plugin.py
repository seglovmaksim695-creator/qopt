"""pytest tests/ -- Qiskit plugin, CLI core and the equality prover."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
from qiskit import QuantumCircuit
from qiskit.transpiler import PassManager
from qiskit.quantum_info import Operator
import qopt
from qopt.qiskit_qopt import QoptPass, QoptMeasPass


def tof3():
    qc = QuantumCircuit(5); qc.ccx(0, 1, 4); qc.ccx(2, 4, 3); qc.ccx(0, 1, 4); return qc


def test_tof3_cx_and_equal():
    qc = tof3(); p = QoptPass(); out = PassManager([p]).run(qc)
    assert p.last["verified"] and out.count_ops().get("cx", 0) <= 12
    assert Operator(out).equiv(Operator(qc))


def test_mod5_4_cli():
    names, ref, best, how = qopt.optimize(str(ROOT / "bench_qc" / "mod5_4.qc"))
    assert how is not None and sum(g == "cx" for g, _ in best) <= 14


def test_canary_rejected():
    names, hl = qopt.load_qc(str(ROOT / "bench_qc" / "tof_3.qc"))
    import evo_lib
    ref = evo_lib.lower_std(hl); bad = list(ref); i = next(k for k, (g, _) in enumerate(bad) if g == "t")
    bad[i] = ("tdg", bad[i][1])  # corruption: one T replaced by T-dagger
    assert qopt.prove(bad, ref, len(names)) is not True


def test_meas_tof3():
    qc = tof3(); p = QoptMeasPass(ancillas=[4]); out = PassManager([p]).run(qc)
    assert p.last["verified"] is True and p.last["t_after"] < p.last["t_before"]
    assert out.count_ops().get("measure", 0) == 1


def test_meas_wrong_ancilla_rejected():
    # wire 4 is not a clean ancilla here (it is flipped first): the pass must not produce an unverified result
    qc = QuantumCircuit(5); qc.x(4); qc.ccx(0, 1, 4); qc.ccx(0, 1, 4)
    p = QoptMeasPass(ancillas=[4]); PassManager([p]).run(qc)
    assert p.last.get("verified") in (True, None) or "skipped" in p.last
