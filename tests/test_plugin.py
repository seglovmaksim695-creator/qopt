"""pytest tests/ -- Qiskit pass and transpiler stage plugin, CLI core and the equality prover."""
import pathlib
import subprocess
import sys

import pytest
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator
from qiskit.transpiler import CouplingMap, PassManager

import qoptct
from qoptct.qiskit_qopt import QoptMeasPass, QoptOptimizationPlugin, QoptPass

ROOT = pathlib.Path(__file__).resolve().parents[1]
BENCH = ROOT / "bench_qc"
CT = ["h", "x", "z", "s", "sdg", "t", "tdg", "cx"]  # Clifford+T basis


def tof3():
    qc = QuantumCircuit(5); qc.ccx(0, 1, 4); qc.ccx(2, 4, 3); qc.ccx(0, 1, 4); return qc


def bench(name):
    """Benchmark .qc file -> QuantumCircuit (Toffoli level)."""
    names, hl = qoptct.load_qc(str(BENCH / f"{name}.qc"))
    qc = QuantumCircuit(len(names))
    for g, q in hl:
        getattr(qc, g)(*q)
    return qc


def cx_t(qc):
    ops = qc.count_ops()
    return ops.get("cx", 0), ops.get("t", 0) + ops.get("tdg", 0)


# -- QoptPass -------------------------------------------------------------------------------------------------
def test_tof3_cx_and_equal():
    qc = tof3(); p = QoptPass(); out = PassManager([p]).run(qc)
    assert p.last["verified"] and out.count_ops().get("cx", 0) <= 12
    assert Operator(out).equiv(Operator(qc))
    assert Operator(out) == Operator(qc)  # global phase is corrected too (<= 10 qubits)


def test_pass_keeps_non_clifford_t_parts():
    # rz(0.3) and the barrier split the circuit; only the Clifford+T runs may change
    qc = tof3(); qc.rz(0.3, 2); qc.barrier(); qc.compose(tof3(), inplace=True)
    p = QoptPass(); out = PassManager([p]).run(qc)
    assert p.last["runs"] == 2 and p.last["replaced"] >= 1
    assert out.count_ops().get("rz", 0) == 1 and out.count_ops().get("barrier", 0) == 1
    assert Operator(out) == Operator(qc)


def test_pass_with_measurements():
    qc = QuantumCircuit(5, 5); qc.compose(tof3(), inplace=True); qc.measure(range(5), range(5))
    p = QoptPass(); out = PassManager([p]).run(qc)
    assert p.last["replaced"] == 1 and out.count_ops()["measure"] == 5
    assert out.data[-1].operation.name == "measure"
    body = out.remove_final_measurements(inplace=False)
    assert Operator(body) == Operator(tof3())


def test_pass_unsupported_only_is_noop():
    qc = QuantumCircuit(2); qc.rx(0.1, 0); qc.rzz(0.2, 0, 1)
    p = QoptPass(); out = PassManager([p]).run(qc)
    assert "skipped" in p.last and out == qc


def test_pass_respects_basis():
    # basis without "t": the optimized run would need t gates, so nothing may be replaced
    qc = tof3()
    p = QoptPass(basis_gates=["h", "x", "s", "sdg", "tdg", "cx", "ccx"]); out = PassManager([p]).run(qc)
    assert p.last["replaced"] == 0 and out == qc


# -- transpiler stage plugin ------------------------------------------------------------------------------------
def test_plugin_registered():
    from qiskit.transpiler.preset_passmanagers.plugin import list_stage_plugins
    assert "qoptct" in list_stage_plugins("optimization")
    assert "qoptct" in list_stage_plugins("init")


@pytest.mark.parametrize("name,max_cx", [("tof_3", 12), ("mod5_4", 14)])
def test_transpile_init_and_optimization(name, max_cx):
    qc = bench(name)
    out = transpile(qc, basis_gates=CT, optimization_level=2, init_method="qoptct",
                    optimization_method="qoptct", seed_transpiler=1)
    assert set(out.count_ops()) <= set(CT)
    assert Operator(out) == Operator(qc)
    assert cx_t(out)[0] <= max_cx


def test_transpile_init_with_coupling_map():
    qc = bench("mod5_4")
    cmap = CouplingMap.from_line(qc.num_qubits)
    out = transpile(qc, basis_gates=CT, coupling_map=cmap, optimization_level=3, init_method="qoptct",
                    optimization_method="qoptct", seed_transpiler=3)
    edges = set(cmap.get_edges())
    assert all(tuple(out.find_bit(q).index for q in i.qubits) in edges for i in out.data if len(i.qubits) == 2)
    assert Operator.from_circuit(out).equiv(Operator(qc))


@pytest.mark.parametrize("name", ["tof_3", "mod5_4"])
@pytest.mark.parametrize("level", [1, 3])
def test_transpile_optimization_method(name, level):
    qc = bench(name)
    ref = transpile(qc, basis_gates=CT, optimization_level=level, seed_transpiler=1)
    out = transpile(qc, basis_gates=CT, optimization_level=level, optimization_method="qoptct", seed_transpiler=1)
    assert set(out.count_ops()) <= set(CT)
    assert Operator(out).equiv(Operator(qc))
    assert Operator(out) == Operator(qc)
    assert cx_t(out) <= cx_t(ref)  # never worse than Qiskit's own stage of the same level (CX first, then T)


def test_transpile_level0_runs_qopt():
    qc = bench("tof_3")
    out = transpile(qc, basis_gates=CT, optimization_level=0, optimization_method="qoptct")
    ref = transpile(qc, basis_gates=CT, optimization_level=0)
    assert Operator(out) == Operator(qc) and cx_t(out) < cx_t(ref)


def test_transpile_with_coupling_map():
    qc = bench("tof_3")
    cmap = CouplingMap.from_line(qc.num_qubits)
    out = transpile(qc, basis_gates=CT, coupling_map=cmap, optimization_level=2,
                    optimization_method="qoptct", seed_transpiler=7)
    edges = set(cmap.get_edges())
    for ins in out.data:
        if len(ins.qubits) == 2:
            assert tuple(out.find_bit(q).index for q in ins.qubits) in edges
    assert Operator.from_circuit(out).equiv(Operator(qc))


def test_transpile_without_basis():
    qc = bench("tof_3")
    out = transpile(qc, optimization_method="qoptct", optimization_level=1)
    assert Operator(out) == Operator(qc) and cx_t(out)[0] <= 12


def test_plugin_object_level_none():
    from qiskit.transpiler import PassManagerConfig
    pm = QoptOptimizationPlugin().pass_manager(PassManagerConfig(basis_gates=CT), None)
    out = pm.run(tof3())
    assert Operator(out) == Operator(tof3())


# -- core ---------------------------------------------------------------------------------------------------------
def test_mod5_4_cli():
    names, ref, best, how = qoptct.optimize(str(BENCH / "mod5_4.qc"))
    assert how is not None and sum(g == "cx" for g, _ in best) <= 14


def test_canary_rejected():
    names, hl = qoptct.load_qc(str(BENCH / "tof_3.qc"))
    import evo_lib
    ref = evo_lib.lower_std(hl); bad = list(ref); i = next(k for k, (g, _) in enumerate(bad) if g == "t")
    bad[i] = ("tdg", bad[i][1])  # corruption: one T replaced by T-dagger
    assert qoptct.prove(bad, ref, len(names)) is not True


def test_meas_tof3():
    qc = tof3(); p = QoptMeasPass(ancillas=[4]); out = PassManager([p]).run(qc)
    assert p.last["verified"] is True and p.last["t_after"] < p.last["t_before"]
    assert out.count_ops().get("measure", 0) == 1


def test_meas_wrong_ancilla_rejected():
    # wire 4 is not a clean ancilla here (it is flipped first): the pass must not produce an unverified result
    qc = QuantumCircuit(5); qc.x(4); qc.ccx(0, 1, 4); qc.ccx(0, 1, 4)
    p = QoptMeasPass(ancillas=[4]); PassManager([p]).run(qc)
    assert p.last.get("verified") in (True, None) or "skipped" in p.last


# -- optional torch -------------------------------------------------------------------------------------------------
def test_default_modes_do_not_import_torch():
    code = ("import sys, qoptct; from qoptct import verify_qasm; from qoptct.qiskit_qopt import QoptPass;"
            f"r = qoptct.optimize({str(BENCH / 'tof_3.qc')!r}); assert r[3] is not None;"
            "assert 'torch' not in sys.modules, 'torch imported'")
    subprocess.run([sys.executable, "-c", code], check=True, cwd=str(ROOT / "tests"))


def test_net_without_torch_gives_clear_error(monkeypatch):
    monkeypatch.setitem(sys.modules, "torch", None)  # makes "import torch" fail
    with pytest.raises(ImportError, match=r"qoptct\[net\]"):
        qoptct.require_torch()
    with pytest.raises(ImportError, match=r"qoptct\[net\]"):
        QoptPass(net=True)


def test_cli_net_without_torch(tmp_path):
    code = ("import sys; sys.modules['torch'] = None; import qoptct;"
            f"sys.argv = ['qoptct', 'optimize', {str(BENCH / 'tof_3.qc')!r}, '--net']; qoptct.main()")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode == 2 and 'pip install "qoptct[net]"' in r.stderr
