"""pytest tests/test_hw.py -- qoptct.hw: device model, router (commutation, SWAP absorption, floating blocks, mirror
rule), finishing to native gates, equivalence check, compile_hw portfolio and the `qoptct hw` command.

Synthetic devices (line / ring / heavy-hex with random calibration) need only Qiskit; the tests on IBM calibration
snapshots (fake_torino, fake_sherbrooke) are skipped when qiskit-ibm-runtime is not installed."""
import json
import math
import random

import numpy as np
import pytest
from qiskit import QuantumCircuit

from qoptct.hw import router as R
from qoptct.hw import benchlib as B
from qoptct.hw.check import verify_compiled
from qoptct.hw.compiler import (VARIANTS, cleanup, compile_hw, finish, flatten, qiskit_o3, relocate, route_hw,
                                valid_for)
from qoptct.hw.device import Device
from qoptct.hw.metrics import analyze

try:
    import qiskit_ibm_runtime  # noqa: F401
    HAVE_IBM = True
except ImportError:
    HAVE_IBM = False
ibm = pytest.mark.skipif(not HAVE_IBM, reason="needs qiskit-ibm-runtime (fake backends)")


@pytest.fixture(scope="module")
def hh():
    return Device("heavyhex:3")


@pytest.fixture(scope="module")
def line8():
    return Device("line:8")


def rand_circuit(n, m, seed, gates=("cx", "cz", "rzz", "cp", "h", "t", "rx", "swap", "rxx")):
    rng = random.Random(seed)
    qc = QuantumCircuit(n)
    for _ in range(m):
        g = rng.choice(gates)
        a, b = rng.sample(range(n), 2)
        if g in ("cx", "cz", "swap"):
            getattr(qc, g)(a, b)
        elif g in ("rzz", "cp", "rxx"):
            getattr(qc, g)(rng.uniform(0.1, 3), a, b)
        elif g == "rx":
            qc.rx(rng.uniform(0.1, 3), a)
        else:
            getattr(qc, g)(a)
    return qc


# ---- device and metrics --------------------------------------------------------------------------------------
def test_device_synthetic(hh):
    assert hh.twoq == "cz" and hh.n == 19 and len(hh.edges) == 20
    assert all(0 <= e < 0.3 for e in hh.err2.values())
    assert hh.hop[0, 0] == 0 and (hh.hop == hh.hop.T).all()
    a, b = next(iter(hh.edges))
    assert hh.has_edge(a, b) and hh.has_edge(b, a) and hh.hop[a, b] == 1
    assert hh.idle_logfid(0, 0.0) == 0.0 and hh.idle_logfid(0, 1e-5) < 0.0


def test_metrics_counts(line8):
    qc = QuantumCircuit(8)
    qc.sx(0); qc.cz(0, 1); qc.rz(0.3, 1); qc.cz(1, 2); qc.cz(0, 1)
    m = analyze(qc, line8, [0, 1, 2])
    assert m["twoq"] == 3 and m["depth2q"] == 3 and m["oneq"] == 1 and m["used_qubits"] == 3 and m["offmap_2q"] == 0
    assert 0 < m["fid"] <= m["esp"] < 1
    want = math.log(1 - line8.err1[0]) + 2 * math.log(1 - line8.err2[(0, 1)]) + math.log(1 - line8.err2[(1, 2)])
    assert abs(m["log_gate"] - want) < 1e-12
    qc.cz(0, 5)
    assert analyze(qc, line8)["offmap_2q"] == 1


# ---- commutation axes and dependencies -------------------------------------------------------------------------
def test_axes():
    from qiskit.circuit.library import CXGate, CZGate, RZZGate, SwapGate, RXXGate, CPhaseGate, HGate, TGate, SXGate
    ax = lambda g: R._axes2(np.asarray(g.to_matrix()))  # noqa: E731
    assert ax(CXGate()) == ["z", "x"] and ax(CZGate()) == ["z", "z"] and ax(RZZGate(0.3)) == ["z", "z"]
    assert ax(CPhaseGate(0.7)) == ["z", "z"] and ax(RXXGate(0.4)) == ["x", "x"] and ax(SwapGate()) == [None, None]
    assert R._axis1(TGate().to_matrix()) == "z" and R._axis1(SXGate().to_matrix()) == "x"
    assert R._axis1(HGate().to_matrix()) is None
    assert R.kak_cost(np.asarray(SwapGate().to_matrix())) == 3 and R.kak_cost(np.asarray(CXGate().to_matrix())) == 1
    assert R.kak_cost(np.asarray(RZZGate(0.3).to_matrix())) == 2


def test_dag_commutation():
    qc = QuantumCircuit(4)
    qc.rzz(0.1, 0, 1); qc.rzz(0.2, 1, 2); qc.rzz(0.3, 2, 3); qc.rzz(0.4, 0, 3)   # all commute
    qc.h(1); qc.cz(1, 2)
    nodes = R.to_nodes(qc)
    npred, succs = R.build_dag(nodes, 4)
    assert npred[:4] == [0, 0, 0, 0]
    assert npred[4] == 2 and set(i for i in range(4) if 4 in succs[i]) == {0, 1}     # h(1) after both rzz on wire 1
    assert npred[5] == 1 and 5 in succs[4]   # cz(1,2) joins the z-group of wire 2 (no edge from rzz(2,3)), waits for h(1)
    npred0, _ = R.build_dag(R.to_nodes(qc, commute=False), 4)
    assert npred0[:4] == [0, 1, 1, 2]


def test_flatten_swap_elision():
    qc = QuantumCircuit(3)
    qc.h(0); qc.swap(0, 2); qc.cx(0, 1); qc.ccx(0, 1, 2)
    flat, wire = flatten(qc)
    assert all(i.operation.num_qubits <= 2 and i.operation.name != "swap" for i in flat.data)
    assert wire == [2, 1, 0]
    # flat acts on relabelled wires: logical i ends on wire[i]
    assert verify_compiled(qc, flat, [0, 1, 2], wire) is True
    assert verify_compiled(qc, flat, [0, 1, 2], [0, 1, 2]) is False


# ---- router -----------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("variant", sorted(VARIANTS))
def test_router_variants_valid_and_equivalent(hh, variant):
    commute, kw = VARIANTS[variant]
    for seed in range(4):
        qc = rand_circuit(6, 40, seed)
        flat, wire = flatten(qc)
        nodes = R.to_nodes(flat, commute=commute)
        dag = R.build_dag(nodes, 6)
        rng = random.Random(seed)
        l0 = R.greedy_layout(R.interaction(nodes, 6), hh, rng)
        assert len(set(l0)) == 6
        ops, lf, cost, nsw = R.route(nodes, dag, hh, l0, R.Params(seed=seed, **kw))
        assert sum(1 for o in ops if o[0] == "g") == len(nodes) and sum(1 for o in ops if o[0] == "s") == nsw
        circ = R.to_circuit(nodes, ops, hh.n)
        for ins in circ.data:
            if ins.operation.num_qubits == 2:
                a, b = (circ.find_bit(q).index for q in ins.qubits)
                assert hh.has_edge(a, b)
        fin = [lf[wire[i]] for i in range(6)]
        assert verify_compiled(qc, circ, l0, fin) is True
        native = finish(circ, hh)
        assert valid_for(native, hh)
        assert verify_compiled(qc, native, l0, fin) is True
        # the router's own estimate is an upper bound of what consolidation finds (noise off -> plain count)
        ops2, _, cost2, _ = R.route(nodes, dag, hh, l0, R.Params(seed=seed, noise=False, **kw))
        n2 = analyze(finish(R.to_circuit(nodes, ops2, hh.n), hh), hh)["twoq"]
        assert n2 <= round(cost2)


def test_floating_block_absorbs_over_commuting_gates():
    """rzz(0,1), rzz(1,2), then a gate that needs 0 next to 2 on a line: the SWAP of (0,1) or (1,2) can still merge with
    the earlier rzz on that pair although another (commuting) rzz was emitted in between."""
    dev = Device("line:3")
    qc = QuantumCircuit(3)
    qc.rzz(0.3, 0, 1); qc.rzz(0.5, 1, 2); qc.rzz(0.7, 0, 2)
    nodes = R.to_nodes(qc)
    dag = R.build_dag(nodes, 3)
    c_strict = R.route(nodes, dag, dev, [0, 1, 2], R.Params(noise=False, floating=False))[2]
    ops, lf, c_float, nsw = R.route(nodes, dag, dev, [0, 1, 2], R.Params(noise=False, floating=True))
    assert nsw == 1 and c_float == 7 and c_strict in (7, 9) and c_float <= c_strict
    circ = finish(R.to_circuit(nodes, ops, 3), dev)
    assert analyze(circ, dev)["twoq"] <= 7
    assert verify_compiled(qc, circ, [0, 1, 2], lf) is True


def test_mirror_rule_gives_linear_qft(line8):
    """QFT on a line in index order: with the mirror rule every controlled phase is merged with a SWAP
    (3 native gates per pair, the known linear-nearest-neighbour construction) or better."""
    qc = B.qft(8)
    cands = route_hw(qc, line8, trials=40, seed=1, variants=["mirror", "mirror+"], noise=False, do_relocate=False,
                     objective="twoq")
    best = cands[0]
    assert best.metrics["twoq"] <= 3 * 28
    assert verify_compiled(qc, best.circuit, best.init_map, best.final_map) is True
    plain = route_hw(qc, line8, trials=40, seed=1, variants=["sabre"], noise=False, do_relocate=False, objective="twoq")[0]
    assert best.metrics["twoq"] < plain.metrics["twoq"]


def test_layout_helpers(hh):
    rng = random.Random(3)
    p = R.find_path(hh, 9, rng)
    assert p and len(set(p)) == 9 and all(hh.has_edge(a, b) for a, b in zip(p, p[1:]))
    adj = R.interaction(R.to_nodes(B.qft(5)), 5)
    assert sorted(R.linear_order(adj, rng)) == list(range(5))
    assert R.line_layout([2, 0, 1], [7, 8, 9]) == [8, 9, 7]


def test_relocate_keeps_equivalence(hh):
    qc = rand_circuit(5, 25, 11)
    c = route_hw(qc, hh, trials=6, seed=2, do_relocate=False)[0]
    circ, init, fin = relocate(c.circuit, c.init_map, c.final_map, hh)
    assert valid_for(circ, hh)
    assert verify_compiled(qc, circ, init, fin) is True
    assert analyze(circ, hh, fin)["esp"] >= c.metrics["esp"] - 1e-12


# ---- verification ---------------------------------------------------------------------------------------------
def test_verify_detects_errors(hh):
    qc = rand_circuit(5, 20, 5)
    c = qiskit_o3(qc, hh, seed=1)
    assert verify_compiled(qc, c.circuit, c.init_map, c.final_map) is True
    bad = c.circuit.copy()
    bad.rz(0.01, c.final_map[0])
    bad.sx(c.final_map[0])
    assert verify_compiled(qc, bad, c.init_map, c.final_map) is False
    fm = list(c.final_map)
    fm[0], fm[1] = fm[1], fm[0]
    assert verify_compiled(qc, c.circuit, c.init_map, fm) is False
    assert verify_compiled(qc, c.circuit, c.init_map, c.final_map, max_qubits=2) is None


# ---- portfolio ------------------------------------------------------------------------------------------------
def test_compile_hw_synthetic(hh):
    qc = B.qaoa(B._graph("full", 5), 5, 1)
    r = compile_hw(qc, hh, seeds=(0, 1), trials=30)
    assert r.verified is True and valid_for(r.circuit, hh)
    assert {"qiskit", "hw"} <= set(r.best)
    assert r.metrics["fid"] >= r.best["qiskit"].metrics["fid"]
    assert all(c.verified is True for c in r.best.values())
    assert r.best["hw"].metrics["twoq"] <= r.best["qiskit"].metrics["twoq"]
    r2 = compile_hw(qc, hh, seeds=(0,), trials=10, objective="twoq", strategies=["qiskit"])
    assert r2.strategy == "qiskit"


def test_compile_hw_uses_qopt_logical(hh):
    """A Toffoli circuit: the proven logical CX/T optimization enters as its own candidates."""
    qc = QuantumCircuit(5)
    qc.ccx(0, 1, 4); qc.ccx(2, 4, 3); qc.ccx(0, 1, 4)
    r = compile_hw(qc, hh, seeds=(0,), trials=20, qopt_seconds=30)
    assert r.verified is True
    assert "qopt+qiskit" in r.best and r.best["qopt+qiskit"].verified is True
    assert r.best["qopt+qiskit"].metrics["twoq"] < r.best["qiskit"].metrics["twoq"]


def test_compile_hw_rejects_measurement_in_the_middle(hh):
    qc = QuantumCircuit(2, 1)
    qc.h(0); qc.measure(0, 0); qc.cx(0, 1)
    with pytest.raises(ValueError):
        compile_hw(qc, hh, strategies=["hw"], trials=2)


@ibm
def test_compile_hw_fake_torino_tiny():
    dev = Device("fake_torino")
    assert dev.n == 133 and dev.twoq == "cz" and dev.bad_qubits
    for b in B.suite("tiny"):
        r = compile_hw(b.circuit, dev, seeds=(0, 1), trials=40, qopt_seconds=30)
        assert r.verified is True, b.name
        assert valid_for(r.circuit, dev)
        assert not (set(r.init_map) | set(r.final_map)) & dev.bad_qubits
        assert r.metrics["fid"] >= r.best["qiskit"].metrics["fid"]


@ibm
def test_ecr_backend():
    dev = Device("fake_sherbrooke")
    assert dev.twoq == "ecr"
    qc = B.heisenberg(B._graph("ring", 6), 6, 1)
    r = compile_hw(qc, dev, seeds=(0,), trials=30, strategies=["qiskit", "hw"])
    for c in r.best.values():
        assert valid_for(c.circuit, dev) and c.verified is True
        assert set(c.metrics["twoq_names"]) == {"ecr"}


# ---- command line ---------------------------------------------------------------------------------------------
def test_cli(tmp_path, capsys):
    from qiskit import qasm2, qasm3
    from qoptct.hw.cli import main_hw
    qc = B.qft(4)
    src = tmp_path / "in.qasm"
    src.write_text(qasm2.dumps(qc))
    out, rep = tmp_path / "out.qasm", tmp_path / "rep.json"
    rc = main_hw([str(src), "--backend", "heavyhex:3", "-o", str(out), "--report", str(rep), "--seeds", "2",
                  "--seconds", "2", "--measure"])
    assert rc == 0
    text = capsys.readouterr().out
    assert "selected:" in text and "init_map" in text
    j = json.loads(rep.read_text())
    assert j["verified"] is True and len(j["init_map"]) == 4 and j["metrics"]["twoq"] > 0
    c = qasm3.loads(out.read_text())
    assert c.num_clbits == 4 and c.count_ops().get("measure") == 4
    assert main_hw([str(tmp_path / "missing.qasm")]) == 3
    assert main_hw([str(src), "--backend", "no_such_backend"]) == 3


def test_cli_entry_point(tmp_path):
    import subprocess
    import sys
    r = subprocess.run([sys.executable, "-m", "qoptct.qopt_core", "hw", "--list-backends"], capture_output=True, text=True)
    assert r.returncode == 0 and "heavyhex:D" in r.stdout
