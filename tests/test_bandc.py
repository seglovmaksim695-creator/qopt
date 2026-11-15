"""pytest tests/test_bandc.py -- scalable logical-AND compiler (qoptct.contrib.bandc) and its built-in checks."""
import os
import pathlib
import random
import re
import tempfile
import time

import pytest

from qoptct import verify_qasm
from qoptct.contrib import bandc
from qoptct.contrib.bandc import compiler as M
from qoptct.contrib.bandc import lcheck, vpt, worker

BENCH = pathlib.Path(os.environ["QOPT_BENCH"])
STD = BENCH / "std"


def need(name):
    p = STD / f"{name}.qc"
    if not p.exists():
        pytest.skip(f"benchmark {p} not available (set QOPT_BENCH)")
    return str(p)


def counts(txt):
    T = len(re.findall(r"\b(t|tdg) q", txt)); CX = len(re.findall(r"\b(cx|cz) q", txt))
    N = int(re.search(r"qubit\[(\d+)\]", txt).group(1))
    return T, CX, N


def cuccaro(path, n):
    """Cuccaro ripple-carry adder (MAJ/UMA), 2n+2 wires, 2n Toffolis: b := a + b, z ^= carry."""
    names = ["c0"] + [x for i in range(n) for x in (f"b{i}", f"a{i}")] + ["z"]
    g = []
    def MAJ(c, b, a): g.extend([(a, b), (a, c), (c, b, a)])
    def UMA(c, b, a): g.extend([(c, b, a), (a, c), (c, b)])
    MAJ("c0", "b0", "a0")
    for i in range(1, n): MAJ(f"a{i-1}", f"b{i}", f"a{i}")
    g.append((f"a{n-1}", "z"))
    for i in range(n - 1, 0, -1): UMA(f"a{i-1}", f"b{i}", f"a{i}")
    UMA("c0", "b0", "a0")
    with open(path, "w") as f:
        f.write(".v " + " ".join(names) + "\n.i " + " ".join(names) + "\n\nBEGIN\n")
        for x in g: f.write("tof " + " ".join(x) + "\n")
        f.write("END\n")
    return str(path)


def compile_ops(qc, clean, **kw):
    names, hl = M.load_qc(qc); n = len(names)
    M.reset_monos()
    C = M.Compiler(n, M.logical(hl, n), kw.pop("anc", 16), random.Random(0),
                   inputs=M.inputs_idx(qc, names) if clean else None, **kw)
    ops, N = M.compact(C.run(), n)
    return ops, N, n


def test_small_statevector_verified():
    """<= 16 qubits: candidates pass the statevector checker independently; T below the unitary 4 Toffolis (28 T)."""
    qc = need("barenco_tof_4")
    c = bandc.candidates(qc, clean=True, budget_sec=12)
    assert c, "no candidate"
    for d in c:
        T, CX, N = counts(d["qasm"])
        assert N <= 16 and "statevector" in d["note"]
        fd, tp = tempfile.mkstemp(suffix=".qasm"); os.write(fd, d["qasm"].encode()); os.close(fd)
        try:
            assert verify_qasm.check(tp, qc, trials=6, clean=True) > 1 - 1e-6
        finally:
            os.unlink(tp)
    assert min(counts(d["qasm"])[0] for d in c) < 28


def test_csum_mux_9_strict_reproduces_lab():
    """lab (agent b1): csum_mux_9 strict 56/112 (+1); unitary base 84/112."""
    qc = need("csum_mux_9")
    c = bandc.candidates(qc, clean=False, budget_sec=8)
    pts = [counts(d["qasm"]) for d in c]
    assert any(T <= 56 and CX <= 112 for T, CX, N in pts), pts
    assert all("lcheck(exact" in d["note"] and "canaries caught" in d["note"] for d in c)


@pytest.mark.slow
def test_large_cuccaro_vpt(tmp_path):
    """66-wire Cuccaro adder (beyond the 22-qubit statevector limit): verified by lcheck + vpt, T = 4n (Gidney)."""
    qc = cuccaro(tmp_path / "cuccaro_32.qc", 32)
    c = bandc.candidates(qc, clean=False, budget_sec=25)
    assert c
    T, CX, N = counts(c[0]["qasm"])
    assert N > 22 and T <= 4 * 32 and "vpt" in c[0]["note"]


def test_lcheck_canaries():
    qc = need("csum_mux_9")
    ops, N, n = compile_ops(qc, False, lazy=True, lctl=True)
    assert lcheck.check_ops(ops, N, qc)["ok"]
    muts = worker._ops_canaries(ops, random.Random(3))
    assert {m for m, _ in muts} == {"drop_corr_cz", "drop_cx"}
    for kind in ("drop_corr_cz", "drop_cx"):
        assert not all(lcheck.check_ops(m, N, qc)["ok"] for nm, m in muts if nm == kind)
    ok, note = worker.check_ops(ops, N, qc, False, time.time() + 30)
    assert ok, note
    bad = list(ops); i = next(i for i, (g, q) in enumerate(bad) if g == "AND"); del bad[i]
    assert not lcheck.check_ops(bad, N, qc)["ok"]           # dropping a logical AND is always caught


def test_vpt_canaries():
    qc = need("csum_mux_9")
    ops, N, n = compile_ops(qc, False, lazy=True, lctl=True)
    tr = M.to_qasm(ops, N, ndata=n)
    assert vpt.check_text(tr, qc, x_trials=3, m_trials=2)["ok"]
    bad = [vpt.check_text(t, qc, x_trials=3, m_trials=2)["ok"] for nm, t in worker.mutants(tr) if nm == "t_to_tdg"]
    assert bad == [False]
    drops = [vpt.check_text(t, qc, x_trials=3, m_trials=2)["ok"] for nm, t in worker.mutants(tr) if nm == "drop_cx"]
    assert drops and not all(drops)
    ok, note = worker.check_qasm(tr, worker.strip_resets(tr), N, qc, False, time.time() + 30, False)
    assert ok, note


def test_strip_resets_guard():
    good = "h q[3];\nc[0] = measure q[3];\nif (c[0]) {cz q[1], q[2]; x q[3];}\nreset q[3];\nx q[1];"
    assert "reset" not in worker.strip_resets(good)
    with pytest.raises(ValueError):
        worker.strip_resets("h q[3];\nc[0] = measure q[3];\nif (c[0]) {cz q[1], q[2];}\nreset q[3];")


def test_non_classical_and_tiny_budget():
    t = time.time()
    assert bandc.candidates(need("qft_4"), clean=False, budget_sec=10) == []
    assert bandc.candidates(need("csum_mux_9"), clean=False, budget_sec=1) == []
    assert bandc.candidates("/nonexistent.qc", clean=False, budget_sec=10) == []
    assert time.time() - t < 15


def test_budget_respected(tmp_path):
    """a 1026-wire adder cannot finish in 6 s: the worker is stopped on time and nothing unverified is returned."""
    qc = cuccaro(tmp_path / "cuccaro_512.qc", 512)
    t = time.time()
    c = bandc.candidates(qc, clean=False, budget_sec=6)
    assert time.time() - t < 9
    assert c == []
