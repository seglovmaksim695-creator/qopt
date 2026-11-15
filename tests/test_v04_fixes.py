"""Release fixes of 0.4.0: `qoptct prove` with .qc circuits, --time-limit over the whole run (also the unitary
portfolio), one-line input errors with exit code 3, OpenQASM output for -o *.qasm."""
import json
import pathlib
import subprocess
import sys
import time
import types

import pytest

import qoptct
from qoptct import qopt_core as core

ROOT = pathlib.Path(__file__).resolve().parents[1]
BENCH = ROOT / "bench_qc"
EXE = [sys.executable, "-c", "import sys; from qoptct.qopt_core import main; sys.argv[0] = 'qoptct'; main()"]


def run_cli(args):
    argv = sys.argv
    sys.argv = ["qoptct"] + [str(a) for a in args]
    try:
        qoptct.main()
    finally:
        sys.argv = argv


def cli(*args, timeout=300):
    return subprocess.run(EXE + [str(a) for a in args], capture_output=True, text=True, timeout=timeout)


def _verify_ready():
    try:
        import qoptct.verify as v
        return core.contrib_available("verify") and callable(getattr(v, "prove", None))
    except Exception:
        return False


def _equal_to_input(out, ref_qc):
    """The written unitary result (OpenQASM 3 or .qc) implements the input (full unitary, <= 10 qubits)."""
    import evo_lib  # noqa: F401  (qoptct/portfolio is on sys.path after an optimize run)
    names, hl = qoptct.load_qc(str(ref_qc))
    ref = evo_lib.lower_std(hl)
    text = pathlib.Path(out).read_text()
    if text.startswith("OPENQASM"):
        n, c = core.qasm_unitary(text)
        assert n == len(names)
    else:
        c = evo_lib.lower_std(qoptct.load_qc(str(out))[1])
    return qoptct.prove(c, ref, len(names)) is True


# -- D1: qoptct prove with a .qc circuit ---------------------------------------------------------------------
@pytest.mark.skipif(not _verify_ready(), reason="qoptct.verify not installed")
def test_prove_qc_circuit_proved_and_refuted(tmp_path):
    out = tmp_path / "o.qc"
    run_cli(["optimize", BENCH / "tof_3.qc", "-o", out, "--mode", "cx", "--no-contrib"])
    r = cli("prove", out, BENCH / "tof_3.qc", "--budget", 60)
    assert r.returncode == 0, r.stdout + r.stderr
    assert json.loads(r.stdout.strip().splitlines()[-1])["status"] == "PROVED"
    assert "Traceback" not in r.stderr
    lines = out.read_text().splitlines()
    i = next(k for k, s in enumerate(lines) if s.startswith("T "))
    lines[i] = "T* " + lines[i].split()[1]  # canary: one T replaced by T-dagger
    bad = tmp_path / "bad.qc"
    bad.write_text("\n".join(lines) + "\n")
    r = cli("prove", bad, BENCH / "tof_3.qc", "--budget", 60)
    assert r.returncode == 1, r.stdout + r.stderr
    assert json.loads(r.stdout.strip().splitlines()[-1])["status"] == "REFUTED"


# -- D2: --time-limit also bounds the unitary portfolio --------------------------------------------------------
def test_time_limit_cuts_hanging_portfolio(tmp_path, monkeypatch):
    """Every portfolio program hangs: the run still ends at the limit; the programs stop early enough to prove the
    phase-folded input (or, at worst, the input itself is written)."""
    import qoptct.portfolio  # noqa: F401
    monkeypatch.setattr(core, "_mod", lambda path: types.SimpleNamespace(run=lambda hl, n: time.sleep(600)))
    out, rep = tmp_path / "o.qasm", tmp_path / "r.json"
    t = time.time()
    run_cli(["optimize", BENCH / "tof_3.qc", "-o", out, "--mode", "t", "--time-limit", 3, "--report", rep])
    assert time.time() - t < 30
    r = json.loads(rep.read_text())
    assert r["search_seconds"] <= 3.5
    fin = r["final"]
    assert fin["format"] == "qc" and fin["file"] == "openqasm3"
    assert r["method"] in ("input+fold", "input (not proven)") and fin["verified"] == (r["method"] == "input+fold")
    assert fin["t"] <= r["t_in"]
    assert out.read_text().startswith("OPENQASM 3.0;")
    assert _equal_to_input(out, BENCH / "tof_3.qc")
    assert core._END is None


def test_time_limit_cuts_hanging_proof(tmp_path, monkeypatch):
    """The candidate proofs of the portfolio hang: cut at the limit; nothing unproven is written as a result."""
    monkeypatch.setattr(core, "prove", lambda *a: time.sleep(600))
    monkeypatch.setattr(core, "_prove_unitary_exact", lambda *a: None)
    out, rep = tmp_path / "o.qc", tmp_path / "r.json"
    t = time.time()
    run_cli(["optimize", BENCH / "tof_3.qc", "-o", out, "--time-limit", 4, "--no-contrib", "--report", rep])
    assert time.time() - t < 30
    r = json.loads(rep.read_text())
    assert r["search_seconds"] <= 4.5 and not r["final"]["verified"] and r["method"] == "input (not proven)"
    assert qoptct.load_qc(str(out))[1] == qoptct.load_qc(str(BENCH / "tof_3.qc"))[1] or _equal_to_input(out, BENCH / "tof_3.qc")


def test_time_limit_writes_best_checked_result(tmp_path):
    """A limit that leaves time for the portfolio: its proven result is written, unitary, as OpenQASM 3."""
    out, rep = tmp_path / "o.qasm", tmp_path / "r.json"
    run_cli(["optimize", BENCH / "tof_3.qc", "-o", out, "--mode", "t", "--time-limit", 20, "--sel-sec", 3,
             "--and-variants", 0, "--no-contrib", "--report", rep])
    r = json.loads(rep.read_text())
    assert r["search_seconds"] <= 21
    fin = r["final"]
    assert fin["verified"]
    if fin["format"] == "qc":
        assert fin["file"] == "openqasm3" and _equal_to_input(out, BENCH / "tof_3.qc")


def test_within_and_limited_respect_the_run_limit(monkeypatch):
    monkeypatch.setattr(core, "_END", time.time() + 0.5)
    assert 0 < core._within(600) <= 0.5 and core._within(0) <= 0.5
    t = time.time()
    with pytest.raises(core._TO):
        core.limited(600, time.sleep, 30)
    assert time.time() - t < 3
    monkeypatch.setattr(core, "_END", time.time() - 1)
    assert core._within(600) is None
    assert core.call_contrib("verify", "prove", ("x", "y", False, 60), seconds=60)[1] in ("skipped (time limit)",
                                                                                          "not installed")
    monkeypatch.setattr(core, "_END", None)
    assert core._within(600) == 600 and core._within(0) == 0


# -- D3: input errors -> one line on stderr, exit code 3 ---------------------------------------------------------
BAD = {
    "garbage.qc": "this is not a circuit\n\x00\x01garbage\n",
    "empty.qc": "",
    "unknown.qc": ".v a b\n.i a b\n\nBEGIN\nH a\nFOO a b\nT a\nEND\n",
    "undeclared.qc": ".v a b\n.i a b\n\nBEGIN\nH a\ntof a q\nT a\nEND\n",
    "wide.qc": ".v a b c d e\n\nBEGIN\ntof a b c d e\nEND\n",
}


@pytest.fixture
def bad_inputs(tmp_path):
    for name, text in BAD.items():
        (tmp_path / name).write_text(text)
    return tmp_path


@pytest.mark.parametrize("name", list(BAD) + ["missing.qc"])
def test_optimize_input_errors(bad_inputs, name):
    r = cli("optimize", bad_inputs / name, "-o", bad_inputs / "o.qc", "--time-limit", 30, timeout=120)
    assert r.returncode == 3, r.stdout + r.stderr
    err = [s for s in r.stderr.splitlines() if s.strip()]
    assert len(err) == 1 and err[0].startswith("error: ") and "Traceback" not in r.stderr, r.stderr
    assert not (bad_inputs / "o.qc").exists()


def test_load_qc_errors_are_value_errors(bad_inputs):
    for name in BAD:
        with pytest.raises(ValueError):
            qoptct.load_qc(str(bad_inputs / name))


def test_optimize_usage_and_output_errors(tmp_path):
    r = cli("optimize", BENCH / "tof_3.qc", "-o", tmp_path / "no_such_dir" / "o.qc", timeout=120)
    assert r.returncode == 3 and r.stderr.startswith("error: output directory does not exist")
    r = cli("optimize", timeout=120)
    assert r.returncode == 3 and "Traceback" not in r.stderr


@pytest.mark.parametrize("circuit,ref", [("missing.qasm", "ok"), ("garbage.qc", "ok"), ("badgate.qasm", "ok"),
                                         ("empty.qc", "ok"), ("ok", "garbage.qc"), ("ok", "missing.qc"),
                                         ("noreg.qasm", "ok")])
def test_prove_input_errors_never_refuted(bad_inputs, circuit, ref):
    (bad_inputs / "badgate.qasm").write_text('OPENQASM 3.0;\ninclude "stdgates.inc";\nqubit[5] q;\nfoo q[0];\n')
    (bad_inputs / "noreg.qasm").write_text('OPENQASM 3.0;\nh q[0];\n')
    path = lambda s: BENCH / "tof_3.qc" if s == "ok" else bad_inputs / s
    r = cli("prove", path(circuit), path(ref), "--budget", 30, timeout=120)
    if circuit in ("badgate.qasm", "noreg.qasm") and not _verify_ready():
        assert r.returncode == 2  # no prover installed: UNKNOWN (the OpenQASM parser belongs to qoptct.verify)
        return
    assert r.returncode == 3, r.stdout + r.stderr
    err = [s for s in r.stderr.splitlines() if s.strip()]
    assert len(err) == 1 and err[0].startswith("error: ") and "Traceback" not in r.stderr, r.stderr


@pytest.mark.skipif(not _verify_ready(), reason="qoptct.verify not installed")
def test_verify_cli_input_errors(bad_inputs):
    from qoptct.verify import main_prove
    assert main_prove([str(bad_inputs / "missing.qasm"), str(BENCH / "tof_3.qc")]) == 3
    assert main_prove([str(BENCH / "tof_3.qc"), str(bad_inputs / "missing.qc")]) == 3


# -- final check: a text proven earlier in the run is not thrown away when the short final re-proof times out -----
GOOD = ROOT / "results" / "best" / "strict" / "tof_3.qasm"
FAKE_SRC = "def candidates(qc_path, clean, budget_sec=60.0):\n    return [{'qasm': open(%r).read(), 'note': 'known good'}]\n"
FAKE_VER = ("import os\nfrom qoptct.verify import check_dense, check_sparse, check_big, quick_check\n"
            "def prove(*a, **k):  # PROVED on the first call, UNKNOWN (out of time) afterwards\n"
            "    m = %r\n"
            "    if os.path.exists(m):\n"
            "        return {'status': 'UNKNOWN', 'method': 'fake', 'detail': 'time budget exceeded'}\n"
            "    open(m, 'w').close()\n"
            "    return {'status': 'PROVED', 'method': 'fake-proof', 'detail': ''}\n")


@pytest.mark.skipif(not (GOOD.exists() and _verify_ready()), reason="results/best or qoptct.verify not present")
def test_final_check_keeps_result_proven_earlier(tmp_path, monkeypatch):
    pkg = tmp_path / "qfixfake"
    pkg.mkdir()
    (pkg / "__init__.py").write_text("")
    (pkg / "src.py").write_text(FAKE_SRC % str(GOOD))
    (pkg / "ver.py").write_text(FAKE_VER % str(tmp_path / "proved_once"))
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setenv("PYTHONPATH", str(tmp_path))
    monkeypatch.setattr(core, "MODULES", {"contrib.arith": "qfixfake.src", "verify": "qfixfake.ver",
                                          "contrib.bandc": "qfixfake.none", "contrib.gf2m": "qfixfake.none",
                                          "contrib.cxpass": "qfixfake.none"})
    core._PROVEN.clear()
    out, rep = tmp_path / "o.qasm", tmp_path / "r.json"
    run_cli(["optimize", BENCH / "tof_3.qc", "-o", out, "--mode", "t", "--sel-sec", 3, "--tgad-sec", 0,
             "--and-variants", 0, "--contrib-sec", 30, "--report", rep])
    fin = json.loads(rep.read_text())["final"]
    assert fin["format"] == "qasm" and (fin["t"], fin["cx"]) == (8, 10), fin
    assert fin["check"].startswith("PROVED earlier in this run (fake-proof)"), fin
    assert out.read_text() == GOOD.read_text()


# -- D6: the written OpenQASM is read by qiskit.qasm3.loads (no ^ / && in conditions) ----------------------------
def _qiskit_loads():
    try:
        import qiskit_qasm3_import  # noqa: F401
        from qiskit import qasm3
        return qasm3.loads
    except ImportError:
        return None


XOR_QASM = ("OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[4] q;\nbit[2] c;\n"
            "h q[2];\ncx q[0], q[2];\nh q[3];\nc[0] = measure q[2];\nif ((c[0])) {x q[2];}\n"
            "c[1] = measure q[3];\nif ((c[0] ^ c[1])) {cz q[0], q[1];}\nif ((c[1])) {x q[3];}\n"
            "if (c[0] && c[1]) { z q[1]; }\nt q[0];\n")


def test_make_qiskit_readable_forms():
    assert not core.qiskit_readable(XOR_QASM)
    for parity in (True, False):
        r = core.make_qiskit_readable(XOR_QASM, parity=parity)
        assert r is not None and core.qiskit_readable(r) and "^" not in r and "&&" not in r
        st0, st = core.qasm_stats(XOR_QASM), core.qasm_stats(r)
        assert st[1] == st0[1]
        # parity form: no extra CX (one more measurement); copies form: the CZ body twice
        assert (st[2], st[3]) == ((st0[2], st0[3] + 1) if parity else (st0[2] + 1, st0[3]))
        loads = _qiskit_loads()
        if loads:
            loads(r)
    # the same operation on every branch: compare the two forms branch by branch with the dense checker of verify
    from qoptct.verify._parse import parse_qasm_text
    assert parse_qasm_text(core.make_qiskit_readable(XOR_QASM))[2]
    assert core.qiskit_readable('OPENQASM 3.0;\nqubit[1] q;\nbit[1] c;\nc[0] = measure q[0];\n'
                                'if (c[0]) { if (!c[0]) { x q[0]; } }\n')


QISKIT_CASES = [("qft_4", "balanced", []), ("mod5_4", "balanced", []), ("tof_3", "balanced", [])]


def _run_readable(tmp_path, circ, mode, extra):
    out, rep = tmp_path / f"{circ}_{mode}.qasm", tmp_path / f"{circ}_{mode}.json"
    run_cli(["optimize", BENCH / f"{circ}.qc", "-o", out, "--mode", mode, "--report", rep] + extra)
    fin = json.loads(rep.read_text())["final"]
    text = out.read_text()
    assert core.qiskit_readable(text) and "^" not in text and "&&" not in text, fin
    loads = _qiskit_loads()
    if loads:
        loads(text)
    if _verify_ready():
        r = cli("prove", out, BENCH / f"{circ}.qc", "--budget", 120)
        assert r.returncode == 0, r.stdout + r.stderr
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted([out.name, rep.name])  # Z11: no stage files left
    return fin, text


@pytest.mark.parametrize("circ,mode,extra", QISKIT_CASES)
def test_output_read_by_qiskit_and_proved(tmp_path, circ, mode, extra):
    fin, text = _run_readable(tmp_path, circ, mode, extra)
    if circ == "qft_4" and _verify_ready():  # the CX post-pass writes an XOR condition here: rewritten, same CX
        assert fin.get("qiskit_rewrite", {}).get("variant", "").startswith("rewritten"), fin
        assert (fin["t"], fin["cx"]) == (62, 52), fin


@pytest.mark.slow
@pytest.mark.parametrize("circ,mode", [("qft_4", "t"), ("qft_4", "t --clean"), ("csla_mux_3", "balanced")])
def test_output_read_by_qiskit_and_proved_slow(tmp_path, circ, mode):
    m = mode.split()
    _run_readable(tmp_path, circ, m[0], m[1:])


@pytest.mark.skipif(not _verify_ready(), reason="qoptct.verify not installed")
def test_rewritten_output_canaries(tmp_path):
    """Mutants of a rewritten (Qiskit-readable) output are REFUTED: a dropped parity X, a dropped correction."""
    out = tmp_path / "o.qasm"
    run_cli(["optimize", BENCH / "qft_4.qc", "-o", out, "--mode", "balanced"])
    lines = out.read_text().splitlines()
    k = max(i for i, s in enumerate(lines) if "measure" in s)  # the measurement of the computed parity
    muts = {"drop_parity_x": lines[:k - 1] + lines[k:],
            "drop_correction": [s for i, s in enumerate(lines) if not (i > k and s.startswith("if") and "cz" in s)]}
    for name, m in muts.items():
        p = tmp_path / f"{name}.qasm"
        p.write_text("\n".join(m) + "\n")
        r = cli("prove", p, BENCH / "qft_4.qc", "--budget", 120)
        assert r.returncode == 1, (name, r.stdout + r.stderr)


def test_keep_temp_and_atomic_output(tmp_path):
    out, rep = tmp_path / "o.qasm", tmp_path / "r.json"
    run_cli(["optimize", BENCH / "tof_3.qc", "-o", out, "--mode", "balanced", "--keep-temp", "--report", rep])
    r = json.loads(rep.read_text())
    work = pathlib.Path(r["temp_dir"])
    assert work.is_dir() and any(work.iterdir())
    assert sorted(p.name for p in tmp_path.iterdir()) == ["o.qasm", "r.json"]
    import shutil
    shutil.rmtree(work)


# -- Z10: malformed OpenQASM / binary files -> exit code 3 ------------------------------------------------------
@pytest.mark.parametrize("text", ['OPENQASM 3.0;\nqubit[5] q;\nh q[0;\n',
                                  'OPENQASM 3.0;\nqubit[5] q;\nbit[1] c;\nif (c[0]) { x q[1]; }\n',
                                  'OPENQASM 3.0;\nqubit[5] q;\ncx q[0], q[9];\n',
                                  b'\xff\xfe\x00binary\x80'])
def test_prove_malformed_qasm_exit_3(tmp_path, text):
    p = tmp_path / "bad.qasm"
    (p.write_bytes if isinstance(text, bytes) else p.write_text)(text)
    r = cli("prove", p, BENCH / "tof_3.qc", "--budget", 30, timeout=120)
    assert r.returncode == 3 and r.stderr.strip().startswith("error: ") and "Traceback" not in r.stderr, r.stderr
    if _verify_ready():
        r = subprocess.run([sys.executable, "-m", "qoptct.verify", str(p), str(BENCH / "tof_3.qc"), "--budget", "30"],
                           capture_output=True, text=True, timeout=120)
        assert r.returncode == 3 and "error: " in r.stderr, r.stdout + r.stderr
