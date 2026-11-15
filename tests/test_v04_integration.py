"""v0.4 integration: optional modules (candidate sources, CX post-pass, exact prover) behind the core.

Fake modules check the safety rules of the core (every module result is checked again; exceptions, hangs and
garbage never break the optimizer); the real modules are exercised when they are installed."""
import json
import pathlib
import re
import subprocess
import sys
import time

import pytest

import qoptct
from qoptct import qopt_core as core
from qoptct import verify_qasm

ROOT = pathlib.Path(__file__).resolve().parents[1]
BENCH = ROOT / "bench_qc"
GOOD = ROOT / "results" / "best" / "strict" / "tof_3.qasm"  # tof_3, 8 T / 10 CX, 2 fresh ancillas, strict model
pytestmark = pytest.mark.skipif(not GOOD.exists(), reason="results/best not present")

FAKE = {
    "src_good": "def candidates(qc_path, clean, budget_sec=60.0):\n"
                "    return [{'qasm': open(%(good)r).read(), 'note': 'known good'}]\n",
    "src_bad": "def candidates(qc_path, clean, budget_sec=60.0):\n"  # canary: one T replaced by T-dagger
               "    t = open(%(good)r).read().replace('\\nt q', '\\ntdg q', 1)\n"
               "    return [{'qasm': t, 'note': 'corrupted'}]\n",
    "src_raise": "def candidates(qc_path, clean, budget_sec=60.0):\n    raise RuntimeError('boom')\n",
    "src_hang": "import time\ndef candidates(qc_path, clean, budget_sec=60.0):\n    time.sleep(3600)\n",
    "src_junk": "def candidates(qc_path, clean, budget_sec=60.0):\n"
                "    return [{'qasm': 'OPENQASM 3.0;\\nqubit[5] q;\\nrz(0.1) q[0];\\n'}, 5, None, {'note': 'x'}]\n",
    "src_exit": "import os\ndef candidates(qc_path, clean, budget_sec=60.0):\n    os._exit(3)\n",
    "ver_unknown": "from qoptct.verify import check_dense, check_sparse, check_big, quick_check\n"  # proof never finishes
                   "def prove(*a, **k):\n    return {'status': 'UNKNOWN', 'method': 'fake', 'detail': 'out of reach'}\n",
    "cx_bad": "import re\ndef improve(qasm, clean, budget_sec=60.0):\n"  # drops a CX: fewer CX, wrong circuit
              "    return re.sub(r'\\ncx [^\\n]*', '', qasm, count=1)\n",
    "cx_zcz": "import re\ndef improve(qasm, clean, budget_sec=60.0, ref=None):\n"  # correct, written with z and cz
              "    L = qasm.rstrip().split('\\n')[:-2]\n"  # drop the redundant trailing cx pair
              "    i = next(k for k, s in enumerate(L) if s.startswith('cx '))\n"
              "    a, b = re.findall(r'q\\[\\d+\\]', L[i])\n"
              "    L[i:i + 1] = [f'h {b};', f'cz {a}, {b};', f'h {b};']\n"
              "    return '\\n'.join(L + ['z q[0];', 'z q[0];']) + '\\n'\n",
}


@pytest.fixture
def fake(tmp_path, monkeypatch):
    pkg = tmp_path / "qv04fake"
    pkg.mkdir()
    (pkg / "__init__.py").write_text("")
    for k, v in FAKE.items():
        (pkg / f"{k}.py").write_text(v % {"good": str(GOOD)})
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setenv("PYTHONPATH", str(tmp_path))
    saved = dict(core.MODULES)

    def use(**kw):  # use(arith="src_good", cxpass="cx_bad", ...)
        core.MODULES.clear()
        for name, mod in kw.items():
            key = name if name == "verify" else f"contrib.{name}"
            core.MODULES[key] = f"qv04fake.{mod}" if mod else "qv04fake.missing_module"
    yield use
    core.MODULES.clear(); core.MODULES.update(saved)


def run_cli(args):
    argv = sys.argv
    sys.argv = ["qoptct"] + [str(a) for a in args]
    try:
        qoptct.main()
    finally:
        sys.argv = argv


# -- module isolation -----------------------------------------------------------------------------------------
def test_call_contrib_isolation(fake):
    fake(arith="src_raise", bandc="src_hang", gf2m="src_exit", cxpass=None)
    r, err = core.call_contrib("contrib.arith", "candidates", ("x.qc", False), seconds=5)
    assert r is None and "boom" in err
    t = time.time()
    r, err = core.call_contrib("contrib.bandc", "candidates", ("x.qc", False), seconds=1, grace=2)
    assert r is None and "timeout" in err and time.time() - t < 20
    r, err = core.call_contrib("contrib.gf2m", "candidates", ("x.qc", False), seconds=5)
    assert r is None and "exit code 3" in err
    assert not core.contrib_available("contrib.cxpass")
    assert core.call_contrib("contrib.cxpass", "improve", ("",), seconds=5) == (None, "not installed")


def test_qasm_stats_subset():
    assert core.qasm_stats(GOOD.read_text()) == (7, 8, 10, 2)
    assert core.qasm_stats("OPENQASM 3.0;\nqubit[2] q;\nccx q[0], q[1], q[2];\n") is None
    assert core.qasm_stats("OPENQASM 3.0;\nqubit[2] q;\nrz(0.1) q[0];\n") is None
    assert core.qasm_stats("OPENQASM 3.0;\nqubit[2] q;\ncx q[0], q[2];\n") is None  # index out of range
    u = core.qasm_unitary("OPENQASM 3.0;\nqubit[2] q;\nh q[0];\ncx q[0], q[1];\nt q[1];\n")
    assert u == (2, [("h", (0,)), ("cx", (0, 1)), ("t", (1,))])


def test_candidates_are_rechecked(fake):
    fake(arith="src_good", bandc="src_bad", gf2m="src_junk")
    log = []
    got = core.contrib_candidates(str(BENCH / "tof_3.qc"), False, ("arith", "bandc", "gf2m"), core.Budget(120), log)
    assert [k for k, _, _ in got] == [(8, 10, 2)]
    assert got[0][2].startswith("arith: known good")
    text = " | ".join(log)
    assert "bandc: rejected candidate" in text and "gf2m: rejected candidate (unsupported qasm)" in text


# -- command line ------------------------------------------------------------------------------------------------
def test_cli_t_mode_fake_modules(fake, tmp_path):
    """Good source wins, corrupted/raising/hanging sources and a wrong CX post-pass are rejected."""
    fake(arith="src_good", bandc="src_hang", gf2m="src_raise", cxpass="cx_bad", verify=None)
    out, rep = tmp_path / "o.qasm", tmp_path / "r.json"
    run_cli(["optimize", BENCH / "tof_3.qc", "-o", out, "--mode", "t", "--sel-sec", 3, "--tgad-sec", 0,
             "--and-variants", 0, "--contrib-sec", 3, "--report", rep])
    r = json.loads(rep.read_text())
    fin = r["final"]
    assert fin["format"] == "qasm" and fin["verified"] and (fin["t"], fin["cx"]) == (8, 10)
    assert fin["how"].startswith("arith: known good")
    log = " | ".join(r["modules"]["log"])
    assert "bandc: timeout" in log and "gf2m: RuntimeError('boom')" in log and "cxpass:" in log
    assert core.qasm_stats(out.read_text())[1:3] == (8, 10)
    assert verify_qasm.check(str(out), str(BENCH / "tof_3.qc"), trials=10) > 1 - 1e-9


@pytest.mark.slow
def test_cli_no_contrib_and_disable_flags(fake, tmp_path):
    fake(arith="src_bad", bandc="src_bad", gf2m="src_bad", cxpass="cx_bad")
    rep = tmp_path / "r.json"
    run_cli(["optimize", BENCH / "tof_3.qc", "--mode", "cx", "--no-contrib", "--report", rep])
    r = json.loads(rep.read_text())
    assert r["modules"]["sources"] == [] and r["modules"]["cxpass"] is False and r["final"]["verified"]
    run_cli(["optimize", BENCH / "tof_3.qc", "--mode", "t", "--sel-sec", 3, "--tgad-sec", 0, "--and-variants", 0,
             "--no-arith", "--no-gf2m", "--no-cxpass", "--report", rep])
    r = json.loads(rep.read_text())
    assert r["modules"]["sources"] == ["bandc"] and r["final"]["verified"]
    assert "bandc: 0 verified of 1" in " | ".join(r["modules"]["log"])


def test_cli_cx_mode_wrong_cxpass_rejected(fake, tmp_path):
    fake(cxpass="cx_bad")
    out, rep = tmp_path / "o.qc", tmp_path / "r.json"
    run_cli(["optimize", BENCH / "tof_3.qc", "-o", out, "--report", rep])
    r = json.loads(rep.read_text())
    assert "cxpass" not in r["method"] and r["final"]["verified"]
    assert any(l.startswith("cxpass:") for l in r["modules"]["log"])
    names, hl = qoptct.load_qc(str(out))
    import evo_lib
    ref = evo_lib.lower_std(qoptct.load_qc(str(BENCH / "tof_3.qc"))[1])
    assert qoptct.prove(evo_lib.lower_std(hl), ref, len(names)) is True


def test_cxpass_unitary_z_cz_output_is_lowered(fake, tmp_path):
    """A post-pass result with z / cz gates (ham15-low in cx mode) must be lowered before the .qc writer sees it."""
    fake(cxpass="cx_zcz")
    import evo_lib
    names, hl = qoptct.load_qc(str(BENCH / "tof_3.qc"))
    ref = evo_lib.lower_std(hl)
    best = ref + [("cx", (0, 1)), ("cx", (0, 1))]
    log = []
    c = core.cxpass_unitary(best, ref, len(names), core.Budget(120), log, path=str(BENCH / "tof_3.qc"))
    assert c is not None, log
    assert {g for g, _ in c} <= {"h", "x", "s", "sdg", "t", "tdg", "cx"}
    assert evo_lib.cx(c) == evo_lib.cx(ref)
    core.save_qc(str(tmp_path / "o.qc"), names, c)
    assert qoptct.prove(evo_lib.lower_std(qoptct.load_qc(str(tmp_path / "o.qc"))[1]), ref, len(names)) is True
    zc = [("z", (0,)), ("cz", (0, 1)), ("h", (1,))]  # the .qc writer itself also accepts z / cz
    core.save_qc(str(tmp_path / "z.qc"), names, zc)
    assert qoptct.load_qc(str(tmp_path / "z.qc"))[1] == zc


def test_time_limit_writes_checked_result(tmp_path):
    """--time-limit: the search stages are skipped or cut, the proven unitary result is written, nothing hangs."""
    out, rep = tmp_path / "o.qasm", tmp_path / "r.json"
    t = time.time()
    run_cli(["optimize", BENCH / "barenco_tof_4.qc", "-o", out, "--mode", "t", "--time-limit", 2, "--report", rep])
    assert time.time() - t < 60
    r = json.loads(rep.read_text())
    assert r["search_seconds"] <= 2.5 and r["final"]["format"] == "qc"
    # the proven portfolio result, or (if the limit already cut the portfolio) the input itself
    assert r["final"]["verified"] or r["method"] == "input (not proven)"
    assert out.read_text().startswith("OPENQASM 3.0;")  # unitary result written as OpenQASM for -o *.qasm
    assert any("time limit" in l for l in r["modules"]["log"])
    assert core._END is None and core.time_left() == float("inf")  # the limit belongs to that run only


def test_and_compiler_deadline():
    """--and-sec: an expired deadline still compiles the unshuffled variant, but no shuffled ones."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("c2_run", core.ROOT / "contrib" / "c2" / "run.py")
    c2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c2)
    path = str(BENCH / "tof_3.qc")
    t = time.time()
    first, _ = c2.candidates(path, 1, None)
    late, _ = c2.candidates(path, 50, None, deadline=time.time() - 1)
    assert late == first and late and time.time() - t < 60


# -- exact prover ----------------------------------------------------------------------------------------------
def _verify_ready():
    if not core.contrib_available("verify"):
        return False
    try:
        import qoptct.verify as v
        return callable(getattr(v, "prove", None))
    except Exception:
        return False


@pytest.mark.skipif(not _verify_ready(), reason="qoptct.verify not installed")
def test_prove_subcommand(tmp_path):
    exe = [sys.executable, "-c", "import sys; from qoptct.qopt_core import main; sys.argv[0] = 'qoptct'; main()"]
    r = subprocess.run(exe + ["prove", str(GOOD), str(BENCH / "tof_3.qc"), "--budget", "60"],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout + r.stderr
    assert json.loads(r.stdout.strip().splitlines()[-1])["status"] == "PROVED"
    bad = tmp_path / "bad.qasm"
    bad.write_text(GOOD.read_text().replace("\nt q", "\ntdg q", 1))
    r = subprocess.run(exe + ["prove", str(bad), str(BENCH / "tof_3.qc"), "--budget", "60"],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode in (1, 2) and json.loads(r.stdout.strip().splitlines()[-1])["status"] != "PROVED"


@pytest.mark.slow
@pytest.mark.skipif(not _verify_ready(), reason="qoptct.verify not installed")
def test_optimize_with_prove(tmp_path):
    out, rep = tmp_path / "o.qasm", tmp_path / "r.json"
    run_cli(["optimize", BENCH / "tof_3.qc", "-o", out, "--mode", "t", "--sel-sec", 3, "--tgad-sec", 0,
             "--prove", "--no-cxpass", "--report", rep])
    fin = json.loads(rep.read_text())["final"]
    assert fin["verified"] and fin["proof"]["status"] in ("PROVED", "UNKNOWN")
    assert fin["t"] <= 8


def _drop_first_conditional(text):
    lines = text.splitlines()
    i = next(k for k, L in enumerate(lines) if L.startswith("if ("))
    return "\n".join(lines[:i] + lines[i + 1:]) + "\n"


@pytest.mark.skipif(not _verify_ready(), reason="qoptct.verify not installed")
def test_strong_check_when_proof_unknown(fake):
    """Proof UNKNOWN: a candidate is accepted only after the strong sampled check; a dropped conditional
    correction is rejected."""
    fake(verify="ver_unknown")
    ref = str(BENCH / "tof_3.qc")
    good, bad = GOOD.read_text(), _drop_first_conditional(GOOD.read_text())
    ok, how = core._also_prove(good, ref, False, core.Budget(120), "simulator")
    assert ok and "prover UNKNOWN" in how and "sampled check dense" in how
    ok, how = core._also_prove(bad, ref, False, core.Budget(120), "simulator")
    assert not ok
    assert core.strong_check(good, ref, True)[0] and not core.strong_check(bad, ref, True)[0]


@pytest.mark.skipif(not _verify_ready(), reason="qoptct.verify not installed")
def test_prover_from_script_without_main_guard(tmp_path):
    """qoptct.verify.prove starts its own worker process; a plain script without a __main__ guard must work."""
    script = tmp_path / "plain.py"
    script.write_text("import qoptct.verify as V\n"
                      f"print(V.prove(open({str(GOOD)!r}).read(), {str(BENCH / 'tof_3.qc')!r}, False, 60)['status'])\n")
    r = subprocess.run([sys.executable, str(script)], capture_output=True, text=True, timeout=300, cwd=tmp_path)
    assert r.returncode == 0 and r.stdout.strip().splitlines()[-1] == "PROVED", r.stdout + r.stderr


@pytest.mark.skipif(not core.contrib_available("contrib.gf2m"), reason="qoptct.contrib.gf2m not installed")
def test_gf2m_tables_are_package_data():
    from qoptct.contrib import gf2m
    for name in ("decomps.json", "recipes.json", "known.json"):
        assert gf2m._data(name), f"{name} missing from the installed package"


def test_qasm_stats_xor_and_reset():
    t = ("OPENQASM 3.0;\nqubit[3] q;\nbit[2] c;\nh q[2];\nc[0] = measure q[2];\nreset q[2];\n"
         "c[1] = measure q[1];\nif ((c[0] ^ c[1])) { cz q[0], q[1]; }\nif (c[0] && c[1]) { x q[0]; }\n")
    assert core.qasm_stats(t) == (3, 0, 1, 2)


# -- Qiskit plugin ----------------------------------------------------------------------------------------------
def test_qiskit_pass_cxpass_option(fake):
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Operator
    from qiskit.transpiler import PassManager
    from qoptct.qiskit_qopt import QoptPass
    fake(cxpass="cx_bad")  # a wrong post-pass must never reach the circuit
    qc = QuantumCircuit(5); qc.ccx(0, 1, 4); qc.ccx(2, 4, 3); qc.ccx(0, 1, 4)
    p = QoptPass(cxpass=True); out = PassManager([p]).run(qc)
    assert p.last["verified"] and "cxpass" not in (p.last["method"] or "")
    assert Operator(out) == Operator(qc)


def test_qiskit_transpile_plugin_still_works():
    from qiskit import QuantumCircuit, transpile
    from qiskit.quantum_info import Operator
    qc = QuantumCircuit(5); qc.ccx(0, 1, 4); qc.ccx(2, 4, 3); qc.ccx(0, 1, 4)
    out = transpile(qc, basis_gates=["h", "x", "z", "s", "sdg", "t", "tdg", "cx"], optimization_method="qoptct",
                    init_method="qoptct", seed_transpiler=1)
    assert Operator(out).equiv(Operator(qc))


# -- real modules (when installed) ------------------------------------------------------------------------------
REAL = [s for s in core.SOURCES if core.contrib_available(f"contrib.{s}")]


@pytest.mark.slow
@pytest.mark.skipif(not REAL, reason="no candidate source module installed")
@pytest.mark.parametrize("clean", [False, True])
def test_real_sources_never_worse_tof3(tmp_path, clean):
    """End to end on tof_3 (t mode): with the installed modules the result is verified and not worse than
    without them (the v0.3 pipeline)."""
    res = {}
    for tag, extra in (("v03", ["--no-contrib"]), ("v04", [])):
        rep = tmp_path / f"{tag}.json"
        run_cli(["optimize", BENCH / "tof_3.qc", "-o", tmp_path / f"{tag}.qasm", "--mode", "t", "--sel-sec", 3,
                 "--tgad-sec", 0, "--contrib-sec", 30, "--report", rep] + (["--clean"] if clean else []) + extra)
        res[tag] = json.loads(rep.read_text())["final"]
        assert res[tag]["verified"]
    assert (res["v04"]["t"], res["v04"]["cx"]) <= (res["v03"]["t"], res["v03"]["cx"])
