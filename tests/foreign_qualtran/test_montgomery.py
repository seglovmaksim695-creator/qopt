"""Montgomery multiplication of Qualtran: positive proofs (small n) and the regression set that must NOT be proved.

Fast tests use n = 16 (seconds).  Tests with the mark ``slow`` use n = 64 or n = 256 (minutes; run heavy ones under
``lab/tools/one``).  Qualtran is optional: without it everything that needs it is skipped.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

qualtran = pytest.importorskip("qualtran")
pytest.importorskip("z3")

from qoptct.foreign.qualtran import montgomery as M


@pytest.fixture(scope="module", autouse=True)
def _cache(tmp_path_factory):
    mp = pytest.MonkeyPatch()
    mp.setenv("QOPT_CACHE", os.environ.get("QOPT_CACHE_TESTS") or str(tmp_path_factory.mktemp("qopt_cache")))
    yield
    mp.undo()


def _pass(v):
    return v["verdict"] == "PROVED" and not v["reasons"]


@pytest.mark.parametrize("w", [1, 2, 4])
def test_proved_n16(w):
    v = M.prove_montgomery(16, w)
    assert _pass(v), v["reasons"]
    assert v["patched"] is True and v["qualtran"]["patch_kind"] == "c"
    assert v["qualtran"]["patch_files_sha256"], "the patch hashes must be recorded"
    assert v["versions"]["z3"] and v["versions"]["qualtran"] == "0.7.0"
    assert v["glue"]["ok"] and v["bookkeeping"]["ok"] and v["telescoping"] is True
    assert v["window_checks"]["all_windows_ok"] and v["window_checks"]["same_hash"] and v["window_checks"]["tiles_bad"] == {}
    assert v["prefix_X_reduced"]["ok"]
    labels = {s["label"] for s in v["segments"]}
    assert {"lt", "add"} <= labels
    g = v["generic_window"]
    assert g["ok"] and [r["kind"] for r in g["subnodes"]] == ["cadd"] * w + ["xor", "qaa"]
    assert all(r["status"] == "PROVED" and r["domain"] == "PASS" for r in g["subnodes"])
    for s in v["segments"]:
        assert len(s["qasm_sha256"]) == 64 and s["status"] == "PROVED"
    assert len(v["program"]["ops_sha256"]) == 64 and len(v["program"]["generic_window_qasm_sha256"]) == 64
    assert v["end_to_end"]["mismatches"] == 0 and v["end_to_end"]["cases_with_reduction"] > 0


def test_installed_qualtran_untouched():
    from qoptct.foreign.qualtran.patches import _installed_root, _MANIFEST
    from hashlib import sha256
    M.prove_montgomery(16, 1)
    f = _installed_root() / _MANIFEST["files"]["mod_multiplication_c.py"]["target"]
    assert sha256(f.read_bytes()).hexdigest() == _MANIFEST["files"]["mod_multiplication_c.py"]["sha256_original"]


# ---------------------------------------------------------------- regression: everything below must be NOT_PROVED
@pytest.mark.parametrize("tag", ["H1", "H2", "M1", "M1b", "M2", "M3", "M4"])
def test_mutant_not_proved_n16_w1(tag):
    v = M.prove_montgomery(16, 1, mutant=tag)
    assert v["verdict"] == "NOT_PROVED" and v["reasons"], tag
    assert not any("stage failed" in r for r in v["reasons"]), v["reasons"]   # refused by the checks, not by a crash


@pytest.mark.parametrize("tag", ["H1", "H2", "M1", "M2", "M3", "M4"])
def test_mutant_not_proved_n16_w2(tag):
    v = M.prove_montgomery(16, 2, mutant=tag)
    assert v["verdict"] == "NOT_PROVED" and v["reasons"], tag
    assert not any("stage failed" in r for r in v["reasons"]), v["reasons"]   # refused by the checks, not by a crash


def test_mutant_M5_prefix_w4():
    v = M.prove_montgomery(16, 4, mutant="M5")
    assert v["verdict"] == "NOT_PROVED"
    assert any("X(reduced)" in r for r in v["reasons"]), v["reasons"]


def test_mutant_reasons_are_the_coverage_checks():
    # H2 / M1: an operation that belongs to no sub-node; M2/M3/H1: a sub-node writes into a wire that is not fresh
    for tag, needle in (("H2", "coverage/freshness"), ("M2", "coverage/freshness"), ("M1", "coverage/freshness")):
        v = M.prove_montgomery(16, 1, mutant=tag)
        assert any(needle in r for r in v["reasons"]), (tag, v["reasons"])


def _canaries(n, w, count, seed):
    from qoptct.foreign.qualtran.patches import patched_env
    root = str(Path(M.__file__).resolve().parents[4])
    env = patched_env("c")
    env["PYTHONPATH"] = root + os.pathsep + env["PYTHONPATH"]
    out = Path(os.environ["QOPT_CACHE"]) / ("canaries_%d_%d_%d.json" % (n, w, seed))
    r = subprocess.run([sys.executable, "-B", "-m", "qoptct.foreign.qualtran.montgomery.regression.canaries", "--n", str(n), "--w", str(w),
                        "--count", str(count), "--seed", str(seed), "--out", str(out)], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-500:]
    return json.loads(out.read_text())


@pytest.mark.parametrize("w", [1, 2])
def test_canaries_sample_n16(w):
    res = _canaries(16, w, 12, 11)
    assert len(res["canaries"]) == 12 and res["targeted"] >= 6
    assert res["false_proved"] == [], res["false_proved"]


def test_original_qualtran_not_proved_with_witness():
    v = M.prove_montgomery(16, 1, patched=False)
    assert v["verdict"] == "NOT_PROVED" and v["patched"] is False and v["qualtran"]["patch_kind"] is None
    wt = v["witness"]
    assert wt and wt["T_N_ge_p"] is True, wt           # the witness is in Bad = {T_N >= p}
    assert wt["got"]["reduced"] == 1 and wt["got"]["target"] != wt["expected"]["target"]
    st = {s["label"]: s["status"] for s in v["segments"]}
    assert st["win0"] == "PROVED" and st["lt"] == "PROVED" and st["add"] != "PROVED"   # the error is only in the "+p" block


@pytest.mark.parametrize("n,w", [(6, 1), (6, 2)])
def test_original_wrong_exactly_on_bad(n, w):
    p = M.prime_for(n)
    orig = M.exhaustive_small(n, w, p, patched=False)
    assert orig["mismatch"] > 0 and orig["mismatch"] == orig["mismatch_with_T_N_ge_p"]
    assert orig["bad_phase"] == 0 and orig["bad_rest"] == 0
    fixed = M.exhaustive_small(n, w, p, patched=True)
    assert fixed["mismatch"] == 0 and fixed["total"] == orig["total"] and fixed["bad_phase"] == 0 and fixed["bad_rest"] == 0


def test_cli_writes_result(tmp_path):
    root = str(Path(M.__file__).resolve().parents[4])
    out = tmp_path / "r.json"
    env = dict(os.environ, PYTHONPATH=root + os.pathsep + os.environ.get("PYTHONPATH", ""))
    r = subprocess.run([sys.executable, "-B", "-m", "qoptct.foreign.qualtran.montgomery", "--n", "16", "--w", "1", "--out", str(out)],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-300:] + r.stderr[-300:]
    assert json.loads(out.read_text())["verdict"] == "PROVED"


# ---------------------------------------------------------------- slow
@pytest.mark.slow
@pytest.mark.parametrize("tag", ["H1", "H2", "M1", "M1b", "M2", "M3", "M4"])
def test_mutant_not_proved_n64_w1(tag):
    v = M.prove_montgomery(64, 1, mutant=tag)
    assert v["verdict"] == "NOT_PROVED" and v["reasons"], tag
    assert not any("stage failed" in r for r in v["reasons"]), v["reasons"]   # refused by the checks, not by a crash


@pytest.mark.slow
def test_proved_n64_w4():
    v = M.prove_montgomery(64, 4)
    assert _pass(v), v["reasons"]


@pytest.mark.slow
def test_proved_n256_w1():
    v = M.prove_montgomery(256, 1)
    assert _pass(v), v["reasons"]
    assert v["p"] == 2 ** 255 - 19
