"""Release checks: independent verifier, model handling, cache guards, no absolute paths."""
import json
import pathlib
import re
import shutil

import pytest

import qopt
from qopt import verify_qasm

ROOT = pathlib.Path(__file__).resolve().parents[1]
BENCH = ROOT / "bench_qc"
PKG = pathlib.Path(qopt.__file__).resolve().parent

# 8 T / 10 CX, 2 fresh ancillas: tof_3 compiled with logical ANDs and measurement-based uncomputation
GOOD = sorted((ROOT / "results" / "best" / "strict").glob("tof_3.qasm"))


def test_no_absolute_project_paths():
    for f in PKG.rglob("*.py"):
        assert "/Users/" not in f.read_text(), f


def test_version():
    assert qopt.__version__ == "0.2.0"


@pytest.mark.skipif(not GOOD, reason="results/best not present")
def test_verify_strict_accepts_and_canary_rejects(tmp_path):
    good = GOOD[0]
    assert verify_qasm.check(str(good), str(BENCH / "tof_3.qc"), trials=10) > 1 - 1e-9
    lines = good.read_text().splitlines()
    # remove the first classically controlled correction -> must fail
    i = next(k for k, l in enumerate(lines) if l.startswith("if ("))
    bad = tmp_path / "bad.qasm"; bad.write_text("\n".join(lines[:i] + lines[i + 1:]) + "\n")
    assert verify_qasm.check(str(bad), str(BENCH / "tof_3.qc"), trials=10) < 0.999
    # replace one T by T-dagger -> must fail
    j = next(k for k, l in enumerate(lines) if re.match(r"t q\[", l))
    bad2 = tmp_path / "bad2.qasm"; bad2.write_text("\n".join(lines[:j] + [lines[j].replace("t q", "tdg q", 1)] + lines[j + 1:]) + "\n")
    assert verify_qasm.check(str(bad2), str(BENCH / "tof_3.qc"), trials=10) < 0.999


def test_unitary_qasm_without_bits(tmp_path):
    names, hl = qopt.load_qc(str(BENCH / "tof_3.qc"))
    import evo_lib
    ref = evo_lib.lower_std(hl)
    f = tmp_path / "u.qasm"
    f.write_text("OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[5] q;\n" +
                 "".join(f"{g} " + ", ".join(f"q[{x}]" for x in q) + ";\n" for g, q in ref))
    assert verify_qasm.check(str(f), str(BENCH / "tof_3.qc"), trials=5) > 1 - 1e-9


def test_ours_cache_guard(tmp_path):
    cached = PKG / "contrib" / "eng1" / "ours" / "tof_3.json"
    assert qopt.ours_file(str(BENCH / "tof_3.qc"), cached) == cached
    # same name, different content -> the shipped cache must not be used
    other = tmp_path / "tof_3.qc"
    other.write_text((BENCH / "tof_4.qc").read_text())
    assert qopt.ours_file(str(other), cached) is None


def test_clean_construction_guard(tmp_path):
    assert qopt.bench_sha_ok(str(BENCH / "mod_red_21.qc"), "mod_red_21")
    other = tmp_path / "mod_red_21.qc"
    shutil.copy(BENCH / "mod_red_21.qc", other)
    other.write_text(other.read_text() + "# modified copy\n")
    assert not qopt.bench_sha_ok(str(other), "mod_red_21")


def test_balanced_mode_tof3(tmp_path):
    out = tmp_path / "o.qasm"; rep = tmp_path / "r.json"
    import sys
    argv = sys.argv
    sys.argv = ["qopt", "optimize", str(BENCH / "tof_3.qc"), "-o", str(out), "--mode", "balanced", "--sel-sec", "20",
                "--report", str(rep)]
    try:
        qopt.main()
    finally:
        sys.argv = argv
    r = json.loads(rep.read_text())
    assert r["verified"] and r["cx_out"] <= 12 and r["t_out"] <= 15
    assert r["and"]["chosen"] and r["and"]["t"] <= 8
    assert verify_qasm.check(str(out), str(BENCH / "tof_3.qc"), trials=10) > 1 - 1e-9
