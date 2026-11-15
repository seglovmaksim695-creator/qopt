"""Tests for qoptct.verify (exact prover, fast checks, c10 verify_big fixes, arith). Runs in well under 90 s.

Benchmarks are read from $QOPT_BENCH (see conftest.py); tests that need them are skipped
when the folder is missing. Every correct circuit is paired with "canaries" (corrupted copies) that must be rejected;
the truth for a canary is taken from the independent dense simulator, so harmless mutations are not miscounted.
"""
import os
import random
import re

import pytest

from qoptct.verify import prove, quick_check, check_dense, check_sparse, check_big, check_qcec, qcec_available
from qoptct.verify import main_prove
from qoptct.verify._parse import parse_qc
from qoptct.verify.pathsum import prove_exact
from qoptct.verify import arith as VA, arith_gens as G

BENCH = os.environ["QOPT_BENCH"]
STD = os.path.join(BENCH, "std")

TOF = lambda a, b, t: ["h q[%d];" % t, "cx q[%d], q[%d];" % (b, t), "tdg q[%d];" % t, "cx q[%d], q[%d];" % (a, t),
                       "t q[%d];" % t, "cx q[%d], q[%d];" % (b, t), "tdg q[%d];" % t, "cx q[%d], q[%d];" % (a, t),
                       "t q[%d];" % b, "t q[%d];" % t, "h q[%d];" % t, "cx q[%d], q[%d];" % (a, b), "t q[%d];" % a,
                       "tdg q[%d];" % b, "cx q[%d], q[%d];" % (a, b)]


def lower_qc(path):
    """Independent textbook lowering of a .qc file to Clifford+T OpenQASM 3 (None if a gate has > 2 controls)."""
    names, ins, gates = parse_qc(path)
    body = []
    for g, q in gates:
        if g in ("h", "x", "y", "z", "s", "sdg", "t", "tdg"):
            body.append("%s q[%d];" % (g, q[0]))
        elif g in ("cx", "cz"):
            body.append("%s q[%d], q[%d];" % (g, q[0], q[1]))
        elif g == "mcx" and len(q) == 3:
            body += TOF(*q)
        elif g == "mcz" and len(q) == 3:
            body += ["h q[%d];" % q[2]] + TOF(*q) + ["h q[%d];" % q[2]]
        else:
            return None
    return "OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[%d] q;\n" % len(names) + "\n".join(body) + "\n"


def canaries(text, k, seed):
    """Corrupted copies: drop/flip T, drop/swap CX, drop an `if`, change a condition bit."""
    rng = random.Random(seed)
    lines = text.splitlines()
    out = []
    for _ in range(200):
        if len(out) >= k:
            break
        L = list(lines)
        kind = rng.choice(["flip_t", "drop_cx", "swap_cx", "drop_if", "cond_bit", "drop_t"])
        idx = {"flip_t": r"(t|tdg) q", "drop_t": r"(t|tdg) q", "drop_cx": r"cx q", "swap_cx": r"cx q",
               "drop_if": r"if", "cond_bit": r"if"}[kind]
        ii = [i for i, l in enumerate(L) if re.match(idx, l)]
        if not ii:
            continue
        i = rng.choice(ii)
        if kind == "flip_t":
            L[i] = ("tdg" + L[i][1:]) if L[i].startswith("t ") else ("t" + L[i][3:])
        elif kind in ("drop_cx", "drop_if", "drop_t"):
            del L[i]
        elif kind == "swap_cx":
            m = re.match(r"cx q\[(\d+)\], q\[(\d+)\];", L[i])
            L[i] = "cx q[%s], q[%s];" % (m.group(2), m.group(1))
        else:
            meas = sorted({int(x) for l in lines for x in re.findall(r"c\[(\d+)\] = measure", l)})
            ks = re.findall(r"c\[(\d+)\]", L[i])
            if len(meas) < 2:
                continue
            old = rng.choice(ks)
            L[i] = L[i].replace("c[%s]" % old, "c[%d]" % rng.choice([k for k in meas if str(k) != old]), 1)
        t = "\n".join(L) + "\n"
        if t != text:
            out.append(t)
    return out


def _std(name):
    p = os.path.join(STD, name + ".qc")
    if not os.path.exists(p):
        pytest.skip("benchmarks not found (set QOPT_BENCH)")
    return p


@pytest.fixture
def adder(tmp_path):
    """Gidney adder with logical-AND gadgets and measurement-based uncomputation (clean model)."""
    n = 3
    q = tmp_path / "add.qasm"
    q.write_text(G.gidney_adder(n))
    r = tmp_path / "add.qc"
    r.write_text(G.adder_qc(n))
    return str(q), str(r)


# ------------------------------------------------------------------------------------------- exact prover

@pytest.mark.parametrize("name", ["tof_3", "barenco_tof_3", "mod5_4", "rc_adder_6"])
def test_prove_unitary_bench_and_canaries(name, tmp_path):
    qc = _std(name)
    text = lower_qc(qc)
    if text is None:
        pytest.skip("gate with more than two controls")
    r = prove(text, qc, budget_sec=60)
    assert r["status"] == "PROVED", r
    assert r["method"].startswith("pathsum"), r
    n_ref = 0
    for t in canaries(text, 6, seed=len(name)):
        truth = check_dense(t, qc)["ok"]
        r = prove(t, qc, budget_sec=60, use_qcec=False)
        assert r["status"] == ("PROVED" if truth else "REFUTED"), (r, truth)
        n_ref += not truth
    assert n_ref >= 3


def test_prove_measurements_feed_forward(adder):
    q, qc = adder
    r = prove(q, qc, clean=True, budget_sec=60)
    assert r["status"] == "PROVED", r
    # the same circuit is not correct in the strict model (carry wires are not clean there)
    assert prove(q, qc, clean=False, budget_sec=60)["status"] == "REFUTED"
    text = open(q).read()
    caught = 0
    for t in canaries(text, 8, seed=3):
        truth = check_dense(t, qc, clean=True)["ok"]
        r = prove(t, qc, clean=True, budget_sec=60)
        assert r["status"] == ("PROVED" if truth else "REFUTED"), (r, truth)
        caught += not truth
        if r["status"] == "REFUTED" and r.get("counterexample") and (r.get("proof") or {}).get("confirmed") is not None:
            assert r["proof"]["confirmed"]
    assert caught >= 5


def test_prove_bad_input_is_unknown(tmp_path):
    qc = _std("tof_3")
    assert prove("OPENQASM 3.0;\nqubit[5] q;\nrz(0.1) q[0];\n", qc)["status"] == "UNKNOWN"
    assert prove(str(tmp_path / "missing.qasm"), qc)["status"] == "UNKNOWN"
    assert prove("OPENQASM 3.0;\nqubit[1] q;\nh q[0];\n", qc)["status"] == "UNKNOWN"   # fewer qubits


def test_prove_time_budget_respected(adder):
    q, qc = adder
    r = prove(q, qc, clean=True, budget_sec=0.0001)
    assert r["status"] in ("PROVED", "UNKNOWN")      # never a wrong verdict, and it returns


# ------------------------------------------------------------------------------------------- fast checks

def test_quick_checks(adder):
    q, qc = adder
    assert check_dense(q, qc, clean=True)["ok"]
    assert check_sparse(q, qc, clean=True)["ok"]
    assert quick_check(q, qc, clean=True)["ok"]
    bad = open(q).read().replace("if (c[2]) {cz", "if (c[0]) {cz", 1)
    assert not check_sparse(bad, qc, clean=True)["ok"]
    assert not check_dense(bad, qc, clean=True)["ok"]


def test_check_big_and_sampled_norm(tmp_path):
    qc = _std("barenco_tof_3")
    text = lower_qc(qc)
    assert check_big(text, qc, x_trials=2, m_trials=1)["ok"]
    # cap=2 forces the sampled mode; the exact branch norm must still reject a T -> Tdg canary
    r = check_big(text, qc, x_trials=2, m_trials=1, cap=2)
    assert r["ok"] and r["sampled"] > 0, r
    bad = text.replace("tdg q", "t q", 1)
    assert not check_dense(bad, qc)["ok"]
    assert not check_big(bad, qc, x_trials=2, m_trials=1, cap=2)["ok"]


def test_check_big_many_hadamards(tmp_path):
    """2^(-nH/2) for thousands of H used to underflow (nan -> vacuous pass). Now: correct passes, canary fails."""
    ref = tmp_path / "id.qc"
    ref.write_text(".v a b\nBEGIN\nEND\n")
    body = "h q[0];\nh q[0];\n" * 1200
    good = "OPENQASM 3.0;\nqubit[2] q;\n" + body
    r = check_big(good, str(ref), x_trials=1, m_trials=1)
    assert r["ok"], r
    bad = "OPENQASM 3.0;\nqubit[2] q;\n" + body + "h q[0];\nt q[0];\nh q[0];\n"
    assert not check_big(bad, str(ref), x_trials=1, m_trials=1)["ok"]


@pytest.mark.skipif(not qcec_available(), reason="mqt.qcec not installed")
def test_qcec_optional():
    qc = _std("tof_3")
    text = lower_qc(qc)
    assert check_qcec(text, qc, budget_sec=30)["status"] == "PROVED"
    assert check_qcec(text.replace("tdg q", "t q", 1), qc, budget_sec=30)["status"] in ("REFUTED", "UNKNOWN")


# ------------------------------------------------------------------------------------------- arith

def test_arith_adder_all_methods(tmp_path):
    n = 4
    p = tmp_path / "a.qasm"
    p.write_text(G.gidney_adder(n))
    q = tmp_path / "a.qc"
    q.write_text(G.adder_qc(n))
    for m in ("sim", "pathsum"):
        r = VA.check(str(p), G.gidney_adder_ref(n), method=m, m_trials=4)
        assert r["ok"] and r["proved_all_x"], r
    assert VA.check_proof(str(p), VA.spec_from_qc(str(q), clean=True), m_trials=4, clean=True)["ok"]
    assert VA.check(str(p), VA.ref_from_qc(str(q), clean=True), anc="zero_clean", m_trials=4)["ok"]


def test_arith_clean_outputs_are_checked(tmp_path):
    """Regression (v0.4): in the clean model check_sim used to compare only the input wires."""
    n = 4
    q = tmp_path / "a.qc"
    q.write_text(G.adder_qc(n))
    lines = G.gidney_adder(n).splitlines()
    i = max(k for k, l in enumerate(lines) if l.startswith("cx q[%d], q[%d]" % (0, n)))   # last sum cx a0 -> b0
    del lines[i]
    p = tmp_path / "m.qasm"
    p.write_text("\n".join(lines) + "\n")
    ref = VA.ref_from_qc(str(q), clean=True)
    for m in ("sim", "pathsum"):
        assert not VA.check(str(p), ref, method=m, anc="zero_clean", m_trials=4)["ok"]
    assert prove(str(p), str(q), clean=True, budget_sec=60)["status"] == "REFUTED"


def test_arith_large(tmp_path):
    n = 256
    p = tmp_path / "a.qasm"
    p.write_text(G.gidney_adder(n))
    assert VA.check(str(p), G.gidney_adder_ref(n), method="sim", B=16, m_trials=2)["ok"]
    bad = p.read_text().replace("tdg q", "t q", 1)
    p.write_text(bad)
    assert not VA.check(str(p), G.gidney_adder_ref(n), method="sim", B=16, m_trials=2)["ok"]


# ------------------------------------------------------------------------------------------- CLI

def test_cli(adder, tmp_path, capsys):
    q, qc = adder
    assert main_prove([q, qc, "--clean", "--budget", "60"]) == 0
    bad = tmp_path / "bad.qasm"
    bad.write_text(open(q).read().replace("if (c[2]) {cz", "if (c[0]) {cz", 1))
    assert main_prove([str(bad), qc, "--clean", "--json"]) == 1
    assert '"REFUTED"' in capsys.readouterr().out
    assert main_prove([]) == 3
