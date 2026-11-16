"""Tests for qoptct.verify.kprove (BDD path-sum prover).

* regression of the garbage-collection bug fixed in 0.5.0: input / outcome literals were cached as BDD edges and
  pointed to foreign nodes after a collection -> false REFUTED. Every case runs with a collection at every budget check
  (KPROVE_FORCE_GC=1) and without one; the verdicts must agree with the known truth;
* a fixed subset of the validation set (tests/data/kprove: correct circuits and mutants of small benchmarks, truth from
  the dense simulator): no false PROVED, no false REFUTED;
* memory guard: a tiny KPROVE_MAX_MB gives UNKNOWN (never PROVED).
"""
import json
import os

import pytest

from qoptct.verify import arith_gens as G
from qoptct.verify.kprove import kprove, main as kprove_main

DATA = os.path.join(os.path.dirname(__file__), "data", "kprove")


def _adder_files(tmp_path, n, mutate=False):
    q = tmp_path / ("add%d%s.qasm" % (n, "_m" if mutate else ""))
    txt = G.gidney_adder(n)
    if mutate:
        L = txt.splitlines()
        i = [k for k, l in enumerate(L) if l.startswith("if")][n // 2]
        del L[i]
        txt = "\n".join(L) + "\n"
    q.write_text(txt)
    c = tmp_path / ("add%d.qc" % n)
    c.write_text(G.adder_qc(n))
    return str(q), str(c)


def _spec(n):
    def spec(S):
        a = S.reg(range(n))
        b = S.reg(range(n, 2 * n))
        s = a + b
        o = S.out(range(n, 2 * n), s)
        o[2 * n] = S.B.XOR(S.x(2 * n), s.bits[n])
        return o
    return spec


@pytest.mark.parametrize("force_gc", ["0", "1"])
@pytest.mark.parametrize("n", [6, 24])
def test_gc_regression(tmp_path, monkeypatch, n, force_gc):
    monkeypatch.setenv("KPROVE_FORCE_GC", force_gc)
    q, c = _adder_files(tmp_path, n)
    m, _ = _adder_files(tmp_path, n, mutate=True)
    assert kprove(q, c, clean=True, engine="bdd", sample=False)["status"] == "PROVED"
    kw = dict(spec=_spec(n), n=3 * n + 1, clean_wires=range(2 * n + 1, 3 * n + 1))
    assert kprove(q, **kw)["status"] == "PROVED"
    assert kprove(m, **kw)["status"] == "REFUTED"


def _cases():
    with open(os.path.join(DATA, "cases.json")) as f:
        return json.load(f)


@pytest.mark.parametrize("case", _cases(), ids=lambda c: c["qasm"])
def test_validation_subset(case, monkeypatch):
    monkeypatch.setenv("KPROVE_FORCE_GC", "3")      # collect often: the bug class fixed in 0.5.0
    r = kprove(os.path.join(DATA, case["qasm"]), os.path.join(DATA, case["qc"]), clean=case["clean"],
               budget_sec=60, engine="bdd", sample=False)
    if case["truth"] == "ok":
        assert r["status"] != "REFUTED", r
        assert r["status"] == "PROVED", r
    else:
        assert r["status"] != "PROVED", r
        assert r["status"] == "REFUTED", r


def test_large_adder_proved(tmp_path):
    q, c = _adder_files(tmp_path, 256)
    r = kprove(q, c, clean=True, engine="bdd", sample=False, budget_sec=120)
    assert r["status"] == "PROVED", r


def test_memory_guard_unknown(tmp_path, monkeypatch):
    monkeypatch.setenv("KPROVE_MAX_MB", "1")         # any real process is above 1 MB: must stop, never PROVED
    q, c = _adder_files(tmp_path, 64)
    r = kprove(q, c, clean=True, engine="bdd", sample=False)
    assert r["status"] == "UNKNOWN" and "memory" in r["detail"], r


def test_memory_budget_is_growth(tmp_path, monkeypatch):
    """max_mb counts the growth of the RSS during the call, not the earlier peak of the process."""
    from qoptct.verify.kprove import MemGuard, current_rss_mb
    monkeypatch.delenv("KPROVE_MAX_MB", raising=False)
    ballast = bytearray(150 * 1000 * 1000)            # the process is now >= 150 MB above its start
    ballast[::4096] = b"\1" * len(ballast[::4096])
    q, c = _adder_files(tmp_path, 64)
    r = kprove(q, c, clean=True, engine="bdd", sample=False, max_mb=100)
    assert r["status"] == "PROVED", r
    assert r["mem"]["max_mb"] == 100 and r["mem"]["peak_growth_mb"] < 100
    if current_rss_mb() is not None:
        assert r["mem"]["base_mb"] >= 150 and r["mem"]["measure"] == "current RSS"
    del ballast
    if current_rss_mb() is not None:
        g = MemGuard(5)
        big = bytearray(60 * 1000 * 1000)
        big[::4096] = b"\1" * len(big[::4096])
        with pytest.raises(Exception, match="memory"):
            g.check()
        del big
    assert MemGuard().max_mb == 900
    monkeypatch.setenv("KPROVE_MAX_MB", "1234")
    assert MemGuard().max_mb == 1234 and MemGuard(77).max_mb == 77


def test_memory_guard_param_unknown(tmp_path):
    q, c = _adder_files(tmp_path, 64)
    r = kprove(q, c, clean=True, engine="bdd", sample=False, max_mb=1)
    assert r["status"] == "UNKNOWN" and "memory" in r["detail"], r


def test_exact2_engine(tmp_path):
    """engine 'exact2' (degree <= 2) on a tiny bilinear circuit: c ^= a & b via one logical AND + erasure."""
    and_body = ["h q[2];", "t q[2];", "cx q[1], q[2];", "tdg q[2];", "cx q[0], q[2];", "t q[2];", "cx q[1], q[2];",
                "tdg q[2];", "h q[2];", "sdg q[2];"]
    qasm = ["OPENQASM 3;", 'include "stdgates.inc";', "qubit[4] q;", "bit[1] c;"] + and_body + [
        "cx q[2], q[3];", "h q[2];", "c[0] = measure q[2];", "if (c[0]) {cz q[0], q[1]; x q[2];}"]
    q = tmp_path / "and.qasm"
    q.write_text("\n".join(qasm) + "\n")
    ref = tmp_path / "and.qc"
    ref.write_text(".v a b w c\n.i a b c\nBEGIN\ntof a b c\nEND\n")
    assert kprove(str(q), str(ref), clean=True, engine="exact2")["status"] == "PROVED"
    r = kprove(str(q), str(ref), clean=True)                   # auto: bilinear reference -> exact2 first
    assert r["status"] == "PROVED" and r["method"] == "exact2", r
    q.write_text("\n".join(qasm[:-4] + qasm[-3:]) + "\n")      # the fan-out CX removed
    assert kprove(str(q), str(ref), clean=True, engine="exact2")["status"] == "REFUTED"
    assert kprove_main([str(q), str(ref), "--clean", "--engine", "exact2", "--max-mb", "500"]) == 1


def test_cli_exit_codes(tmp_path, capsys):
    q, c = _adder_files(tmp_path, 4)
    m, _ = _adder_files(tmp_path, 4, mutate=True)
    assert kprove_main([q, c, "--clean", "--engine", "bdd"]) == 0
    assert kprove_main([m, c, "--clean", "--engine", "bdd"]) == 1
    assert kprove_main([str(tmp_path / "missing.qasm"), c]) == 3
