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


def test_cli_exit_codes(tmp_path, capsys):
    q, c = _adder_files(tmp_path, 4)
    m, _ = _adder_files(tmp_path, 4, mutate=True)
    assert kprove_main([q, c, "--clean", "--engine", "bdd"]) == 0
    assert kprove_main([m, c, "--clean", "--engine", "bdd"]) == 1
    assert kprove_main([str(tmp_path / "missing.qasm"), c]) == 3
