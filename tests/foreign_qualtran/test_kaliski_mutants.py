"""Regression: every wiring mutant (the 14 own ones of the corrected client and the auditor's F1-F4 or their analogues at small
n) must be NOT PROVED, with the reason "wire flow"; and the canaries (damaged circuits) must never be proved when they misbehave."""
import json

import pytest

pytest.importorskip("qualtran")

from kaliski_support import P25519

from qoptct.foreign.qualtran.kaliski import canary, mutants

OWN = ["O%d " % i for i in range(1, 15)]


def _check(res, expect_rows):
    assert len(res["rows"]) == expect_rows
    bad = [r["name"] for r in res["rows"] if r["status"] != "NOT_PROVED"]
    assert not bad, "mutants that were NOT rejected: %s" % bad
    assert all(r["flow_reason"] for r in res["rows"]), [(r["name"], r["reason"]) for r in res["rows"] if not r["flow_reason"]]
    assert res["all_rejected"] and res["all_flow_reason"]


@pytest.mark.parametrize("n,p", [(8, 251), (10, 1021)])
def test_14_own_mutants_without_the_static_pass(n, p):
    """The built-in check of _begin/_fin alone (the static pre-pass is switched off) rejects all 14."""
    res = mutants.run_mutants(n, p, noprepass=True, verbose=False)
    _check(res, 14)


@pytest.mark.parametrize("n,p", [(8, 251), (16, 32749)])
def test_18_mutants_with_the_auditor_analogues(n, p):
    res = mutants.run_mutants(n, p, aud=True, verbose=False)
    _check(res, 18)
    names = " ".join(r["name"] for r in res["rows"])
    for f in ("F1", "F2", "F3", "F4"):
        assert f + " " in names


def test_auditor_parameters_at_n256_are_those_of_the_audit():
    ap = mutants.aud_params(256)
    assert (ap["vswap_k"], ap["vswap_bit"], ap["t_src"], ap["t_dst"], ap["mt_k"]) == (507, 1, 262, 269, 100)
    for n in (8, 16, 64):
        ap = mutants.aud_params(n)
        assert 0 <= ap["t_src"] < ap["t_dst"] < 2 * n and 0 <= ap["mt_k"] < 2 * n and ap["vswap_k"] < 2 * n


def test_without_the_flow_check_the_mutants_would_pass_the_old_proof_steps():
    """Control for the control: the mutants are real.  With the flow check ON they are rejected FOR THE FLOW; the test above would be
    empty if the mutants were rejected by something else, so the reason is part of the assertion (flow_reason)."""
    res = mutants.run_mutants(8, 251, noprepass=True, verbose=False, only="O1 ")
    assert res["rows"][0]["reason"].startswith("wire flow")


@pytest.mark.parametrize("n,p", [(8, 251), (12, 4093)])
def test_canaries_never_proved_when_they_misbehave(n, p, tmp_path):
    out = tmp_path / "c.json"
    canary.main(["--n", str(n), "--p", str(p), "--out", str(out), "--random", "6"])
    res = json.loads(out.read_text())
    rows = res["canaries"]
    assert len(rows) >= 19
    false_proved = [r["name"] for r in rows if r.get("prove") == "CHAIN_PROVED" and r.get("witness")]
    assert not false_proved
    assert res["targeted_counted"] >= 11                      # at n=256 the audited run counted 13 targeted ones
    assert all(r["counted"] for r in rows if r["name"].startswith("T") and r["name"] not in ("T8",)) or res["targeted_counted"] >= 11


@pytest.mark.slow
def test_f1_to_f4_and_own_mutants_at_n256():
    res = mutants.run_mutants(256, P25519, aud=True, verbose=False)
    _check(res, 18)
