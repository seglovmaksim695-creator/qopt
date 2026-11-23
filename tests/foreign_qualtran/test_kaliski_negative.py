"""Negative tests of the verdict: damaged or foreign evidence is NOT PROVED (empty or foreign lemma file, substituted query,
no wire-flow check, another client, dirty simulation ...)."""
import copy

import pytest

pytest.importorskip("qualtran")

from kaliski_support import lemma_result

from qoptct.foreign.qualtran.kaliski.primes import certify_prime
from qoptct.foreign.qualtran.kaliski.run import run_pair
from qoptct.foreign.qualtran.kaliski.verdict import judge, lemmas_ok


@pytest.fixture(scope="module")
def good_lemmas():
    return lemma_result()


@pytest.fixture(scope="module")
def evidence():
    return run_pair(4, 13, verbose=False)


def verdict(rep, lemmas, **kw):
    return judge(rep, lemmas=lemmas, prime=kw.pop("prime", certify_prime(13)), **kw)["status"]


def test_control_real_files_proved(evidence, good_lemmas):
    assert verdict(evidence, good_lemmas) == "PROVED"


def test_empty_lemma_file(evidence):
    empty = dict(ok=True, implications=0, proved=0, controls=0, controls_ok=0, rows=[], z3="x", source_sha256="0")
    assert verdict(evidence, empty) == "NOT PROVED (lemmas)"
    assert verdict(evidence, None) == "NOT PROVED (lemmas)"
    assert verdict(evidence, {}) == "NOT PROVED (lemmas)"


def test_rows_removed_counters_kept(evidence, good_lemmas):
    assert verdict(evidence, dict(good_lemmas, rows=[])) == "NOT PROVED (lemmas)"


def test_foreign_lemma_file(evidence, good_lemmas):
    foreign = copy.deepcopy(good_lemmas)
    foreign["source_sha256"] = "0" * 64                  # a file made for another lemma source
    assert verdict(evidence, foreign) == "NOT PROVED (lemmas)"
    short = copy.deepcopy(good_lemmas)                    # an older attempt with fewer queries
    short["rows"] = short["rows"][:-30]
    assert verdict(evidence, short) == "NOT PROVED (lemmas)"


def test_substituted_query(evidence, good_lemmas):
    bad = copy.deepcopy(good_lemmas)
    bad["rows"][5]["sha256"] = "0" * 64
    assert verdict(evidence, bad) == "NOT PROVED (lemmas)"
    bad = copy.deepcopy(good_lemmas)
    bad["rows"][7]["name"] = "something else"
    assert not lemmas_ok(bad)
    bad = copy.deepcopy(good_lemmas)
    bad["rows"][0], bad["rows"][1] = bad["rows"][1], bad["rows"][0]        # same queries, other order
    assert not lemmas_ok(bad)


def test_a_failed_lemma_row(evidence, good_lemmas):
    bad = copy.deepcopy(good_lemmas)
    bad["rows"][3]["ok"] = False
    bad["rows"][3]["answer"] = "sat"
    assert verdict(evidence, bad) == "NOT PROVED (lemmas)"
    bad = copy.deepcopy(good_lemmas)
    bad["proved"] -= 1
    assert verdict(evidence, bad) == "NOT PROVED (lemmas)"
    bad = copy.deepcopy(good_lemmas)
    bad["z3"] = None                                       # no Z3 version recorded
    assert verdict(evidence, bad) == "NOT PROVED (lemmas)"


def test_no_wire_flow_check(evidence, good_lemmas):
    rep = copy.deepcopy(evidence)
    del rep["prove"]["flow"]
    assert verdict(rep, good_lemmas) == "NOT PROVED (flow)"
    rep = copy.deepcopy(evidence)
    del rep["run"]["flow"]
    assert verdict(rep, good_lemmas) == "NOT PROVED (flow)"


def test_old_result_without_flow_version(evidence, good_lemmas):
    rep = copy.deepcopy(evidence)
    rep["prove"]["flow"].pop("version")
    assert verdict(rep, good_lemmas) == "NOT PROVED (flow)"


def test_too_few_flow_checks(evidence, good_lemmas):
    rep = copy.deepcopy(evidence)
    rep["prove"]["flow"]["eq_checks"] -= 1
    assert verdict(rep, good_lemmas) == "NOT PROVED (flow)"
    rep = copy.deepcopy(evidence)
    rep["run"]["flow"]["junk_ok"] = False
    assert verdict(rep, good_lemmas) == "NOT PROVED (flow)"


def test_result_made_by_another_client(evidence, good_lemmas):
    rep = copy.deepcopy(evidence)
    rep["manifest_before"]["src"]["each"]["kaliski/chain.py"] = "0" * 64
    assert verdict(rep, good_lemmas) == "NOT PROVED (client_hash)"
    rep = copy.deepcopy(evidence)
    rep["manifest_before"]["src"]["each"].pop("kaliski/lemmas_z3.py")
    assert verdict(rep, good_lemmas) == "NOT PROVED (client_hash)"


def test_failed_chain_proof(evidence, good_lemmas):
    rep = copy.deepcopy(evidence)
    rep["prove"]["status"] = "CHAIN_FAILED"
    assert verdict(rep, good_lemmas).startswith("NOT PROVED (chain")
    rep = copy.deepcopy(evidence)
    rep["prove"]["admitted"] = ["something admitted"]
    assert verdict(rep, good_lemmas).startswith("NOT PROVED (chain")
    rep = copy.deepcopy(evidence)
    rep["prove"]["assumed"] = dict(rep["prove"]["assumed"], **{"M-new": 1})      # an assumption no lemma covers
    assert verdict(rep, good_lemmas) == "NOT PROVED (assumptions_covered)"


def test_run_mode_failures(evidence, good_lemmas):
    for key in ("violations", "n_wrong_output", "n_differs", "n_pq_mismatch"):
        rep = copy.deepcopy(evidence)
        rep["run"][key] = 1
        assert verdict(rep, good_lemmas) == "NOT PROVED (run)", key
    rep = copy.deepcopy(evidence)
    rep["run"] = dict(error="wire flow: something")                # the client refused in mode run
    assert verdict(rep, good_lemmas).startswith("NOT PROVED (run")


def test_dirty_or_wrong_simulation(evidence, good_lemmas):
    for key, val in (("dirty_wires", 3), ("n_mismatch", 1), ("gadget_error", ["x", 1])):
        rep = copy.deepcopy(evidence)
        rep["sim"][key] = val
        assert verdict(rep, good_lemmas) == "NOT PROVED (sim)", key
    rep = copy.deepcopy(evidence)
    rep["sim"]["inputs"] = 2
    assert verdict(rep, good_lemmas) == "NOT PROVED (sim)"


def test_sources_changed_during_run(evidence, good_lemmas):
    rep = copy.deepcopy(evidence)
    rep["sources_unchanged_during_run"] = False
    assert verdict(rep, good_lemmas) == "NOT PROVED (sources_unchanged)"


def test_too_few_sample_inputs(evidence, good_lemmas):
    rep = copy.deepcopy(evidence)
    rep["run"]["inputs"] = 1
    assert verdict(rep, good_lemmas) == "NOT PROVED (samples)"


def test_modulus_without_certificate(evidence, good_lemmas):
    assert verdict(evidence, good_lemmas, prime=None) == "NOT PROVED (prime)"
    assert verdict(evidence, good_lemmas, prime=dict(verdict="PROBABLE PRIME (x)")) == "NOT PROVED (prime)"
    assert verdict(evidence, good_lemmas, prime=certify_prime(21)) == "NOT PROVED (prime)"


def test_canaries_foreign_or_weak(evidence, good_lemmas):
    weak = dict(n=256, canaries=[dict(name="T1", counted=True, prove="NOT_PROVED", witness=dict(x=1))] * 3)
    assert verdict(evidence, good_lemmas, canaries=weak) == "NOT PROVED (canaries)"
    false_proved = dict(n=256, canaries=[dict(name="T%d" % i, counted=True, prove="NOT_PROVED", witness=dict(x=1)) for i in range(30)]
                        + [dict(name="R01", counted=False, prove="CHAIN_PROVED", witness=dict(x=1))])
    assert verdict(evidence, good_lemmas, canaries=false_proved) == "NOT PROVED (canaries)"
    assert verdict(evidence, good_lemmas, canaries_required=True) == "NOT PROVED (canaries)"
