"""prove_kaliski: positive results at small n, the n=16 pair against the reference of the audited round, the n=256 pairs (slow)."""
import copy

import pytest

pytest.importorskip("qualtran")

from kaliski_support import P25519, P256, SMALL_PAIRS, canary_result, lemma_result, reference

from qoptct.foreign.qualtran.kaliski import prove_kaliski
from qoptct.foreign.qualtran.kaliski.primes import certify_prime
from qoptct.foreign.qualtran.kaliski.verdict import flow_expected


@pytest.fixture(scope="module")
def lemmas():
    return lemma_result()


@pytest.mark.parametrize("n,p", SMALL_PAIRS)
def test_small_pairs_proved(n, p, lemmas):
    v = prove_kaliski(n, p, lemmas=lemmas)
    assert v["status"] == "PROVED", v["checks"]
    assert all(v["checks"].values())
    assert v["lemmas"]["implications"] == 243 and v["lemmas"]["proved"] == 243
    assert v["lemmas"]["z3"]
    assert v["program"]["top_ops_sha256"] and v["prime"]["verdict"].startswith("PROVEN PRIME")


def test_verdict_has_the_documented_fields(lemmas, tmp_path):
    out = tmp_path / "v.json"
    v = prove_kaliski(4, 13, lemmas=lemmas, out=str(out))
    for k in ("status", "checks", "program", "prime", "lemmas", "canaries", "client", "qualtran", "manifest_sha256", "prove"):
        assert k in v
    assert v["qualtran"] == "0.7.0"
    assert v["client"]["flow_version"] == "flow-1" and len(v["client"]["chain_sha256"]) == 64
    assert out.exists()


def test_composite_modulus_is_refused(lemmas):
    from qoptct.foreign.qualtran.kaliski.run import run_pair
    from qoptct.foreign.qualtran.kaliski.verdict import judge
    rep = run_pair(4, 13, verbose=False)            # real evidence, but the claimed modulus is declared composite
    row = certify_prime(21)
    assert row["verdict"] == "NOT A PRIME"
    assert judge(rep, lemmas=lemmas, prime=row)["status"] == "NOT PROVED (prime)"


def test_probable_prime_alone_is_not_accepted(lemmas):
    v = prove_kaliski(4, 13, lemmas=lemmas, prime=dict(verdict="PROBABLE PRIME (probabilistic test only)"))
    assert v["status"] == "NOT PROVED (prime)"


def test_no_lemmas_never_proved():
    v = prove_kaliski(4, 13, lemmas=dict(ok=False, error="no Z3"))
    assert v["status"] == "NOT PROVED (lemmas)"


def test_n256_needs_canaries(lemmas):
    # the evidence is not run at n = 256 here: only the rule is checked, on the small evidence with the rule switched on
    from qoptct.foreign.qualtran.kaliski.run import run_pair
    from qoptct.foreign.qualtran.kaliski.verdict import judge
    rep = run_pair(4, 13, verbose=False)
    assert judge(rep, lemmas=lemmas, prime=certify_prime(13), canaries_required=True)["status"] == "NOT PROVED (canaries)"
    assert judge(rep, lemmas=lemmas, prime=certify_prime(13), canaries=dict(canaries=[]))["status"] == "NOT PROVED (canaries)"


def test_n16_matches_the_audited_reference(lemmas):
    ref = [r for r in reference() if r["n"] == 16][0]
    v = prove_kaliski(16, int(ref["p"]), lemmas=lemmas)
    assert v["status"] == ref["status"] == "PROVED"
    assert v["program"]["ops"] == ref["ops"] and v["program"]["wires"] == ref["wires"]
    assert v["program"]["top_ops_sha256"] == ref["top_ops_sha256"]
    assert v["program"]["piece_qasm_sha256"] == ref["piece_qasm_sha256"]
    assert v["prove"]["segments"] == ref["segments"] and v["prove"]["maxnodes"] == ref["maxnodes"]


def test_flow_formula_is_written_in_the_verdict_not_taken_from_the_client():
    from qoptct.foreign.qualtran.kaliski.chain import expected_flow_counts
    for n in (4, 16, 256):
        assert flow_expected(n) == expected_flow_counts(n)


@pytest.mark.slow
@pytest.mark.parametrize("name", ["2^255-19", "secp256k1", "P-256"])
def test_n256_pairs_proved_and_equal_to_the_reference(name, lemmas):
    can = canary_result()
    if can is None:
        pytest.skip("no canary evidence (results/qualtran/kaliski_canaries_n256_p25519.json)")
    ref = [r for r in reference() if r["name"] == name][0]
    v = prove_kaliski(256, int(ref["p"]), lemmas=lemmas, canaries=can)
    assert v["status"] == "PROVED", v["checks"]
    assert v["program"]["ops"] == ref["ops"] and v["program"]["wires"] == ref["wires"]
    assert v["program"]["top_ops_sha256"] == ref["top_ops_sha256"]
    assert v["prove"]["segments"] == ref["segments"] and v["prove"]["maxnodes"] == ref["maxnodes"]
