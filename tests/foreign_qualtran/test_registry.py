"""The registry of proved results: consistent with the result files, with the contract table and with the shipped evidence.

These tests need neither Qualtran nor Z3 (they read data files only)."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from qoptct.foreign.qualtran import registry
from qoptct.foreign.qualtran.ecadd.contracts import CONTRACTS

REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "results" / "qualtran"
needs_results = pytest.mark.skipif(not RESULTS.is_dir(), reason="results/qualtran is not part of this tree")


def test_registry_shape():
    reg = registry.load_registry()
    assert reg["schema"] == "qoptct.qualtran.registry/1" and reg["qualtran_version"] == "0.7.0"
    es = reg["entries"]
    assert es and len({e["id"] for e in es}) == len(es)
    assert {e["name"] for e in es} == set(registry.NAMES)
    for e in es:
        assert e["status"] in ("PROVED", "PROVED_OUTSIDE_G1_G3"), e["id"]
        if e["status"] != "PROVED":
            assert e["name"] == "ecadd" and "original" in e["id"]            # the only non-PROVED result is the original scheme
        assert e["assumptions"] and e["claim"] and e["block"]
        assert e["reproduce"].startswith("qoptct qualtran ")
        assert len(e["program"]["value"]) == 64 and int(e["program"]["value"], 16) >= 0
        assert e["qualtran"]["version"] == "0.7.0"
        if e["qualtran"]["variant"] != "original":
            assert e["qualtran"]["patch_files_sha256"]
        assert e["id"].split("/")[0] == e["name"] and e["source"]["file"].endswith(".json")


def test_table_and_long_forms():
    es = registry.entries()
    txt = registry.format_table(es)
    assert all(name in txt for name in registry.NAMES)
    assert registry.entries("ecadd", 16) and not registry.entries("ecadd", 17)
    long = registry.format_entry(registry.entries("kaliski_inverse", 16)[0])
    assert "assumptions:" in long and "reproduce:" in long


@needs_results
def test_registry_matches_result_files():
    assert registry.verify_sources(RESULTS) == []


@needs_results
def test_registry_is_up_to_date():
    spec = importlib.util.spec_from_file_location("build_qualtran_registry", REPO / "scripts" / "build_qualtran_registry.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.main(["--check"]) == 0


def test_contract_evidence_exists_and_is_proved():
    by_id = {e["id"]: e for e in registry.entries()}
    seen = 0
    for box, c in CONTRACTS.items():
        for which, name in (("forward_evidence", c["forward_proved_by"]), ("pair_evidence", c["pair_proved_by"])):
            for (n, p), rid in c[which].items():
                e = by_id[rid]
                assert e["name"] == name and e["status"] == "PROVED", (box, n, rid)
                assert (e["n"], int(e["p"])) == (n, p), (box, n, rid)
                seen += 1
    assert seen >= 18
    # S3: the Kaliski contract and its pair at n = 16, p = 32749 are results, not assumptions
    assert "kaliski_inverse/kaliski_n16_p32749" in by_id and "pair_kaliski/pairs_kaliski_n16_p32749" in by_id


def test_every_result_is_either_in_the_registry_or_supporting():
    if not RESULTS.is_dir():
        pytest.skip("results/qualtran is not part of this tree")
    reg = registry.load_registry()
    known = {e["source"]["file"] for e in reg["entries"]} | {s["file"] for s in reg["support_files"]}
    assert {p.name for p in RESULTS.glob("*.json")} <= known


def test_shipped_canaries_are_the_evidence_file_and_come_from_this_client():
    path = registry.packaged_canaries()
    assert path, "the canary evidence of Kaliski n = 256 must be shipped with the package"
    data = json.loads(Path(path).read_text())
    if RESULTS.is_dir():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == \
            hashlib.sha256((RESULTS / "kaliski_canaries_n256_p25519.json").read_bytes()).hexdigest()
    from qoptct.foreign.qualtran.kaliski import verdict
    ok, summary = verdict.canaries_summary(data)
    assert ok and summary["made_by_this_client"] and summary["false_proved"] == [], summary
    assert data["n"] == 256 and summary["counted"] >= 20 and summary["targeted"] >= 10
