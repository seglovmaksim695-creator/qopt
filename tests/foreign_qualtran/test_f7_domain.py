"""Domain of the proofs (primality of p, 2p <= 2^n, p < 2^n - 1) and the status PROVED_CONDITIONAL (audit AUD7, remark 1).

The refusals are decided before any heavy work, so these tests take seconds."""
import json

import pytest

from qoptct.foreign.qualtran import cli
from qoptct.foreign.qualtran.domain import check_domain
from qoptct.foreign.qualtran.ecadd.contracts import CONTRACTS, contracts_status, instance_covered
from qoptct.foreign.qualtran import registry

pytest.importorskip("qualtran")


def run(argv, capsys):
    code = cli.main(argv)
    cap = capsys.readouterr()
    return code, cap.out, cap.err


# ---------------------------------------------------------------- the domain function
def test_domain_composite_and_even():
    for kind in ("ecadd", "ecwindowaddr", "kaliski"):
        d = check_domain(kind, 16, 32767)                      # 32767 = 7 * 31 * 151
        assert not d["ok"] and any("prime" in x for x in d["problems"]), d
    assert not check_domain("montgomery", 16, 32748)["ok"]     # even
    assert check_domain("montgomery", 16, 32767)["ok"]         # the multiplier needs only an odd p (contract C)
    assert not check_domain("pair-kaliski", 16, 32766)["ok"]


def test_domain_range_conditions():
    assert not check_domain("montgomery", 16, 32771)["ok"]     # 2p > 2^16
    assert not check_domain("pair-montgomery", 16, 65521)["ok"]
    assert not check_domain("ecadd", 16, 65521)["ok"]          # prime but 2p > 2^n
    for p, n in ((255, 8), (257, 8), (127, 7)):                # the bridge refuses these (audit, A5)
        assert not check_domain("pair-kaliski", n, p)["ok"], (n, p)
    assert not check_domain("kaliski", 7, 127)["ok"]           # p = 2^n - 1
    assert check_domain("pair-kaliski", 16, 65521)["ok"]       # registry result: 2p > 2^n is fine for the inverse
    assert check_domain("kaliski", 16, 65521)["ok"]


def test_domain_accepts_every_registry_instance():
    for e in registry.entries():
        kind = {"montgomery_mul": "montgomery", "pair_montgomery": "pair-montgomery", "kaliski_inverse": "kaliski",
                "pair_kaliski": "pair-kaliski", "ecadd": "ecadd", "ecwindowaddr": "ecwindowaddr"}[e["name"]]
        d = check_domain(kind, e["n"], int(e["p"]))
        assert d["ok"], (e["id"], d["problems"])


# ---------------------------------------------------------------- the commands refuse with a reason and exit 1
@pytest.mark.parametrize("argv", [
    ["prove", "ecadd", "--n", "16", "--p", "32767"],
    ["prove", "ecwindowaddr", "--n", "16", "--p", "32767"],
    ["prove", "kaliski", "--n", "16", "--p", "32767"],
    ["prove", "kaliski", "--n", "8", "--p", "257"],
    ["prove", "montgomery", "--n", "16", "--w", "1", "--p", "32748"],
    ["prove", "montgomery", "--n", "16", "--w", "1", "--p", "40001"],
    ["pair", "montgomery", "--n", "16", "--w", "1", "--p", "40001"],
    ["pair", "kaliski", "--n", "8", "--p", "257"],
    ["pair", "kaliski", "--n", "8", "--p", "255"],
    ["pair", "kaliski", "--n", "7", "--p", "127"],
])
def test_command_refuses_outside_domain(argv, capsys):
    code, out, err = run(argv, capsys)
    assert code == cli.EXIT_NOT_PROVED == 1, (code, out, err)
    assert out.startswith("NOT_PROVED") and "outside the domain" in out and "Traceback" not in err


def test_command_refusal_json(capsys):
    code, out, _ = run(["prove", "ecadd", "--n", "16", "--p", "32767", "--json"], capsys)
    v = json.loads(out)
    assert code == 1 and v["status"] == "NOT_PROVED" and v["domain"]["ok"] is False


def test_api_refusals_without_running_children():
    from qoptct.foreign.qualtran import ecadd, pairs
    from qoptct.foreign.qualtran.kaliski import prove_kaliski
    for v in (ecadd.prove_ecadd(16, 32767), ecadd.prove_ecwindowaddr(16, 32767), pairs.check_pair_kaliski(8, 257),
              pairs.check_pair_montgomery(16, 1, 40001)):
        assert v["status"] == "NOT_PROVED" and v["domain"]["ok"] is False and "domain" in v["detail"], v
    assert prove_kaliski(8, 255)["status"].startswith("NOT PROVED")


# ---------------------------------------------------------------- the registry holds only fully proved results
def test_registry_has_no_conditional_result():
    for e in registry.entries():
        assert e["status"] in ("PROVED", "PROVED_OUTSIDE_G1_G3"), e["id"]
    txt = registry.format_table(registry.entries())
    assert "PROVED_CONDITIONAL" not in txt


def test_registry_point_addition_results_have_all_contracts_proved():
    for e in registry.entries("ecadd") + registry.entries("ecwindowaddr"):
        n, p = e["n"], int(e["p"])
        cs = contracts_status(n, p, ("MM", "KINV") if e["name"] == "ecadd" else ("MM", "KINV", "ECADD"))
        for b, st in cs.items():
            for part in ("forward", "pair"):
                assert st[part] is None or st[part]["status"] == "PROVED_FOR_THIS_INSTANCE", (e["id"], b, part)


def test_contract_evidence_ids_exist_and_are_proved():
    ids = {e["id"]: e for e in registry.entries()}
    for box, c in CONTRACTS.items():
        for key in ("forward_evidence", "pair_evidence"):
            for (n, p), rid in c[key].items():
                assert rid in ids and ids[rid]["status"] == "PROVED" and ids[rid]["n"] == n and int(ids[rid]["p"]) == p, (box, rid)
                assert instance_covered(box, n, p, "forward" if key == "forward_evidence" else "pair")


# ---------------------------------------------------------------- PROVED_CONDITIONAL (n = 16, a prime outside the registry)
@pytest.fixture(scope="module")
def conditional(tmp_path_factory):
    pytest.importorskip("z3")
    import os
    os.environ["QOPT_CACHE"] = str(tmp_path_factory.mktemp("qopt_cache"))
    from qoptct.foreign.qualtran import ecadd
    from qoptct.foreign.qualtran.ecadd import _z3
    if _z3.version() is None:
        pytest.skip("Z3 is not available")
    return ecadd.prove_ecadd(16, 32719)             # prime, 2p <= 2^16, but no registry result: contracts only ASSUMED


def test_ecadd_outside_registry_is_conditional_not_proved(conditional, capsys):
    v = conditional
    assert v["status"] == "PROVED_CONDITIONAL", v.get("detail")
    assert v["status"] != "PROVED" and set(v["conditional_on"]) == {"MM forward", "MM pair", "KINV forward", "KINV pair"}
    assert all(s["status"] == "PROVED" for s in v["steps"].values()) and v["glue"]["status"] == "PROVED"
    assert "box contracts" in v["detail"]


def test_conditional_exit_code_is_3(conditional, capsys, monkeypatch):
    from qoptct.foreign.qualtran import ecadd
    monkeypatch.setattr(ecadd, "prove_ecadd", lambda *a, **k: conditional)
    code, out, _ = run(["prove", "ecadd", "--n", "16", "--p", "32719"], capsys)
    assert code == cli.EXIT_CONDITIONAL == 3
    assert out.startswith("PROVED_CONDITIONAL") and "CONDITIONAL on box contracts" in out and "--prove-boxes" in out


def test_ecwindowaddr_without_ecadd_result_is_conditional(capsys):
    pytest.importorskip("z3")
    from qoptct.foreign.qualtran.ecadd import _z3
    if _z3.version() is None:
        pytest.skip("Z3 is not available")
    from qoptct.foreign.qualtran import ecadd
    v = ecadd.prove_ecwindowaddr(16, 32719, 1)
    assert v["status"] == "PROVED_CONDITIONAL" and v["conditional_on"] == ["ECADD forward"], v.get("detail")


def test_prove_boxes_marks_only_proved_parts(monkeypatch):
    """prove_boxes upgrades a part only when its own proof says exactly PROVED; anything else stays ASSUMED."""
    from qoptct.foreign.qualtran.ecadd import api
    import qoptct.foreign.qualtran.kaliski as K
    import qoptct.foreign.qualtran.montgomery as M
    import qoptct.foreign.qualtran.pairs as P
    calls = []
    monkeypatch.setattr(M, "prove_montgomery", lambda n, w, p, **k: calls.append("mm") or {"verdict": "PROVED"})
    monkeypatch.setattr(P, "check_pair_montgomery", lambda n, w, p, **k: calls.append("mmp") or {"status": "NOT_PROVED"})
    monkeypatch.setattr(K, "prove_kaliski", lambda n, p, **k: calls.append("k") or {"status": "NOT PROVED (lemmas)"})
    monkeypatch.setattr(P, "check_pair_kaliski", lambda n, p, **k: calls.append("kp") or {"status": "PROVED"})
    cs = api.prove_missing_boxes(16, 32719, contracts_status(16, 32719))
    assert cs["MM"]["forward"]["status"] == "PROVED_FOR_THIS_INSTANCE"
    assert cs["MM"]["pair"]["status"] == "ASSUMED_FOR_THIS_INSTANCE" and cs["MM"]["pair"]["attempt"] == "NOT_PROVED"
    assert cs["KINV"]["forward"]["status"] == "ASSUMED_FOR_THIS_INSTANCE"
    assert cs["KINV"]["pair"]["status"] == "PROVED_FOR_THIS_INSTANCE"
    assert sorted(calls) == ["k", "kp", "mm", "mmp"]
    # a crash is never a proof
    monkeypatch.setattr(M, "prove_montgomery", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")))
    cs = api.prove_missing_boxes(16, 32719, contracts_status(16, 32719))
    assert cs["MM"]["forward"]["status"] == "ASSUMED_FOR_THIS_INSTANCE" and cs["MM"]["forward"]["attempt"].startswith("ERROR")
    # a covered part is not run again
    calls.clear()
    api.prove_missing_boxes(16, 32749, contracts_status(16, 32749))
    assert calls == []
