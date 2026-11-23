"""Point addition (ECAdd, ECWindowAddR) through the package: positive proofs at n=16 and the verdict logic.

The n=256 reruns are marked ``slow`` (about an hour, ~3 GB; the registry results are produced by the same API)."""
import pytest

pytest.importorskip("qualtran")
pytest.importorskip("qoptct.verify.bdd")

from qoptct.foreign.qualtran.ecadd import _z3
from qoptct.foreign.qualtran.ecadd import api, contracts_status
from qoptct.foreign.qualtran.ecadd.contracts import CONTRACTS, P255, instance_covered

needs_z3 = pytest.mark.skipif(_z3.version() is None, reason="Z3 is not available (verdicts without Z3 are never PROVED)")
P16 = 32749


@pytest.fixture(scope="session", autouse=True)
def _cache(tmp_path_factory):
    import os
    old = os.environ.get("QOPT_CACHE")
    os.environ["QOPT_CACHE"] = str(tmp_path_factory.mktemp("qopt_cache"))
    yield
    if old is None:
        os.environ.pop("QOPT_CACHE", None)
    else:
        os.environ["QOPT_CACHE"] = old


@pytest.fixture(scope="session")
def ecadd16():
    return api.prove_ecadd(16, P16)


@needs_z3
def test_ecadd_fixed_n16_proved(ecadd16):
    r = ecadd16
    assert r["status"] == "PROVED", r
    assert r["qualtran_kind"] == "e"
    assert sorted(r["steps"]) == ["1", "2", "3", "4", "5", "6"]
    assert all(s["status"] == "PROVED" and s["kind"] == "e" for s in r["steps"].values())
    assert r["wiring"]["status"] == "OK" and r["wiring"]["order"] == [1, 2, 3, 4, 5, 6]
    g = r["glue"]
    assert g["status"] == "PROVED" and g["wiring_equals_reference"] is True
    assert g["wiring_source"] == "live ECAdd.decompose_bloq()"
    assert g["answers"] == {"unsat": 49} and g["failed"] == []
    assert "Z3" in g["z3_version"] and r["z3"]


@needs_z3
def test_ecadd_verdict_names_assumed_contracts(ecadd16):
    cs = ecadd16["contracts_status"]
    assert cs["MM"]["forward"]["status"] == "PROVED_FOR_THIS_INSTANCE"
    # n=16, p=32749: the Kaliski contract and its pair are an assumption (goal D and H1 did 65521)
    assert cs["KINV"]["forward"]["status"] == "ASSUMED_FOR_THIS_INSTANCE"
    assert cs["KINV"]["pair"]["status"] == "ASSUMED_FOR_THIS_INSTANCE"


@needs_z3
def test_ecadd_original_is_correct_exactly_outside_groups():
    r = api.prove_ecadd(16, P16, patched=False)
    assert r["status"] == "PROVED_OUTSIDE_G1_G3", r
    assert r["qualtran_kind"] == "original"
    assert r["glue"]["mode"] == "E3"
    assert all(s["kind"] == "original" and s["status"] == "PROVED" for s in r["steps"].values())


@needs_z3
def test_partial_run_is_never_proved():
    r = api.prove_ecadd(16, P16, steps=(1, 2))
    assert r["status"] == "NOT_PROVED" and "partial" in r["detail"]


def test_wrong_qualtran_variant_is_refused():
    r = api._check_variant(dict(status="PROVED", kind="original"), "e")
    assert r["status"] == "NOT_PROVED" and "variant" in r["detail"]
    assert api._check_variant(dict(status="PROVED", kind="e"), "e")["status"] == "PROVED"


@needs_z3
@pytest.mark.parametrize("window", [1, 2, 4])
def test_ecwindowaddr_n16_proved(window):
    r = api.prove_ecwindowaddr(16, P16, window)
    assert r["status"] == "PROVED", r
    assert r["e4"]["kind"] == "e4" and r["e4"]["status"] == "PROVED"
    assert "1:QROAMClean" in r["e4"]["regions"] and "2:ECAdd" in r["e4"]["regions"]
    assert r["depends_on"]["proved_by"] == "ecadd"


def test_contract_tables():
    assert instance_covered("MM", 256, P255) and instance_covered("MM", 16, P16)
    assert instance_covered("KINV", 256, P255) and not instance_covered("KINV", 16, P16)
    assert instance_covered("KINV", 16, 65521, "pair") and not instance_covered("KINV", 16, P16, "pair")
    cs = contracts_status(256, P255)
    for box in ("MM", "KINV"):
        assert cs[box]["forward"]["status"] == "PROVED_FOR_THIS_INSTANCE"
        assert cs[box]["pair"]["status"] == "PROVED_FOR_THIS_INSTANCE"
    # every contract names the result that backs it
    for name, c in CONTRACTS.items():
        assert c["forward_proved_by"] and c["statement"] and c["forward_instances"]


@needs_z3
def test_no_z3_is_not_proved(monkeypatch):
    monkeypatch.setattr(_z3, "version", lambda: None)
    monkeypatch.setattr(api._z3, "version", lambda: None)
    r = api.prove_ecadd(16, P16)
    assert r["status"] == "NOT_PROVED" and "Z3" in r["detail"]


@pytest.mark.slow
@needs_z3
def test_ecadd_n256_rerun():
    r = api.prove_ecadd(256, P255)
    assert r["status"] == "PROVED", r
    assert r["contracts_status"]["MM"]["forward"]["status"] == "PROVED_FOR_THIS_INSTANCE"
    assert r["contracts_status"]["KINV"]["pair"]["status"] == "PROVED_FOR_THIS_INSTANCE"


@pytest.mark.slow
@needs_z3
def test_ecwindowaddr_n256_rerun():
    r = api.prove_ecwindowaddr(256, P255, 1)
    assert r["status"] == "PROVED", r
