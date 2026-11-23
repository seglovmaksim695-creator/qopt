"""Regression for the point-addition prover: every corrupted scheme, glue and wiring must come out NOT proved.

Sources: the glue and wiring corruptions of goal E, the step canaries (n=16), the ECWindowAddR canaries, and the
auditor's mutants (box interfaces, glue contracts, XOR properties, the brute force of the glue lemmas)."""
import json

import pytest

pytest.importorskip("qualtran")
pytest.importorskip("qoptct.verify.bdd")

from qoptct.foreign.qualtran._variant import run_json
from qoptct.foreign.qualtran.ecadd import _z3
from qoptct.foreign.qualtran.ecadd.zglue import LEMMAS, WMUT
from qoptct.foreign.qualtran.ecadd.zmutants import MUTANTS

needs_z3 = pytest.mark.skipif(_z3.version() is None, reason="Z3 is not available")
PKG = "qoptct.foreign.qualtran.ecadd"
P16, P255 = 32749, 2 ** 255 - 19


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


def glue(n, p, *args):
    return run_json(None, f"{PKG}.zglue", ["--n", n, "--p", p] + list(args))


# ---------------------------------------------------------------- glue: control and corruptions
@needs_z3
@pytest.mark.parametrize("n,p", [(16, str(P16)), (256, "2**255-19")])
def test_glue_control_is_proved(n, p):
    r = glue(n, p, "--mode", "E2")
    assert r["status"] == "PROVED" and r["answers"] == {"unsat": 49}
    assert r["wiring_equals_reference"] is True


@needs_z3
@pytest.mark.parametrize("n,p", [(16, str(P16)), (256, "2**255-19")])
def test_glue_corruptions_not_proved(n, p):
    bad = []
    for m in MUTANTS:                                    # corrupted step contracts
        r = glue(n, p, "--mode", "E2", "--mutant", m)
        if r["status"] == "PROVED":
            bad.append(("mutant", m))
    for m in WMUT:                                       # corrupted wiring between the steps
        r = glue(n, p, "--mode", "E2", "--wmut", m)
        if r["status"] == "PROVED":
            bad.append(("wiring", m))
    for lem in LEMMAS:                                   # a missing lemma
        r = glue(n, p, "--mode", "E2", "--drop", lem)
        if r["status"] == "PROVED":
            bad.append(("drop", lem))
    for o in ("1,2,6", "1", "2"):                        # original steps 1 / 2 do not glue to the point sum
        r = glue(n, p, "--mode", "E2", "--orig", o)
        if r["status"] == "PROVED":
            bad.append(("orig", o))
    assert not bad, bad
    assert len(MUTANTS) >= 10 and len(WMUT) >= 4


@needs_z3
def test_glue_refuses_wiring_that_differs_from_the_reference(tmp_path):
    from qoptct.foreign.qualtran.ecadd.zglue import WIRING
    w = json.loads(json.dumps(WIRING))
    w["steps"]["5"]["x"] = [3, "x"]
    f = tmp_path / "w.json"
    f.write_text(json.dumps(w))
    r = glue(16, str(P16), "--mode", "E2", "--wiring", str(f))
    assert r["status"] != "PROVED" and r["wiring_equals_reference"] is False


@needs_z3
def test_glue_original_groups_are_not_correct():
    # on the groups G1-G3 the original scheme leaves a dirty flag: "always fails" is proved there, at n=16
    for g in (1, 2, 3):
        r = glue(16, str(P16), "--mode", "E3fail", "--orig", "1,2,6", "--group", g)
        assert r["status"] == "PROVED_ALWAYS_FAILS", (g, r)


# ---------------------------------------------------------------- auditor's mutants
def test_audit_lemmas_bruteforce_no_violation():
    r = run_json(None, f"{PKG}.audit_mutants", ["lemmas"])
    assert r["status"] == "OK"
    for p, row in r["lemmas_own"].items():
        assert row["cases"] > 1000 and set(row["violations"].values()) == {0}, (p, row)


@needs_z3
def test_audit_glue_mutants_not_proved():
    r = run_json(None, f"{PKG}.audit_mutants", ["glue"])
    assert r["status"] == "OK"
    for n, rows in r["glue_canaries"].items():
        for name, row in rows.items():
            if name.startswith("G0"):
                assert row["status"] == "PROVED", (n, name)
            elif name.startswith("G8.3"):
                continue                                  # a property that is not needed (reference only)
            else:
                assert row["status"] != "PROVED", (n, name, row)


def test_audit_box_interface_mutants_not_proved():
    r = run_json("e", f"{PKG}.audit_mutants", ["box", "--n", "16"])
    assert r["status"] == "OK"
    rows = r["box_canaries_n16"]
    assert len(rows) == 6
    for name, row in rows.items():
        if name.startswith("B0"):
            assert row["status"] == "PROVED"
        else:
            assert row["status"] != "PROVED", (name, row)


def test_audit_identity_blocks_in_ecadd_are_the_proved_steps():
    r = run_json("e", f"{PKG}.audit_mutants", ["identity"])
    assert r["status"] == "OK"
    for n, row in r["identity"].items():
        assert [s["step"] for s in row["steps"]] == [1, 2, 3, 4, 5, 6]
        assert all(s["same_as_E1"] for s in row["steps"]), n
        boxes = [b for s in row["steps"] for b in s["boxes"]]
        assert boxes and all(b["mod"] and b["bitsize"] == int(n[1:]) for b in boxes)
        assert all(b["window"] == 1 for b in boxes if b["type"] == "DirtyOutOfPlaceMontgomeryModMul")


# ---------------------------------------------------------------- step and ECWindowAddR canaries (n=16)
@pytest.mark.parametrize("step", [1, 2, 3, 4, 5, 6])
def test_step_canaries_n16(step):
    r = run_json("e", f"{PKG}.canary", ["--step", step, "--n", 16, "--p", P16, "--count", 6])
    s = r["summary"]
    assert s["holes"] == 0 and s["errors"] == 0, s
    assert s["refuted_with_witness"] >= 6 and s["targeted_refuted"] >= 3, s
    assert not any(row["status"] == "HOLE" for row in r["rows"])
    # a corruption that is still PROVED must have been shown neutral on the whole precondition by simulation
    assert all(row["neutral_check"]["violations"] == 0 for row in r["rows"] if row["status"] == "PROVED")


@needs_z3
def test_ecwindowaddr_canaries_n16():
    r = run_json("e4", f"{PKG}.e4", ["canary", "--n", 16, "--p", P16, "--w", 2, "--count", 10])
    s = r["summary"]
    assert r["status"] == "OK" and s["holes"] == 0 and s["errors"] == 0, s
    assert s["refuted_with_witness"] >= 8 and s["targeted"] >= 5, s


@pytest.mark.slow
@pytest.mark.parametrize("step", [1, 2, 3, 4, 5, 6])
def test_step_canaries_n16_full_sample(step):
    r = run_json("e", f"{PKG}.canary", ["--step", step, "--n", 16, "--p", P16, "--count", 24])
    s = r["summary"]
    assert s["holes"] == 0 and s["errors"] == 0 and s["refuted_with_witness"] >= 20, s
