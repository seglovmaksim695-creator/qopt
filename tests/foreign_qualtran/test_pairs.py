"""Pairs "forward + backward = identity" (Kaliski, Montgomery): positive checks, canaries, and the auditor's mutants."""
import pytest

pytest.importorskip("qualtran")
pytest.importorskip("qoptct.verify.compose.lift")

from qoptct.foreign.qualtran._variant import run_json
from qoptct.foreign.qualtran.pairs import check_pair_kaliski, check_pair_montgomery

PKG = "qoptct.foreign.qualtran.pairs"


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


# ---------------------------------------------------------------- positive
def test_pair_kaliski_n8_with_full_lowering():
    r = check_pair_kaliski(8, 251)
    assert r["status"] == "PROVED", r
    assert r["compare_full"] is True and r["equals_full_lowering"] is True
    assert r["qualtran_kind"] == "original"
    assert r["pi_id_on_join"] and r["pi_xin_is_xout"] and r["pi_injective"]


def test_pair_kaliski_n16():
    r = check_pair_kaliski(16, 65521, compare_full=False)
    assert r["status"] == "PROVED", r
    assert r["macros"] == 64 and r["top_parts"]["F_post"] == r["top_parts"]["B_pre"]


def test_pair_montgomery_n8_with_full_lowering():
    r = check_pair_montgomery(8, 1, 127)
    assert r["status"] == "PROVED", r
    assert r["equals_full_lowering"] is True and r["qualtran_kind"] == "c"
    assert r["macros"] == 16                              # 2 * (n / w) window nodes


def test_pair_montgomery_n16():
    r = check_pair_montgomery(16, 1, 32749)
    assert r["status"] == "PROVED", r
    assert r["qualtran_kind"] == "c" and r["pi_id_on_join"] and r["pi_xin_is_xout"]


def test_pair_verdict_is_refused_on_the_wrong_qualtran(monkeypatch):
    # the Montgomery proof is about the FIXED multiplier: a child that ran on another variant must not count
    from qoptct.foreign.qualtran import pairs
    r = pairs._verdict(dict(mirror_ok=True, qualtran_kind="original"), "c")
    assert r["status"] == "NOT_PROVED" and "variant" in r["detail"]
    assert pairs._verdict(dict(mirror_ok=True, qualtran_kind="c"), "c")["status"] == "PROVED"
    assert pairs._verdict(dict(mirror_ok=False, qualtran_kind="c"), "c")["status"] == "NOT_PROVED"
    assert pairs._verdict(dict(mirror_ok=True, qualtran_kind="c"), "c", ("equals_full_lowering",))["status"] == "NOT_PROVED"


# ---------------------------------------------------------------- canaries (every corruption must NOT pass)
def _check_canaries(r):
    assert r["status"] == "OK" and r["base_mirror_ok"] is True, r.get("detail")
    assert r["false_pass"] == []
    assert r["n_fail_canaries"] >= 20 and r["n_caught"] == r["n_fail_canaries"]
    rows = r["rows"]
    assert all(x["verdict"] == "NOT PASS" for x in rows if x["expect"] == "FAIL")
    assert all(x["verdict"] == "PASS" for x in rows if x["expect"] == "EQUIV")      # equivalent permutations are accepted
    names = " ".join(x["name"] for x in rows)
    assert "macro:" in names and "top " in names and "piece:" in names                # all three levels are hit


def test_canaries_kaliski_n8():
    _check_canaries(run_json(None, PKG, ["canary-kaliski", "--n", 8, "--p", 251]))


def test_canaries_montgomery_n8():
    r = run_json("c", PKG, ["canary-montgomery", "--n", 8, "--w", 1, "--p", 127])
    _check_canaries(r)
    assert "token:" in " ".join(x["name"] for x in r["rows"])                         # LD / CLR abstraction canaries


def test_canaries_montgomery_n16():
    _check_canaries(run_json("c", PKG, ["canary-montgomery", "--n", 16, "--w", 1, "--p", 32749]))


@pytest.mark.slow
def test_canaries_kaliski_n16():
    _check_canaries(run_json(None, PKG, ["canary-kaliski", "--n", 16, "--p", 65521]))


# ---------------------------------------------------------------- auditor's mutants of the mirror and of the LD/CLR abstraction
def _setup():
    from qoptct.foreign.qualtran.pairs import mirror as M
    from qoptct.foreign.qualtran.pairs.mont_compose import extract_token
    from qoptct.foreign.qualtran.pairs.pair_compose import CErr
    return M, extract_token, CErr


def test_audit_mirror_match_mutants():
    M, _, _ = _setup()
    K_X, K_CX, K_AND, K_ERASE = M.K_X, M.K_CX, M.K_AND, M.K_ERASE
    X = lambda a: (K_X, a, 0, 0, -1)
    CX = lambda a, b: (K_CX, a, b, 0, -1)
    AND = lambda a, b, t: (K_AND, a, b, t, -1)
    ER = lambda a, b, t, l=-1: (K_ERASE, a, b, t, l)
    F = [CX(0, 2), AND(0, 1, 3), CX(3, 4), ER(0, 1, 3), X(2)]
    JOIN = [0, 1, 2, 4]

    def ren(prog, mp):
        return [(k, mp.get(a, a), mp.get(b, b) if k != K_X else b, mp.get(t, t) if k >= K_AND else t, mp.get(l, l) if l >= 0 else l)
                for (k, a, b, t, l) in prog]
    Bok = ren(M.mirror_prog(F), {3: 13})
    cases = {
        "M0 control: backward = mirror with a fresh ancilla": (Bok, True),
        "M1 two non-commuting gates swapped": ([Bok[0], Bok[2], Bok[1]] + Bok[3:], False),
        "M2 backward takes joint wire 2 instead of a fresh ancilla": (ren(M.mirror_prog(F), {3: 2}), False),
        "M3 target of a cx moved to another joint wire": (Bok[:-1] + [CX(0, 4)], False),
        "M4 a and b of an AND swapped (equivalent)": ([Bok[0], AND(1, 0, 13)] + Bok[2:], True),
        "M5 commuting X and CX swapped (equivalent)": ([Bok[-1]] + Bok[1:-1] + [Bok[0]], True),
        "M6 erase with an added linear correction": (Bok[:3] + [ER(0, 1, 13, 2)] + Bok[4:], False),
        "M7 backward ancilla = an input wire (joint)": (ren(M.mirror_prog(F), {3: 1}), False),
    }
    for name, (B, want) in cases.items():
        r = M.match(F, B, JOIN)
        assert bool(r["ok"]) == want, (name, r["msg"])


def test_audit_ld_clr_abstraction_mutants():
    M, extract_token, CErr = _setup()
    K_X, K_CX, K_AND, K_ERASE = M.K_X, M.K_CX, M.K_AND, M.K_ERASE
    X = lambda a: (K_X, a, 0, 0, -1)
    CX = lambda a, b: (K_CX, a, b, 0, -1)
    AND = lambda a, b, t: (K_AND, a, b, t, -1)
    ER = lambda a, b, t, l=-1: (K_ERASE, a, b, t, l)

    def accepted(kind, ops):
        try:
            extract_token(kind, ops)
            return True
        except CErr:
            return False
    QC = [CX(0, 5), CX(0, 6), X(0), CX(0, 6), X(0)]
    S = [CX(0, 8)]
    assert accepted("QC", QC)                                         # A0 control
    assert not accepted("QC", QC + [AND(0, 5, 7)])                    # A1 AND inside a QROAMClean block
    assert not accepted("QC", QC + [CX(5, 6)])                        # A2 cx not controlled by sel
    assert accepted("QA", S + [ER(0, 8, 5)] + S)                      # A3c control
    assert not accepted("QA", S + [ER(0, 8, 5, 6)] + S)               # A3 erase with a correction
    assert not accepted("QA", S + [ER(0, 8, 5), CX(5, 6)] + S)        # A4 erased wire read as a control
    assert not accepted("QA", S + [X(9), ER(0, 8, 5)] + S)            # A5 internal wire not back to zero


@pytest.mark.slow
def test_pairs_n256_rerun():
    r = check_pair_kaliski(256, 2 ** 255 - 19)
    assert r["status"] == "PROVED", r
    r = check_pair_montgomery(256, 1, 2 ** 255 - 19)
    assert r["status"] == "PROVED", r
