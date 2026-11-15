"""qoptct.estimate.physical: each model reproduces the published numbers of its paper; monotonicity; input errors.

Tolerances (relative):
* ge2019_rsa_reference (re-implementation of GE19's exact layout formulas): 10 % -- the paper rounds to 2 digits
  and Table 3 quotes 20 Mqubits / 5.1 h / 31 % risk;
* ge2019 generic (routing factor + factories sized to the reaction rate, not GE's exact layout): qubits and time
  per run 20 %; expected time 30 % (its optimum is d = 29 at 3 % risk, GE's is d = 27 at 31 % risk), 15 % at d = 27;
* gidney2025: qubits 2 % (the formula is the paper's own sum, only the 1399 vs 1537 split of hot/cold differs),
  time 10 % (we use a fixed 2.5 slack over the CCZ supply, the paper sums per-operation durations);
* litinski2023: 5 % (closed formulas of the paper; only rounding)."""
import json
import math

import pytest

from qoptct.estimate import physical as ph
from qoptct.estimate.cli import main_estimate


def rel(a, b):
    return abs(a / b - 1)


@pytest.fixture(scope="module")
def ge_ref():
    return ph.ge2019_rsa_reference()


def test_ge2019_reference_layout_reproduces_paper(ge_ref):
    P = ph.PAPER["ge2019"]
    assert ge_ref["d1"] == 15 and ge_ref["d2"] == 27          # Table 3: d1 = 15, d2 = 27
    assert rel(ge_ref["physical_qubits"], P["physical_qubits"]) < 0.10
    assert rel(ge_ref["hours"], P["runtime_h"]) < 0.10
    assert rel(ge_ref["expected_hours"], P["expected_runtime_h"]) < 0.10
    assert rel(ge_ref["toffoli"], P["toffoli"]) < 0.10
    assert abs(ge_ref["error"] - P["retry_risk"]) < 0.05


def test_ge2019_generic_reproduces_paper(ge_ref):
    P = ph.PAPER["ge2019"]
    r = ph.ge2019(ge_ref["toffoli"], P["qubits"], ge_ref["reaction_depth"])
    assert r.code_distance in (27, 29)
    assert rel(r.physical_qubits, P["physical_qubits"]) < 0.20
    assert rel(r.runtime_s / 3600, P["runtime_h"]) < 0.20
    # the generic model sizes factories smoothly, so its volume optimum is d = 29 with ~3 % risk instead of GE's
    # d = 27 with 31 % risk; the expected time is therefore shorter: 30 %
    assert rel(r.expected_runtime_s / 3600, P["expected_runtime_h"]) < 0.30
    # at GE's own distances d1 = 15, d2 = 27, with GE's approximate-arithmetic (deviation) error added, the risk
    # and the expected time match closely
    r27 = ph.ge2019(ge_ref["toffoli"], P["qubits"], ge_ref["reaction_depth"], code_distance=27, d1_options=(15,),
                    algorithm_error=ge_ref["deviation_error"])
    assert rel(r27.physical_qubits, P["physical_qubits"]) < 0.20
    assert abs(r27.failure_prob - P["retry_risk"]) < 0.10
    assert rel(r27.expected_runtime_s / 3600, P["expected_runtime_h"]) < 0.15
    # default depth (GE's depth/Toffoli ratio, Eq. 3 / Eq. 2) gives the same order
    r2 = ph.ge2019(P["toffoli"], P["qubits"])
    assert rel(r2.runtime_s / 3600, P["runtime_h"]) < 0.20
    assert r2.notes


def test_gidney2025_reproduces_paper():
    P = ph.PAPER["gidney2025"]
    r = ph.gidney2025(P["toffoli"], P["qubits"], shots=P["shots"])
    assert r.code_distance == 25
    assert rel(r.physical_qubits, P["physical_qubits"]) < 0.02
    assert r.physical_qubits < 1e6
    assert rel(r.runtime_s / 3600, P["shot_hours"]) < 0.10
    assert rel(r.expected_runtime_s / 86400, P["expected_runtime_days"]) < 0.10
    assert r.expected_runtime_s < 7 * 86400
    assert abs((1 - r.failure_prob) - P["success"]) < 0.01


def test_litinski2023_reproduces_paper():
    P = ph.PAPER["litinski2023"]
    r10 = ph.litinski2023(P["toffoli"], P["qubits"], cycle_us=10)
    assert r10.code_distance == P["code_distance"]
    assert rel(r10.physical_qubits, P["physical_qubits"]) < 0.05
    assert rel(r10.runtime_s / 60, P["runtime_min_10us"]) < 0.05
    r1 = ph.litinski2023(P["toffoli"], P["qubits"], cycle_us=1)
    assert rel(r1.runtime_s, P["runtime_s_1us"]) < 0.05
    B = ph.PAPER["litinski2023_baseline"]
    b = ph.litinski2023(B["toffoli"], B["qubits"], architecture="baseline")
    assert b.code_distance == B["code_distance"]
    assert rel(b.physical_qubits, B["physical_qubits"]) < 0.05
    assert rel(b.runtime_s / 3600, B["runtime_h_1us"]) < 0.05
    # reaction limit of Sec. 2.3: depth 2T at 10 us -> 36 min
    rr = ph.litinski2023(P["toffoli"], P["qubits"], 2 * P["toffoli"], cycle_us=10)
    assert rr.reaction_limited and rel(rr.runtime_s / 60, 36.3) < 0.05


@pytest.mark.parametrize("model,q,kw", [("ge2019", 6189, {}), ("gidney2025", 1399, {"shots": 9.2}),
                                         ("litinski2023", 3000, {})])
def test_monotone_in_toffoli_and_qubits(model, q, kw):
    T = [1e7, 1e8, 1e9, 4e9, 1e10]
    rs = [ph.estimate(model, t, q, **kw) for t in T]
    for a, b in zip(rs, rs[1:]):
        assert b.expected_runtime_s > a.expected_runtime_s
        assert b.physical_qubits >= a.physical_qubits
        assert b.code_distance >= a.code_distance
    Q = [q // 2, q, 2 * q, 4 * q]
    rq = [ph.estimate(model, 1e9, x, **kw) for x in Q]
    for a, b in zip(rq, rq[1:]):
        assert b.physical_qubits > a.physical_qubits


@pytest.mark.parametrize("model", sorted(ph.MODELS))
def test_monotone_in_depth_and_p(model):
    base = ph.estimate(model, 1e9, 2000, 1e9)
    deep = ph.estimate(model, 1e9, 2000, 1e10)
    assert deep.expected_runtime_s >= base.expected_runtime_s
    good = ph.estimate(model, 1e9, 2000, p=1e-4)
    bad = ph.estimate(model, 1e9, 2000, p=2e-3)
    mid = ph.estimate(model, 1e9, 2000)
    assert good.code_distance <= mid.code_distance <= bad.code_distance
    assert good.physical_qubits <= mid.physical_qubits <= bad.physical_qubits


def test_reaction_limit_gidney2025():
    r = ph.gidney2025(1e9, 1400)
    slow = ph.gidney2025(1e9, 1400, 1e9, reaction_us=100)   # 100 us per reaction > 62.5 us per Toffoli
    assert not r.reaction_limited and slow.reaction_limited
    assert rel(slow.runtime_s, 1e9 * 100e-6) < 1e-9


@pytest.mark.parametrize("bad", [dict(toffoli=0, qubits=10), dict(toffoli=-1, qubits=10),
                                 dict(toffoli=float("nan"), qubits=10), dict(toffoli=float("inf"), qubits=10),
                                 dict(toffoli=1e9, qubits=0), dict(toffoli=1e9, qubits=10.5),
                                 dict(toffoli=1e9, qubits=10, reaction_depth=-5),
                                 dict(toffoli=1e9, qubits=10, p=0.02), dict(toffoli=1e9, qubits=10, p=0),
                                 dict(toffoli="1e9", qubits=10), dict(toffoli=True, qubits=10)])
@pytest.mark.parametrize("model", sorted(ph.MODELS))
def test_input_errors(model, bad):
    with pytest.raises(ValueError):
        ph.estimate(model, **bad)


def test_unknown_model_and_bad_options():
    with pytest.raises(ValueError):
        ph.estimate("nope", 1e9, 10)
    with pytest.raises(ValueError):
        ph.gidney2025(1e9, 10, shots=0.5)
    with pytest.raises(ValueError):
        ph.litinski2023(1e9, 10, architecture="x")
    with pytest.raises(ValueError):
        ph.ge2019(1e9, 10, max_error=1e-30)


def test_cli(capsys):
    assert main_estimate(["--toffoli", "4.91e9", "--qubits", "1400", "--shots", "9.2", "--model", "gidney2025",
                          "--json"]) == 0
    d = json.loads(capsys.readouterr().out)
    assert d[0]["model"] == "gidney2025" and d[0]["physical_qubits"] < 1e6
    assert main_estimate(["--toffoli", "1e9", "--qubits", "2000"]) == 0
    out = capsys.readouterr().out
    assert all(m in out for m in ph.MODELS)
    assert main_estimate(["--qubits", "10"]) == 3
    assert main_estimate(["--toffoli", "-1", "--qubits", "10"]) == 3
    assert main_estimate(["--toffoli", "1e9", "--qubits", "10.5"]) == 3
    assert main_estimate(["--bogus"]) == 3


def test_cli_check(capsys):
    assert main_estimate(["--check", "--json"]) == 0
    rows = json.loads(capsys.readouterr().out)
    assert len(rows) >= 10
    for r in rows:
        assert rel(r["ours"], r["published"]) < 0.30, r
