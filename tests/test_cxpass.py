# SPDX-License-Identifier: Apache-2.0
"""pytest tests/test_cxpass.py -- CX post-passes for measurement-based circuits (qoptct.contrib.cxpass)."""
import importlib.util
import os
import pathlib
import random
import re
import sys

import pytest

import qoptct
from qoptct.contrib import cxpass
from qoptct.contrib.cxpass import check, ir

BENCH = pathlib.Path(os.environ["QOPT_BENCH"])
PKG = pathlib.Path(qoptct.__file__).resolve().parent


def bench(name):
    p = BENCH / "std" / f"{name}.qc"
    if not p.exists():
        pytest.skip(f"benchmark {p} not found (set QOPT_BENCH)")
    return str(p)


_CACHE = {}


def c2_circuit(name, clean=False):
    """Lowest-T output of the logical-AND compiler (contrib c2) for a std benchmark (a typical --mode t circuit)."""
    key = (name, clean)
    if key not in _CACHE:
        p = bench(name)
        sys.path.insert(0, str(PKG / "contrib" / "c2"))
        spec = importlib.util.spec_from_file_location("c2_run_t", PKG / "contrib" / "c2" / "run.py")
        c2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c2)
        cands, n = c2.candidates(p, 3, None, clean=clean)
        k = min(cands, key=lambda k: (k[0], k[1], k[2]))
        _CACHE[key] = (p, cands[k])
    return _CACHE[key]


def test_ir_normalises_bits_and_conditions():
    txt = """OPENQASM 3.0;
include "stdgates.inc";
qubit[3] q;
bit[1] c;
h q[2];
c[0] = measure q[2];
if (c[0]) {cz q[0], q[1]; x q[2];}
h q[2];
c[0] = measure q[2];
if (c[0] ^ c[0]) {z q[0];}
if ((c[0]) && (c[0])) {cz q[0], q[1]; x q[2];}
reset q[2];
"""
    N, items = ir.parse(txt)
    assert N == 3
    assert [it[2] for it in items if it[0] == "M"] == [0, 1]           # one fresh bit per measurement
    ifs = [it for it in items if it[0] == "IF"]
    assert ifs[0][1] == "(c[0])" and ifs[1][1] == "(c[1]) && (c[1])"    # c[0]^c[0] == 0 dropped
    assert ("R", 2) in items
    out = ir.emit(N, items)
    assert ir.parse(out)[1] == items and "bit[2] c;" in out
    for bad in ("if (!c[0]) {x q[0];}", "ccx q[0], q[1], q[2];", "if (c[0] == 1) {x q[0];}"):
        with pytest.raises(ir.ParseError):
            ir.parse("qubit[3] q;\nbit[1] c;\nc[0] = measure q[0];\n" + bad)


def test_improve_c2_strict_with_reference():
    p, txt = c2_circuit("mod_mult_55")
    T0, CX0 = ir.counts_text(txt)
    out, rep = cxpass.improve_report(txt, False, 60, ref=p)
    T1, CX1 = ir.counts_text(out)
    assert T1 <= T0 and CX1 <= CX0
    assert CX1 < CX0, rep          # the conditional-phase passes save at least one CZ here
    assert rep["verified"]
    f = check.check_vs_ref(out, p, clean=False, trials=20, seed=123)
    assert f > 1 - 1e-9
    try:  # independent dense simulator of the verify team, if present
        from qoptct.verify import check_dense
        assert check_dense(out, p, clean=False, trials=10)["ok"]
    except ImportError:
        pass


def test_improve_without_reference_checks_against_input():
    p, txt = c2_circuit("mod_mult_55")
    out = cxpass.improve(txt, False, 60)
    T0, CX0 = ir.counts_text(txt); T1, CX1 = ir.counts_text(out)
    assert T1 <= T0 and CX1 <= CX0
    assert check.check_vs_ref(out, p, trials=12) > 1 - 1e-9   # still equal to the specification


def test_improve_clean_model():
    p, txt = c2_circuit("tof_4", clean=True)
    out = cxpass.improve(txt, True, 30, ref=p)
    T0, CX0 = ir.counts_text(txt); T1, CX1 = ir.counts_text(out)
    assert T1 <= T0 and CX1 <= CX0
    assert check.check_vs_ref(out, p, clean=True, trials=12) > 1 - 1e-9


def _mutants(txt, rng):
    L = txt.splitlines()
    out = []
    ifcz = [i for i, l in enumerate(L) if l.startswith("if") and "cz" in l]
    if ifcz:
        i = rng.choice(ifcz); M = list(L); del M[i]; out.append(("drop conditional cz", M))
    xr = [i for i, l in enumerate(L) if l.startswith("if") and "^" in l]
    if xr:
        i = rng.choice(xr); M = list(L); M[i] = re.sub(r" \^ c\[\d+\]", "", M[i], count=1); out.append(("drop XOR bit", M))
    cx = [i for i, l in enumerate(L) if l.startswith("cx ")]
    i = rng.choice(cx); M = list(L); del M[i]; out.append(("drop cx", M))
    ts = [i for i, l in enumerate(L) if l.startswith("t ")]
    if ts:
        i = rng.choice(ts); M = list(L); M[i] = "tdg" + M[i][1:]; out.append(("t -> tdg", M))
    return [(k, "\n".join(M) + "\n") for k, M in out]


def test_canaries_rejected():
    p, txt = c2_circuit("mod_mult_55")
    out = cxpass.improve(txt, False, 60, ref=p)
    assert check.verify(out, False, ref_qc=p, official=False)
    assert check.verify(out, False, ref_qc=p, official=False, trials=30, strong=True)
    rng = random.Random(5)
    muts = _mutants(out, rng) + _mutants(txt, rng)
    assert len(muts) >= 6
    for kind, m in muts:
        assert not check.verify(m, False, ref_qc=p, official=False), kind
        assert not check.verify(m, False, ref_qc=p, official=False, trials=30, strong=True), kind
        assert not check.verify(m, False, ref_qc=p, official=True), kind          # qoptct.verify.prove (REFUTED)
        # reference-free mode: a mutant is never accepted as equal to the input program
        N, items = ir.parse(txt)
        D = check.deterministic_wires(N, items, range(N))
        assert not check.verify(m, False, ref_text=ir.normalize(txt), D=D), kind


def test_every_lost_conditional_cz_rejected():
    """Dropping ANY single conditional CZ (incl. XOR-conditioned ones) of the improved circuit is caught by the final check."""
    p, txt = c2_circuit("mod_mult_55")
    out = cxpass.improve(txt, False, 60, ref=p)
    L = out.splitlines()
    idx = [i for i, l in enumerate(L) if l.startswith("if") and "cz" in l]
    assert idx
    for i in idx:
        m = "\n".join(L[:i] + L[i + 1:]) + "\n"
        assert not check.verify(m, False, ref_qc=p, official=False, trials=30, strong=True), L[i]


def test_wrong_pass_never_returned(monkeypatch):
    """A pass that returns a cheaper but wrong circuit must be rejected (the input comes back unchanged)."""
    p, txt = c2_circuit("mod_mult_55")

    def bad_pass(name, text, n_data, deadline, seed=1):
        L = text.splitlines()
        i = max(k for k, l in enumerate(L) if l.startswith("cx "))
        return "\n".join(L[:i] + L[i + 1:]) + "\n"
    monkeypatch.setattr(cxpass, "_run_pass", bad_pass)
    assert cxpass.improve(txt, False, 20, ref=p) == txt
    assert cxpass.improve(txt, False, 20) == txt


def test_bad_input_returned_unchanged():
    for t in ("", "garbage", "OPENQASM 3.0;\nqubit[2] q;\nccx q[0], q[1], q[1];\n",
              "OPENQASM 3.0;\nqubit[2] q;\nbit[1] c;\nc[0] = measure q[0];\nif (!c[0]) {x q[1];}\n"):
        assert cxpass.improve(t, False, 5) == t
    u = "OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[2] q;\nh q[0];\ncx q[0], q[1];\n"
    assert cxpass.improve(u, False, 5) == u          # nothing to gain: the very same text comes back
