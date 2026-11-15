"""Tests of qoptct.contrib.gf2m (bilinear GF(2^m) multipliers and known constructions)."""
import os
import re
import tempfile

import pytest

from qoptct.contrib import gf2m
from qoptct.contrib.gf2m import anf, check, compile as comp, field, flip

BENCH = os.environ["QOPT_BENCH"]


def bench(rel):
    p = os.path.join(BENCH, rel)
    if not os.path.exists(p):
        pytest.skip("benchmark not available: " + rel)
    return p


def write_tmp(text, suffix):
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    return path


def numeric(qasm, qc, clean):
    return check.numeric_check(qasm, qc, clean, trials=3)


def gf_qc(m, f, order=None, inputs=True):
    """schoolbook GF(2^m) multiplier c ^= a*b mod f as a .qc text (optionally with permuted wire order)."""
    names = [f"a{i}" for i in range(m)] + [f"b{i}" for i in range(m)] + [f"c{i}" for i in range(m)]
    vs = list(names) if order is None else [names[k] for k in order]
    lines = [".v " + " ".join(vs)]
    if inputs:
        lines.append(".i " + " ".join(f"a{i}" for i in range(m)) + " " + " ".join(f"b{i}" for i in range(m)))
    lines.append("BEGIN")
    for i in range(m):
        for j in range(m):
            v = field.pmod(1 << (i + j), f)
            for k in range(m):
                if v >> k & 1:
                    lines.append(f"tof a{i} b{j} c{k}")
    lines.append("END")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------ gadget lemmas (numerically)

def test_and_and_erase_gadgets_exact():
    qc = write_tmp(".v a b c\n.i a b\nBEGIN\ntof a b c\nEND\n", ".qc")
    E = comp.Emit(3)
    E.AND(0, 1, 2)
    assert numeric(E.text(), qc, clean=True) > 1 - 1e-12
    E = comp.Emit(4)  # strict: c ^= ab through an extra ancilla, erased by measurement
    E.AND(0, 1, 3)
    E.cx(3, 2)
    E.erase(0, 1, 3)
    assert numeric(E.text(), qc, clean=False) > 1 - 1e-12
    os.unlink(qc)


# ------------------------------------------------------------------------------ recognition + synthesis

@pytest.mark.parametrize("clean,cx_lab", [(False, 76), (True, 70)])
def test_gf2_4(clean, cx_lab):
    qc = bench("std/gf2^4_mult.qc")
    res = gf2m.candidates(qc, clean, budget_sec=4)
    assert res
    T, CX = check.counts(res[0]["qasm"])
    assert T == 36 and CX <= cx_lab
    assert numeric(res[0]["qasm"], qc, clean) > 1 - 1e-9


def test_gf2_5_clean():
    qc = bench("std/gf2^5_mult.qc")
    res = gf2m.candidates(qc, True, budget_sec=4)
    assert res and check.counts(res[0]["qasm"])[0] == 52


def test_gf2_10_clean_with_input_target():
    # the benchmark lists c9 in .i: in the clean model c9 is an input and must be updated in place
    qc = bench("large/gf2^10_mult.qc")
    info = gf2m.analyse(qc, True)
    assert info["st"] is not None and not info["st"]["zero"]
    res = gf2m.candidates(qc, True, budget_sec=4)
    assert res and check.counts(res[0]["qasm"])[0] == 132


def test_ham15_known_construction():
    qc = bench("std/ham15-low.qc")
    res = gf2m.candidates(qc, True, budget_sec=5)
    assert res and check.counts(res[0]["qasm"]) == (60, 179)
    assert gf2m.candidates(qc, False, budget_sec=2) == []  # strict ham15 is not handled


def test_permuted_wires_and_other_modulus():
    # GF(2^4) with f = x^4 + x^3 + 1 (not the benchmark modulus), wires shuffled
    order = [5, 9, 0, 11, 2, 7, 4, 1, 10, 3, 8, 6]
    qc = write_tmp(gf_qc(4, 0b11001, order), ".qc")
    info = gf2m.analyse(qc, True)
    A, B, T, terms = info["tensor"]
    m, f = field.identify(A, B, T, terms)[:2]
    assert (m, f) == (4, 0b11001)
    res = gf2m.candidates(qc, True, budget_sec=3)
    assert res and check.counts(res[0]["qasm"])[0] == 36
    assert numeric(res[0]["qasm"], qc, True) > 1 - 1e-9
    os.unlink(qc)


def test_generic_bilinear_form():
    # not a field multiplication (2 x 3 -> 3): strict, with target-linear and source-linear parts
    txt = ".v a0 a1 b0 b1 b2 c0 c1 c2\nBEGIN\ntof a0 b1 c0\ntof a1 b0 c0\ntof a1 b2 c0\ntof a0 b0 c1\n" \
          "tof a1 b1 c1\ntof a0 b2 c2\ntof a1 b2 c2\ntof a0 b1 c2\ntof c0 c1\ntof a0 c1\nX c2\nEND\n"
    qc = write_tmp(txt, ".qc")
    info = gf2m.analyse(qc, False)
    A, B, T, terms = info["tensor"]
    assert field.identify(A, B, T, terms) is None
    res = gf2m.candidates(qc, False, budget_sec=3)
    assert res and "flip" in res[0]["note"]
    assert numeric(res[0]["qasm"], qc, False) > 1 - 1e-9
    os.unlink(qc)


def test_unrecognised_returns_empty():
    assert gf2m.candidates(bench("std/tof_4.qc"), True, budget_sec=2) == []
    assert gf2m.candidates(bench("std/qft_4.qc"), False, budget_sec=2) == []
    assert gf2m.candidates("/nonexistent/file.qc", True, budget_sec=1) == []


def test_flip_graph_keeps_tensor():
    t = field.field_tensor(4, 0b10011)
    naive = [(1 << i, 1 << j, c) for (i, j), c in t.items()]
    import time
    d = flip.search(naive, time.time() + 1.5, seed=3)
    assert anf.tensor_of(d) == t and len(d) < len(naive)


def test_change_basis_and_isotope():
    tab = gf2m._data("decomps.json")
    ent = next(e for e in tab["fields"] if e["m"] == 5)
    d0 = [tuple(x) for x in ent["decs"][0]["dec"]]
    assert field.check_dec(d0, 5, ent["decs"][0]["f"])
    for f in (0b100101, 0b111101):
        d1 = field.change_basis(d0, 5, ent["decs"][0]["f"], f)
        assert field.check_dec(d1, 5, f) and len(d1) == 13
        assert field.check_dec(field.isotope(d1, 5, f, 2, 7, 19), 5, f)
    assert field.check_dec(field.karatsuba_field(7, 0b10000011), 7, 0b10000011)
    # large field: root finding by trace splitting (no exhaustive search)
    e16 = next(e for e in tab["fields"] if e["m"] == 16)["decs"][0]
    f2 = next(f for f in range((1 << 16) | 1, 1 << 17, 2) if field.is_irreducible(f) and f != e16["f"])
    d16 = field.change_basis([tuple(x) for x in e16["dec"]], 16, e16["f"], f2)
    assert field.check_dec(d16, 16, f2) and len(d16) == 60
    assert [f for f in range(16, 64) if field.is_irreducible(f)] == [19, 25, 31, 37, 41, 47, 55, 59, 61]


def test_decomposition_table_ranks():
    tab = gf2m._data("decomps.json")
    ranks = {e["m"]: e["rank"] for e in tab["fields"]}
    assert ranks[4] == 9 and ranks[8] == 24 and ranks[10] == 33 and ranks[16] == 60 and ranks[32] == 135
    for e in tab["fields"]:
        for d in e["decs"][:1]:
            assert field.check_dec([tuple(x) for x in d["dec"]], e["m"], d["f"])


# ------------------------------------------------------------------------------ canaries

def _mutants(q):
    out = []
    i = q.find("if (c[")
    if i >= 0:  # drop the CZ of one erasure
        out.append(re.sub(r"if \(c\[(\d+)\]\) \{cz q\[\d+\], q\[\d+\]; x q\[(\d+)\];\}",
                          r"if (c[\1]) {x q[\2];}", q, count=1))
    lines = q.split("\n")
    cx = [k for k, L in enumerate(lines) if L.startswith("cx ")]
    out.append("\n".join(lines[:cx[len(cx) // 2]] + lines[cx[len(cx) // 2] + 1:]))  # drop a CX
    out.append(q.replace("tdg q[", "t q[", 1))  # T -> T
    k = next(k for k, L in enumerate(lines) if L.startswith("h q["))
    out.append("\n".join(lines[:k] + lines[k + 10:]))  # drop a whole AND
    return out


def test_canaries_rejected():
    qc = bench("std/gf2^4_mult.qc")
    info = gf2m.analyse(qc, True)
    res = gf2m.candidates(qc, True, budget_sec=2)
    q = res[0]["qasm"]
    assert check.structural_check(q, info["P"], info["n"], info["inp"], True)[0]
    for bad in _mutants(q):
        assert not check.structural_check(bad, info["P"], info["n"], info["inp"], True)[0]
        # a missing conditional CZ is only visible when that measurement gives 1: enough trials (miss prob 2^-16)
        assert check.numeric_check(bad, qc, True, trials=16) < 0.99
    # a correct circuit for the clean model is rejected in the strict model (c must be updated in place)
    infos = gf2m.analyse(qc, False)
    assert not check.structural_check(q, infos["P"], infos["n"], infos["inp"], False)[0]


def test_ham15_canary():
    qc = bench("std/ham15-low.qc")
    info = gf2m.analyse(qc, True)
    q = gf2m.candidates(qc, True, budget_sec=3)[0]["qasm"]
    for bad in _mutants(q):
        assert not check.structural_check(bad, info["P"], info["n"], info["inp"], True)[0]
