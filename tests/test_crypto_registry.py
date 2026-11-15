"""results/crypto: every entry has its circuit and reference, the counts in its .json match the circuit, and the
circuit is proved again by qoptct.verify.kprove (entries whose recorded proof took long are marked slow)."""
import glob
import json
import os
import re

import pytest

from qoptct.verify.kprove import kprove

ROOT = os.path.join(os.path.dirname(__file__), "..", "results", "crypto")
ENTRIES = sorted(glob.glob(os.path.join(ROOT, "*.json")))
ENTRIES = [e for e in ENTRIES if os.path.basename(e) != "index.json"]

pytestmark = pytest.mark.skipif(not ENTRIES, reason="results/crypto not present")


def _counts(txt):
    T = len(re.findall(r'(?:^|[;{]\s*)(?:t|tdg) q\[', txt, re.M))
    CX = len(re.findall(r'(?:^|[;{]\s*)(?:cx|cz) q\[', txt, re.M))
    Q = int(re.search(r'qubit\[(\d+)\]', txt).group(1))
    return T, CX, Q


def _params():
    out = []
    for e in ENTRIES:
        rec = json.load(open(e))
        marks = [pytest.mark.slow] if rec.get("prove_time_s", 0) > 60 else []
        out.append(pytest.param(rec, marks=marks, id=rec["name"]))
    return out


@pytest.mark.parametrize("rec", _params())
def test_entry(rec):
    qasm = os.path.join(ROOT, rec["name"] + ".qasm")
    qc = os.path.join(ROOT, rec["reference"])
    assert os.path.isfile(qasm) and os.path.isfile(qc)
    T, CX, Q = _counts(open(qasm).read())
    assert (T, CX, Q) == (rec["T"], rec["CX"], rec["qubits"])
    assert rec["AND"] == T // 4 and rec["verdict"] == "PROVED"
    r = kprove(qasm, qc, clean=rec["model"] == "clean", budget_sec=max(120, 3 * rec.get("prove_time_s", 0)),
               sample=False)
    assert r["status"] == "PROVED", r


def test_index_lists_all():
    idx = json.load(open(os.path.join(ROOT, "index.json")))
    assert sorted(i["name"] for i in idx["circuits"]) == sorted(json.load(open(e))["name"] for e in ENTRIES)
