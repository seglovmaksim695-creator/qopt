"""Canary mutants for short circuits (a single wide AND, a row of X gates) and the evenly spread sample for huge ones."""
import warnings

import pytest

pytest.importorskip("qiskit")
warnings.filterwarnings("ignore")

from qoptct.foreign.qiskit.library import FAMILIES
from qoptct.foreign.qiskit.mutants import make_mutants


@pytest.mark.parametrize("fam,n", [("and", 8), ("or", 6), ("xor-const", 8)])
def test_short_circuits_get_at_least_twenty_mutants(fam, n):
    c = FAMILIES[fam](n)
    q = next(q for q in c.circuit.qregs if len(q) >= 3)
    ms = make_mutants(c.circuit, rare_wires=[c.circuit.find_bit(b).index for b in q])
    assert len(ms) >= 20, len(ms)


def test_sample_is_lazy_and_spread():
    c = FAMILIES["cdkm-half"](6)
    full = make_mutants(c.circuit)
    got = list(make_mutants(c.circuit, sample=10))
    assert len(got) == 10
    labels = [l for l, _ in full]
    pos = [labels.index(l) for l, _ in got]
    assert pos == sorted(pos) and pos[0] == 0 and pos[-1] == len(labels) - 1
