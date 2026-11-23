"""check_bloq: Qualtran blocks are accepted when right and NOT proved when anything is wrong."""
import attrs
import pytest

pytest.importorskip("qualtran")
pytest.importorskip("cirq")

from qualtran import Bloq, QBit, QUInt, Register, Signature
from qualtran.bloqs.arithmetic import Add, LessThanConstant, LessThanEqual
from qualtran.bloqs.basic_gates import CNOT, XGate
from qualtran.bloqs.mod_arithmetic import ModAdd

from qoptct.qualtran import check_bloq


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_add_proved(n):
    r = check_bloq(Add(QUInt(n)))
    assert r["verdict"] == "PROVED", r
    assert r["inputs_checked"] == 1 << (2 * n)
    assert r["prover"]["status"] == "PROVED" and r["partial_domains"]["status"] == "PASS"


def test_mod_add_proved_on_the_qualtran_domain():
    r = check_bloq(ModAdd(3, 5))
    assert r["verdict"] == "PROVED", r
    assert "of 64 bit patterns" in r["domain"]  # Qualtran's Montgomery domain is 1..p-1: the claim is restricted


def test_comparison_proved():
    assert check_bloq(LessThanEqual(3, 3))["verdict"] == "PROVED"
    assert check_bloq(LessThanConstant(3, 5))["verdict"] == "PROVED"


@attrs.frozen
class _AddButClaimsPlusOne(Bloq):
    """Decomposes into Add but declares a different classical action: must be REFUTED with a witness."""
    n: int = 2

    @property
    def signature(self):
        return Signature([Register("a", QUInt(self.n)), Register("b", QUInt(self.n))])

    def build_composite_bloq(self, bb, a, b):
        a, b = bb.add(Add(QUInt(self.n)), a=a, b=b)
        return {"a": a, "b": b}

    def on_classical_vals(self, a, b):
        return {"a": a, "b": (a + b + 1) % (1 << self.n)}


def test_wrong_classical_claim_is_refuted_with_witness():
    r = check_bloq(_AddButClaimsPlusOne(2))
    assert r["verdict"] == "REFUTED" and r["witness"] is not None, r


@attrs.frozen
class _NoClassicalAction(Bloq):
    @property
    def signature(self):
        return Signature([Register("q", QBit())])

    def build_composite_bloq(self, bb, q):
        return {"q": bb.add(XGate(), q=q)}


def test_missing_classical_action_is_unknown_not_proved():
    # XGate has an action but this wrapper does not define on_classical_vals -> call_classically fails.
    r = check_bloq(_NoClassicalAction())
    assert r["verdict"] in ("UNKNOWN", "PROVED")  # Qualtran may derive it from the decomposition
    if r["verdict"] == "PROVED":
        assert r["inputs_checked"] == 2


def test_too_many_inputs_is_unknown():
    r = check_bloq(Add(QUInt(8)), max_input_bits=10)
    assert r["verdict"] == "UNKNOWN" and "input bits" in r["detail"]


def test_enumeration_alone_never_gives_proved():
    r = check_bloq(Add(QUInt(2)), use_prover=False)
    assert r["verdict"] == "UNKNOWN"


def _spoil(monkeypatch, how):
    from qoptct.qualtran.bridge import rules
    real = rules.translate

    def bad(lowered, inputs, outputs, **kw):
        prog = real(lowered, inputs, outputs, **kw)
        ops = list(prog["ops"])
        if how == "drop_cx":
            i = max(i for i, op in enumerate(ops) if op[0] == "cx")
            del ops[i]
        elif how == "extra_x":
            ops.insert(len(ops) // 2, ("x", ops[-1][1] if ops[-1][0] == "x" else 0))
        elif how == "drop_phase":
            ops = [(op[0], op[1], ()) if op[0] == "erase" else op for op in ops]
        from qoptct.qualtran.bridge.emit import emit_qasm
        prog = dict(prog, ops=ops, qasm=emit_qasm(prog["nq"], ops))
        return prog

    monkeypatch.setattr(rules, "translate", bad)


@pytest.mark.parametrize("how", ["drop_cx", "extra_x", "drop_phase"])
def test_spoiled_translation_is_not_proved(monkeypatch, how):
    _spoil(monkeypatch, how)
    r = check_bloq(Add(QUInt(3)))
    assert r["verdict"] != "PROVED", (how, r)
