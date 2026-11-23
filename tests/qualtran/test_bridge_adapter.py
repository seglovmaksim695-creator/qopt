"""Checks of the strict Qiskit adapter.
"""

from __future__ import annotations
import pytest
pytest.importorskip("qiskit")

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import Gate, Parameter
from qiskit.circuit.library import CCXGate, CXGate
from qiskit.quantum_info import Operator

from qoptct.qualtran.bridge.from_qiskit import TranslationError, lower


def selftest() -> dict:
    records: list[dict] = []

    def accept(label: str, circuit: QuantumCircuit, expected: list[tuple[str, list[int]]] | None = None) -> None:
        record = {"case": label, "expected": "ACCEPT"}
        try:
            result = lower(circuit)
            got = [(leaf["name"], leaf["qubits"]) for leaf in result["leaves"]]
            deviation = float(np.max(np.abs(Operator(circuit).data - Operator(result["flat_source"]).data)))
            ok = deviation < 1e-12 and (expected is None or expected == got)
            record.update(observed="ACCEPT", passed=ok, matrix_deviation=deviation, leaves=got)
        except Exception as exc:
            record.update(observed=type(exc).__name__, passed=False, detail=str(exc))
        records.append(record)

    def reject(label: str, circuit: QuantumCircuit, **kwargs) -> None:
        record = {"case": label, "expected": "REJECT"}
        try:
            result = lower(circuit, **kwargs)
            record.update(observed="ACCEPT", passed=False, leaf_count=len(result["leaves"]))
        except TranslationError as exc:
            record.update(observed="REJECT", passed=True, detail=str(exc))
        except Exception as exc:
            record.update(observed=type(exc).__name__, passed=False, detail=str(exc))
        records.append(record)

    q = QuantumCircuit(3)
    q.x(2)
    q.cx(2, 0)
    q.ccx(0, 2, 1)
    q.swap(1, 2)
    accept("четыре_примитива", q, [("x", [2]), ("cx", [2, 0]), ("ccx", [0, 2, 1]), ("swap", [1, 2])])

    inner = QuantumCircuit(2, name="inner")
    inner.x(0)
    inner.cx(0, 1)
    outer = QuantumCircuit(3)
    outer.append(inner.to_gate(), [2, 0])
    accept("вложенные_провода", outer, [("x", [2]), ("cx", [2, 0])])

    q = QuantumCircuit(2)
    q.append(CXGate(ctrl_state=0), [0, 1])
    accept("управление_CX_по_нулю", q, [("x", [0]), ("cx", [0, 1]), ("x", [0])])
    q = QuantumCircuit(3)
    q.append(CCXGate(ctrl_state=1), [2, 0, 1])
    accept("смешанное_управление_CCX", q, [("x", [0]), ("ccx", [2, 0, 1]), ("x", [0])])

    q = QuantumCircuit(2)
    q.cx(0, 1)
    q.cx(0, 1)
    accept("нет_сокращения_обратных_вентилей", q, [("cx", [0, 1]), ("cx", [0, 1])])

    forged = Gate("x", 1, [])
    definition = QuantumCircuit(1)
    definition.h(0)
    forged.definition = definition
    q = QuantumCircuit(1)
    q.append(forged, [0])
    reject("поддельное_имя_x_с_H", q)

    q = QuantumCircuit(1)
    q.append(Gate("x", 1, []), [0])
    reject("непрозрачный_gate_с_именем_x", q)

    phased = QuantumCircuit(1, global_phase=0.125, name="phase_inside")
    phased.x(0)
    q = QuantumCircuit(1)
    q.append(phased.to_gate(), [0])
    reject("вложенная_глобальная_фаза", q)
    q = QuantumCircuit(1, global_phase=0.25)
    q.x(0)
    reject("глобальная_фаза_корня", q)

    q = QuantumCircuit(1)
    q.rx(np.pi / 2, 0)
    reject("дробный_поворот", q)
    q = QuantumCircuit(1)
    q.rx(Parameter("theta"), 0)
    reject("свободный_параметр", q)

    for name in ("h", "z", "s", "t", "id", "reset"):
        q = QuantumCircuit(1)
        getattr(q, name)(0)
        reject(f"запрещённый_{name}", q)
    q = QuantumCircuit(1)
    q.barrier()
    reject("barrier_не_исчезает", q)

    q = QuantumCircuit(1, 1)
    q.measure(0, 0)
    reject("измерение", q)
    q = QuantumCircuit(1, 1)
    with q.if_test((q.clbits[0], 1)):
        q.x(0)
    reject("классическое_управление", q)

    q = QuantumCircuit(2)
    q.append(Gate("unknown_atomic", 2, []), [0, 1])
    reject("неизвестный_атомарный_gate", q)

    empty = Gate("empty_custom", 1, [])
    empty.definition = QuantumCircuit(1)
    q = QuantumCircuit(1)
    q.append(empty, [0])
    reject("пустое_определение_gate", q)

    q = QuantumCircuit(2)
    q.x(0)
    q.cx(0, 1)
    reject("лимит_операций", q, max_operations=1)

    a = QuantumCircuit(1, name="a")
    a.x(0)
    b = QuantumCircuit(1, name="b")
    b.append(a.to_gate(), [0])
    c = QuantumCircuit(1)
    c.append(b.to_gate(), [0])
    reject("лимит_глубины", c, max_depth=1)

    return {
        "utc": datetime.now(timezone.utc).isoformat(),
        "total": len(records),
        "passed": sum(record["passed"] for record in records),
        "failed": sum(not record["passed"] for record in records),
        "cases": records,
    }



def test_qiskit_adapter():
    report = selftest()
    assert report["failed"] == 0 and report["total"] > 0, report
