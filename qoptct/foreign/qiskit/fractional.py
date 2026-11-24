# SPDX-License-Identifier: Apache-2.0
"""Q1c: which Qiskit library circuits translate (no rotations) and which need fractional phases (Q2)."""
from __future__ import annotations

import warnings

from qiskit import QuantumCircuit
from qiskit.circuit.library import (CDKMRippleCarryAdder, DraperQFTAdder, ExactReciprocalGate, FullAdderGate,
                                    HRSCumulativeMultiplier, HalfAdderGate, IntegerComparator, LinearPauliRotations,
                                    ModularAdderGate, MultiplierGate, PhaseEstimation, PiecewiseLinearPauliRotations,
                                    PolynomialPauliRotations, QFTGate, QuadraticForm, RGQFTMultiplier,
                                    VBERippleCarryAdder, WeightedAdder)

from .translate import NeedsFractionalPhases, QiskitTranslationError, translate


def _gate_circuit(g):
    qc = QuantumCircuit(g.num_qubits)
    qc.append(g, range(g.num_qubits))
    return qc


def _native_control(c):
    g = c.to_gate().control(1)
    return _gate_circuit(g)


def candidates(n=4):
    from qiskit.circuit.library import IntegerComparatorGate
    return [
        ('DraperQFTAdder(half)', lambda: DraperQFTAdder(n, kind='half')),
        ('DraperQFTAdder(fixed)', lambda: DraperQFTAdder(n, kind='fixed')),
        ('RGQFTMultiplier', lambda: RGQFTMultiplier(n)),
        ('HRSCumulativeMultiplier', lambda: HRSCumulativeMultiplier(n)),
        ('MultiplierGate', lambda: _gate_circuit(MultiplierGate(n, n))),
        ('QFTGate', lambda: _gate_circuit(QFTGate(n))),
        ('QuadraticForm', lambda: QuadraticForm(n, quadratic=[[1, 0, 0, 0]] * 0 or None, linear=[1] * n, offset=1)),
        ('WeightedAdder', lambda: WeightedAdder(n, [1] * n)),
        ('LinearPauliRotations', lambda: LinearPauliRotations(n, 0.5, 0.1)),
        ('PiecewiseLinearPauliRotations', lambda: PiecewiseLinearPauliRotations(n, [0, 4], [1, 1], [0, 0])),
        ('PolynomialPauliRotations', lambda: PolynomialPauliRotations(n, [0, 1, 1])),
        ('PhaseEstimation', lambda: PhaseEstimation(3, QuantumCircuit(1))),
        ('ExactReciprocalGate', lambda: _gate_circuit(ExactReciprocalGate(n, 0.5))),
        ('CDKMRippleCarryAdder.control(1) (Qiskit native)', lambda: _native_control(CDKMRippleCarryAdder(n))),
        ('VBERippleCarryAdder.control(1) (Qiskit native)', lambda: _native_control(VBERippleCarryAdder(n))),
        ('HalfAdderGate.control(1)', lambda: _gate_circuit(HalfAdderGate(n).control(1))),
        ('FullAdderGate.control(1)', lambda: _gate_circuit(FullAdderGate(n).control(1))),
        ('ModularAdderGate.control(1)', lambda: _gate_circuit(ModularAdderGate(n).control(1))),
        ('MultiplierGate.control(1)', lambda: _gate_circuit(MultiplierGate(n, n).control(1))),
        ('IntegerComparator.control(1)', lambda: _native_control(IntegerComparator(n, 3))),
        ('IntegerComparatorGate.control(1)', lambda: _gate_circuit(IntegerComparatorGate(n, 3).control(1))),
    ]


def classify(n=4):
    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for name, build in candidates(n):
            try:
                tr = translate(build())
                rows.append((name, 'translates', 'ops=%s' % tr.counts()))
            except NeedsFractionalPhases as e:
                rows.append((name, 'needs-fractional-phases', str(e)))
            except QiskitTranslationError as e:
                rows.append((name, 'unsupported', str(e)[:160]))
            except Exception as e:
                rows.append((name, 'build-error', '%s: %s' % (type(e).__name__, str(e)[:120])))
    return rows
