# SPDX-License-Identifier: Apache-2.0
"""Qiskit -> qoptct bridge for circuits without non-Clifford+T rotations.

    from qoptct.foreign.qiskit import check_qiskit
    r = check_qiskit(circuit, spec)            # spec: Spec or classical function
    r = check_qiskit(circuit, reference=ref)   # miter against another Qiskit circuit
    r['status']  # PROVED | REFUTED (with witness) | UNKNOWN (r['needs'] == 'fractional_phases' for Q2 circuits)
"""
from .check import Spec, check_qiskit, classical_validate, run_classical
from .translate import NeedsFractionalPhases, QiskitTranslationError, Translation, translate

__all__ = ['check_qiskit', 'Spec', 'translate', 'Translation', 'NeedsFractionalPhases', 'QiskitTranslationError',
           'classical_validate', 'run_classical']
