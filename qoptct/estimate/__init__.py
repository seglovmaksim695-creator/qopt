# SPDX-License-Identifier: Apache-2.0
"""qoptct.estimate: physical (fault-tolerant) cost estimates from logical counts.

    from qoptct.estimate import estimate
    r = estimate('gidney2025', toffoli=6.5e9, qubits=1399, shots=9.2)
    print(r.physical_qubits, r.expected_runtime_s / 86400)

See qoptct.estimate.physical for the models, their formulas and the paper references."""
from .physical import (MODELS, PAPER, PhysicalEstimate, estimate, ge2019, ge2019_rsa_reference, gidney2025,
                       litinski2023)

__all__ = ["MODELS", "PAPER", "PhysicalEstimate", "estimate", "ge2019", "ge2019_rsa_reference", "gidney2025",
           "litinski2023", "main_estimate"]


def main_estimate(argv=None):
    from .cli import main_estimate as _m
    return _m(argv)
