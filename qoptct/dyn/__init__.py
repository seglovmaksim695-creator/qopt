# SPDX-License-Identifier: Apache-2.0
"""qoptct.dyn: measurement-based erasure (dynamic circuits) versus ordinary reverse computation on noisy hardware.

    gen      arithmetic blocks in the AND / MZ form (adders, comparators, modular adders, lookups, multipliers)
    lower    the two lowerings to Qiskit + the ideal reference + an exact gate-level check
    hw       hardware cost model from a backend calibration: schedule (with measurement and feed-forward delays),
             estimated fidelity
    sim      noisy trajectory simulation (qiskit-aer) of a scheduled circuit
    choose   choose(block, backend): compile the variants, estimate, pick the best; DynErasePass for Qiskit circuits

CLI: python -m qoptct.dyn --help
"""
