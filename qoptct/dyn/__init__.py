# SPDX-License-Identifier: Apache-2.0
"""qoptct.dyn: measurement-based erasure (dynamic circuits) versus ordinary reverse computation on noisy hardware.

    gen      arithmetic blocks in the AND / MZ form (adders, comparators, modular adders, lookups, multipliers)
    lower    the two lowerings to Qiskit + the ideal reference + an exact gate-level check
    hw       hardware cost model from a backend calibration: schedule (with measurement and feed-forward delays),
             estimated fidelity
    sim      noisy trajectory simulation (qiskit-aer) of a scheduled circuit
    choose   choose(block, backend): compile the variants, estimate, pick the best; DynErasePass for Qiskit circuits

Default: ordinary (static) uncomputation, always. Measurement-based erasure (M) and the hybrid (H) are available only
by an explicit option (allow_dynamic=True / --allow-dynamic), which prints a warning.

Check on a real chip (ibm_kingston, 20.11.2026, one job per experiment): measurement-based erasure lost to ordinary
uncomputation on all blocks checked (adder2: model 0.918, measured 0.525; the model badly underestimates the cost of a
measurement in the middle of a circuit), and the hybrid lost too (lookadd2x2, X test: 0.682 against 0.772). The hybrid
won only in the model. See docs/real_runs/2026-11-20_ibm_kingston/.

CLI: qoptct dyn --help  (or python -m qoptct.dyn --help)
"""
