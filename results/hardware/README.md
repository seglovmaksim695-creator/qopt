# Hardware run: IBM Quantum ibm_marrakesh (2026-11-13)

One Sampler V2 job (`dbr5s71abdes73fpchmg`), 60 circuits × 4000 shots, 67 s QPU time.
For each benchmark, four computational-basis inputs (two that the circuit changes, two that it leaves unchanged) were
run in three variants with identical transpilation (`optimization_level=1`, seed 7, same initial layout for the data
qubits):

* **A** — the input `.qc` circuit decomposed to Clifford+T,
* **B** — `qopt` unitary result (`results/best/unitary`),
* **C** — `qopt` measurement-based result (`results/best/strict`; mid-circuit measurements + `if_else`).

Metric: fraction of shots that return the correct output bit string (mean over the four inputs; statistical error ≈ ±0.004).

| circuit | success A → B → C | two-qubit gates after transpilation A → B → C |
|---|---|---|
| tof_3 | 0.769 → **0.789** → 0.554 | 30 → 21 → 28 |
| barenco_tof_3 | 0.682 → **0.766** → 0.576 | 42 → 26 → 27 |
| mod5_4 | 0.730 → **0.834** → 0.745 | 52 → 22 → 32 |
| tof_4 | 0.616 → **0.681** → 0.467 | 51 → 27 → 42 |
| barenco_tof_4 | 0.567 → **0.665** → 0.410 | 90 → 42 → 40 |

* The unitary `qopt` circuits (B) give higher success on all five benchmarks (+2 to +10 percentage points; better on
  17 of 20 individual inputs), consistent with 30–58 % fewer two-qubit gates after transpilation.
* The measurement-based circuits (C) were on par with B in a noisy simulation (Aer with the device noise model), but
  worse than A on the device. A likely reason (not separately tested) is that mid-circuit measurement and classical
  feed-forward take microseconds during which the other qubits decohere — this is not captured by the noise model —
  plus extra routing for the additional ancillas. On current NISQ hardware the measurement-based modes are therefore
  not advantageous; they target fault-tolerant cost models, where T gates dominate.

Raw counts (simulated and measured), layouts and gate counts: `ibm_marrakesh_2026-11-13.json`; experiment script: `exp.py`.
