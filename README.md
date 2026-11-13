# qopt — a verified Clifford+T circuit optimizer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24182059.svg)](https://doi.org/10.5281/zenodo.24182059)

`qopt` reduces the **T count** and the **CNOT (CX) count** of Clifford+T circuits given as Toffoli networks
(Feynman `.qc` format). It has a unitary mode and measurement-based modes, and it **never reports a circuit that
has not been checked**: unitary outputs are proven equal to the input (or the input is returned), and every
measurement-based output is re-checked by an independent simulator before it is written.

Repository: <https://github.com/seglovmaksim695-creator/qopt>  · Version 0.2.0 · License Apache-2.0 (benchmark files: see `bench_qc/SOURCES.md`) · Author: Maksym Shchehlov

## Installation

```bash
git clone https://github.com/seglovmaksim695-creator/qopt.git && cd qopt
python3 -m venv .venv && . .venv/bin/activate        # Python 3.12+ (tested with 3.14.6 on macOS arm64)
pip install -r requirements.txt                      # exact versions used for the results
pip install -e .                                     # provides the `qopt` command
pytest -q tests                                      # ~1 minute
```

`requirements-lock.txt` pins the full environment used for testing. Caches are written to `$QOPT_CACHE`
(default `~/.cache/qopt`).

## Usage

```bash
qopt optimize bench_qc/tof_3.qc -o tof_3_opt.qc                          # unitary, minimum CX (default --mode cx)
qopt optimize bench_qc/tof_3.qc -o tof_3_opt.qc --net --post             # + neural block resynthesis + CX post-pass
qopt optimize bench_qc/tof_3.qc -o tof_3.qasm --mode balanced            # fewer T, CX <= 1.5x unitary CX
qopt optimize bench_qc/mod_red_21.qc -o m.qasm --mode t --clean          # minimum T, clean-ancilla model
qopt optimize X.qc --report report.json                                  # JSON report (also printed to stdout)
```

Python: `names, ref, best, how = qopt.optimize("X.qc")`. Qiskit: `from qopt.qiskit_qopt import QoptPass`
(a `TransformationPass`; h x z s sdg t tdg cx cz ccx ccz).

| option | meaning |
|---|---|
| `--mode cx` (default) | unitary output (`.qc`), minimum CX; proven equal to the input, otherwise the input is returned |
| `--mode balanced` | measurement-based output (OpenQASM 3) with lower T and CX ≤ `--cx-factor` (1.5) × unitary CX, if one is found and verified |
| `--mode t` | lowest T found (measurement-based), CX may grow |
| `--clean` | clean-ancilla model (see below) for the measurement-based modes |
| `--net` | neural block resynthesis (models in `qopt/models/`) |
| `--post` | CX post-pass on H-free regions (slow) |
| `--meas` | older Gidney transformation (clean-ancilla model), kept for compatibility |
| `--tgad`, `--vartodd` | Hadamard gadgets + TODD; `--vartodd` also tries an external VarTODD binary given by `QOPT_VARTODD_BIN` |

Measurement-based outputs are written only if they pass `verify_qasm` (≤ 22 qubits in total); otherwise the unitary
result is written. The JSON report says which candidate was used.

## Cost models — read this before comparing numbers

The same benchmark has different "best" numbers in different models. All tables say which model they use.

| model | allowed in the output | input wires | ancillas | checked by |
|---|---|---|---|---|
| **unitary** | Clifford+T gates only, same wires | arbitrary state | none | full unitary / PyZX / numerical |
| **strict** (measurement-based) | + fresh ancillas in \|0>, mid-circuit measurements, classically controlled Clifford gates | **all** `.qc` wires are data in an arbitrary state, including wires the benchmark intends as ancillas | fresh ones must end in a basis state, unentangled | `verify_qasm` (strict), `verify_big` |
| **clean** (clean ancillas) | as strict | wires listed in the `.i` line are data; **the other `.qc` wires start in \|0> and must end in \|0>** (the intended use of e.g. `tof_n`, `mod_red_21`) | as strict | `verify_qasm(clean=True)` |

* *unitary* is the model of most circuit-optimization papers (Nam et al., PyZX, GUOQ, Quarl, …).
* *strict* is stronger than *clean*: a strict circuit is valid in the clean model, not vice versa.
* In both measurement-based models the logical AND costs 4 T and its uncomputation 0 T (Gidney 2018); that is
  why T counts there are far below unitary ones. They must not be compared with unitary T-count records.
* Published H-gadget T-count records (TODD/FastTODD, AlphaTensor-Quantum) use yet another model (gadget ancillas,
  no classical feed-forward on the output). They are listed for orientation only.

## Results (23 benchmarks)

| circuit | qubits | Toffolis | input T/CX | unitary T/CX | strict T/CX (+extra) | clean T/CX (+extra) |
|---|---|---|---|---|---|---|
| adder_8 | 24 | 57 | 399/409 | 193/223 | — | — |
| barenco_tof_10 | 19 | 32 | 224/192 | 72/72 | 59/48 (+3) | 59/48 (+3) |
| barenco_tof_3 | 5 | 4 | 28/24 | 16/14 | 8/9 (+2) | 8/9 (+2) |
| barenco_tof_4 | 7 | 8 | 56/48 | 24/24 | 12/13 (+3) | 12/13 (+3) |
| barenco_tof_5 | 9 | 12 | 84/72 | 32/32 | 16/17 (+4) | 16/17 (+4) |
| csla_mux_3 | 15 | 10 | 70/80 | 64/64 | 40/70 (+1) | 28/84 (+3) |
| csum_mux_9 | 30 | 28 | 196/168 | 84/112 | — | — |
| gf2^4_mult | 12 | 16 | 112/99 | 68/86 | 58/115 (+0) | 36/79 (+1) |
| gf2^5_mult | 15 | 25 | 175/154 | 115/132 | 90/183 (+0) | 56/124 (+1) |
| grover_5 | 9 | 48 | 336/288 | 166/160 | 124/155 (+4) | 124/155 (+4) |
| ham15-low | 17 | 23 | 161/236 | 97/179 | — | 68/177 (+1) |
| mod5_4 | 5 | 4 | 28/28 | 7/10 | 4/11 (+1) | 4/11 (+0) |
| mod_mult_55 | 9 | 7 | 49/48 | 35/39 | 20/35 (+2) | 12/24 (+1) |
| mod_red_21 | 11 | 17 | 119/105 | 73/61 | 44/58 (+5) | 28/48 (+0) |
| qcla_com_7 | 24 | 29 | 203/186 | 95/111 | 59/565 (+37) | 59/565 (+37) |
| qcla_mod_7 | 26 | 59 | 413/382 | 237/256 | — | — |
| qft_4 | 5 | 2 | 69/46 | 67/37 | 53/265 (+42) | 53/265 (+42) |
| rc_adder_6 | 14 | 11 | 77/93 | 47/65 | 24/54 (+6) | 24/53 (+6) |
| tof_10 | 19 | 17 | 119/102 | 71/54 | 59/71 (+3) | 36/37 (+1) |
| tof_3 | 5 | 3 | 21/18 | 15/12 | 8/10 (+2) | 8/9 (+1) |
| tof_4 | 7 | 5 | 35/30 | 23/18 | 12/15 (+3) | 12/13 (+1) |
| tof_5 | 9 | 7 | 49/42 | 31/24 | 16/20 (+4) | 16/17 (+1) |
| vbe_adder_3 | 10 | 10 | 70/70 | 24/28 | 12/33 (+3) | 12/33 (+3) |

Best point per model: unitary = minimum CX; strict and clean = minimum T (then CX, then extra qubits). Other Pareto points (e.g. fewer extra qubits) are in results/TABLE.md.

Totals over the circuits that have a result in the model: unitary: 23 circuits, T 3093→1656, CX 2920→1813; strict: 19 circuits, T 1924→718, CX 1725→1747; clean: 20 circuits, T 2085→671, CX 1961→1781

Full tables with Pareto fronts and provenance: `results/TABLE.md`, machine-readable `results/table.csv`;
circuits: `results/best/<model>/<name>.qasm` and `results/pareto/<model>/`.

## Comparison with nine tools (summary; details in `results/COMPARISON.md`)

Tools run on the same files: Qiskit 2.6.1, PyZX 0.10.8, pytket 2.18.5, Cirq 1.7, Q# Resource Estimator 1.31,
Qualtran 0.7.0, pyLIQTR 1.0.1, Qrisp 0.9.10, PennyLane 0.45.1. All produced circuits were verified with the same
checker; estimates without circuits are marked as such.

| circuit | best tool circuit by T (as-is or MCX rewrite) | best tool circuit by CX | qopt unitary | qopt strict | qopt clean |
|---|---|---|---|---|---|
| barenco_tof_3 | 8/9 +2 (PennyLane 0.45.1, and_ladder; strict) | 8/9 +2 (PennyLane 0.45.1, and_ladder; strict) | 16/14 | 8/9 (+2) | 8/9 (+2) |
| barenco_tof_4 | 12/13 +3 (PennyLane 0.45.1, and_ladder; strict) | 12/13 +3 (PennyLane 0.45.1, and_ladder; strict) | 24/24 | 12/13 (+3) | 12/13 (+3) |
| barenco_tof_5 | 16/17 +4 (PennyLane 0.45.1, and_ladder; strict) | 16/17 +4 (PennyLane 0.45.1, and_ladder; strict) | 32/32 | 16/17 (+4) | 16/17 (+4) |
| csla_mux_3 | 62/161 +0 (PyZX 0.10.8, full_reduce; strict) | 64/72 +0 (Qiskit 2.6.1, as-is; strict) | 64/64 | 40/70 (+1) | 28/84 (+3) |
| gf2^4_mult | 68/274 +0 (PyZX 0.10.8, full_reduce; strict) | 96/99 +0 (Qiskit 2.6.1, as-is; strict) | 68/86 | 58/115 (+0) | 36/79 (+1) |
| gf2^5_mult | 115/509 +0 (pytket 2.18.5, zxperm_fpo; strict) | 155/154 +0 (Qiskit 2.6.1, as-is; strict) | 115/132 | 90/183 (+0) | 56/124 (+1) |
| grover_5 | 166/301 +0 (PyZX 0.10.8, full_reduce; strict) | 296/288 +0 (Qiskit 2.6.1, as-is; strict) | 166/160 | 124/155 (+4) | 124/155 (+4) |
| ham15-low | 97/322 +0 (PyZX 0.10.8, full_reduce; strict) | 147/236 +0 (Qiskit 2.6.1, as-is; strict) | 97/179 | — | 68/177 (+1) |
| mod5_4 | 8/20 +0 (PyZX 0.10.8, full_reduce; strict) | 8/20 +0 (PyZX 0.10.8, full_reduce; strict) | 7/10 | 4/11 (+1) | 4/11 (+0) |
| mod_mult_55 | 35/89 +0 (PyZX 0.10.8, full_reduce; strict) | 45/48 +0 (Qiskit 2.6.1, as-is; strict) | 35/39 | 20/35 (+2) | 12/24 (+1) |
| mod_red_21 | 73/149 +0 (PyZX 0.10.8, full_reduce; strict) | 107/105 +0 (Qiskit 2.6.1, as-is; strict) | 73/61 | 44/58 (+5) | 28/48 (+0) |
| qft_4 | 67/55 +0 (PyZX 0.10.8, full_reduce; strict) | 409/45 +0 (Qiskit 2.6.1, as-is; strict) | 67/37 | 53/265 (+42) | 53/265 (+42) |
| rc_adder_6 | 47/103 +0 (PyZX 0.10.8, full_reduce; strict) | 67/83 +0 (Qiskit 2.6.1, as-is; strict) | 47/65 | 24/54 (+6) | 24/53 (+6) |
| tof_10 | 36/37 +1 (PennyLane 0.45.1, and_ladder; clean) | 36/37 +1 (PennyLane 0.45.1, and_ladder; clean) | 71/54 | 59/71 (+3) | 36/37 (+1) |
| tof_3 | 8/9 +1 (PennyLane 0.45.1, and_ladder; clean) | 8/9 +1 (PennyLane 0.45.1, and_ladder; clean) | 15/12 | 8/10 (+2) | 8/9 (+1) |
| tof_4 | 12/13 +1 (PennyLane 0.45.1, and_ladder; clean) | 12/13 +1 (PennyLane 0.45.1, and_ladder; clean) | 23/18 | 12/15 (+3) | 12/13 (+1) |
| tof_5 | 16/17 +1 (PennyLane 0.45.1, and_ladder; clean) | 16/17 +1 (PennyLane 0.45.1, and_ladder; clean) | 31/24 | 16/20 (+4) | 16/17 (+1) |
| vbe_adder_3 | 24/71 +0 (PyZX 0.10.8, full_reduce; strict) | 56/58 +0 (Qiskit 2.6.1, as-is; strict) | 24/28 | 12/33 (+3) | 12/33 (+3) |

Each tool cell: T/CX extra-qubits (tool, variant; model in which the circuit passed verification).

Honest caveats:

* **Given the gate list "as is", no tool replaces Toffolis by logical ANDs with measurement-based uncomputation
  automatically.** Unitary tools give ~7 T per Toffoli minus phase folding; Q#/Qualtran give only estimates
  (4 T per Toffoli) without a circuit.
* **When a human rewrites the circuit as a multi-controlled X** (tof_n, barenco_tof_n), tools reach our numbers or
  beat them: PennyLane's `TemporaryAND` ladder equals our clean-model tof_n points; on barenco_tof_n Qiskit's
  Khattar–Gidney syntheses with one ancilla were better than our general compiler until the MCX recognition
  (`contrib/c12`) was added. Our advantage is automatic discovery in arbitrary networks, not better gadgets.
* Third-party pitfalls and one wrong circuit of a tool are documented in `docs/FINDINGS.md`.

## How to check the results yourself

```bash
python scripts/verify_all.py --jobs 2          # re-verifies all 200 rows of results/table.csv (~40 min)
python scripts/verify_all.py --only tof_3       # one circuit
sh scripts/reproduce.sh cx out                  # re-run the optimizer on all benchmarks (see script header)
```

Every qasm in `results/` is self-contained OpenQASM 3 using `h x z s sdg t tdg cx cz`, `measure`, and `if` on
classical bits; data qubits are `q[0..n-1]` in the order of the `.v` line of the `.qc` file, extra qubits follow
and start in |0>.

## Limitations

* The checks are randomized: `verify_qasm` uses 30 random input states (10 for 20–22-qubit circuits) and sampled
  measurement outcomes; `verify_big` samples basis inputs and branches. A wrong circuit could pass only if it acts
  correctly on all samples. Unitary results ≤ 10 qubits are proven exactly; larger ones rely on PyZX or a
  numerical check.
* Measurement-based modes handle at most 22 qubits in total (state-vector check). Larger outputs are discarded.
* The neural models were trained on random circuits of 2–4 qubits; `--net` gives small gains (≈ 1–2 % CX).
* Runtime: unitary mode takes seconds for small circuits and up to an hour for the largest benchmarks; the
  measurement-based modes take minutes. Several rewriting programs have internal time limits, so results can
  depend on machine speed.
* Some numbers in `results/` come from longer development searches than a single default run (column `source`
  in `results/TABLE.md`); a default run may give slightly worse points.
* The hand-made clean-model constructions (`contrib/c1/best_clean`) apply only to the exact benchmark files.
* hwb6 is not included (license of the benchmark file).

## Citation

See `CITATION.cff`. A preprint draft is in `paper/preprint.md`.

## Acknowledgements

Benchmarks from the Feynman repository (M. Amy). Logical AND construction: C. Gidney, *Halving the cost of quantum
addition*, Quantum 2, 74 (2018). Sum-over-paths verification: M. Amy, arXiv:1805.06908.
