# Changelog

## 0.2.0 — 2026-11-12 (first public release candidate; not yet published)

### Added
* Measurement-based optimization modes `--mode balanced` and `--mode t` (OpenQASM 3 output with mid-circuit
  measurements and classically controlled Clifford corrections):
  * logical-AND compiler with measurement-based uncomputation and exact GF(2) tracking of wire values;
  * clean-ancilla model `--clean` (wires not listed in `.i` start in |0>);
  * CX polishing of AND circuits (re-synthesis of unconditioned pieces using |+> and measured-wire freedoms);
  * hybrid compilation for circuits with genuine H gates (grover_5, qft_4);
  * recognition of multi-controlled Toffolis (AND ladders, Khattar–Gidney syntheses with AND-based outer pair);
  * Hadamard gadgets + TODD with deterministic corrections (`--tgad`), partial gadgetization;
  * optional external VarTODD via `--vartodd` and `QOPT_VARTODD_BIN`.
* Independent checkers: `verify_qasm` (strict and clean models, no renormalization) and `verify_big`
  (sum-over-paths, > 22 qubits).
* `results/`: verified best and Pareto-optimal circuits for 23 benchmarks in three models, `TABLE.md`,
  `table.csv`, comparison with nine tools (`COMPARISON.md`).
* `scripts/verify_all.py`, `scripts/reproduce.sh`, `scripts/fetch_benchmarks.sh`.
* Documentation: `docs/METHOD.md`, `docs/FINDINGS.md`, `paper/preprint.md`.

### Changed
* Installable package `qopt` with the `qopt` console command; all paths relative to the package; caches in
  `$QOPT_CACHE` (default `~/.cache/qopt`).
* Cached data that is keyed by benchmark name (unitary start circuits, hand-made clean-ancilla constructions) is
  used only if the SHA-256 of the input file matches the benchmark.
* CLI messages and help in English. The `bench` subcommand of 0.1 was replaced by `scripts/reproduce.sh`.

### Removed / not included
* `hwb6` (and its optimized circuits): the benchmark file forbids redistribution; `scripts/fetch_benchmarks.sh`
  downloads it for local use.

### Fixed (found by independent verification, see docs/FINDINGS.md)
* A 51-T hwb6 circuit of an earlier gadget pipeline was wrong; it was never part of a release.

## 0.1 — 2026-11-10 (internal)
* Unitary CX/T optimizer: portfolio of evolved rewriting programs, phase folding, neural block resynthesis,
  proof of equality (unitary / PyZX / numerical) with fallback to the input.
