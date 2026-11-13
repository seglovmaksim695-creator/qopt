# Changelog

## 0.3.0 — 2026-11-13 (first PyPI release)

### Changed — package renamed to `qoptct`
* **The distribution, the import package and the command are now called `qoptct`** (the PyPI name `qopt`
  belongs to an unrelated project with the import name `qopt`). The project itself is still called qopt.
  Migration: `pip install qoptct`; `import qopt` → `import qoptct`; `from qopt.qiskit_qopt import QoptPass` →
  `from qoptct.qiskit_qopt import QoptPass`; command `qopt optimize ...` → `qoptct optimize ...`.
  Environment variables (`QOPT_CACHE`, `QOPT_PROG_SEC`, ...) and the cache directory `~/.cache/qopt` are unchanged.
* PyTorch is now optional: `pip install "qoptct[net]"`. Only `--net` (and `QoptPass(net=True)`) needs it and
  stops with an installation hint if it is missing; the other modes no longer import torch at all.
* Python ≥ 3.11 (was 3.12); tested with 3.11, 3.12 and 3.14.

### Added
* Qiskit transpiler stage plugins, both named `qoptct`:
  `transpile(qc, optimization_method="qoptct")` (Qiskit's optimization stage of the chosen level, then
  `QoptPass` restricted to the target / basis / coupling map) and `init_method="qoptct"` (`QoptPass` at the
  Toffoli level, before decomposition).
* `QoptPass` now works on any circuit: it optimizes each maximal Clifford+T run separately (measurements,
  barriers, rotations and other gates are kept), replaces a run only if the result is proven equal and strictly
  cheaper, respects basis gates / coupling map / target, keeps the global phase (runs ≤ 10 qubits), and reports
  `runs`, `replaced`, CX/T before and after in `pass.last`. Messages and documentation in English.
* `.github/workflows/tests.yml` (pytest on Ubuntu and macOS against the built wheel) and
  `.github/workflows/publish.yml` (PyPI Trusted Publishing on a published GitHub release).
* `MANIFEST.in` (the sdist contains the tests and the redistributable benchmark inputs, never hwb6).

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
