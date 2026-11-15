# Changelog

## 0.4.0 — 2026-11-15

### Added
* **Exact prover `qoptct.verify`** and the subcommand **`qoptct prove CIRCUIT REF.qc [--clean] [--budget S]`**
  (JSON output; exit code 0 PROVED, 1 REFUTED, 2 UNKNOWN, 3 usage or input error); CIRCUIT is OpenQASM 3 or a
  `.qc` file. The proof covers every basis input and every
  combination of measurement outcomes (sum-over-paths normal form or trace miter over Z[ω]) in the unitary, strict
  and clean models; counterexamples are confirmed by an independent sparse simulator. MQT QCEC is used for unitary
  circuits when installed (`pip install "qoptct[qcec]"`). Fast sampled checkers: `check_dense` (also XOR
  conditions and resets), `check_sparse`, `check_big`, `quick_check`. The prover runs in its own process with a hard
  time limit, so it works from scripts without a `__main__` guard.
* `qoptct optimize --prove`: prove the written result exactly (status in the report).
* Optional candidate sources for `--mode balanced` / `--mode t` (all included in the package):
  * `qoptct.contrib.arith` — recognized classical arithmetic (adders, comparators, addition mod 2^n−1, hidden
    weighted bit) rebuilt on logical ANDs;
  * `qoptct.contrib.bandc` — logical-AND compiler for large classical circuits (up to 1600 wires);
  * `qoptct.contrib.gf2m` — GF(2^m) multipliers from low-rank bilinear decompositions (tables shipped as package
    data).
* CX post-pass `qoptct.contrib.cxpass` for unitary and measurement-based results (block re-synthesis, merged and
  XOR-conditioned corrections); `QoptPass(cxpass=True)` in the Qiskit plugin.
* Options `--no-arith`, `--no-bandc`, `--no-gf2m`, `--no-cxpass`, `--no-contrib`, `--budget`, `--contrib-sec`,
  `--source-sec NAME=SEC`, `--require-proof`. The JSON report has a `final` section (format, T, CX, extra qubits,
  measurements, how the result was obtained and checked) and a `modules` log.
* Measurement-based modes now also produce results above 22 qubits (from the modules; proven, or sampled-checked when
  the proof is out of reach and `--require-proof` is not given).
* Extra `qcec` (`mqt.qcec`); pytest marker `slow` (CI runs `-m "not slow"`).

### Changed
* Every module call runs in a child process with a hard wall-clock limit; an exception, a crash, a hang or garbage
  output of a module only skips it. Every circuit a module returns is checked again by the core.
* Measurement-based results are re-checked more strongly before they are written: exact proof, and if it ends
  UNKNOWN a sampled check with 64–256 random inputs and fresh measurement outcomes (the 30-sample check alone could
  miss a dropped conditional correction behind rarely sampled outcomes). A REFUTED result is never written; the
  optimizer falls back to the proven unitary result.
* The 0.3 pipeline is unchanged and the modules only add candidates (the minimum over all is written). One exception: the variant
  search of the logical-AND compiler now stops after `--and-sec` seconds (default 1800, `$QOPT_AND_SEC`; the
  unshuffled variant is always compiled completely). Before, it had no time limit and could run for hours on
  inputs such as ham15-low (0.3.0: balanced and t modes did not finish within 90 minutes).
* New overall time limit `--time-limit` (default 3600 s, `$QOPT_TIME_LIMIT`, 0 = none). All searches stop there —
  the unitary portfolio (rewriting programs, phase folding and their proofs), the module stages and the
  measurement-based searches — and the best result checked so far is written (at worst the input itself); only the
  final check of a written measurement-based result (or `--prove`) may follow. Under a limit the rewriting programs
  stop a quarter of the limit (at most 600 s) early, so the best of them can still be proven. The optional modules now run before the open-ended AND-compiler search, so their
  results survive the limit. 0.3.0 ran ham15-high (balanced/t) for more than 90 minutes.

### Comparison with 0.3.0
29 benchmarks (the 23 of `bench_qc/` and six larger ones) in the modes cx, balanced, t and t --clean, same time limit
for both versions: 0.4.0 is never worse by the objective of the mode; 55 results are better, 55 equal, and in 6 runs
(ham15-low and ham15-high, measurement-based modes) only 0.4.0 finished. Examples (T/CX): gf2^10_mult t 410/537 →
132/506, mod_adder_1024 t --clean 1011/979 → 36/78, qcla_adder_10 t 162/163 → 100/194, adder_8 cx 211/248 → 211/245.

How to read the comparison: in the measurement-based modes 0.3.0 has no result on the inputs above 22 qubits (all its
measurement stages need the 22-qubit simulator; its balanced/t runs there end in a timeout or return the unitary
circuit), so for those rows (adder_8, csum_mux_9, gf2^8/gf2^10_mult, mod_adder_1024, qcla_adder_10, qcla_com_7,
qcla_mod_7) the baseline is the unitary (cx-mode) result of 0.3.0. Exact proofs: every 0.4.0 result with measurements
in the comparison (81 of 85) was proven exactly by `qoptct prove`, except four large module results
(mod_adder_1024 and qcla_mod_7, modes balanced and t, 37–38 qubits): there the proof ran out of time and they were
accepted on sampled checks only (sparse simulator, 32–64 random inputs with random measurement outcomes) — checked, not
proven. Unitary results above 19 qubits were proven by the optimizer itself (PyZX, numerical check or path-sum prover).

### Fixed
* `qoptct prove CIRCUIT.qc REF.qc` failed with `ModuleNotFoundError` (exit code 1, the code of REFUTED).
* Input errors (missing or unreadable file, malformed or empty circuit, unknown gate, undeclared wire, missing output
  directory) are reported as one `error: ...` line on stderr with exit code 3 instead of a traceback with exit code 1;
  for `qoptct prove` and `python -m qoptct.verify` an input error is never reported as REFUTED. Usage errors of
  `qoptct` also exit with 3 (2 is UNKNOWN for `qoptct prove`).
* A unitary result written to `-o *.qasm` is now OpenQASM 3 (before, `.qc` text was written into the `.qasm` file);
  the report's `final.file` names the format. A measurement-based result written to a non-`.qasm` name gives a
  warning on stderr.
* The written OpenQASM 3 is readable by `qiskit.qasm3.loads`: XOR (`^`) conditions of the CX post-pass and AND
  (`&&`) conditions of older stages are rewritten into single-bit (possibly nested) `if`s — the condition is computed
  into an idle measured qubit and measured once more (same T and CX, one more measurement) or the correction is
  repeated per term — and the rewritten circuit is proven again (otherwise the previous readable result is written).
  On the regression set this changed qft_4 (balanced, t, t --clean: same T/CX, +1 measurement) and csla_mux_3
  (balanced, t: 36/80 → 36/81, the measured-parity form could not be proven there).
* `qoptct prove` / `python -m qoptct.verify`: syntactically malformed OpenQASM (also binary files, out-of-range
  operands, conditions on unmeasured bits) is an input error, exit code 3 (was 2, UNKNOWN).
* Intermediate files of the stages (`*.and.qasm`, `*.sel.qasm`, `*_tg3.qasm`, ...) are no longer left next to `-o`
  (new option `--keep-temp`); `-o` is written atomically and never holds an intermediate circuit during the run.
* Unitary results of more than 26 qubits whose PyZX proof does not finish in time are now proven with the exact
  path-sum prover of `qoptct.verify`, instead of falling back to the input (the cx result of large inputs depended on
  the machine load).
* `verify_big` (sampled path-sum checker for > 22 qubits): normalization overflow on circuits with many H gates
  (the check could become empty), false rejections of large correct circuits, and a false acceptance in the sampled
  mode (the branch norm is now computed exactly).
* sdist: `global-exclude *hwb*` also dropped the module `contrib/arith/hwb.py`; only benchmark files are excluded now.

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
