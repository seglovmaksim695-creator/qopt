# Changelog

## 0.6.0 — unreleased

### Added
* **`qoptct.verify.algebra`** and the subcommand **`qoptct aprove CIRCUIT.qasm SPEC.json [--budget S] [--json]`**
  (also `qoptct kprove CIRCUIT.qasm SPEC.json --engine algebra`): a word-level algebraic exact prover for integer
  arithmetic on logical ANDs with measurement-based erasure, the class the BDD engine cannot handle (integer
  multipliers). The circuit is read by a strict parser, every carry chain gives an exact word identity over the
  integers, dropped top carries are resolved by interval analysis, and the claims of a JSON specification
  (`"out == x*y"`, `"r + p*(H + thi1 - 1 + g) == L + 2**n*H"`, `"r <= p-1"`) are proved by backward rewriting to the
  zero polynomial. PROVED is a statement for every basis input and every combination of measurement outcomes
  (classical function, |0> ancillas, no input-dependent phase). Proved at full size (`scripts/algebra_proofs.py`):
  `multiply` and `square` at n = 256 with Karatsuba levels 0..3 (4-9 s, at most 0.8 GB), the fold reductions for
  Curve25519, secp256k1, P-256, P-384, P-521 (under 1 s), and whole modular multiplications / squarings
  (multiplier followed by the reduction, `algebra.compose_modmul`) for Curve25519, secp256k1, P-256 (4-18 s).
  Limits: not a general equivalence checker (carry-chain arithmetic only; anything else ends in UNKNOWN or ERROR);
  some correct circuits are not proved (Karatsuba on very narrow or odd-width leaves, `square(256, 3, legacy=True)`);
  garbage registers are only shown to be classical functions of the input.
* **`qoptct.hw`** and the subcommand **`qoptct hw --backend fake_torino IN.qasm [-o OUT.qasm] [--report R.json]`**
  (experimental; extra `qoptct[hw]`): compilation for noisy IBM-like devices (heavy-hex, CZ/ECR + RZ, SX, X) with
  calibration data. A verified portfolio of Qiskit level 3, qoptct logical CX/T optimization + Qiskit, and a new
  router (commutation-aware dependencies, exact SWAP-absorption cost from the KAK decomposition, floating blocks,
  mirror rule, line layouts, calibration-aware relocation). Metrics: two-qubit gates, depth, expected fidelity from
  calibration; equivalence check with qubit maps (`qoptct.hw.check.verify_compiled`); benchmark driver
  `python -m qoptct.hw.bench`, tket comparison `scripts/hw_tket.py`, real-run preparation `scripts/hw_real_run.py`
  (never submits by itself). Results and limitations: `docs/HW.md`.
* **`qoptct compose`** (`python -m qoptct.verify.compose`, `qoptct.verify.compose`): compositional exact proof of the
  Kaliski inversion at n = 256 on the real circuit text (10 configurations: secp256k1, P-256, a random prime; generic /
  special reduction; forward / clean; also Curve25519, P-384, P-521 clean + special), 40-250 s, about 0.5 GB. A generic
  segment engine (`walker.py`, contracts on ROBDDs, untrusted labels from the generator) and Z3 verification conditions
  for the lemmas about the integer algorithm (`lemmas_z3.py`, extra `qoptct[z3]`, needs the `z3` binary). Machine
  checked: segments, gadget lemmas, the Z3 lemmas; written argument (`PROOF.md`): the gluing of segments and the link
  between the Z3 statements and the contracts. `kaliski.py` only gained segment labels (the circuit text is unchanged)
  and a streamed `write_qasm`.
* **`qoptct dyn`** (`qoptct.dyn`, extra `qoptct[dyn]`, experimental): measurement-based erasure versus ordinary reverse
  computation on noisy hardware. The default is the ordinary (static) uncomputation; another form is taken only when its
  estimated fidelity is strictly better (hybrid for "lookup + add" blocks). With `--route SEC` the static variant is
  compiled by `qoptct.hw`. On today's chips dynamic erasure does not pay in our models. Model estimates, not measured on
  hardware.
* Subcommands `qoptct compose` and `qoptct dyn`; `PROOF.md` is shipped in the package.
* **`qoptct.crypto.ecpoint`** (`ec_point_add`, `ec_point_add_window`; core in `qoptct/crypto/_ecp/`): exact complete addition
  of a classical point to a quantum point on a prime-field curve (also P = O, P = Q, P = -Q, P = -2Q), windowed form.
  P-256 / secp256k1 at n = 256: 3 041 421 / 3 033 954 ANDs. An exact scheme, not the cheapest known one: the approximate schemes of
  2026 are cheaper. Checked on small curves exhaustively and on edge points at n = 256; no machine proof of the whole circuit.
* **`qoptct.crypto.residue`** (accumulator compression modulo q, 106 / 138 / 148 ANDs against 682 by the formula of the
  paper; against the paper's exact circuit count 352-418 the gain is smaller; junction with the next step not checked on
  a whole circuit), **`qoptct.crypto.modadd`** (modular adders, 4n-2 ANDs), `lookup` option `fanout='gray'`
  (up to -68 % CX on large tables, same ANDs; the default is unchanged).
* **`qoptct.crypto.gf2ec`**: point addition on binary curves (m = 163: 18 160 ANDs, m = 409: 55 006); exceptional cases are
  not handled (as in BBHL / Taguchi-Takayasu). New registry entries `gf2m_square_m163_798cx`, `gf2m_square_m283_1386cx`.
  `ErasurePhaseChecker` no longer puts wires into the zero assumptions after the last phase of an erasure group.
* **`qoptct.estimate.logical`** and `qoptct estimate --attack NAME [--variant L3|L2|L0] [--source ...] [--no-reaction-limit]`,
  `--list-attacks`, `--check-attacks`: models of RSA / ECDLP / AES-Grover from published algorithms and from checked
  qoptct blocks (statuses published / qoptct / estimate); not end-to-end circuits. Kaliski count 12n^2 - 8n - 3 is used.

### Changed
* `qoptct hw`: all fidelity figures are model estimates from calibration snapshots; the portfolio keeps the Qiskit
  candidate where the qoptct one is worse.

## 0.5.0 — unreleased (release candidate)

### Added
* **`qoptct.verify.kprove`** and the subcommand **`qoptct kprove CIRCUIT.qasm REF.qc [--clean] [--budget S]
  [--engine auto|bdd|anf|sample|exact2] [--pre lt:REG:N] [--max-mb MB] [--json]`**: exact verification of large Clifford+T circuits with
  mid-circuit measurements and classically controlled corrections against a *classical* (reversible, possibly with
  diagonal phases) reference. Sum-over-paths with symbolic measurement outcomes, wire functions as reduced ordered
  BDDs with complement edges (`qoptct.verify.bdd`), so ripple carries stay linear in size; e.g. a 256-bit Gidney adder
  is proved in about a second. PROVED means for every basis input and every combination of measurement outcomes;
  REFUTED comes with a witness confirmed by the independent sparse simulator; SAMPLED / UNKNOWN are never a proof.
  Word-level specifications (`qoptct.verify.ksym`), input promises (`--pre`), garbage wires, phase checked on the
  promised domain only. Exit codes: 0 PROVED, 1 REFUTED, 2 not proved, 3 usage or input error.
  Memory budget per call: `kprove(..., max_mb=MB)` / `--max-mb MB` / `$KPROVE_MAX_MB` (default 900), counted as the
  growth of the resident set during the call (plus a floor of 200 bytes per live BDD node, for memory reused from
  earlier work), not as the peak of the whole process — a test process that already used 1 GB is not charged for it;
  over budget the BDD engine stops with UNKNOWN and reports the share of circuit ops processed (`result['mem']` gives
  the budget, the base and the peak growth). The budget limits the BDD engine only: with the default `--engine auto`
  the ANF engine (a separate worker process, not limited by `max_mb`) and the sampled check run next, so a small
  budget usually ends in SAMPLED (exit code 2); `--engine bdd` keeps the whole call within the budget.
  Engine **`exact2`** (`qoptct.verify.exact2.prove_degree2`, tried first by `auto` for a bilinear reference):
  circuits of cx / x / logical AND / AND erasure whose wires stay of GF(2) degree ≤ 2 are decided exactly by the
  1 + n + n(n−1)/2 inputs of Hamming weight ≤ 2 (a polynomial of degree ≤ 2 vanishing there vanishes everywhere; the
  phases of the AND gadget and of the erasure are checked separately), at any size: the GF(2^163) multiplier
  (489 input bits, 119 806 inputs) in 0.4 s. A mismatch with all preconditions held is a genuine counterexample.
  Limits: references must be classical; multiplications (x·y, x·c mod N, high-degree GF(2^m)) and reductions modulo
  primes with a large c usually exceed the memory budget (UNKNOWN), no induction over the width. The BDD engine on the
  GF(2^163) multiplier reaches 20–27 % of the circuit within 1.7 GB.
* **`results/crypto/`**: circuits of cryptographic building blocks with their references, each one proved by
  `qoptct kprove` (only PROVED entries are kept; see `results/crypto/index.json` and the README section), including
  the GF(2^163) multiplier `gf2m_mult_m163_908and` (strict, 908 ANDs, proved by the exact2 engine against a
  schoolbook reference built from the NIST polynomial).
* `qoptct.verify.branchsim.run_branch_scaled`: branch simulation that returns the state with its scale as log2.
* **`qoptct.crypto`** — circuit generators for cryptographic arithmetic on logical ANDs (pure functions, OpenQASM 3
  output, bit-sliced checkers, exact proofs of small instances in the tests; one shared `CheckError`):
  * `kaliski_inverse(p, n=None, reduction='generic'|'special', clean=False, legacy=False)` — Kaliski modular
    inversion: 12n²−8n−3 ANDs forward for any odd prime, 11n−2 qubits with the generic reduction (the special
    reduction needs more: secp256k1 3 107, P-256 5 513 qubits) (Litinski 2023, arXiv:2306.08585: 26n²+2n in
    the same model); the count was checked for every n tried (up to 256), not proved for all n. Special-prime
    reduction (P-256: 664 966, secp256k1: 657 503 ANDs); clean compute-copy-uncompute mode; exact proof of the emitted
    circuit at n ≤ 7 in the tests (qoptct.verify.prove and kprove). The default omits 9n+3 AND/erase pairs whose
    product is identically constant (or which only feed a constant) for every odd prime — five classes, each proved in
    the module docstring, found by an exhaustive run over all inputs at n ≤ 8; −0.3 % ANDs, −0.3 % CX at n = 256.
    `legacy=True` keeps the previous circuit (12n²+n ANDs, 11n qubits, byte-identical OpenQASM; P-256 667 273,
    secp256k1 659 810).
  * `multiply` / `square` (`qoptct.crypto.mulsq`): out-of-place integer multiplication (x, y, 0) -> (x, y, xy) and
    squaring (x, 0) -> (x, x²) on conditionally inverted additions (the idea is Häner et al., arXiv:2609.05625):
    base n²+4n−2 and (n²+7n−16)/2 ANDs (identically-0 ANDs removed; `legacy=True` keeps the published n²+4n,
    (n²+7n−6)/2 circuits), optional Karatsuba levels with garbage (n = 256, L = 3: 37 826 and 22 145 ANDs, legacy
    37 880 and 22 332; the squaring at L = 3 is −33 % ANDs against the base formula ½n(n+3)−1 = 33 151 of Häner et
    al., Sec. VII.A, at the cost of qubits and garbage; their App. A has a one-level Karatsuba squaring, Algorithm
    7, with no published count — our L = 1 is 27 049, −18.4 % against 33 151), Bennett
    version `clean=True`, independent schoolbook reference `.qc`, `prove(engine='prove'|'kprove')`.
  * `qoptct.crypto.fold`: `reduce_special` (fold for p = 2^n − c, Solinas P-256, Mersenne P-521; 2^255−19: 1 304,
    secp256k1: 2 175, P-256: 2 744 ANDs), `add_const` / `controlled_add_const`; small engine `qoptct.crypto.andcirc`.
    Not a record: special-prime arithmetic is published (arXiv:2606.02235, 2609.05625); the fair comparison is the
    whole modular multiplication (≈ −11 %, special primes only).
  * `qoptct.crypto.lookup`: narrow table lookups z ⊕= T[a] — bilinear address split (6→32: 40 ANDs vs 2^c−c−1 = 57),
    QROAM s = 1, Davio fallback, measurement-based unlookup, lookup pairs and lookup-add (acc += T[a], 79 ANDs at
    6→32); exact build-time checks, independent `.qc` references. Prior art: Luongo et al. 2025 (arXiv:2502.17325)
    for the narrow lookup in the Gidney–Ekerå framework; QROAM: Berry et al. 2019 / Low et al. 2018.
  * GF(2^m): `gf2m_multiply` (CRT bilinear multiplier, 908 / 1 343 / 1 669 / 2 501 / 3 569 ANDs for the NIST m =
    163…571; same model as Kim et al. 2022 and Garn–Kan, smaller bilinear complexities are published in
    arXiv:1303.7082), `gf2m_square` (in place, CNOT only), `gf2m_inverse` / `gf2m_divide` (Itoh–Tsujii with
    measurement-based erasure of product registers, (L+1)·M ANDs; prior art for that erasure: Khattar et al. 2025,
    arXiv:2510.10967), exact counts without building, gate-level stream check (`gf2m_itoh_gate_check`). New package
    data `qoptct/contrib/gf2m/data/crt_ranks.json`. Default register erasure `erase='terms'` (X measurement + one
    classically controlled CZ per (bilinear term, output bit): 0.67 M CZ instead of 5.4 M for the m = 163 division;
    `erase='monomials'` keeps the per-bit quadratic forms): its phases are checked exactly for any m by
    `qoptct.contrib.gf2m.qcirc.ErasurePhaseChecker` (symbolic phase polynomial of the emitted gates; all 668 520
    conditional phases of the m = 163 division; exact for the phase part only, given its classical preconditions —
    the register content and the |0> wires — which are checked by simulation on sample inputs, 16 at m = 163, so
    the whole m = 163 circuit is not proved), and whole circuits are proved by kprove against an independent
    truth-table reference for inversion m ≤ 14 and division m ≤ 12 (the tests go to 12 and 10; a removed
    conditional phase is REFUTED);
    the argument for every m is in the docstring of `qoptct.contrib.gf2m.itoh`.
* **`qoptct.estimate`** and the subcommand **`qoptct estimate --toffoli T --qubits Q [--depth D] [--shots S]
  [--model all|ge2019|gidney2025|litinski2023] [--json]`** / **`qoptct estimate --check`**: physical cost of a
  surface-code run (physical qubits, runtime, expected runtime with retries, failure probability, code distance,
  megaqubit-days) from logical counts. Formula models calibrated on Gidney–Ekerå 2019 (arXiv:1905.09749), Gidney 2025
  (arXiv:2505.15917) and Litinski 2023 (arXiv:2306.08585); not a layout tool or a compiler, and each model assumes
  the architecture of its paper. `--check` reproduces the published figures (see the module docstring for the
  tolerances and the assumptions that are ours).

### Fixed
* kprove (research version, before inclusion): after a BDD garbage collection the cached input and outcome literals
  pointed to foreign nodes, which gave false REFUTED (and "list index out of range" errors) on long runs. Literals are
  now re-created from their variable level; regression tests run with a collection at every check
  (`KPROVE_FORCE_GC`) on adders and on a fixed subset of the validation set (correct circuits and mutants of the
  benchmarks, truth from the dense simulator).
* `qoptct.verify.branchsim`: a fixed branch was simulated without renormalisation, every X-basis measurement halved
  the norm and terms below 1e-14 were dropped, so after about 50 measurements the branch vanished and
  `confirm_branch` silently reported no deviation (a real counterexample on a long circuit could not be confirmed).
  The state is now rescaled after every measurement and the scale kept as log2; `confirm_branch` compares the
  factors in log form. Test: a 120-measurement chain.

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
