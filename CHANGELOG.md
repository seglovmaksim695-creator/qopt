# Changelog

Versions 0.5.0, 0.6.0, 0.7.0 and 0.8.0 are published together as one release, 0.8.0; 0.5.0 to 0.7.0 were not released separately.

## 0.8.0 — unreleased

### Added
* **`qoptct.dyadic`** and **`qoptct dyadic prove|results`**: exact proofs for circuits with rotations pi*p/2^k (QFT, Draper adders, QFT
  multipliers), a sum over paths whose phases are integers modulo 2^K (Amy, arXiv:1805.06908v2; rules Elim, HH and omega generalised from
  Z8, applied only when their condition is checked exactly; the Clifford+T prover `qoptct.verify.pathsum` is unchanged).  API: `prove_dyadic`,
  `prove_equal`, `prove_spec`, `prove_qasm`; specifications `QFTSpec`, `ArithSpec` (`adder_spec`, `const_adder_spec`, `multiplier_spec`; the
  fourier-shift lemma is reported under `assumptions`), `ClassicalSpec` / `PermutationSpec`.  Registry: Qiskit QFT, inverse QFT, QFT·QFT^-1 and
  `DraperQFTAdder` (fixed, half) up to n = 256, `RGQFTMultiplier` up to n = 32, controlled Draper adder up to n = 10,
  `HRSCumulativeMultiplier` up to n = 7; 20-22 canaries at every largest size, all REFUTED; 4000 random circuits against a dense simulator.
* **`qoptct qiskit`** (`qoptct.foreign.qiskit`): Qiskit circuits against a classical or symbolic specification or a reference circuit:
  CDKM, VBE, `HalfAdderGate`, `FullAdderGate`, `ModularAdderGate`, `IntegerComparator`, `AndGate`, `OrGate`, `BitwiseXorGate` and controlled
  adders at n = 256 for all inputs (`IntegerComparatorGate` up to n = 16); 33-34 canaries per family.  Rotations that are not multiples of
  pi/4 are proved by `qoptct.dyadic` when the call gives a dyadic specification or a reference circuit (otherwise `UNKNOWN`, `needs ==
  'fractional_phases'`).
* **`qoptct projectq`** (`qoptct.foreign.projectq`, extra `qoptct[projectq]`: `projectq==0.8.0`): ProjectQ programs against a classical
  specification (`check_projectq`) in NCT and Clifford+T modes; Binary_ECC multipliers m = 8, 16, 127, 163, 233, 283 and point addition m = 8
  (the source is not part of the package: `--source` or `$QOPT_BINARY_ECC`; the garbage wires are pinned to the original program, which
  closes a hole of the first method: one mutant at m = 163 was proved); `libs.math` and `libs.revkit` blocks; 24 canaries per size.
  `check_projectq_dyadic` proves programs with rotations through `qoptct.dyadic`; an angle is accepted only if it is exactly what ProjectQ
  stores for p*pi/2^k with k <= 38, otherwise `UNKNOWN` (QFT up to 39 qubits, `AddConstant` up to 38).
* **`qoptct auto`** (`qoptct.foreign.auto`, extra `qoptct[qualtran]`): a generic proof of Qualtran blocks from their decomposition, without
  code per block: leaves by the existing exact prover, nodes by composing the contracts of the children in Z3 with wire-flow, coverage,
  freshness and domain checks; loops need an invariant.  14 results at n = 256 (Montgomery w = 1 and 2 and the first step of ECAdd with the same
  program hashes as in 0.7; `ModAdd`, `CModAdd`, `ModSub`, `CModSub`, `ModNeg`, `CModNeg`, `ModDbl`, `CModAddK`, `Subtract`, `SubtractFrom`,
  `Negate`), 416 canaries, none proved.  `CtrlScaleModAdd` and `CModMulK` are proved at n = 8 only.
* **`qoptct loop run|verify|table|status`** (`qoptct.loop`): propose, prove, accumulate.  Proposers suggest circuits for a target; only a
  candidate that is PROVED, agrees with an independent simulation and survives at least 20 canaries enters the journal (JSONL, hash chain);
  the table of the best known circuits is built from the journal only; `verify` re-checks the chain, the file hashes, the costs and every
  proof, and fails on a missing or empty journal.  A run of 4263 proposals gave 116 proved records (`results/loop/state`).
* **`qoptct results`**: one table over all registries (Qualtran, Qiskit, ProjectQ, dyadic, auto, loop), `--long`, `--json`, `--library`.
* **Documentation**: `docs/FOREIGN.md` (overview), `docs/QISKIT.md`, `docs/PROJECTQ.md`, `docs/DYADIC.md`, `docs/AUTO.md`, `docs/LOOP.md`.

### Fixed (independent audit of 0.8)
* **Qiskit -> dyadic**: the registers of a reference circuit are matched to those of the circuit by NAME (as the Clifford+T path does;
  qubits without a unique register or an unknown `ref_clean` register give `UNKNOWN`); a gate is taken by its name only when it is an
  instance of a standard Qiskit gate class with the matching matrix, any other gate (for example a user gate called `h` whose body is S)
  is unrolled by its definition.  Before, such circuits could be reported `PROVED` although not equal.  Same rule for ProjectQ gates.
* **Float angles**: an angle above 4 pi is first reduced exactly modulo 4 pi, then compared with the tolerance (nearly always `UNKNOWN`).
* **Canaries** for the registry rows that had none (Qiskit `and`, `or`, `xor-const`, `c-cdkm-fixed`; auto `cmodmulk` and `ctrlscalemodadd`
  at n = 8); two Qiskit rows are marked in the registry instead (`comparator-gate`: none, `modular-adder-gate`: at n = 32); `qoptct loop`: records carry a hash of the target flags and `verify` checks it; `qoptct.foreign.auto`
  canaries write the candidates dropped by the filters D13 / D19 to `<result>.excluded.jsonl`.

### Findings (described in the documents of the parts)
* Qiskit 2.6.2: `VBERippleCarryAdder(kind='fixed')` is wrong for n >= 2 (an open pull request exists); `control()` of `CPhase(theta)` with
  |theta| <= pi/2^51 gives the angles (theta, theta, theta) instead of (theta/2, -theta/2, theta/2); `transpile` with the default
  `qubits_initially_zero=True` keeps only the action on |0...0>.
* ProjectQ 0.8.0: `ComparatorQuantum` is wrong for n = 1 and n >= 3, `MultiplyQuantum` for b >= 2; angles are stored with 12 decimals and
  rotations below 5e-13 become 0.
* Qualtran: the findings of 0.7 (see below) are published with this release; nothing was reported to the developers separately.

### Changed
* `qoptct/qopt_core.py`: the subcommands `qiskit`, `projectq`, `dyadic`, `auto`, `results`, `loop`.  `qoptct/verify/*` and
  `qoptct/foreign/qualtran/*` are unchanged (the proofs and the Kaliski canary evidence of 0.7 stay valid).

## 0.7.0 — published together with 0.8.0 (not released separately)

### Added
* **`qoptct qualtran`** (`qoptct.foreign.qualtran`, optional extra `qoptct[qualtran]`: `qualtran==0.7.0`, `z3-solver`): exact
  proofs about blocks of Google's Qualtran 0.7.0, gate by gate, with the unchanged qoptct provers. Subcommands `check` (a small
  bloq against its `call_classically`), `prove kaliski|montgomery|ecadd|ecwindowaddr`, `pair kaliski|montgomery`, `results`;
  exit code 0 only for PROVED (3: `PROVED_CONDITIONAL`; 4: usage error). Without the extra qoptct installs and all earlier tests pass; the new tests are skipped.
  * Strict bridge Qualtran -> qoptct (`bridge`, rules A1 controlled-And and C1 QROAMCleanAdjoint) and `check_bloq`
    (PROVED / REFUTED with a witness / UNKNOWN).
  * Kaliski inversion (`kaliski.prove_kaliski`), original Qualtran: n = 16 ... 256, three 256-bit primes; includes the
    wire-flow check of the audit (every step reads the wires the previous one wrote; damaged wirings are not proved).
  * Montgomery multiplication (`montgomery.prove_montgomery`), patched Qualtran: n = 256 windows 1 and 2, n = 64 window 4;
    includes the coverage and freshness checks of the window pieces of the audit. The original Qualtran is not proved
    (witness in `T_N >= p`).
  * Point addition (`ecadd.prove_ecadd`, `prove_ecwindowaddr`): fixed `ECAdd` at n = 256, original `ECAdd` correct exactly
    outside the groups G1-G3 of finding 05, `ECWindowAddR` with windows 1, 2, 4; the boxes of the multiplier and the inverse are
    contracts stored as data (`ecadd/contracts.py`), each tied to a registry result, and every verdict says which are only assumed.
  * Pairs "forward + `uncompute=True` = identity" (`pairs.check_pair_kaliski`, `check_pair_montgomery`) by a gate-level mirror check.
  * Registry of the proved results (`qoptct qualtran results`, `qoptct.foreign.qualtran.registry`, data in the package):
    block, Qualtran variant with the SHA-256 of the patch files, n, p, w, status, program hash, assumptions, command, source file.
  * Qualtran patches (PR 1997 sign, tables without int64, flag fix of `ECAdd`) only as files with SHA-256
    (`qoptct/foreign/qualtran/patches`), applied to a copy of Qualtran in `$QOPT_CACHE`; the installed Qualtran is never changed.
  * Documentation: `docs/QUALTRAN.md` (overview, table, assumptions), `QUALTRAN_BRIDGE.md`, `QUALTRAN_KALISKI.md`,
    `QUALTRAN_MONTGOMERY.md`, `QUALTRAN_ECADD.md`; `NOTICE` lists the modified Qualtran files (Apache-2.0).
* The Kaliski pair and the Kaliski contract at n = 16, p = 32749 are now results (`pairs_kaliski_n16_p32749`,
  `kaliski_n16_p32749`); the contract status of `ECAdd` at n = 16 no longer reports an assumption for them.

### Changed
* Kaliski lowering/macro-leaf routine, the child-process environment of the proof clients and the variant check are single
  implementations shared by the Kaliski, Montgomery, point-addition and pair clients (behaviour unchanged; pair hashes and mirror
  statistics are equal to the earlier results).
* `qoptct/qopt_core.py`: the subcommand `qualtran`. `qoptct/verify/*` is unchanged.

## 0.6.0 — published together with 0.8.0 (not released separately)

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
  estimated fidelity is strictly better. With `--route SEC` the static variant is
  compiled by `qoptct.hw`. After the real-chip check (below) the dynamic variants (measurement M, hybrid H) are **not
  selected automatically**: `DynErasePass`, `choose_circuit` and `qoptct dyn choose` return the static form unless
  `allow_dynamic=True` / `--allow-dynamic` is given (a warning is printed).
* Subcommands `qoptct compose` and `qoptct dyn`; `PROOF.md` is shipped in the package.
* **`qoptct.crypto.ecpoint`** (`ec_point_add`, `ec_point_add_window`; core in `qoptct/crypto/_ecp/`): exact complete addition
  of a classical point to a quantum point on a prime-field curve (also P = O, P = Q, P = -Q, P = -2Q), windowed form.
  P-256 / secp256k1 at n = 256: 3 031 380 / 3 023 913 ANDs in the default construction on the library blocks (Kaliski core
  12n^2 - 8n - 3, `mulsq`, `fold`; segment marks, `hints()`, `write_qasm()`); `legacy=True`: the lab construction,
  3 041 421 / 3 033 954 ANDs. An exact scheme, not the cheapest known one: the approximate schemes of 2026 are cheaper.
  Checked on small curves exhaustively and on edge points at n = 256.
* **`qoptct.verify.compose.ecadd`** (`verify_ecadd`, `wec_test`; `qoptct compose ecgen | ecverify | ecwec`; `compose.ecmath`):
  compositional proof of the point addition with a classical Q at n = 256 (secp256k1, P-256: PROVED, about 2 minutes and
  1.7 GB each): ROBDD segments, algebra sub-goals (products, fold reductions), mirror segments (uncompute halves),
  the Kaliski sub-walk; the word-level layer (chord formula, exceptional cases, safe point) by a written elementary
  argument with machine-checked polynomial identities and enumeration over all curves with p <= 31 (`PROOF.md`, part
  II). The windowed form (`ec_point_add_window`, also rebuilt on the library blocks: w = 4, 3 027 017 / 3 034 484 ANDs)
  is proved for every index k by `verify_ecadd_window` (`--window W`). Canaries (thousands of
  small mutants against exhaustive truth; full-size mutants). PROVED is valid for the default full check (chain plus the
  final strict run of the circuit); without the final run (`--no-test`, `test_samples=0`) the result is never PROVED
  (UNKNOWN, `chain_status`). Not proved: `legacy=True`.
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
  For ECDLP on P-256 / secp256k1 the point addition is the real, proved scheme (`POINTADD_PROVED`: 3 031 380 / 3 023 913
  ANDs at 6 724 / 6 916 logical qubits at L3, 3 284 271 / 3 276 804 at 4 365 / 4 557 at L0; a slow test compares the
  constants with `ec_point_add(...).counts()`), instead of the former estimate (3 041 421 - 4 x 2307 = 3 032 193). Totals
  (13 windows of 16 bits): P-256 L3 40 269 892 Toffoli (was 40 280 461), L0 43 557 475 (was 43 557 735); secp256k1 L3
  40 172 821, L0 43 460 404; physical, P-256: L3 13.0 million qubits, 0.121 megaqubit-days (was 0.122), L0 8.45 million,
  0.085. The lab numbers remain under `kaliski='legacy'`.

### Fixed (audit of the provers)
* `verify.algebra`: (1) the mirror lemma of `clean=True` now requires that the copying part writes to no wire that the
  computing part mentions in any role (before: only to the wires it changes), otherwise a false PROVED was possible;
  (2) the virtual carry of an adder (an integer in [-2, 1]) is no longer treated as a bit: its square gives UNKNOWN;
  (3) strict parsing of wire indices (no leading zeros).
* `verify.compose`: the lifter rejects negative indices `q[-j]`, indices beyond the register and classical bit numbers
  outside `bit[M]`; `verify_kaliski` checks the primality of p: a composite p gives ERROR, for p >= 2^64 outside the named
  curves PROVED carries the note "under the condition that p is prime (probabilistic test)" (key `p_primality`).
* `verify.algebra`: (4) an erasure whose condition (wire xor correction) is the constant 1 (a global sign on the outcome-1
  branch) is no longer accepted, only a proved 0 (`Engine.op_erase`, and the simulation pre-check). Before, a circuit with
  such an erasure in the compute half of a "compute - copy - uncompute" form got a false PROVED from `aprove`, and the
  chain of the point addition without its final run (`--no-test`) was fooled the same way (the mirror lemma assumes exact
  erasures). Found by an independent audit; regression tests `tests/test_verify_const1.py`. `ECChain.seg_algebra` now
  records the wires of algebra segments for the frame condition of the mirror lemma. `verify_ecadd` /
  `verify_ecadd_window` with `test_samples=0` never answer PROVED.
* An independent audit found false PROVED results on three classes of specially constructed circuits; they are closed.
  Regression tests on the audit examples: `tests/data/audit`. Some correct circuits with a virtual carry may now end in
  UNKNOWN.

### Real-chip check (ibm_kingston, 20.11.2026, one job per experiment)
Protocol and results: `docs/real_runs/2026-11-20_ibm_kingston/`.
* `qoptct hw` against Qiskit O3 on 4 small circuits: +5 % (QFT-5) and +4 % (adder) in Hellinger fidelity with fewer
  two-qubit gates (25 against 30; 27 against 43); QAOA-5 and the Heisenberg chain within noise; mean ratio 1.02.
  By the pre-set criterion: a draw, not a success. The model gains ("-14...-41 % two-qubit gates, x1.2-1.85 expected
  fidelity") are model figures; the measured gain in fidelity on small circuits was 0-5 %.
* Measurement-based erasure lost to ordinary uncomputation on every block checked; the hybrid lost too (lookadd2x2:
  0.682 against 0.772). It won only in the model (adder2: predicted 0.918, measured 0.525 for measurement erasure).
* The calibration estimate stays optimistic (adder: 0.877 estimated, 0.781 measured).

### Changed
* `qoptct dyn`: measurement-based and hybrid erasure are chosen only with an explicit option (`allow_dynamic=True`,
  `--allow-dynamic`); by default the static form is always returned. The `choose()` result gets `best_model`
  (best of all variants by the model) next to `best` (best static).
* `qoptct hw`: all fidelity figures are model estimates from calibration snapshots; the portfolio keeps the Qiskit
  candidate where the qoptct one is worse.

## 0.5.0 — published together with 0.8.0 (not released separately)

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
