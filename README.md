# qopt — a verified Clifford+T circuit optimizer

[![PyPI](https://img.shields.io/pypi/v/qoptct.svg)](https://pypi.org/project/qoptct/)
[![tests](https://github.com/seglovmaksim695-creator/qopt/actions/workflows/tests.yml/badge.svg)](https://github.com/seglovmaksim695-creator/qopt/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24182059.svg)](https://doi.org/10.5281/zenodo.24182059)

> **Package name:** the project is called *qopt*, but on PyPI the name `qopt` belongs to an unrelated project.
> Since version 0.3.0 the package is therefore installed and imported as **`qoptct`**:
> `pip install qoptct`, `import qoptct`, command `qoptct`.

`qopt` reduces the **T count** and the **CNOT (CX) count** of Clifford+T circuits given as Toffoli networks
(Feynman `.qc` format). It has a unitary mode and measurement-based modes, and it **never reports a circuit that
has not been checked**: unitary outputs are proven equal to the input (or the input is returned), and every
measurement-based output is proven equal by the exact prover `qoptct.verify` or, when a proof is out of reach,
re-checked by independent simulators over many sampled inputs and measurement outcomes before it is written
(the JSON report says which).

Repository: <https://github.com/seglovmaksim695-creator/qopt>  · Version 0.6.0 (unreleased) · License Apache-2.0 (benchmark files: see `bench_qc/SOURCES.md`) · Author: Maksym Shchehlov

## Installation

```bash
pip install qoptct              # Python 3.11+, macOS or Linux; provides the `qoptct` command and the Qiskit plugin
pip install "qoptct[net]"       # optional: + PyTorch for the neural mode --net
pip install "qoptct[qcec]"      # optional: + MQT QCEC as an extra equivalence checker for unitary circuits
```

Dependencies: numpy, qiskit ≥ 2.0, pyzx. PyTorch is needed only for `--net`; without it every other mode works
and `--net` stops with a message that says what to install.

From the repository (to run the tests, the benchmarks and the scripts):

```bash
git clone https://github.com/seglovmaksim695-creator/qopt.git && cd qopt
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt                      # exact versions used for the results (incl. torch)
pip install -e ".[test]"
pytest -q tests                                      # a few minutes
```

`requirements-lock.txt` pins the full environment used for testing. Caches are written to `$QOPT_CACHE`
(default `~/.cache/qopt`).

## Usage

```bash
qoptct optimize bench_qc/tof_3.qc -o tof_3_opt.qc                          # unitary, minimum CX (default --mode cx)
qoptct optimize bench_qc/tof_3.qc -o tof_3_opt.qc --net --post             # + neural block resynthesis + CX post-pass
qoptct optimize bench_qc/tof_3.qc -o tof_3.qasm --mode balanced            # fewer T, CX <= 1.5x unitary CX
qoptct optimize bench_qc/mod_red_21.qc -o m.qasm --mode t --clean          # minimum T, clean-ancilla model
qoptct optimize X.qc --report report.json                                  # JSON report (also printed to stdout)
qoptct optimize X.qc -o x.qasm --mode t --prove                            # + exact proof of the written result
qoptct prove x.qasm X.qc                                                   # exact proof of any output against its input
```

Python: `import qoptct; names, ref, best, how = qoptct.optimize("X.qc")`.

| option | meaning |
|---|---|
| `--mode cx` (default) | unitary output, minimum CX; proven equal to the input, otherwise the input is returned. Written as `.qc`, or as OpenQASM 3 if the `-o` file name ends in `.qasm` |
| `--mode balanced` | measurement-based output (OpenQASM 3) with lower T and CX ≤ `--cx-factor` (1.5) × unitary CX, if one is found and verified |
| `--mode t` | lowest T found (measurement-based), CX may grow |
| `--clean` | clean-ancilla model (see below) for the measurement-based modes |
| `--net` | neural block resynthesis (models in `qoptct/models/`; needs `pip install "qoptct[net]"`) |
| `--post` | CX post-pass on H-free regions (slow) |
| `--meas` | older Gidney transformation (clean-ancilla model), kept for compatibility |
| `--tgad`, `--vartodd` | Hadamard gadgets + TODD; `--vartodd` also tries an external VarTODD binary given by `QOPT_VARTODD_BIN` |
| `--prove` | prove the written result exactly with `qoptct.verify` (also for unitary results); the status is in the report |
| `--no-arith`, `--no-bandc`, `--no-gf2m`, `--no-cxpass` | switch off one optional module (see below) |
| `--no-contrib` | switch off all optional modules (0.3 behaviour) |
| `--budget S` | total seconds for all optional-module stages (default 1500; 0 = unlimited) |
| `--time-limit S` | overall seconds of a run (default 3600; 0 = none): all searches — the unitary portfolio, its proofs, the module stages and the measurement-based searches — stop there and the best result checked so far is written (in the worst case the input itself). Only the final check of a written measurement-based result (or `--prove`) may follow: usually under a minute, at most ~15 min |
| `--and-sec S` | seconds for the variant search of the logical-AND compiler in balanced/t modes (default 1800; 0 = unlimited) |
| `--contrib-sec S`, `--source-sec NAME=S` | time limit of one module call (default 240 s; `bandc` gets more on large circuits) |
| `--require-proof` | accept module results above 22 qubits only with an exact proof |
| `--keep-temp` | keep the work directory with the intermediate files of the stages (path in the report: `temp_dir`) |

Measurement-based outputs are written only if they pass the checks described below; otherwise the unitary result is
written. The JSON report (`final`) says which candidate was used, how it was checked and its T/CX/qubit counts
(`final.file`: `qc` or `openqasm3`; `search_seconds`: time until the searches ended). The written OpenQASM 3 is
readable by `qiskit.qasm3.loads` (with `qiskit-qasm3-import`): conditions are single measured bits or their negation,
possibly nested. Where a stage produced an XOR (`^`) or AND (`&&`) condition, the result is rewritten — the condition
is computed into an idle measured qubit and measured once more (same T and CX, one more measurement), or the
correction is repeated per term — and proven again; `final.qiskit_rewrite` gives the counts before and after. If no
rewritten form can be proven, the previous readable result is written. Intermediate files are kept in a private
temporary directory (removed at the end unless `--keep-temp`), and `-o` is written atomically (temporary file +
rename), so it never holds a partial or intermediate circuit.

Exit codes of `qoptct optimize`: 0 — a result was written; 3 — usage or input error (missing or unreadable file,
malformed circuit, unknown gate, undeclared wire, missing output directory), reported as one `error: ...` line on
stderr. `qoptct prove`: 0 PROVED, 1 REFUTED, 2 UNKNOWN, 3 usage or input error (an input error is never reported as
REFUTED). `python -m qoptct.verify` uses the same codes.

## Optional modules and exact proofs (new in 0.4)

All modules are part of the package; each one can be switched off, and none of them can break a run: every call
runs in a child process with a hard time limit, an error or a timeout only skips the module, and **every circuit a
module returns is checked again by the core** before it is compared with the other candidates.

| module | what it does | used in |
|---|---|---|
| `qoptct.contrib.arith` | recognizes classical arithmetic (adders, comparators, addition mod 2^n−1, hidden weighted bit) by its truth table and rebuilds it on logical ANDs | `--mode balanced`, `--mode t` |
| `qoptct.contrib.bandc` | logical-AND compiler for large classical circuits (hundreds to 1600 wires) with exact wire-value tracking | `--mode balanced`, `--mode t` |
| `qoptct.contrib.gf2m` | GF(2^m) multipliers from low-rank bilinear decompositions (tables shipped with the package) | `--mode balanced`, `--mode t` |
| `qoptct.contrib.cxpass` | CX post-pass: block re-synthesis of CX/phase regions, merged and XOR-conditioned corrections | all modes |
| `qoptct.verify` | exact prover and fast checkers (below) | checks of all measurement-based results |

How a candidate is accepted:

* ≤ 22 qubits in total: the state-vector checker (30 random inputs), then the exact proof; if the proof does not
  finish in time, a stronger sampled check with 64–256 random inputs and freshly sampled measurement outcomes.
* > 22 qubits: the exact proof; if it ends UNKNOWN (time or width), the sparse-state / path-sum sampled checks of
  `qoptct.verify` (not a proof; the report says so). `--require-proof` accepts only proven results.
* A REFUTED result is never written; the optimizer falls back to the proven unitary result.

`qoptct prove CIRCUIT REF.qc [--clean] [--budget S]` proves that an output (OpenQASM 3 as written by qoptct, or a
`.qc` file) implements the input in the unitary, strict or clean model. It prints JSON and exits with 0 (PROVED),
1 (REFUTED), 2 (UNKNOWN) or 3 (usage or input error). The proof covers every basis input and every combination of measurement outcomes
(sum-over-paths normal form or a trace miter over Z[ω]); unitary circuits can also be handed to MQT QCEC when
`qoptct[qcec]` is installed. From Python:

```python
from qoptct.verify import prove, check_dense, check_sparse, quick_check
r = prove("x.qasm", "X.qc", clean=False, budget_sec=120)    # {"status": "PROVED"|"REFUTED"|"UNKNOWN", "method", "detail", ...}
```

The prover runs in its own process (`python -m qoptct.verify._worker`), so it can be called from scripts without an
`if __name__ == "__main__":` guard and from the Qiskit plugin. `QoptPass(cxpass=True)` also applies the CX post-pass
to every optimized run (a post-pass result is used only if it is proven equal again).

`qoptct kprove CIRCUIT.qasm REF.qc [--clean] [--budget S] [--engine auto|bdd|anf|sample|exact2] [--max-mb MB] [--json]`
(new in 0.5) is an exact prover for large
circuits with a *classical* reference (arithmetic, table lookups, S-boxes): the same sum-over-paths with symbolic
measurement outcomes, but wire functions are BDDs, so adders of thousands of bits are proved in seconds. Exit codes
as above (2 also for SAMPLED, a sampled check that is not a proof). The BDD engine stops with UNKNOWN when the memory
it uses during the call (growth of the resident set) exceeds `--max-mb` / `kprove(max_mb=...)` / `$KPROVE_MAX_MB`
(default 900 MB); multiplications and reductions modulo general primes usually do. In the default mode
`--engine auto` the run then continues with the ANF engine, which runs in a separate process and is not limited by
this budget, and with a sampled check, so the result is then usually SAMPLED (exit code 2). `--engine bdd` uses only
the BDD engine and keeps the whole call within the budget. For circuits whose wires stay of GF(2)
degree ≤ 2 (GF(2^m) multipliers) the engine `exact2` decides exactly at any size from the inputs of Hamming weight
≤ 2; `auto` tries it first when the reference is bilinear.
From Python: `from qoptct.verify.kprove import kprove; kprove("c.qasm", "ref.qc", clean=True, max_mb=1500)`.

`qoptct aprove CIRCUIT.qasm SPEC.json [--budget S] [--json]` (new in 0.6; the same engine is
`qoptct kprove CIRCUIT.qasm SPEC.json --engine algebra`) is a word-level algebraic prover for the integer arithmetic
that BDDs cannot handle: multipliers, squarers, reductions modulo special primes, and their compositions, built from
logical ANDs with measurement-based erasure. The specification is a JSON file with input words, output words and
integer claims, for example
`{"inputs": {"x": {"range": [0, 256]}, "y": {"range": [256, 512]}}, "outputs": {"out": {"range": [512, 1024]}},
"claims": ["out == x*y"], "garbage": []}`.
PROVED means: for every basis input and every combination of measurement outcomes the circuit computes a classical
function for which the claims hold, the input wires are unchanged, all other wires end in |0> (garbage wires: any
classical function), and no branch has an input-dependent phase. `multiply` / `square` at n = 256 with Karatsuba
levels 0..3 take 4-9 s and at most 0.8 GB; the fold reductions take under a second; a full multiplication modulo
Curve25519 / secp256k1 / P-256 takes 9-18 s (`python scripts/algebra_proofs.py proofs`). The prover is not a general
equivalence checker: circuits that are not carry-chain arithmetic end in UNKNOWN (exit code 2) or, when they use
gates outside its strict language, in an input error (exit code 3); a few correct Karatsuba circuits on very narrow
or odd-width leaves are also left UNKNOWN. At n = 256 the following is proved: multiplication and squaring (Karatsuba levels 0..3), the fold reductions for
Curve25519, secp256k1, P-256, P-384, P-521 and whole modular multiplications for three of them. Left UNKNOWN: squaring
with Karatsuba level >= 2 on odd or very narrow widths, `square(256, 3, legacy=True)`, and every circuit that is not
carry-chain arithmetic (controlled swaps, general Boolean logic; e.g. the Kaliski inversion, see `qoptct compose`).
From Python:
`from qoptct.verify.algebra import aprove, spec_mulsq; c = multiply(256, 3); aprove(c.to_qasm(), spec_mulsq(c))`.

`qoptct compose` (new in 0.6; `python -m qoptct.verify.compose`) proves the Kaliski inversion `kaliski_inverse` at
n = 256 (secp256k1, P-256, a random 256-bit prime; generic and special reduction; forward and clean; Curve25519, P-384 and
P-521 in the clean mode with special reduction) on the real OpenQASM text (250-520 MB), 40-250 s and about 0.5 GB per circuit: PROVED for every
x in [1, p) and every combination of measurement outcomes. The proof has layers, and they are not equally strong:
the 2 061-2 302 segments of the circuit are checked against their contracts on BDDs, and the gadget lemmas are exact
arithmetic, all by machine; the lemmas about the integer algorithm (bounds, parity, the invariant u*s + v*r = p) are
checked by Z3 for all n and p (`pip install "qoptct[z3]"`, needs the `z3` binary); the gluing of segments and of the Z3 formulation to the
contracts is a written argument (`qoptct/verify/compose/PROOF.md`), not machine-checked. Only the default
construction (`legacy=False`) is covered.

Limits of these proofs (after an independent audit of the provers): garbage registers (the upper half of a product, the
quotient, Karatsuba garbage) are not specified, such products are "with garbage"; the `clean` variant was not proved by
`aprove` at n = 256 (only up to n = 64; the clean Kaliski inversion at n = 256 is proved by `qoptct compose`). The Kaliski
statement is for x in [1, p) and prime p: x = 0 and x >= p are not covered; the primality of the named curves is a known
fact, for any other p < 2^64 it is checked deterministically and for larger p by a probabilistic test, and the result is
marked (`p_primality`); a composite p gives ERROR. The provers give no certificate for an independent checker. The audit
found false PROVED results on three classes of specially constructed circuits (a flaw in the mirror lemma, a virtual
carry treated as a bit, negative wire indices in the lifter); they are fixed and the audit examples are regression tests
(`tests/data/audit`). Some correct circuits with a virtual carry may now end in UNKNOWN (conservative).

## Qiskit plugin

`pip install qoptct` registers two stage plugins for the Qiskit transpiler, both named `qoptct`:

* `optimization_method="qoptct"` — Qiskit's own optimization stage of the chosen level, then `QoptPass` on the
  Clifford+T parts of the (already translated and routed) circuit;
* `init_method="qoptct"` (optional, recommended for Toffoli circuits) — `QoptPass` on the circuit as written,
  while the Toffolis are still Toffolis, before Qiskit's own init stage. Most of the savings come from here.

```python
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator

qc = QuantumCircuit(5)                                                 # tof_3 from bench_qc/
qc.ccx(0, 1, 4); qc.ccx(2, 4, 3); qc.ccx(0, 1, 4)

basis = ["h", "x", "z", "s", "sdg", "t", "tdg", "cx"]                 # Clifford+T
plain = transpile(qc, basis_gates=basis, optimization_level=3)
ours = transpile(qc, basis_gates=basis, optimization_level=3,
                 init_method="qoptct", optimization_method="qoptct")
print(plain.count_ops(), ours.count_ops())
assert Operator(ours).equiv(Operator(qc))                              # same operator up to a global phase
```

CX/T after `transpile` with a Clifford+T basis (Qiskit 2.6.2, levels 1–3 gave the same numbers):

| circuit | Qiskit level 3 | `optimization_method="qoptct"` | `+ init_method="qoptct"` |
|---|---|---|---|
| tof_3 | 18 / 19 | 16 / 15 | 12 / 15 |
| mod5_4 | 28 / 22 | 25 / 16 | 11 / 8 |
| barenco_tof_3 | 24 / 24 | 20 / 16 | 16 / 16 |

A part of the circuit is replaced only if the new version is **proven** equal, strictly cheaper (CX, then T),
and uses only gates and qubit pairs allowed by the basis / coupling map / target; everything else
(measurements, rotations, other gates) is left as it is. With a coupling map, compare with
`Operator.from_circuit(out)`. On circuits with up to 10 qubits per part the global phase is kept exactly.
The "strictly cheaper" rule holds for each replaced part; the later stages of `transpile` (for example the
synthesis of rotations) can still change the totals, so on circuits with rotations compare the final counts
with and without the plugin.

The pass can also be used directly:

```python
from qiskit.transpiler import PassManager
from qoptct.qiskit_qopt import QoptPass
p = QoptPass(); out = PassManager([p]).run(qc); print(p.last)   # runs, replaced, CX/T before and after
```

`QoptMeasPass(ancillas=[...])` is the measurement-based variant (logical AND, needs qubits known to start in |0>).

## Compiling for noisy hardware: `qoptct hw` (new in 0.6, experimental)

`qoptct hw` compiles a unitary circuit for an IBM-like device (heavy-hex coupling, native CZ or ECR + RZ, SX, X) and
minimizes what matters on NISQ hardware: the number of two-qubit gates after routing, the depth, and the expected
fidelity computed from the calibration data of the backend.

```
pip install "qoptct[hw]"                      # qiskit-ibm-runtime (fake backends = calibration snapshots), qiskit-aer
qoptct hw --backend fake_torino in.qasm -o out.qasm --report out.json
```
```python
from qoptct.hw import compile_hw
r = compile_hw(qc, "fake_torino", seeds=range(5), seconds=20)
r.circuit, r.init_map, r.final_map, r.metrics["twoq"], r.metrics["fid"], r.strategy, r.verified
```

It is a portfolio; the candidate with the best expected fidelity that passes the equivalence check is returned:
Qiskit `optimization_level=3` (one run per seed, always included, so the result is never worse than Qiskit by the
metric used), qoptct's proven CX/T optimization of the logical circuit followed by Qiskit, and the qoptct router
(`qoptct/hw/router.py`) with and without the logical optimization. The router is SABRE-like with four additions:
commutation-aware dependencies (all gates that may legally come next are in the front layer), exact SWAP cost
(a SWAP next to a gate on the same pair is merged with it; the cost is taken from the KAK decomposition of the
merged block: +0, +1 or +3 native gates), "floating" blocks (a block stays mergeable while only gates that commute
with it pass over its wires) and a mirror rule (right after a gate, swap its pair if that helps), plus line layouts on
low-error paths and a final move of the used region to the best-calibrated isomorphic region of the device.
Every returned circuit is checked against the input with the qubit maps taken into account (state-vector test on
random product states, up to 26 used qubits). Nothing is ever submitted to IBM; `scripts/hw_real_run.py` prepares one
short real run and refuses to submit without `--submit` and `QOPT_ALLOW_REAL=1`.

The portfolio keeps the Qiskit candidate wherever the qoptct one is worse by the metric, so the result is never worse
than Qiskit level 3 by that metric. All fidelity numbers are model estimates from the chip model and a calibration
snapshot (fake backends); only the small real-chip check below was measured. Check against real runs: on 8 real IBM jobs (ibm_fez, ibm_marrakesh; 286 result rows, one calibration snapshot per chip) the calibration estimate understated the error by a factor of about 1.6 (1.6-1.7 for the models compared). As a tool for a preliminary choice between variants it is useful: it ordered 93.7 % of the distinguishable pairs correctly (100 % on Fez, 76.4 % on Marrakesh, where only 6 tasks were available). The absolute fidelity numbers must not be trusted. A fitted correction factor exists but did not pass the accuracy criterion and is not built into the code.

Results on 76 algorithmic circuits (4-30 qubits), methods, the comparison with Qiskit level 3 and tket, the noisy
simulation and the limitations are in [`docs/HW.md`](docs/HW.md). Those gains are model figures (fake backends); on a
real chip the measured gain in fidelity was 0-5 % on small circuits:

### Check on a real chip (ibm_kingston, 20.11.2026, one job per experiment)

Two experiments, one Sampler job each, success criteria written down before submission
(`docs/real_runs/2026-11-20_ibm_kingston/`). Statistical error: about 0.01 (2000 shots) and 0.025 (400 shots).

* **`qoptct hw` against Qiskit O3 (best of 5 seeds), 4 small circuits.** qoptct is better by +5 % (QFT-5) and +4 %
  (adder) in Hellinger fidelity, with fewer two-qubit gates (25 against 30; 27 against 43); on QAOA-5 and the Heisenberg
  chain the difference is within noise. Mean ratio 1.02. By the pre-set criterion (success needed at least 3 of 4
  circuits and a mean ratio of at least 1.05) the result is a **draw**, not a success.
* **Measurement-based erasure on the chip loses to ordinary uncomputation on every block checked** (adder2, adder3,
  cmp2); the hybrid lost as well (lookadd2x2, X test: 0.682 against 0.772). The model badly underestimates the cost of a
  measurement in the middle of a circuit (adder2: predicted 0.918, measured 0.525). The hybrid won only in the model.
  Therefore dynamic variants are no longer chosen by the model: only by the explicit option `--allow-dynamic`
  (`allow_dynamic=True`), which prints a warning.
* **The calibration estimate stays optimistic** (adder after `qoptct hw`: estimate 0.877, measured 0.781). Use it to rank
  variants, not as an absolute number.

## Measurement-based erasure on hardware: `qoptct dyn` (new in 0.6, experimental)

On IBM chips the default is the ordinary (static) uncomputation of a logical AND (reverse gadget, 3 CX); measurement-based
erasure with a classical correction (M) and the hybrid (H: single ANDs reversed, multi-term erasures measured) are **not
selected automatically**: `DynErasePass`, `choose_circuit` and `qoptct dyn choose` always return the static form unless
`--allow-dynamic` (`allow_dynamic=True`) is given, which prints a warning. `qoptct.dyn` compiles the variants of an
arithmetic block for a backend calibration and lists their model estimates. Reason: on a real chip (ibm_kingston,
20.11.2026) the model strongly underestimated the cost of a measurement in the middle of a circuit (adder2: predicted
0.918, measured 0.525), and the hybrid lost to ordinary uncomputation (0.682 against 0.772); it won only in the model
(a measurement with feed-forward idles every live qubit for about 2.8 µs). `qoptct dyn choose adder 4 --backend
FakeBoston`; with `--route SEC` the static variant is compiled by `qoptct.hw`. Extra: `pip install "qoptct[dyn]"`.
Model figures come from chip models and calibration snapshots; the real-chip check is described under `qoptct hw`.


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

## Cryptographic building blocks (new in 0.5)

`results/crypto/` holds circuits of building blocks of quantum attacks on AES and elliptic-curve cryptography, each
with its reference (`.qc`) and a `.json` record (counts, model, wire layout, verdict). **Every circuit there is
proved by `qoptct kprove`** for all inputs and all measurement outcomes; only PROVED circuits are included.
References of the AES and GF(2^m) circuits were checked against independent definitions (AES tables and the
MixColumns test vector of FIPS-197; GF(2^m) arithmetic); linear references were generated from the definition,
independently of the circuit. Model: logical AND = 4 T + 3 CX, uncomputation by X measurement = 0 T (Gidney 2018);
"CX" below counts `cx` and `cz`, including classically controlled ones.

| circuit | what | AND (T) | CX | qubits | compared with |
|---|---|---|---|---|---|
| `aes_sbox_clean_330cx` | \|x,0⟩ → \|x,S(x)⟩ | 32 (128) | 330 | 26 | same AND count as Meuli et al. 2019 (arXiv:1908.01609: 32 AND at 48 qubits); the difference here is 26 qubits and the CX count. T-depth 98 — designs optimized for depth use 34 AND at T-depth 4–6 |
| `aes_sbox_depth4_33and` | \|x,0⟩ → \|x,S(x)⟩, AND-depth 4 | 33 (132) | 552 | 44 | published AND-depth-4 designs: 34 AND |
| `aes_sbox_strict_371cx` | y ⊕= S(x), any y | 32 (128) | 371 | 26 | |
| `aes_isbox_clean_339cx`, `aes_isbox_strict_365cx` | S⁻¹, clean / any y | 32 (128) | 339 / 365 | 26 | |
| `aes_sbox_inplace_64and`, `aes_isbox_inplace_64and` | \|x,0⟩ → \|0,S(x)⟩ (and S⁻¹) | 64 (256) | 704 / 710 | 26 | Toffoli-model in-place results are not directly comparable |
| `aes_mixcolumns_90cx` | one MixColumns column in place | 0 | 90 | 32 | published in-place CNOT circuits: 92 (Xiang et al. 2020), 91. Depth 34 is worse than theirs (28–30), and the output bits end on permuted wires (`out_wires` in the .json); circuits that may use new wires (g-XOR metric) reach fewer XORs |
| `lookup_c6_w32_40and`, `lookup_c6_w21_33and` | table lookup z ⊕= T[a], 6-bit address, 32 / 21-bit word | 40 / 33 | 1477 / 991 | 47 / 40 | 2^c − c − 1 = 57 AND for the standard lookup. The idea of a narrow lookup combined with an adder in the Gidney–Ekerå framework is due to Luongo et al. 2025 (arXiv:2502.17325); with lookups shared between several primes (Sharma–Lau 2026) this construction gives no gain |
| `gf2m_mult_m8_clean` | c = a·b in GF(2^8) | 26 (104) | 347 | 25 | small analogue of the m = 163…571 multipliers |
| `gf2m_mult_m163_908and` | c ⊕= a·b in GF(2^163) (NIST polynomial), any c | 908 (3632) | 75 823 | 490 | same model as Kim et al. 2022 (992 Toffoli) and Garn–Kan (999); smaller bilinear complexity published (906, arXiv:1303.7082). Proved by the `exact2` engine (degree ≤ 2) |
| `gf2m_square_m{163,233,283,409,571}_*` | in place a ← a² mod f (NIST polynomials) | 0 | 875 / 825 / 1534 / 1418 / 3239 | m | compared with the generic in-place bound m² − m used in the literature, not with a published circuit. The earlier lab files are kept as proved; the generator `gf2m_square(m)` builds 798 / 825 / 1386 / 1418 / 3239 CX, and the circuits for m = 163 (`..._798cx`) and m = 283 (`..._1386cx`) are in the registry as well |

What is **not** claimed: no record for whole attacks (RSA, ECDLP); this release contains proved building blocks, not
end-to-end attack circuits. The measurement-based erasure of GF(2^m) products (used for division circuits) was
published before by Khattar et al. 2025 (arXiv:2510.10967).

### Generators: `qoptct.crypto` (new in 0.5)

Pure functions (same arguments, same circuit) that build the arithmetic of these attacks in the logical-AND model,
with OpenQASM 3 output, a bit-sliced checker and, for small sizes, an exact proof in the tests. Counts are logical
ANDs (= Toffoli, 4 T each).

| function | what | count | compared with |
|---|---|---|---|
| `kaliski_inverse(p, n=None, reduction='generic'\|'special', clean=False, legacy=False)` | Kaliski modular inversion, result −x⁻¹·2^{2n} mod p | 12n²−8n−3 AND forward (generic; `legacy=True`: 12n²+n), 11n−2 qubits (generic); P-256 special 664 966 AND / 5 513 qubits, secp256k1 657 503 / 3 107 | Litinski 2023 (arXiv:2306.08585): 26n²+2n in the same model. The count was checked by counting for every n tried up to n = 256; it is not proved for all n. Proved exactly at n = 256 by `qoptct compose` (segments and lemmas by machine, the gluing by a written argument, `qoptct/verify/compose/PROOF.md`) |
| `multiply(n, karatsuba_levels=0, legacy=False)`, `square(n, ...)` | integer x·y and x² (no reduction) on conditionally inverted additions | base n²+4n−2 and (n²+7n−16)/2; n = 256, 3 Karatsuba levels: 37 826 / 22 145 (published construction `legacy=True`: 37 880 / 22 332) | squaring with 3 Karatsuba levels is −33.2 % ANDs against the base formula of Häner et al. (arXiv:2609.05625, Sec. VII.A: ½n(n+3)−1 = 33 151), at the cost of qubits (3 000 vs 1 275 at L = 0) and garbage (clearing it by the inverse pass: 44 277 ANDs in total). Häner et al. also give a one-level Karatsuba squaring (App. A, Algorithm 7) without a published count; our one level (L = 1, 27 049) is −18.4 % against 33 151. The conditionally inverted additions and the one-level Karatsuba squaring on them are their ideas; the multi-level recursion with garbage is ours. The base multiplier is not better than theirs |
| `fold.reduce_special(p, method='auto')` | reduction of a 2n-bit product modulo p = 2^n − c, NIST P-256 (Solinas), Mersenne P-521 | 2^255−19: 1 304, secp256k1: 2 175, P-256: 2 744 | not a record: special-prime arithmetic is already in arXiv:2606.02235 and 2609.05625. The fair comparison is the whole modular multiplication, ≈ −11 % only for special primes |
| `lookup.narrow_lookup(table, c, w)`, `lookup.lookup_add`, `lookup.qroam1_lookup` | table lookup z ⊕= T[a] with a bilinear address split; QROAM s = 1; measurement-based unlookup | 6→32: 40 AND (2^c−c−1 = 57) | the narrow lookup inside lookup-add in the Gidney–Ekerå framework is due to Luongo et al. 2025 (arXiv:2502.17325); QROAM is Berry et al. 2019 / Low et al. 2018. The saving is per lookup, not end-to-end |
| `gf2m_multiply(m)`, `gf2m_square`, `gf2m_inverse`, `gf2m_divide` | GF(2^m) for the NIST binary fields: CRT bilinear multiplier, in-place squaring (CNOT only), Itoh–Tsujii inversion / division | m = 163…571: 908 / 1 343 / 1 669 / 2 501 / 3 569 AND per product | same model as Kim et al. 2022 and Garn–Kan; smaller bilinear complexities are published (arXiv:1303.7082: 906 / 1 340 / 1 668 / 2 495 / 3 566). Measurement-based erasure of product registers: Khattar et al. 2025 (arXiv:2510.10967) |

```python
from qoptct.crypto import kaliski_inverse, multiply, square, gf2m_mult_count, default_modulus
from qoptct.crypto import fold, lookup
import random

kaliski_inverse(13).counts()['AND']                    # 157 = 12·4² − 8·4 − 3 (legacy=True: 196)
multiply(256, karatsuba_levels=3).counts()['AND']      # 37826 (+ garbage, uncomputed after use; legacy=True: 37880)
fold.reduce_special(fold.P25519).n_and                 # 1304
table = [random.getrandbits(32) for _ in range(64)]
lookup.narrow_lookup(table, 6, 32).n_and               # 40 for most random tables (sometimes 41)
gf2m_mult_count(163, default_modulus(163))['AND']      # 908
c = square(8); qasm, ref = c.to_qasm(), c.reference_qc()   # circuit and independent reference for `qoptct kprove`
```

### Physical cost: `qoptct estimate` (new in 0.5)

`qoptct.estimate` turns logical counts (Toffoli/AND, logical qubits, optionally reaction depth) into physical
qubits, runtime and failure probability of a surface-code run. These are **formula models calibrated on the papers**
— `ge2019` (Gidney–Ekerå 2019, arXiv:1905.09749), `gidney2025` (Gidney 2025, arXiv:2505.15917), `litinski2023`
(Litinski 2023, arXiv:2306.08585, active volume and baseline) — **not a layout tool or a compiler**. Each model
assumes the structure of its paper (e.g. `gidney2025`: few hot qubits, the rest in cold storage), so results are
comparable only within one model. `qoptct estimate --check` reproduces the published figures (within 0–4 % for the
exact GE layout, Gidney 2025 and Litinski; the generic GE model differs more, see the module docstring).

```
qoptct estimate --toffoli 6.5e9 --qubits 1399 --shots 9.2 --model gidney2025   # ≈ 0.89 M physical qubits, ≈ 5 days
qoptct estimate --check
```
```python
from qoptct.estimate import estimate
r = estimate("litinski2023", toffoli=1.09e8, qubits=3000)
r.physical_qubits, r.runtime_s, r.failure_prob, r.code_distance
```

### Exact point addition, residues, modular adders, binary curves (new in 0.6)

Building blocks added to `qoptct.crypto`; used as modules (`from qoptct.crypto import ecpoint, residue, modadd, gf2ec`).
All counts are logical ANDs (= Toffoli) in the same model as above.

* `ecpoint.ec_point_add(p, a, b, Q)`, `ec_point_add_window(...)`: **exact** complete addition of a classical point Q to a
  quantum point P on a prime-field curve (also for P = O, P = Q, P = -Q, P = -2Q); the windowed form takes Q from a table.
  P-256 / secp256k1, n = 256: 3 031 380 / 3 023 913 ANDs (6 724 / 6 916 qubits) in the default construction on the library
  blocks (Kaliski core, `mulsq`, `fold`); `legacy=True` keeps the lab construction (3 041 421 / 3 033 954 ANDs,
  6 742 / 6 934 qubits). This is an exact scheme, not the cheapest known one: the approximate schemes of 2026
  (arXiv:2609.05625, arXiv:2606.02235) are cheaper and use fewer qubits (they have a non-zero failure probability).
  Checked exhaustively on small curves, on edge points at n = 256, phases by sparse simulation.
  **The full exact point addition is proved at n = 256** (P-256 and secp256k1; classical Q and the window w = 4) by
  `qoptct compose ecverify` (about 2 minutes and 1.5-2.1 GB per run). What exactly is proved:
  * PROVED is valid for the default full check: the chain of segment proofs plus the final strict run of the circuit on
    concrete points. A run without that final run (`--no-test`, `test_samples=0`) is never reported as PROVED (status
    UNKNOWN, field `chain_status`): an independent audit showed that the chain alone could be fooled (see CHANGELOG, Fixed).
  * By machine: the circuit repeats the word-level program (every segment of the OpenQASM text: ROBDD segments, algebra
    sub-goals for the products and reductions, syntactic mirror segments for the uncompute halves; the Kaliski lemmas by Z3).
    That the word-level program equals P + Q is a written theorem (`PROOF.md`, part II) checked by enumeration over all
    curves with p <= 31 and by runs at n = 256, not a machine proof.
  * Domain: points of the curve and O; p prime; odd group order; order of Q above 4; Karatsuba L = 3; window w = 4
    (every index k). `legacy=True` is not proved.
  * Trusted in addition: the code of the provers (walker, ROBDD, algebra prover, Z3), not checked by an independent tool.
* `residue.fold_accumulator(q, clen, nbits)`: compression of the accumulator modulo q (loop 2 of Gidney 2025):
  106 / 138 / 148 ANDs for q = 1 574 500 against 682 by the formula of the paper. Against the exact count of the
  paper's circuit (352-418) the gain is smaller; the junction with the next step (loop 3) was argued in words and not
  checked on a whole circuit. Proved by `qoptct kprove` for 29 bits / window 4 and 32 bits / window 4 at l = 21.
* `lookup.davio_lookup(..., fanout='gray')`: optional Gray-code fan-out of CX; up to -68 % CX on large tables
  (c = 10, w = 1024: 167 231 instead of 528 201 CX, same ANDs). It is an option, the default is unchanged.
* `modadd`: modular adders / subtractors, doubling, comparisons. Modular addition 4n-2 ANDs (published: 4n), a gain of 2 ANDs.
* `gf2ec.gf2_point_add(m, ...)`: point addition on a binary curve over GF(2^m), second point classical:
  m = 163: 18 160 ANDs (10 registers) / 29 056 (4 registers); m = 409: 55 006. Exceptional cases (point at infinity,
  P1 = P2, P1 = -P2, P1 + P2 = +-P2) are **not** handled, as in Banegas et al. 2021 and Taguchi-Takayasu 2023.
  With the default `lowering='crt'` the 4-register variant for m = 163 needs 1 142 qubits; 981 qubits need
  `lowering='parity'` (about 6 times more CX).

### Models of whole attacks: `qoptct estimate --attack` (new in 0.6)

`qoptct.estimate.logical` reproduces the logical cost (ANDs, logical qubits) of published RSA, ECDLP and AES-Grover
algorithms and gives the same frameworks with the blocks of `qoptct.crypto`. It is a resource-estimation tool: **the
numbers are models built from checked blocks, not end-to-end circuits.** Every line has a status: `published`
(reproduces the number of the paper), `qoptct` (blocks checked, whole circuit not built), `estimate`.

```
qoptct estimate --list-attacks
qoptct estimate --attack ecdlp-p256 --variant L0
qoptct estimate --attack rsa2048-g25
qoptct estimate --check-attacks        # compares with the published tables
```

## Hardware check (IBM Quantum)

On `ibm_marrakesh`, the unitary `qopt` circuits returned the correct answer more often than the input circuits on all
five benchmarks tried (+2 to +10 percentage points), while the measurement-based circuits were worse on current
hardware (mid-circuit measurement latency). Details: `results/hardware/README.md`.

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
classical bits (eight files from 0.3.0 — best and Pareto points of qft_4 and qcla_com_7 — use `&&` in conditions, which
`qiskit.qasm3.loads` does not read; `qoptct prove` and `qoptct.verify` read them); data qubits are `q[0..n-1]` in the order of the `.v` line of the `.qc` file, extra qubits follow
and start in |0>.

## Limitations

* Not every result is proven: when the exact prover ends UNKNOWN (deep classical logic such as hwb, very wide
  circuits), measurement-based results rely on sampled checks (many random inputs and measurement outcomes); the
  report's `final.check` / `final.proof` say which. Use `--require-proof` (modules) or `qoptct prove` to insist on a
  proof. Unitary results ≤ 10 qubits are proven exactly; larger ones rely on PyZX, a numerical check (≤ 26 qubits)
  or the exact path-sum prover of `qoptct.verify`.
* Above 22 qubits in total only the module results (arith, bandc, gf2m, cxpass) are available in the
  measurement-based modes; the older H-gadget and AND stages need the state-vector check and are skipped there.
* The neural models were trained on random circuits of 2–4 qubits; `--net` gives small gains (≈ 1–2 % CX).
* Runtime: unitary mode takes seconds for small circuits and up to an hour for the largest benchmarks; the
  measurement-based modes take minutes (with the optional modules up to `--budget` seconds more). Several rewriting programs have internal time limits, so results can
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
