# qopt — methods

This document describes what `qoptct optimize` does internally, in the order the stages run, and how every
output is checked. Code locations are relative to the `qoptct/` package directory (the package of the project qopt; `pip install qoptct`). Directories under
`qoptct/contrib/` keep the names of the development work packages (c1, c2, c8, c10, c12, eng1, eng2, tg2, tg3).

Terminology used throughout:

* **T** — number of `t`/`tdg` gates in the output text (conditional gates count too).
* **CX** — number of two-qubit gates `cx` + `cz` in the output text (conditional ones count; CZ = 1 CX).
* **extra qubits (A)** — qubits in the output minus wires declared in the input `.qc` file.
* **Cost models** (see README, "Models"):
  * *unitary* — a plain Clifford+T circuit on the same wires, equal to the input unitary up to global phase;
  * *strict (measurement-based)* — the output may use fresh ancillas starting in |0>, mid-circuit measurements and
    classically controlled Clifford corrections; **every wire of the input is data in an arbitrary state**, and
    every fresh ancilla must end in a computational basis state, unentangled with the data;
  * *clean ancillas* — as strict, but wires of the `.qc` file that are **not** listed in its `.i` line are assumed
    to start in |0> (the intended use of these benchmark circuits) and must be returned to |0>.

## 1. Unitary mode (`--mode cx`, default)

1. **Parsing** (`qopt_core.load_qc`): Feynman `.qc` → list of Toffoli-level gates (h, x, t, tdg, s, sdg, cx, ccx, z, cz, ccz).
2. **Phase folding + global cancellation** (`portfolio/evo_lib.fold_cancel`): phase-polynomial merging of T/S rotations
   acting on equal parities, followed by commutation-aware cancellation of inverse pairs.
3. **Portfolio of rewriting programs** (`portfolio/programs/c00*.py`): 13 small programs discovered by an
   LLM-driven evolutionary search (AlphaEvolve-style loop over Python programs, scored on a development set of
   generated circuits disjoint from the published benchmarks): relative-phase Toffoli pairing, CCZ target choice,
   phase-polynomial resynthesis with Gray-code parity networks, carry/majority recognition, ESOP resynthesis of
   multi-controlled gates, etc. Each program receives the Toffoli-level circuit and returns a Clifford+T circuit;
   each run is time-limited (`QOPT_PROG_SEC`, default 120 s) and followed by phase folding.
4. **Proof** (`qopt_core.prove`): the three cheapest candidates by (CX, T) are checked in order; the first proven wins:
   * n ≤ 10: full unitary comparison (`qiskit.quantum_info.Operator.equiv`);
   * otherwise PyZX `verify_equality` (time-limited), then a numerical comparison on random states in complex128
     with tolerance 1e-11 (n ≤ 26; `proto2/verify.numeric_equal`);
   * if nothing is proven, **the input circuit is returned unchanged**.
5. **Neural block resynthesis** (`--net`, optional; `bench_real.py`, `block_synth.py`): the proven circuit is cut into
   3–4-qubit blocks; each block is resynthesized by a beam search guided by graph-neural-network policy/value models
   trained by expert iteration (`models/mx_4_4q.pt`, `models/rt_3_si2_3q.pt`, `models/u_2_best_2q.pt`,
   value model `models/value_u4q.pt`). Block replacements are exact (checked as unitaries) and the whole circuit
   is proven again.
6. **CX post-pass** (`--post`, optional; `contrib/eng2/improve.py`): regions without H are re-synthesized as a phase
   polynomial + linear map; accepted only if CX decreases and the result is proven equal.

## 2. Measurement-based modes (`--mode balanced`, `--mode t`, `--meas`)

Outputs of these modes are OpenQASM 3 with `measure` and `if (c[k]) { ... }` blocks. They are written only after
the **independent** check of Section 4 (≤ 22 qubits in total; larger candidates are discarded).

### 2.1 Logical AND with measurement-based uncomputation (`contrib/c2/andc.py`)

All classical reversible parts of a circuit (after removing the H "frames" around CCZ, i.e. Toffolis) are compiled
with Gidney's temporary logical AND (arXiv:1709.06648): a product `a·b` is computed into a fresh |0> ancilla with
4 T and 3 CX, and erased later by an X-basis measurement followed by a classically controlled CZ/Z correction
(0 T). The compiler tracks the value of **every physical wire exactly, as a polynomial over GF(2)** in the input
variables:

* a product already present on some wire, or expressible as an XOR of live wires, is reused linearly (no T);
* *carrier transfer*: if the old value of a wire is needed again (compute–uncompute patterns), the new logical value
  lives in the ancilla and the wire is left untouched (restoring it is then free);
* *lazy target accumulation*: several Toffolis into the same target are summed into one polynomial whose quadratic
  part is decomposed symplectically into r products of linear forms — r ANDs instead of the number of Toffolis
  (mod5_4: 4 Toffolis → 1 AND);
* a garbage ancilla is erased by measurement as soon as its value is a quadratic form of live wires (the correction
  is a set of CZ/Z; the number of CZ is minimized over the kernel);
* if the ancilla budget is exhausted, a unitary Toffoli is used and the unitary pieces between measurements are
  phase-folded.

The search enumerates random orders of commuting logical gates × laziness levels × ancilla limits and keeps the
Pareto front in (T, CX, extra qubits). `--and-variants` sets the number of variants (default 120).

### 2.2 Clean-ancilla model (`--clean`)

Wires not in the `.i` line start with value 0. The compiler then (a) computes an AND directly into a wire whose
current value is the constant 0 (no fresh ancilla, no erasure), and (b) erases a wire itself by X-measurement when
its logical value becomes constant and its old value is quadratic in the remaining wires. For some benchmarks there
are hand-derived constructions (`contrib/c1/best_clean/*.qasm`: Karatsuba-style GF(2^m) multipliers, a 7-AND
`mod_red_21`, ripple adders with MAJ via one AND, …). They are offered as candidates **only if the input file is
byte-identical to that benchmark** (SHA-256 in `contrib/c1/best_clean/MANIFEST.json`) and are verified like
everything else.

### 2.3 CX polishing of AND circuits (`contrib/c1/run.py`, `resyn*.py`, `csyn.py`)

Each unconditioned unitary piece between H gates/measurements (a phase polynomial + linear map) is re-synthesized.
Two freedoms are exploited: a wire in |+> is invariant under CX into it (free initial row), and a CX into a wire
that is about to be X-measured equals a classically controlled Z on the control after the measurement (free final
row). Even-weight Clifford remainders are absorbed with S gates on parities that occur somewhere in the network.
A piece is replaced only if CX strictly decreases, T does not increase and the (phase polynomial, linear map) pair
is identical; the whole circuit is then re-verified. Typical gain 0–20 % CX.

### 2.4 Hybrid circuits with genuine H (`contrib/c8/hyb.py`)

For circuits that are not classical reversible (grover_5, qft_4) the H gates that break the frame analysis are
declared genuine; the circuit is cut at genuine H/T/S gates into classical blocks, each compiled by 2.1 with ancillas
returned to |0> (and reused between blocks), and a Pareto DP chooses one variant per block.

### 2.5 Multi-controlled Toffoli recognition (`contrib/c12/mcx.py`)

The exact algebraic normal form (ANF) of each output wire is computed. If the outputs have the form `y_i ⊕ P_i(x)`
(a sum/chain of multi-controlled Toffolis), candidates are synthesized (i) from a prefix tree of monomials as a ladder
of logical ANDs (T = 4(k−1), CX = 4(k−1)+1 with k−1 ancillas), (ii) from Qiskit's conditionally-clean-ancilla MCX
syntheses (Khattar & Gidney, arXiv:2407.17966; taken from `qiskit.synthesis` at run time) with the outer RCCX pair
replaced by an AND and a measurement-based erasure (−4 T, −2 CX per ancilla), and (iii) hybrids of the two.

### 2.6 Hadamard gadgets + TODD (`--mode t`, `--tgad`; `contrib/eng1`, `contrib/tg2`, `contrib/tg3`)

Internal H gates are replaced by Hadamard gadgets (fresh |+> ancilla, coupling, X-measurement, Clifford correction);
the resulting phase polynomial is T-reduced with a re-implementation of FastTODD (Vandaele, arXiv:2407.08695;
TODD: Heyfron & Campbell, arXiv:1712.01557) with deterministic classically controlled Clifford corrections.
`contrib/tg2/sel.py` gadgetizes only a subset of the H gates (windows + Pareto DP) for `--mode balanced`;
`contrib/tg3` re-synthesizes CX using the |+>/measured-wire freedoms. `--vartodd` additionally tries an external
VarTODD binary if the user provides one (`QOPT_VARTODD_BIN`); its output is used only if the signature tensor
matches. On all classical benchmarks the AND compiler (2.1) dominates the gadget route; the gadget route matters
for circuits with genuinely non-classical parts.

### 2.7 Mode selection

* `--mode cx`: unitary result only.
* `--mode balanced`: lowest T subject to CX ≤ `--cx-factor` × (unitary CX), among partial gadgetization (2.6) and
  the AND compiler (2.1–2.5); otherwise the unitary result is kept.
* `--mode t`: lowest T (ties broken by CX) among full gadgetization + TODD, partial gadgetization and the AND compiler.
* `--meas`: an older Gidney transformation (`gidney.py`, clean-ancilla model), kept for compatibility.

## 3. What is not claimed

* No optimality claims. Lower bounds are mentioned only where trivially available (multiplicative complexity
  ≥ deg − 1 per output bit).
* The AND/measurement constructions are standard (Jones 2013, Gidney 2018). The contribution is their **automatic**
  placement in arbitrary Toffoli networks via exact GF(2) tracking of wire values, combined with CX polishing,
  MCX recognition, and independent verification of every reported circuit.

## 4. Verification

Every circuit in `results/` and every output of the measurement-based modes is checked by code that does not share
the optimizer's data structures:

* `qoptct/verify_qasm.py` (≤ 22 qubits): state-vector simulation of an OpenQASM 3 subset with measurements (outcomes
  sampled) and `if` blocks on ANDs of classical bits. For 30 random input states (complex Gaussian vectors on the
  data wires, fresh ancillas in |0>) it computes `|<U_ref ψ | out>|` **without renormalization**, where `out` is the
  output state projected onto the most likely basis state of every ancilla. A value of 1 therefore certifies that
  the ancillas end in a basis state unentangled with the data and that the data part equals the reference (the
  input `.qc` lowered with standard 7-T Toffolis). The worst value over the trials is reported; acceptance threshold
  F > 1 − 1e−9. `clean=True` sets to zero the amplitudes in which a non-`.i` wire is 1. This is a randomized test,
  not a proof; a wrong circuit passes only if it acts correctly on all sampled states and branches.
* `qoptct/contrib/c10/verify_big.py` + `pathsum.py` (> 22 qubits, deterministic measurement-based circuits):
  sum-over-paths evaluation (Amy, arXiv:1805.06908) with exact variable elimination. For sampled basis inputs x and
  sampled measurement branches m it checks that the output equals c_m · U_ref|x> with c_m independent of x, and that
  live ancillas end in one basis state. Randomized over (x, m). The engine was validated against brute force on
  300 random path sums and against simulation on 80 mutants of small circuits (80/80 agree).
* Unitary circuits: full operator comparison (≤ 10 qubits), state-vector check with `verify_qasm` (≤ 22 qubits),
  PyZX `verify_equality` or a numerical check (larger).
* **Canaries**: during development every accepted point was also checked to *fail* after removing or altering one
  essential correction, CX or T gate (fidelity typically 0.0–0.9). `tests/test_release.py` keeps two such canaries.
