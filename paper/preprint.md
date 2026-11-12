---
title: "Automatic logical-AND compilation with measurement-based uncomputation and independent verification for T- and CNOT-count reduction of Clifford+T circuits"
author: "Maksym Shchehlov"
date: "Draft, 2026-11-12 (not submitted)"
---

> **Status: draft.** Numbers are taken from `results/` of qopt 0.2.0 and were re-verified for this draft.
> References marked [unverified] were not checked against the original source.

## Abstract

Fault-tolerant cost estimates are dominated by the number of T gates, while the number of two-qubit (CNOT)
gates drives routing and latency. For circuits given as Toffoli networks — arithmetic, oracles, reversible
benchmarks — the cheapest known way to implement a Toffoli whose target is uncomputed later is Gidney's
temporary logical AND (4 T) with measurement-based uncomputation (0 T). Existing toolchains apply this
construction only where a programmer writes it explicitly. We describe qopt, an open-source optimizer that places
logical ANDs automatically in arbitrary Toffoli networks by tracking the value of every wire exactly as a polynomial
over GF(2), reuses products, erases garbage by X-measurement when it is a quadratic form of live wires, polishes
the CNOT count of the result, and recognizes multi-controlled Toffolis. Every output is re-checked by checkers that
are independent of the compiler: a state-vector simulator that rejects ancillas left entangled with the data, and
a sum-over-paths checker for larger circuits. On 23 standard benchmarks we report verified circuits in three
explicitly separated cost models — unitary, strict measurement-based (all wires are data), and clean ancillas —
together with a comparison against nine quantum software tools run on the same inputs. Independent verification
also exposed a wrong 51-T circuit produced by an earlier gadget pipeline of ours and a wrong circuit produced by a
third-party tool; both are reported with reproductions.

## 1. Introduction

T-count is the standard proxy for the cost of fault-tolerant Clifford+T circuits because each T gate consumes a
distilled magic state. Unitary T-count optimization has been studied extensively: phase folding and matroid
partitioning [AMM14], ZX-calculus rewriting [KvdW20], TODD and its successors [HC18, Van25], and reinforcement or
tensor-decomposition methods [Ruiz24]. Unitary CNOT optimization is addressed by rewriting systems and search
[Nam18, Xu24, Li23, Riu23]. When mid-circuit measurements and classical feed-forward are available — as they are
in any fault-tolerant architecture — a Toffoli whose result is later uncomputed costs only 4 T [Jon13, Gid18],
and multi-controlled Toffolis with clean or conditionally clean ancillas become much cheaper [KG25, Kole26a, Kole26b].

These constructions are available as library primitives (e.g. `And`/`And†` in Qualtran [Har24], `TemporaryAND` in
PennyLane [Ber18], `AND` in Q# [Bev22]), but a user must know where to apply them. Given a Toffoli network "as is",
none of the nine tools we tested inserts them automatically (Section 6). Moreover, measurement-based circuits are
harder to check than unitary ones: an ancilla left entangled with the data, or a wrong classically controlled
correction, is invisible to checks that renormalize the output or test only basis states.

**Contributions.**
1. A compiler that places logical ANDs and measurement-based uncomputation automatically in arbitrary Toffoli
   networks, based on exact GF(2) polynomial tracking of all wires (Section 4.1), with CNOT polishing (4.3), hybrid
   handling of genuine H gates (4.4) and multi-controlled-Toffoli recognition (4.5).
2. A clear separation of three cost models (Section 3) and verified results in each (Section 5).
3. Independent verification of every reported circuit (Section 4.6), which found real errors (Section 7).
4. A reproducible package: code, circuits, tables and verification scripts (Apache-2.0).

## 2. Background

*Clifford+T and counts.* We count T as the number of `t`/`tdg` gates and CX as the number of `cx` + `cz` gates in the
output text, including classically controlled ones. A Toffoli in the standard decomposition costs 7 T and 6 CX.

*Logical AND [Gid18].* For a fresh |0> ancilla `a`, the AND computes `a = x·y` with 4 T and 3 CX (up to a Clifford
phase). Uncomputation: measure `a` in the X basis; if the outcome is 1, apply CZ(x, y). No T gates are needed.

*Benchmarks.* We use the 24-circuit suite of the Feynman toolkit [Amy19; AMM14]; hwb6 is excluded from the
release because its file forbids redistribution, leaving 23 circuits with 5–30 qubits.

## 3. Cost models

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

The models are nested: a circuit valid in the strict model is valid in the clean model. Numbers from different
models must not be compared directly. Published H-gadget T-count records [HC18, Van25, Ruiz24] are yet another
model (gadget ancillas, no feed-forward on the output).

## 4. Methods

### 4.1 Logical-AND compiler with exact GF(2) tracking

After removing the H "frames" that turn CCZ into Toffolis, a benchmark circuit is a classical reversible network.
The compiler simulates it symbolically: each physical wire holds a polynomial over GF(2) in the input variables
(bitset monomials). For each logical Toffoli `t ^= x·y` it chooses between: (i) linear reuse — the product, or the
new value, is already an XOR of live wires; (ii) a new AND into a fresh ancilla (4 T, 3 CX); (iii) carrier transfer —
if the old value of `t` is needed later, the new value lives in the ancilla and `t` is untouched. Toffolis into the
same target can be accumulated lazily: the quadratic part of the accumulated polynomial is decomposed symplectically
into r products of linear forms, giving r ANDs (mod5_4: four Toffolis become one AND). A garbage ancilla is erased
by X-measurement as soon as its value is a quadratic form q of live wires; the correction is the CZ/Z
pattern of q, with the number of CZ minimized over the kernel. If the ancilla budget is exhausted, a unitary Toffoli
is emitted and unitary pieces are phase-folded. A randomized search over gate orders, laziness and ancilla limits
yields a Pareto front in (T, CX, extra qubits).

In the clean model the compiler additionally knows that non-input wires start at 0: it computes an AND directly
into such a wire, and erases a wire itself by measurement once its logical value is constant.

### 4.2 Unitary optimizer

The unitary mode uses a portfolio of 13 small rewriting programs (relative-phase Toffolis [Mas16], phase-polynomial
resynthesis with Gray-code parity networks [AAM18], carry/majority recognition, ESOP resynthesis, …) found by an
LLM-driven evolutionary program search in the spirit of AlphaEvolve [Nov25], followed by phase folding and
cancellation, optional neural block resynthesis (graph-neural-network policy trained by expert iteration on random
2–4-qubit blocks) and a proof of equality (full unitary ≤ 10 qubits, otherwise PyZX [KvdW19] or a numerical check).
If equality is not proven, the input is returned.

### 4.3 CNOT polishing

Each unconditioned piece between H gates and measurements is a phase polynomial with a linear map. It is
re-synthesized using two freedoms specific to measurement-based circuits: a wire in |+> is invariant under CX into
it, and a CX into a wire that will be X-measured equals a classically controlled Z after the measurement. Even-weight
Clifford remainders are absorbed with S gates on parities that occur in the network. A replacement is accepted
only if CX strictly decreases, T does not increase and the algebraic data are identical.

### 4.4 Genuine Hadamards

For circuits that are not classical reversible (grover_5, qft_4), H gates that break the frame analysis are kept;
the circuit is cut into classical blocks compiled by 4.1 (ancillas returned to |0> and reused), and a dynamic program
chooses one Pareto point per block.

### 4.5 Multi-controlled Toffoli recognition

From the exact algebraic normal form of the outputs, sums/chains of multi-controlled Toffolis are recognized and
synthesized as AND ladders (T = 4(k−1)), as Khattar–Gidney syntheses [KG25] with the outer RCCX pair replaced by
an AND and measurement-based erasure, or as hybrids. This closes the gap to hand-written MCX constructions on
barenco_tof_n.

### 4.6 Independent verification

*verify_qasm* (≤ 22 qubits) simulates the OpenQASM 3 output with sampled measurement outcomes on random input states
(30; 10 for 20–22 qubits) and computes |⟨U_ref ψ|out⟩| **without renormalization** after projecting each ancilla onto
its most likely basis state. A value of 1 certifies that ancillas end unentangled in a basis state and that the data
part equals the reference unitary (the input with 7-T Toffolis). The clean variant zeroes amplitudes with a non-input
wire set. *verify_big* (> 22 qubits) evaluates sum-over-paths expressions [Amy18] with exact variable elimination for
sampled basis inputs and measurement branches, requiring the output to equal c_m·U_ref|x> with c_m independent of x.
It was validated against brute force (300 random path sums) and against simulation on 80 mutants (80/80 agree).
During development every accepted circuit was also checked to fail after removing or altering one correction, CX or
T gate ("canaries").

## 5. Results

Table 1 lists the best verified point per benchmark and model (files in `results/best/`). Input counts use 7 T and 6 CX per Toffoli/CCZ.

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

**Observations.**

* *Unitary model.* The unitary optimizer is our baseline; per-circuit numbers are competitive with published
  unitary CNOT results on this suite but this paper does not claim new unitary records (an earlier internal comparison
  with GUOQ, PyZX and Qiskit is part of the repository history, not of this release).
* *Strict model.* 19 of the 23 circuits have verified strict-model results (state-vector check up to 22 qubits; the
  minimum-T points of qcla_com_7 and qft_4 are large H-gadget circuits with 61 and 47 qubits, checked by
  sum-over-paths). T drops to 4 per necessary AND:
  tof_n and barenco_tof_n reach 4(k−1) for k controls (tof_3 8 T, tof_10 59 T in the strict model where the
  benchmark's ancilla wires are data), mod5_4 4 T (one AND), mod_red_21 44 T vs 73 unitary, grover_5 124 T vs 166.
  At the minimum-T point CX is often *also* lower than the unitary result (tof_3 10 vs 12, barenco_tof_4 13 vs 24,
  grover_5 155 vs 160, mod_red_21 58 vs 61, rc_adder_6 54 vs 65) but not always (vbe_adder_3 33 vs 28,
  gf2^4_mult 115 vs 86 at its 0-ancilla minimum-T point); the Pareto fronts in `results/TABLE.md` show the trade-off
  with the number of extra qubits.
* *Clean model.* Knowing that the benchmark's ancilla wires start in |0> gives further gains where the circuits
  compute into those wires: tof_10 36/37 (+1) vs 59/71 (+3) strict, mod_mult_55 12/24 (+1), and with hand-derived
  constructions mod_red_21 28/48 with no extra qubit, gf2^4_mult 36/79 and gf2^5_mult 56/124 (Karatsuba-style
  bilinear multiplication). The hand constructions are labelled as such and are not produced by the general compiler.
* *Not covered.* adder_8, csum_mux_9 and qcla_mod_7 (24–30 qubits) have no measurement-based result because the
  optimizer only emits circuits it can check with the state-vector checker (≤ 22 qubits); ham15-low has a clean-model
  result only.
* *Lower bounds.* For tof_n, barenco_tof_3, mod5_4 and vbe_adder_3 the strict-model AND count equals deg−1 of the
  highest-degree output (a lower bound on ANDs per output), so these T counts cannot be improved within the AND
  construction; for mod_red_21 the bound is 6 ANDs vs 11 used.


## 6. Comparison with other tools

We ran Qiskit 2.6.1, PyZX 0.10.8, pytket 2.18.5, Cirq 1.7, the Q# Resource Estimator (qsharp 1.31), Qualtran 0.7.0,
pyLIQTR 1.0.1, Qrisp 0.9.10 and PennyLane 0.45.1 on the same `.qc` files, each in its own environment, in two ways:
*as is* (the gate list of the file) and, for tof_n and barenco_tof_n only, *by meaning* (we rewrote the circuit by
hand as a multi-controlled X or a ladder of the tool's own logical-AND primitives). Every circuit a tool produced was
verified with the same checker in the strict and in the clean model; resource estimates without a circuit (Q#,
Qualtran, cirq_ft) are reported as estimates only. The full table is `results/COMPARISON.md`.

Findings:

1. *As is*, no tool replaces Toffolis by logical ANDs with measurement-based uncomputation. Unitary tools give
   about 7 T per Toffoli minus phase folding (best: PyZX and tket's ZX passes, with T equal to our unitary mode and
   1.5–3× more CX). Q# and Qualtran estimate 4 T per Toffoli but output no circuit for these inputs. Our strict-model
   circuits are at or below that estimate in T on every circuit, with a concrete verified circuit and CX count.
2. *By meaning*, tools match or beat a general-purpose compiler: PennyLane's `TemporaryAND` ladder (3 CX per AND)
   gives exactly our clean-model tof_n numbers (8/9, 12/13, 16/17, 36/37 with one ancilla); for barenco_tof_n with
   dirty ancillas, Qiskit's Khattar–Gidney syntheses dominated our earlier points with one ancilla (barenco_tof_4
   23/18 vs 27/29) — our MCX recognition (Section 4.5) now gives 19/16. The advantage of qopt is that it needs no
   knowledge of what the circuit computes.
3. Recent measurement-based MCX constructions [Kole26a, Kole26b] report, e.g., C³X with 11 T / 10 CX (worst-case
   branch) and, with one clean ancilla, C⁴X 19 T / 17 CX; our tof_3 strict point is 8 T / 10 CX (+2 ancillas) and our
   barenco_tof_4 point with one ancilla is 19 T / 16 CX. These numbers are taken from the papers, not re-run.


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

## 7. Errors found by independent verification

1. **Our own 51-T hwb6 circuit was wrong.** An earlier pipeline (H gadgets + TODD + deterministic corrections + CX
   re-synthesis) reported hwb6 at 51 T / 535 CX on 33 qubits — equal to the best published H-gadget T-count. It had
   passed the pipeline's internal checks. The sum-over-paths checker rejects it in every sampled branch: besides
   outcome-dependent X flips (repairable for free) an input-dependent phase of degree up to 5 remains, i.e. some gadget
   corrections are wrong. Other circuits of the same pipeline (mod_red_21, qcla_com_7, rc_adder_6, qft_4) pass.
2. **Qrisp 0.9.10, `mcx(method="jones")` with `compile_mcm=True`** produces a circuit with an extra phase on the
   controls; a three-qubit test in Qrisp's own simulator gives uniformly random outcomes where a deterministic result
   is expected (with `compile_mcm=False` the result is correct).
3. **Pitfalls:** pytket's ZX extraction leaves an implicit qubit permutation that a gate-by-gate export loses;
   Qualtran costs an imported `cirq.CCZ` at 7 T but `cirq.CCX` at 4 T; Qrisp assumes MCX work qubits are |0>.

Reproductions are in `docs/FINDINGS.md`. The lesson is that measurement-based circuits need an independent checker
that does not renormalize and that inspects ancillas.

## 8. Limitations

* Verification is randomized (inputs and measurement branches are sampled); only unitary results ≤ 10 qubits are
  proven exactly.
* Measurement-based outputs are limited to 22 qubits in total by the state-vector checker; larger candidates are
  discarded by the optimizer (verify_big covers some deterministic gadget circuits up to 61 qubits).
* Several numbers come from longer development searches than one default run; hand-derived clean-model constructions
  (e.g. Karatsuba-style GF(2^m) multipliers, a 7-AND mod_red_21) are used only for the exact benchmark files and are
  marked as such.
* We make no optimality claims. For tof_n, barenco_tof_n, mod5_4 and vbe_adder_3 the strict-model T equals 4 × (a
  lower bound deg − 1 on the number of ANDs per output), but for hwb6-type permutations the gap is unknown.
* The benchmark suite is small and dominated by classical reversible circuits, where the AND construction applies
  everywhere; circuits with many genuine non-Clifford rotations are not addressed.
* The neural component contributes little (≈ 1–2 % CX in unitary mode).

## 9. Conclusion

Automatic placement of logical ANDs with measurement-based uncomputation, driven by exact GF(2) tracking of wire
values, gives large T reductions on Toffoli networks at equal or lower CNOT counts than unitary optimizers, without
hand-written MCX knowledge. Independent, non-renormalizing verification is essential: it caught errors both in our
pipeline and in a third-party tool. Code, circuits and checkers are released for reproduction.

## References

All references below were checked on arXiv (abstract pages) on 2026-11-12 unless marked [unverified].

* [AMM14] M. Amy, D. Maslov, M. Mosca. Polynomial-time T-depth optimization of Clifford+T circuits via matroid partitioning. IEEE TCAD (2014), doi:10.1109/TCAD.2014.2341953. arXiv:1303.2042.
* [AMMR13] M. Amy, D. Maslov, M. Mosca, M. Roetteler. A meet-in-the-middle algorithm for fast synthesis of depth-optimal quantum circuits. IEEE TCAD 32(6):818–830 (2013), doi:10.1109/TCAD.2013.2244643. arXiv:1206.0758.
* [AAM18] M. Amy, P. Azimzadeh, M. Mosca. On the CNOT-complexity of CNOT-PHASE circuits. Quantum Sci. Technol. (2018), doi:10.1088/2058-9565/aad8ca. arXiv:1712.01859.
* [Amy18] M. Amy. Towards large-scale functional verification of universal quantum circuits. EPTCS 287, doi:10.4204/EPTCS.287.1. arXiv:1805.06908.
* [Amy19] M. Amy. Formal Methods in Quantum Circuit Design. PhD thesis, University of Waterloo, 2019 (referenced by the Feynman repository README; thesis itself not checked) [unverified].
* [Bar95] A. Barenco et al. Elementary gates for quantum computation. Phys. Rev. A 52, 3457 (1995), doi:10.1103/PhysRevA.52.3457. arXiv:quant-ph/9503016.
* [Ber18] V. Bergholm et al. PennyLane: Automatic differentiation of hybrid quantum-classical computations. arXiv:1811.04968.
* [Bev22] M. E. Beverland et al. Assessing requirements to scale to practical quantum advantage. arXiv:2211.07629.
* [Gid18] C. Gidney. Halving the cost of quantum addition. Quantum 2, 74 (2018), doi:10.22331/q-2018-06-18-74. arXiv:1709.06648.
* [HC18] L. E. Heyfron, E. T. Campbell. An efficient quantum compiler that reduces T count. arXiv:1712.01557.
* [Har24] M. P. Harrigan et al. Expressing and analyzing quantum algorithms with Qualtran. arXiv:2409.04643.
* [Jav24] A. Javadi-Abhari et al. Quantum computing with Qiskit. arXiv:2405.08810.
* [Jon13] C. Jones. Novel constructions for the fault-tolerant Toffoli gate. Phys. Rev. A 87, 022328 (2013), doi:10.1103/PhysRevA.87.022328. arXiv:1212.5069.
* [KG25] T. Khattar, C. Gidney. Rise of conditionally clean ancillae for efficient quantum circuit constructions. Quantum 9, 1752 (2025). arXiv:2407.17966.
* [Kole26a] A. Kole, T. Schnittka, R. Drechsler. Measurement-driven adaptive low-overhead implementation of multi-controlled Toffoli gates. arXiv:2605.18159 (2026).
* [Kole26b] A. Kole, M. Assaad, T. Schnittka, R. Drechsler. Adaptive Clifford+T decomposition of large Toffoli gates with one clean ancilla. arXiv:2605.18169 (2026).
* [KvdW19] A. Kissinger, J. van de Wetering. PyZX: Large scale automated diagrammatic reasoning. EPTCS 318 (2020), doi:10.4204/EPTCS.318.14. arXiv:1904.04735.
* [KvdW20] A. Kissinger, J. van de Wetering. Reducing T-count with the ZX-calculus. Phys. Rev. A 102, 022406 (2020), doi:10.1103/PhysRevA.102.022406. arXiv:1903.10477.
* [Li23] Z. Li et al. Quarl: A learning-based quantum circuit optimizer. arXiv:2307.10120.
* [Mas16] D. Maslov. On the advantages of using relative phase Toffolis with an application to multiple control Toffoli optimization. Phys. Rev. A 93, 022311 (2016), doi:10.1103/PhysRevA.93.022311. arXiv:1508.03273.
* [MSC19] G. Meuli, M. Soeken, E. Campbell, M. Roetteler, G. De Micheli. The role of multiplicative complexity in compiling low T-count oracle circuits. arXiv:1908.01609.
* [Nam18] Y. Nam, N. J. Ross, Y. Su, A. M. Childs, D. Maslov. Automated optimization of large quantum circuits with continuous parameters. npj Quantum Inf. (2018), doi:10.1038/s41534-018-0072-4. arXiv:1710.07345.
* [Nov25] A. Novikov et al. AlphaEvolve: A coding agent for scientific and algorithmic discovery. arXiv:2506.13131.
* [Riu23] J. Riu, J. Nogué, G. Vilaplana, A. Garcia-Saez, M. P. Estarellas. Reinforcement learning based quantum circuit optimization via ZX-calculus. arXiv:2312.11597.
* [Ruiz24] F. J. R. Ruiz et al. Quantum circuit optimization with AlphaTensor. arXiv:2402.14396.
* [Sae10] M. Saeedi, M. Saheb Zamani, M. Sedighi, Z. Sasanian. Reversible circuit synthesis using a cycle-based approach. arXiv:1004.4320.
* [Sei24] R. Seidel et al. Qrisp: A framework for compilable high-level programming of gate-based quantum computers. arXiv:2406.14792.
* [Siv20] S. Sivarajah et al. t|ket⟩: A retargetable compiler for NISQ devices. arXiv:2003.10611.
* [Van25] V. Vandaele. Lower T-count with faster algorithms. Quantum (2025), doi:10.22331/q-2025-09-16-1860. arXiv:2407.08695.
* [Xu24] A. Xu, A. Molavi, S. Tannu, A. Albarghouthi. Optimizing quantum circuits, fast and slow (GUOQ). arXiv:2411.04104.
* Cirq (Google Quantum AI), pyLIQTR (MIT Lincoln Laboratory): software, cited by version only [unverified: no paper reference checked].
