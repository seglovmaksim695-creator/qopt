# Findings: incorrect circuits found during this work

Independent verification caught several incorrect circuits: one produced by our own earlier pipeline, others
by, or through pitfalls of, third-party tools. They matter for how published numbers should be compared.
Reproduction scripts are in `docs/findings/` (re-run on 2026-11-12 with the versions below).

Tool versions: Qiskit 2.6.1, PyZX 0.10.8, pytket 2.18.5, Cirq 1.7, Qualtran 0.7.0, Qrisp 0.9.10, PennyLane 0.45.1,
qsharp 1.31.0, pyLIQTR 1.0.1. Each third-party tool was installed in its own virtual environment; none of them
is a dependency of qopt.

## F1. Our own "hwb6: 51 T" measurement-based circuit is wrong

* **What.** During development a full-gadgetization pipeline (Hadamard gadgets + TODD + deterministic Clifford
  corrections, followed by CX re-synthesis) reported `hwb6` at **51 T / 535 CX** (33 qubits) — equal to the best
  published T-count for hwb6 in the H-gadget model. The circuit was too large for the state-vector checker and had
  been accepted on the pipeline's internal checks only.
* **How it was found.** The sum-over-paths checker (`qopt/contrib/c10/verify_big.py`, METHOD.md §4) rejects it in
  every sampled measurement branch, and also a 51 T / 587 CX variant. The same checker confirms the other large
  circuits of that pipeline (mod_red_21 51/358, qcla_com_7 59/565, rc_adder_6 37/311, qft_4 53/265).
* **Diagnosis.** Two error types: (a) X flips of output bits that are affine in the measurement outcomes and
  independent of the input (fixable for free with conditional X gates); (b) after fixing (a), an input-dependent
  phase (−1)^{f_m(x)} with deg f_m up to 5 remains, i.e. some gadget Z-type corrections are wrong. This cannot be
  repaired with Clifford corrections. The checker itself was validated against brute force (300 random path sums)
  and on 80 mutants of small circuits (80/80 verdicts agree with simulation).
* **Consequence.** No 51-T hwb6 circuit is claimed. Pipelines that output measurement-based circuits need an
  independent check; internal consistency checks of the same pipeline were not sufficient.
* **Reproduce.** hwb6 and circuits derived from it are not redistributed (bench_qc/SOURCES.md). With a local copy:
  `sh scripts/fetch_benchmarks.sh bench_extra && python docs/findings/f1_hwb6_51T.py <circuit.qasm> bench_extra/hwb6.qc`
  → `{'ok': False, 'F': 0.0}`.

## F2. Qrisp 0.9.10: `mcx(method="jones")` with `compile(compile_mcm=True)` produces a wrong circuit

* **What.** The Jones Toffoli (4 T, one ancilla, mid-circuit measurement) compiled with mid-circuit-measurement
  support leaves an extra phase on the two controls.
* **Minimal reproduction** (`docs/findings/f2_qrisp_jones.py`, Qrisp's own simulator): prepare |+,+,0>, apply
  `mcx(method="jones")`, then an ordinary Toffoli (which undoes it), then H on both controls. Expected: the all-zero
  outcome with certainty. Observed (2000 shots):

  | method | compile_mcm | outcomes |
  |---|---|---|
  | jones | False | `000`: 2000 |
  | jones | **True** | 8 outcomes, each ≈ 230–275 shots (uniform) |
  | gray | False / True | `000`: 2000 |

* **Effect on benchmarks.** Our strict check gives F = 0.15–0.45 on 18 of 19 benchmark circuits compiled with this
  method. On a circuit whose wires all start in |0> (grover_5) the extra phase happens to be invisible, so that point
  passes a clean-ancilla check by accident; we treat all Qrisp-jones numbers as invalid.

## F3. Qrisp: MCX work qubits are assumed to be |0>

Qrisp places MCX work qubits on wires it has not touched yet and assumes they are |0>. On benchmarks whose
"ancilla" wires carry data (barenco_tof_n) the methods `auto`, `maslov`, `balauca`, `khattar` and the gidney ladder
fail the strict check; in the clean-ancilla model (tof_n) they pass. `khattar` with 4 controls was also observed to
fail in the clean model (not investigated further, no minimal reproduction). This is a modelling assumption rather
than a bug, but it must be accounted for when comparing numbers.

## F4. pytket 2.18.5: `ZXGraphlikeOptimisation` leaves an implicit wire permutation (pitfall)

The optimized circuit is correct only together with its implicit qubit permutation. A gate-by-gate export (e.g. to
QASM) without `replace_implicit_wire_swaps()` is a different unitary. `docs/findings/f4_tket_implicit_swaps.py` on a
3-Toffoli circuit prints:

```
implicit permutation after ZX pass: {q[0]: q[1], q[1]: q[4], q[2]: q[3], q[3]: q[2], q[4]: q[0]}
gate-by-gate export equals the input: False
after replace_implicit_wire_swaps():   True
```

This is documented tket behaviour, not a bug, but our first export made exactly this mistake (F ≈ 0.05). In the
comparison we always call `replace_implicit_wire_swaps()` (each SWAP then costs 3 CX). After the ZX pass,
`FullPeepholeOptimise` often produces rotation angles that are not multiples of π/4 (outside Clifford+T); such
variants were discarded.

## F5. Other observations (not bugs)

* Qiskit: MCX syntheses based on rotations (b95, noaux_v24/hp24, gray code) become 1,900–10,000 T after exact
  Clifford+T approximation; on qft_4, `transpile(optimization_level=3)` merges single-qubit gates into arbitrary
  rotations and the re-synthesis inflates T to 409.
* Q#: `Controlled X` decomposes into a unitary V-chain (2n−3 CCZ) without measurement-based uncomputation; the
  n−1 CCZ ladder requires writing `AND`/`Adjoint AND` explicitly.
* Q# Resource Estimator and Qualtran return resource *estimates* for these inputs, not circuits; their numbers are
  listed as estimates, never as verified circuits.
* Qualtran 0.7.0: the cost of an imported Cirq circuit depends on how the Toffoli is written
  (`docs/findings/f5_qualtran_ccz.py`): `cirq.CCX`/`cirq.TOFFOLI` become Qualtran's `Toffoli` (1 Toffoli = 4 T in
  its cost model), while `cirq.CCZ` is decomposed by Cirq into 7 T + 8 Cliffords. The benchmark files write
  Toffolis as H·CCZ·H, so an "as-is" import is costed at 7 T per Toffoli; the comparison table therefore lists
  both numbers ("as-is" and "Toffoli bloq, 4 T").

  ```
  cirq.CCX      -> toffoli: 1               total_t_count = 4
  cirq.TOFFOLI  -> toffoli: 1               total_t_count = 4
  cirq.CCZ      -> t: 7, clifford: 8        total_t_count = 7
  ```
