# qoptct.dyadic — exact proofs for circuits with rotations π·p/2^k

`qoptct.verify.pathsum` (the Clifford+T prover of 0.4–0.7) keeps phases modulo 8 (ω = e^{iπ/4}).  Circuits with
rotations π/2^k — the QFT, Draper's adders, QFT multipliers in Qiskit and ProjectQ — need finer phases.
`qoptct.dyadic` is a separate prover (the old one is not changed) whose phases are integers modulo 2^K, with K
chosen per problem.

```python
from qoptct.dyadic import DCircuit, Pi, prove_dyadic, QFTSpec, adder_spec, multiplier_spec, ClassicalSpec
from qoptct.dyadic.qiskit_io import from_qiskit          # optional, needs qiskit

c = from_qiskit(adder_qft_d00(64, kind='fixed'))                # or DCircuit(n); c.add('cp', (i, j), Pi(1, 8))
r = prove_dyadic(c, adder_spec(a=range(64), b=range(64, 128)))
r['verdict']        # 'PROVED' | 'REFUTED' | 'UNKNOWN' | 'ERROR' | 'INTERNAL'
r['global_phase']   # exact, e.g. 'pi*1/4' (circuit = e^{i g} * spec)
r['assumptions']    # ['fourier-shift lemma: ...'] for arithmetic specifications (Lemma 2 below)
r['counterexample'] # for REFUTED: a basis input where the circuit differs
```

| function | claim checked |
|---|---|
| `prove_dyadic(circuit, target, **kw)` | dispatches to the two below (target: a DCircuit or a specification) |
| `prove_equal(circuit, reference, phase='global', ancillas=None)` | `circuit |x,0> = e^{ig} (reference ⊗ I) |x,0>` for all basis x; ancillas (default: the extra wires) start in 0, the extra ones must end in 0 |
| `prove_spec(circuit, spec, phase='global', ancillas=None)` | the same with a specification: `QFTSpec`, `ArithSpec` (`adder_spec`, `const_adder_spec`, `multiplier_spec`), `ClassicalSpec` / `PermutationSpec` |
| `prove_qasm(qasm, reference, clean=False)` | OpenQASM 2/3 text or path against a DCircuit, OpenQASM or a qoptct `.qc` file (clean model as in `qoptct prove`) |

None of them raises: bad input gives `ERROR`, an internal inconsistency `INTERNAL`, anything not decided `UNKNOWN`;
only `PROVED` means proved.  `phase='exact'` also demands g = 0.  Options: `budget_sec` (default 120),
`max_terms` (400 000), `selfcheck` (dense re-check up to 10 qubits, default on).

Supported gates (Qiskit names and conventions, controls first): `id x y z h s sdg t tdg sx sxdg p u1 phase rz rx ry u
u3 u2 gphase cx cy cz cs csdg csx cp cu1 cphase crz swap dcx rzz ccx ccz cswap mcx mcz mcp` (also `toffoli`, `mct`,
`c3x`, `c4x`, `mcphase`, `mcu1`, `fredkin`).  Measurements, resets, classical control, gates outside the list and
non-dyadic angles: `UNKNOWN`.

## 1. Representation

Following Amy [1, Definition 2.1], an operator acts on basis states as

    |x>  ->  √2^e · Σ_{y ∈ {0,1}^m} exp(2πi F(x,y) / 2^K) · |W(x,y)>

* x: free input variables; y: path variables (summed);
* F: multilinear polynomial in the Boolean variables with coefficients in **Z_{2^K}**.  Amy takes coefficients in
  the dyadic fractions modulo 1, D/Z [1, Def. 2.1]; all phases of a circuit with angles π·p/2^k lie in
  2^{-K}Z/Z ≅ Z_{2^K}, so everything is exact integer arithmetic (Python integers — no overflow; Amy reports the
  QFT up to 31 qubits, stopped by overflow of fixed-width dyadic arithmetic [1, Sec. 5.2]);
* W: one Boolean polynomial per wire in algebraic normal form (ANF);
* e: integer.

Both F and W are **canonical**: a function {0,1}^v → Z_{2^K} has exactly one multilinear representation (Möbius
inversion is an integer formula, so it holds over any commutative ring), a Boolean function exactly one ANF.
A Boolean polynomial P enters F through its integer lift [P] — the multilinear integer polynomial equal to P on
every Boolean point — built by [A ⊕ m] = [A] + m − 2[A]m [1, Sec. 2.1]; with a coefficient c = 2^j·odd only the lift
modulo 2^{K−j} matters.  K = max(3, largest k of the phases π·p/2^k of the expanded circuit and reference, plus 1,
width of a QFT / arithmetic specification).

### Gates and angles

Each gate is expanded into six exact primitives (`circuit.expand`): Hadamard (new path variable y,
F += 2^{K−1}[W_q]·y, W_q := y, e −= 1), XOR W_t ⊕= ΠW_c (X, CX, Toffoli, MCX), swap, controlled swap, the phase
e^{iπr[ΠW_q]} (Z, S, T, P, CP, MCP, CZ, CCZ; the global phase for the empty product) and e^{iπr[⊕W_q]} (RZZ).
Identities used: RZ(t) = e^{−it/2}P(t); RX(t) = H·RZ(t)·H; RY(t) = S·RX(t)·S†; U(θ,φ,λ) = P(φ)RY(θ)P(λ);
SX = HSH; Y = iXZ; CY = S_t·CX·S_t†; CRZ(t) = P_c(−t/2)·CP(t); CSX = H_t·CS·H_t; DCX = CX(a,b)·CX(b,a).
The inverse of a circuit is its reversed primitive list with every phase negated (H, XOR, swap, cswap are
self-inverse), so no gate-level inverse table can be wrong.

Angles are exact multiples of π (`Pi(p, q)`, or strings such as `'pi/4'`, `'-3*pi/2**5'` read without eval) and are
kept **modulo 4π**: the half-angle gates RZ, RX, RY, U, CRZ, RZZ have period 4π (RZ(t+2π) = −RZ(t); CRZ(t+2π) has a
relative phase −1 on the control).  A float θ (radians, as Qiskit and ProjectQ give them) is read as π·p/2^k only if
θ/π, reduced to (−2, 2], equals o + p/2^k (o ∈ {0, ±2}, 0 < |p| ≤ 2^20) within 32 units in the last place of θ/π;
otherwise the circuit is UNKNOWN.  Every float is itself a dyadic rational, so without the bound on p every angle
would be "recognised"; with it a random non-dyadic angle matches with probability about 10^{−6}.  For a source that
rounds its angles, `DCircuit(n, float_abs_tol=…)` gives an absolute tolerance (ProjectQ keeps 12 decimals,
`ANGLE_PRECISION = 12`, and re-rounds after every merge and decomposition: 2·10^{−11} rad is used); then a nonzero
angle below 64 × tolerance is refused as ambiguous.  The proof is about the circuit with the recognised exact angles;
the count and the largest deviation are returned in `float_angles`.

## 2. Rewrite rules (Amy [1, Fig. 2, Sec. 3.2], generalised from Z_8 to Z_{2^K})

In every rule y is a path variable that occurs on no wire (and is not frozen, see §4).  The side conditions are
integer tests; if a condition does not hold the rule is not applied.

**[Elim]** y does not occur in F: Σ_y g = 2g.  Remove y, e += 2.

**[HH]** every term of F that contains y has coefficient 2^{K−1}.  Then F = 2^{K−1}·y·Q_int + R with Q_int the
integer sum of the quotients, and e^{2πiF/2^K} = (−1)^{y·Q}·e^{2πiR/2^K} with Q = Q_int mod 2, the XOR of the
quotients (an F2 polynomial without y).  If Q = z ⊕ Q′ where z is a path variable that occurs in Q only as the linear
monomial z (so z ∉ Q′):

    Σ_{y,z} (−1)^{y(z ⊕ Q′)} g(z) = Σ_z 2·δ(z = Q′) g(z) = 2·g(Q′).

Drop the y terms, substitute z := Q′ in F (through the lift: c·z·m → c·[Q′]·m) and in W; remove y, z; e += 2.
If Q = 1 the whole sum is 0 (impossible for a unitary; reported `INTERNAL`).
*The condition is also necessary*: Σ_y e^{2πi yG/2^K} = 1 + e^{2πiG/2^K} is 0 or 2 at every point iff
2G ≡ 0 (mod 2^K) pointwise, iff (by canonicity) every coefficient of G is 0 or 2^{K−1}.

**[ω]** the linear term y has coefficient c = 2^{K−2} or 3·2^{K−2}, all other terms with y have 2^{K−1}:

    Σ_y e^{2πi y(c + 2^{K−1}Q)/2^K} = 1 ± i(−1)^Q = √2 · e^{±iπ/4} · e^{∓iπQ/2}

(check Q = 0, 1).  Drop the y terms, F += ±2^{K−3} ∓ 2^{K−2}[Q], e += 1.

Every rule removes at least one path variable, so reduction terminates.  Soundness is a direct calculation [1,
Prop. 3.1]; the implementation is tested against an exact enumeration of the operator before and after reduction on
600 random path sums (`test_rules_sound_on_random_path_sums`, every rule fires many times).  Amy's fourth rule
[Case] is not implemented (not needed here); circuits that would need it stay UNKNOWN.

## 3. Verdicts

The circuit is composed with the adjoint of the reference / specification (adjoint first: circuit ∘ ref† = I ⟺
circuit = ref; for ancillas inside the reference — the qoptct clean model — adjoint last, so that exactly the
inputs with ancillas 0 are checked) and the sum is reduced.

* **Normal form** (no path variable left): the operator is |x⟩ ↦ √2^e e^{2πiF(x)/2^K}|W(x)⟩ with canonical F, W.
  PROVED iff e = 0, every wire equals the expected polynomial (the input variable, 0 for an ancilla, the ANF of the
  map for a ClassicalSpec) and F is a constant F_0 (global phase 2πF_0/2^K).  Otherwise REFUTED with an exact
  witness: a minimal-degree monomial of the nonzero difference polynomial, set to 1, is an input where it differs
  (for the phase: compared with the all-zero input).  e ≠ 0 with no variables left is `INTERNAL`.
* **Exact enumeration** (≤ 16 inputs + path variables left): amplitudes are summed in Z[ζ], ζ = e^{2πi/2^K}, in the
  power basis 1, ζ, …, ζ^{2^{K−1}−1} (ζ^{2^{K−1}} = −1; the 2^K-th cyclotomic polynomial x^{2^{K−1}}+1 is
  irreducible, so the representation is unique).  PROVED iff all off-diagonal amplitudes are 0 and all diagonal ones
  are the same element.
* **Witness search** (refutation only; 7 inputs before the symbolic strategies, 24 after the first one fails): for a
  concrete basis input x the composition maps |x⟩ to a unit vector ψ(x) given by a sum over path variables only.
  The output wires are fixed to the expected basis state, the rules reduce the then internal variables, a remainder of
  ≤ 18 variables is enumerated in Z[ζ]: this gives the amplitude A(x) = ⟨expected(x)|ψ(x)⟩ exactly.
  |A(x)|² ≠ 1 (exact test in Z[ζ]) or A(x) ≠ A(x′) proves that the circuit differs.  Inputs: all-0, all-1,
  alternating, single bits, seeded random.  It uses only unitarity (|A| = 1 ⟺ ψ = A|expected⟩).
* Anything else is UNKNOWN (rules stuck on a large sum, time, terms, ANF size).

With `selfcheck=True` the independent dense simulator (`sim.py`: gate matrices from their definitions, no code shared
with the path sum) re-checks every answer up to 10 qubits; a disagreement turns the answer into `INTERNAL`.

## 4. Specifications

**QFTSpec(wires, swaps=True, inverse=False)** — F|x⟩ = 2^{−m/2} Σ_y e^{2πi XY/2^m}|y⟩, X = Σ_j 2^j x_j
(little endian, wire j = bit j, as in Qiskit); Amy's QFT specification [1, Example 2.2].  `swaps=False`: output bit
y_k on wire m−1−k (Qiskit `QFT(do_swaps=False)`; ProjectQ's `QFT` is this one).  As a path sum the phase is
Σ_{j+k<m} 2^{K−m+j+k}·x_j·y_k.  No lemma is involved.

**ArithSpec(target, terms)** — |x, s⟩ ↦ |x, (s + G(x)) mod 2^m⟩ on m target wires, for every s;
G(x) = Σ_i w_i Π_{q∈Q_i} x_q (integer weights; term wires disjoint from the target).
`adder_spec(a, b, carry)`: G = Σ_j 2^j a_j (Qiskit DraperQFTAdder; `kind='half'`: the carry wire is the top target
bit); `const_adder_spec(c, b, controls)`: G = c·Πcontrols (ProjectQ AddConstant); `multiplier_spec(a, b, out)`:
G = Σ_{i,j} 2^{i+j} a_i b_j (Qiskit RGQFTMultiplier).  The prover uses its exact Fourier path sum

    A |x, s⟩ = 2^{−m} Σ_{Y,Z ∈ [0,2^m)} e^{2πi (S + G(x) − Z)·Y / 2^m} |x, Z⟩

(Y = Σ 2^k y_k, Z = Σ 2^j z_j, S the integer of the target input; only the terms of total weight < 2^m survive
modulo 1).

> **Lemma 2 (Fourier shift).** A is the classical map |x,s⟩ ↦ |x,(s+G(x)) mod 2^m⟩.
>
> *Proof.* For an integer L, Σ_{Y=0}^{2^m−1} e^{2πi LY/2^m} = 2^m if 2^m | L and 0 otherwise (geometric series
> with ratio ρ = e^{2πiL/2^m}: (ρ^{2^m} − 1)/(ρ − 1) = 0 for ρ ≠ 1).  As (y_0,…,y_{m−1}) runs over {0,1}^m, Y runs
> over [0, 2^m) once.  With L = S + G(x) − Z: A|x,s⟩ = Σ_Z [Z ≡ S + G(x) mod 2^m] |x,Z⟩, and exactly one Z in
> [0, 2^m) satisfies the congruence.  ∎

So `prove_spec(c, ArithSpec)` = PROVED means: the circuit equals the path sum A exactly (machine-checked), and A is
the classical map by Lemma 2; the answer lists the lemma under `assumptions`.  For small sizes the same circuits are
also proved against the classical map directly (`PermutationSpec`, no lemma) and re-checked by the dense simulator.

The adjoint A† |x,z⟩ = 2^{−m} Σ_{Y,S′} e^{2πi(Z − G − S′)Y/2^m}|x,S′⟩ is applied before the circuit.  Strategy
**freeze**: the frequency variables Y are excluded from the rules while the circuit is added (and reduced gate by
gate), then released.  Without it, [HH] on Y would first expand the carries of the classical subtraction
(exponential ANF); with it, the circuit's QFT cancels against the specification first (substitutions go towards the
frozen variables), its rotations cancel the −G·Y terms, and finally its inverse QFT cancels the remaining sum.

**ClassicalSpec(n, func, zero_wires=())** / **PermutationSpec(n, func)** — |b⟩ ↦ |func(b)⟩ for the basis states whose
zero_wires are 0 (at most 20 free inputs; func injective there, checked).  Expected wires: the ANF of each output bit
over the free inputs (fast Möbius transform of the truth table).  No lemma.  Example: Qiskit's HRSCumulativeMultiplier
adds the partial products with half adders whose carry-out must start in 0, so its claim is
|a, b, 0⟩ ↦ |a, b, a·b⟩: `ClassicalSpec(4n, f, zero_wires=out)` (against `multiplier_spec` for every `out` it is
REFUTED for n ≥ 2, which is correct).

## 5. Using it from a bridge (Q1: Qiskit, P: ProjectQ)

1. Build a `DCircuit(n)` with `c.add(name, qubits, *angles)` (controls first, Qiskit conventions) and
   `c.add_global_phase(angle)`.  Prefer exact angles `Pi(p, 2**k)`; floats in radians are accepted (see §1).
   Qiskit: `qoptct.dyadic.qiskit_io.from_qiskit(qc)` unrolls library blocks by their definitions (and checks
   control states).  Do not `transpile` with the default `qubits_initially_zero=True`: it may use idle qubits as
   clean ancillas and keeps only the action on |0…0⟩; use `qubits_initially_zero=False`.
   ProjectQ: rotation angles are rounded to 12 decimals by ProjectQ itself.  The bridge of 0.8
   (`qoptct.foreign.projectq.dyadic`, `check_projectq_dyadic`) accepts an angle only if it is exactly what ProjectQ stores for
   p·π/2^k with k ≤ 38 (the multiples of π/2^k are then further apart than the rounding step, so the reading is unique) and
   refuses 0, angles within 1e-11 of 0 or of the period, and everything else (`UNKNOWN`), and does not use the tolerance
   `float_abs_tol` that the dump scripts of the first round used.  Rotations below 5·10^{−13} rad are rounded to 0 by ProjectQ, so its
   QFT and AddConstant on 40 or more qubits are `UNKNOWN` (QFT on 64 qubits: 231 of 2016 controlled R have angle 0).
   Qiskit: `check_qiskit(circuit, spec_or_reference)` calls the prover itself (`qoptct/foreign/qiskit/dyadic_bridge.py`).
2. Call `prove_dyadic(c, spec_or_reference)`.  Only `verdict == 'PROVED'` is a proof; report `assumptions`
   (the Fourier lemma for ArithSpec), `global_phase` and `float_angles` with it.
3. A circuit outside the gate set or with measurements: UNKNOWN with the reason; convert first or use the Clifford+T
   provers.

## 6. Results (registry `results/dyadic/*.json` from `scripts/dyadic_runs.py` and `scripts/projectq_dyadic_runs.py`, lab Mac, one process)

Command line (0.8): `qoptct dyadic prove CIRCUIT.qasm --reference REF.qasm` or `--qft`, `qoptct dyadic results`; exit code 0 only for PROVED.

| circuit (library) | specification | proved up to | qubits / gates at max | prove time, peak memory | canaries at max | first size not proved |
|---|---|---|---|---|---|---|
| Qiskit `QFTGate` (with swaps) | `QFTSpec` | n = 256 | 256 / 33 024 | 1.0 s, 102 MB | 22 / 22 REFUTED | — |
| Qiskit `synth_qft_full(do_swaps=False)` | `QFTSpec(swaps=False)` | n = 256 | 256 / 32 896 | 1.0 s, 102 MB | 21 / 21 | — |
| Qiskit `QFTGate.inverse()` | `QFTSpec(inverse=True)` | n = 256 | 256 / 33 024 | 1.0 s, 127 MB | 22 / 22 | — |
| Qiskit QFT · QFT† | identity | n = 256 | 256 / 66 048 | 1.7 s, 141 MB | 22 / 22 | — |
| Qiskit `DraperQFTAdder` fixed | `adder_spec` (+ lemma) | n = 256 | 512 / 98 944 | 9.9 s, 191 MB | 22 / 22 | — |
| Qiskit `DraperQFTAdder` half | `adder_spec` with carry (+ lemma) | n = 256 | 513 / 99 714 | 9.9 s, 191 MB | 22 / 22 | — |
| Qiskit `DraperQFTAdder(n).control(1)` | controlled `ArithSpec` (+ lemma) | n = 10 | 21 / 895 | 48 s, 125 MB | 20 / 20 | REFUTED at n = 64 (see below) |
| Qiskit `RGQFTMultiplier` | `multiplier_spec` (+ lemma) | n = 32 | 128 / 69 760 | 1.8 s, 163 MB | 22 / 22 | — |
| Qiskit `HRSCumulativeMultiplier` (p(±π/8)) | `ClassicalSpec`, out = 0 | n = 7 | 28 / 4 634 | 166 s, 92 MB | 20 / 20 | n = 8 UNKNOWN (time) |
| ProjectQ `QFT`, exact-angle bridge (`pq_exact_qft`) | `QFTSpec(swaps=False)` | n = 39 | 39 / 780 | 0.2 s | 21 / 21 | n = 40 UNKNOWN (angle 6e-12 is ambiguous) |
| ProjectQ `QFT` (native, first round, tolerance 2e-11) | `QFTSpec(swaps=False)` | n = 32 | 32 / 528 | 0.04 s | 21 / 21 | n = 48 UNKNOWN (angles) |
| ProjectQ `AddConstant(c)`, exact-angle bridge (`pq_exact_addconst`) | `const_adder_spec` (+ lemma) | n = 38 | 38 / 1 961 | 0.3 s | 21 / 21 | n = 40 UNKNOWN |
| ProjectQ `AddConstant(c)` (native, first round) | `const_adder_spec` (+ lemma) | n = 32 | 32 / 1 351 | 0.07 s | 21 / 21 | n = 48 UNKNOWN (angles) |
| ProjectQ `QFT` (basic: CNOT + 1-qubit) | `QFTSpec(swaps=False)` | n = 31 | 31 / 2 356 | 0.08 s | 20 / 20 | n = 32 UNKNOWN (angles) |
| ProjectQ `AddConstant(c)` (basic) | `const_adder_spec` (+ lemma) | n = 16 | 16 / 1 292 | 0.06 s | 20 / 20 | n = 24 UNKNOWN (angle numerator > 2^20) |

Independent checks: every size ≤ 10 qubits re-checked by the dense simulator (`dense_check: agrees`) and, for Qiskit
circuits, the translation compared with Qiskit's own `Operator`; at the largest size the product-state simulator
(`prodsim.py`) ran 20 random basis inputs with 0 errors (it cannot follow HRS, whose Toffoli decompositions
entangle).  `results/dyadic/fuzz.json`: 3000 random circuit pairs and 1000 random circuit/specification pairs up to
10 qubits, 0 disagreements with the dense simulator.

Findings.  (1) Qiskit 2.6.2 `control()` of `CPhase(θ)` gives the angles (θ, θ, θ) instead of (θ/2, −θ/2, θ/2) when
|θ| ≤ π/2^51, so `DraperQFTAdder(n).control(1)` is not exactly a controlled adder for n ≥ 52 (REFUTED at n = 64 with
the all-ones witness; the error is below 10^{−15} rad, invisible to numeric checks;
`results/dyadic/cdraper_n64_finding.json`, `test_controlled_cp_pattern_tiny_angles`).  (2) ProjectQ keeps 12 decimals
of every angle and rounds rotations below 5·10^{−13} rad to 0 (its QFT on 64 qubits has 231 of 2016 controlled R with
angle 0): from n ≈ 43 on its QFT and AddConstant are not exact, and from n = 33 the float angles cannot be read back
uniquely.  (3) Qiskit `transpile` with the default `qubits_initially_zero=True` may borrow idle qubits as clean
ancillas; the result is equal only on |0…0⟩ (HRS n = 2: 4 inputs with out = 0 differ).

## 7. Limits

* Unitary circuits only (plus clean ancillas).  No measurements or classical control.
* The rules are incomplete (as in [1]): circuits whose path sum does not cancel (e.g. a ripple-carry multiplier
  against the Fourier specification) stay UNKNOWN; classical maps are expanded to ANF (exponential in the carry
  chain), so ClassicalSpec / HRS go only to small n.
* Limits: `max_terms` (terms of F and monomials on all wires), ANF size per wire 100 000, time budget.

## 8. Reproduce

```
python -m pytest -q tests/test_dyadic.py -m "not slow"     # seconds
python -m pytest -q tests/test_dyadic.py                   # + 1500 random circuits up to 10 qubits, n = 64 canaries
python scripts/dyadic_projectq_dump.py                     # ProjectQ gate lists (needs projectq), also --mode basic,
                                                           # --optimizer on
python scripts/dyadic_runs.py registry                     # results/dyadic/<family>.json, one process per size
python scripts/dyadic_runs.py run draper_fixed 64 --canaries
```

## References

1. M. Amy, *Towards Large-scale Functional Verification of Universal Quantum Circuits*, Proc. QPL 2018, EPTCS 287
   (2019), doi:10.4204/EPTCS.287.1, arXiv:1805.06908v2 — Definition 2.1 (path sums), Sec. 2.1 (lifting),
   Fig. 2 / Sec. 3.2 (rules [Elim], [HH], [ω], [Case]), Prop. 3.1 (soundness), Example 2.2 (QFT specification),
   Sec. 5.2 (QFT verified up to 31 bits).
2. T. G. Draper, *Addition on a Quantum Computer*, arXiv:quant-ph/0008033 (2000).
3. L. Ruiz-Perez, J. C. Garcia-Escartin, *Quantum arithmetic with the quantum Fourier transform*, Quantum Inf.
   Process. 16, 152 (2017), arXiv:1411.5949.
4. S. Beauregard, *Circuit for Shor's algorithm using 2n+3 qubits*, arXiv:quant-ph/0205095 (2003) (constant
   addition in Fourier space, as in ProjectQ's `add_constant`).
5. T. Häner, M. Roetteler, K. M. Svore, *Factoring using 2n+2 qubits with Toffoli based modular multiplication*,
   arXiv:1611.07995 (2017) (Qiskit HRSCumulativeMultiplier).
