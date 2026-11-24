# Qiskit bridge (qoptct 0.8, subpackage `qoptct.foreign.qiskit`)

Proves that a Qiskit `QuantumCircuit` without non-Clifford+T rotations computes a stated function (or equals a
reference circuit) for **all** inputs, using the unchanged qoptct provers (`qoptct.verify.kprove`, BDD sum over paths
with symbolic measurement outcomes).  Qiskit is the only extra dependency (already required by qoptct).

## API

```python
from qoptct.foreign.qiskit import check_qiskit, Spec, translate
r = check_qiskit(circuit, spec, clean=("help",))        # spec: Spec or classical function f(vals) -> dict
r = check_qiskit(circuit, reference=other_circuit, clean=(...), ref_clean=(...))   # miter circuit * reference^-1
r["status"]    # "PROVED" | "REFUTED" (r["witness"]) | "UNKNOWN"
r["needs"]     # "fractional_phases" when a rotation that is not a multiple of pi/4 stopped the translation
```

* `translate(circuit)` expands every instruction **only from its own `definition`** down to leaves whose exact
  matrix was compared with a rule table: `x cx ccx swap h s sdg t tdg z cz`, `p/u1(k*pi/4)` (as `t^k`), `mcx`
  (C3X, C4X, `MCXGate`, all-ones control state; matrix checked up to 9 wires, attributes above), `measure`, and
  `if_else` on one bit.  Nothing is optimised or recognised.  Any other rotation (`p(theta)` with theta not a multiple
  of pi/4, `u`, `r`, `rx`, `ry`, `rz`) gives `NeedsFractionalPhases` and the verdict `UNKNOWN`, never `PROVED`.
* Global phase of the Qiskit circuit is recorded; every verdict is "up to a global phase", the contract of the
  qoptct provers (the phase must not depend on the input or on measurement outcomes).
* `Spec(name, classical, symbolic)`: `classical(vals)` maps register name to integer (clean registers are 0) and
  returns the registers that change; `symbolic(S, R)` is the same function over `qoptct.verify.ksym.SymInt` words
  (`R[name]`), used by the BDD prover.  Unlisted input registers must stay unchanged, clean registers must return
  to 0.
* Independent checks before and besides the proof: the translated leaves are replayed by a separate classical
  simulator on all inputs (up to 18 free bits: that alone is an exhaustive proof for a classical circuit) or on
  64+ edge and random inputs; a mismatch is a `REFUTED` with a replayed witness.  The tests also compare the
  translation with Qiskit's own state-vector evolution.
* `mcx='expand'` expands multi-controlled X from Qiskit's definition instead of keeping it as a primitive.

## Proved (Qiskit 2.6.2), registry `qoptct/foreign/qiskit/data/registry.json`

Reproduce: `python -m qoptct.foreign.qiskit prove FAMILY N` (`... list` shows the families; `--canaries` runs the
mutants).  All rows below are PROVED by the BDD prover for **every** input and every value of the dirty outputs:

| family | contract | max n proved |
|---|---|---|
| `cdkm-half/full/fixed` (CDKMRippleCarryAdder) | `b := (cin + a + b) mod 2^n`, `cout ^= carry` | 256 |
| `vbe-half/full` (VBERippleCarryAdder) | same | 256 |
| `half-adder-gate`, `full-adder-gate`, `modular-adder-gate` | sum of the gate, `ModularAdderGate` is Clifford+T (h, t, tdg) | 256 (`modular-adder-gate`: canaries at n = 32 only, mutants at n = 256 exceed the prover's memory) |
| `comparator-geq`, `comparator-lt` (IntegerComparator) | `compare ^= [state >= v]` / `< v` | 256 |
| `comparator-gate` (IntegerComparatorGate) | same | 16 (the BDD exceeds 3 GB at 64); **no canaries** (one proof of its 415 073 gates takes tens of minutes; `prove comparator-gate 16 --canaries --sample 26` produces them) |
| `and`, `or`, `xor-const` (AndGate, OrGate, BitwiseXorGate) | boolean function | 256 |
| `c-cdkm-*`, `c-vbe-half/full` | controlled versions, see below | 256 |

Controlled versions: Qiskit's generic `Gate.control(1)` of a whole adder inserts fractional multi-controlled
phases (see the list below).  The registry versions are built by taking Qiskit's own `.control(1)` of every primitive
of the flattened adder (x -> cx, cx -> ccx, ccx -> C3X, swap -> cswap) and are proved against `ctrl ? f : identity`.

### Defect found in Qiskit

`VBERippleCarryAdder(n, kind="fixed")` (Qiskit 2.6.2) does not compute `(a + b) mod 2^n` for n >= 2: 8 of 16 inputs
at n=2, 16 of 64 at n=3, 64 of 256 at n=4 are wrong (Qiskit's own state-vector simulation; example n=3: a=2, b=2).
`check_qiskit` returns `REFUTED` with this witness (family `vbe-fixed`, kept as a negative test).

## Mutants (canaries)

`qoptct.foreign.qiskit.mutants.make_mutants` rebuilds the flattened circuit with one defect: a dropped gate at 13
positions (seams, first, last, quartiles), a moved or swapped gate, a wrong control or wrong target at five places,
an extra X / CX, a phase-only defect (z, t, s) and a defect on one rare input (all-ones multi-controlled flip).  A
mutant that the prover accepts must also be accepted by the reference miter (genuinely equivalent, for example
`s` on a clean ancilla at the start, or the last CNOT of the 2-CNOT UMA when the carry-in is 0); everything else must
be `REFUTED`.  Canaries at the largest size are in the registry (`canaries` of each row; rows without them say so).  For short circuits (a single wide AND, a row of X gates) extra mutants are added (dropped control, shifted target, extra X / Z); `--sample N` takes N mutants spread over the whole list and builds them one at a time (huge circuits).

## Qiskit circuits with fractional phases: the dyadic prover (0.8)

`qoptct qiskit fractional` prints the list below.  Such a circuit stays `UNKNOWN` (`needs == 'fractional_phases'`) while the call
gives no target; with a dyadic specification (`QFTSpec`, `ArithSpec` via `adder_spec` / `multiplier_spec`, `ClassicalSpec` from
`qoptct.dyadic`) or a reference `QuantumCircuit`, `check_qiskit` hands the circuit to `qoptct.dyadic.prove_dyadic`
(`qoptct/foreign/qiskit/dyadic_bridge.py`):

```python
from qoptct.dyadic import QFTSpec, adder_spec
check_qiskit(QFTGate(8) in a circuit, QFTSpec(range(8)))["status"]                 # PROVED
check_qiskit(DraperQFTAdder(4, kind="fixed").decompose(), adder_spec(range(4), range(4, 8)))   # PROVED, assumptions: fourier-shift lemma
```

The answer carries `method == 'dyadic:...'`, the answer of the prover under `r['dyadic']` and its `assumptions`.  Clean registers are
ancillas.  A float angle that is not (within 32 ulp, after an exact reduction modulo 4 pi) a multiple of pi/2^k, an unsupported gate or a measurement give `UNKNOWN`.
The registers of a reference circuit are matched to those of the circuit by name (the same rule as in the Clifford+T path); a gate is taken by its name only when it is a standard Qiskit gate class with the matching matrix, a user gate is unrolled by its definition.
Overview: [FOREIGN.md](FOREIGN.md); results: `qoptct dyadic results`.  At n=4 the following give `UNKNOWN` without a target ("needs
fractional phases"): `DraperQFTAdder` (half, fixed), `RGQFTMultiplier`, `HRSCumulativeMultiplier`, `MultiplierGate`,
`QFTGate`, `QuadraticForm`, `LinearPauliRotations`, `PiecewiseLinearPauliRotations`, `PolynomialPauliRotations`,
`ExactReciprocalGate`, and every native `Gate.control(1)` of `CDKMRippleCarryAdder`, `VBERippleCarryAdder`,
`HalfAdderGate`, `FullAdderGate`, `ModularAdderGate`, `MultiplierGate`, `IntegerComparator`,
`IntegerComparatorGate`.  `PhaseEstimation` (power of a unitary) is unsupported for another reason.
`WeightedAdder` translates (Clifford+T) but has no specification in the registry yet.

## Limits

* No hand-written rotations, no `reset`, no real-time variables, `if_else` only on one bit with an x/z/cx/cz body.
* A classical spec function alone proves only by exhaustive replay (<= 18 free bits); without `Spec.symbolic` the
  result for larger circuits is `UNKNOWN` with a sampled check, never `PROVED`.
