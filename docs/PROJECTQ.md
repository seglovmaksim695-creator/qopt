# ProjectQ bridge (`qoptct.foreign.projectq`)

Proofs about ProjectQ programs with the qoptct provers. ProjectQ is an optional dependency:

    pip install "qoptct[projectq]"          # ProjectQ 0.8.0; libs.revkit additionally needs the external `revkit` package

Importing `qoptct.foreign.projectq` never imports ProjectQ; using it without ProjectQ raises a clear `ImportError`.

## What it does

1. The program (a function `program(eng)`) runs on a ProjectQ `MainEngine` whose last engine records commands. ProjectQ's own
   decomposition rules (`projectq.setups.decompositions`, `projectq.libs.math`) bring every command into
   * `mode="nct"`: X with at most two controls (X, CNOT, Toffoli), or
   * `mode="clifford_t"`: X with at most one control, H, S, Sdg, T, Tdg, Z, CZ.
2. Every command is translated by exactly one rule (`_logical.py`): X, CNOT; a gate with k >= 2 controls becomes a chain of
   logical ANDs on fresh scratch wires (Gidney's AND and its measurement-based erasure, as in the Qualtran bridge); the
   Clifford+T atoms stay atoms. A single-qubit gate is accepted only if its 2x2 matrix is exactly (no global phase) the matrix of
   an atom. Measurements, dirty qubits, rotations and unknown gates give `UNKNOWN`, never `PROVED`.
3. Provers:
   * `ClassicalSpec(fn, promise=None)`: all 2^n inputs, bit-parallel (the OpenQASM text that is emitted is parsed back and
     simulated); needs n <= `max_input_bits` (20).
   * `GF2MulSpec(m, poly)`: `qoptct.verify.exact2`, for all inputs at any size (degree <= 2 circuits); the polynomial is an
     explicit argument and is tested for irreducibility with Rabin's exact test (any degree).
   * `mode="clifford_t"`: the NCT version of the program is proved against the specification, and the Clifford+T version is proved
     equal to it with the exact path-sum prover `qoptct.verify.prove` (clean-ancilla model).
4. Contract: the outputs equal the specification; every other wire (including freed ancillas and scratch wires) ends in |0>.
   Unspecified garbage is only allowed through a **pin** (below).
5. Independent check: the genuine program is run on ProjectQ's own `ClassicalSimulator` (inputs injected after allocation; when
   the simulator cannot run a gate, after ProjectQ's decomposition) on all inputs (<= 10 input bits) or on corner and random inputs.
   A disagreement turns `PROVED` into `UNKNOWN`.

```python
from qoptct.foreign.projectq import check_projectq, ClassicalSpec

def program(eng):
    from projectq.libs.math import AddQuantum
    a = eng.allocate_qureg(4); b = eng.allocate_qureg(4)
    AddQuantum | (a, b)
    return dict(inputs=dict(a=a, b=b), outputs=dict(a=a, b=b))

r = check_projectq(program, ClassicalSpec(lambda v: dict(a=v['a'], b=(v['a'] + v['b']) % 16)))
r['status']            # 'PROVED'; REFUTED comes with r['witness']
```

## Garbage pins

The Binary_ECC multipliers and point addition leave non-zero values on many wires (the upper bits of the Karatsuba product and
intermediate sums). `pin_garbage(program, spec)` records the final value of every such wire of the golden circuit; with
`pin=...` a candidate must end in exactly these values. The values are compared, not proved meaningful. (In the first campaign
the garbage of the circuit under test was taken as its own reference, and one mutant at m = 163 was PROVED; with a pin it is
REFUTED.)

## Canaries

`qoptct.foreign.projectq.canaries` makes deliberately wrong versions: dropped / retargeted / reversed / duplicated gates, seams
(the first or last use of a wire), a wrong or missing erasure correction, an AND replaced by a CNOT, a wrong phase in the AND
gadget (Sdg to S, Tdg to T), changes visible on a small part of the inputs, flipped output or garbage bits. A mutant that is
PROVED is checked for equivalence with the golden program (identical value on every wire); an equivalent mutant (e.g. a Toffoli
control that is always 1) is reported as `EQUIVALENT` and is not counted. A result is registered only if at least 20 canaries
were decided, none was PROVED and at least half were targeted (seam, erasure, phase, rare).

## Binary_ECC

The sources (Eclipse Public License 2.0) are not part of this repository. Point `--source` or `$QOPT_BINARY_ECC` at the
`Binary_ECC` folder; the text before the first `global resource_check` is executed unchanged with the authors' flags
(`NCT=1`), and the SHA-256 of the file is in every result.

    python -m qoptct.foreign.projectq prove multiplier --m 163 --source /path/to/Binary_ECC
    python -m qoptct.foreign.projectq prove point-addition --m 8 --source /path/to/Binary_ECC
    python -m qoptct.foreign.projectq prove library --block add_quantum --n 6
    python -m qoptct.foreign.projectq results [--name NAME] [--size N] [--long] [--json]

Point addition: the claim is for q = 0 (unchanged) and for q = 1 with x1 != x2; at x1 = x2 the circuit inverts by Fermat's
little theorem and the value is not a point addition. The registry is `data/registry.json`, built by
`scripts/build_projectq_registry.py` from `results/projectq/*.json`.

## Library blocks (`library.py`)

AddQuantum (plain, with carry, controlled), SubtractQuantum, ComparatorQuantum, MultiplyQuantum, DivideQuantum, multi-controlled
X, controlled swap and `libs.revkit.PermutationOracle`, each against plain Python arithmetic. Findings about ProjectQ 0.8.0
itself (decomposition rules, not our code; reproduction in `FINDINGS/` of the round folder, nothing was sent upstream):
`ComparatorQuantum` decomposes wrongly for n = 1 and for n >= 3 (for n >= 3 a and b are changed; n = 2 is correct), `MultiplyQuantum` decomposes wrongly when b >= 2.
The gate itself (its `math_fun`) is correct; the disagreement is between the gate and its decomposition.

## Blocks with rotations (not here)

AddConstant / SubConstant, AddConstantModN / SubConstantModN, MultiplyByConstantModN (Draper QFT adder with phases R(pi/2^k)),
QFT, Rx / Ry / Rz / R / Ph with angles that are not multiples of pi/4. They need the fractional-phase prover. `PhaseOracle`
of `libs.revkit` produces phase gates and is not covered either.

## Assumptions

Exact atom matrices; one explicit rule per command; the AND gadget lemma (`exact2.gadget_check`); ProjectQ's decomposition rules are
not audited (what they produce is what is proved); a freed qubit must end in |0>.
