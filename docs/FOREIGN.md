# Proofs about circuits from other libraries (0.8)

qoptct proves, with its own exact provers, that circuits built by other quantum libraries do what they claim.  Nothing here
changes the provers of `qoptct.verify`; every new part is a new module, and every optional library is an extra (`qoptct[qualtran]`,
`qoptct[projectq]`; Qiskit is a dependency of qoptct itself).  Without an extra qoptct installs and works as before; the tests of
that part are skipped.

| Part | Command | Module | What it proves | Details |
|---|---|---|---|---|
| Qualtran | `qoptct qualtran ...` | `qoptct.foreign.qualtran` | Kaliski inversion, Montgomery multiplier, point addition, forward + `uncompute` pairs (n up to 256) | [QUALTRAN.md](QUALTRAN.md) |
| Qiskit | `qoptct qiskit ...` | `qoptct.foreign.qiskit` | adders, comparators, `AndGate`, `OrGate`, `BitwiseXorGate` and controlled adders at n = 256 (all inputs); rotations through `dyadic` | [QISKIT.md](QISKIT.md) |
| ProjectQ | `qoptct projectq ...` | `qoptct.foreign.projectq` | Binary_ECC multipliers m = 8 ... 283, point addition m = 8, `libs.math` / `libs.revkit` blocks; rotations through `dyadic` | [PROJECTQ.md](PROJECTQ.md) |
| Rotations pi*p/2^k | `qoptct dyadic ...` | `qoptct.dyadic` | QFT, Draper adders, QFT multipliers (Qiskit up to n = 256, ProjectQ up to n = 38 / 39) | [DYADIC.md](DYADIC.md) |
| Generic parse | `qoptct auto ...` | `qoptct.foreign.auto` | Qualtran blocks from their decomposition, without code per block (14 results at n = 256) | [AUTO.md](AUTO.md) |
| The circle | `qoptct loop ...` | `qoptct.loop` | propose, prove, accumulate: a journal of proved improvements | [LOOP.md](LOOP.md) |
| All results | `qoptct results` | `qoptct.results` | one table over all registries shipped with the package | below |

Every verdict is PROVED, REFUTED (with a witness) or UNKNOWN; the exit code of a command is 0 only for PROVED (or PASS).  A
result holds for the stated size, library version and assumptions only.  Every PROVED at the largest size comes with at least 20
canaries: damaged variants of the same circuit (seams between blocks, erasure corrections, phases, rare inputs) that must not be proved.

## Fractional phases: Qiskit and ProjectQ through `qoptct.dyadic`

The Clifford+T provers read phases that are multiples of pi/4.  The QFT, Draper's adder, QFT multipliers and the constant adders
of ProjectQ use rotations pi*p/2^k.  `qoptct.dyadic` is a path-sum prover whose phases are integers modulo 2^K (Amy, arXiv:1805.06908;
the rules are generalised from Z8 and applied only when their condition is checked exactly).  The two bridges hand such circuits to it:

* **Qiskit.**  `check_qiskit(circuit, spec_or_reference)`.  A rotation that is not a multiple of pi/4 used to end in `UNKNOWN`
  (`needs == 'fractional_phases'`).  Now, when the call gives a dyadic specification (`QFTSpec`, `ArithSpec` with `adder_spec` /
  `multiplier_spec`, `ClassicalSpec`) or a reference `QuantumCircuit`, the circuit is unrolled to its Qiskit definition and proved by
  `prove_dyadic`.  Without such a target the answer stays `UNKNOWN` with `needs == 'fractional_phases'`.  Clean registers
  (`clean=[...]`) are ancillas that start and end in 0.  The verdict is up to a global phase.  A float angle is read as pi*p/2^k within a few
  units in the last place; an angle that is not close to such a value gives `UNKNOWN`, never `PROVED` (a corrupted angle
  that is another dyadic value gives `REFUTED`).
* **ProjectQ.**  `qoptct.foreign.projectq.dyadic.check_projectq_dyadic(program, target)`.  The program runs on a ProjectQ engine with
  ProjectQ's own decomposition rules and no optimizer.  ProjectQ stores every angle as a float rounded to 12 decimals and rounds an
  angle below 5e-13 to 0 (a rotation by pi/2^43 is a gate with angle 0).  An angle is therefore accepted only if it is **exactly** what
  ProjectQ stores for p*pi/2^k with k <= 38 (or for the inverse gate, which ProjectQ builds from the rounded angle); an angle that is 0,
  within 1e-11 of 0 or of the period, or matches no such value (the sum of two rounded angles, a corrupted angle) is `UNKNOWN`.  So the
  QFT is proved up to 39 qubits and `AddConstant` up to 38; on 40 qubits or more the program is `UNKNOWN`.
  Allocated-and-freed wires are clean ancillas.

Assumptions that stay in the answer (`assumptions`): for `ArithSpec` (addition, multiplication in the Fourier basis) the
fourier-shift lemma (a geometric sum, proved in `DYADIC.md`, Lemma 2; the circuit's equality with the exact path sum of the
arithmetic is checked by the machine); for ProjectQ the exact angle reading above.  The proofs are about the circuit with the exact
dyadic angles; the real floats of the libraries differ from them by their rounding.

Findings on the way (described in the documents of the parts; minimal reproductions without qoptct: `docs/findings/f6_projectq_decomposition_bugs.py`,
`docs/findings/f7_qiskit_vbe_fixed.py`): Qiskit 2.6.2 `VBERippleCarryAdder(kind='fixed')` is wrong for n >= 2 (an
open Qiskit pull request exists); Qiskit 2.6.2 `control()` of a `CPhase` with an angle of at most pi/2^51 gives wrong angles
(`DraperQFTAdder(64).control(1)` is REFUTED); ProjectQ 0.8.0 `ComparatorQuantum` is wrong for n = 1 and n >= 3 and `MultiplyQuantum` for b >= 2;
Qiskit `transpile` with the default `qubits_initially_zero=True` keeps only the action on |0...0>.

## Generic parse: `qoptct auto`

`qoptct.foreign.auto` descends into the decomposition of a Qualtran block by itself.  Leaves go to the existing exact prover
(ROBDD with exact AND / erasure domains); nodes are proved by composing the contracts of the children in Z3, with the checks of the
audits of round 3: the input of every piece is the output of the previous one (wire flow), the pieces cover all operations and foreign
wires are fresh, simplicity and the parameter domain are checked before any computation, and a contract that is not proved gives
`PROVED_CONDITIONAL`, not `PROVED`.  A loop needs an invariant formula, otherwise it is not PROVED.  Per block only data is written:
the specification of the root, the contracts of the primitives (88 instances compared with `call_classically`) and the invariant.
It repeats the results of 0.7 (Montgomery w = 1, 2, the first step of ECAdd: the same program hashes) and adds eleven blocks at n = 256
(`ModAdd`, `CModAdd`, `ModSub`, `CModSub`, `ModNeg`, `CModNeg`, `ModDbl`, `CModAddK`, `Subtract`, `SubtractFrom`, `Negate`).
Not done: `CtrlScaleModAdd` and `CModMulK` above n = 8, steps 2-6 of ECAdd, `ECWindowAddR`.  See [AUTO.md](AUTO.md).

## The circle: `qoptct loop`

`qoptct loop run` lets proposers (parameter search for lookup and lookup-add, passes of the optimizer) suggest a circuit; the provers
of qoptct check it against the reference of its target; only a proved candidate with an independent simulation and at least 20 canaries
enters the journal (`journal.jsonl`, a chain of hashes).  The table of the best known circuits is built from the journal only.
`qoptct loop verify --state DIR` re-checks the chain, the file hashes, the costs and every proof.  The shipped journal
(`results/loop/state`, 116 records) re-checks with `qoptct loop verify --state results/loop/state`.  See [LOOP.md](LOOP.md).

## All results: `qoptct results`

```
qoptct results                       # the table of everything proved (registries of the installed package)
qoptct results --library dyadic      # one library; repeatable
qoptct results --long                # one block per result with the command that reproduces it
qoptct results --json                # the same as data
qoptct results --all                 # also the rows of the registries that are not proved
```

`qoptct results --all` lists 125 rows in the 6 registries (33 Qualtran, 31 ProjectQ, 21 Qiskit, 18 auto, 11 loop, 11 dyadic); the table without `--all` has only the proved ones.

The registries are data files of the package (`foreign/qualtran`, `foreign/qiskit`, `foreign/projectq`, `foreign/auto`, `dyadic`, `loop`),
built from the result files of the repository by the `scripts/build_*_registry.py` scripts (each has `--check`); the table needs no
optional library.  Exit code 0 when the table was printed.

## Reproduce

```
qoptct qiskit prove cdkm-half 256
qoptct projectq prove library --block add_quantum --n 8          # needs the extra `projectq`
qoptct dyadic prove circuit.qasm --qft                           # circuit.qasm: OpenQASM with rotations pi*p/2^k
qoptct auto prove modadd --n 256 --p "2**255-19"                 # needs the extra `qualtran`
qoptct loop verify --state results/loop/state
python scripts/dyadic_runs.py run draper_fixed 64 --canaries     # in a checkout
python scripts/projectq_dyadic_runs.py run qft 39 --canaries     # in a checkout, needs ProjectQ
```
