# Generic hierarchical prover for foreign circuits (`qoptct.foreign.auto`, new in 0.8)

The clients of 0.7 (`qoptct.foreign.qualtran.kaliski`, `.montgomery`, `.ecadd`, `.pairs`) each carry code and lemmas
written for one Qualtran block.  `qoptct.foreign.auto` is one mechanism for all blocks: it walks down the
decomposition of a block, proves the small pieces directly and every composite node from the contracts of its
children, with no code for the block itself.  What a block needs is **data**: the specification of the root (a
formula), contracts of the primitives it is built from (a registry, each entry checked against Qualtran's classical
simulation), and, for a loop, its invariant (a formula).  A false `PROVED` is worse than none, so every verdict is
positive-only (`PROVED` or `NOT_PROVED`; a time-out, a missing contract or a refused bloq is never a proof), carries
an independent end-to-end simulation, and every result of the registry comes with canaries (damaged copies that must
not be proved).

```
pip install "qoptct[qualtran]"                                   # qualtran==0.7.0 and z3-solver, as for 0.7
python -m qoptct.foreign.auto prove cmodadd --n 256              # PROVED (seconds)
python -m qoptct.foreign.auto prove montgomery --n 256 --w 1     # PROVED (5 minutes, 1.8 GB); same program hash as 0.7
python -m qoptct.foreign.auto canaries montgomery --n 256        # 30+ damaged copies, none proved (about 40 minutes)
python -m qoptct.foreign.auto results [--long]                   # the registry of results
```

## Results

All at n = 256, p = 2^255 - 19, Qualtran 0.7.0 (variant `original` unless stated: `c` is the Montgomery fix of
PR 1997 with exact-integer tables, `e` the fixed `ECAdd`, both the patched copies of 0.7).  The hash is `ops_sha256` of
the whole lowered program (the same definition as in 0.7).  `python -m qoptct.foreign.auto results` prints the current
table from `qoptct/foreign/auto/data/registry.json`; this one is a copy.

| block | n | parameter | status | criterion | canaries: not proved / total (joints) | operations | program hash | = 0.7 |
|---|---|---|---|---|---|---|---|---|
| DirtyOutOfPlaceMontgomeryModMul (Qualtran kind c) | 256 | w = 1 | PROVED | met | 34 / 34 (20) | 1443063 | `90e665cd297649d8` | yes |
| DirtyOutOfPlaceMontgomeryModMul (Qualtran kind c) | 256 | w = 2 | PROVED | met | 34 / 34 (20) | 1249271 | `ffd71399723c972a` | yes |
| CModAdd | 256 | - | PROVED | met | 30 / 30 (16) | 11272 | `fc5c3e00362e71f3` | - |
| CModAddK | 256 | k = 3^101 mod p | PROVED | met | 30 / 30 (16) | 11434 | `d583358058a742f0` | - |
| CModNeg | 256 | - | PROVED | met | 29 / 29 (15) | 6128 | `1e80b93370558c73` | - |
| CModSub | 256 | - | PROVED | met | 31 / 31 (17) | 16870 | `1bd7c64477350139` | - |
| _ECAddStepOne of ECAdd (Qualtran kind e) | 256 | - | PROVED | met | 31 / 31 (17) | 37839 | `ded17ef5f27f898a` | yes |
| ModAdd | 256 | - | PROVED | met | 30 / 30 (16) | 9207 | `7ace01b8437bcb57` | - |
| ModDbl | 256 | - | PROVED | met | 29 / 29 (15) | 4619 | `16a33bced9a4b3bb` | - |
| ModNeg | 256 | - | PROVED | met | 28 / 28 (14) | 6126 | `140ed1a588a031e5` | - |
| ModSub | 256 | - | PROVED | met | 31 / 31 (17) | 14805 | `6bd686ff22d6b653` | - |
| Negate | 256 | - | PROVED | met | 26 / 26 (12) | 2295 | `4af597993328eede` | - |
| Subtract | 256 | - | PROVED | met | 27 / 27 (13) | 2805 | `e3f31459c0d3acfc` | - |
| SubtractFrom | 256 | - | PROVED | met | 26 / 26 (12) | 2549 | `654ec5ddec0cb107` | - |

"= 0.7": the program that is proved is, operation by operation, the one the 0.7 client proved (Montgomery: the program
hash of the 0.7 registry entry; the QASM hashes of the window, of `CAdd`, `Xor` and `LessThanConstant` also occur in
the 0.7 result file; ECAdd step 1: both lowerings compared one operation after the other by `compare`).
Apart from Montgomery, none of these blocks was proved on its own before (several of the modular blocks occur inside
the ECAdd steps, which 0.7 proved as whole programs).

Not proved: `CtrlScaleModAdd` and `CModMulK` are proved at n = 8 but not at n = 16 (Z3 cannot follow a chain of n different modular additions; their loop has iterations of three
children with different constants, which `Fold` does not cover); the ECAdd steps 2 to 6 and `ECWindowAddR` need the
multiplier and the inverse as boxes with relational contracts (the uncompute box returns the junk of its forward box),
which the contract language does not have; pairs (forward + uncompute) are not checked by the mirror of `pairs`.

## How a block is proved

```
bloq --(adapter: lower)--> Program (operations + tree of nodes) --(Engine.prove(root, spec))--> verdict
                                                                 \--(end_to_end: simulation)--/
```

1. **Lowering** (`qualtran_adapter.lower`).  The strict bridge of 0.7 (`lower_bloqgraph` with the real connections of
   `decompose_bloq` at every level, fresh wires for every allocation, the atom rules A1 controlled-And, C1
   `QROAMCleanAdjoint`, basis states, and `TwoBitCSwap` by the opt-in Fredkin rule; then `translate(light=True)`)
   produces one flat program over four operations `x`, `cx`, `and` (into a fresh zero wire) and `erase` (the
   measurement-based uncompute of an `and`, with its phase correction).  The hash of this list is the program hash.
   The lowering records the calls of the root; a call inside a call is lowered alone and aligned with the program
   operation by operation (an injective renaming of wires that must reproduce every operation), so memory stays at the
   size of the program (1.8 GB for the 1.44 million operations of Montgomery n = 256).
2. **Tree** (`tree`).  Every call is a node with its span of operations, its input and output wires per register
   (packing: a scalar register is an unsigned integer, least significant wire first; a bit array of shape (m,) has its
   element 0 as the most significant bit, as Qualtran; other shapes are separate elements `name[i]`), its bloq and a
   key.  The core (`expr`, `contract`, `tree`, `compose`, `leaf`) does not import Qualtran.
3. **Contracts** (`contract`, `expr`).  A contract has input and output registers with widths, a precondition `pre`,
   definitions `defs`, and one expression per output register.  Expressions: integers, `+ - * div mod`, comparisons,
   `and or not ite`, bitwise `xor and or` (expanded bit by bit up to 64 bits, or rewritten: XOR with all ones is the
   complement), `bit`, `slice`, shifts, table lookup `lut`, sums.  Validation rejects unknown names and outputs that do
   not fit their width on the domain.  Kinds: `primitive` (proved on its gates), `spec` (proved from its children),
   `body` (loop body), `derived` (computed by the engine).
4. **Leaves** (`leaf`).  A node whose contract is a primitive is proved on its own operations by the unchanged qoptct
   prover (`bridge.strict_proof.prove_program`: kprove ROBDD path sum, every input and every measurement outcome, with
   the exact AND/ERASE domain guard), and gets an independent simulation receipt (24 inputs of its domain, every wire,
   phases of every outcome).  A leaf with a precondition whose domain has at most 256 points (enumerated exactly by
   Z3: the precondition with every found point excluded is unsat) is proved by simulating every point with the phase
   of every outcome instead (the table uncompute `QROAMCleanAdjoint`, whose precondition `target = data[sel]` makes
   the BDD of the phase explode); so is a node without inputs or without outputs (states, effects).
5. **Composition** (`compose.Engine`).  A composite node is proved by applying the contracts of its children in program
   order to a symbolic state of the wires (every wire is 0, 1 or bit j of a symbolic value) and asking Z3 (linear
   integer arithmetic with `div`/`mod` by constants) the obligations below.  A child without a registry contract is
   proved recursively and gets a **derived** contract (its outputs are what its children compute; the preconditions
   of its children are passed to the parent as obligations).  Identical children (same QASM, same layout, same
   contract) are proved once and reused; the 256 windows of the Montgomery multiplier cost one proof.
6. **Loops** (`contract.Fold`).  A run of children of one type is accepted only with a fold in the parent's
   specification: body contract, carried registers, initial values, the inputs of every iteration and an invariant,
   all as formulas.  The engine proves the body once, checks that every iteration really reads the stated inputs, and
   proves the invariant at iteration 0 and by one step lemma for a symbolic iteration.  Without a fold the children are
   composed one by one; if Z3 cannot follow them, the verdict is `NOT_PROVED` (never `PROVED` by default).  The test
   `test_loop_without_invariant_is_not_proved` removes the invariant of the Montgomery fold: not proved.
7. **End to end** (`api.end_to_end`).  The whole lowered program is simulated on 16 inputs of the domain (edge values,
   Z3 models, random) and compared with the specification on every wire, phases included.  A mismatch is a reason
   against `PROVED` and is reported with the input.
8. **Specification check** (`api.spec_selfcheck`).  The specification of a named target is compared with
   `call_classically` of the same Qualtran on every input of a small instance (n = 4, Montgomery n = 6); a mismatch is
   a reason against `PROVED`.

## The checks at the joints

The audits of round 3 found the holes of composite proofs at the joints: wires that do not flow where the proof
assumes, operations that belong to no proved piece.  Every composition here checks, and every check has canaries:

| check | rule | canaries |
|---|---|---|
| coverage | the spans of the children tile the span of the parent without gaps or overlaps; an operation owned by no child, or by a child that does not contain it contiguously, refuses the node | J4 (X between two children), J4b (before the first), J4c (after the last), J6 (between grandchildren), J7 (last operation of a child attributed to the next one) |
| wire flow | a child reads its inputs from the wires the program really connects (the events of the lowering), never from the register names; outputs are written to its real output wires | J1, J1b (two wires of a child exchanged), J2 (two output wires exchanged), J3 (wires exchanged between two children), J13 (an iteration of a loop reads its carried register from other wires) |
| freshness | every wire a child takes as an ancilla or allocation must be 0 in the symbolic state when the child starts; freed wires must be provably 0; at the end every non-output wire is 0 | J5 (X inside a child), gate canaries G5, G7 |
| contracts of children | a leaf is proved on its own operations against its contract; a wrong contract cannot be used | J8 (an output +1), J8b (a child that claims to change nothing), J12 (wrong loop step) |
| domain | every child is applied only where its precondition provably holds (Z3 obligation with the facts of the state); the parameter domain of a specification is checked before anything is computed (`SpecDomainError`: p odd, p < 2^n, 2p <= 2^n where needed, window divides n, 0 <= k < p ...) | J9 (root precondition weakened), J10 (child domain narrowed to x < 2), J11 (wrong loop invariant) |
| gates | the leaf provers see every operation | G1 (phase correction of an erase removed), G2 (moved), G3 (carry CX reversed), G4 (AND of a high carry weakened: wrong only on rare inputs), G5 (X removed), G6 (two non-commuting CX exchanged), G7 (last operation removed), G8 (random single-operation changes) |

A canary that changes the program is used only if the change is real: an exchange of wires is chosen where the
child's contract gives a different function, a gate change only where the damaged leaf differs from its contract in
simulation, a moved boundary only between children of a node proved by composition, a weakened precondition first
where the whole program then differs from the specification on the added inputs; the damaged leaves are spread over
all leaf types (at most four of one type).  Every gate canary is counted only
with a witness (an input on which the damaged program differs, found by simulation).  The canary set of a result fails
if any canary is proved.

## Registry of primitive contracts and specifications

`registry.contract_for(bloq)` gives the contract of a Qualtran bloq or `None`:

* gates and states: `XGate`, `CNOT`, `Toffoli`, `And` (all control values, and its uncompute with precondition),
  `MultiAnd`, `MultiControlX`, `ZeroState/OneState/ZeroEffect/OneEffect`, `IntState`, `IntEffect` (precondition
  value = k), `Swap`, `CSwap`, `TwoBitSwap`, `TwoBitCSwap`, `MultiTargetCNOT`, `OnEach(n, X)`;
* arithmetic: `Add` (unsigned, Montgomery and two's complement words: the same bits), `AddK`, `CAdd`,
  `OutOfPlaceAdder`, `LessThanConstant`, `EqualsAConstant`, `LessThanEqual`, `GreaterThan`, `LinearDepthGreaterThan`,
  `CLinearDepthGreaterThan`, `Equals`, `Xor`, `XorK`, `BitwiseNot`;
* table lookup: `QROAMClean` (one target, no junk) and its adjoint (rule C1, precondition target = data[sel]);
* wrappers: `Controlled`/`ControlledViaAnd` of a primitive (any control value, multi-bit controls), `Adjoint` of a
  primitive (self-inverse ones keep their contract; additions get the subtraction; a block that computes into fresh
  registers gets the uncompute with the precondition "the registers hold the computed value"), `Always`;
* specifications (`specs`): `ModAdd`, `CModAdd` (domain x <= p, y < p: at y = p the circuit leaves an ancilla dirty),
  `ModSub`, `CModSub`, `ModNeg`, `CModNeg`, `ModDbl`, `CModAddK`, `CtrlScaleModAdd`, `CModMulK`, `Subtract`,
  `SubtractFrom`, `Negate`, `DirtyOutOfPlaceMontgomeryModMul` (recurrence of the windows with fold and invariant
  T < y + p), step 1 of `ECAdd`.

`python -m qoptct.foreign.auto selfcheck` compares every entry with `call_classically` on every input of small
instances (88 instances, all agree; `tests/foreign_auto/test_auto_qualtran.py` runs it).

## Results registry and reproduction

`results/auto/*.json` are the result files: verdicts (`schema qoptct.auto/1`: target, parameters, versions, Qualtran
variant with the SHA-256 of the patched files, SHA-256 of the sources the verdict depends on (this package without
the canary generator, the results registry and the cross-check; the bridge; `qoptct/verify`; `montgomery/_core.py`),
specification and its SHA-256, program hash, every proved piece with its QASM hash and how it was proved, statistics
of the composition, end-to-end simulation, specification check, assumptions), canaries (`qoptct.auto.canaries/1`:
every canary with its mutation, verdict and witness, the SHA-256 of the verdict sources and of `canary.py`) and the
cross-check with 0.7 (`compare`).
`scripts/build_auto_registry.py` builds `qoptct/foreign/auto/data/registry.json` from them (`--check` compares);
`tests/foreign_auto/test_auto_registry.py` checks that the registry is a fresh build and that the files have the
recorded SHA-256.  The criterion of a result is shown separately from its status: `PROVED`, end-to-end simulation
without a mismatch, at least 20 counted canaries with at least 10 at the joints, none proved, made by the same sources.

Each entry has its command (`reproduce`); the full set was made by

```
python -m qoptct.foreign.auto prove montgomery --n 256 --w 1 --out results/auto/montgomery_n256_w1.json   # kind c, 5 min
python -m qoptct.foreign.auto prove montgomery --n 256 --w 2 --out results/auto/montgomery_n256_w2.json
python -m qoptct.foreign.auto prove TARGET --n 256 --out results/auto/TARGET_n256.json      # modadd ... ecadd_step1
python -m qoptct.foreign.auto compare ecadd_step1 --n 256 --out results/auto/compare_ecadd_step1_n256.json
python -m qoptct.foreign.auto canaries TARGET --n 256 [--w W] --out results/auto/canaries_TARGET_n256[_wW].json
python scripts/build_auto_registry.py
```

`prove`, `canaries` and `compare` start a child process with the right Qualtran (`qualtran._variant`, the patched
copies live in `$QOPT_CACHE`).  Exit code 0 only for `PROVED` (1 not proved, 4 usage error or Qualtran missing).
Tests: `pytest tests/foreign_auto` (fast, under a minute); `-m slow` reruns n = 256.

## Trusted base and assumptions

* the strict bridge of 0.7 and its atom rules (A1, C1 with its assumptions A1-A4, basis states, Fredkin), as in 0.7;
* the qoptct prover kprove with the exact AND/ERASE domain guard and the gadget lemma for AND/erase, as in 0.7;
* Z3 for the composition (its answer `unsat` is trusted; `unknown` is never a proof);
* the code of `compose` (the rules above) and of `leaf` (the bridge between a leaf and kprove);
* the specification of the root is the statement.  It is compared with `call_classically` on small instances, not
  proved equal to Qualtran's intent.  For Montgomery the statement is the recurrence of the windows; that its result
  equals x*y*2^-n mod p is the telescoping identity of Montgomery multiplication (as in 0.7, not re-proved here);
* the registry contracts are checked against `call_classically` only on small instances; at n = 256 a primitive is
  proved against its contract on its own gates, so a wrong contract makes a proof fail, not succeed.

## Other libraries (Qiskit, ProjectQ)

The core only needs a `tree.Program`: a list of operations in the four-operation language, the number of wires, and an
adapter that answers `node(path)`, `children(node)`, `contract(node)` and `bloq_of(node)`.  A node has its span of
operations, its registers with their wires, an opaque object and a key (equal keys = same block).  To plug in another
library:

* **Qiskit**: the hierarchy is the instruction tree.  `QuantumCircuit.data` gives `CircuitInstruction`s; an
  instruction with a `definition` is a composite node whose children are the instructions of its definition, mapped to
  the parent's qubits by position; a gate without definition is an atom.  Registers: the `qargs` of the instruction,
  grouped as the library groups them (`QuantumRegister`s of the definition, or one register per argument).  Lowering:
  every atom must become `x`, `cx`, `and` or `erase` (a Toffoli onto a fresh ancilla is `and`; its measurement-based
  uncompute or a Toffoli onto the same ancilla with the same controls is `erase`; other gates are refused).  The key is
  the instruction name with its parameters; contracts of the primitives (for example `qiskit.circuit.library`
  adders and comparators) go into a registry like `registry.py` and are checked against a classical simulation of the
  definition on small widths.
* **ProjectQ**: record the commands of a `MainEngine` with a backend that keeps the decomposition levels
  (`DecompositionRule`s applied by the `AutoReplacer`), and the meta blocks `Compute`/`Uncompute`, `Control`, `Loop`,
  `Dagger` as nodes; a `Loop` is a natural fold.  The rest as for Qiskit.

The engine, the leaf prover, the canaries of the joints and the end-to-end check stay as they are.

## Modules

| module | role |
|---|---|
| `expr` | expressions: validation, evaluation, bounds, Z3 terms, BDD words |
| `contract` | `Contract`, `Fold`, validation and evaluation |
| `tree` | `Program`, `Node` (library independent) |
| `leaf` | proof of a leaf (kprove, small domains, simulation receipts) |
| `compose` | `Engine`: recursive descent, composition in Z3, folds, derived contracts |
| `qualtran_adapter` | Qualtran bloq -> `Program` (strict bridge, alignment of nested calls) |
| `registry`, `specs` | contracts of Qualtran primitives and specifications of composite blocks |
| `api`, `cli` | `prove_bloq`, named targets, end-to-end check, specification check, command line |
| `canary` | damaged copies at the joints and in the gates |
| `results`, `crosscheck` | registry of results, comparison with the 0.7 clients |
