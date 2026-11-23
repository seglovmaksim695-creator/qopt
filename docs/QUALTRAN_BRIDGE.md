# Qualtran bridge (`qoptct.foreign.qualtran`)

Qualtran is an optional dependency: `pip install "qoptct[qualtran]"` (Qualtran 0.7.0, z3-solver).
Without it qoptct installs as before and the Qualtran tests are skipped.

## What it does

`qoptct.foreign.qualtran.bridge` lowers a Qualtran bloq **gate by gate** (the actual `decompose_bloq` connections, no
optimizer, no re-synthesis) to the logical language of qoptct (`x`, `cx`, `and`, `erase`) and hands the result to
the **unchanged** qoptct provers (`qoptct.verify.kprove`, BDD).  Every source gate has exactly one explicit rule;
anything else is refused.

| module | role |
|---|---|
| `qualtran_layout` | `lower_bloqgraph`: physical wires of the real graph, registers, `ClassicalContract` |
| `from_qualtran` | `lower_graph` (BloqGraph) and the old Cirq path with the mandatory THRU guard |
| `rules` | `translate(lowered, inputs, outputs)`: leaf -> logical operations, with history |
| `controlled_and` | rule A1: singly controlled And / And-dagger |
| `qroam_clean_adjoint` | rule C1: `QROAMCleanAdjoint` as a primitive with a contract (assumptions A1-A4 in its docstring) |
| `fredkin_rule`, `basis_vectors` | rules for CSwap-like and basis-state atoms |
| `strict_proof` | `prove_program`: full output contract over kprove (no wire excluded, junk is an output) |
| `domain_guard` | exact domains of the partial AND/ERASE gates |
| `sim`, `sim_parse` | independent simulator with the phase of every measurement outcome |
| `gf_measure` | GF2MulMBUC measurement flow |
| `targeted_canaries` | mutants (rare inputs, erase corrections, joints) that must NOT be proved |

## `check_bloq`

```python
from qualtran import QUInt
from qualtran.bloqs.arithmetic import Add
from qoptct.foreign.qualtran import check_bloq

r = check_bloq(Add(QUInt(4)))
r['verdict']   # 'PROVED' | 'REFUTED' (r['witness']) | 'UNKNOWN'
```

The bloq is checked against `bloq.call_classically` on the whole classical domain (all inputs, every measurement
outcome, phase, all non-output wires zero) by exhaustive simulation **and** by the qoptct BDD prover.  `UNKNOWN`
is never a proof: unsupported gate, `call_classically` undefined on part of the domain, too many input bits
(`max_input_bits`, default 12), budget exceeded.  `domain=` states a promise explicitly; the result field
`domain` always says which inputs the claim covers (Qualtran's own domain of e.g. `QMontgomeryUInt` is `1..p-1`).

## Patched Qualtran

The proofs of Montgomery multiplication and point addition used Qualtran 0.7.0 with small fixes.
`qoptct.foreign.qualtran.patches.patched_qualtran(kind)` (`'c'`, `'e'`, `'e4'`) builds a copy of the installed Qualtran in
`$QOPT_CACHE` (default `~/.cache/qopt`) with exactly these files; the installed package is never touched.
Use `patched_env(kind)` / `run_with_patched_qualtran(kind, argv)` to run a fresh interpreter with it.
See the module docstring for what each kind contains; `patches/diff/*.diff` shows the changes.

## Tests

`pytest tests/foreign_qualtran` (seconds).  They skip without Qualtran (`pytest.importorskip`).
