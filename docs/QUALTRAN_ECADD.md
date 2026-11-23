# Point addition and "forward + backward" pairs for Qualtran (`qoptct.foreign.qualtran.ecadd`, `.pairs`)

Optional (`pip install "qoptct[qualtran]"`: Qualtran 0.7.0, z3-solver).  Everything here is positive-only: a
verdict is `PROVED` or it is not, and "skipped" or "unknown" is never a proof.

## What is proved

| API | claim | where it is proved |
|---|---|---|
| `ecadd.prove_ecadd(n, p)` | the fixed `ECAdd(n, p, 1)` is the point addition of the Montgomery form (all cases: generic, doubling, `P = -Q`, `P = O`, `Q = O`), `a`, `b`, `lam_r` unchanged, every other wire clean, for all inputs with every word `< p` | six steps by the ROBDD prover (`ecadd.core`, qoptct's own `bdd`/`lift`/gadget lemmas) + glue in Z3 |
| `ecadd.prove_ecadd(n, p, patched=False)` | the original Qualtran 0.7.0 scheme is the point addition outside the groups G1-G3 of finding 05 (status `PROVED_OUTSIDE_G1_G3`); that the scheme leaves a dirty flag ON the groups was shown in round 3 (finding 05), not by this package | same, glue mode `E3` |
| `ecadd.prove_ecwindowaddr(n, p, window)` | `ECWindowAddR` (fixed): `(c, Q) -> (c, Q + c*R)`, lookup and its removal gate by gate (rule C1), inner `ECAdd` as a box with the contract above | `ecadd.e4` |
| `pairs.check_pair_kaliski(n, p)` | `KaliskiModInverse` followed by its `uncompute=True` twin is the identity on all wires (exact phase, all measurement outcomes) | gate mirror (`pairs.mirror`, `pair_compose`) |
| `pairs.check_pair_montgomery(n, w, p)` | the same for `DirtyOutOfPlaceMontgomeryModMul` (fixed Qualtran, kind `'c'`) | gate mirror with `LD`/`CLR` abstraction of the table blocks (`mont_compose`) |

The steps use the multiplier and the Kaliski inverse as **boxes with contracts**
(`ecadd.contracts.CONTRACTS`: statement, the registry result that proves it, and the (n, p) instances for which
that result exists).  `contracts_status(n, p)` is part of every verdict: for an instance without a proved
contract the point-addition proof is conditional on it, and the verdict says so.  For n=16, p=32749 and for n=256,
p=2^255-19 every contract (multiplier, inverse, and the pair of each) is proved by a registry result.

### Status of a point-addition verdict

* `PROVED` (exit code 0): the six steps and their glue are proved AND every box contract is proved for exactly this (n, p, w=1):
  by a registry result, or, with `--prove-boxes` (API: `prove_boxes=True`), by a proof made in the same run.
* `PROVED_CONDITIONAL` (exit code 3): the steps and the glue are proved, but at least one box contract is only assumed for this
  (n, p).  The verdict lists them (`conditional_on`, e.g. `MM forward`, `KINV pair`).  It is NOT a proof of the point addition and
  is never stored in the registry.  `ecwindowaddr` is conditional the same way on the result for `ECAdd` for the same (n, p).
* `NOT_PROVED` (exit code 1): anything else, including a modulus outside the domain (below).

### Domain of p

Checked before any work (`qoptct.foreign.qualtran.domain`); outside it the status is `NOT_PROVED` with the reason, for `prove` and
`pair` alike.  Montgomery (and its pair): p odd, 2p <= 2^n, p < 2^n - 1 (primality is not needed).  Kaliski: p an odd prime
(certificate, or a named standard prime) with p < 2^n - 1.  Pair of the inverse: p odd with p < 2^n - 1 (the mirror does not use
primality).  `ecadd`, `ecwindowaddr`: all of these together (an odd prime p with 2p <= 2^n), because the glue divides modulo p.

## Qualtran variants

The fixed scheme is Qualtran 0.7.0 plus the patches of `qoptct.foreign.qualtran.patches` (flag fix E05, PR 1997,
exact-integer tables).  The proof clients run in a child process whose `PYTHONPATH` starts with a COPY of Qualtran
in `$QOPT_CACHE` (the installed one is never touched).  Every child reports which variant it really imported
(SHA-256 of the four touched files, `_variant.detect_variant`); a result from another variant is `NOT_PROVED`.

| client | variant |
|---|---|
| `prove_ecadd(patched=True)` | `e` |
| `prove_ecadd(patched=False)`, `check_pair_kaliski` | `original` |
| `prove_ecwindowaddr` | `e4` |
| `check_pair_montgomery` | `c` |

The wiring between the six steps is read from the real `ECAdd.decompose_bloq().connections` on every call and must be
equal to the reference wiring shipped in `ecadd/data/wiring_ecadd.json` (checked for both variants at n=16 and 256);
the glue refuses anything else.

Z3: the `z3` binary (`$QOPT_Z3`, `PATH`) via SMT-LIB, else the `z3-solver` package; the version is in the result.
Without Z3 no verdict is `PROVED`.

## Reproduce

```
python -m qoptct.foreign.qualtran.ecadd ecadd  --n 16  --p 32749 --out results/qualtran/ecadd_n16.json
python -m qoptct.foreign.qualtran.ecadd ecadd  --n 256 --p '2**255-19' --out results/qualtran/ecadd_n256.json
python -m qoptct.foreign.qualtran.ecadd window --n 256 --p '2**255-19' --window 1 --out results/qualtran/ecadd_window_n256_w1.json
python -m qoptct.foreign.qualtran.pairs kaliski    --n 256 --p '2**255-19' --out results/qualtran/pairs_kaliski_n256.json
python -m qoptct.foreign.qualtran.pairs montgomery --n 256 --w 1 --p '2**255-19' --out results/qualtran/pairs_montgomery_n256.json
```

(the pair commands expect the right Qualtran variant on `PYTHONPATH`; the Python API sets it up, or use
`patches.run_with_patched_qualtran('c', [...])`).  n=256 takes about an hour and 3 GB for `ecadd` (step 3 is the
longest), a few minutes for `window` and the pairs.

## Tests

`tests/foreign_qualtran/test_ecadd_prove.py`, `test_ecadd_regression.py`, `test_pairs.py`.  Fast tests (about 1.5 min) cover
n=16 proofs and all regression sets: corrupted glue contracts and wiring, missing lemmas, original steps 1/2 in the
fixed glue, step canaries (n=16), ECWindowAddR canaries, the pair canaries (n=8, 16), and the independent
mutants of the audit (box interfaces, XOR properties, mirror matching, LD/CLR abstraction, brute force of the eight
glue lemmas on small curves).  Every one of them must come out NOT proved.  `-m slow` adds the n=256 reruns and
the large canary samples.

## Assumptions and limits

* Boxes: Montgomery multiplier (window 1) and Kaliski inverse are taken with contracts (see above), not re-proved here.
* The three records of the step formulas (`specs` for the ROBDD prover, `zspecs` for Z3, `ispecs` for the simulator)
  agree on sample points (n=16: 3000 per step, n=256: 600 per step); equality is not proved formally.
* The lemmas about the curve (Z3 glue) were checked by brute force over all curves and point pairs for p <= 23 and
  by sympy; they are not proved for p = 2^255-19.
* `ECWindowAddR`: window sizes 1, 2, 4 (`QROAMClean` with K=1); curve `y^2 = x^3 + 2x + 3`, a fixed point R.
* The pair check is a structural gate-level argument: it relies on the gadget lemma for AND/erase of qoptct and on the
  forward block leaving its internal wires clean (contract of the forward block).
