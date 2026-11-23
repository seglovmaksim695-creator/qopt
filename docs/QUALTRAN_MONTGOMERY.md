# Montgomery multiplication of Qualtran (`qoptct.foreign.qualtran.montgomery`)

Proof client for `qualtran.bloqs.mod_arithmetic.DirtyOutOfPlaceMontgomeryModMul` (bit size `n`, window `w`, odd modulus
`p` with `2p <= 2^n`).  It needs the optional extra `qoptct[qualtran]` (Qualtran 0.7.0, z3-solver).

## API

```python
from qoptct.foreign.qualtran.montgomery import prove_montgomery, exhaustive_small, prime_for

v = prove_montgomery(16, 2)                  # p = largest prime <= 2^15; for n = 256: p = 2^255 - 19
v = prove_montgomery(256, 1, out="results/qualtran/montgomery_n256_w1.json")
v = prove_montgomery(16, 1, patched=False)   # the ORIGINAL installed Qualtran: NOT_PROVED, with a witness
v['verdict']                                 # 'PROVED' or 'NOT_PROVED' (nothing else; UNKNOWN is NOT_PROVED)
```

Command line (heavy sizes: run under `lab/tools/one`, one process, 3 GB):

```
python -m qoptct.foreign.qualtran.montgomery --n 256 --w 1 --out results/qualtran/montgomery_n256_w1.json
python -m qoptct.foreign.qualtran.montgomery --n 16 --w 1 --original        # installed Qualtran without the patches
python -m qoptct.foreign.qualtran.montgomery --n 16 --w 1 --mutant H2       # a damaged program: must be NOT_PROVED
```

`prove_montgomery(n, w, p=None, *, patched=True, mutant=None, budget=3600, max_mb=2800, seed=7, out=None, fresh=False)`.
The stages (lowering with all structural checks, invariants, proofs) run in separate processes so that the memory of
the lowering is freed; work files are in `$QOPT_CACHE/montgomery/<key>/` (`run.log` is the log).

`patched=True` (default) runs with a COPY of Qualtran carrying the patches the proofs were made with
(`patches.patched_qualtran('c')`: PR 1997 -- the final reduction subtracts `p` instead of adding it -- and exact integer
tables, no int64 overflow above 63 bits).  The installed Qualtran is never modified.  The verdict records the version
and the SHA-256 of every patched file.

## What is proved

For all `y < p` and all n-bit `x`, the lowered program is the exact Montgomery recurrence with window `w`: it outputs `x`, `y`,
`qrom_indices` (the digits `m_i`), `reduced = [T_N >= p]` and `target = T_N - p*reduced`; every other wire returns to
zero; no measurement outcome leaves a phase.  With the telescoping identity `2^n T_N = x y + p sum c_i 2^(w i)` and the
invariant `T_i < y + p` this is `x*y*2^-n mod p`, below `p`.

The multiplier is cut at the boundaries of the Qualtran nodes: `X(reduced)` + window 0, windows 1..N-1,
`LessThanConstant`, controlled `AddK`.  Each piece goes through the unchanged prover `kprove` (BDD path sum) and the
exact partial-domain guard, and is cross-checked by the independent simulator.  The generic window is proved by its
sub-nodes (`CAdd` x w, `Xor`, `QROAMClean + Add + QROAMCleanAdjoint`) and their composition (BDD; for `w >= 4`
by wire positions, interval arithmetic and exact integer runs).  Windows 2..N-1 are accepted only if their program
(QASM and layout hash) equals that of window 1.

Checks that make the sub-node route sound and that are never skipped (a failure gives NOT_PROVED):

* wire glue between the windows and the program registers; live-wire bookkeeping (every input of a segment is alive,
  every non-input wire is fresh);
* coverage: every operation of every window belongs to exactly one sub-node and the children of the window tile its span
  with no gaps or overlaps;
* freshness: the wires that a sub-node assumes to be zero are not inputs of the window and were not touched before;
* the operation `X(reduced)` before window 0 is exactly one `X` on the `reduced` wire and window 0 does not touch it;
* the end-to-end simulation of the whole program against the exact recurrence (n <= 64);
* interval invariants and the telescoping identity;
* z3 must be importable (its version is recorded); a missing z3 never gives PROVED.

## Regression (all must be NOT_PROVED)

`montgomery.mutants`: H1, H2 (holes found by the audit) and M1, M1b, M2, M3, M4 (internal wires, extra operations,
operations counted twice) in every generic window, M5 (wrong wire in `X(reduced)`).  `montgomery.regression.canaries`:
single-operation mutants of the real program (window 0, the sub-nodes, `LessThanConstant`, the controlled `AddK`) that
the independent simulator separates from the contract; `kprove` with a negative control and the domain guard must not
accept them.  The original Qualtran (`patched=False`) is NOT_PROVED against the fixed contract with a witness in
`Bad = {T_N >= p}` (about a quarter of all inputs); the error is only in the final `+p` block.

```
python -m pytest tests/foreign_qualtran/test_montgomery_pure.py tests/foreign_qualtran/test_montgomery.py -m "not slow"
python -m pytest tests/foreign_qualtran/test_montgomery.py -m slow      # n = 64 and n = 256; run under lab/tools/one
```

## Assumptions

Qualtran 0.7.0 with the patches above; `p` odd, `2p <= 2^n`; the rules of the bridge (A1 controlled-And,
C1 QROAMCleanAdjoint, see `QUALTRAN_BRIDGE.md`) and the prover `kprove` with the exact domain guard; the interval
invariant of the generic window is analytic (not BDD) for `n > 192` and for `w >= 4`; for `w >= 4` the composition of
the sub-nodes is by wire positions + interval arithmetic + exact integer runs, not by BDD.  The pair "forward + inverse"
is not covered here.
