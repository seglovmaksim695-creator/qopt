# qoptct loop: propose, prove, accumulate

`qoptct loop` closes the circle of the optimizer and the provers. Proposers suggest improved circuits for a target,
the provers of `qoptct.verify` (unchanged) check each candidate against the target's reference `.qc`, and only proved
improvements are written to an append-only journal. The "best known" table is built from the journal alone.

## Commands

    qoptct loop run    --state DIR --target SPEC [--target SPEC ...] [--hours H] [--proposers lookup,optimize]
                       [--canaries 20] [--budget SEC] [--max-steps N] [--prove-all]
    qoptct loop verify --state DIR [--no-reprove] [--budget SEC]     # exit 0: journal holds; 1: problems listed
    qoptct loop table  --state DIR
    qoptct loop status --state DIR

Target specs: `crypto:<name>` (a circuit of `results/crypto`), `lookup:c6w16s1` and `lookupadd:c5w12s1` (random
table, xor lookup or lookup-add), `qc:<path>[:clean]` (any `.qc`; the first proved candidate becomes its baseline).
Long runs: start through the lab limiter (`lab/tools/one`); stop at any time, the next `run` continues.

## How a step works

1. A proposer yields the next parameter set (deterministic order; sets already in `attempts.jsonl` are skipped).
2. A child process builds the candidate and computes its cost vector (AND = T/4, T, CX = cx + cz, qubits used).
   The child is killed on a time limit; a crash is an attempt with outcome `error`, never an acceptance.
3. A candidate that cannot enter the Pareto front (not better in any of AND, T, CX, qubits) is `no_gain` and is not
   proved (`--prove-all` proves it anyway, for statistics).
4. Otherwise it is proved against the stored reference: `qoptct.verify.kprove`, then `qoptct.verify.prove`.
   Only PROVED passes; SAMPLED, UNKNOWN and ERROR never do.
5. Independent check: the sampled simulators of `qoptct.verify.quick_check` must not disagree.
6. Canaries: at least 20 mutants of the candidate (phase, erasure correction, seam, flip) must be REFUTED by the same
   prover on the same reference, at least half of them in the classes phase / erasure / seam. A mutant that the prover
   calls PROVED but the simulators refute is a false PROVED: the candidate is rejected and the case is written to
   `HOLES/` of the state directory.
7. The file is hashed before and after the proof; the journal record is written only if the hashes agree.

## Journal record

`journal.jsonl`, one JSON object per line: `seq`, `prev`, `hash` (SHA-256 chain over the canonical record),
`role` (`baseline` | `improvement`), `kind` (`baseline` | `dominates` | `tradeoff`), `target`, `ref_sha`, `cand_sha`,
`cost_before`, `cost_after`, `prover` (name, version = qoptct version + hash of `qoptct/verify/*.py`, status, method,
time), `sim`, `canaries`, `proposer` (parameters), `repro` (command; the files are in `store/`).
The baseline of every target is proved and recorded like any candidate. `dominates` means better than the best known
circuit by the order (AND, CX, qubits, T) without being worse in any key; `tradeoff` is Pareto-incomparable.

## Guard against false accumulation

`qoptct loop verify` believes nothing: it checks the hash chain and the replay rules (sequence, first record of a
target is its baseline, reference hash pinned by the baseline, `cost_before` equals the best known cost at that point,
`kind` matches the costs, status PROVED), the stored files against their SHA-256, the costs recomputed from the stored
circuit, the proof itself again on the stored files plus the simulators, and that `best.json` / `BEST.md` equal the
table derived from the journal. `tests/test_loop.py` covers 18 kinds of corruption (replaced or edited candidate,
wrong hashes, swapped reference, deleted / reordered / duplicated / truncated records, forged costs and status, proof
of another circuit, hand-edited table, ...) and a prover that proves everything.
Limits: the proofs are those of `qoptct.verify` and carry its assumptions; the state directory is not signed, so an
attacker who also reruns every proof and rewrites every hash is not stopped by the chain but only by the re-proof,
which is why `verify` re-proves.

## Proposers

* `lookup`: parameter search for table lookups (`qoptct.crypto.lookup`): seeds, splits and tries of the bilinear plan
  search, Gray-code fan-out batch sizes, QROAM with one select bit, and the same for lookup-add.
* `optimize`: the optimizer portfolio, `qoptct optimize` on the reference with several modes under a time limit.

New proposers: add `(params(target), build(target, params, workdir))` to `PROPOSERS` in `qoptct/loop/proposers.py`.

## Not covered

Qualtran blocks (Montgomery, Kaliski, ECAdd) are proved compositionally from Qualtran bloqs; they are not given as a
candidate `.qasm` plus a reference `.qc`, so the loop does not propose for them.
