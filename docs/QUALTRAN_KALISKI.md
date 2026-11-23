# Kaliski inverse of Qualtran: proof client (`qoptct.foreign.qualtran.kaliski`)

Proves, gate by gate, that Qualtran's `KaliskiModInverse(bitsize=n, mod=p)` (Qualtran 0.7.0, unpatched) is correct for
**every** input, with the unchanged qoptct provers and the strict bridge (`docs/QUALTRAN_BRIDGE.md`).
Needs `pip install "qoptct[qualtran]"` and Z3 (the `z3` binary, or the `z3-solver` package).

```python
from qoptct.foreign.qualtran.kaliski import prove_kaliski

v = prove_kaliski(16, 32749)          # n = register width, p = odd prime, p < 2^n - 1
v['status']                           # 'PROVED' or 'NOT PROVED (<failed checks>)'
v['checks']                           # every item of the criterion: chain, run, flow, sim, lemmas, prime, ...
v['program']['top_ops_sha256']        # hash of the lowered program the proof is about
```

Command line: `python -m qoptct.foreign.qualtran.kaliski prove N P [--out FILE]`, also `lemmas`, `canaries`, `mutants`.

## What PROVED means

For every `0 <= x < p` the circuit outputs `x^-1 * 2^(2n) mod p` (`0` for `x = 0`) and the junk register `M*2^(2n) + T`
(`T = 2^(2n-1-K)`, one bit, `K` = the iteration that sees `v = 0`; `M` = the step-type flags `m_k`, `k <= K`); every other wire is zero at the end; phase and measurement flow are those of the bridge rules.
All of this holds **under the mathematical assumptions M-carry, M-dbl, M-a, M-s, M-end, M-end-vf**, which are closed
for all `n, p, x` by 243 Z3 implications (`lemmas_z3`, 25 controls that must fail), **provided p is prime**.
Not covered: the translation (bridge + composition of the pieces) is the trusted part; it is checked against the full expansion
for small `n` (`lower.check_against_full`) and by an independent bit-parallel simulation.

## The checks of the verdict (`verdict.judge`)

All must hold: `chain` (Walker proof of 2n iterations of 6 steps, nothing admitted), `assumptions_covered`, `run` (contracts executed
on 256 sample inputs equal Qualtran's `call_classically` and the integer recursion `pq`), `flow` (**wire flow**: every step reads
exactly the wires the previous step wrote, new registers sit on fresh wires, junk and output wire lists are in the documented
order; run in the proof and in the execution, with the expected number of checks), `client_hash`, `sim` (127 inputs),
`sources_unchanged`, `lemmas` (the Z3 file is for this lemma source and the full query list, hashes included), `prime`
(trial division, deterministic Miller-Rabin, Lucas-Lehmer or Pocklington; a named standard prime without a certificate
is labelled `KNOWN PRIME`; a probabilistic test alone is refused), `samples`, and at `n = 256` also `canaries`.
No Z3 means no `lemmas` and therefore never PROVED. The Z3 version is written into the verdict.

## Canaries and mutants (regression)

* `mutants.run_mutants(n, p)`: 14 wiring mutants of the corrected client (+ the auditor's F1-F4 with `aud=True`) built on the
  real circuit; all must be rejected with the reason "wire flow".  At `n = 256` F1-F4 use the parameters of the audit.
* `canary.py`: damaged copies of the real circuit (one operation replaced); a canary counts when the prover does not prove it and a
  witness input is found.  `python -m qoptct.foreign.qualtran.kaliski canaries --n 256 --p P --out FILE` (about 45 minutes) makes the
  evidence that `prove_kaliski(256, ...)` needs (`canaries=FILE`; default: `$QOPT_CACHE/kaliski/canaries_n256.json`).
* The lemma result is cached in `$QOPT_CACHE/kaliski` (about 150 s the first time); a copy for the current source is in
  `results/qualtran/kaliski_lemmas.json`.

## Reproduce

```
pytest tests/foreign_qualtran/test_kaliski_*.py               # seconds to a minute: small n, mutants, negative tests, primes
pytest tests/foreign_qualtran/test_kaliski_*.py -m slow       # n = 256 (needs the canary file), minutes
python -m qoptct.foreign.qualtran.kaliski prove 16 32749 --lemmas results/qualtran/kaliski_lemmas.json
```

Results of the re-run through this package are in `results/qualtran/kaliski_*.json` (verdict + evidence, with the hash of the program);
they agree with the audited round (same program hashes, operation counts, segments, BDD node counts, flow counts, statuses).
