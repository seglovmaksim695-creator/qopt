#!/usr/bin/env python
"""D2 evidence for qoptct.dyadic: random circuits with rotations pi*p/2^k against the dense simulator.

    python scripts/dyadic_fuzz.py [--cases 3000] [--spec-cases 1000] [--out results/dyadic/fuzz.json]

Part 1: circuit = reference (prove_equal) on random circuits up to 10 qubits, equivalent rewrites and mutants.
Part 2: circuit = specification (prove_spec): QFT variants, Draper adders, constant adders, RGQFT multipliers,
correct and mutated, up to 10 qubits.  Every verdict is compared with the independent dense simulator: a PROVED must
be dense-equal (with the reported global phase), a REFUTED dense-different.  Exit code 1 if any case disagrees.
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from qoptct.dyadic.fuzz import run_fuzz, run_spec_fuzz  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--cases', type=int, default=3000)
    ap.add_argument('--spec-cases', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=8)
    ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                                  'results', 'dyadic', 'fuzz.json'))
    a = ap.parse_args(argv)
    t0 = time.time()
    eq = run_fuzz(a.seed, a.cases, nmax=10, depth=(3, 30))
    t1 = time.time()
    sp = run_spec_fuzz(a.seed + 1, a.spec_cases)
    out = dict(equal_vs_dense=eq, spec_vs_dense=sp, time_equal=round(t1 - t0, 1), time_spec=round(time.time() - t1, 1))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, 'w') as f:
        json.dump(out, f, indent=1)
    for k, st in (('equal', eq), ('spec', sp)):
        print(k, {x: st[x] for x in ('cases', 'PROVED', 'REFUTED', 'UNKNOWN', 'equal', 'different')},
              'disagreements:', len(st['bad']))
    return 1 if eq['bad'] or sp['bad'] else 0


if __name__ == '__main__':
    sys.exit(main())
