#!/usr/bin/env python
"""Registry runs of the ProjectQ -> dyadic bridge (exact angles; needs ProjectQ).

    python scripts/projectq_dyadic_runs.py run FAMILY N [--canaries] [--budget S]     FAMILY: qft | addconst
    python scripts/projectq_dyadic_runs.py registry [--out results/dyadic]            all sizes, canaries at the largest

Every size is recorded from the real ProjectQ program (AutoReplacer with ProjectQ's own decomposition rules, no
optimizer) by ``qoptct.foreign.projectq.dyadic.record_dyadic`` and proved by ``qoptct.dyadic``.  Angles are accepted only
when they are exactly what ProjectQ stores for p*pi/2^k (k <= 38); QFT on 40 qubits or more (R(pi/2^39) is stored as
6e-12, R(pi/2^43) as 0) is UNKNOWN.  At the largest proved size: 20 or more wrong variants (scripts/dyadic_runs.py
``canaries``) that must all be REFUTED, and a numeric product-state check of the recorded circuit.
Output: results/dyadic/pq_exact_<family>.json (same row format as the Qiskit families).
"""
import argparse
import json
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import dyadic_runs as dr                                                      # noqa: E402
from qoptct.dyadic import QFTSpec, const_adder_spec, prodsim                  # noqa: E402
from qoptct.foreign.projectq.dyadic import check_projectq_dyadic, record_dyadic   # noqa: E402

SIZES = {'qft': [2, 4, 8, 16, 32, 38, 39, 40, 43], 'addconst': [2, 4, 8, 16, 32, 38, 40]}
CANARY_AT = {'qft': 39, 'addconst': 38}


def constant(n):
    return (0xB5C3A1F7D94E62B8 * 0x9E3779B97F4A7C15 + n) % (1 << n) | 1       # mixed bits, odd


def program(family, n):
    def prog(eng):
        q = eng.allocate_qureg(n)
        if family == 'qft':
            from projectq.ops import QFT
            QFT | q
        else:
            from projectq.libs.math import AddConstant
            AddConstant(constant(n)) | q
        return dict(registers=dict(q=q))
    return prog


def spec_of(family, n):
    if family == 'qft':
        return lambda regs: QFTSpec(regs['q'], swaps=False)      # ProjectQ's QFT has no final reversal
    return lambda regs: const_adder_spec(constant(n), regs['q'])


def run_one(family, n, do_canaries, budget, seed=1):
    t0 = time.time()
    r = check_projectq_dyadic(program(family, n), spec_of(family, n), budget_sec=budget)
    d = r.get('dyadic') or {}
    res = dict(family='pq_exact_' + family, n=n, qubits=(r.get('circuit') or {}).get('wires'),
               gates=(r.get('circuit') or {}).get('ops'), verdict=r['status'], method=d.get('method'),
               reason=r['detail'] if r['status'] != 'PROVED' else None, assumptions=r['assumptions'],
               angles=r.get('angles'), projectq=(r.get('circuit') or {}).get('projectq'),
               global_phase=d.get('global_phase'), dense_check=d.get('dense_check'), time=round(time.time() - t0, 2),
               command='python scripts/projectq_dyadic_runs.py run %s %d --canaries' % (family, n))
    if r['status'] == 'PROVED' and do_canaries:
        c, _ = record_dyadic(program(family, n))
        spec = spec_of(family, n)({'q': list(range(n))})
        rng = random.Random(seed + n)
        cres = []
        for name, e, ph in dr.canaries(c, spec, family, n, rng):
            t2 = time.time()
            rr = dr.prove(e, spec, min(budget, 120.0), phase=ph)
            cres.append(dict(name=name, verdict=rr['verdict'], method=rr.get('method'), time=round(time.time() - t2, 2)))
        res['canaries'] = cres
        res['canaries_refuted'] = sum(1 for x in cres if x['verdict'] == 'REFUTED')
        res['canaries_proved'] = sum(1 for x in cres if x['verdict'] == 'PROVED')
        chk = prodsim.check_qft if family == 'qft' else prodsim.check_arith
        t, bad, ent = chk(c, spec, 20, random.Random(seed))
        res['product_check'] = dict(inputs=t, bad=bad, entangling=ent)
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    p1 = sub.add_parser('run')
    p1.add_argument('family', choices=sorted(SIZES))
    p1.add_argument('n', type=int)
    p1.add_argument('--canaries', action='store_true')
    p1.add_argument('--budget', type=float, default=300.0)
    p2 = sub.add_parser('registry')
    p2.add_argument('--out', default=os.path.join(os.path.dirname(HERE), 'results', 'dyadic'))
    p2.add_argument('--budget', type=float, default=300.0)
    a = ap.parse_args(argv)
    if a.cmd == 'run':
        print(json.dumps(run_one(a.family, a.n, a.canaries, a.budget), indent=1, default=str))
        return 0
    for fam, sizes in SIZES.items():
        rows = []
        for n in sizes:
            row = run_one(fam, n, n == CANARY_AT[fam], a.budget)
            print(fam, n, row['verdict'], row.get('canaries_refuted'), flush=True)
            rows.append(row)
        with open(os.path.join(a.out, 'pq_exact_%s.json' % fam), 'w') as f:
            json.dump(dict(family='pq_exact_' + fam, rows=rows, python=sys.version.split()[0]), f, indent=1, default=str)
    return 0


if __name__ == '__main__':
    sys.exit(main())
