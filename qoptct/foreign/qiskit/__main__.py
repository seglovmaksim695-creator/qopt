# SPDX-License-Identifier: Apache-2.0
"""python -m qoptct.foreign.qiskit prove FAMILY N [--canaries] [--budget S]   |   list   |   fractional

Reproduces every row of the registry (qoptct/foreign/qiskit/data/registry.json) and docs/QISKIT.md.
"""
import argparse
import json
import sys

from .library import FAMILIES, prove_case


def run_canaries(family, n, budget=120, count=None):
    from .check import check_qiskit
    from .mutants import make_mutants
    case = FAMILIES[family](n)
    base = prove_case(case, budget_sec=budget)
    clean = tuple(r for r in case.clean if r in {q.name for q in case.circuit.qregs})
    rare = next((q for q in case.circuit.qregs if q.name not in clean and q.name != 'ctrl' and len(q) >= 3), None)
    wires = [case.circuit.find_bit(q).index for q in rare] if rare is not None else None
    rows = []
    for label, mut in make_mutants(case.circuit, rare_wires=wires, count=count):
        r = check_qiskit(mut, case.spec, clean=clean, budget_sec=budget)
        row = dict(mutant=label, status=r['status'], method=r['method'], time=r['time'])
        if r['status'] == 'PROVED':       # equivalent mutant? settle it with the reference miter
            q = check_qiskit(mut, reference=case.circuit, clean=clean, ref_clean=clean, budget_sec=budget)
            row['equivalent_to_original'] = q['status'] == 'PROVED'
        rows.append(row)
    return base, rows


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m qoptct.foreign.qiskit')
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('prove')
    p.add_argument('family', choices=sorted(FAMILIES))
    p.add_argument('n', type=int)
    p.add_argument('--canaries', action='store_true')
    p.add_argument('--budget', type=float, default=300)
    sub.add_parser('list')
    a = ap.parse_args(argv)
    if a.cmd == 'list':
        print('\n'.join(sorted(FAMILIES)))
        return 0
    if a.canaries:
        base, rows = run_canaries(a.family, a.n, a.budget)
        bad = [r for r in rows if r['status'] != 'REFUTED' and not r.get('equivalent_to_original')]
        print(json.dumps(dict(base=dict(status=base['status'], time=base['time']), canaries=rows, not_refuted=len(bad)),
                         indent=1))
        return 0 if base['status'] == 'PROVED' and not bad else 1
    r = prove_case(FAMILIES[a.family](a.n), budget_sec=a.budget)
    print(json.dumps({k: r[k] for k in ('case', 'n', 'status', 'method', 'detail', 'time') if k in r}, indent=1))
    return 0 if r['status'] == 'PROVED' else 1


if __name__ == '__main__':
    sys.exit(main())
