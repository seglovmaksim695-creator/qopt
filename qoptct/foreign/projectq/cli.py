"""``python -m qoptct.foreign.projectq``: proofs about ProjectQ programs.

    python -m qoptct.foreign.projectq prove multiplier --m 163 --source /path/to/Binary_ECC [--canaries 24] [--out FILE.json]
    python -m qoptct.foreign.projectq prove point-addition --m 8 --source /path/to/Binary_ECC [--canaries 24] [--out FILE.json]
    python -m qoptct.foreign.projectq prove library --block add_quantum --n 4 [--canaries 24] [--out FILE.json]
    python -m qoptct.foreign.projectq results [--name NAME] [--size N] [--long] [--json]

Exit code 0 only when the verdict is PROVED and all canaries behaved (``results``: the registry was printed).  1: REFUTED or a
canary was PROVED, 2: UNKNOWN, 4: usage error or a missing optional dependency / source.  Nothing is sent anywhere.
The Binary_ECC sources (Eclipse Public License 2.0) are not part of this package: ``--source`` or $QOPT_BINARY_ECC.
"""
from __future__ import annotations

import argparse
import json
import sys

__all__ = ['main']


def _parser():
    ap = argparse.ArgumentParser(prog='qoptct projectq', description=__doc__.split('\n')[0])
    sp = ap.add_subparsers(dest='cmd', required=True)
    pr = sp.add_parser('prove', help='prove one block and run its canaries')
    ps = pr.add_subparsers(dest='what', required=True)
    for name in ('multiplier', 'point-addition'):
        p = ps.add_parser(name)
        p.add_argument('--m', type=int, required=True, help='field width')
        p.add_argument('--source', help='Binary_ECC folder (default: $QOPT_BINARY_ECC)')
    ps.choices['point-addition'].add_argument('--variant', default='FLT_in', choices=['FLT_in', 'FLT_out'])
    p = ps.add_parser('library')
    p.add_argument('--block', required=True)
    p.add_argument('--n', type=int, required=True)
    for p in ps.choices.values():
        p.add_argument('--canaries', type=int, default=24)
        p.add_argument('--seed', type=int, default=1)
        p.add_argument('--out', help='write the evidence as JSON')
    rs = sp.add_parser('results', help='the registry of proved results')
    rs.add_argument('--name')
    rs.add_argument('--size', type=int)
    rs.add_argument('--long', action='store_true')
    rs.add_argument('--json', action='store_true')
    return ap


def main(argv=None):
    args = _parser().parse_args(argv)
    if args.cmd == 'results':
        from . import registry
        es = registry.entries(args.name, args.size)
        if args.json:
            print(json.dumps(es, indent=1))
        elif args.long:
            print('\n\n'.join(registry.format_entry(e) for e in es))
        else:
            print(registry.format_table(es))
        return 0 if es else 1
    from ._record import ProjectQUnavailable, need_projectq
    try:
        need_projectq()
        from . import results
        if args.what == 'multiplier':
            r = results.prove_multiplier(args.m, args.source, args.canaries, args.seed)
        elif args.what == 'point-addition':
            r = results.prove_point_addition(args.m, args.source, args.canaries, args.seed, args.variant)
        else:
            r = results.prove_library(args.block, args.n, args.canaries, args.seed)
    except (ProjectQUnavailable, FileNotFoundError, KeyError) as e:
        print(f'error: {e}', file=sys.stderr)
        return 4
    if args.out:
        with open(args.out, 'w') as f:
            json.dump(r, f, indent=1, default=str)
    c = r.get('canaries')
    print(f"{r['id']}: {r['status']}  {r['detail']}")
    if c:
        print(f"canaries: {c['not_proved']} of {c['decided']} not proved ({c['equivalent']} equivalent), groups {c['by_group']}"
              + ('' if c['ok'] else '  CANARY CRITERION NOT MET'))
    if r['status'] == 'REFUTED':
        return 1
    if r['status'] != 'PROVED':
        return 2
    return 0 if c and c['ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
