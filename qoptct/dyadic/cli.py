# SPDX-License-Identifier: Apache-2.0
"""Command line of qoptct.dyadic (``qoptct dyadic ...``).

    qoptct dyadic prove CIRCUIT --reference REF [--exact-phase] [--budget S] [--json]
    qoptct dyadic prove CIRCUIT --qft [--no-swaps] [--inverse] [--wires 0,1,2] [--budget S] [--json]
    qoptct dyadic results [--name NAME] [--long] [--json]

CIRCUIT and REF are OpenQASM 2/3 files (rotations pi*p/2^k, float or written as pi/2^k) or, for REF, a .qc file.
Exit code: 0 only when the verdict is PROVED; 1 REFUTED; 2 UNKNOWN (including unsupported input); 3 internal error;
4 usage error.  The verdict is up to a global phase unless --exact-phase is given.  Nothing is sent anywhere.
"""
from __future__ import annotations

import argparse
import json
import sys

EXIT = {'PROVED': 0, 'REFUTED': 1, 'UNKNOWN': 2, 'ERROR': 4, 'INTERNAL': 3}


def _parser():
    ap = argparse.ArgumentParser(prog='qoptct dyadic', description=__doc__.split('\n')[0])
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('prove', help='prove a circuit against a reference circuit or the QFT')
    p.add_argument('circuit')
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--reference')
    g.add_argument('--qft', action='store_true', help='the specification is the QFT on --wires (default: all wires)')
    p.add_argument('--no-swaps', action='store_true', help='QFT without the final reversal of the wires')
    p.add_argument('--inverse', action='store_true')
    p.add_argument('--wires', help='comma separated wire numbers of the QFT (default: all)')
    p.add_argument('--exact-phase', action='store_true', help='also require that the global phase is 0')
    p.add_argument('--budget', type=float, default=120.0)
    p.add_argument('--json', action='store_true')
    r = sp.add_parser('results', help='the registry of proved results')
    r.add_argument('--name')
    r.add_argument('--long', action='store_true')
    r.add_argument('--json', action='store_true')
    return ap


def main(argv=None):
    a = _parser().parse_args(argv)
    if a.cmd == 'results':
        from . import registry
        es = registry.entries(a.name)
        if a.json:
            print(json.dumps(es, indent=1))
        elif a.long:
            print('\n\n'.join(registry.format_entry(e) for e in es))
        else:
            print(registry.format_table(es))
        return 0 if es else 4
    from . import NotDyadic, QFTSpec, Unsupported, parse_qasm, prove_qasm
    from .api import prove_spec
    phase = 'exact' if a.exact_phase else 'global'
    try:
        if a.reference:
            r = prove_qasm(a.circuit, a.reference, phase=phase, budget_sec=a.budget)
        else:
            c = parse_qasm(open(a.circuit, encoding='utf8').read() if '\n' not in a.circuit else a.circuit)
            wires = [int(x) for x in a.wires.split(',')] if a.wires else list(range(c.n))
            r = prove_spec(c, QFTSpec(wires, swaps=not a.no_swaps, inverse=a.inverse), phase=phase, budget_sec=a.budget)
    except (NotDyadic, Unsupported) as e:           # an angle that is not pi*p/2^k, a gate that is not supported
        print('UNKNOWN  (input)\nreason: %s' % e)
        return 2
    except (OSError, ValueError, TypeError, KeyError, IndexError) as e:
        print('error: %s: %s' % (type(e).__name__, e), file=sys.stderr)
        return 4
    v = r.get('verdict', 'ERROR')
    if a.json:
        print(json.dumps(r, indent=1, default=str))
    else:
        print('%s  (%s)' % (v, r.get('method')))
        if r.get('reason'):
            print('reason: %s' % r['reason'])
        for x in r.get('assumptions') or []:
            print('assumption: %s' % x)
        if r.get('counterexample'):
            print('counterexample: %s' % json.dumps(r['counterexample']))
    return EXIT.get(v, 3)


if __name__ == '__main__':
    sys.exit(main())
