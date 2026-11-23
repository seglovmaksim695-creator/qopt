"""Command line: prove | lemmas | canaries | mutants.   python -m qoptct.foreign.qualtran.kaliski <command> --help"""
import argparse
import json
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m qoptct.foreign.qualtran.kaliski', description=__doc__)
    sub = ap.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('prove', help='prove KaliskiModInverse(n, p) and print the verdict')
    a.add_argument('n', type=int)
    a.add_argument('p', type=int)
    a.add_argument('--out', default=None, help='write verdict + evidence to this JSON file')
    a.add_argument('--lemmas', default=None, help='lemma result file (default: computed once, cached in $QOPT_CACHE)')
    a.add_argument('--canaries', default='auto', help='canary evidence file (required at n=256)')
    a.add_argument('--run-inputs', type=int, default=256)
    a.add_argument('--sim-inputs', type=int, default=127)
    b = sub.add_parser('lemmas', help='run the Z3 lemmas and write the result')
    b.add_argument('--out', default=None)
    b.add_argument('--timeout', type=int, default=120)
    c = sub.add_parser('canaries', help='run the canaries (damaged circuits that must NOT be proved); n=256 takes about 45 minutes')
    c.add_argument('--n', type=int, required=True)
    c.add_argument('--p', type=int, required=True)
    c.add_argument('--out', required=True)
    c.add_argument('--only', default='')
    c.add_argument('--random', type=int, default=12)
    d = sub.add_parser('mutants', help='wiring mutants: all must be rejected with the reason "wire flow"')
    d.add_argument('n', type=int)
    d.add_argument('p', type=int)
    d.add_argument('--out', default=None)
    d.add_argument('--noprepass', action='store_true')
    d.add_argument('--aud', action='store_true')
    args = ap.parse_args(argv)
    if args.cmd == 'prove':
        from .verdict import prove_kaliski
        v = prove_kaliski(args.n, args.p, lemmas=args.lemmas, canaries=args.canaries, run_inputs=args.run_inputs,
                          sim_inputs=args.sim_inputs, out=args.out, verbose=True)
        print(json.dumps({k: v[k] for k in ('status', 'n', 'p', 'checks', 'program', 'sec')}, indent=1, default=str))
        return 0 if v['status'] == 'PROVED' else 1
    if args.cmd == 'lemmas':
        from . import lemmas_z3
        res = lemmas_z3.ensure_lemmas(path=args.out, timeout=args.timeout, refresh=True, verbose=True)
        print('%d of %d implications unsat; controls ok %d of %d; ok=%s; %s' % (
            res['proved'], res['implications'], res['controls_ok'], res['controls'], res['ok'], res['z3']))
        return 0 if res['ok'] else 1
    if args.cmd == 'canaries':
        from . import canary
        argv2 = ['--n', str(args.n), '--p', str(args.p), '--out', args.out, '--random', str(args.random)] + (['--only', args.only] if args.only else [])
        canary.main(argv2)
        return 0
    from . import mutants
    res = mutants.run_mutants(args.n, args.p, noprepass=args.noprepass, aud=args.aud, outp=args.out)
    print('rejected %d of %d; reason "wire flow" for all: %s' % (sum(r['status'] == 'NOT_PROVED' for r in res['rows']), len(res['rows']), res['all_flow_reason']))
    return 0 if res['all_rejected'] else 1


if __name__ == '__main__':
    sys.exit(main())
