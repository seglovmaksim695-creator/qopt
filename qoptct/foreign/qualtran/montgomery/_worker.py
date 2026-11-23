"""Entry point of one stage in its own process: ``python -m qoptct.foreign.qualtran.montgomery._worker <stage> ...``"""
import argparse
import sys

from . import _stages as S


def main(argv=None):
    ap = argparse.ArgumentParser(prog='montgomery._worker')
    sub = ap.add_subparsers(dest='stage', required=True)
    for name in ('lower', 'inv', 'prove'):
        sp = sub.add_parser(name)
        sp.add_argument('--n', type=int, required=True)
        sp.add_argument('--w', type=int, required=True)
        sp.add_argument('--p', type=int, required=True)
        sp.add_argument('--out', required=True)
        if name == 'lower':
            sp.add_argument('--kind', default='')
            sp.add_argument('--mutant', default='')
            sp.add_argument('--e2e-max-n', type=int, default=64)
            sp.add_argument('--seed', type=int, default=1)
        if name == 'inv':
            sp.add_argument('--part', choices=('first', 'generic'), required=True)
        if name == 'prove':
            sp.add_argument('--inp', required=True)
            sp.add_argument('--budget', type=int, default=3600)
            sp.add_argument('--max-mb', type=int, default=2800)
            sp.add_argument('--seed', type=int, default=7)
    a = ap.parse_args(argv)
    if a.stage == 'lower':
        S.stage_lower(a.n, a.w, a.p, a.kind or None, a.out, mutant=a.mutant or None, e2e_max_n=a.e2e_max_n, seed=a.seed)
    elif a.stage == 'inv':
        S.stage_inv(a.n, a.w, a.p, a.part, a.out)
    else:
        S.stage_prove(a.inp, a.out, a.budget, a.max_mb, a.seed)
    return 0


if __name__ == '__main__':
    sys.exit(main())
