"""``python -m qoptct.foreign.qualtran.montgomery --n 256 --w 1 [--original] [--mutant H1] --out results/qualtran/montgomery_n256_w1.json``"""
import argparse
import json
import sys

from . import prove_montgomery


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m qoptct.foreign.qualtran.montgomery',
                                 description='Prove the Montgomery multiplier of Qualtran (heavy sizes: run under lab/tools/one).')
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--w', type=int, required=True)
    ap.add_argument('--p', type=int, default=0, help='modulus (default: 2^255-19 for n=256, else the largest prime <= 2^(n-1))')
    ap.add_argument('--original', action='store_true', help='use the installed Qualtran without the patches')
    ap.add_argument('--mutant', default='', help='run a deliberately damaged program (H1 H2 M1 M1b M2 M3 M4 M5)')
    ap.add_argument('--fresh', action='store_true', help='ignore cached lowering')
    ap.add_argument('--budget', type=int, default=3600)
    ap.add_argument('--max-mb', type=int, default=2800)
    ap.add_argument('--out', default='')
    a = ap.parse_args(argv)
    v = prove_montgomery(a.n, a.w, a.p or None, patched=not a.original, mutant=a.mutant or None, budget=a.budget, max_mb=a.max_mb,
                         out=a.out or None, fresh=a.fresh)
    print(json.dumps({k: v.get(k) for k in ('verdict', 'n', 'w', 'patched', 'mutant', 'reasons', 'seconds')}, indent=1, default=str))
    return 0 if v['verdict'] == 'PROVED' else 1


if __name__ == '__main__':
    sys.exit(main())
