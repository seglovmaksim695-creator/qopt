"""Command line: ``python -m qoptct.foreign.qualtran.ecadd {ecadd,window} --n N --p P [--out FILE]``.

Examples (the registry results ``ecadd_*`` are reproduced by these commands; n=256 needs about 3 GB of memory
and runs for about an hour, so start it through the lab wrapper / with ``ulimit`` as you see fit)::

    python -m qoptct.foreign.qualtran.ecadd ecadd  --n 16  --p 32749 --out results/qualtran/ecadd_n16.json
    python -m qoptct.foreign.qualtran.ecadd ecadd  --n 256 --p '2**255-19' --out results/qualtran/ecadd_n256.json
    python -m qoptct.foreign.qualtran.ecadd window --n 256 --p '2**255-19' --window 1 --out results/qualtran/ecadd_window_n256_w1.json
    python -m qoptct.foreign.qualtran.ecadd ecadd  --n 16  --p 32749 --original     # original 0.7.0: correct outside G1-G3

Exit code 0 only for a PROVED verdict.
"""
import argparse
import json
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m qoptct.foreign.qualtran.ecadd', description=__doc__.split('\n')[0])
    ap.add_argument('cmd', choices=['ecadd', 'window'])
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', required=True, help="modulus, e.g. 32749 or '2**255-19'")
    ap.add_argument('--window', type=int, default=1)
    ap.add_argument('--original', action='store_true', help='ecadd: the original Qualtran 0.7.0 scheme (not the fixed one)')
    ap.add_argument('--steps', default='1,2,3,4,5,6', help='ecadd: subset of steps (a partial run is never PROVED)')
    ap.add_argument('--out', default=None)
    a = ap.parse_args(argv)
    p = eval(a.p, {'__builtins__': {}}, {})
    from . import api
    log = lambda m: print(m, file=sys.stderr, flush=True)
    if a.cmd == 'ecadd':
        res = api.prove_ecadd(a.n, p, patched=not a.original, steps=tuple(int(s) for s in a.steps.split(',')), log=log)
    else:
        res = api.prove_ecwindowaddr(a.n, p, a.window, log=log)
    txt = json.dumps(res, ensure_ascii=False, indent=1, default=str)
    if a.out:
        with open(a.out, 'w') as f:
            f.write(txt + '\n')
    print(json.dumps({k: res[k] for k in ('api', 'n', 'status', 'wall_seconds') if k in res}))
    return 0 if res['status'] in ('PROVED', 'PROVED_OUTSIDE_G1_G3') else 1


if __name__ == '__main__':
    sys.exit(main())
