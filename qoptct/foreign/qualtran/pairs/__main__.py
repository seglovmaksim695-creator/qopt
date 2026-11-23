"""Command line of the pair checks (this is what the child processes of the API run)::

    python -m qoptct.foreign.qualtran.pairs kaliski    --n 256 --p '2**255-19' [--full] [--out FILE]
    python -m qoptct.foreign.qualtran.pairs montgomery --n 256 --w 1 --p '2**255-19' [--full] [--out FILE]
    (public commands: they call the API, which starts the right Qualtran variant; exit code 0 only for PROVED)
    python -m qoptct.foreign.qualtran.pairs canary-kaliski    --n 16 --p 65521 [--seed 1] [--only 0,3,7] [--out FILE]
    python -m qoptct.foreign.qualtran.pairs canary-montgomery --n 16 --w 1 --p 32749 [--seed 1] [--only ...] [--out FILE]

The canary commands, and the internal ``run-kaliski`` / ``run-montgomery`` that the API starts, run on the Qualtran that is
imported: Kaliski on the installed 0.7.0, Montgomery on the patched copy (kind ``'c'``, see
``qoptct.foreign.qualtran.patches.run_with_patched_qualtran``).
"""
import argparse
import json
import resource
import sys
import time


def _p(s):
    return int(eval(s, {'__builtins__': {}}, {}))


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m qoptct.foreign.qualtran.pairs')
    ap.add_argument('cmd', choices=['kaliski', 'montgomery', 'run-kaliski', 'run-montgomery', 'canary-kaliski', 'canary-montgomery'])
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--w', type=int, default=1)
    ap.add_argument('--p', required=True)
    ap.add_argument('--full', action='store_true', help='also compare the assembly with the full lowering (small n)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--only', default=None, help='canaries: only these job indexes (comma separated)')
    ap.add_argument('--out', default=None)
    a = ap.parse_args(argv)
    n, p = a.n, _p(a.p)
    if a.cmd in ('kaliski', 'montgomery'):
        from . import check_pair_kaliski, check_pair_montgomery
        res = (check_pair_kaliski(n, p, compare_full=a.full or None) if a.cmd == 'kaliski'
               else check_pair_montgomery(n, a.w, p, compare_full=a.full or None))
        if a.out:
            with open(a.out, 'w') as f:
                json.dump(res, f, ensure_ascii=False, indent=1, default=str)
        print(json.dumps({k: res.get(k) for k in ('api', 'n', 'status', 'mirror_ok', 'sec')}, default=str))
        return 0 if res['status'] == 'PROVED' else 1
    from qoptct.foreign.qualtran._variant import detect_variant
    var = detect_variant()
    t0 = time.time()
    res = dict(cmd=a.cmd, qualtran_kind=var['kind'], qualtran_file=var['qualtran_file'], files_sha256_16=var['files_sha256_16'])
    try:
        if a.cmd == 'run-kaliski':
            from .pair_compose import check_pair, CErr
            try:
                res.update(check_pair(n, p, compare_full=a.full))
            except CErr as e:
                res.update(n=n, p=p, mirror_ok=False, error=str(e))
            res['p'] = str(p)
        elif a.cmd == 'run-montgomery':
            from .mont_compose import check_pair_mont
            from .pair_compose import CErr
            try:
                res.update(check_pair_mont(n, a.w, p, compare_full=a.full))
            except CErr as e:
                res.update(n=n, w=a.w, p=p, mirror_ok=False, error=str(e))
            res['p'] = str(p)
        else:
            only = {int(t) for t in a.only.split(',')} if a.only else None
            if a.cmd == 'canary-kaliski':
                from .pair_canary import main as cmain
                r = cmain(n, p, a.seed, None, only)
            else:
                from .mont_canary import main as cmain
                r = cmain(n, a.w, p, a.seed, None, only)
            res.update(r)
            res['p'] = str(p)
            res['status'] = 'OK' if res.get('base_mirror_ok') and not res.get('false_pass') else 'CHECK'
    except Exception as e:  # noqa: BLE001
        import traceback
        res.update(status='ERROR', detail=f'{type(e).__name__}: {e}', tb=traceback.format_exc()[-2500:])
    res['wall_seconds'] = round(time.time() - t0, 2)
    res['peak_rss_mb'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 1)
    if a.out:
        with open(a.out, 'w') as f:
            json.dump(res, f, ensure_ascii=False, indent=1, default=str)
    short = {k: res.get(k) for k in ('cmd', 'n', 'mirror_ok', 'status', 'false_pass', 'wall_seconds')}
    print(json.dumps(short, default=str))
    return 0 if res.get('status') != 'ERROR' else 1


if __name__ == '__main__':
    sys.exit(main())
