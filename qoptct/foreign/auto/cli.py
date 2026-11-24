# SPDX-License-Identifier: Apache-2.0
"""Command line of ``qoptct.foreign.auto``.

    python -m qoptct.foreign.auto prove TARGET --n N [--p P] [--w W] [--cv CV] [--out FILE]
    python -m qoptct.foreign.auto selfcheck [--out FILE]
    python -m qoptct.foreign.auto canaries TARGET --n N [--p P] [--w W] [--out FILE] [--only NAME,...]

TARGET: montgomery (DirtyOutOfPlaceMontgomeryModMul, fixed Qualtran kind 'c'), modadd, cmodadd, modsub, cmodsub,
modneg, cmodneg, moddbl (installed Qualtran 0.7.0).  ``prove`` runs in a child process whose Qualtran is the right
variant (``qualtran._variant``); the exit code is 0 only for PROVED (1 not proved, 4 usage error or Qualtran missing).
Numbers are read without eval (``2**255-19`` is accepted).  Heavy sizes (n = 256): run under lab/tools/one.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

__all__ = ['main', 'prove_target', 'variant_of']


def variant_of(target):
    return 'c' if target == 'montgomery' else None


def _num(s):
    from ..qualtran._intexpr import int_expr
    return int_expr(s)


def prove_target(target, n, p=None, *, w=1, cv=1, out=None, timeout=None, e2e_cases=16, extra=()):
    """Prove a named target in a child process with the right Qualtran; returns the verdict dict."""
    from ..qualtran._variant import run_json
    args = ['prove', target, '--n', str(n), '--w', str(w), '--cv', str(cv), '--inprocess',
            '--e2e-cases', str(e2e_cases)] + (['--p', str(p)] if p else []) + list(extra)
    res = run_json(variant_of(target), 'qoptct.foreign.auto', args, timeout=timeout)
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text(json.dumps(res, indent=1, default=str))
    return res


def _print_summary(v):
    print(f"{v.get('target')} {v.get('params')}: {v.get('verdict', v.get('status'))}")
    for r in v.get('reasons', [])[:5]:
        print('  reason:', str(r)[:300])
    prog = v.get('program') or {}
    if prog:
        print(f"  program: {prog.get('nq')} wires, {prog.get('ops')} operations, sha256 {prog.get('ops_sha256')}")
    for pc in (v.get('pieces') or [])[:40]:
        print(f"  {pc['status']:10s} {pc['mode']:7s} {pc['type']:32s} x{pc['count']:<4d} ops {pc['ops']:<8} "
              f"{(pc.get('qasm_sha256') or '')[:16]}  {pc.get('contract') or ''}")
    e2e = v.get('end_to_end')
    if e2e:
        print(f'  end-to-end simulation: {e2e}')


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m qoptct.foreign.auto', description=__doc__.split('\n\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    pp = sub.add_parser('prove', help='prove a named target')
    pp.add_argument('target')
    pp.add_argument('--n', required=True)
    pp.add_argument('--p', default=None)
    pp.add_argument('--w', default='1')
    pp.add_argument('--cv', default='1')
    pp.add_argument('--out', default=None)
    pp.add_argument('--e2e-cases', default='16')
    pp.add_argument('--inprocess', action='store_true', help='run here (the Qualtran on the path must be the right one)')
    pp.add_argument('-v', '--verbose', action='store_true')
    sc = sub.add_parser('selfcheck', help='check every registry contract against call_classically (small instances)')
    sc.add_argument('--out', default=None)
    cn = sub.add_parser('canaries', help='damaged copies of a target that must NOT be proved')
    cn.add_argument('target')
    cn.add_argument('--n', required=True)
    cn.add_argument('--p', default=None)
    cn.add_argument('--w', default='1')
    cn.add_argument('--out', default=None)
    cn.add_argument('--only', default=None)
    cn.add_argument('--inprocess', action='store_true')
    a = ap.parse_args(argv)
    try:
        if a.cmd == 'selfcheck':
            from .registry import selfcheck
            rows = selfcheck()
            ok = all(r['ok'] for r in rows)
            for r in rows:
                print('ok ' if r['ok'] else 'BAD', r['checked'], r['contract'], '|', r['bloq'][:80])
            if a.out:
                Path(a.out).write_text(json.dumps(dict(ok=ok, rows=rows), indent=1, default=str))
            return 0 if ok else 1
        n = _num(a.n)
        p = _num(a.p) if a.p else None
        w = _num(a.w)
        if a.cmd == 'canaries':
            from .canary import run_canaries, canaries_child
            if not a.inprocess:
                res = canaries_child(a.target, n, p, w=w, only=a.only, out=a.out)
            else:
                res = run_canaries(a.target, n, p, w=w, only=(a.only.split(',') if a.only else None), out=a.out)
            print(json.dumps({k: res.get(k) for k in ('target', 'params', 'total', 'not_proved', 'proved', 'ok')},
                             default=str))
            return 0 if res.get('ok') else 1
        cv = _num(a.cv)
        if not a.inprocess:
            v = prove_target(a.target, n, p, w=w, cv=cv, out=a.out, e2e_cases=_num(a.e2e_cases))
            _print_summary(v)
            return 0 if v.get('verdict') == 'PROVED' else 1
        from .api import build_target, prove_bloq
        bloq, kind, params = build_target(a.target, n, p, w=w, cv=cv)
        log = (lambda m: print('  ', m, flush=True)) if a.verbose else None
        v = prove_bloq(bloq, target=a.target, params=params, expect_variant=kind or 'original', log=log,
                       e2e_cases=_num(a.e2e_cases), out=a.out)
        if not a.out:
            _print_summary(v)
        return 0 if v.get('verdict') == 'PROVED' else 1
    except (ValueError, ImportError) as e:
        print(f'error: {type(e).__name__}: {e}', file=sys.stderr)
        return 4


if __name__ == '__main__':
    sys.exit(main())
