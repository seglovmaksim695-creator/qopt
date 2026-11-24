"""``qoptct qualtran``: one command line over the Qualtran proof clients.

    qoptct qualtran check "qualtran.bloqs.arithmetic.addition.Add(QUInt(3))"
    qoptct qualtran prove kaliski|montgomery|ecadd|ecwindowaddr --n N [--p P] [--w W | --window W]
    qoptct qualtran pair kaliski|montgomery --n N [--p P] [--w W]
    qoptct qualtran results [--name NAME] [--n N] [--long] [--json]

Exit code 0 only when the verdict is PROVED (``check``: PROVED; ``results``: the registry was printed).  1: not proved or
refuted (the status says which part failed; also a p outside the domain of the proof), 2: UNKNOWN (``check``),
3: PROVED_CONDITIONAL (ecadd / ecwindowaddr: the point addition is proved, but a box contract is only assumed for this (n, p);
``--prove-boxes`` proves the missing ones), 4: usage error or a missing optional dependency.
Heavy sizes (n = 256) take minutes to an hour and up to 3 GB; start them under a memory limit.  Nothing is sent anywhere.
"""
from __future__ import annotations

import argparse
import json
import re
import sys

__all__ = ['main']

P255 = 2 ** 255 - 19
EXIT_OK, EXIT_NOT_PROVED, EXIT_UNKNOWN, EXIT_CONDITIONAL, EXIT_USAGE = 0, 1, 2, 3, 4
DEFAULT_P = {'kaliski': {16: 32749, 256: P255}, 'ecadd': {16: 32749, 256: P255}, 'ecwindowaddr': {16: 32749, 256: P255},
             'pair-kaliski': {16: 32749, 256: P255}}


class UsageError(Exception):
    pass


def _int_expr(text):
    """An integer given as ``2**255-19`` or as a number (no names, no calls)."""
    if len(text or '') > 80 or not re.fullmatch(r'[0-9+\-*() ]+', text or ''):
        raise UsageError(f'cannot read the number {text!r} (digits and + - * ** ( ) only)')
    try:
        from qoptct.foreign.qualtran._intexpr import int_expr
        return int_expr(text)
    except Exception as e:  # noqa: BLE001
        raise UsageError(f'cannot read the number {text!r}: {e}') from e


def _need_qualtran():
    try:
        import qualtran  # noqa: F401
    except ImportError:
        raise UsageError("Qualtran is not installed: pip install 'qoptct[qualtran]'")


def _modulus(a, kind):
    if a.p is not None:
        return _int_expr(a.p)
    if kind == 'montgomery':
        _need_qualtran()
        from .montgomery import prime_for
        return prime_for(a.n)
    p = DEFAULT_P.get(kind, {}).get(a.n)
    if p is None:
        raise UsageError(f'--p is required for n = {a.n} (defaults exist only for n in {sorted(DEFAULT_P.get(kind, {}))})')
    return p


def _write(path, obj):
    if path:
        with open(path, 'w', encoding='utf8') as f:
            json.dump(obj, f, ensure_ascii=False, indent=1, default=str)
            f.write('\n')


def _refuse_domain(kind, n, p, title, as_json, out=None):
    """Outside the domain of the proof: print NOT_PROVED with the reasons (exit 1); ``None`` when the domain holds."""
    from .domain import check_domain
    d = check_domain(kind, n, p)
    if d['ok']:
        return None
    res = dict(status='NOT_PROVED', detail='outside the domain of the proof: ' + '; '.join(d['problems']), domain=d, n=n, p=str(p))
    _write(out, res)
    return _report(title, 'NOT_PROVED', False, ['outside the domain of the proof:'] + ['- ' + x for x in d['problems']], as_json, res)


def _report(title, status, ok, lines, as_json, res, code=None):
    if as_json:
        print(json.dumps(res, ensure_ascii=False, default=str))
    else:
        print(f'{status}  {title}')
        for ln in lines:
            if ln:
                print('  ' + ln)
    if code is not None:
        return code
    return EXIT_OK if ok else EXIT_NOT_PROVED


# ---------------------------------------------------------------- check
def _eval_bloq(expr):
    _need_qualtran()
    import importlib
    import qualtran
    names = {k: getattr(qualtran, k) for k in dir(qualtran) if not k.startswith('_')}
    names['qualtran'] = qualtran
    for dotted in sorted(set(re.findall(r'qualtran(?:\.\w+)+', expr)), key=len):
        parts = dotted.split('.')
        for i in range(len(parts), 1, -1):          # import the longest module prefix that exists
            try:
                importlib.import_module('.'.join(parts[:i]))
                break
            except ImportError:
                continue
    try:
        return eval(expr, {'__builtins__': {'range': range, 'int': int, 'tuple': tuple, 'list': list}}, names)
    except Exception as e:  # noqa: BLE001
        raise UsageError(f'cannot build the bloq from {expr!r}: {type(e).__name__}: {e}') from e


def cmd_check(a):
    bloq = _eval_bloq(a.bloq)
    from .check import check_bloq
    res = check_bloq(bloq, max_input_bits=a.max_input_bits)
    _write(a.out, res)
    v = res['verdict']
    if a.json:
        print(json.dumps(res, ensure_ascii=False, default=str))
    else:
        print(f'{v}  {bloq}')
        for k in ('detail', 'witness', 'inputs_checked', 'domain', 'qasm_sha256'):
            if res.get(k) not in (None, '', [], {}):
                print(f'  {k}: {res[k]}')
    return {'PROVED': EXIT_OK, 'REFUTED': EXIT_NOT_PROVED}.get(v, EXIT_UNKNOWN)


# ---------------------------------------------------------------- prove
def cmd_prove(a):
    _need_qualtran()
    kind = a.block
    n = a.n
    if kind == 'kaliski':
        p = _modulus(a, 'kaliski')
        r = _refuse_domain('kaliski', n, p, f'KaliskiModInverse n={n} p={p}', a.json, a.out)
        if r is not None:
            return r
        from .kaliski import prove_kaliski
        v = prove_kaliski(n, p, lemmas=a.lemmas, canaries=a.canaries, out=a.out, verbose=a.verbose)
        ok = v['status'] == 'PROVED'
        prog = v.get('program') or {}
        return _report(f'KaliskiModInverse n={n} p={p}', v['status'], ok,
                       [f"program hash (top ops): {prog.get('top_ops_sha256')}", f"ops {prog.get('ops')}, wires {prog.get('wires')}",
                        f"checks failed: {[k for k, c in v['checks'].items() if not c]}" if not ok else '',
                        f"Qualtran {v.get('qualtran')}, Z3 {(v.get('lemmas') or {}).get('z3')}"], a.json, v)
    if kind == 'montgomery':
        p = _modulus(a, 'montgomery')
        r = _refuse_domain('montgomery', n, p, f'Montgomery multiplier n={n} w={a.w} p={p}', a.json, a.out)
        if r is not None:
            return r
        from .montgomery import prove_montgomery
        v = prove_montgomery(n, a.w, p, patched=not a.original, out=a.out)
        ok = v['verdict'] == 'PROVED'
        prog = v.get('program') or {}
        return _report(f"Montgomery multiplier n={n} w={a.w} p={p}{' (original Qualtran)' if a.original else ''}", v['verdict'], ok,
                       [f"program hash (ops): {prog.get('ops_sha256')}", f"reasons: {v.get('reasons')}" if not ok else '',
                        f"Qualtran {(v.get('qualtran') or {}).get('version')} patch kind {(v.get('qualtran') or {}).get('patch_kind')}, "
                        f"Z3 {(v.get('versions') or {}).get('z3')}"], a.json, v)
    p = _modulus(a, kind)
    title = (f"ECAdd n={n} p={p}{' (original scheme)' if a.original else ''}" if kind == 'ecadd'
             else f'ECWindowAddR n={n} p={p} window={a.window}')
    r = _refuse_domain(kind, n, p, title, a.json, a.out)
    if r is not None:
        return r
    from . import ecadd
    if kind == 'ecadd':
        v = ecadd.prove_ecadd(n, p, patched=not a.original, prove_boxes=a.prove_boxes)
    else:
        v = ecadd.prove_ecwindowaddr(n, p, a.window, prove_boxes=a.prove_boxes)
    _write(a.out, v)
    # PROVED_OUTSIDE_G1_G3 (original scheme) is a different statement: not PROVED.  PROVED_CONDITIONAL: exit 3.
    ok = v['status'] == 'PROVED'
    assumed = v.get('conditional_on') or []
    lines = [f"claim: {v.get('claim')}",
             ('CONDITIONAL on box contracts that are not proved for this (n, p, w=1): ' + ', '.join(assumed)
              + ' (run with --prove-boxes to prove them here)') if assumed else '',
             f"detail: {v['detail']}" if v.get('detail') and v['status'] == 'NOT_PROVED' else '',
             f"Qualtran variant {v.get('qualtran_kind')}, Z3 {v.get('z3')}"]
    return _report(title, v['status'], ok, lines, a.json, v, code=EXIT_CONDITIONAL if v['status'] == 'PROVED_CONDITIONAL' else None)


def cmd_pair(a):
    _need_qualtran()
    from . import pairs
    n = a.n
    if a.block == 'kaliski':
        p = _modulus(a, 'pair-kaliski')
        r = _refuse_domain('pair-kaliski', n, p, f'pair KaliskiModInverse + uncompute n={n} p={p}', a.json, a.out)
        if r is not None:
            return r
        v = pairs.check_pair_kaliski(n, p)
        title = f'pair KaliskiModInverse + uncompute n={n} p={p}'
    else:
        p = _modulus(a, 'montgomery')
        r = _refuse_domain('pair-montgomery', n, p, f'pair Montgomery multiplier + uncompute n={n} w={a.w} p={p}', a.json, a.out)
        if r is not None:
            return r
        v = pairs.check_pair_montgomery(n, a.w, p)
        title = f'pair Montgomery multiplier + uncompute n={n} w={a.w} p={p}'
    _write(a.out, v)
    ok = v['status'] == 'PROVED'
    return _report(title, v['status'], ok, [f"mirror_ok: {v.get('mirror_ok')}", f"detail: {v.get('detail')}" if v.get('detail') else '',
                                          f"Qualtran variant {v.get('qualtran_kind')}"], a.json, v)


# ---------------------------------------------------------------- results
def cmd_results(a):
    from . import registry
    es = registry.entries(a.name, a.n)
    if not es:
        print('no such result in the registry', file=sys.stderr)
        return EXIT_USAGE
    if a.json:
        print(json.dumps(es, ensure_ascii=False, indent=1))
    elif a.long:
        print('\n\n'.join(registry.format_entry(e) for e in es))
    else:
        print(registry.format_table(es))
        print(f'\n{len(es)} results; Qualtran {registry.load_registry()["qualtran_version"]}; '
              'a result holds for the stated n, p, w and the stated assumptions only (qoptct qualtran results --long)')
    return EXIT_OK


def build_parser():
    ap = argparse.ArgumentParser(prog='qoptct qualtran', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    c = sp.add_parser('check', help='prove a small Qualtran bloq against its call_classically (check_bloq)')
    c.add_argument('bloq', help='Python expression of a bloq, e.g. "qualtran.bloqs.arithmetic.addition.Add(QUInt(3))"')
    c.add_argument('--max-input-bits', type=int, default=12)
    c.add_argument('--out', help='write the full result as JSON')
    c.add_argument('--json', action='store_true', help='print the full result as JSON')
    c.set_defaults(fn=cmd_check)
    pr = sp.add_parser('prove', help='prove one block of Qualtran')
    pr.add_argument('block', choices=['kaliski', 'montgomery', 'ecadd', 'ecwindowaddr'])
    pr.add_argument('--n', type=int, required=True)
    pr.add_argument('--p', help='modulus, e.g. 32749 or 2**255-19 (default: 32749 at n=16, 2**255-19 at n=256; Montgomery: the largest prime <= 2**(n-1))')
    pr.add_argument('--w', type=int, default=1, help='montgomery: window size')
    pr.add_argument('--window', type=int, default=1, help='ecwindowaddr: window size')
    pr.add_argument('--original', action='store_true', help='montgomery / ecadd: the unpatched Qualtran 0.7.0 (it is not proved correct)')
    pr.add_argument('--prove-boxes', action='store_true',
                    help='ecadd / ecwindowaddr: also prove here the box contracts (multiplier, inverse, pairs) that no registry result covers; '
                         'without it a missing contract gives PROVED_CONDITIONAL (exit 3)')
    pr.add_argument('--lemmas', help='kaliski: Z3 lemma result file (default: computed once, cached in $QOPT_CACHE)')
    pr.add_argument('--canaries', default='auto', help='kaliski n=256: canary evidence file (default: cache, else the file shipped with the package)')
    pr.add_argument('--out', help='write the full result as JSON')
    pr.add_argument('--json', action='store_true')
    pr.add_argument('--verbose', action='store_true')
    pr.set_defaults(fn=cmd_prove)
    pa = sp.add_parser('pair', help='check that a block followed by its uncompute twin is the identity')
    pa.add_argument('block', choices=['kaliski', 'montgomery'])
    pa.add_argument('--n', type=int, required=True)
    pa.add_argument('--p', help='modulus (default as for prove)')
    pa.add_argument('--w', type=int, default=1, help='montgomery: window size')
    pa.add_argument('--out')
    pa.add_argument('--json', action='store_true')
    pa.set_defaults(fn=cmd_pair)
    rs = sp.add_parser('results', help='table of the proved results (registry shipped with the package)')
    rs.add_argument('--name', choices=['montgomery_mul', 'kaliski_inverse', 'pair_montgomery', 'pair_kaliski', 'ecadd', 'ecwindowaddr'])
    rs.add_argument('--n', type=int)
    rs.add_argument('--long', action='store_true', help='claim, assumptions and the reproduce command of every result')
    rs.add_argument('--json', action='store_true')
    rs.set_defaults(fn=cmd_results)
    return ap


def main(argv=None):
    ap = build_parser()
    a = ap.parse_args(argv)
    try:
        return a.fn(a)
    except UsageError as e:
        print(f'error: {e}', file=sys.stderr)
        return EXIT_USAGE


if __name__ == '__main__':
    sys.exit(main())
