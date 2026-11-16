# SPDX-License-Identifier: Apache-2.0
"""python -m qoptct.verify.compose gen|verify|pstd ...   (see qoptct.verify.compose)

    gen    --prime P --reduction R [--clean] -o DIR      write DIR/<name>.qasm (streamed) and DIR/<name>.hints.json
    verify --prime P --reduction R [--clean] FILE.qasm [--hints FILE.json] [--json]
    pstd   --prime P --reduction R [--count N | --all]   execute P_std on inputs: test of the lemmas M1..M6

P: a decimal / 0x number or one of: secp256k1, p256, curve25519, p384, p521, rand256.  Exit code of verify: 0 PROVED, 1 REFUTED, 2 UNKNOWN.
"""
import argparse
import json
import os
import sys
import time

from . import PRIMES, named_prime, verify_kaliski, pstd_test


def _prime(s):
    if s in PRIMES or s == 'rand256':
        return named_prime(s), s
    p = int(s, 0)
    return p, 'p%d' % p if p < 10 ** 6 else 'p%dbit_%x' % (p.bit_length(), p & 0xffffffff)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m qoptct.verify.compose')
    sub = ap.add_subparsers(dest='cmd', required=True)
    for name in ('gen', 'verify', 'pstd'):
        a = sub.add_parser(name)
        a.add_argument('--prime', required=True)
        a.add_argument('--reduction', default='generic', choices=('generic', 'special'))
        a.add_argument('--clean', action='store_true')
        if name == 'gen':
            a.add_argument('-o', '--out', required=True)
        if name == 'verify':
            a.add_argument('qasm')
            a.add_argument('--hints')
            a.add_argument('--json', action='store_true')
            a.add_argument('--no-test', action='store_true', help='skip the sampled execution of P_std')
            a.add_argument('--quiet', action='store_true')
        if name == 'pstd':
            a.add_argument('--count', type=int, default=1024)
            a.add_argument('--all', action='store_true')
            a.add_argument('--hints')
    a = ap.parse_args(argv)
    p, pname = _prime(a.prime)
    if a.cmd == 'gen':
        from ...crypto.kaliski import kaliski_inverse
        t0 = time.time()
        C = kaliski_inverse(p, reduction=a.reduction, clean=a.clean)
        base = os.path.join(a.out, 'kaliski_%s_%s_%s' % (pname, a.reduction, 'clean' if a.clean else 'fwd'))
        os.makedirs(a.out, exist_ok=True)
        with open(base + '.hints.json', 'w') as f:
            json.dump(C.hints(), f)
        C.write_qasm(base + '.qasm')
        print('%s.qasm  %s  %.1f s' % (base, C.counts(), time.time() - t0))
        return 0
    if a.cmd == 'pstd':
        hints = json.load(open(a.hints)) if a.hints else None
        t0 = time.time()
        r = pstd_test(p, a.reduction, xs=range(1, p) if a.all else None, hints=hints, count=a.count)
        print('P_std on %d inputs: %s  (%.1f s)' % (r['inputs'], 'all assumptions and the result hold' if r['ok']
                                                    else 'VIOLATED %r' % r['violations'][:3], time.time() - t0))
        return 0 if r['ok'] else 1
    hints = json.load(open(a.hints)) if a.hints else None
    prog = None if (a.quiet or a.json) else (lambda s: print('  ' + s, flush=True))
    r = verify_kaliski(a.qasm, p, a.reduction, a.clean, hints=hints, progress=prog,
                       test_samples=0 if a.no_test else 128)
    if a.json:
        print(json.dumps(r, default=str))
    else:
        print('%s  (n=%d, %s, %s, %d ops, %.1f s)' % (r['status'], r['n'], r['reduction'],
                                                       'clean' if r['clean'] else 'forward', r.get('ops', 0), r['time']))
        for k in ('detail', 'where', 'counterexample', 'assumptions', 'segments', 'lemma_facts', 'bdd_nodes_max',
                  'pstd_test'):
            if r.get(k) is not None:
                print('  %s: %s' % (k, r[k]))
    return {'PROVED': 0, 'REFUTED': 1}.get(r['status'], 2)


if __name__ == '__main__':
    sys.exit(main())
