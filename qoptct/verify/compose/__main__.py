# SPDX-License-Identifier: Apache-2.0
"""python -m qoptct.verify.compose gen|verify|pstd ...   (see qoptct.verify.compose)

    gen    --prime P --reduction R [--clean] -o DIR      write DIR/<name>.qasm (streamed) and DIR/<name>.hints.json
    verify --prime P --reduction R [--clean] FILE.qasm [--hints FILE.json] [--json]
    pstd   --prime P --reduction R [--count N | --all]   execute P_std on inputs: test of the lemmas M1..M6

P: a decimal / 0x number or one of: secp256k1, p256, curve25519, p384, p521, rand256.  Exit code of verify: 0 PROVED, 1 REFUTED, 2 UNKNOWN.

Elliptic-curve point addition P -> P + Q (qoptct.crypto.ecpoint.ec_point_add, compose.ecadd):

    ecgen    --curve C [--k K | --qx X --qy Y] [--levels L] -o DIR     write DIR/<name>.qasm (streamed) and .hints.json
    ecverify --curve C [--k K | --qx X --qy Y] FILE.qasm [--hints FILE.json] [--json] [--no-test]
    ecwec    --curve C [--k K | --qx X --qy Y] [--levels L] [--count N]  execute the word-level program W_ec: test
                                                                        of the assumptions G1..G4, M1..M6

C: secp256k1, p256 (Q = K * G, default K = 123456789) or "p,a,b" (then --qx --qy).  Without a Q option ecverify takes
Q from the certificate and prints it.  Exit code of ecverify: 0 PROVED, 1 REFUTED, 2 UNKNOWN, 3 bad parameters.
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


def _ec_curve(a):
    from ...crypto import ecpoint
    from .ecadd import ec_mul
    names = {'secp256k1': 'secp256k1', 'p256': 'P-256', 'P-256': 'P-256'}
    if a.curve in names:
        cv = ecpoint.CURVES[names[a.curve]]
        p, ca, cb = cv['p'], cv['a'], cv['b']
        Q = None
        if a.qx is None and (a.k is not None or a.cmd != 'ecverify'):
            Q = ec_mul(123456789 if a.k is None else int(a.k, 0), (cv['gx'], cv['gy']), p, ca)
        name = a.curve.lower().replace('-', '')
    else:
        p, ca, cb = (int(v, 0) for v in a.curve.split(','))
        Q = None
        name = 'p%d' % p if p < 10 ** 6 else 'p%dbit' % p.bit_length()
    if a.qx is not None:
        Q = (int(a.qx, 0), int(a.qy, 0))
    return p, ca, cb, Q, name


def ec_main(argv):
    ap = argparse.ArgumentParser(prog='python -m qoptct.verify.compose')
    sub = ap.add_subparsers(dest='cmd', required=True)
    for name in ('ecgen', 'ecverify', 'ecwec'):
        a = sub.add_parser(name)
        a.add_argument('--curve', required=True)
        a.add_argument('--k')
        a.add_argument('--qx')
        a.add_argument('--qy')
        a.add_argument('--levels', type=int, default=3)
        if name == 'ecgen':
            a.add_argument('-o', '--out', required=True)
        if name == 'ecverify':
            a.add_argument('qasm')
            a.add_argument('--hints')
            a.add_argument('--json', action='store_true')
            a.add_argument('--no-test', action='store_true', help='skip the sampled execution of W_ec')
            a.add_argument('--quiet', action='store_true')
        if name == 'ecwec':
            a.add_argument('--count', type=int, default=64)
    a = ap.parse_args(argv)
    from .ecadd import verify_ecadd, wec_test
    p, ca, cb, Q, cname = _ec_curve(a)
    if a.cmd == 'ecgen':
        from ...crypto import ecpoint
        if Q is None:
            ap.error('give Q: --k or --qx/--qy')
        t0 = time.time()
        C = ecpoint.ec_point_add(p, ca, cb, Q, karatsuba_levels=a.levels)
        base = os.path.join(a.out, 'ecadd_%s_L%d' % (cname, a.levels))
        os.makedirs(a.out, exist_ok=True)
        with open(base + '.hints.json', 'w') as f:
            json.dump(C.hints(), f)
        C.write_qasm(base + '.qasm')
        print('%s.qasm  %s  %.1f s' % (base, C.counts(), time.time() - t0))
        return 0
    if a.cmd == 'ecwec':
        if Q is None:
            ap.error('give Q: --k or --qx/--qy')
        t0 = time.time()
        r = wec_test(p, ca, cb, Q, count=a.count, karatsuba_levels=a.levels)
        print('W_ec on %d points: %s  (%.1f s)' % (r['inputs'], 'all assumptions and the result hold' if r['ok']
                                                    else 'VIOLATED %r' % r['violations'][:3], time.time() - t0))
        return 0 if r['ok'] else 1
    hints = json.load(open(a.hints)) if a.hints else None
    if Q is None:
        if hints is None:
            ap.error('give Q (--k or --qx/--qy) or --hints')
        Q = tuple(hints['Q'])
    prog = None if (a.quiet or a.json) else (lambda s: None if s.startswith('step') else print('  ' + s, flush=True))
    r = verify_ecadd(a.qasm, p, ca, cb, Q, hints=hints, progress=prog, test_samples=0 if a.no_test else 64,
                     karatsuba_levels=a.levels)
    if a.json:
        print(json.dumps(r, default=str))
    else:
        print('%s  (point addition P -> P + Q, n=%d, %d ops, %.1f s)' % (r['status'], r['n'], r.get('ops', 0),
                                                                          r['time']))
        print('  Q: (%#x, %#x)' % tuple(Q))
        for k in ('detail', 'where', 'counterexample', 'proved_modulo', 'assumptions', 'segments', 'mirrors',
                  'lemma_facts', 'bdd_nodes_max', 'wec_test'):
            if r.get(k) is not None:
                print('  %s: %s' % (k, r[k]))
        if r.get('algebra'):
            print('  algebra sub-goals: %d, all PROVED: %s' % (
                len(r['algebra']), all(x['status'] == 'PROVED' for x in r['algebra'])))
    return {'PROVED': 0, 'REFUTED': 1, 'ERROR': 3}.get(r['status'], 2)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] in ('ecgen', 'ecverify', 'ecwec'):
        return ec_main(argv)
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
