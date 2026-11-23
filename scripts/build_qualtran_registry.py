#!/usr/bin/env python3
"""Build qoptct/foreign/qualtran/data/registry.json from results/qualtran/*.json.

    python scripts/build_qualtran_registry.py            # (re)write the registry
    python scripts/build_qualtran_registry.py --check    # exit 1 if the registry on disk differs from a fresh build

Every entry copies its facts from ONE result file (name, SHA-256 of the file are in ``source``); nothing is taken from memory.
The assumptions are the ones listed in docs/QUALTRAN.md for the block, plus those the result file itself records.
"""
import argparse
import glob
import hashlib
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(REPO, 'results', 'qualtran')
PKG = os.path.join(REPO, 'qoptct', 'foreign', 'qualtran')
OUT = os.path.join(PKG, 'data', 'registry.json')
P255 = 2 ** 255 - 19
SCHEMA = 'qoptct.qualtran.registry/1'

COMMON = ['rules of the strict bridge (A1 controlled-And, C1 QROAMCleanAdjoint where used) and the gadget lemma for AND/erase of qoptct',
          'the qoptct ROBDD prover and lift are the unchanged ones of qoptct.verify']

ASSUME = {
    'montgomery_mul': COMMON + [
        'Qualtran 0.7.0 with patch kind c (PR 1997 + exact-integer tables) is the program that is proved; the unpatched 0.7.0 is wrong for T_N >= p',
        'the final claim x*y*2^-n mod p with result < p follows from the proved window contracts, the invariant T_i < y + p and the '
        'telescoping sum (the last two are checked by code, not by a separate proof)',
        'invariant of the generic window: interval arithmetic where the BDD does not fit in memory (n > 192 or w >= 4)',
        'windows 1..N-1 are the same program with their own wires (equal QASM and layout hash, checked)'],
    'kaliski_inverse': COMMON + [
        'the conditions M-carry, M-dbl, M-a, M-s, M-end, M-end-vf are closed by Z3 lemmas for all n, p, x (kaliski_lemmas.json)',
        'p is prime: 2^255-19 and the secp256k1 prime are accepted as known primes without a certificate',
        'the assembly of the program from the real pieces is trusted; it was compared with the full expansion at small n and with an '
        'independent simulation',
        'the canary evidence (damaged circuits that must not be proved) is part of the criterion at n = 256'],
    'pair_kaliski': COMMON + [
        'the forward block is NOT verified here; its contract is the result kaliski_inverse',
        'structural argument: mirror of the forward program up to a renaming of wires; relies on the forward block leaving its wires clean'],
    'pair_montgomery': COMMON + [
        'the forward block is NOT verified here; its contract is the result montgomery_mul',
        'structural argument: mirror of the forward program up to a renaming of wires; the table blocks are abstracted to LD/CLR '
        '(the abstraction is checked block by block)'],
    'ecadd': COMMON + [
        'the Montgomery multiplier (window 1) and the Kaliski inverse are used as boxes with contracts (ecadd/contracts.py); an instance '
        'without a proved contract is reported as ASSUMED in the verdict',
        'the three records of the step formulas (BDD, Z3, integer) agree on sample points only, not by proof',
        'lemmas about the curve (MM0, K0, C1, C2, EC_DBL, EC_TAN, EC_RNOTO, EC_RNOTP) are used in Z3 glue; checked by brute force for p <= 23, '
        'not proved for p = 2^255-19',
        'primality of the modulus; Z3 as the solver of the glue'],
    'ecwindowaddr': COMMON + [
        'rule C1 QROAMCleanAdjoint (assumptions A1-A4 in bridge/qroam_clean_adjoint.py); QROAMClean with K = 1 expanded by the bridge',
        'the inner ECAdd is a box with the contract of the result ecadd for the same (n, p)',
        'curve y^2 = x^3 + 2x + 3 and a fixed point R; window sizes 1, 2, 4'],
}
CLAIM = {
    'montgomery_mul': 'DirtyOutOfPlaceMontgomeryModMul(n, w, p) computes target = x*y*2^-n mod p (< p) with the garbage registers as stated, '
                      'all other wires clean, no phase; for all x, y < p',
    'kaliski_inverse': 'KaliskiModInverse(n, p): x_out = x^-1 * 2^(2n) mod p for 0 < x < p, 0 for x = 0; all other wires clean; for all x',
    'pair_montgomery': 'DirtyOutOfPlaceMontgomeryModMul followed by its uncompute=True twin is the identity on all wires (exact phase, every outcome)',
    'pair_kaliski': 'KaliskiModInverse followed by its uncompute=True twin is the identity on all wires (exact phase, every outcome)',
    'ecadd': 'ECAdd(n, p, 1) is the point addition in Montgomery form (all special cases), a, b, lam_r unchanged, other wires clean',
    'ecwindowaddr': '(c, Q) -> (c, Q + c*R) for all c < 2^window and x, y < p, given the ECAdd contract',
}
BLOCK = {
    'montgomery_mul': 'DirtyOutOfPlaceMontgomeryModMul', 'kaliski_inverse': 'KaliskiModInverse',
    'pair_montgomery': 'DirtyOutOfPlaceMontgomeryModMul + uncompute', 'pair_kaliski': 'KaliskiModInverse + uncompute',
    'ecadd': 'ECAdd', 'ecwindowaddr': 'ECWindowAddR'}


def sha_file(p):
    with open(p, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def canon_sha(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()


def patch_info(kind, manifest):
    if kind in (None, 'original'):
        return dict(variant='original', version=manifest['qualtran_version'], patch_files_sha256={})
    files = {}
    for name in manifest['kinds'][kind]['files']:
        info = manifest['files'][name]
        files[info['target']] = info['sha256_patched']
    return dict(variant=kind, version=manifest['qualtran_version'], patch_files_sha256=files)


def pstr(p):
    return str(p)


def pexpr(p):
    return "'2**255-19'" if int(p) == P255 else str(p)


def build():
    manifest = json.load(open(os.path.join(PKG, 'patches', 'manifest.json')))
    es = []

    def add(name, fn, d, *, n, p, w, status, kind, z3, prog, cmd, extra_assume=(), claim_extra=''):
        # the registry holds fully proved instances only: PROVED_CONDITIONAL (a box contract only assumed) is never stored
        assert status in ('PROVED', 'PROVED_OUTSIDE_G1_G3'), (fn, status)
        es.append(dict(
            id='%s/%s' % (name, os.path.splitext(fn)[0]), name=name, block=BLOCK[name], claim=CLAIM[name] + claim_extra,
            status=status, n=n, p=pstr(p), w=w, qualtran=patch_info(kind, manifest), z3=z3, program=prog,
            assumptions=list(ASSUME[name]) + list(extra_assume), reproduce=cmd,
            source=dict(file=fn, sha256=sha_file(os.path.join(RES, fn)))))

    for path in sorted(glob.glob(os.path.join(RES, '*.json'))):
        fn = os.path.basename(path)
        d = json.load(open(path))
        if fn.startswith('montgomery_n') and 'verdict' in d:
            n, w, p = d['n'], d['w'], d['p']
            cmd = 'qoptct qualtran prove montgomery --n %d --w %d' % (n, w) + ('' if int(p) == P255 and n == 256 else ' --p %s' % pexpr(p))
            add('montgomery_mul', fn, d, n=n, p=p, w=w, status=d['verdict'], kind=d['qualtran'].get('patch_kind'),
                z3=d['versions'].get('z3'), prog=dict(name='ops_sha256', value=d['program']['ops_sha256']), cmd=cmd)
        elif fn.startswith('kaliski_n') and 'verdict' in d:
            v = d['verdict']
            n, p = v['n'], v['p']
            cmd = 'qoptct qualtran prove kaliski --n %d --p %s' % (n, pexpr(p))
            add('kaliski_inverse', fn, d, n=n, p=p, w=None, status=v['status'], kind='original', z3=v['lemmas'].get('z3'),
                prog=dict(name='top_ops_sha256', value=v['program']['top_ops_sha256']), cmd=cmd)
        elif fn.startswith('pairs_kaliski'):
            n, p = d['n'], d['p']
            cmd = 'qoptct qualtran pair kaliski --n %d --p %s' % (n, pexpr(p))
            add('pair_kaliski', fn, d, n=n, p=p, w=None, status=d['status'], kind=d['qualtran_kind'], z3=None,
                prog=dict(name='qasm_sha_fwd_piece', value=d['qasm_sha_fwd_piece']), cmd=cmd)
        elif fn.startswith('pairs_montgomery'):
            n, p, w = d['n'], d['p'], d['w']
            cmd = 'qoptct qualtran pair montgomery --n %d --w %d --p %s' % (n, w, pexpr(p))
            stat = {k: d.get(k) for k in ('n', 'w', 'p', 'piece_fwd_ops', 'piece_bwd_ops', 'top_ops', 'macros', 'tokens_fwd', 'tokens_bwd')}
            add('pair_montgomery', fn, d, n=n, p=p, w=w, status=d['status'], kind=d['qualtran_kind'], z3=None,
                prog=dict(name='mirror_stats_sha256', value=canon_sha(stat)), cmd=cmd)
        elif fn.startswith('ecadd_window_'):
            n, p, w = d['n'], d['p'], d['window']
            cmd = 'qoptct qualtran prove ecwindowaddr --n %d --window %d' % (n, w) + ('' if int(p) == P255 and n == 256 else ' --p %s' % pexpr(p))
            e4 = {k: d['e4'].get(k) for k in ('logical_ops', 'regions', 'R', 'curve', 'stats')}
            e4['stats'] = {k: v for k, v in (e4['stats'] or {}).items() if k != 'seconds'}
            add('ecwindowaddr', fn, d, n=n, p=p, w=w, status=d['status'], kind=d['qualtran_kind'], z3=d.get('z3'),
                prog=dict(name='checked_stats_sha256', value=canon_sha(e4)), cmd=cmd)
        elif fn.startswith('ecadd_') and 'steps' in d:
            n, p = d['n'], d['p']
            orig = not d['patched']
            cmd = 'qoptct qualtran prove ecadd --n %d' % n + ('' if int(p) == P255 and n == 256 else ' --p %s' % pexpr(p)) + (' --original' if orig else '')
            steps = {k: {'logical_ops': v.get('logical_ops'), 'stats': {a: b for a, b in (v.get('stats') or {}).items() if a != 'seconds'}}
                     for k, v in d['steps'].items()}
            cs = d.get('contracts_status', {})
            extra = []
            for box, st in cs.items():
                for part in ('forward', 'pair'):
                    if st.get(part) and st[part]['status'] != 'PROVED_FOR_THIS_INSTANCE':
                        extra.append('contract of the box %s (%s, %s) is ASSUMED for this instance' % (box, part, st[part]['proved_by']))
            if orig:
                extra.append('the original Qualtran 0.7.0 scheme: the point addition holds outside the groups G1-G3 of finding 05; that the scheme is dirty on G1-G3 was shown in round 3 (finding 05), not by this package')
            add('ecadd', fn, d, n=n, p=p, w=None, status=d['status'], kind=d['qualtran_kind'], z3=d.get('z3'),
                prog=dict(name='steps_stats_sha256', value=canon_sha(steps)), cmd=cmd, extra_assume=extra,
                claim_extra=' (original scheme: correct outside G1-G3)' if orig else '')
    support = []
    for fn in ('kaliski_lemmas.json', 'kaliski_canaries_n256_p25519.json', 'kaliski_canaries_T8b_n256_p25519.json', 'kaliski_mutants_n256.json',
               'montgomery_regression.json', 'montgomery_compare_C_CFIX.json'):
        if os.path.isfile(os.path.join(RES, fn)):
            support.append(dict(file=fn, sha256=sha_file(os.path.join(RES, fn))))
    order = {'montgomery_mul': 0, 'kaliski_inverse': 1, 'pair_montgomery': 2, 'pair_kaliski': 3, 'ecadd': 4, 'ecwindowaddr': 5}
    es.sort(key=lambda e: (order[e['name']], e['n'], e['w'] or 0, e['p'], e['id']))
    return dict(schema=SCHEMA, qualtran_version=manifest['qualtran_version'],
                note='generated by scripts/build_qualtran_registry.py from results/qualtran/*.json; only PROVED-type statuses are recorded as results',
                entries=es, support_files=support)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args(argv)
    reg = build()
    txt = json.dumps(reg, indent=1, ensure_ascii=False) + '\n'
    if a.check:
        cur = open(OUT, encoding='utf8').read() if os.path.isfile(OUT) else ''
        if cur != txt:
            print('registry.json differs from a fresh build', file=sys.stderr)
            return 1
        print('registry.json is up to date (%d entries)' % len(reg['entries']))
        return 0
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w', encoding='utf8') as f:
        f.write(txt)
    print('wrote %s (%d entries)' % (OUT, len(reg['entries'])))
    return 0


if __name__ == '__main__':
    sys.exit(main())
