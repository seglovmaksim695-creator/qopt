#!/usr/bin/env python3
"""Build qoptct/foreign/auto/data/registry.json from results/auto/*.json.

    python scripts/build_auto_registry.py            # (re)write the registry
    python scripts/build_auto_registry.py --check    # exit 1 if the registry on disk differs from a fresh build

Every entry copies its facts from result files only (verdict, canaries, cross-check with 0.7; each file is named
with its SHA-256).  The comparison with 0.7 is done here on the recorded hashes: the program hash of the verdict
against the entry of the 0.7 registry (qoptct/foreign/qualtran/data/registry.json) for the same block and instance,
and the QASM hashes of the pieces against the hashes that occur in the 0.7 result file.  No Qualtran is imported.
"""
import argparse
import hashlib
import json
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(REPO, 'results', 'auto')
RES07 = os.path.join(REPO, 'results', 'qualtran')
REG07 = os.path.join(REPO, 'qoptct', 'foreign', 'qualtran', 'data', 'registry.json')
OUT = os.path.join(REPO, 'qoptct', 'foreign', 'auto', 'data', 'registry.json')
SCHEMA = 'qoptct.auto.registry/1'
P255 = 2 ** 255 - 19

BLOCK = {
    'montgomery': 'DirtyOutOfPlaceMontgomeryModMul (Qualtran kind c)',
    'modadd': 'ModAdd', 'cmodadd': 'CModAdd', 'modsub': 'ModSub', 'cmodsub': 'CModSub', 'modneg': 'ModNeg',
    'cmodneg': 'CModNeg', 'moddbl': 'ModDbl', 'ecadd_step1': '_ECAddStepOne of ECAdd (Qualtran kind e)'}
CLAIM = {
    'montgomery': 'for every n-bit x and every y < p: target = T_N - p*[T_N >= p], reduced = [T_N >= p], qrom_indices = '
                  'the digits m_i of the windows, x and y unchanged, every other wire clean, no phase; T_0 = 0, '
                  'T_{i+1} = (T_i + x_i*y + c(m_i)*p) / 2^w with m_i = (T_i + x_i*y) mod 2^w (Montgomery recurrence, '
                  'PR 1997); the invariant T_i < y + p is proved for every window',
    'modadd': 'y = (x + y) mod p, x unchanged, other wires clean, no phase; for x <= p, y < p',
    'cmodadd': 'y = (x + y) mod p if ctrl = cv else y; ctrl, x unchanged; other wires clean, no phase; for x <= p, y < p',
    'modsub': 'y = (y - x) mod p, x unchanged, other wires clean, no phase; for x, y < p',
    'cmodsub': 'y = (y - x) mod p if ctrl = cv else y; ctrl, x unchanged; other wires clean, no phase; for x, y < p',
    'modneg': 'x = (-x) mod p, other wires clean, no phase; for x < p',
    'cmodneg': 'x = (-x) mod p if ctrl = cv else x; ctrl unchanged; other wires clean, no phase; for x < p',
    'moddbl': 'x = 2x mod p, other wires clean, no phase; for x < p',
    'ecadd_step1': 'f1 = [a = x], f3 = [a = b = 0], f4 = [x = y = 0], ctrl = [not f3 and not f4 and not (f1 and '
                   'b = -y mod p)], f2 = [b = -y mod p] and not ctrl; a, b, x, y unchanged; other wires clean, no '
                   'phase; for a, b, x, y < p (statement S1 of the 0.7 ECAdd proof)'}
EXTRA = {
    'montgomery': ['Qualtran 0.7.0 with patch kind c (PR 1997 + exact-integer tables), as in 0.7',
                   'the statement is the recurrence; that T_N = x*y*2^-n (mod p) is the telescoping identity of '
                   'Montgomery multiplication, not proved here (as in 0.7)'],
    'ecadd_step1': ['Qualtran 0.7.0 with patch kind e (fixed ECAdd), as in 0.7'],
}
NAME07 = {'montgomery': 'montgomery_mul'}
HEX = re.compile(r'\b[0-9a-f]{64}\b')


def sha_file(p):
    with open(p, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def load(p):
    with open(p, encoding='utf8') as f:
        return json.load(f)


def key(params, target):
    return (target, int(params['n']), str(params.get('p')), int(params.get('w', 1)) if target == 'montgomery' else None)


def reproduce(target, params):
    n = int(params['n'])
    cmd = f'python -m qoptct.foreign.auto prove {target} --n {n}'
    if target == 'montgomery':
        cmd += f" --w {params.get('w', 1)}"
    p = int(params['p'])
    default = P255 if n == 256 else None
    if p != default:
        cmd += f' --p {p}'
    if params.get('cv', 1) != 1:
        cmd += f" --cv {params['cv']}"
    return cmd


def compare_07(target, v, compares, reg07):
    params = v['params']
    if target in NAME07:
        for e in reg07.get('entries', []):
            if e['name'] == NAME07[target] and e['n'] == params['n'] and str(e['p']) == str(params['p']) and \
                    e.get('w') == params.get('w'):
                src = os.path.join(RES07, e['source']['file'])
                pieces = {pc['type']: pc.get('qasm_sha256') for pc in v.get('pieces', []) if pc.get('qasm_sha256')}
                found = set(HEX.findall(open(src, encoding='utf8').read())) if os.path.isfile(src) else set()
                return dict(entry_07=e['id'], program_07=e['program']['value'],
                            program_equal=(e['program']['value'] == v['program']['ops_sha256']),
                            pieces_with_hash_in_07=sorted(t for t, h in pieces.items() if h in found),
                            pieces_without=sorted(t for t, h in pieces.items() if h not in found))
        return None
    c = compares.get(key(params, target))
    if c is None:
        return None
    d, fname, sha = c
    return dict(file=fname, sha256=sha, program_07=d['digest_07'], program_equal=bool(d['equal']) and
                d['digest_07'] == v['program']['ops_sha256'], same_inputs=d.get('same_inputs'),
                ops_07=d['ops_07'], ops_auto=d['ops_auto'])


def build():
    files = sorted(f for f in os.listdir(RES) if f.endswith('.json'))
    verdicts, canaries, compares = [], {}, {}
    for f in files:
        d = load(os.path.join(RES, f))
        sha = sha_file(os.path.join(RES, f))
        if d.get('schema') == 'qoptct.auto/1':
            verdicts.append((f, sha, d))
        elif d.get('schema') == 'qoptct.auto.canaries/1':
            k = key(d['params'], d['target'])
            canaries[k] = (d, f, sha)
        elif d.get('kind') == 'compare_ecadd_step':
            k = key(dict(n=d['n'], p=d['p']), f"ecadd_step{d['step']}")
            compares[k] = (d, f, sha)
    reg07 = load(REG07) if os.path.isfile(REG07) else {}
    out = []
    for f, sha, v in verdicts:
        t, params = v['target'], v['params']
        if params.get('drop_invariant'):
            continue
        c = canaries.get(key(params, t))
        can = None
        if c is not None:
            d, cf, csha = c
            can = dict(file=cf, sha256=csha, total=d['total'], counted=d['counted'], joints=d['joints'],
                       gates=d['gates'], not_proved=d['not_proved'], proved=d['proved'], ok=d['ok'],
                       source_sha256=d.get('source_sha256'), program_ops_sha256=d.get('program_ops_sha256'))
        qv = v.get('qualtran_variant') or {}
        e = dict(id=f"auto/{f[:-5]}", target=t, block=BLOCK.get(t, t), claim=CLAIM.get(t, v['spec']['name']),
                 status=v['verdict'], reasons=v.get('reasons', []), n=params['n'], p=str(params['p']),
                 w=params.get('w') if t == 'montgomery' else None, cv=params.get('cv'),
                 qualtran=dict(version=(v.get('versions') or {}).get('qualtran'), variant=qv.get('kind'),
                               files_sha256_16=qv.get('files_sha256_16')),
                 z3=(v.get('versions') or {}).get('z3-solver'),
                 program={k: v['program'].get(k) for k in ('nq', 'ops', 'ops_sha256')},
                 spec=v['spec'], source_sha256=v.get('source_sha256'),
                 end_to_end={k: (v.get('end_to_end') or {}).get(k) for k in ('cases', 'mismatches', 'skipped')},
                 root={k: (v.get('root') or {}).get(k) for k in ('status', 'detail', 'obligations', 'z3_calls')},
                 pieces=[{k: pc.get(k) for k in ('type', 'mode', 'status', 'count', 'ops', 'qasm_sha256', 'contract')}
                         for pc in v.get('pieces', [])],
                 canaries=can, compare_07=compare_07(t, v, compares, reg07),
                 assumptions=list(v.get('assumptions', [])) + EXTRA.get(t, []),
                 reproduce=reproduce(t, params) + f' --out results/auto/{f}',
                 source=dict(file=f, sha256=sha), seconds=v.get('seconds'))
        out.append(e)
    out.sort(key=lambda e: (e['target'] != 'montgomery', e['target'], e['n'], e.get('w') or 0))
    return dict(schema=SCHEMA, note='generated by scripts/build_auto_registry.py from results/auto; do not edit',
                entries=out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    reg = build()
    text = json.dumps(reg, indent=1, sort_keys=False) + '\n'
    if a.check:
        old = open(OUT, encoding='utf8').read() if os.path.isfile(OUT) else ''
        if old != text:
            print('registry differs from a fresh build of results/auto', file=sys.stderr)
            return 1
        print(f"registry up to date ({len(reg['entries'])} entries)")
        return 0
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w', encoding='utf8') as f:
        f.write(text)
    print(f"wrote {OUT} ({len(reg['entries'])} entries)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
