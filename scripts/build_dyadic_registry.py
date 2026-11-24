#!/usr/bin/env python3
"""Build qoptct/dyadic/data/registry.json from results/dyadic/*.json.

    python scripts/build_dyadic_registry.py            # (re)write the registry
    python scripts/build_dyadic_registry.py --check    # exit 1 if the registry on disk differs from a fresh build

One entry per family: the largest proved size, with its canaries (at least 20, all REFUTED) and the SHA-256 of the result
file.  The ProjectQ rows come from the exact-angle bridge (results/dyadic/pq_exact_*.json); the older dump-based files
results/dyadic/pq_*.json used a coarser angle tolerance and are not part of the registry.
"""
import argparse
import hashlib
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(REPO, 'results', 'dyadic')
OUT = os.path.join(REPO, 'qoptct', 'dyadic', 'data', 'registry.json')
SCHEMA = 'qoptct.dyadic.registry/1'

FAMILIES = {
    'qft': ('Qiskit QFTGate', 'circuit = QFT (with the final reversal of the wires), up to a global phase', 'run'),
    'qft_noswap': ('Qiskit synth_qft_full(do_swaps=False)', 'circuit = QFT without the final reversal, up to a global phase', 'run'),
    'iqft': ('Qiskit QFTGate.inverse()', 'circuit = inverse QFT, up to a global phase', 'run'),
    'qft_qftdg': ('Qiskit QFTGate then its inverse', 'circuit = identity, up to a global phase', 'run'),
    'draper_fixed': ('Qiskit DraperQFTAdder kind=fixed', 'b += a mod 2^n', 'run'),
    'draper_half': ('Qiskit DraperQFTAdder kind=half', '(b, carry) += a', 'run'),
    'cdraper': ('Qiskit DraperQFTAdder(n).control(1)', 'b += c*a mod 2^n', 'run'),
    'rgqft': ('Qiskit RGQFTMultiplier', 'out += a*b mod 2^(2n)', 'run'),
    'hrs': ('Qiskit HRSCumulativeMultiplier', '|a, b, 0> -> |a, b, a*b>', 'run'),
    'pq_exact_qft': ('ProjectQ QFT (exact-angle bridge)', 'program = QFT without the final reversal, up to a global phase', 'pq'),
    'pq_exact_addconst': ('ProjectQ AddConstant (exact-angle bridge)', 'b += c mod 2^n for a fixed mixed-bit constant c', 'pq'),
}
ASSUMPTION_LEMMA = ('fourier-shift lemma: the exact Fourier path sum of an ArithSpec equals the classical map '
                    '(geometric sum; docs/DYADIC.md, Lemma 2)')


def sha256(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def build():
    entries = []
    for fam, (block, claim, kind) in FAMILIES.items():
        path = os.path.join(RES, fam + '.json')
        rows = json.load(open(path))['rows']
        ok = [r for r in rows if r.get('verdict') == 'PROVED' and r.get('canaries')]
        if not ok:
            raise SystemExit('%s: no PROVED row with canaries' % fam)
        m = max(ok, key=lambda r: r['n'])
        cans = m['canaries']
        refuted = sum(1 for c in cans if c['verdict'] == 'REFUTED')
        if len(cans) < 20 or refuted != len(cans) or m.get('canaries_proved'):
            raise SystemExit('%s n=%s: canary criterion not met (%d of %d)' % (fam, m['n'], refuted, len(cans)))
        bigger = [r for r in rows if r['n'] > m['n']]
        nxt = min(bigger, key=lambda r: r['n']) if bigger else None
        pc = m.get('product_check')
        cmd = ('python scripts/dyadic_runs.py run %s %d --canaries' if kind == 'run' else
               'python scripts/projectq_dyadic_runs.py run %s %d --canaries')
        fam_arg = fam if kind == 'run' else fam[len('pq_exact_'):]
        entries.append(dict(
            id='dyadic/%s_n%d' % (fam, m['n']), name=fam, block=block, claim=claim, status='PROVED', n=m['n'],
            qubits=m.get('qubits'), gates=m.get('gates'), method=m.get('method'),
            global_phase=m.get('global_phase'),
            canaries=dict(total=len(cans), refuted=refuted, proved=0),
            product_check=pc, assumptions=list(m.get('assumptions') or []) +
            ([ASSUMPTION_LEMMA] if fam in ('draper_fixed', 'draper_half', 'cdraper', 'rgqft') and not m.get('assumptions') else []),
            first_not_proved=(dict(n=nxt['n'], verdict=nxt['verdict'], reason=(nxt.get('reason') or '')[:200]) if nxt else None),
            versions=dict(projectq=m.get('projectq')) if kind == 'pq' else {},
            reproduce=cmd % (fam_arg, m['n']),
            evidence=dict(file='results/dyadic/%s.json' % fam, sha256=sha256(path))))
    return dict(schema=SCHEMA, note='each entry: the largest proved size of a family, with its canaries; see docs/DYADIC.md',
                entries=entries)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args(argv)
    text = json.dumps(build(), indent=1, ensure_ascii=False) + '\n'
    if a.check:
        cur = open(OUT).read() if os.path.exists(OUT) else ''
        print('registry matches' if cur == text else 'registry DIFFERS from the result files')
        return 0 if cur == text else 1
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    open(OUT, 'w').write(text)
    print('wrote %s (%d entries)' % (os.path.relpath(OUT, REPO), len(json.loads(text)['entries'])))
    return 0


if __name__ == '__main__':
    sys.exit(main())
