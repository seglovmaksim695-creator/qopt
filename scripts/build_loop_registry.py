#!/usr/bin/env python3
"""Build qoptct/loop/data/registry.json (the best-known table of the circle) from results/loop/state/journal.jsonl.

    python scripts/build_loop_registry.py            # (re)write the registry
    python scripts/build_loop_registry.py --check    # exit 1 if the registry on disk differs from a fresh build

The table is computed from the journal only (qoptct.loop.journal.table_from); the entry carries the SHA-256 of the journal and
the command that re-checks every record (chain, file hashes, costs, proofs, simulation): ``qoptct loop verify --state
results/loop/state``.
"""
import argparse
import hashlib
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
STATE = os.path.join(REPO, 'results', 'loop', 'state')
OUT = os.path.join(REPO, 'qoptct', 'loop', 'data', 'registry.json')
SCHEMA = 'qoptct.loop.registry/1'


def build():
    from qoptct.loop.journal import State, table_from
    recs, probs = State(STATE).read()
    if probs:
        raise SystemExit('journal problems: %s' % probs[:3])
    tab = table_from(recs)
    with open(os.path.join(STATE, 'journal.jsonl'), 'rb') as f:
        jsha = hashlib.sha256(f.read()).hexdigest()
    entries = []
    for t in sorted(tab):
        b, a = tab[t]['best'], tab[t]['baseline']
        entries.append(dict(
            id='loop/%s' % t, name='loop', target=t, status='PROVED', n=None, baseline=a, best=b,
            records=tab[t]['records'], front=[e['cost'] for e in tab[t]['front']],
            improves_baseline=all(b[k] <= a[k] for k in a) and b != a,
            claim='every record is a candidate circuit proved equal to the reference of the target by the provers of '
                  'qoptct.verify, with a simulation check and at least 20 canaries (all not proved)',
            assumptions=['for qc: targets the baseline is the first proved result of the optimizer, not the original file'
                         if t.startswith('qc_') else 'baseline = the reference circuit of the target'],
            reproduce='qoptct loop verify --state results/loop/state'))
    return dict(schema=SCHEMA, journal=dict(file='results/loop/state/journal.jsonl', sha256=jsha, records=len(recs)),
                entries=entries)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args(argv)
    text = json.dumps(build(), indent=1, ensure_ascii=False, sort_keys=True) + '\n'
    if a.check:
        cur = open(OUT).read() if os.path.exists(OUT) else ''
        print('registry matches' if cur == text else 'registry DIFFERS from the journal')
        return 0 if cur == text else 1
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    open(OUT, 'w').write(text)
    print('wrote %s (%d targets)' % (os.path.relpath(OUT, REPO), len(json.loads(text)['entries'])))
    return 0


if __name__ == '__main__':
    sys.exit(main())
