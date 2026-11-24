# SPDX-License-Identifier: Apache-2.0
"""Registry of the proved results of qoptct.dyadic (the file ``data/registry.json`` of the installed package).

Built from ``results/dyadic/*.json`` by ``scripts/build_dyadic_registry.py``.  Every entry is the largest proved size of a
family, with its canaries (all REFUTED), its assumptions, the command that reproduces it and the SHA-256 of its evidence.
"""
import hashlib
import json
from pathlib import Path

__all__ = ['load_registry', 'entries', 'format_table', 'format_entry', 'verify_sources']

_DATA = Path(__file__).resolve().parent / 'data'


def load_registry():
    with open(_DATA / 'registry.json', encoding='utf8') as f:
        return json.load(f)


def entries(name=None, n=None):
    out = load_registry()['entries']
    if name is not None:
        out = [e for e in out if e['name'] == name]
    if n is not None:
        out = [e for e in out if e['n'] == n]
    return out


def format_table(es):
    rows = [('name', 'n', 'status', 'qubits', 'canaries (refuted/total)', 'assumptions')]
    for e in es:
        c = e['canaries']
        rows.append((e['name'], str(e['n']), e['status'], str(e.get('qubits')), '%d/%d' % (c['refuted'], c['total']),
                     'fourier-shift lemma' if any('fourier-shift' in a for a in e['assumptions']) else '-'))
    width = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    lines = ['  '.join(c.ljust(width[i]) for i, c in enumerate(r)).rstrip() for r in rows]
    lines.insert(1, '  '.join('-' * w for w in width))
    return '\n'.join(lines)


def format_entry(e):
    lines = [e['id'], '  block:      %s' % e['block'], '  claim:      %s' % e['claim'], '  status:     %s' % e['status'],
             '  size:       n = %s (%s qubits, %s gates)' % (e['n'], e.get('qubits'), e.get('gates')),
             '  canaries:   %d of %d refuted' % (e['canaries']['refuted'], e['canaries']['total']),
             '  evidence:   %s (sha256 %s)' % (e['evidence']['file'], e['evidence']['sha256'][:16]),
             '  reproduce:  %s' % e['reproduce']]
    if e.get('first_not_proved'):
        f = e['first_not_proved']
        lines.append('  next size:  n = %s is %s' % (f['n'], f['verdict']))
    lines.append('  assumptions:' if e['assumptions'] else '  assumptions: none beyond the contract of qoptct.dyadic')
    lines += ['    - %s' % a for a in e['assumptions']]
    return '\n'.join(lines)


def verify_sources(results_dir):
    """Compare the evidence files in ``results_dir`` with the recorded SHA-256; returns a list of problems."""
    problems = []
    for e in entries():
        p = Path(results_dir) / Path(e['evidence']['file']).name
        if not p.is_file():
            problems.append('%s: missing %s' % (e['id'], p))
        elif hashlib.sha256(p.read_bytes()).hexdigest() != e['evidence']['sha256']:
            problems.append('%s: %s differs from the registry' % (e['id'], p.name))
    return problems
