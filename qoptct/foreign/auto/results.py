# SPDX-License-Identifier: Apache-2.0
"""Registry of the results of the generic mechanism.

The registry is ``data/registry.json`` of this package.  ``scripts/build_auto_registry.py`` generates it from the result
files of the repository (``results/auto/*.json``: verdicts, canaries, cross-checks with 0.7); every entry names its
source files and their SHA-256, so a checkout can verify that the registry still matches the evidence
(``verify_sources``).  Nothing is rounded up: an entry is ``PROVED`` only if its verdict file says so, and the
criterion of the round (``criterion``) is reported separately: independent end-to-end simulation without a mismatch,
at least 20 canaries none of which is proved, at least 10 of them at the joints of the composition.

    python -m qoptct.foreign.auto results            # table
    python -m qoptct.foreign.auto results --long     # every entry with its command, hashes and assumptions
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

__all__ = ['load_registry', 'entries', 'format_table', 'format_entry', 'verify_sources', 'criterion']

_DATA = Path(__file__).resolve().parent / 'data'


def load_registry():
    with open(_DATA / 'registry.json', encoding='utf8') as f:
        return json.load(f)


def entries(target=None, n=None):
    out = load_registry()['entries']
    if target is not None:
        out = [e for e in out if e['target'] == target]
    if n is not None:
        out = [e for e in out if e['n'] == n]
    return out


def criterion(e):
    """(met, list of what is missing) for the criterion of a new PROVED"""
    miss = []
    if e['status'] != 'PROVED':
        miss.append(f"verdict {e['status']}")
    e2e = e.get('end_to_end') or {}
    if not e2e.get('cases') or e2e.get('mismatches'):
        miss.append('independent end-to-end simulation')
    c = e.get('canaries')
    if not c:
        miss.append('canaries at this size')
    else:
        if c.get('proved'):
            miss.append(f"canaries proved: {c['proved']}")
        if (c.get('counted') or 0) < 20:
            miss.append(f"only {c.get('counted')} counted canaries (< 20)")
        if (c.get('joints') or 0) < 10:
            miss.append(f"only {c.get('joints')} joint canaries (< 10)")
        if c.get('source_sha256') != e.get('source_sha256'):
            miss.append('canaries were made by other sources than the verdict')
    return not miss, miss


def _short(s, k=12):
    s = '' if s is None else str(s)
    return s if len(s) <= k else s[:k - 1] + '~'


_NAMED_P = {2 ** 255 - 19: '2^255-19', 2 ** 256 - 2 ** 32 - 977: 'secp256k1', 2 ** 127 - 1: '2^127-1'}


def _pname(p):
    try:
        return _NAMED_P.get(int(p), _short(p, 14))
    except (TypeError, ValueError):
        return _short(p, 14)


def format_table(es):
    rows = [('target', 'n', 'p', 'w', 'status', 'canaries', 'criterion', 'program hash', '= 0.7')]
    for e in es:
        c = e.get('canaries') or {}
        met, _ = criterion(e)
        same = e.get('compare_07') or {}
        rows.append((e['target'], str(e['n']), _pname(e['p']), '-' if e.get('w') is None else str(e['w']), e['status'],
                     f"{c.get('not_proved', '-')}/{c.get('total', '-')}" if c else '-', 'met' if met else 'NOT met',
                     _short(e['program']['ops_sha256'], 12),
                     '-' if not same else ('yes' if same.get('program_equal') else 'no')))
    width = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    lines = ['  '.join(c.ljust(width[i]) for i, c in enumerate(r)).rstrip() for r in rows]
    lines.insert(1, '  '.join('-' * w for w in width))
    return '\n'.join(lines)


def format_entry(e):
    met, miss = criterion(e)
    q = e['qualtran']
    lines = [e['id'], f"  block:      {e['block']}", f"  claim:      {e['claim']}", f"  status:     {e['status']}",
             f"  criterion:  {'met' if met else 'NOT met: ' + '; '.join(miss)}",
             f"  n, p, w:    {e['n']}, {e['p']}, {e.get('w')}",
             f"  qualtran:   {q.get('version')} ({q.get('variant')})",
             f"  z3:         {e.get('z3')}",
             f"  program:    {e['program']['nq']} wires, {e['program']['ops']} operations, "
             f"ops_sha256 = {e['program']['ops_sha256']}",
             f"  spec:       {e['spec']['name']} (sha256 {e['spec']['sha256'][:16]})",
             f"  sources:    sha256 {e.get('source_sha256')}"]
    if e.get('compare_07'):
        lines.append(f"  vs 0.7:     {json.dumps(e['compare_07'])}")
    c = e.get('canaries')
    if c:
        lines.append(f"  canaries:   {c['not_proved']} of {c['total']} not proved, {c['counted']} counted "
                     f"({c['joints']} joints, {c['gates']} gates); proved: {c['proved'] or 'none'}; {c['file']}")
    lines += [f"  reproduce:  {e['reproduce']}", f"  source:     {e['source']['file']} (sha256 {e['source']['sha256'][:16]})",
              '  pieces:']
    lines += [f"    {p['status']:10s} {p['mode']:7s} {p['type']:32s} x{p['count']:<4d} {(p.get('qasm_sha256') or '')[:16]}"
              for p in e.get('pieces', [])]
    lines.append('  assumptions:')
    lines += [f'    - {a}' for a in e['assumptions']]
    return '\n'.join(lines)


def verify_sources(results_dir):
    """compare the result files in ``results_dir`` with the SHA-256 recorded by the registry; returns the problems"""
    problems = []
    for e in entries():
        files = [e['source']] + ([e['canaries']] if e.get('canaries') else []) + \
            ([e['compare_07']] if (e.get('compare_07') or {}).get('file') else [])
        for s in files:
            p = Path(results_dir) / s['file']
            if not p.is_file():
                problems.append(f"{e['id']}: missing {p}")
            elif hashlib.sha256(p.read_bytes()).hexdigest() != s['sha256']:
                problems.append(f"{e['id']}: {p.name} differs from the registry")
    return problems
