"""Registry of the proved results about Qualtran blocks.

The registry is the file ``data/registry.json`` of the installed package.  It is generated from the result files of the
repository (``results/qualtran/*.json``) by ``scripts/build_qualtran_registry.py``; every entry names the source file and its
SHA-256, so a checkout can verify that the registry still matches the evidence (``verify_sources``).

Entry fields: ``id``, ``name`` (``montgomery_mul``, ``kaliski_inverse``, ``pair_montgomery``, ``pair_kaliski``, ``ecadd``,
``ecwindowaddr``), ``block``, ``claim``, ``status``, ``n``, ``p``, ``w``, ``qualtran`` (variant, version, SHA-256 of the
patch files), ``z3``, ``program`` (name and value of the hash that identifies the checked program), ``assumptions``,
``reproduce`` (the command) and ``source`` (file, SHA-256).

A status other than ``PROVED`` is shown as it is (``PROVED_OUTSIDE_G1_G3`` for the original ``ECAdd``); it is never rounded up.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

__all__ = ['load_registry', 'entries', 'format_table', 'format_entry', 'verify_sources', 'packaged_canaries', 'NAMES']

NAMES = ('montgomery_mul', 'kaliski_inverse', 'pair_montgomery', 'pair_kaliski', 'ecadd', 'ecwindowaddr')
_DATA = Path(__file__).resolve().parent / 'data'


def load_registry():
    """The registry as a dictionary (``entries`` is the list of results)."""
    with open(_DATA / 'registry.json', encoding='utf8') as f:
        return json.load(f)


def entries(name=None, n=None):
    """Entries of the registry, optionally only one ``name`` and/or one ``n``."""
    out = load_registry()['entries']
    if name is not None:
        out = [e for e in out if e['name'] == name]
    if n is not None:
        out = [e for e in out if e['n'] == n]
    return out


def _short(s, k=12):
    s = str(s)
    return s if len(s) <= k else s[:k - 1] + '~'


_NAMED_P = {2 ** 255 - 19: '2^255-19', 2 ** 256 - 2 ** 32 - 977: 'secp256k1', 2 ** 127 - 1: '2^127-1', 2 ** 63 - 25: '2^63-25',
            2 ** 256 - 2 ** 224 + 2 ** 192 + 2 ** 96 - 1: 'P-256'}


def _pname(p):
    try:
        return _NAMED_P.get(int(p), _short(p, 14))
    except ValueError:
        return _short(p, 14)


def format_table(es):
    """Plain-text table of the entries."""
    rows = [('name', 'n', 'p', 'w', 'qualtran', 'status', 'program hash', 'assumptions')]
    for e in es:
        rows.append((e['name'], str(e['n']), _pname(e['p']), '-' if e.get('w') is None else str(e['w']),
                     e['qualtran']['variant'], e['status'], _short(e['program']['value'], 12), str(len(e['assumptions']))))
    width = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    lines = ['  '.join(c.ljust(width[i]) for i, c in enumerate(r)).rstrip() for r in rows]
    lines.insert(1, '  '.join('-' * w for w in width))
    return '\n'.join(lines)


def format_entry(e):
    """Long text of one entry: claim, variant, assumptions, command."""
    q = e['qualtran']
    lines = [f"{e['id']}", f"  block:      {e['block']}", f"  claim:      {e['claim']}", f"  status:     {e['status']}",
             f"  n, p, w:    {e['n']}, {e['p']}, {e.get('w')}",
             f"  qualtran:   {q['version']} ({q['variant']})" + ('; patch files: ' + ', '.join(f"{k} {v[:12]}" for k, v in q['patch_files_sha256'].items()) if q.get('patch_files_sha256') else ''),
             f"  z3:         {e.get('z3')}", f"  program:    {e['program']['name']} = {e['program']['value']}",
             f"  source:     {e['source']['file']} (sha256 {e['source']['sha256'][:16]})", f"  reproduce:  {e['reproduce']}",
             '  assumptions:']
    lines += [f'    - {a}' for a in e['assumptions']]
    return '\n'.join(lines)


def verify_sources(results_dir):
    """Compare the source files of the registry in ``results_dir`` with the recorded SHA-256; returns a list of problems."""
    problems = []
    for e in entries():
        p = Path(results_dir) / e['source']['file']
        if not p.is_file():
            problems.append(f"{e['id']}: missing {p}")
        elif hashlib.sha256(p.read_bytes()).hexdigest() != e['source']['sha256']:
            problems.append(f"{e['id']}: {p.name} differs from the registry")
    return problems


def packaged_canaries():
    """Path of the canary evidence for Kaliski at n=256 shipped with the package (``None`` if it is not there).

    The file is accepted by ``kaliski.verdict`` only if it was made by exactly this client (SHA-256 of the six source files)."""
    p = _DATA / 'kaliski_canaries_n256.json'
    return str(p) if p.is_file() else None
