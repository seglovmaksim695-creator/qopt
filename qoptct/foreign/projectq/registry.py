"""Registry of the proved results about ProjectQ blocks.

The registry is the file ``data/registry.json`` of the installed package.  It is generated from the result files of the
repository (``results/projectq/*.json``) by ``scripts/build_projectq_registry.py``; every entry names its source file and the
SHA-256 of it, so a checkout can verify that the registry still matches the evidence (``verify_sources``).

Entry fields: ``id``, ``name`` (``binary_ecc_multiplier``, ``binary_ecc_point_addition``, ``projectq_library``), ``block``,
``claim``, ``status`` (always PROVED here), ``size`` (m or n), ``projectq`` (version), ``source`` (the Binary_ECC file and its
SHA-256, or None), ``program`` (hash of the logical program), ``garbage`` (number of pinned garbage wires and the pin hash),
``canaries`` (decided / not proved / equivalent), ``assumptions``, ``reproduce`` and ``evidence`` (result file, SHA-256).
Results that are not PROVED (for example the wrong decomposition of ComparatorQuantum for n >= 3) are in ``findings``.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

__all__ = ['load_registry', 'entries', 'format_table', 'format_entry', 'verify_sources', 'NAMES']

NAMES = ('binary_ecc_multiplier', 'binary_ecc_point_addition', 'projectq_library')
_DATA = Path(__file__).resolve().parent / 'data'


def load_registry():
    with open(_DATA / 'registry.json', encoding='utf8') as f:
        return json.load(f)


def entries(name=None, size=None):
    out = load_registry()['entries']
    if name is not None:
        out = [e for e in out if e['name'] == name]
    if size is not None:
        out = [e for e in out if e['size'] == size]
    return out


def format_table(es):
    rows = [('id', 'size', 'status', 'garbage', 'canaries (not proved/decided)', 'program hash')]
    for e in es:
        c = e['canaries']
        rows.append((e['id'], str(e['size']), e['status'], str(e['garbage']['wires']), f"{c['not_proved']}/{c['decided']}",
                     e['program']['value'][:12]))
    width = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    lines = ['  '.join(c.ljust(width[i]) for i, c in enumerate(r)).rstrip() for r in rows]
    lines.insert(1, '  '.join('-' * w for w in width))
    return '\n'.join(lines)


def format_entry(e):
    lines = [e['id'], f"  block:      {e['block']}", f"  claim:      {e['claim']}", f"  status:     {e['status']}",
             f"  projectq:   {e['projectq']}", f"  program:    {e['program']['name']} = {e['program']['value']}",
             f"  garbage:    {e['garbage']['wires']} wires, pin {e['garbage'].get('pin_sha256')}",
             f"  canaries:   {e['canaries']['not_proved']} of {e['canaries']['decided']} not proved "
             f"({e['canaries']['equivalent']} equivalent, groups {e['canaries']['by_group']})",
             f"  evidence:   {e['evidence']['file']} (sha256 {e['evidence']['sha256'][:16]})", f"  reproduce:  {e['reproduce']}"]
    if e.get('source'):
        lines.insert(4, f"  source:     {e['source']['file']} (sha256 {e['source']['sha256'][:16]})")
    lines.append('  assumptions:')
    lines += [f'    - {a}' for a in e['assumptions']]
    return '\n'.join(lines)


def verify_sources(results_dir):
    """Compare the evidence files of the registry in ``results_dir`` with the recorded SHA-256; returns a list of problems."""
    problems = []
    for e in entries():
        p = Path(results_dir) / e['evidence']['file']
        if not p.is_file():
            problems.append(f"{e['id']}: missing {p}")
        elif hashlib.sha256(p.read_bytes()).hexdigest() != e['evidence']['sha256']:
            problems.append(f"{e['id']}: {p.name} differs from the registry")
    return problems
