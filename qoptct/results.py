# SPDX-License-Identifier: Apache-2.0
"""``qoptct results``: one table of everything that is proved, over all registries of the package.

    qoptct results [--library qualtran|qiskit|projectq|dyadic|auto|loop] [--json] [--long]

The registries are the data files of the installed package (``foreign/qualtran``, ``foreign/qiskit``, ``foreign/projectq``,
``foreign/auto``, ``dyadic``, ``loop``); nothing is recomputed and no optional dependency is needed.  A row says what was
proved, up to which size, and with how many canaries (wrong variants that must NOT be proved); the command in the
``reproduce`` column of ``--json`` re-runs it.  A result holds for the stated size and the stated assumptions only.

Only proved results are listed (status PROVED, or PROVED_CONDITIONAL / PROVED_OUTSIDE_G1_G3, which say what is assumed);
``--all`` adds the rows of the registries that are not proved (for example NOT_PROVED sizes recorded by ``auto``).
Exit code: 0 when the table was printed, 1 when nothing matched.
"""
from __future__ import annotations

import argparse
import json
import sys

LIBRARIES = ('qualtran', 'qiskit', 'projectq', 'dyadic', 'auto', 'loop')
_OK = ('PROVED', 'PROVED_CONDITIONAL', 'PROVED_OUTSIDE_G1_G3')


def _row(library, name, size, status, canaries, assumptions, reproduce, detail=''):
    return dict(library=library, name=name, size=size, status=status, canaries=canaries, assumptions=assumptions,
                reproduce=reproduce, detail=detail)


def _qualtran():
    from .foreign.qualtran import registry
    out = []
    for e in registry.entries():
        size = 'n=%s' % e['n'] + (', w=%s' % e['w'] if e.get('w') else '')
        out.append(_row('qualtran', e['name'], size, e['status'], '-', len(e.get('assumptions') or []), e.get('reproduce', ''),
                        e.get('block', '')))
    return out


def _qiskit():
    import json as _j
    from importlib import resources
    data = _j.loads(resources.files('qoptct.foreign.qiskit').joinpath('data/registry.json').read_text(encoding='utf8'))
    out = []
    for r in data['rows']:
        c = r.get('canaries')
        cs = '%d/%d' % (c['refuted'] + c['equivalent'], c['mutants']) if c else ('none' if r.get('canaries_note') else '-')
        if c and r.get('canaries_note'):
            cs += ' (at n=%s)' % c['n']
        out.append(_row('qiskit', r['family'], 'n=%s' % r['max_n_proved'], r['status'], cs, 1, r.get('command', ''),
                        'Qiskit %s%s' % (data.get('qiskit', ''), ('; ' + r['canaries_note']) if r.get('canaries_note') else '')))
    return out


def _projectq():
    from .foreign.projectq import registry
    out = []
    for e in registry.entries():
        c = e['canaries']
        out.append(_row('projectq', e['id'], 'm/n=%s' % e['size'], e['status'], '%d/%d' % (c['not_proved'], c['decided']),
                        len(e.get('assumptions') or []), e.get('reproduce', ''), e.get('block', '')))
    return out


def _dyadic():
    from .dyadic import registry
    out = []
    for e in registry.entries():
        c = e['canaries']
        out.append(_row('dyadic', e['name'], 'n=%s' % e['n'], e['status'], '%d/%d' % (c['refuted'], c['total']),
                        len(e.get('assumptions') or []), e.get('reproduce', ''), e.get('block', '')))
    return out


def _auto():
    from .foreign.auto import results as registry
    out = []
    for e in registry.entries():
        c = e.get('canaries') or {}
        cs = '%d/%d' % (c.get('total', 0) - len(c.get('proved') or []), c['total']) if c.get('total') else '-'
        out.append(_row('auto', e['id'], 'n=%s' % e['n'] if e.get('n') else '-', e['status'], cs,
                        len(e.get('assumptions') or []), e.get('reproduce', ''), e.get('block', '')))
    return out


def _loop():
    from importlib import resources
    data = json.loads(resources.files('qoptct.loop').joinpath('data/registry.json').read_text(encoding='utf8'))
    out = []
    for e in data['entries']:
        b = e['best']
        out.append(_row('loop', e['target'], 'AND %d, CX %d' % (b['AND'], b['CX']), e['status'], '>=20 per record', 1,
                        e['reproduce'], '%d records%s' % (e['records'], ', improves the baseline' if e['improves_baseline'] else '')))
    return out


_LOADERS = dict(qualtran=_qualtran, qiskit=_qiskit, projectq=_projectq, dyadic=_dyadic, auto=_auto, loop=_loop)


def collect(libraries=LIBRARIES):
    rows = []
    for lib in libraries:
        rows += _LOADERS[lib]()
    return rows


def format_table(rows):
    head = ('library', 'result', 'size', 'status', 'canaries', 'assumptions')
    body = [(r['library'], r['name'], r['size'], r['status'], r['canaries'], str(r['assumptions'])) for r in rows]
    width = [max(len(x[i]) for x in [head] + body) for i in range(len(head))]
    lines = ['  '.join(c.ljust(width[i]) for i, c in enumerate(x)).rstrip() for x in [head] + body]
    lines.insert(1, '  '.join('-' * w for w in width))
    return '\n'.join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='qoptct results', description='everything that is proved, over all registries')
    ap.add_argument('--library', action='append', choices=LIBRARIES, help='only this library (repeatable)')
    ap.add_argument('--all', action='store_true', help='also list rows that are not proved')
    ap.add_argument('--json', action='store_true')
    ap.add_argument('--long', action='store_true', help='one block per row with the command that reproduces it')
    a = ap.parse_args(argv)
    rows = collect(tuple(a.library) if a.library else LIBRARIES)
    omitted = 0
    if not a.all:
        omitted = sum(1 for r in rows if r['status'] not in _OK)
        rows = [r for r in rows if r['status'] in _OK]
    if a.json:
        print(json.dumps(rows, indent=1))
    elif a.long:
        for r in rows:
            print('%s / %s\n  size: %s   status: %s   canaries: %s   assumptions: %s\n  %s\n  reproduce: %s\n'
                  % (r['library'], r['name'], r['size'], r['status'], r['canaries'], r['assumptions'], r['detail'],
                     r['reproduce']))
    else:
        print(format_table(rows))
        libs = sorted({r['library'] for r in rows})
        print('\n%d results in %d libraries (%s).  Each holds for the stated size and assumptions only; '
              '`qoptct results --long` shows the command that reproduces it.' % (len(rows), len(libs), ', '.join(libs)))
        if omitted:
            print('%d registry rows are not proved and not listed (--all shows them).' % omitted)
    return 0 if rows else 1


if __name__ == '__main__':
    sys.exit(main())
