#!/usr/bin/env python3
"""Карта точки перехода из runs/sens_*.json -> markdown."""
import json, sys

for p in sys.argv[1:]:
    d = json.load(open(p)); m = d['model']
    print('\n### %s (сейчас: измерение %.0f нс + задержка %.0f нс = %.0f нс; T1 %.0f мкс, T2 %.0f мкс, CZ %.2f %%, считывание %.2f %%)\n' % (
        m['name'], m['meas_ns'], m['ff_ns'], m['meas_ns'] + m['ff_ns'], m['t1_us'], m['t2_us'], 100 * m['e2q'], 100 * m['ro']))
    grid = d['rows'][0]['grid']
    print('| блок | ошибка считывания | ' + ' | '.join('T x%g, CZ x%g' % tuple(g) for g in grid) + ' |')
    print('|---|---|' + '---|' * len(grid))
    for r in d['rows']:
        cells = ['уже выигрывает' if t is None else 'никогда' if t == 0 else '%.0f нс' % t for t in r['breakeven_ns']]
        print('| %s | x%g | %s |' % (r['block'], r['ro'], ' | '.join(cells)))
