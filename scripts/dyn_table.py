#!/usr/bin/env python3
"""Таблицы из runs/bench_*.json: python scripts/dyn_table.py runs/bench_FakeBoston.json ... > table.md"""
import json, sys

NAMES = {'U': 'обычное (U)', 'T': 'Тоффоли (T)', 'M': 'измерение (M)', 'Md': 'измерение, прямой cz (Md)', 'H': 'гибрид (H)'}


def main(paths):
    summ = []
    for p in paths:
        d = json.load(open(p)); rows = d['rows']
        be = rows[0]['backend']
        print('\n### %s\n' % be)
        print('| блок | вариант | кубиты | 2Q (+усл.) | измер. | длит., мкс | F оценка | F симул. Z | F симул. X |')
        print('|---|---|---|---|---|---|---|---|---|')
        order = []
        for r in rows:
            k = (r['family'], r['size'])
            if k not in order:
                order.append(k)
        for k in order:
            rs = [r for r in rows if (r['family'], r['size']) == k]
            best = rs[0]['best']
            for r in rs:
                z = '%.3f' % r['pZ'] if 'pZ' in r else '-'; x = '%.3f' % r['pX'] if 'pX' in r else '-'
                mark = ' **<- лучший**' if r['variant'] == best else ''
                print('| %s | %s%s | %d | %d (+%d) | %d | %.1f | %.3f | %s | %s |' % (
                    r['circuit'], r['variant'], mark, r['qubits'], r['n2q'], r['n2q_cond'], r['n_meas'],
                    r['duration_us'], r['f_est'], z, x))
            static = max(r['f_est'] for r in rs if r['variant'] in ('U', 'T'))
            dyn = [r for r in rs if r['variant'] not in ('U', 'T')]
            bd = max(dyn, key=lambda r: r['f_est']) if dyn else None
            summ.append((be, rs[0]['circuit'], best, static, bd['variant'] if bd else '-', bd['f_est'] if bd else 0.0))
        # does the simulation agree with the estimate about who wins (static vs dynamic), where the gap > 0.02?
        agree = tot = 0
        for k in order:
            rs = [r for r in rows if (r['family'], r['size']) == k and 'pZ' in r]
            st = [r for r in rs if r['variant'] in ('U', 'T')]; dy = [r for r in rs if r['variant'] not in ('U', 'T')]
            if st and dy:
                ge = max(r['f_est'] for r in st) - max(r['f_est'] for r in dy)
                gs = max(r['pZ'] for r in st) - max(r['pZ'] for r in dy)
                if abs(ge) > 0.02:
                    tot += 1; agree += (ge > 0) == (gs > 0)
        print('\nСовпадение победителя (статика/динамика) по оценке и по симуляции Z при разрыве > 0,02: %d из %d' % (agree, tot))
        # accuracy of the estimate
        errs = [abs(r['pZ'] - r['f_est']) for r in rows if 'pZ' in r and r['f_est'] > 0.05]
        if errs:
            print('\nОценка против симуляции (Z, F_est>0.05): среднее |разница| = %.3f, максимум = %.3f, n=%d' % (
                sum(errs) / len(errs), max(errs), len(errs)))
    print('\n### Итог: кто выигрывает\n')
    print('| бэкенд | блок | победитель | F лучшего статического | лучший динамический | F динамического |')
    print('|---|---|---|---|---|---|')
    for s in summ:
        print('| %s | %s | %s | %.3f | %s | %.3f |' % s)
    wins = [s for s in summ if s[2] not in ('U', 'T')]
    print('\nДинамический/гибрид победил в %d из %d случаев:' % (len(wins), len(summ)))
    for s in wins:
        print('- %s %s: %s (%.3f против %.3f)' % (s[0], s[1], s[2], s[5], s[3]))


if __name__ == '__main__':
    main(sys.argv[1:])
