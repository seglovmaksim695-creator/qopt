"""Пакет склеек Z3: E2 (исправленная), канарейки склейки (без каждой леммы; исходные шаги), E3, E3fail."""
import json, subprocess, sys, os
from .zglue import LEMMAS

def run(args):
    r = subprocess.run([sys.executable, '-B', '-s', '-m', 'qoptct.foreign.qualtran.ecadd.zglue'] + args, capture_output=True, text=True)
    try:
        return json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        return dict(status='ERROR', stderr=(r.stderr or '')[-300:])

def main():
    out = sys.argv[1]
    rows = []
    for n, p in ((16, '32749'), (256, '2**255-19')):
        base = ['--n', str(n), '--p', p]
        jobs = [('E2', ['--mode', 'E2'])]
        jobs += [(f'E2-drop-{L}', ['--mode', 'E2', '--drop', L]) for L in LEMMAS]
        jobs += [('E2-orig126', ['--mode', 'E2', '--orig', '1,2,6']), ('E2-orig1', ['--mode', 'E2', '--orig', '1']),
                 ('E2-orig2', ['--mode', 'E2', '--orig', '2']), ('E2-orig6', ['--mode', 'E2', '--orig', '6'])]
        jobs += [('E3', ['--mode', 'E3', '--orig', '1,2,6'])]
        for g in (1, 2, 3):
            jobs += [(f'E3fail-G{g}-orig', ['--mode', 'E3fail', '--orig', '1,2,6', '--group', str(g)]),
                     (f'E3fail-G{g}-patched', ['--mode', 'E3fail', '--group', str(g)])]
        from .zmutants import MUTANTS
        jobs += [(f'E2-mutant-{m}', ['--mode', 'E2', '--mutant', m]) for m in MUTANTS]
        from .zglue import WMUT
        jobs += [(f'E2-wiring-{m}', ['--mode', 'E2', '--wmut', m]) for m in WMUT]
        for name, extra in jobs:
            d = run(base + extra)
            row = dict(n=n, job=name, status=d.get('status'), answers=d.get('answers'), failed=d.get('failed', [])[:3],
                       z3_seconds=d.get('z3_seconds'), err=d.get('stderr', '')[-200:] if d.get('status') == 'ERROR' else '')
            rows.append(row)
            print(json.dumps(row, ensure_ascii=False), flush=True)
    json.dump(rows, open(out, 'w'), ensure_ascii=False, indent=1)

if __name__ == '__main__':
    main()
