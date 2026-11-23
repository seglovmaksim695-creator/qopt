"""Проводка ECAdd между шагами — берётся из настоящего составного блока Qualtran (cbloq.connections), а не
из рукописной схемы. Пишет results/wiring_ecadd_<версия>_n<n>.json; ew.wiring --merge сверяет все файлы и пишет
results/wiring_ecadd.json, который читает zglue (склейка E2/E3). Пример:
  python -m ew.wiring --n 16 --p 32749          (в окружении env_e.sh или env_orig.sh)
  python -m ew.wiring --merge"""
import argparse, json, os

NAMES = {'_ECAddStepOne': 1, '_ECAddStepTwo': 2, '_ECAddStepThree': 3, '_ECAddStepFour': 4,
         '_ECAddStepFive': 5, '_ECAddStepSix': 6}


def extract(n, p, window_size=1):
    import qualtran
    from qualtran import LeftDangle, RightDangle
    from qualtran.bloqs.cryptography.ecc.ec_add import ECAdd
    cb = ECAdd(n, p, window_size).decompose_bloq()
    kinds = {}
    for bi in cb.bloq_instances:
        nm = type(bi.bloq).__name__
        if nm not in NAMES:
            raise SystemExit(f'неожиданный блок в ECAdd: {nm}')
        k = NAMES[nm]
        if k in kinds.values():
            raise SystemExit(f'шаг {k} встречается дважды')
        kinds[bi] = k
    if sorted(kinds.values()) != [1, 2, 3, 4, 5, 6]:
        raise SystemExit(f'не все шаги: {sorted(kinds.values())}')

    def end(soq, side):
        if soq.idx != ():
            raise SystemExit(f'регистр с индексом {soq.idx} — схема рассчитана на целые регистры')
        if soq.binst is (LeftDangle if side == 'L' else RightDangle):
            return ['in' if side == 'L' else 'out', soq.reg.name]
        return [kinds[soq.binst], soq.reg.name]

    steps, out, used = {str(k): {} for k in range(1, 7)}, {}, set()
    for c in cb.connections:
        s, d = end(c.left, 'L'), end(c.right, 'R')
        if s[0] == 'out' or d[0] == 'in':
            raise SystemExit(f'странное соединение {s} -> {d}')
        if tuple(s) in used:
            raise SystemExit(f'выход {s} использован дважды')
        used.add(tuple(s))
        if d[0] == 'out':
            out[d[1]] = s
        else:
            if d[1] in steps[str(d[0])]:
                raise SystemExit(f'вход {d} подключён дважды')
            steps[str(d[0])][d[1]] = s
    # каждый шаг получает все свои входы, порядок шагов — топологический
    order = []
    for bi, k in sorted(kinds.items(), key=lambda t: t[1]):
        lefts = {r.name for r in bi.bloq.signature.lefts()}
        if lefts != set(steps[str(k)]):
            raise SystemExit(f'шаг {k}: входы {sorted(lefts)} != подключённые {sorted(steps[str(k)])}')
        for r, s in steps[str(k)].items():
            if s[0] != 'in' and s[0] >= k:
                raise SystemExit(f'шаг {k} берёт {r} у шага {s[0]} (не раньше)')
        order.append(k)
    rights = {r.name for r in ECAdd(n, p, window_size).signature.rights()}
    if rights != set(out):
        raise SystemExit(f'выходы ECAdd {sorted(rights)} != {sorted(out)}')
    from qoptct.foreign.qualtran._variant import detect_variant
    var = detect_variant()
    if var['kind'] not in ('original', 'e', 'e4'):
        raise SystemExit(f"unexpected Qualtran variant {var['kind']}")
    return dict(n=n, p=str(p), window_size=window_size, version='original' if var['kind'] == 'original' else 'patched', kind=var['kind'],
                qualtran_file=qualtran.__file__, order=order, steps=steps, out=out,
                connections=len(cb.connections))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', required=True)
    ap.add_argument('--out', default=None)
    a = ap.parse_args(argv)
    n, p = a.n, eval(a.p, {}, {})
    d = extract(n, p)
    d['status'] = 'OK'
    if a.out:
        json.dump(d, open(a.out, 'w'), ensure_ascii=False, indent=1)
    print(json.dumps(dict(version=d['version'], kind=d['kind'], connections=d['connections'], steps=d['steps'], out=d['out']),
                     ensure_ascii=False))


if __name__ == '__main__':
    main()
