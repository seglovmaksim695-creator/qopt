"""Поиск: случайные перестановки коммутирующих логических вентилей x параметры компилятора x лимит анцилл.
python sweep2.py имена число_перестановок"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, json, os, itertools, random, signal
class _TO(Exception): pass
def _h(*a): raise _TO()
signal.signal(signal.SIGALRM, _h)
sys.path.insert(0, (_QPKG + '/contrib/c2'))
import andc, qopt
B = "tof_3 tof_4 tof_5 barenco_tof_3 barenco_tof_4 mod5_4 vbe_adder_3 hwb6 mod_red_21".split()
OUT = (_QPKG + '/contrib/c2/out'); os.makedirs(OUT, exist_ok=True)


def supp(g):
    nm, q = g
    if nm in ('cx', 'ccx'): return set(q[:-1]), {q[-1]}
    if nm == 'x': return set(), set(q)
    return set(q), set()   # фазовые: «читают» все


def commute(g, h):
    rg, wg = supp(g); rh, wh = supp(h)
    if g[0] == 'x' and h[0] == 'x': return True
    if wg & (rh | wh) - (wg & wh) or wh & (rg | wg) - (wg & wh): return False
    if wg & wh and not (g[0] in ('cx', 'ccx', 'x') and h[0] in ('cx', 'ccx', 'x')): return False
    return not (wg & rh) and not (wh & rg)


def shuffle(lg, rng, steps=None):
    lg = list(lg); steps = steps or 4 * len(lg)
    for _ in range(steps):
        i = rng.randrange(len(lg) - 1)
        if commute(lg[i], lg[i + 1]): lg[i], lg[i + 1] = lg[i + 1], lg[i]
    return lg


def search(nm, NV, log=True, clean=False, out=None):
    names_, hl = qopt.load_qc(f'{_QPKG}/../bench_qc/{nm}.qc'); n = len(names_)
    lg0 = andc.logical(hl, n); rng = random.Random(1)
    ins = andc.inputs_idx(f'{_QPKG}/../bench_qc/{nm}.qc', names_) if clean else None
    pts = {}
    for v in range(NV):
        lg = lg0 if v == 0 else shuffle(lg0, rng)
        for lazy, lm, ma in itertools.product((True, False, 0.2, 0.4, 0.6), (3, 1), (22 - n, 4, 3, 2, 1, 0)):
            if ma > 22 - n: continue
            signal.alarm(8)
            try:
                C = andc.Compiler(n, lg, ma, random.Random(v), lazy=lazy, linmax=lm, inputs=ins)
                ops = C.run(); ops, N = andc.compact(ops, n)
            except BaseException as e:
                signal.alarm(0)
                continue
            signal.alarm(0)
            for fold in (False, True):
                txt = andc.to_qasm(ops, N, fold=fold); T, CX = andc.count(txt); key = (T, CX, N - n)
                if key not in pts: pts[key] = ((v, lazy, lm, ma, fold), txt)
    keys = sorted(pts)
    par = [k for k in keys if not any(o != k and o[0] <= k[0] and o[1] <= k[1] and o[2] <= k[2] for o in keys)]
    od = out or (OUT + ('_clean' if clean else '')); os.makedirs(od, exist_ok=True)
    for k in par: open(f'{od}/{nm}_T{k[0]}_CX{k[1]}_A{k[2]}.qasm', 'w').write(pts[k][1])
    if log: print(nm, 'pareto', [(k, pts[k][0]) for k in par], flush=True)
    return par


if __name__ == '__main__':
    names = sys.argv[1].split(',') if len(sys.argv) > 1 else B
    NV = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    clean = len(sys.argv) > 3 and sys.argv[3] == 'clean'
    res = {nm: search(nm, NV, clean=clean) for nm in names}
    json.dump(res, open(f'sweep2_{names[0]}_{len(names)}{"_clean" if clean else ""}.json', 'w'))
