"""c2: improve(path, qasm_out, clean=False) — логическое И + стирание измерением (см. andc.py, REPORT.md).
clean=False (строгая модель): все провода .qc — данные в любом состоянии, |0> только у новых кубитов.
clean=True: провода .qc не из строки .i стартуют в |0> (verify_qasm.check(..., clean=True)); тогда компилятор
знает их значения с начала (И прямо в чистый провод, стирание провода измерением), а строгие кандидаты
тоже годятся (схема, верная на всех входах, верна и на чистых) — берётся лучшее из обоих.
Перебор: перестановки коммутирующих логических вентилей × ленивость × лимит анцилл (по 8 с на компиляцию);
минимум T (затем CX) среди N <= 22 (или max_extra). python run.py схема.qc out.qasm [max_extra] [clean]"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, random, itertools, signal
sys.path[:0] = [_QPKG, (_QPKG + '/contrib/c2')]
import andc, qoptct, sweep2, verify_qasm


def candidates(path, variants=20, max_extra=None, clean=False):
    names, hl = qoptct.load_qc(path); n = len(names)
    lim = 22 - n if max_extra is None else min(max_extra, 22 - n)
    lg0 = andc.logical(hl, n); rng = random.Random(1); out = {}
    modes = [None] + ([andc.inputs_idx(path, names)] if clean else [])
    for v in range(variants):
        lg = lg0 if v == 0 else sweep2.shuffle(lg0, rng)
        for ins, lazy, lm, ma in itertools.product(modes, (True, False, 0.2, 0.4, 0.6), (3, 1), sorted({lim, 4, 3, 2, 1, 0})):
            if ma > lim or ma < 0: continue
            signal.alarm(8)
            try:
                C = andc.Compiler(n, lg, ma, random.Random(v), lazy=lazy, linmax=lm, inputs=ins)
                ops, N = andc.compact(C.run(), n)
            except BaseException:
                signal.alarm(0); continue
            signal.alarm(0)
            for fold in (False, True):
                txt = andc.to_qasm(ops, N, fold=fold); T, CX = andc.count(txt)
                k = (T, CX, N - n)
                if k not in out: out[k] = txt
    return out, n


def improve(path, qasm_out, max_extra=None, variants=120, trials=30, objective='T', clean=False):
    """objective: 'T' — минимум T, затем CX; 'CX' — минимум CX, затем T."""
    cands, n = candidates(path, variants, max_extra, clean)
    if not cands: return None
    key = min(cands, key=(lambda k: (k[0], k[1], k[2])) if objective == 'T' else (lambda k: (k[1], k[0], k[2])))
    open(qasm_out, 'w').write(cands[key])
    F = verify_qasm.check(qasm_out, path, trials=trials, clean=clean) if n + key[2] <= 22 else None
    return dict(T=key[0], CX=key[1], extra_qubits=key[2], fidelity=None if F is None else float(F))


if __name__ == '__main__':
    me = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != 'clean' else None
    print(improve(sys.argv[1], sys.argv[2], max_extra=me, clean='clean' in sys.argv[3:]))
