"""Независимая проверка схем OpenQASM 3 с измерениями и условными командами (подмножество). Мой код, не агента."""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[0])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import re, sys, numpy as np
sys.path[:0] = [_QPKG, (_QPKG + '/portfolio')]
import gidney, qopt, evo_lib
def parse(path):
    prog = []
    for L in open(path):
        L = L.strip()
        if not L or L.startswith(("OPENQASM", "include", "qubit", "bit")): continue
        m = re.match(r"c\[(\d+)\] = measure q\[(\d+)\];", L)
        if m: prog.append(("M", int(m.group(2)), int(m.group(1)))); continue
        m = re.match(r"if \((.*)\) \{(.*)\}$", L)
        if m:
            cond = m.group(1)
            assert re.fullmatch(r"[\s()c\[\]\d&]+", cond), "неизвестное условие: " + cond
            idx = [int(x) for x in re.findall(r"c\[(\d+)\]", cond)]
            body = [b.strip() for b in m.group(2).split(";") if b.strip()]
            prog.append(("IF", idx, [gate(b) for b in body])); continue
        prog.append(("G",) + gate(L.rstrip(";")))
    return prog
def gate(s):
    s = s.rstrip(";"); name, args = s.split(None, 1)
    return (name, tuple(int(x) for x in re.findall(r"q\[(\d+)\]", args)))
def run(prog, N, psi, rng, nbits):
    c = [0] * nbits
    for op in prog:
        if op[0] == "G": psi = gidney.apply(psi, N, op[1], op[2])
        elif op[0] == "IF":
            if all(c[i] for i in op[1]):
                for g, q in op[2]: psi = gidney.apply(psi, N, g, q)
        else:
            _, q, k = op; p = psi.reshape([2] * N)
            i1 = [slice(None)] * N; i1[q] = 1; p1 = float(np.sum(np.abs(p[tuple(i1)]) ** 2))
            m = int(rng.random() < p1); i0 = [slice(None)] * N; i0[q] = 1 - m
            p = p.copy(); p[tuple(i0)] = 0; psi = (p / np.sqrt(p1 if m else 1 - p1)).reshape(-1); c[k] = m
    return psi, c
def check(qasm, qc_path, trials=30, seed=1, strict=True, clean=False):
    """clean=True: провода .qc, не указанные в строке .i, — чистые анциллы (старт |0>), как задумано в бенчмарках."""
    names, hl = qopt.load_qc(qc_path); n = len(names); ref = evo_lib.lower_std(hl)
    prog = parse(qasm); N = int(re.search(r"qubit\[(\d+)\]", open(qasm).read()).group(1))
    _mb = re.search(r"(?<!qu)bit\[(\d+)\]", open(qasm).read()); nb = int(_mb.group(1)) if _mb else 0
    rng = np.random.default_rng(seed); meas_q = [op[1] for op in prog if op[0] == "M"]; worst = 1
    for _ in range(trials):
        v = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n)
        if clean:
            ins = gidney.inputs_of(qc_path)
            anc = [i for i, nm in enumerate(names) if ins is not None and nm not in ins]
            vv = v.reshape([2] * n)
            for a in anc:
                sl = [slice(None)] * n; sl[a] = 1; vv[tuple(sl)] = 0
            v = vv.reshape(-1)
        v /= np.linalg.norm(v)
        want = v.copy()
        for g, q in ref: want = gidney.apply(want, n, g, q)
        psi = np.zeros([2] * N, complex); psi[(slice(None),) * n + (0,) * (N - n)] = v.reshape([2] * n); psi = psi.reshape(-1)
        out, c = run(prog, N, psi, rng, nb)
        o = out.reshape([2] * N); idx = [slice(None)] * N
        for q in range(n, N): idx[q] = int(np.argmax([np.sum(np.abs(np.take(o, b, axis=q)) ** 2) for b in (0, 1)]))
        got = o[tuple(idx)].reshape(-1)
        if strict:  # без перенормировки: анциллы обязаны быть в базисном состоянии (не запутаны с данными)
            worst = min(worst, abs(np.vdot(want, got)) / np.linalg.norm(out))
        else:
            got = got / np.linalg.norm(got); worst = min(worst, abs(np.vdot(want, got)))
    return worst
if __name__ == "__main__":
    print(sys.argv[1], "худшая точность", check(sys.argv[1], sys.argv[2]))
