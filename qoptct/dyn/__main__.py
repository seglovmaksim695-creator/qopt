# SPDX-License-Identifier: Apache-2.0
"""CLI of qoptct.dyn (needs qiskit-ibm-runtime for the fake backends and qiskit-aer for --sim):

    python -m qoptct.dyn choose adder 4 --backend FakeBoston [--sim] [--qasm best.qasm]
    python -m qoptct.dyn bench --backend FakeBoston --families adder,cmp --max-qubits 16 --sim --out res.json
    python -m qoptct.dyn threshold adder 4 --backend FakeBoston
    python -m qoptct.dyn convert in.qasm out.qasm --backend FakeBoston --mode auto
"""
import argparse
import json
import sys
import time
import warnings


def _fmt(row):
    s = '%-14s %-3s q=%-2d 2Q=%-4d(+%-3d cond) meas=%-2d depth=%-4d T=%7.2f us  F_est=%.3f' % (
        row['circuit'], row['variant'], row['qubits'], row['n2q'], row['n2q_cond'], row['n_meas'], row['depth'],
        row['duration_us'], row['f_est'])
    if 'pZ' in row:
        s += '  sim: Z=%.3f+-%.3f' % (row['pZ'], row['pZ_err'])
    if 'pX' in row:
        s += ' X=%.3f+-%.3f' % (row['pX'], row['pX_err'])
    return s


def _model(a):
    from . import choose, hw
    return hw.HwModel.from_backend(choose.resolve_backend(a.backend), ff_latency=a.ff_ns * 1e-9, fast_meas=a.fast_meas)


def cmd_choose(a):
    from . import choose
    M = _model(a)
    r = choose.choose(a.family, a.size, M, seeds=tuple(range(a.seeds)), simulate=a.sim, shots=a.shots,
                      n_inputs=a.inputs)
    print('model:', json.dumps(r['model']))
    for row in r['rows']:
        print(_fmt(row))
    print('best by estimate:', r['best'])
    if a.qasm:
        from qiskit import qasm3
        with open(a.qasm, 'w') as f:
            f.write(qasm3.dumps(r['comps'][r['best']]['logical']))
        print('wrote', a.qasm)


def cmd_bench(a):
    from . import choose, gen
    M = _model(a); rows = []; t0 = time.time()
    fams = a.families.split(',')
    for fam in fams:
        for size in gen.FAMILIES[fam]:
            if max(B.n for B in gen.build(fam, size)) > a.max_qubits:
                continue
            t = time.time()
            r = choose.choose(fam, size, M, seeds=tuple(range(a.seeds)), simulate=a.sim, shots=a.shots,
                              n_inputs=a.inputs, sim_max_qubits=a.sim_max_qubits)
            for row in r['rows']:
                row.update(family=fam, size=size, backend=M.name, best=r['best']); rows.append(row)
                print(_fmt(row), flush=True)
            print('  -> best %s  (%.0f s)' % (r['best'], time.time() - t), flush=True)
            if a.out:
                with open(a.out, 'w') as f:
                    json.dump(dict(model=r['model'], rows=rows, elapsed=time.time() - t0), f, indent=1)


def cmd_threshold(a):
    from . import choose
    M = _model(a)
    r = choose.choose(a.family, a.size, M, seeds=tuple(range(a.seeds)))
    s = r['model']; t_chain = s['meas_ns'] + s['ff_ns']
    print('model:', json.dumps(s))
    print('estimated fidelity of the best static (U/T) and the best dynamic (M/Md/H) variant when the measurement '
          'chain is improved')
    print('%-22s' % 'readout error \\ meas+ff', ''.join('%14s' % ('%.0f ns' % t) for t in a.times))
    for ro in a.ro:
        line = '%-22s' % ('%.2g' % ro)
        for t in a.times:
            f = t / t_chain
            M2 = M.scaled(ro=ro / s['ro'], meas_dur=f, ff_latency=M.ff_latency * f)
            sc = choose.rescore(r, M2)
            st = max(v for k, v in sc.items() if k in ('U', 'T')); dy = max(v for k, v in sc.items() if k not in ('U', 'T'))
            line += '%14s' % ('%.3f/%.3f%s' % (st, dy, '*' if dy > st else ' '))
        print(line)
    print('(static/dynamic; * = the dynamic variant wins)')


def cmd_convert(a):
    from qiskit import qasm3
    from . import choose
    qc = qasm3.load(a.inp)
    out, rep = choose.choose_circuit(qc, _model(a), seeds=tuple(range(a.seeds)), mode=a.mode)
    with open(a.out, 'w') as f:
        f.write(qasm3.dumps(out))
    print(json.dumps({k: v for k, v in rep.items() if not k.startswith('stats')}))


def main(argv=None):
    warnings.filterwarnings('ignore')
    p = argparse.ArgumentParser(prog='python -m qoptct.dyn', description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    def common(s):
        s.add_argument('--backend', default='FakeBoston'); s.add_argument('--ff-ns', type=float, default=600.0)
        s.add_argument('--seeds', type=int, default=4)
        s.add_argument('--fast-meas', action='store_true', help='model mid-circuit measurements with measure_2')
    s = sub.add_parser('choose'); s.add_argument('family'); s.add_argument('size', type=int); common(s)
    s.add_argument('--sim', action='store_true'); s.add_argument('--shots', type=int, default=300)
    s.add_argument('--inputs', type=int, default=8); s.add_argument('--qasm'); s.set_defaults(fn=cmd_choose)
    s = sub.add_parser('bench'); common(s)
    s.add_argument('--families', default='adder,cmp,modadd,lookup,lookadd,mult')
    s.add_argument('--max-qubits', type=int, default=25); s.add_argument('--sim', action='store_true')
    s.add_argument('--sim-max-qubits', type=int, default=22)
    s.add_argument('--shots', type=int, default=300); s.add_argument('--inputs', type=int, default=8)
    s.add_argument('--out'); s.set_defaults(fn=cmd_bench)
    s = sub.add_parser('threshold'); s.add_argument('family'); s.add_argument('size', type=int); common(s)
    s.add_argument('--times', type=float, nargs='+', default=[2800, 1000, 300, 100, 30])
    s.add_argument('--ro', type=float, nargs='+', default=[5e-3, 2e-3, 1e-3, 3e-4, 1e-4])
    s.set_defaults(fn=cmd_threshold)
    s = sub.add_parser('convert'); s.add_argument('inp'); s.add_argument('out'); common(s)
    s.add_argument('--mode', default='auto', choices=['auto', 'rule', 'unitary', 'keep']); s.set_defaults(fn=cmd_convert)
    a = p.parse_args(argv)
    a.fn(a)


if __name__ == '__main__':
    main()
