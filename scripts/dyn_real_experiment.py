#!/usr/bin/env python
# SPDX-License-Identifier: Apache-2.0
"""One short real-hardware experiment: measurement-based erasure (dynamic circuits) vs reverse computation.

NOTHING is sent to IBM unless the command is `run` AND the flag --submit-to-real-hardware is given.

    plan  [--backend FakeBoston]            build and compile the circuits, print the plan and the estimates (default)
    run   --backend FakeBoston              run the same job locally (fake backend, Aer) to test the whole path
    run   --backend ibm_xxx --submit-to-real-hardware
                                            ONE Sampler job on the real device, max_execution_time = 60 s
    fetch dyn_real_<job id>.json            collect the result of a submitted job and print the table

Circuits (44 by default, 400 shots each, about 18 000 shots, a few seconds of QPU time): adder2, adder3, cmp2,
lookadd2x2, each in the static variant(s) (U: reverse computation; T: Toffoli adder) and the dynamic one(s)
(M: measurement + feed-forward; H: hybrid). Two tests per variant:
    Z: 3 basis inputs, success = every data wire right and every ancilla 0;
    X: |+> on all data wires, the block (a permutation of basis states keeps |+...+>), H, success = all zeros;
       this one is sensitive to phase errors, i.e. to a wrong or late fix-up.
Account: ~/.qiskit/quantum-start.json (name "quantum-start").
"""
import argparse
import json
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

PLAN = [('adder', 2, ('U', 'M*', 'T')), ('adder', 3, ('U', 'M*', 'T')), ('cmp', 2, ('U', 'M*')),
        ('lookadd', 2, ('U', 'H', 'M*'))]


def build_circuits(model, shots, n_z=3, seeds=(0, 1, 2, 3)):
    from qiskit import QuantumCircuit, ClassicalRegister
    from qoptct.dyn import choose, hw, lower
    out = []
    for fam, size, want in PLAN:
        vs = {lab: (B, e, vt) for lab, B, e, vt in choose.variants(fam, size)}
        if 'M*' in want:                                    # the better of the two fix-up styles, by estimate
            cand = {lab: choose.evaluate(*vs[lab], model, seeds)['sched'].fid for lab in ('M', 'Md')}
            vs['M*'] = vs[max(cand, key=cand.get)]
        for lab in want:
            B, erase, vt = vs[lab]
            assert set(B.xfree) == set(B.data), 'the hardware X test needs a permutation on all data wires'
            ok, why, _, _ = B.check()
            assert ok, why
            body = lower.to_qiskit(B.C, erase, vt); layout = None
            tests = [('Z', d) for d in B.inputs(n_z, seed=7)] + [('X', None)]
            for kind, d in tests:
                qc = QuantumCircuit(*body.qregs, *body.cregs, ClassicalRegister(B.C.N, 'out'))
                expect = 0
                if kind == 'Z':
                    e = B.fn(d)
                    for k, ws in B.regs.items():
                        for j, w in enumerate(ws):
                            if (d[k] >> j) & 1:
                                qc.x(w)
                            if (e[k] >> j) & 1:
                                expect |= 1 << w
                else:
                    for w in B.data:
                        qc.h(w)
                qc.compose(body, inplace=True)
                if kind == 'X':
                    for w in B.data:
                        qc.h(w)
                qc.measure(range(B.C.N), qc.cregs[-1])
                comp = hw.compile_circuit(qc, model, seeds=seeds if layout is None else seeds[:1],
                                          live_virtual=B.data, initial_layouts=[layout] if layout else [])
                layout = layout or comp['layout_init']
                st = comp['sched'].stats()
                out.append(dict(circuit=B.name, variant=lab, test=kind, input=d, expect=expect, tqc=comp['tqc'],
                                f_est=st['f_est'], n2q=st['n2q'], n2q_cond=st['n2q_cond'],
                                duration_us=st['duration_us'], n_wires=B.C.N))
    return out


def summarize(items, counts_list, shots):
    rows = {}
    for it, counts in zip(items, counts_list):
        good = counts.get(format(it['expect'], '0%db' % it['n_wires']), 0)
        r = rows.setdefault((it['circuit'], it['variant']), dict(Z=[], X=[], est=[], n2q=it['n2q'],
                                                                 cond=it['n2q_cond'], dur=it['duration_us']))
        r[it['test']].append(good / shots); r['est'].append(it['f_est'])
    print('%-12s %-3s %5s %5s %9s %8s %8s %8s' % ('circuit', 'var', '2Q', 'cond', 'T, us', 'F_est', 'Z', 'X'))
    for (c, v), r in rows.items():
        z = sum(r['Z']) / len(r['Z']); x = sum(r['X']) / len(r['X']); e = sum(r['est']) / len(r['est'])
        print('%-12s %-3s %5d %5d %9.2f %8.3f %8.3f %8.3f' % (c, v, r['n2q'], r['cond'], r['dur'], e, z, x))
    print('(F_est includes the final readout here; statistical error about +-%.3f per entry)'
          % (0.5 / (shots ** 0.5)))


def get_service():
    from qiskit_ibm_runtime import QiskitRuntimeService
    return QiskitRuntimeService(filename=str(Path.home() / '.qiskit' / 'quantum-start.json'), name='quantum-start')


def main():
    warnings.filterwarnings('ignore')
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', nargs='?', default='plan', choices=['plan', 'run', 'fetch'])
    p.add_argument('job', nargs='?')
    p.add_argument('--backend', default='FakeBoston')
    p.add_argument('--shots', type=int, default=400)
    p.add_argument('--submit-to-real-hardware', action='store_true')
    a = p.parse_args()
    from qoptct.dyn import hw

    if a.cmd == 'fetch':
        meta = json.load(open(a.job))
        job = get_service().job(meta['job_id'])
        print('status:', job.status())
        res = job.result()
        counts = [res[i].data.out.get_counts() for i in range(len(meta['items']))]
        summarize(meta['items'], counts, meta['shots'])
        try:
            print('QPU usage, s:', job.usage())
        except Exception:
            pass
        return

    real = not a.backend.startswith('Fake')
    if real and not (a.cmd == 'run' and a.submit_to_real_hardware):
        sys.exit('a real backend is only contacted by: run --backend %s --submit-to-real-hardware' % a.backend)
    if real:
        backend = get_service().backend(a.backend)
    else:
        import qiskit_ibm_runtime.fake_provider as fp
        backend = getattr(fp, a.backend)()
    if 'if_else' not in backend.target.operation_names:
        sys.exit('the backend does not support dynamic circuits (if_else)')
    model = hw.HwModel.from_backend(backend)
    items = build_circuits(model, a.shots)
    n = len(items); tot = n * a.shots
    longest = max(it['duration_us'] for it in items)
    print('backend %s: %d circuits x %d shots = %d shots; longest circuit %.1f us' % (backend.name, n, a.shots, tot,
                                                                                       longest))
    print('rough QPU time: %.0f s (0.3 ms per shot incl. the repetition delay) + job overhead; hard limit 60 s'
          % (tot * 3e-4))
    if a.cmd == 'plan':
        summarize(items, [{} for _ in items], a.shots)
        print('plan only: nothing was submitted (columns Z, X are empty).')
        return
    from qiskit_ibm_runtime import SamplerV2
    opts = {'max_execution_time': 60} if real else {}
    sampler = SamplerV2(mode=backend, options=opts)
    job = sampler.run([it['tqc'] for it in items], shots=a.shots)
    meta = dict(job_id=job.job_id(), backend=backend.name, shots=a.shots, time=time.strftime('%Y-%m-%d %H:%M:%S'),
                items=[{k: v for k, v in it.items() if k != 'tqc'} for it in items])
    if real:
        fn = 'dyn_real_%s.json' % job.job_id()
        json.dump(meta, open(fn, 'w'), indent=1)
        print('submitted job %s; later: python %s fetch %s' % (job.job_id(), sys.argv[0], fn))
        return
    res = job.result()
    summarize(items, [res[i].data.out.get_counts() for i in range(n)], a.shots)
    print('local run on the fake backend (Aer noise from its calibration, without idle errors): a test of the '
          'pipeline, not a measurement.')


if __name__ == '__main__':
    main()
