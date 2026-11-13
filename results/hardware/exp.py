"""Опыт: исходная схема (A) vs qopt унитарная (B) vs qopt с измерениями (C) на реальном IBM.
Этапы: prep (сборка, транспиляция, шумовая симуляция) -> submit (ОДНО задание) -> fetch.
Запуск: .venv/bin/python agents/ibm/exp.py prep|submit|fetch
"""
import itertools, json, sys
from pathlib import Path
import numpy as np
sys.path[:0] = ['..', '../cloud_v2/evolve']
import qopt, evo_lib
from qiskit import QuantumCircuit, ClassicalRegister, qasm3, qpy
from qiskit.quantum_info import Statevector
from qiskit.transpiler import generate_preset_pass_manager

R = Path('../release_qopt')
OUT = Path('../agents/ibm')
SCHEMES = ['tof_3', 'barenco_tof_3', 'mod5_4', 'tof_4', 'barenco_tof_4']
N_IN = 4
SHOTS = 4000
SEED = 7


def service():
    from qiskit_ibm_runtime import QiskitRuntimeService
    return QiskitRuntimeService(filename=str(Path.home() / '.qiskit/quantum-start.json'), name='quantum-start')


def variant_A(nm):
    names, hl = qopt.load_qc(R / 'bench_qc' / f'{nm}.qc')
    n = len(names)
    qc = QuantumCircuit(n, name=f'{nm}_A')
    for g, q in evo_lib.lower_std(hl):
        getattr(qc, g)(*q)
    return n, qc


def variant_qasm(nm, kind):
    return qasm3.loads(open(R / 'results/best' / kind / f'{nm}.qasm').read())


def perm_out(qcA, n, x):
    """Классический выход исходной схемы для базисного входа x (список битов)."""
    c = QuantumCircuit(n)
    for i, b in enumerate(x):
        if b: c.x(i)
    c.compose(qcA, inplace=True)
    p = Statevector(c).probabilities()
    k = int(np.argmax(p)); assert p[k] > 0.999
    return [(k >> i) & 1 for i in range(n)]


def with_io(core, n, x):
    nq = core.num_qubits
    c = QuantumCircuit(*core.qregs, *core.cregs) if core.qregs else QuantumCircuit(nq)
    out = ClassicalRegister(n, 'out'); c.add_register(out)
    for i, b in enumerate(x):
        if b: c.x(i)
    c.compose(core, qubits=range(nq), clbits=range(core.num_clbits), inplace=True)
    c.measure(range(n), out)
    return c


def n2q(c):
    cnt = 0
    for inst in c.data:
        op = inst.operation
        if op.name == 'if_else':
            for b in op.blocks:
                if b is not None: cnt += n2q(b)
        elif op.num_qubits == 2:
            cnt += 1
    return cnt


def prep():
    s = service()
    bk = s.least_busy(operational=True, simulator=False, dynamic_circuits=True, min_num_qubits=10)
    print('backend', bk.name, flush=True)
    rng = np.random.default_rng(2026)
    pm = lambda lay: generate_preset_pass_manager(backend=bk, optimization_level=1, seed_transpiler=SEED, initial_layout=lay)
    cmap = bk.coupling_map
    items, circs = [], []
    for nm in SCHEMES:
        n, A = variant_A(nm)
        B = variant_qasm(nm, 'unitary'); C = variant_qasm(nm, 'strict')
        # входы: половина «срабатывающих» (выход != вход), половина нет; все провода .qc — данные
        allx = [[(k >> i) & 1 for i in range(n)] for k in range(2 ** n)]
        fire = [x for x in allx if perm_out(A, n, x) != x]; idle = [x for x in allx if perm_out(A, n, x) == x]
        xs = [fire[i] for i in rng.choice(len(fire), N_IN // 2, replace=False)] + \
             [idle[i] for i in rng.choice(len(idle), N_IN - N_IN // 2, replace=False)]
        # раскладка: свободно для исходной A; её же для данных B и C
        tA = pm(None).run(with_io(A, n, xs[0]))
        init = tA.layout.initial_index_layout(filter_ancillas=True)[:n]
        anc_n = C.num_qubits - n
        lay_C = None
        if anc_n:
            near = sorted({nb for p in init for nb in cmap.neighbors(p)} - set(init))
            near2 = near + sorted({nb for p in near for nb in cmap.neighbors(p)} - set(init) - set(near))
            near2 = near2[:9]
            best = None
            for combo in itertools.permutations(near2, anc_n):
                if len(set(combo)) < anc_n: continue
                t = pm(init + list(combo)).run(with_io(C, n, xs[0]))
                k = n2q(t)
                if best is None or k < best[0]: best = (k, list(combo))
            lay_C = init + best[1]
        else:
            lay_C = init
        for var, core, lay in [('A', A, init), ('B', B, init), ('C', C, lay_C)]:
            for x in xs:
                y = perm_out(A, n, x)
                t = pm(lay).run(with_io(core, n, x))
                items.append({'scheme': nm, 'variant': var, 'n': n, 'qubits_logical': core.num_qubits,
                              'x': x, 'expect': ''.join(str(b) for b in reversed(y)),
                              'cx_logical': n2q(core), 'twoq': n2q(t), 'depth': t.depth(),
                              'layout': lay})
                circs.append(t)
            print(nm, var, 'logical2q', n2q(core), '2q', items[-1]['twoq'], 'depth', items[-1]['depth'], flush=True)
    with open(OUT / 'circuits.qpy', 'wb') as f:
        qpy.dump(circs, f)
    # шумовая симуляция с моделью шума этого устройства
    from qiskit_aer import AerSimulator
    from qiskit_aer.primitives import SamplerV2 as AerSampler
    sim = AerSimulator.from_backend(bk)
    sampler = AerSampler(options={'backend_options': {'noise_model': sim.options.noise_model}, 'run_options': {'seed_simulator': SEED}})
    r = sampler.run(circs, shots=SHOTS).result()
    for it, pr in zip(items, r):
        cnt = pr.data.out.get_counts()
        it['sim_counts'] = cnt; it['sim_p'] = cnt.get(it['expect'], 0) / SHOTS
    json.dump({'backend': bk.name, 'shots': SHOTS, 'items': items}, open(OUT / 'prep.json', 'w'), indent=1)
    summary(items, 'sim_p')


def summary(items, key):
    from collections import defaultdict
    g = defaultdict(list)
    for it in items:
        if key in it: g[(it['scheme'], it['variant'])].append(it)
    for (nm, v), L in g.items():
        print(f"{nm:15s} {v} 2q={L[0]['twoq']:3d} depth={L[0]['depth']:4d} {key}={np.mean([i[key] for i in L]):.3f}", flush=True)


def submit():
    st = OUT / 'job.json'
    if st.exists():
        print('уже отправлено', st.read_text()); return
    P = json.load(open(OUT / 'prep.json'))
    circs = qpy.load(open(OUT / 'circuits.qpy', 'rb'))
    s = service(); bk = s.backend(P['backend'])
    u = s.usage(); print('остаток секунд', u.get('usage_remaining_seconds'), flush=True)
    from qiskit_ibm_runtime import SamplerV2
    job = SamplerV2(mode=bk, options={'max_execution_time': 230}).run(circs, shots=SHOTS)
    st.write_text(json.dumps({'job_id': job.job_id(), 'backend': bk.name, 'n_circuits': len(circs), 'shots': SHOTS}, indent=1))
    print('job', job.job_id(), flush=True)


def fetch():
    J = json.load(open(OUT / 'job.json')); P = json.load(open(OUT / 'prep.json'))
    s = service(); job = s.job(J['job_id'])
    stt = str(job.status()); print('status', stt, flush=True)
    if stt.upper() != 'DONE':
        print(job.error_message() if 'ERROR' in stt.upper() else ''); return
    r = job.result()
    for it, pr in zip(P['items'], r):
        cnt = pr.data.out.get_counts(); it['real_counts'] = cnt; it['real_p'] = cnt.get(it['expect'], 0) / SHOTS
    m = job.metrics()
    P['job_id'] = J['job_id']; P['metrics'] = m
    try: P['usage_seconds'] = job.usage()
    except Exception as e: P['usage_seconds'] = repr(e)
    json.dump(P, open(OUT / 'results.json', 'w'), indent=1, default=str)
    print('usage', P['usage_seconds'], json.dumps(m.get('timestamps', {}), default=str))
    summary(P['items'], 'real_p')


if __name__ == '__main__':
    {'prep': prep, 'submit': submit, 'fetch': fetch}[sys.argv[1]]()
