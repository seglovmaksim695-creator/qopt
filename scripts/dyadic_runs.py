#!/usr/bin/env python
"""Registry runs of qoptct.dyadic on library circuits with rotations pi/2^k (Qiskit; ProjectQ from JSON dumps).

    python scripts/dyadic_runs.py run FAMILY N [--canaries] [--budget S]     # one result as a JSON line
    python scripts/dyadic_runs.py registry [--families f1,f2] [--sizes 2,4,8] [--out results/dyadic]

Families (Qiskit 2.x): qft, qft_noswap, iqft, qft_qftdg, draper_fixed, draper_half, cdraper (DraperQFTAdder
.control(1)), rgqft, hrs;
ProjectQ (needs dumps from scripts/dyadic_projectq_dump.py in results/dyadic/projectq/): pq_addconst, pq_qft,
pq_addconst_basic, pq_qft_basic (decomposed further by ProjectQ), pq_addconst_opt (with ProjectQ's LocalOptimizer).
Each run: build, translate to a DCircuit, check the translation against Qiskit's own Operator (<= 10 qubits), prove,
and with --canaries prove >= 20 wrong variants (every one must be REFUTED).  `registry` runs every (family, n) in its
own process (the peak memory is that of one proof) and writes <out>/<family>.json.
"""
import argparse
import json
import os
import random
import resource
import subprocess
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from qoptct.dyadic import (DCircuit, Pi, prove_equal, prove_spec, QFTSpec, adder_spec, const_adder_spec,  # noqa
                                  multiplier_spec, ClassicalSpec)

BASIS = ['h', 'x', 'y', 'z', 's', 'sdg', 't', 'tdg', 'sx', 'sxdg', 'cx', 'cz', 'cy', 'ccx', 'ccz', 'p', 'cp', 'rz',
         'rx', 'ry', 'u', 'swap', 'cswap', 'crz', 'cs', 'csdg']
SIZES = {'qft': [2, 4, 8, 16, 32, 64, 128, 256], 'qft_noswap': [2, 4, 8, 16, 32, 64, 128, 256],
         'iqft': [2, 4, 8, 16, 32, 64, 128, 256], 'qft_qftdg': [2, 4, 8, 16, 32, 64, 128, 256],
         'draper_fixed': [1, 2, 4, 8, 16, 32, 64, 128, 256], 'draper_half': [1, 2, 4, 8, 16, 32, 64, 128, 256],
         'rgqft': [1, 2, 3, 4, 6, 8, 12, 16, 24, 32], 'hrs': [1, 2, 3, 4, 5, 6, 7, 8],
         'cdraper': [1, 2, 4, 6, 8, 9, 10],
         'pq_addconst': [2, 4, 8, 16, 32, 48, 64], 'pq_qft': [2, 4, 8, 16, 32, 48, 64],
         'pq_addconst_basic': [2, 4, 8, 16, 24, 28, 31, 32, 48, 64],
         'pq_qft_basic': [2, 4, 8, 16, 24, 28, 31, 32, 48, 64],
         'pq_addconst_opt': [2, 4, 8, 16, 24, 32]}
PQ_DIR = os.path.join(os.path.dirname(HERE), 'results', 'dyadic', 'projectq')


def peak_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(r / (1 << 20) if sys.platform == 'darwin' else r / 1024, 1)


# ----------------------------------------------------------------------------- circuits and their specifications
def qiskit_case(family, n):
    """-> (qiskit circuit or None, DCircuit, kind, spec or reference, extra)."""
    from qiskit import QuantumCircuit, transpile
    from qiskit.circuit.library import QFTGate
    from qiskit.synthesis.qft import synth_qft_full
    from qiskit.synthesis.arithmetic import adder_qft_d00, multiplier_qft_r17, multiplier_cumulative_h18
    from qoptct.dyadic.qiskit_io import from_qiskit
    if family in ('qft', 'iqft', 'qft_qftdg'):
        qc = QuantumCircuit(n)
        g = QFTGate(n)
        if family == 'qft':
            qc.append(g, range(n))
        elif family == 'iqft':
            qc.append(g.inverse(), range(n))
        else:
            qc.append(g, range(n))
            qc.append(g.inverse(), range(n))
    elif family == 'qft_noswap':
        qc = synth_qft_full(n, do_swaps=False)
    elif family in ('draper_fixed', 'draper_half'):
        qc = adder_qft_d00(n, kind=family.split('_')[1])
    elif family == 'cdraper':
        qc = adder_qft_d00(n, kind='fixed').control(1, annotated=False)     # control = qubit 0
    elif family == 'rgqft':
        qc = multiplier_qft_r17(n)
    elif family == 'hrs':
        qc = multiplier_cumulative_h18(n)
    else:
        raise ValueError(family)
    if family == 'hrs':
        # qubits_initially_zero=False: by default transpile may use idle qubits as clean ancillas, which keeps only
        # the action on |0...0> - not the unitary (found here: 4 of 16 inputs with out = 0 differ at n = 2)
        c = from_qiskit(transpile(qc, basis_gates=BASIS, optimization_level=0, qubits_initially_zero=False))
    else:
        c = from_qiskit(qc)
    if family == 'qft':
        spec = QFTSpec(range(n))
    elif family == 'iqft':
        spec = QFTSpec(range(n), inverse=True)
    elif family == 'qft_noswap':
        spec = QFTSpec(range(n), swaps=False)
    elif family == 'qft_qftdg':
        spec = DCircuit(n)
    elif family == 'draper_fixed':
        spec = adder_spec(range(n), range(n, 2 * n))
    elif family == 'draper_half':
        spec = adder_spec(range(n), range(n, 2 * n), 2 * n)
    elif family == 'cdraper':
        from qoptct.dyadic import ArithSpec
        spec = ArithSpec(range(1 + n, 1 + 2 * n), [(1 << j, (1 + j, 0)) for j in range(n)], name='controlled add')
    elif family == 'hrs':
        # HRS adds the partial products with half adders whose carry-out qubit must start in 0: the claim is
        # |a, b, 0> -> |a, b, a*b> (out register zero), checked against the classical map (no lemma)
        mask = (1 << n) - 1
        spec = ClassicalSpec(4 * n, lambda v: v | (((v & mask) * ((v >> n) & mask)) << (2 * n)),
                             zero_wires=range(2 * n, 4 * n), name='hrs')
    else:
        spec = multiplier_spec(range(n), range(n, 2 * n), range(2 * n, 4 * n))
    return qc, c, spec


def projectq_case(family, n):
    """ProjectQ circuits from the JSON dumps of scripts/dyadic_projectq_dump.py."""
    path = os.path.join(PQ_DIR, '%s_%d.json' % (family, n))
    d = json.load(open(path))
    # ProjectQ rounds rotation angles to 12 decimals (ANGLE_PRECISION) after every merge and decomposition step,
    # so errors add up (7e-13 seen after one halving); tolerance 2e-11 rad
    c = DCircuit(d["n"], float_abs_tol=2e-11)
    for nm, qs, angles in d['ops']:
        c.add(nm, qs, *angles)                    # float radians, as ProjectQ has them
    if d.get('global_phase'):
        c.add_global_phase(d['global_phase'])
    if family.startswith('pq_addconst'):
        spec = const_adder_spec(d['constant'], d['b'])
    else:
        spec = QFTSpec(d['wires'], swaps=d.get('swaps', False))
    return None, c, spec, d


def prove(c, spec, budget, phase='global'):
    if isinstance(spec, DCircuit):
        return prove_equal(c, spec, phase=phase, budget_sec=budget)
    return prove_spec(c, spec, phase=phase, budget_sec=budget)


def translation_check(qc, c):
    """Qiskit's own Operator of the original circuit vs the dense unitary of the translated DCircuit (<= 10 qubits)."""
    if qc is None or c.n > 10:
        return None
    import numpy as np
    from qiskit.quantum_info import Operator
    from qoptct.dyadic import sim
    U = Operator(qc).data
    V = sim.unitary(c)
    return bool(np.max(np.abs(U - V)) < 1e-8)


# ----------------------------------------------------------------------------- canaries
def canaries(c, spec, family, n, rng):
    """>= 20 wrong variants; about half at targeted places (smallest rotation, block junctions, global phase, rare
    inputs).  -> list of (name, circuit, phase mode)."""
    out = []
    # rotations with a nonzero angle only: a rotation by a multiple of 2*pi is the identity (RGQFT has many), so
    # changing it would not be a canary
    rot = [i for i, o in enumerate(c.ops) if o[2] and o[2][0] % 2 and o[0] in ('p', 'cp', 'mcp', 'rz', 'crz', 'u')]
    if not rot:
        rot = [i for i, o in enumerate(c.ops) if o[2] and any(r % 2 for r in o[2])]
    def ang(i):
        r = c.ops[i][2][0] % 2
        return min(r, 2 - r) if r else Fraction(2)
    nz = [i for i in rot if ang(i) != 2]
    small = min(nz, key=ang) if nz else None
    if small is not None:
        k = ang(small).denominator.bit_length() - 1
        for d in (Fraction(1, 2 ** k), Fraction(1, 2 ** (k + 1)), Fraction(-1, 2 ** (k + 1))):
            e = c.copy()
            nm, qs, rs = e.ops[small]
            e.ops[small] = (nm, qs, (rs[0] + d,) + rs[1:])
            out.append(('smallest_rotation%+s' % d, e, 'global'))
        e = c.copy()
        del e.ops[small]
        out.append(('lost_smallest_rotation', e, 'global'))
    for i in rng.sample(rot, min(6, len(rot))):
        e = c.copy()
        nm, qs, rs = e.ops[i]
        e.ops[i] = (nm, qs, (rs[0] + Fraction(1, 2 ** rng.randint(0, max(1, n))),) + rs[1:])
        out.append(('shifted_rotation', e, 'global'))
    for i in rng.sample(rot, min(3, len(rot))):
        e = c.copy()
        del e.ops[i]
        out.append(('lost_rotation', e, 'global'))
    two = [i for i in rot if len(c.ops[i][1]) >= 2]
    for i in rng.sample(two, min(3, len(two))):
        e = c.copy()
        nm, qs, rs = e.ops[i]
        free = [q for q in range(c.n) if q not in qs]
        if free:
            qs = list(qs)
            qs[rng.randrange(len(qs))] = rng.choice(free)
            e.ops[i] = (nm, tuple(qs), rs)
            out.append(('wrong_qubit', e, 'global'))
    hs = [i for i, o in enumerate(c.ops) if o[0] == 'h']
    if hs:
        e = c.copy()
        del e.ops[rng.choice(hs)]
        out.append(('lost_h', e, 'global'))
    sw = [i for i, o in enumerate(c.ops) if o[0] == 'swap']
    if sw:
        e = c.copy()
        del e.ops[rng.choice(sw)]
        out.append(('lost_swap', e, 'global'))
    # junction: a gate moved across the boundary of two blocks (QFT | rotations | inverse QFT)
    first_rot_after_h = None
    for i in range(1, len(c.ops)):
        if c.ops[i][0] in ('cp', 'p', 'mcp') and c.ops[i - 1][0] == 'h' and set(c.ops[i][1]) & set(c.ops[i - 1][1]):
            first_rot_after_h = i
            break
    if first_rot_after_h is not None:
        e = c.copy()
        e.ops[first_rot_after_h - 1], e.ops[first_rot_after_h] = e.ops[first_rot_after_h], e.ops[first_rot_after_h - 1]
        out.append(('junction_swap_order', e, 'global'))
    e = c.copy()
    e.gphase = (e.gphase + Fraction(1, 2 ** min(n + 1, 60))) % 2
    out.append(('global_phase_exact', e, 'exact'))
    # rare inputs: a phase or a flip that acts only when many inputs are 1
    zero_w = set(getattr(spec, 'zero_wires', ()))
    qs_all = [q for q in range(c.n) if q not in zero_w][:8]      # reachable: free inputs only
    if len(qs_all) >= 2:
        e = c.copy()
        e.ops.append(('mcp', tuple(qs_all), (Fraction(1, 2 ** min(n, 30)),)))
        out.append(('rare_input_phase', e, 'global'))
        e = c.copy()
        zero = list(getattr(spec, 'zero_wires', ()))
        if zero:
            # domain with zero wires: flip a zero wire when all free inputs are 1 (an output state that the circuit
            # really reaches; a flip conditioned on the whole register could be unreachable and then no canary)
            free = [q for q in range(c.n) if q not in zero]
            e.ops.append(('mcx', tuple(free) + (zero[0],), ()))
        else:
            e.ops.append(('mcx', tuple(range(c.n)), ()))
        out.append(('rare_input_flip_all_ones', e, 'global'))
    guard = 0
    while len(out) < 20 and rot and guard < 200:
        guard += 1
        i = rng.choice(rot)
        e = c.copy()
        nm, qs, rs = e.ops[i]
        e.ops[i] = (nm, qs, (-rs[0],) + rs[1:])
        if rs[0] % 2 not in (0, 1):         # P(-t) != P(t), and not only a global phase for RZ
            out.append(('negated_rotation', e, 'global'))
    return out


# ----------------------------------------------------------------------------- one run
def run_one(family, n, do_canaries, budget, seed=1, canary_budget=120.0):
    t0 = time.time()
    if family.startswith('pq_'):
        from qoptct.dyadic import NotDyadic, Unsupported
        try:
            qc, c, spec, meta = projectq_case(family, n)
        except (NotDyadic, Unsupported) as e:      # the gate list cannot be read exactly: not proved
            return dict(family=family, n=n, verdict='UNKNOWN', method='input', reason=str(e)[:300],
                        time_build=round(time.time() - t0, 2), peak_mb=peak_mb())
    else:
        qc, c, spec = qiskit_case(family, n)
        meta = {}
    t_build = time.time() - t0
    tr = translation_check(qc, c)
    t1 = time.time()
    r = prove(c, spec, budget)
    t_prove = time.time() - t1
    res = dict(family=family, n=n, qubits=c.n, gates=len(c.ops), verdict=r['verdict'], method=r.get('method'),
               strategy=r.get('strategy'), K=r.get('K'), global_phase=r.get('global_phase'),
               assumptions=r.get('assumptions'), dense_check=r.get('dense_check'), translation_check=tr,
               float_angles=r.get('float_angles'), reason=r.get('reason'), peak_terms=(r.get('stats') or {}).get(
                   'peak_terms'), time_build=round(t_build, 2), time_prove=round(t_prove, 2), peak_mb=peak_mb())
    if meta:
        res['source'] = {k: meta[k] for k in ('generator', 'projectq_version', 'constant', 'engine') if k in meta}
    if do_canaries:
        rng = random.Random(seed + n)
        cres = []
        for name, e, ph in canaries(c, spec, family, n, rng):
            t2 = time.time()
            rr = prove(e, spec, min(budget, canary_budget), phase=ph)
            cres.append(dict(name=name, verdict=rr['verdict'], method=rr.get('method'),
                             time=round(time.time() - t2, 2)))
        res['canaries'] = cres
        res['canaries_refuted'] = sum(1 for x in cres if x['verdict'] == 'REFUTED')
        res['canaries_proved'] = sum(1 for x in cres if x['verdict'] == 'PROVED')
        # independent numeric check at the largest size (product-state simulator, random basis inputs; not a proof)
        if not isinstance(spec, DCircuit):
            from qoptct.dyadic import prodsim
            chk = prodsim.check_qft if spec.kind == 'qft' else prodsim.check_arith
            t, bad, ent = chk(c, spec, 20, random.Random(seed))
            res['product_check'] = dict(inputs=t, bad=bad, entangling=ent)
    return res


def table(out):
    """Markdown summary of <out>/<family>.json: largest proved size, its time / memory / canaries / checks."""
    lines = ['| family | proved up to n | qubits | gates | prove s | peak MB | canaries REFUTED | product check | '
             'first size not proved |', '|' + '---|' * 9]
    for fam in SIZES:
        p = os.path.join(out, fam + '.json')
        if not os.path.exists(p):
            continue
        rows = json.load(open(p))['rows']
        ok = [r for r in rows if r.get('verdict') == 'PROVED']
        bad = [r for r in rows if r.get('verdict') != 'PROVED']
        if not ok:
            lines.append('| %s | - | | | | | | | %s %s |' % (fam, bad[0].get('n'), bad[0].get('verdict')))
            continue
        m = max(ok, key=lambda r: r['n'])
        cans = m.get('canaries') or []
        pc = m.get('product_check')
        lines.append('| %s | %d | %s | %s | %s | %s | %s | %s | %s |' % (
            fam, m['n'], m.get('qubits'), m.get('gates'), m.get('time_prove'), m.get('peak_mb'),
            '%d / %d' % (m.get('canaries_refuted', 0), len(cans)) if cans else '-',
            '%d inputs, %d bad' % (pc['inputs'] - pc['entangling'], pc['bad']) if pc else '-',
            ('%s %s' % (bad[0].get('n'), bad[0].get('verdict'))) if bad else '-'))
    return '\n'.join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    p1 = sub.add_parser('run')
    p1.add_argument('family')
    p1.add_argument('n', type=int)
    p1.add_argument('--canaries', action='store_true')
    p1.add_argument('--budget', type=float, default=600.0)
    p1.add_argument('--canary-budget', type=float, default=120.0)
    p3 = sub.add_parser('table')
    p3.add_argument('--out', default=os.path.join(os.path.dirname(HERE), 'results', 'dyadic'))
    p2 = sub.add_parser('registry')
    p2.add_argument('--families', default=','.join(f for f in SIZES if not f.startswith('pq_')))
    p2.add_argument('--sizes', default=None)
    p2.add_argument('--budget', type=float, default=600.0)
    p2.add_argument('--out', default=os.path.join(os.path.dirname(HERE), 'results', 'dyadic'))
    a = ap.parse_args(argv)
    if a.cmd == 'run':
        print(json.dumps(run_one(a.family, a.n, a.canaries, a.budget, canary_budget=a.canary_budget)), flush=True)
        return 0
    if a.cmd == 'table':
        print(table(a.out))
        return 0
    os.makedirs(a.out, exist_ok=True)
    for fam in a.families.split(','):
        sizes = [int(s) for s in a.sizes.split(',')] if a.sizes else SIZES[fam]
        rows = []

        def one(n, can):
            cmd = [sys.executable, os.path.abspath(__file__), 'run', fam, str(n), '--budget', str(a.budget)]
            if can:
                cmd.append('--canaries')
            try:
                p = subprocess.run(cmd, capture_output=True, text=True, timeout=a.budget * 40 + 120)
                line = [x for x in p.stdout.splitlines() if x.startswith('{')]
                return json.loads(line[-1]) if line else dict(family=fam, n=n, verdict='ERROR',
                                                              reason=(p.stderr or '')[-400:])
            except subprocess.TimeoutExpired:
                return dict(family=fam, n=n, verdict='UNKNOWN', reason='driver timeout')
        for n in sizes:
            row = one(n, False)
            rows.append(row)
            print(fam, n, row.get('verdict'), row.get('time_prove'), row.get('peak_mb'), flush=True)
            if row.get('verdict') != 'PROVED' and not fam.startswith('pq_'):
                break
        ok = [i for i, r in enumerate(rows) if r.get('verdict') == 'PROVED']
        if ok:
            i = ok[-1]
            rows[i] = one(rows[i]['n'], True)
            print(fam, rows[i]['n'], 'canaries', rows[i].get('canaries_refuted'), '/',
                  len(rows[i].get('canaries', [])), 'proved', rows[i].get('canaries_proved'), flush=True)
        with open(os.path.join(a.out, fam + '.json'), 'w') as f:
            json.dump(dict(family=fam, rows=rows, python=sys.version.split()[0]), f, indent=1)
    return 0


if __name__ == '__main__':
    sys.exit(main())
