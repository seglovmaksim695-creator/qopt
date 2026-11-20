#!/usr/bin/env python
# SPDX-License-Identifier: Apache-2.0
"""Word-level algebraic proofs (qoptct.verify.algebra) of the qoptct.crypto arithmetic at full size, canaries, and a
comparison with the BDD engine of kprove.  Every job runs in its own process (time and peak memory are its own).

    python scripts/algebra_proofs.py proofs   [--out results.json]     n = 256 multiply / square / fold / modmul
    python scripts/algebra_proofs.py canaries [--out canaries.json]    broken n = 256 circuits: never PROVED
    python scripts/algebra_proofs.py compare  [--out compare.json]     algebra vs kprove (BDD) on small multipliers
    python scripts/algebra_proofs.py job '<json>'                      one job (used internally)
"""
import json
import random
import subprocess
import sys
import time


def build(job):
    from qoptct.crypto import fold, mulsq
    from qoptct.verify import algebra as A
    k = job['kind']
    if k in ('mul', 'sq'):
        c = (mulsq.multiply if k == 'mul' else mulsq.square)(job['n'], karatsuba_levels=job.get('L', 0),
                                                              clean=job.get('clean', False),
                                                              legacy=job.get('legacy', False))
        N, ops = A.ops_of_mulsq(c)
        return N, ops, A.spec_mulsq(c)
    p = fold.CURVES[job['curve']]
    R = fold.reduce_special(p)
    if k == 'fold':
        N, ops = A.ops_of_andcirc(R.circ, R.order)
        return N, ops, A.spec_fold(R)
    m = (mulsq.multiply if job.get('op', 'mul') == 'mul' else mulsq.square)(R.n, karatsuba_levels=job.get('L', 0))
    return A.compose_modmul(m, R)


def mutate(N, ops, spec, mut):
    """one deliberate error; -> (ops, spec)."""
    kind, frac = mut['type'], mut['at']
    rng = random.Random(mut.get('seed', 1))

    def pick(pred):
        idx = [i for i, o in enumerate(ops) if pred(o)]
        return idx[min(len(idx) - 1, int(frac * len(idx)))]
    ops = list(ops)
    if kind == 'drop_cx':
        del ops[pick(lambda o: o[0] == 'cx')]
    elif kind == 'drop_x':
        del ops[pick(lambda o: o[0] == 'x')]
    elif kind == 'cx_control':
        i = pick(lambda o: o[0] == 'cx')
        c = (ops[i][1] + 1 + rng.randrange(N - 2)) % N
        ops[i] = ('cx', c if c != ops[i][2] else (c + 1) % N, ops[i][2])
    elif kind == 'and_operand':
        i = pick(lambda o: o[0] == 'and')
        a = next(w for w in range(ops[i][1] + 1, ops[i][1] + N) if w % N not in ops[i][1:]) % N
        ops[i] = ('and', a, ops[i][2], ops[i][3])
    elif kind == 'drop_z':
        i = pick(lambda o: o[0] == 'erase' and any(c[0] == 'z' for c in o[2]))
        ops[i] = ('erase', ops[i][1], tuple(c for c in ops[i][2] if c[0] != 'z'))
    elif kind == 'drop_cz':
        i = pick(lambda o: o[0] == 'erase' and any(c[0] == 'cz' for c in o[2]))
        ops[i] = ('erase', ops[i][1], tuple(c for c in ops[i][2] if c[0] != 'cz'))
    elif kind == 'drop_erase':
        del ops[pick(lambda o: o[0] == 'erase')]
    elif kind == 'claim':
        spec = dict(spec, claims=[c.replace('==', '== 1 +', 1) if '==' in c else c for c in spec['claims']])
    elif kind == 'rare':
        # an error on a 2^-40 fraction of the inputs (first input word: 20 low bits 1, the next 20 bits 0): random
        # simulation does not see it; the output bit 0 is flipped there.  New clean wires hold the AND chain.
        ins = list(spec['inputs'].values())[0]
        out = list(spec['outputs'].values())[0]
        ones, zeros = ins[:20], ins[20:40]
        ch = [('x', w) for w in zeros]
        lits = ones + zeros
        prev = lits[0]
        anc = []
        for j, w in enumerate(lits[1:]):
            t = N + j
            anc.append((prev, w, t))
            ch.append(('and', prev, w, t))
            prev = t
        ops += ch + [('cx', prev, out[0])]
        ops += [('erase', t, (('cz', a, b),)) for a, b, t in reversed(anc)] + [('x', w) for w in zeros]
        N += len(anc)
    else:
        raise ValueError(kind)
    return N, ops, spec


def run_job(job):
    from qoptct.verify import algebra as A
    t0 = time.time()
    N, ops, spec = build(job)
    if job.get('mut'):
        try:
            N, ops, spec = mutate(N, ops, spec, job['mut'])
        except IndexError:                           # the circuit has no gate of this kind (e.g. no Z corrections)
            return dict(status='N/A', detail='no such gate in this circuit', job=job, time=0.0)
    tb = time.time() - t0
    if job.get('engine') == 'kprove':
        import os
        import tempfile
        from qoptct.crypto import mulsq
        from qoptct.verify.kprove import kprove
        c = (mulsq.multiply if job['kind'] == 'mul' else mulsq.square)(job['n'])
        d = tempfile.mkdtemp()
        qa, qc = os.path.join(d, 'c.qasm'), os.path.join(d, 'r.qc')
        open(qa, 'w').write(c.to_qasm(num_qubits=c.reference_width()))
        open(qc, 'w').write(c.reference_qc())
        t1 = time.time()
        r = kprove(qa, qc, clean=True, budget_sec=job.get('budget', 300), engine='bdd', sample=False,
                   max_mb=job.get('max_mb', 1500))
        r = dict(status=r['status'], detail=str(r.get('detail'))[:160], time=round(time.time() - t1, 2),
                 counts=dict(AND=c.counts()['AND'], qubits=c.num_qubits))
    else:
        text = A.ops_to_qasm(N, ops)                 # through the text: the strict parser is part of the proof
        t1 = time.time()
        r = A.aprove(text, spec, budget_sec=job.get('budget', 1800))
        r['time'] = round(time.time() - t1, 2)
        r['detail'] = r.get('detail', '')[:200]
    r['build_time'] = round(tb, 2)
    r['peak_mb'] = A._peak_mb()
    r['job'] = job
    return r


def child(job, timeout=3600):
    try:
        p = subprocess.run([sys.executable, __file__, 'job', json.dumps(job)], capture_output=True, text=True,
                           timeout=timeout)
        return json.loads(p.stdout.strip().splitlines()[-1])
    except Exception as ex:                          # noqa
        return dict(status='ERROR', detail='job failed: %r' % ex, job=job)


def name(job):
    s = job['kind'] + (' ' + job['curve'] if 'curve' in job else '') + (' n=%d' % job['n'] if 'n' in job else '')
    if job['kind'] != 'fold':
        s += ' L=%d' % job.get('L', 0)
    for f in ('legacy', 'clean'):
        if job.get(f):
            s += ' ' + f
    if job.get('op'):
        s += ' ' + job['op']
    if job.get('engine'):
        s += ' [kprove bdd]'
    if job.get('mut'):
        s += ' !%s@%.2f' % (job['mut']['type'], job['mut']['at'])
    return s


def main(argv):
    mode = argv[0]
    if mode == 'job':
        print(json.dumps(run_job(json.loads(argv[1])), default=str))
        return 0
    out = argv[argv.index('--out') + 1] if '--out' in argv else None
    jobs = []
    if mode == 'proofs':
        for k in ('mul', 'sq'):
            jobs += [dict(kind=k, n=256, L=L) for L in (0, 1, 2, 3)]
            jobs += [dict(kind=k, n=256, L=L, legacy=True) for L in (0, 3)]
        jobs += [dict(kind='fold', curve=c) for c in ('curve25519', 'secp256k1', 'p256', 'p384', 'p521')]
        for c in ('curve25519', 'secp256k1', 'p256'):
            jobs += [dict(kind='modmul', curve=c, L=L, op='mul') for L in (0, 3)]
            jobs += [dict(kind='modmul', curve=c, L=L, op='sq') for L in (0, 2)]
    elif mode == 'canaries':
        base = [dict(kind='mul', n=256, L=0), dict(kind='mul', n=256, L=3), dict(kind='sq', n=256, L=3),
                dict(kind='fold', curve='secp256k1'), dict(kind='modmul', curve='curve25519', L=3, op='mul')]
        for b in base:
            for t in ('drop_cx', 'drop_x', 'cx_control', 'and_operand', 'drop_z', 'drop_cz', 'drop_erase'):
                for at in (0.03, 0.5, 0.97):
                    jobs.append(dict(b, mut=dict(type=t, at=at)))
            jobs.append(dict(b, mut=dict(type='claim', at=0)))
            jobs.append(dict(b, mut=dict(type='rare', at=0)))
    elif mode == 'compare':
        for n in (4, 6, 8, 10, 12):
            jobs.append(dict(kind='mul', n=n))
            jobs.append(dict(kind='mul', n=n, engine='kprove', budget=300))
    else:
        print(__doc__)
        return 3
    res = []
    for job in jobs:
        r = child(job)
        res.append(r)
        print('%-52s %-8s %7.1f s %7s MB  AND=%s' % (name(job), r['status'], r.get('time', 0.0), r.get('peak_mb'),
                                                  (r.get('counts') or {}).get('AND')), flush=True)
        if out:
            json.dump(res, open(out, 'w'), indent=1, default=str)
    if mode == 'canaries':
        res = [r for r in res if r['status'] != 'N/A']
        bad = [r for r in res if r['status'] == 'PROVED']
        cnt = {}
        for r in res:
            cnt[r['status']] = cnt.get(r['status'], 0) + 1
        print('canaries: %d, by status: %s, PROVED (must be 0): %d' % (len(res), cnt, len(bad)))
        return 1 if bad else 0
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
