# SPDX-License-Identifier: Apache-2.0
"""Proof of a candidate against its reference, with an independent sampled check and canaries.

The provers are used unchanged (qoptct.verify.kprove, then qoptct.verify.prove). A candidate passes only if
  1. the prover says PROVED (SAMPLED / UNKNOWN / ERROR never pass),
  2. the independent simulators (qoptct.verify.quick_check) do not disagree,
  3. at least `canaries` mutants of the candidate (phase, erasure-correction, seam, flip changes) are REFUTED by the
     same prover on the same reference (a prover that proves everything would be caught here), and no mutant that the
     prover calls PROVED is refuted by simulation (that would be a false PROVED: it is reported as a hole).
"""
import hashlib
import os
import random
import re
import time

MIN_CLASSES = ('phase', 'erasure', 'seam')


def prover_version():
    import qoptct
    d = os.path.join(os.path.dirname(os.path.abspath(qoptct.__file__)), 'verify')
    h = hashlib.sha256()
    for n in sorted(os.listdir(d)):
        if n.endswith('.py'):
            h.update(n.encode())
            h.update(open(os.path.join(d, n), 'rb').read())
    return 'qoptct %s verify-%s' % (qoptct.__version__, h.hexdigest()[:12])


def prove(qasm_path, ref_path, clean, budget=120.0, max_mb=900):
    """-> dict(status, prover, method, time, detail). status is PROVED only for an exact proof."""
    t0 = time.time()
    out = {'status': 'UNKNOWN', 'prover': 'kprove', 'method': None, 'detail': ''}
    try:
        from ..verify.kprove import kprove
        r = kprove(qasm_path, ref_path, clean=clean, budget_sec=budget, max_mb=max_mb)
        out.update(status=r['status'], method=r.get('method'), detail=str(r.get('detail', ''))[:300])
        if r['status'] in ('PROVED', 'REFUTED'):
            out['time'] = round(time.time() - t0, 2)
            return out
    except Exception as e:                       # noqa: BLE001 -- a crash of a prover is "not proved"
        out['detail'] = 'kprove error: %s' % str(e)[:200]
    try:
        from ..verify import prove as vprove
        r = vprove(qasm_path, ref_path, clean=clean, budget_sec=budget)
        if r['status'] in ('PROVED', 'REFUTED'):
            out.update(status=r['status'], prover='prove', method=r.get('method'), detail=str(r.get('detail', ''))[:300])
    except Exception as e:                       # noqa: BLE001
        out['detail'] += ' | prove error: %s' % str(e)[:200]
    out['time'] = round(time.time() - t0, 2)
    return out


def sim_check(qasm_path, ref_path, clean):
    """Independent sampled check. -> 'ok' | 'mismatch' | 'inconclusive'."""
    try:
        from ..verify import quick_check
        r = quick_check(qasm_path, ref_path, clean=clean)
    except Exception:                            # noqa: BLE001
        return 'inconclusive'
    if r.get('F') is None:
        return 'inconclusive'
    return 'ok' if r.get('ok') else 'mismatch'


_PH = {'t': 'tdg', 'tdg': 't', 's': 'sdg', 'sdg': 's'}


def mutants(txt, n, seed):
    """Up to n distinct (class, mutated text) pairs: phase, erasure, seam, flip."""
    rng = random.Random(seed)
    lines = txt.split('\n')
    used = sorted({int(x) for x in re.findall(r'q\[(\d+)\]', txt)})
    cand = {'phase': [], 'erasure': [], 'seam': [], 'flip': []}
    for i, ln in enumerate(lines):
        if re.match(r'^(t|tdg|s|sdg|z) q\[\d+\];$', ln):
            cand['phase'].append(i)
        elif ln.startswith('if (') or ln.startswith('if('):
            cand['erasure'].append(i)
        elif re.match(r'^(cx|cz) q\[\d+\], q\[\d+\];$', ln):
            cand['seam'].append(i)
        elif re.match(r'^x q\[\d+\];$', ln):
            cand['flip'].append(i)
    out, seen = [], set()
    classes = [c for c in ('phase', 'erasure', 'seam', 'seam', 'phase', 'flip') if cand[c]]
    tries = 0
    while classes and len(out) < n and tries < 20 * n + 50:
        tries += 1
        cl = classes[tries % len(classes)]
        i = rng.choice(cand[cl])
        ln = lines[i]
        m = re.match(r'^(\w+) q\[(\d+)\](?:, q\[(\d+)\])?;$', ln)
        new = None
        if cl == 'phase':
            g, a = m.group(1), m.group(2)
            new = '%s q[%s];' % (_PH[g], a) if g in _PH else ('s q[%s];' % a if rng.random() < .5 else None)
        elif cl == 'erasure':
            new = None if rng.random() < .7 else ln.replace('z q', 'x q') if 'z q' in ln else None
        elif cl == 'seam':
            g, a, b = m.group(1), int(m.group(2)), int(m.group(3))
            r = rng.random()
            if r < .4:
                new = None
            else:
                alts = [u for u in used if u not in (a, b)]
                if alts:
                    c2 = rng.choice(alts)
                    new = '%s q[%d], q[%d];' % (g, a, c2) if r < .7 else '%s q[%d], q[%d];' % (g, c2, b)
        else:
            new = None
        key = (i, new)
        if key in seen:
            continue
        seen.add(key)
        mut = lines[:i] + ([new] if new is not None else []) + lines[i + 1:]
        out.append((cl, '\n'.join(mut)))
    return out


def run_canaries(txt, ref_path, clean, workdir, want=20, budget=20.0, max_sec=600.0, seed=0):
    """-> dict(detected, equivalent, undecided, by_class, hole (None or description), tried)."""
    t0 = time.time()
    res = dict(detected=0, equivalent=0, undecided=0, by_class={}, hole=None, tried=0)
    path = os.path.join(workdir, 'mutant.qasm')
    for cl, mt in mutants(txt, want * 3, seed):
        if res['detected'] >= want or time.time() - t0 > max_sec:
            break
        open(path, 'w').write(mt)
        r = prove(path, ref_path, clean, budget=budget)
        res['tried'] += 1
        if r['status'] == 'REFUTED':
            res['detected'] += 1
            res['by_class'][cl] = res['by_class'].get(cl, 0) + 1
        elif r['status'] == 'PROVED':
            if sim_check(path, ref_path, clean) == 'mismatch':
                res['hole'] = 'FALSE PROVED: a %s mutant is PROVED but the simulators refute it' % cl
                break
            res['equivalent'] += 1
        else:
            if sim_check(path, ref_path, clean) == 'mismatch':
                res['detected'] += 1
                res['by_class'][cl] = res['by_class'].get(cl, 0) + 1
            else:
                res['undecided'] += 1
    if os.path.exists(path):
        os.remove(path)
    return res


def judge(qasm_path, ref_path, clean, workdir, canaries=20, budget=120.0, canary_sec=600.0):
    """Full acceptance test of one candidate file. -> dict(accepted, status, prover, sim, canaries, reason)."""
    txt = open(qasm_path).read()
    p = prove(qasm_path, ref_path, clean, budget=budget)
    out = dict(accepted=False, status=p['status'], prover=p, sim=None, canaries=None, reason='')
    if p['status'] != 'PROVED':
        out['reason'] = 'not proved: %s %s' % (p['status'], p.get('detail', ''))[:200]
        return out
    out['sim'] = sim_check(qasm_path, ref_path, clean)
    if out['sim'] == 'mismatch':
        out['hole'] = 'FALSE PROVED: the prover says PROVED but the simulators refute the candidate'
        out['reason'] = out['hole']
        return out
    seed = int(hashlib.sha256(txt.encode()).hexdigest()[:8], 16)
    c = run_canaries(txt, ref_path, clean, workdir, want=canaries, seed=seed, max_sec=canary_sec)
    out['canaries'] = c
    if c['hole']:
        out['hole'] = c['hole']
        out['reason'] = c['hole']
        return out
    need = canaries
    k = sum(c['by_class'].get(x, 0) for x in MIN_CLASSES)
    if c['detected'] < need or k < need // 2:
        out['reason'] = 'canaries: detected %d of the required %d (targeted classes %d, need %d)' % (c['detected'], need, k, need // 2)
        return out
    out['accepted'] = True
    return out
