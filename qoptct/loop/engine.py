# SPDX-License-Identifier: Apache-2.0
"""The loop: propose -> prove -> accumulate. Resumable (attempts.jsonl + journal.jsonl), time-limited, one child process
per step (killed on a timeout)."""
import json
import os
import shutil
import subprocess
import sys
import time
import uuid

from . import proposers
from .costs import best_of, dominates, front_insert
from .journal import State, flags_hash, replay, sha256_file, utc_now, write_table
from .targets import register


def _fronts(state):
    recs, probs = state.read()
    fronts, pinned, rp = replay(recs)
    if probs or rp:
        raise RuntimeError('journal is not consistent, run `qoptct loop verify`: ' + '; '.join((probs + rp)[:3]))
    return fronts, pinned


def _run_child(state, job, timeout):
    jid = uuid.uuid4().hex[:10]
    wd = os.path.join(state.work, jid)
    os.makedirs(wd, exist_ok=True)
    job = dict(job, state=state.root, workdir=wd)
    jf, rf = os.path.join(wd, 'job.json'), os.path.join(wd, 'result.json')
    json.dump(job, open(jf, 'w'))
    try:
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(
            [os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))] + [os.environ.get('PYTHONPATH', '')]))
        subprocess.run([sys.executable, '-m', 'qoptct.loop.worker', jf, rf], env=env, stdin=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=timeout)
    except subprocess.TimeoutExpired:
        shutil.rmtree(wd, ignore_errors=True)
        return {'outcome': 'timeout'}, None
    try:
        res = json.load(open(rf))
    except Exception:                              # noqa: BLE001
        res = {'outcome': 'error', 'detail': 'child produced no result (killed or crashed)'}
    return res, wd


def _accept(state, tid, tdef, res, wd, role, params, before):
    cand_file = os.path.join(wd, 'cand.qasm')
    if sha256_file(cand_file) != res['cand_sha']:
        return None
    data = open(cand_file, 'rb').read()
    sha = state.put(data, '.qasm')
    assert sha == res['cand_sha']
    j = res['judge']
    cl = j['canaries']
    clean = ' --clean' if tdef['clean'] else ''
    kind = 'baseline' if before is None else ('dominates' if dominates(res['cost'], before) else 'tradeoff')
    rec = dict(role=role, kind=kind, target=tid, ref_sha=tdef['ref_sha'], flags_sha=flags_hash(tdef), cand_sha=sha, cost_before=before, cost_after=res['cost'],
               prover={'name': j['prover']['prover'], 'version': res['prover_version'], 'status': j['status'],
                       'method': j['prover'].get('method'), 'time': j['prover'].get('time')},
               sim=j['sim'], canaries={k: cl[k] for k in ('detected', 'equivalent', 'undecided', 'by_class', 'tried')},
               proposer=params,
               repro='qoptct kprove store/%s.qasm store/%s.qc%s' % (sha, tdef['ref_sha'], clean))
    return state.append(rec)


def _hole(state, tid, res, wd, params):
    d = os.path.join(state.root, 'HOLES')
    os.makedirs(d, exist_ok=True)
    name = '%s_%s' % (tid, utc_now().replace(':', ''))
    if wd and os.path.exists(os.path.join(wd, 'cand.qasm')):
        shutil.copy(os.path.join(wd, 'cand.qasm'), os.path.join(d, name + '.qasm'))
    json.dump({'target': tid, 'params': params, 'result': res}, open(os.path.join(d, name + '.json'), 'w'), indent=1,
              default=str)


def run(state_dir, specs, names=('lookup', 'optimize'), hours=1.0, canaries=20, budget=120.0, canary_sec=600.0,
        job_timeout=1800, prove_all=False, max_steps=None, log=print, repo_root=None):
    """Run the loop until the time limit, `max_steps` or until every proposer is exhausted. Returns the counters."""
    state = State(state_dir).create()
    lock = state.lock()
    t_end = time.time() + hours * 3600
    targets = []
    for s in specs:
        d, base = register(state, s, repo_root)
        targets.append((d, base))
    done = {a['key'] for a in state.read_attempts() if a.get('key') and a.get('outcome') not in ('error', 'timeout')}
    fronts, pinned = _fronts(state)
    # baselines (proved like everything else)
    for d, base in targets:
        tid = d['id']
        if tid in fronts or base is None:
            continue
        key = proposers.key_of(tid, {'p': 'baseline'})
        bf = os.path.join(state.work, 'baseline_%s.qasm' % tid)
        open(bf, 'w').write(base)
        res, wd = _run_child(state, dict(target=tid, qasm_file=bf, front=[], canaries=canaries, budget=budget,
                                         canary_sec=canary_sec, prove_all=True), job_timeout)
        out = res['outcome']
        if out == 'accepted':
            rec = _accept(state, tid, d, res, wd, 'baseline', {'p': 'baseline'}, None)
            log('baseline %s: accepted (AND %d, CX %d)' % (tid, res['cost']['AND'], res['cost']['CX']))
        else:
            log('baseline %s: %s %s' % (tid, out, res.get('judge', {}).get('reason', res.get('detail', ''))))
            if out == 'hole':
                _hole(state, tid, res, wd, {'p': 'baseline'})
        state.log_attempt(dict(key=key, target=tid, outcome=out, params={'p': 'baseline'}, secs=res.get('secs')))
        if wd:
            shutil.rmtree(wd, ignore_errors=True)
        fronts, pinned = _fronts(state)
    streams = {}
    for d, _ in targets:
        streams[d['id']] = (d, proposers.stream(dict(d, ref_text=open(state.path(d['ref_sha'], '.qc')).read()), names))
    steps = 0
    live = list(streams)
    while live and time.time() < t_end and (max_steps is None or steps < max_steps):
        nxt = []
        for tid in live:
            if time.time() >= t_end or (max_steps is not None and steps >= max_steps):
                nxt.append(tid)
                continue
            d, gen = streams[tid]
            params = None
            for p in gen:
                if proposers.key_of(tid, p) not in done:
                    params = p
                    break
            if params is None:
                continue
            nxt.append(tid)
            if tid not in fronts and not d.get('has_baseline') and params['p'] != 'optimize':
                pass                                # a lookup candidate may become the baseline of a qc target
            key = proposers.key_of(tid, params)
            front = fronts.get(tid, [])
            res, wd = _run_child(state, dict(target=tid, params=params, front=front, canaries=canaries, budget=budget,
                                             canary_sec=canary_sec, prove_all=prove_all), job_timeout)
            out = res['outcome']
            steps += 1
            if out == 'accepted':
                role = 'improvement' if tid in fronts else 'baseline'
                rec = _accept(state, tid, d, res, wd, role, params, best_of(front) if front else None)
                if rec is None:
                    out = 'error'
                else:
                    log('ACCEPTED %s %s: AND %d CX %d qubits %d' % (tid, key.split('|', 1)[1], res['cost']['AND'],
                                                                   res['cost']['CX'], res['cost']['qubits']))
                    fronts, pinned = _fronts(state)
            elif out == 'hole':
                _hole(state, tid, res, wd, params)
                log('HOLE %s %s' % (tid, res['judge'].get('reason')))
            att = dict(key=key, target=tid, outcome=out, params=params, secs=round(res.get('secs', 0), 2))
            if 'cost' in res:
                att['cost'] = res['cost']
            if out in ('rejected', 'unproved', 'refuted', 'error', 'build_error'):
                att['detail'] = (res.get('judge', {}).get('reason') or res.get('detail') or '')[:200]
            state.log_attempt(att)
            done.add(key)
            if wd:
                shutil.rmtree(wd, ignore_errors=True)
        live = nxt
    write_table(state)
    lock.close()
    return stats(state)


def stats(state):
    att = state.read_attempts()
    c = {}
    for a in att:
        c[a['outcome']] = c.get(a['outcome'], 0) + 1
    recs, _ = state.read()
    return {'attempts': len(att), 'outcomes': c, 'journal_records': len(recs),
            'improvements': sum(1 for r in recs if r['role'] == 'improvement')}
