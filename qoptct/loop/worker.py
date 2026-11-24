# SPDX-License-Identifier: Apache-2.0
"""One loop step in a child process: build a candidate, cost it, and (if it can enter the front) prove it.

    python -m qoptct.loop.worker JOB.json RESULT.json

The child is killed by the parent on a timeout; a crash or a kill is an attempt with outcome "error" and nothing is
accepted. The parent re-hashes the candidate file after the child has finished.
"""
import json
import os
import sys
import time
import traceback

from .costs import front_insert, qasm_cost
from .journal import State, sha256_bytes, sha256_file


def run(job):
    t0 = time.time()
    st = State(job['state'])
    tdef = st.load_target(job['target'])
    ref_path = st.path(tdef['ref_sha'], '.qc')
    if sha256_file(ref_path) != tdef['ref_sha']:
        return {'outcome': 'error', 'detail': 'reference copy does not match its hash'}
    target = dict(tdef, ref_text=open(ref_path).read())
    work = job['workdir']
    os.makedirs(work, exist_ok=True)
    if job.get('qasm_file'):
        txt = open(job['qasm_file']).read()
    else:
        from .proposers import build
        try:
            txt = build(target, job['params'], work)
        except Exception as e:                    # noqa: BLE001
            return {'outcome': 'build_error', 'detail': ('%s: %s' % (type(e).__name__, e))[:300], 'secs': time.time() - t0}
    if txt is None:
        return {'outcome': 'no_candidate', 'secs': time.time() - t0}
    cost = qasm_cost(txt)
    _, gain = front_insert(job['front'], cost)
    if not gain and not job.get('prove_all'):
        return {'outcome': 'no_gain', 'cost': cost, 'secs': time.time() - t0}
    cand = os.path.join(work, 'cand.qasm')
    with open(cand, 'w') as f:
        f.write(txt)
    sha = sha256_bytes(txt.encode())
    from .prover import judge, prover_version
    j = judge(cand, ref_path, bool(tdef['clean']), work, canaries=job['canaries'], budget=job['budget'],
              canary_sec=job['canary_sec'])
    if sha256_file(cand) != sha:
        return {'outcome': 'error', 'detail': 'candidate file changed during the proof'}
    out = {'cost': cost, 'cand_sha': sha, 'judge': j, 'prover_version': prover_version(), 'gain': gain,
           'secs': time.time() - t0}
    if j.get('hole'):
        out['outcome'] = 'hole'
    elif j['accepted'] and gain:
        out['outcome'] = 'accepted'
    elif j['accepted']:
        out['outcome'] = 'proved_no_gain'
    elif j['status'] == 'REFUTED':
        out['outcome'] = 'refuted'
    elif j['status'] != 'PROVED':
        out['outcome'] = 'unproved'
    else:
        out['outcome'] = 'rejected'
    return out


def main():
    job = json.load(open(sys.argv[1]))
    try:
        res = run(job)
    except Exception:                             # noqa: BLE001
        res = {'outcome': 'error', 'detail': traceback.format_exc()[-600:]}
    with open(sys.argv[2], 'w') as f:
        json.dump(res, f, default=str)


if __name__ == '__main__':
    main()
