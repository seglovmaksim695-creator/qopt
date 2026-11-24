# SPDX-License-Identifier: Apache-2.0
"""`qoptct loop verify`: re-check the whole journal. Nothing is believed because it is written in a file.

For every record: hash chain and replay rules, the stored files against their SHA-256, the reference against the hash
pinned by the baseline and against the target definition, the costs recomputed from the stored circuit, and (unless
--no-reprove) the proof itself, run again on the stored files together with the independent simulators. The table
files are regenerated from the journal and compared with the ones on disk.
"""
import json
import os
import tempfile

from .costs import qasm_cost
from .journal import State, flags_hash, replay, sha256_file, table_from
from .prover import prove, sim_check


def verify(state_dir, reprove=True, budget=120.0, log=None, allow_empty=False):
    """-> (ok, problems, n_records). A missing or empty journal is a problem unless allow_empty (an empty check that
    says OK is a false success); the number of records that went through the full check must equal the number of
    records in the journal."""
    st = State(state_dir)
    if not os.path.isdir(st.root):
        return False, ['state directory not found: %s' % st.root], 0
    if not os.path.exists(st.journal):
        if not allow_empty:
            return False, ['journal not found: %s' % st.journal], 0
    recs, probs = st.read()
    if not recs and not allow_empty and not probs:
        probs.append('journal is empty: nothing to verify (use --allow-empty if that is intended)')
    checked = 0
    _, _, rp = replay(recs)
    probs = list(probs) + rp
    tdefs = {}
    for i, r in enumerate(recs):
        tag = 'record %d' % (i + 1)
        try:
            t, csha, rsha = r['target'], r['cand_sha'], r['ref_sha']
        except KeyError as e:
            probs.append('%s: missing field %s' % (tag, e))
            continue
        cp, rpth = st.path(csha, '.qasm'), st.path(rsha, '.qc')
        if not os.path.exists(cp):
            probs.append('%s: candidate file is missing from the store' % tag)
            continue
        if not os.path.exists(rpth):
            probs.append('%s: reference file is missing from the store' % tag)
            continue
        bad = False
        if sha256_file(cp) != csha:
            probs.append('%s: candidate file does not match its SHA-256 (replaced after the proof)' % tag)
            bad = True
        if sha256_file(rpth) != rsha:
            probs.append('%s: reference file does not match its SHA-256' % tag)
            bad = True
        if t not in tdefs:
            try:
                tdefs[t] = st.load_target(t)
            except Exception:                      # noqa: BLE001
                tdefs[t] = None
        td = tdefs[t]
        if td is None:
            probs.append('%s: no target definition for %s' % (tag, t))
            continue
        if td['ref_sha'] != rsha:
            probs.append('%s: reference hash differs from the target definition of %s' % (tag, t))
        if 'flags_sha' in r and r['flags_sha'] != flags_hash(td):
            probs.append('%s: the flags of target %s (clean / spec / generator) differ from those the proofs were made with'
                         % (tag, t))
        if not bad:
            try:
                real = qasm_cost(open(cp).read())
            except Exception as e:                 # noqa: BLE001
                real = None
                probs.append('%s: cannot read candidate (%s)' % (tag, e))
            if real is not None and real != r.get('cost_after'):
                probs.append('%s: cost_after %s differs from the cost of the stored circuit %s' % (tag, r.get('cost_after'), real))
        if reprove and not bad:
            with tempfile.TemporaryDirectory() as _:
                p = prove(cp, rpth, bool(td['clean']), budget=budget)
            if p['status'] != 'PROVED':
                probs.append('%s: the proof does not hold on re-check: %s %s' % (tag, p['status'], p.get('detail', '')[:120]))
            elif sim_check(cp, rpth, bool(td['clean'])) == 'mismatch':
                probs.append('%s: the independent simulation refutes the candidate' % tag)
        checked += 1
        if log:
            log('checked record %d of %d' % (i + 1, len(recs)))
    if checked != len(recs) and not probs:
        probs.append('only %d of %d records were checked' % (checked, len(recs)))
    # the table is a view of the journal
    try:
        want = table_from(recs)
    except Exception:                              # noqa: BLE001
        want = None
        probs.append('the table cannot be derived from the journal (malformed record)')
    bp = os.path.join(st.root, 'best.json')
    if os.path.exists(bp) and want is not None:
        try:
            have = json.load(open(bp))
        except ValueError:
            have = None
        try:
            if have != json.loads(json.dumps(want)):
                probs.append('best.json differs from the table derived from the journal (edited by hand?)')
        except Exception as e:                     # noqa: BLE001
            probs.append('best.json: ' + str(e))
    return (not probs), probs, len(recs)
