# SPDX-License-Identifier: Apache-2.0
"""The journal of proved improvements (JSONL, one record per line, hash chain) and the state directory.

State directory layout::

    journal.jsonl     accepted records only (append only, fsync after every line, hash chain)
    attempts.jsonl    every proposal with its outcome (used to resume and for the statistics)
    store/<sha256>.qasm | .qc     content-addressed copies of the proved candidates and of the references
    targets/<id>.json             target definitions (reference hash, model, generator parameters)
    best.json, BEST.md            the "best known" table; always regenerated from the journal

A record is valid only if its hash chain, its files, its costs and its proof can be re-derived: see verify.py.
"""
import fcntl
import hashlib
import json
import os
import time

from .costs import best_of, dominates, front_insert

GENESIS = '0' * 64


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def canon(rec):
    return json.dumps(rec, sort_keys=True, separators=(',', ':'), ensure_ascii=True)


def record_hash(rec):
    body = {k: v for k, v in rec.items() if k != 'hash'}
    return sha256_bytes(canon(body).encode())


def flags_hash(tdef):
    """SHA-256 of the flags that give the proofs of a target their meaning: the model (``clean``), the target spec string,
    the generator parameters and the reference hash.  Every new journal record carries it; ``verify`` compares it with the
    target definition on disk, so that flipping ``clean`` (or the spec) after the proofs is caught."""
    body = {'clean': bool(tdef.get('clean')), 'spec': tdef.get('spec'), 'gen': tdef.get('gen'), 'ref_sha': tdef.get('ref_sha')}
    return sha256_bytes(canon(body).encode())


def utc_now():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


class State:
    def __init__(self, root):
        self.root = os.path.abspath(root)
        self.journal = os.path.join(self.root, 'journal.jsonl')
        self.attempts = os.path.join(self.root, 'attempts.jsonl')
        self.store = os.path.join(self.root, 'store')
        self.targets = os.path.join(self.root, 'targets')
        self.work = os.path.join(self.root, 'work')

    def create(self):
        for d in (self.root, self.store, self.targets, self.work):
            os.makedirs(d, exist_ok=True)
        return self

    # ---- files in the store
    def put(self, data, suffix):
        """Store bytes under their sha256; returns the hash."""
        h = sha256_bytes(data)
        p = os.path.join(self.store, h + suffix)
        if not os.path.exists(p):
            tmp = p + '.tmp%d' % os.getpid()
            with open(tmp, 'wb') as f:
                f.write(data)
            os.replace(tmp, p)
        return h

    def path(self, h, suffix):
        return os.path.join(self.store, h + suffix)

    # ---- journal
    def read(self):
        """-> (records, problems). Unparsable lines are reported, not skipped silently."""
        recs, probs = [], []
        if not os.path.exists(self.journal):
            return recs, probs
        with open(self.journal, 'rb') as f:
            raw = f.read()
        lines = raw.split(b'\n')
        if lines and lines[-1] == b'':
            lines.pop()
        elif lines:
            probs.append('journal: the last line is not terminated (truncated write)')
        for i, ln in enumerate(lines):
            try:
                r = json.loads(ln.decode('utf8'))
                if not isinstance(r, dict):
                    raise ValueError('not an object')
                recs.append(r)
            except Exception as e:
                probs.append('journal line %d: unparsable (%s)' % (i + 1, e))
        return recs, probs

    def append(self, rec):
        """Add the next record (seq, prev, hash are set here). Caller holds the run lock."""
        recs, probs = self.read()
        if probs:
            raise RuntimeError('journal is damaged, refusing to append: ' + '; '.join(probs[:3]))
        prev = recs[-1]['hash'] if recs else GENESIS
        rec = dict(rec, seq=len(recs) + 1, prev=prev, ts=utc_now())
        rec['hash'] = record_hash(rec)
        with open(self.journal, 'ab') as f:
            f.write((canon(rec) + '\n').encode())
            f.flush()
            os.fsync(f.fileno())
        return rec

    def log_attempt(self, rec):
        rec = dict(rec, ts=utc_now())
        with open(self.attempts, 'ab') as f:
            f.write((json.dumps(rec, sort_keys=True) + '\n').encode())
            f.flush()
            os.fsync(f.fileno())

    def read_attempts(self):
        out = []
        if os.path.exists(self.attempts):
            for ln in open(self.attempts, 'rb').read().split(b'\n'):
                if ln.strip():
                    try:
                        out.append(json.loads(ln))
                    except ValueError:
                        pass          # a half-written last line after a kill: that attempt is simply repeated
        return out

    def lock(self):
        os.makedirs(self.root, exist_ok=True)
        f = open(os.path.join(self.root, 'run.lock'), 'a+')
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            f.close()
            raise RuntimeError('another loop run holds the lock of %s' % self.root)
        return f

    # ---- target definitions
    def save_target(self, tid, d):
        with open(os.path.join(self.targets, tid + '.json'), 'w') as f:
            json.dump(d, f, indent=1, sort_keys=True)

    def load_target(self, tid):
        with open(os.path.join(self.targets, tid + '.json')) as f:
            return json.load(f)

    def target_ids(self):
        if not os.path.isdir(self.targets):
            return []
        return sorted(x[:-5] for x in os.listdir(self.targets) if x.endswith('.json'))


def replay(records):
    """Re-derive the state from the records alone, checking every rule. -> (fronts, pinned, problems).
    fronts: target -> list of cost dicts; pinned: target -> ref sha256."""
    probs, fronts, pinned, flags = [], {}, {}, {}
    prev = GENESIS
    for i, r in enumerate(records):
        tag = 'record %d' % (i + 1)
        if r.get('seq') != i + 1:
            probs.append('%s: seq is %r (a record was removed, duplicated or reordered)' % (tag, r.get('seq')))
        if r.get('prev') != prev:
            probs.append('%s: hash chain broken (prev does not match the previous record)' % tag)
        if r.get('hash') != record_hash(r):
            probs.append('%s: record hash does not match its content (record was edited)' % tag)
        prev = r.get('hash', prev)
        t = r.get('target')
        role = r.get('role')
        after = r.get('cost_after')
        if not t or not isinstance(after, dict):
            probs.append('%s: no target or cost_after' % tag)
            continue
        fl = r.get('flags_sha')
        if t in flags:
            if flags[t] != fl:
                probs.append('%s: the flags hash of %s differs from the one of the earlier records of that target' % (tag, t))
        else:
            flags[t] = fl
        pv = r.get('prover') or {}
        if pv.get('status') != 'PROVED':
            probs.append('%s: not proved (prover status %r): a candidate without a proof is never accepted' % (tag, pv.get('status')))
        if t not in fronts:
            if role != 'baseline':
                probs.append('%s: first record of %s is not a baseline' % (tag, t))
            fronts[t] = []
            pinned[t] = r.get('ref_sha')
            if r.get('cost_before') is not None:
                probs.append('%s: baseline with a cost_before' % tag)
        else:
            if role != 'improvement':
                probs.append('%s: second baseline for %s' % (tag, t))
            if r.get('ref_sha') != pinned[t]:
                probs.append('%s: reference hash differs from the one pinned by the baseline of %s (proof of another circuit)' % (tag, t))
            if r.get('cost_before') != best_of(fronts[t]):
                probs.append('%s: cost_before %s is not the best known cost at that point (%s)'
                             % (tag, r.get('cost_before'), best_of(fronts[t])))
            want_kind = 'dominates' if r.get('cost_before') and dominates(after, r['cost_before']) else 'tradeoff'
            if r.get('kind') != want_kind:
                probs.append('%s: kind %r does not match the costs (%s)' % (tag, r.get('kind'), want_kind))
        fronts[t], ok = front_insert(fronts[t], after)
        if not ok:
            probs.append('%s: cost_after %s does not improve the front of %s' % (tag, after, t))
    return fronts, pinned, probs


def table_from(records):
    """The best-known table, from the journal only. -> dict target -> {best, baseline, front, records}."""
    out = {}
    for r in records:
        t = r['target']
        row = out.setdefault(t, {'baseline': r['cost_after'], 'records': 0, '_f': []})
        row['records'] += 1
        new, ok = front_insert([e['cost'] for e in row['_f']], r['cost_after'])
        if ok:
            row['_f'] = [e for e in row['_f'] if e['cost'] in new] + [{'cost': r['cost_after'], 'seq': r['seq'],
                                                                      'cand_sha': r['cand_sha']}]
    for row in out.values():
        row['front'] = row.pop('_f')
        row['best'] = best_of([e['cost'] for e in row['front']])
    return out


def write_table(state):
    recs, _ = state.read()
    tab = table_from(recs)
    with open(os.path.join(state.root, 'best.json'), 'w') as f:
        json.dump(tab, f, indent=1, sort_keys=True)
    lines = ['# Best known (generated from journal.jsonl; do not edit)', '',
             '| target | AND | T | CX | qubits | baseline AND / CX | journal records |', '|---|---|---|---|---|---|---|']
    for t in sorted(tab):
        b, a = tab[t]['best'], tab[t]['baseline']
        lines.append('| %s | %d | %d | %d | %d | %d / %d | %d |' % (t, b['AND'], b['T'], b['CX'], b['qubits'], a['AND'], a['CX'],
                                                                 tab[t]['records']))
    with open(os.path.join(state.root, 'BEST.md'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return tab
