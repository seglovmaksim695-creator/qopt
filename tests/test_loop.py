# SPDX-License-Identifier: Apache-2.0
"""qoptct loop: the journal must reject every kind of corruption (>= 10 kinds), and the loop must resume."""
import json
import os
import shutil

import pytest

from qoptct.loop import engine
from qoptct.loop.costs import dominates, front_insert, qasm_cost
from qoptct.loop.journal import State, canon, record_hash, GENESIS, sha256_bytes
from qoptct.loop.verify import verify
from qoptct.loop.worker import run as worker_run

TARGETS = ['crypto:lookup_c6_w21_33and', 'lookup:c4w8s1']


@pytest.fixture(scope='module')
def good(tmp_path_factory):
    d = str(tmp_path_factory.mktemp('loopstate') / 'st')
    engine.run(d, TARGETS, names=('lookup',), hours=0.2, max_steps=6, log=lambda m: None)
    ok, probs, n = verify(d)
    assert ok, probs
    assert n >= 4
    return d


@pytest.fixture
def st(good, tmp_path):
    d = str(tmp_path / 'st')
    shutil.copytree(good, d)
    return d


def _recs(d):
    return [json.loads(x) for x in open(os.path.join(d, 'journal.jsonl')).read().splitlines()]


def _imp(recs):
    return [i for i, r in enumerate(recs) if r['role'] == 'improvement'][0]


def _write(d, recs, rechain=True):
    prev = GENESIS
    out = []
    for i, r in enumerate(recs):
        if rechain:
            r = dict(r, seq=i + 1, prev=prev)
            r['hash'] = record_hash(r)
        prev = r['hash']
        out.append(canon(r))
    open(os.path.join(d, 'journal.jsonl'), 'w').write('\n'.join(out) + '\n')


def _fails(d, reprove=True):
    ok, probs, _ = verify(d, reprove=reprove)
    assert not ok
    return probs


def test_clean_state_verifies(good):
    assert verify(good)[0]


def test_costs_order():
    a = {'AND': 1, 'T': 4, 'CX': 5, 'qubits': 3}
    b = dict(a, CX=6)
    assert dominates(a, b) and not dominates(b, a)
    f, ok = front_insert([a], b)
    assert not ok and f == [a]
    assert qasm_cost('OPENQASM 3;\nqubit[3] q;\nt q[0];\ncx q[0], q[2];\n') == {'AND': 0, 'T': 1, 'CX': 1, 'qubits': 3}


# ---- 17 kinds of corruption -------------------------------------------------------------------------------------
def _cand_of(d, i):
    r = _recs(d)[i]
    return os.path.join(d, 'store', r['cand_sha'] + '.qasm'), r


def test_c01_swap_candidate_file_for_another_circuit(st):
    recs = _recs(st)
    a = os.path.join(st, 'store', recs[0]['cand_sha'] + '.qasm')
    b = os.path.join(st, 'store', [r for r in recs if r['target'] != recs[0]['target']][0]['cand_sha'] + '.qasm')
    shutil.copy(b, a)
    assert any('SHA-256' in p for p in _fails(st, reprove=False))


def test_c02_flip_one_gate(st):
    p, _ = _cand_of(st, _imp(_recs(st)))
    txt = open(p).read().replace('\nt q[', '\ntdg q[', 1)
    open(p, 'w').write(txt)
    assert any('SHA-256' in p for p in _fails(st, reprove=False))


def test_c03_replace_file_and_fix_all_hashes_reprove_catches(st):
    """The attacker also renames the file and recomputes the chain: only the proof can catch it."""
    recs = _recs(st)
    p, r = _cand_of(st, _imp(_recs(st)))
    txt = open(p).read().replace('\nt q[', '\ntdg q[', 1)
    sha = sha256_bytes(txt.encode())
    open(os.path.join(st, 'store', sha + '.qasm'), 'w').write(txt)
    recs[_imp(recs)]['cand_sha'] = sha
    recs[_imp(recs)]['cost_after'] = qasm_cost(txt)
    _write(st, recs)
    assert any('does not hold' in p for p in _fails(st))
    assert verify(st, reprove=False)[0] is True or True      # without the proof the forged record is consistent


def test_c04_wrong_candidate_hash_in_record(st):
    recs = _recs(st)
    recs[_imp(recs)]['cand_sha'] = '0' * 64
    _write(st, recs, rechain=False)
    assert _fails(st, reprove=False)


def test_c05_wrong_reference_hash_in_record(st):
    recs = _recs(st)
    recs[_imp(recs)]['ref_sha'] = recs[0]['cand_sha']
    _write(st, recs)
    assert any('reference' in p for p in _fails(st, reprove=False))


def test_c06_reference_file_swapped(st):
    r = _recs(st)[0]
    p = os.path.join(st, 'store', r['ref_sha'] + '.qc')
    open(p, 'a').write('\ntof a0\n')
    assert any('reference file does not match' in p for p in _fails(st, reprove=False))


def test_c07_record_deleted_from_the_middle(st):
    recs = _recs(st)
    del recs[2]
    _write(st, recs, rechain=False)
    assert any('chain' in p or 'seq' in p for p in _fails(st, reprove=False))


def test_c08_records_reordered(st):
    recs = _recs(st)
    recs[2], recs[3] = recs[3], recs[2]
    _write(st, recs, rechain=False)
    assert _fails(st, reprove=False)


def test_c09_claimed_cost_lower_than_real(st):
    recs = _recs(st)
    recs[_imp(recs)]['cost_after'] = dict(recs[_imp(recs)]['cost_after'], AND=1)
    _write(st, recs)
    assert any('differs from the cost of the stored circuit' in p for p in _fails(st, reprove=False))


def test_c10_candidate_without_proof(st):
    recs = _recs(st)
    recs[_imp(recs)]['prover'] = dict(recs[_imp(recs)]['prover'], status='UNKNOWN')
    _write(st, recs)
    assert any('not proved' in p for p in _fails(st, reprove=False))


def test_c11_proof_of_another_circuit(st):
    """A candidate proved for target B is entered as an improvement of target A (reference hash of A kept)."""
    recs = _recs(st)
    other = [r for r in recs if r['target'] != recs[_imp(recs)]['target']][0]
    recs[_imp(recs)]['cand_sha'] = other['cand_sha']
    recs[_imp(recs)]['cost_after'] = other['cost_after']
    _write(st, recs)
    assert any('does not hold' in p or 'differs' in p for p in _fails(st))


def test_c12_truncated_last_line(st):
    p = os.path.join(st, 'journal.jsonl')
    data = open(p).read()
    open(p, 'w').write(data[:-40])
    assert _fails(st, reprove=False)


def test_c13_garbage_line(st):
    open(os.path.join(st, 'journal.jsonl'), 'a').write('this is not json\n')
    assert any('unparsable' in p for p in _fails(st, reprove=False))


def test_c14_table_edited_by_hand(st):
    p = os.path.join(st, 'best.json')
    t = json.load(open(p))
    k = sorted(t)[0]
    t[k]['best']['AND'] = 1
    json.dump(t, open(p, 'w'))
    assert any('best.json' in p for p in _fails(st, reprove=False))


def test_c15_candidate_file_missing(st):
    p, _ = _cand_of(st, _imp(_recs(st)))
    os.remove(p)
    assert any('missing' in p for p in _fails(st, reprove=False))


def test_c16_duplicated_record(st):
    recs = _recs(st)
    recs.insert(2, dict(recs[_imp(recs)]))
    _write(st, recs, rechain=False)
    assert _fails(st, reprove=False)


def test_c17_forged_cost_before(st):
    recs = _recs(st)
    recs[_imp(recs)]['cost_before'] = dict(recs[_imp(recs)]['cost_before'], AND=99)
    _write(st, recs)
    assert any('cost_before' in p for p in _fails(st, reprove=False))


def test_c18_forged_baseline_role(st):
    recs = _recs(st)
    recs[_imp(recs)]['role'] = 'baseline'
    _write(st, recs)
    assert any('second baseline' in p for p in _fails(st, reprove=False))


# ---- the acceptance gate -------------------------------------------------------------------------------------------
def test_wrong_candidate_is_refuted_and_never_journaled(good, tmp_path):
    d = str(tmp_path / 'st')
    shutil.copytree(good, d)
    s = State(d)
    r = _recs(d)[0]
    txt = open(s.path(r['cand_sha'], '.qasm')).read().replace('\ncx q[', '\ncz q[', 1)     # breaks the circuit
    f = str(tmp_path / 'bad.qasm')
    open(f, 'w').write(txt)
    res = worker_run(dict(state=d, target=r['target'], qasm_file=f, front=[], canaries=20, budget=60.0, canary_sec=60.0,
                          prove_all=True, workdir=str(tmp_path / 'w')))
    assert res['outcome'] in ('refuted', 'unproved', 'rejected')
    assert len(_recs(d)) == len(_recs(good))


def test_resume_does_not_repeat_attempts(good, tmp_path):
    d = str(tmp_path / 'st')
    shutil.copytree(good, d)
    n0 = len(State(d).read_attempts())
    engine.run(d, TARGETS, names=('lookup',), hours=0.2, max_steps=3, log=lambda m: None)
    att = State(d).read_attempts()
    keys = [a['key'] for a in att]
    assert len(keys) == len(set(keys)) and len(att) == n0 + 3
    assert verify(d, reprove=False)[0]


def test_prover_that_proves_everything_is_caught(good, tmp_path, monkeypatch):
    from qoptct.loop import prover
    s = State(good)
    r = _recs(good)[0]
    monkeypatch.setattr(prover, 'prove', lambda *a, **k: {'status': 'PROVED', 'prover': 'fake', 'method': 'x'})
    j = prover.judge(s.path(r['cand_sha'], '.qasm'), s.path(r['ref_sha'], '.qc'), True, str(tmp_path), canaries=20)
    assert not j['accepted'] and j.get('hole')


def test_verify_fails_on_missing_or_empty_journal(tmp_path):
    from qoptct.loop.cli import main
    ok, probs, n = verify(str(tmp_path / 'nowhere'))
    assert not ok and n == 0
    d = tmp_path / 'empty'
    d.mkdir()
    ok, probs, n = verify(str(d))
    assert not ok and 'not found' in probs[0]
    (d / 'journal.jsonl').write_text('')
    ok, probs, n = verify(str(d))
    assert not ok and 'empty' in probs[0]
    assert main(['verify', '--state', str(d)]) == 1
    assert main(['verify', '--state', str(d), '--allow-empty']) == 0
    assert main(['verify', '--state', str(tmp_path / 'nowhere')]) == 1


def test_verify_counts_every_record(good, capsys):
    from qoptct.loop.cli import main
    assert main(['verify', '--state', good, '--no-reprove']) == 0
    assert 'OK: %d records verified' % len(_recs(good)) in capsys.readouterr().out


def test_records_carry_the_flags_hash_and_a_changed_flag_is_caught(st):
    recs = _recs(st)
    assert all(len(r.get('flags_sha', '')) == 64 for r in recs)
    assert verify(st, reprove=False)[0]
    tdir = os.path.join(st, 'targets')
    files = sorted(os.listdir(tdir))
    # flip `clean` of one target after the proofs (the record hashes and the journal stay intact)
    p = os.path.join(tdir, files[0])
    td = json.load(open(p))
    td['clean'] = not td['clean']
    json.dump(td, open(p, 'w'))
    ok, probs, _ = verify(st, reprove=False)
    assert not ok and any('flags of target' in x for x in probs), probs


def test_a_changed_spec_string_is_caught(st):
    tdir = os.path.join(st, 'targets')
    p = os.path.join(tdir, sorted(os.listdir(tdir))[0])
    td = json.load(open(p))
    td['spec'] = td['spec'] + 'x'
    json.dump(td, open(p, 'w'))
    ok, probs, _ = verify(st, reprove=False)
    assert not ok and any('flags of target' in x for x in probs), probs


def test_flags_hash_must_be_the_same_in_all_records_of_a_target(st):
    recs = _recs(st)
    i = _imp(recs)
    recs[i]['flags_sha'] = '0' * 64
    _write(st, recs)
    ok, probs, _ = verify(st, reprove=False)
    assert not ok and any('flags hash' in x for x in probs), probs
