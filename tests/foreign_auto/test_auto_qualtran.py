"""qoptct.foreign.auto on real Qualtran blocks (small sizes, seconds); n = 256 reruns are marked slow."""
import json
from pathlib import Path

import pytest

pytest.importorskip('qualtran')
pytest.importorskip('z3')

ROOT = Path(__file__).resolve().parents[2]


def test_registry_agrees_with_call_classically():
    from qoptct.foreign.auto.registry import selfcheck
    rows = selfcheck()
    bad = [r for r in rows if not r['ok']]
    assert not bad, bad[:3]
    assert len(rows) >= 60


@pytest.mark.parametrize('target', ['modadd', 'cmodadd', 'modsub', 'cmodsub', 'modneg', 'cmodneg', 'moddbl'])
def test_modular_blocks_small(target):
    from qoptct.foreign.auto.cli import prove_target
    v = prove_target(target, 6, timeout=600)
    assert v.get('verdict') == 'PROVED', v.get('reasons')
    assert v['end_to_end']['mismatches'] == 0
    assert v['qualtran_variant']['kind'] == 'original'


def test_montgomery_n16_has_the_programs_of_07():
    """the generic mechanism proves the same programs as the 0.7 client (fixed Qualtran, kind c)"""
    from qoptct.foreign.auto.cli import prove_target
    v = prove_target('montgomery', 16, timeout=900)
    assert v.get('verdict') == 'PROVED', v.get('reasons')
    old = json.loads((ROOT / 'results' / 'qualtran' / 'montgomery_n16_w1.json').read_text())
    assert v['program']['ops_sha256'] == old['program']['ops_sha256']
    hashes = {p['type']: p['qasm_sha256'] for p in v['pieces']}
    assert hashes['SingleWindowModMul'] == old['program']['generic_window_qasm_sha256']
    sub = {s['label']: s['qasm_sha256'] for s in old['generic_window']['subnodes']}
    assert hashes['CAdd'] == sub['cadd0'] and hashes['Xor'] == sub['xor']
    seg = {s['label']: s['qasm_sha256'] for s in old['segments']}
    assert hashes['LessThanConstant'] == seg['lt']
    assert v['qualtran_variant']['kind'] == 'c'


def test_loop_without_invariant_is_not_proved():
    from qoptct.foreign.auto.cli import prove_target
    v = prove_target('montgomery', 8, timeout=900, extra=['--drop-invariant'])
    assert v.get('verdict') != 'PROVED'
    assert any('loop' in r or 'precondition of the loop body' in r for r in v['reasons']), v['reasons']


@pytest.mark.parametrize('target', ['montgomery', 'cmodadd'])
def test_canaries_small(target):
    from qoptct.foreign.auto.canary import canaries_child
    r = canaries_child(target, 8, timeout=1200)
    assert r.get('proved') == [], r.get('proved')
    assert r['counted'] >= 20 and r['joints'] >= 10, (r['counted'], r['joints'])
    assert r['ok']


@pytest.mark.slow
def test_montgomery_n256_has_the_program_of_07():
    from qoptct.foreign.auto.cli import prove_target
    v = prove_target('montgomery', 256, timeout=3600)
    assert v.get('verdict') == 'PROVED', v.get('reasons')
    old = json.loads((ROOT / 'results' / 'qualtran' / 'montgomery_n256_w1.json').read_text())
    assert v['program']['ops_sha256'] == old['program']['ops_sha256']
    hashes = {p['type']: p['qasm_sha256'] for p in v['pieces']}
    assert hashes['SingleWindowModMul'] == old['program']['generic_window_qasm_sha256']


@pytest.mark.slow
@pytest.mark.parametrize('target', ['modadd', 'cmodadd', 'modsub', 'cmodsub', 'modneg', 'cmodneg', 'moddbl'])
def test_modular_blocks_n256(target):
    from qoptct.foreign.auto.cli import prove_target
    v = prove_target(target, 256, timeout=3600)
    assert v.get('verdict') == 'PROVED', v.get('reasons')
