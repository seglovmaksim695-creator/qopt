"""Binary_ECC circuits (need a checkout: $QOPT_BINARY_ECC) and the recorded evidence."""
import json

import pytest

from .conftest import RESULTS, binary_ecc_source, needs_binary_ecc, needs_projectq


@needs_projectq
@needs_binary_ecc
def test_multiplier_m8_proved_with_pin_and_canaries():
    from qoptct.foreign.projectq import results
    r = results.prove_multiplier(8, binary_ecc_source(), canaries=20)
    assert r['status'] == 'PROVED' and r['garbage']['wires'] == 19
    c = r['canaries']
    assert c['ok'] and c['false_proved'] == [] and c['decided'] >= 20


@needs_projectq
@needs_binary_ecc
def test_point_addition_m8_whole_domain():
    from qoptct.foreign.projectq import results
    r = results.prove_point_addition(8, binary_ecc_source(), canaries=20)
    assert r['status'] == 'PROVED' and r['garbage']['wires'] == 364 and r['canaries']['ok']
    assert 'q = 0 or x1 != x2' in r['claim'] or 'x1 != x2' in r['claim']


@pytest.mark.skipif(not RESULTS.is_dir(), reason='results/projectq is not part of this tree')
def test_recorded_evidence_has_clean_canaries():
    files = sorted(RESULTS.glob('*.json'))
    assert files
    for f in files:
        d = json.loads(f.read_text())
        if d['status'] == 'PROVED':
            c = d['canaries']
            assert c and c['ok'] and not c['false_proved'] and c['decided'] >= 20, f.name
