"""The registry of the results of qoptct.foreign.auto matches the result files (no Qualtran needed)."""
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


def test_registry_is_a_fresh_build_of_the_results():
    if not (ROOT / 'qoptct' / 'foreign' / 'auto' / 'data' / 'registry.json').exists():
        pytest.skip('the source tree is not here (tests run against an installed wheel)')
    r = subprocess.run([sys.executable, str(ROOT / 'scripts' / 'build_auto_registry.py'), '--check'],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_registry_sources_and_commands():
    from qoptct.foreign.auto import results as RR
    es = RR.entries()
    assert es
    assert RR.verify_sources(ROOT / 'results' / 'auto') == []
    for e in es:
        assert e['reproduce'].startswith('python -m qoptct.foreign.auto prove ')
        assert e['status'] in ('PROVED', 'NOT_PROVED')
    assert RR.format_table(es).count('\n') == len(es) + 1


def test_criterion_is_not_rounded_up():
    from qoptct.foreign.auto.results import criterion
    e = dict(status='PROVED', end_to_end=dict(cases=16, mismatches=0), source_sha256='s',
             canaries=dict(proved=[], counted=30, joints=16, source_sha256='s'))
    assert criterion(e) == (True, [])
    assert not criterion(dict(e, canaries=None))[0]
    assert not criterion(dict(e, canaries=dict(e['canaries'], joints=9)))[0]
    assert not criterion(dict(e, canaries=dict(e['canaries'], proved=['J1'])))[0]
    assert not criterion(dict(e, canaries=dict(e['canaries'], source_sha256='other')))[0]
    assert not criterion(dict(e, end_to_end=dict(cases=16, mismatches=1)))[0]
    assert not criterion(dict(e, status='NOT_PROVED'))[0]
