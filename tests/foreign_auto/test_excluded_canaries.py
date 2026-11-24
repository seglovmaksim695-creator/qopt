"""D13 / D19: a candidate canary that a filter drops is recorded with its reason, never dropped silently."""
import json
import random
from types import SimpleNamespace

import pytest

pytest.importorskip('z3')
from qoptct.foreign.auto import canary, registry


def test_swap_filter_records_why_it_dropped_the_candidate(monkeypatch):
    monkeypatch.setattr(registry, 'contract_for', lambda obj: None)
    del canary.EXCLUDED[:]
    k = SimpleNamespace(obj=object(), ins={'x': [0, 1]})
    assert canary._swap_effective(k, 0, random.Random(1)) is False
    assert len(canary.EXCLUDED) == 1
    e = canary.EXCLUDED[0]
    assert e['filter'] == 'D13' and 'no registry contract' in e['reason'] and e['candidate']


def test_leaf_filter_records_why_it_dropped_the_candidate(monkeypatch):
    monkeypatch.setattr(registry, 'contract_for', lambda obj: None)
    del canary.EXCLUDED[:]
    lf = SimpleNamespace(obj=object(), type='Leaf', span=(3, 9))
    assert canary._leaf_differs(None, lf, lambda o: None, random.Random(1)) is False
    assert canary.EXCLUDED[0]['filter'] == 'D19' and canary.EXCLUDED[0]['site']


def test_sidecar_file_is_written(tmp_path):
    out = tmp_path / 'c.json'
    canary._write_excluded(out, {'excluded': [dict(filter='D13', candidate='c', reason='r', site='s:1')]})
    lines = (tmp_path / 'c.json.excluded.jsonl').read_text().splitlines()
    assert json.loads(lines[0])['reason'] == 'r'
