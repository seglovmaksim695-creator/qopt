"""The registry of proved ProjectQ results: consistent with the result files. No ProjectQ needed."""
import importlib.util

import pytest

from qoptct.foreign.projectq import registry

from .conftest import REPO, RESULTS

needs_results = pytest.mark.skipif(not RESULTS.is_dir(), reason='results/projectq is not part of this tree')


def test_registry_shape():
    reg = registry.load_registry()
    assert reg['schema'] == 'qoptct.projectq.registry/1'
    es = reg['entries']
    assert es and len({e['id'] for e in es}) == len(es)
    assert {e['name'] for e in es} <= set(registry.NAMES)
    for e in es:
        assert e['status'] == 'PROVED' and e['assumptions'] and e['claim']
        c = e['canaries']
        assert c['not_proved'] == c['decided'] >= 20
        assert len(e['program']['value']) == 64
        assert e['reproduce'].startswith('python -m qoptct.foreign.projectq prove ')
    for f in reg['findings']:
        assert f['status'] == 'REFUTED'


def test_binary_ecc_results_cover_the_target_sizes():
    ids = {e['id'] for e in registry.entries()}
    for m in (8, 16, 127, 163, 233):
        assert f'binary_ecc_multiplier/m{m}' in ids
    assert 'binary_ecc_point_addition/FLT_in_m8' in ids
    for e in registry.entries('binary_ecc_multiplier'):
        assert e['source']['sha256'] and e['garbage']['pin_sha256']


def test_table_and_long_forms():
    es = registry.entries()
    assert 'binary_ecc_multiplier/m8' in registry.format_table(es)
    assert 'assumptions:' in registry.format_entry(es[0])


@needs_results
def test_registry_matches_result_files():
    assert registry.verify_sources(RESULTS) == []


@needs_results
def test_registry_is_up_to_date():
    spec = importlib.util.spec_from_file_location('build_projectq_registry', REPO / 'scripts' / 'build_projectq_registry.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.main(['--check']) == 0


@needs_results
def test_every_result_is_registered_or_a_finding():
    reg = registry.load_registry()
    known = {e['evidence']['file'] for e in reg['entries']} | {f['evidence']['file'] for f in reg['findings']}
    unregistered = {p.name for p in RESULTS.glob('*.json')} - known
    assert not unregistered, unregistered      # e.g. UNKNOWN results or failed canary criteria must be looked at
