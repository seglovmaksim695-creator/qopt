"""Unified command line of 0.8: qoptct results / dyadic / qiskit / projectq / auto / loop dispatch, exit codes, registries."""
import json
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]

QFT3 = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
h q[2];
cp(pi/2) q[1],q[2];
cp(pi/4) q[0],q[2];
h q[1];
cp(pi/2) q[0],q[1];
h q[0];
swap q[0],q[2];
"""


def run(*args):
    """The command line through the module entry point of qoptct.qopt_core (what the `qoptct` script calls)."""
    code = 'import sys; from qoptct.qopt_core import main; sys.argv[0] = "qoptct"; main()'
    return subprocess.run([sys.executable, '-c', code, *args], capture_output=True, text=True, cwd=str(ROOT), timeout=300)


def test_results_table_has_every_library():
    r = run('results', '--json')
    assert r.returncode == 0, r.stderr
    rows = json.loads(r.stdout)
    libs = {x['library'] for x in rows}
    assert libs == {'qualtran', 'qiskit', 'projectq', 'dyadic', 'auto', 'loop'}
    assert all(x['status'].startswith('PROVED') for x in rows)
    assert all(x['reproduce'] for x in rows)


def test_results_library_filter_and_text():
    r = run('results', '--library', 'dyadic')
    assert r.returncode == 0 and 'qft' in r.stdout and 'qualtran' not in r.stdout.split('\n\n')[0]


@pytest.mark.parametrize('lib,cmd', [('dyadic', 'results'), ('projectq', 'results'), ('auto', 'results'),
                                     ('qualtran', 'results')])
def test_library_results_commands(lib, cmd):
    r = run(lib, cmd)
    assert r.returncode == 0, r.stderr
    assert 'PROVED' in r.stdout


def test_qiskit_list_and_fractional():
    pytest.importorskip('qiskit')
    assert run('qiskit', 'list').returncode == 0
    r = run('qiskit', 'fractional')
    assert r.returncode == 0 and 'QFT' in r.stdout


def test_dyadic_prove_exit_codes(tmp_path):
    good = tmp_path / 'qft3.qasm'
    good.write_text(QFT3)
    bad = tmp_path / 'bad.qasm'
    bad.write_text(QFT3.replace('cp(pi/4)', 'cp(pi/8)'))
    assert run('dyadic', 'prove', str(good), '--qft').returncode == 0
    r = run('dyadic', 'prove', str(bad), '--qft')
    assert r.returncode == 1 and 'REFUTED' in r.stdout
    nd = tmp_path / 'nd.qasm'
    nd.write_text(QFT3.replace('cp(pi/4)', 'cp(0.3)'))
    r = run('dyadic', 'prove', str(nd), '--qft')
    assert r.returncode == 2, r.stdout          # not a multiple of pi/2^k: UNKNOWN, never PROVED
    assert run('dyadic', 'prove', str(tmp_path / 'missing.qasm'), '--qft').returncode == 4


def test_help_lists_the_new_commands():
    r = run('--help')
    for name in ('qiskit', 'projectq', 'dyadic', 'auto', 'results', 'loop', 'qualtran'):
        assert name in r.stdout, name


def test_loop_verify_of_the_shipped_journal():
    state = ROOT / 'results' / 'loop' / 'state'
    if not (state / 'journal.jsonl').exists():
        pytest.skip('results/ is not part of the installed package')
    r = run('loop', 'verify', '--state', str(state), '--no-reprove')
    assert r.returncode == 0 and 'OK: 116 records verified' in r.stdout, r.stdout + r.stderr


@pytest.mark.parametrize('script', ['build_dyadic_registry.py', 'build_loop_registry.py'])
def test_registry_matches_results(script):
    if not (ROOT / 'results' / 'dyadic').is_dir() or not (ROOT / 'qoptct' / 'dyadic' / 'data').is_dir():
        pytest.skip('the source tree is not here (tests run against an installed wheel)')
    r = subprocess.run([sys.executable, str(ROOT / 'scripts' / script), '--check'], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
