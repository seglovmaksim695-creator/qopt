"""Import guard, recorder atoms, classical checks."""
import subprocess
import sys

from .conftest import needs_projectq


def test_import_without_projectq_is_clean():
    code = ("import sys; sys.modules['projectq'] = None\n"
            "from qoptct.foreign.projectq import check_projectq, ClassicalSpec\n"
            "try:\n    check_projectq(lambda eng: None, ClassicalSpec(lambda v: {}))\n    print('NO ERROR')\n"
            "except ImportError as e:\n    print('ERR', e)\n")
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, timeout=120)
    assert out.stdout.startswith('ERR') and 'pip install' in out.stdout, (out.stdout, out.stderr)


def _adder(n):
    def program(eng):
        from projectq.libs.math import AddQuantum
        a = eng.allocate_qureg(n)
        b = eng.allocate_qureg(n)
        AddQuantum | (a, b)
        return dict(inputs=dict(a=a, b=b), outputs=dict(a=a, b=b))
    return program


def _add_spec(n, extra=0):
    from qoptct.foreign.projectq import ClassicalSpec
    return ClassicalSpec(lambda v: dict(a=v['a'], b=(v['a'] + v['b'] + extra) % (1 << n)))


@needs_projectq
def test_adder_proved_and_independent_run_agrees():
    from qoptct.foreign.projectq import check_projectq
    r = check_projectq(_adder(3), _add_spec(3))
    assert r['status'] == 'PROVED' and r['native']['ok'] is True and r['native']['tested'] == 64
    assert r['circuit']['mode'] == 'nct' and r['assumptions']


@needs_projectq
def test_wrong_specification_is_refuted_with_a_witness():
    from qoptct.foreign.projectq import check_projectq
    r = check_projectq(_adder(3), _add_spec(3, extra=1))
    assert r['status'] == 'REFUTED' and set(r['witness']) == {'a', 'b'}


@needs_projectq
def test_atoms_are_exact_matrices():
    import math
    from projectq.ops import H, Ph, R, Rz, S, Sdag, T, Tdag, X, Z
    from qoptct.foreign.projectq._record import atom_of
    assert [atom_of(g) for g in (X, Z, H, S, Sdag, T, Tdag)] == ['x', 'z', 'h', 's', 'sdg', 't', 'tdg']
    assert atom_of(R(math.pi / 4)) == 't' and atom_of(R(math.pi / 2)) == 's'
    assert atom_of(Rz(math.pi / 4)) is None            # differs by a global phase: not accepted
    assert atom_of(Rz(0.3)) is None and atom_of(R(math.pi / 8)) is None and atom_of(Ph(math.pi)) is None


@needs_projectq
def test_rotations_and_measurements_are_unknown_never_proved():
    from qoptct.foreign.projectq import ClassicalSpec, check_projectq

    def rot(eng):
        from projectq.ops import Rz
        q = eng.allocate_qureg(1)
        Rz(0.3) | q[0]
        return dict(inputs=dict(q=q), outputs=dict(q=q))

    def meas(eng):
        from projectq.ops import Measure
        q = eng.allocate_qureg(1)
        Measure | q[0]
        return dict(inputs=dict(q=q), outputs=dict(q=q))
    for prog in (rot, meas):
        r = check_projectq(prog, ClassicalSpec(lambda v: dict(q=v['q'])))
        assert r['status'] == 'UNKNOWN' and 'unsupported' in r['detail']


@needs_projectq
def test_zero_control_state_and_unclean_ancilla():
    from qoptct.foreign.projectq import ClassicalSpec, check_projectq

    def zero_control(eng):
        from projectq.meta import Control
        from projectq.ops import X
        q = eng.allocate_qureg(1)
        t = eng.allocate_qureg(1)
        with Control(eng, q, ctrl_state='0'):
            X | t[0]
        return dict(inputs=dict(q=q, t=t), outputs=dict(q=q, t=t))
    r = check_projectq(zero_control, ClassicalSpec(lambda v: dict(q=v['q'], t=v['t'] ^ (1 - v['q']))), native=False)
    assert r['status'] == 'PROVED'

    def garbage(eng):
        from projectq.ops import CNOT
        a = eng.allocate_qureg(1)
        g = eng.allocate_qureg(1)
        CNOT | (a[0], g[0])
        return dict(inputs=dict(a=a), outputs=dict(a=a))
    r = check_projectq(garbage, ClassicalSpec(lambda v: dict(a=v['a'])), native=False)
    assert r['status'] == 'REFUTED' and 'must end in 0' in r['detail']


@needs_projectq
def test_logical_rules_and_text():
    from qoptct.foreign.projectq import record, to_logical
    lg = to_logical(record(_adder(2)))
    assert {op[0] for op in lg.ops} <= {'x', 'cx', 'and', 'erase'}
    q = lg.qasm()
    assert q.startswith('OPENQASM 3;') and q.count('measure') == sum(op[0] == 'erase' for op in lg.ops)
    assert sum(op[0] == 'and' for op in lg.ops) == sum(op[0] == 'erase' for op in lg.ops) > 0
    assert all(op[3] not in lg.free_wires for op in lg.ops if op[0] == 'and')
