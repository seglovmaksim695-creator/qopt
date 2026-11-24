"""Library blocks, Clifford+T comparison, canaries, GF(2^m) multipliers and the garbage pin."""
import pytest

from .conftest import needs_projectq
from .test_record_check import _add_spec, _adder


@needs_projectq
@pytest.mark.parametrize('name,n', [('add_quantum', 4), ('add_quantum_carry', 3), ('add_quantum_controlled', 3),
                                   ('subtract_quantum', 4), ('divide_quantum', 2), ('multi_controlled_x', 4),
                                   ('controlled_swap', 1), ('comparator_quantum', 2)])
def test_library_blocks_proved(name, n):
    from qoptct.foreign.projectq import check_projectq
    from qoptct.foreign.projectq.library import block
    r = check_projectq(*block(name, n))
    assert r['status'] == 'PROVED', r['detail']
    assert r['native']['ok'] is True


@needs_projectq
def test_projectq_comparator_decomposition_is_wrong_except_n2():
    from qoptct.foreign.projectq import check_projectq
    from qoptct.foreign.projectq.library import block
    assert check_projectq(*block('comparator_quantum', 1))['status'] == 'REFUTED'
    r = check_projectq(*block('comparator_quantum', 3))
    assert r['status'] == 'REFUTED'
    assert r['native']['decomposed_run'] != r['native']['specification']


@needs_projectq
def test_clifford_t_adder_is_proved_equal_to_nct_and_mutants_are_not():
    from qoptct.foreign.projectq import check_projectq, record, to_logical
    from qoptct.foreign.projectq._ct import check_clifford_t
    r = check_projectq(_adder(2), _add_spec(2), mode='clifford_t')
    assert r['status'] == 'PROVED', r['detail']
    ref = record(_adder(2), 'nct', None)
    lg = to_logical(record(_adder(2), 'clifford_t'))
    ts = [i for i, op in enumerate(lg.ops) if op[0] in ('t', 'tdg')]
    assert len(ts) >= 8
    for i in ts[:6]:
        for new in ([], [('tdg' if lg.ops[i][0] == 't' else 't',) + lg.ops[i][1:]], [('s',) + lg.ops[i][1:]]):
            bad = lg.copy(lg.ops[:i] + new + lg.ops[i + 1:])
            assert check_clifford_t(bad, ref)['status'] != 'PROVED', (i, new)


@needs_projectq
def test_canaries_of_an_adder_are_all_not_proved():
    from qoptct.foreign.projectq import record
    from qoptct.foreign.projectq.canaries import run_canaries
    rows = run_canaries(record(_adder(4)), _add_spec(4), None, count=20)
    decided = [r for r in rows if r['status'] != 'EQUIVALENT']
    assert len(decided) >= 20 and all(r['status'] != 'PROVED' for r in decided)
    assert {'junction', 'erasure', 'phase', 'rare'} <= {r['group'] for r in decided}


def _schoolbook(m, poly):
    def program(eng):
        from projectq.ops import Toffoli
        from qoptct.foreign.projectq._gf2 import poly_mod
        x = eng.allocate_qureg(m)
        y = eng.allocate_qureg(m)
        r = eng.allocate_qureg(m)
        for i in range(m):
            for j in range(m):
                rem = poly_mod(1 << (i + j), poly)
                for k in range(m):
                    if (rem >> k) & 1:
                        Toffoli | (x[i], y[j], r[k])
        return dict(inputs=dict(x=x, y=y), outputs=dict(x=x, y=y, result=r))
    return program


@needs_projectq
def test_gf2_multiplier_exact2_proof_and_wrong_polynomial():
    from qoptct.foreign.projectq import GF2MulSpec, check_projectq
    r = check_projectq(_schoolbook(5, 0b100101), GF2MulSpec(5, 0b100101))
    assert r['status'] == 'PROVED' and r['exact2']['status'] == 'PROVED' and r['independent']['ok'] and r['native']['ok']
    assert check_projectq(_schoolbook(5, 0b100101), GF2MulSpec(5, 0b101001))['status'] == 'REFUTED'
    r = check_projectq(_schoolbook(5, 0b100101), GF2MulSpec(5, 0b100001))     # x^5 + 1 is reducible
    assert r['status'] == 'UNKNOWN' and 'irreducible' in r['detail']


def test_rabin_test():
    from qoptct.foreign.projectq._gf2 import rabin_irreducible

    def ex(*e):
        return sum(1 << i for i in e)
    assert rabin_irreducible(0x11b) and rabin_irreducible(ex(163, 7, 6, 3, 0)) and rabin_irreducible(ex(233, 74, 0))
    assert rabin_irreducible(ex(283, 12, 7, 5, 0))
    assert not rabin_irreducible(0x11a) and not rabin_irreducible(ex(8, 0)) and not rabin_irreducible(ex(6, 0))


@needs_projectq
def test_garbage_pin_refutes_a_changed_garbage_wire():
    from projectq.ops import CNOT, Toffoli
    from qoptct.foreign.projectq import GF2MulSpec, check_projectq, pin_garbage, record

    def program(eng):
        x = eng.allocate_qureg(2)
        y = eng.allocate_qureg(2)
        r = eng.allocate_qureg(2)
        g = eng.allocate_qureg(1)
        for (i, j), tg in {(0, 0): [0], (0, 1): [1], (1, 0): [1], (1, 1): [0, 1]}.items():     # x^2 = x + 1 over GF(4)
            for k in tg:
                Toffoli | (x[i], y[j], r[k])
        CNOT | (x[0], g[0])
        return dict(inputs=dict(x=x, y=y), outputs=dict(x=x, y=y, result=r))
    spec = GF2MulSpec(2, 0b111)
    assert check_projectq(program, spec)['status'] == 'REFUTED'
    pin = pin_garbage(program, spec)
    assert len(pin) == 1 and check_projectq(program, spec, pin=pin)['status'] == 'PROVED'
    rec = record(program)
    ops = list(rec.ops) + [('x', (rec.inputs['y'][0],), rec.n_wires - 1)]
    assert check_projectq(rec.copy(ops), spec, pin=pin)['status'] == 'REFUTED'
