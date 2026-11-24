"""ProjectQ -> qoptct.dyadic: exact angles, QFT / constant adder proofs, corrupted angles.  Needs ProjectQ."""
import math

import pytest

from .conftest import needs_projectq

pytestmark = needs_projectq


def _pq():
    from qoptct.foreign.projectq.dyadic import check_projectq_dyadic, exact_angle
    return check_projectq_dyadic, exact_angle


def _qft_prog(n, tweak=None):
    def program(eng):
        from projectq.ops import QFT
        q = eng.allocate_qureg(n)
        QFT | q
        if tweak:
            tweak(eng, q)
        return dict(registers=dict(q=q))
    return program


def _qft_spec(regs):
    from qoptct.dyadic import QFTSpec
    return QFTSpec(regs['q'], swaps=False)


def test_exact_angle_reading():
    _, exact_angle = _pq()
    from fractions import Fraction
    for k in (1, 2, 3, 10, 30, 38):
        assert exact_angle(round(math.pi / 2 ** k, 12), 2 * math.pi) == Fraction(1, 2 ** k)
    assert exact_angle(round(3 * math.pi / 2 ** 7, 12), 2 * math.pi) == Fraction(3, 2 ** 7)
    assert exact_angle(round(-math.pi / 8 % (2 * math.pi), 12), 2 * math.pi) == Fraction(15, 8)


@pytest.mark.parametrize('bad', [0.0, 1e-12, round(math.pi / 2 ** 43, 12), round(math.pi / 2 ** 39, 12),
                                 round(math.pi / 8 + 1e-9, 12), round(math.pi / 8 + 2e-12, 12), 0.3,
                                 2 * math.pi - 1e-12])
def test_inexact_or_ambiguous_angle_is_refused(bad):
    from qoptct.foreign.projectq._record import Unsupported
    _, exact_angle = _pq()
    with pytest.raises(Unsupported):
        exact_angle(bad, 2 * math.pi)


@pytest.mark.parametrize('n', [1, 2, 3, 5, 8])
def test_qft_proved(n):
    check, _ = _pq()
    r = check(_qft_prog(n), _qft_spec)
    assert r['status'] == 'PROVED', r
    assert r['angles']['count'] == n * (n - 1) // 2 or n == 1


def test_qft_reference_circuit_proved_and_inverse_is_not_qft():
    check, _ = _pq()
    from qoptct.dyadic import QFTSpec
    r = check(_qft_prog(4), lambda regs: QFTSpec(regs['q'], swaps=False, inverse=True))
    assert r['status'] == 'REFUTED', r
    assert r.get('witness') is not None


def test_add_constant_proved_and_wrong_constant_refuted():
    check, _ = _pq()
    from qoptct.dyadic import const_adder_spec

    def prog(c):
        def program(eng):
            from projectq.libs.math import AddConstant
            q = eng.allocate_qureg(4)
            AddConstant(c) | q
            return dict(registers=dict(b=q))
        return program
    r = check(prog(5), lambda regs: const_adder_spec(5, regs['b']))
    assert r['status'] == 'PROVED', r
    assert r['dyadic'].get('assumptions') or r['assumptions']
    r2 = check(prog(5), lambda regs: const_adder_spec(6, regs['b']))
    assert r2['status'] == 'REFUTED', r2


def test_corrupted_angle_to_another_dyadic_angle_is_refuted():
    check, _ = _pq()
    import projectq.ops as ops

    def tweak(eng, q):
        ops.R(math.pi / 2 ** 5) | q[0]         # an extra, exactly dyadic rotation
    r = check(_qft_prog(4, tweak), _qft_spec)
    assert r['status'] == 'REFUTED', r


@pytest.mark.parametrize('eps', [1e-9, 3e-12, 1e-11, 1e-6])
def test_corrupted_angle_not_dyadic_is_unknown(eps):
    check, _ = _pq()
    import projectq.ops as ops

    def tweak(eng, q):
        ops.R(math.pi / 8 + eps) | q[0]
    r = check(_qft_prog(4, tweak), _qft_spec)
    assert r['status'] == 'UNKNOWN', r
    assert 'not exactly' in r['detail'] or 'ambiguous' in r['detail'] or 'within' in r['detail'], r['detail']


def test_rounded_away_rotation_is_unknown():
    check, _ = _pq()
    r = check(_qft_prog(46), _qft_spec, budget_sec=5)     # QFT(46) has R(pi/2^45) which ProjectQ rounds to 0
    assert r['status'] == 'UNKNOWN', r
    assert r['status'] != 'PROVED'


def test_unsupported_measurement_is_unknown():
    check, _ = _pq()

    def program(eng):
        from projectq.ops import Measure
        q = eng.allocate_qureg(2)
        Measure | q
        return dict(registers=dict(q=q))
    r = check(program, _qft_spec)
    assert r['status'] == 'UNKNOWN', r


@pytest.mark.slow
def test_largest_exact_sizes_and_first_ambiguous_size():
    check, _ = _pq()
    from qoptct.dyadic import const_adder_spec
    assert check(_qft_prog(39), _qft_spec)['status'] == 'PROVED'      # R(pi/2^38) is the smallest unique angle
    assert check(_qft_prog(40), _qft_spec)['status'] == 'UNKNOWN'     # R(pi/2^39) = 6e-12: ambiguous

    def program(eng):
        from projectq.libs.math import AddConstant
        q = eng.allocate_qureg(38)
        AddConstant(0x2AAAAAAAAB) | q
        return dict(registers=dict(b=q))
    assert check(program, lambda regs: const_adder_spec(0x2AAAAAAAAB, regs['b']))['status'] == 'PROVED'
