import pytest

from qoptct.foreign.qualtran._intexpr import int_expr


def test_values():
    assert int_expr('2**255-19') == 2**255 - 19
    assert int_expr(' 65521 ') == 65521
    assert int_expr('(1<<64)-59') == 2**64 - 59
    assert int_expr(32749) == 32749


@pytest.mark.parametrize('bad', ['__import__("os")', 'p', 'abs(3)', '9**9**9', '2**(10**9)', '1/2', '[1]', 'True', ''])
def test_refused(bad):
    with pytest.raises(ValueError):
        int_expr(bad)


def test_zglue_refuses_a_modulus_outside_the_domain():
    pytest.importorskip('qualtran')
    from qoptct.foreign.qualtran.ecadd import zglue
    assert zglue.main(['--n', '16', '--p', '32767'])['status'] == 'NOT_PROVED'   # 32767 = 7*31*151
