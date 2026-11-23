"""Montgomery proof client: checks that need neither Qualtran nor a prover run (API, arguments, identities, mutant tags)."""
import pytest

from qoptct.foreign.qualtran import montgomery as M
from qoptct.foreign.qualtran.montgomery import mutants as MUT


def test_prime_for():
    assert M.prime_for(256) == 2 ** 255 - 19
    for n in (8, 16, 32, 64):
        p = M.prime_for(n)
        assert p % 2 == 1 and 2 * p <= 1 << n and p > 1 << (n - 2)


@pytest.mark.parametrize("n,w", [(16, 1), (16, 2), (16, 4), (64, 4), (256, 1), (256, 2)])
def test_telescoping_identity(n, w):
    assert M.telescoping(n, w) is True


def test_bad_arguments():
    with pytest.raises(ValueError):
        M.prove_montgomery(15, 2)           # w does not divide n
    with pytest.raises(ValueError):
        M.prove_montgomery(16, 2, 2 ** 15)  # even modulus
    with pytest.raises(ValueError):
        M.prove_montgomery(16, 2, 40000)    # 2p > 2^n


def test_mutant_tags_known():
    assert set(MUT.TAGS) == {"H1", "H2", "M1", "M1b", "M2", "M3", "M4", "M5"}
    with pytest.raises(ValueError):
        MUT.apply("nope", None, [], [], 0, 0)
