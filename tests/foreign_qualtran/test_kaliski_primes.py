"""Primality evidence of the modulus."""
import pytest

from qoptct.foreign.qualtran.kaliski.primes import certify_prime, lucas_lehmer, prime_accepted

P256 = 2 ** 256 - 2 ** 224 + 2 ** 192 + 2 ** 96 - 1


@pytest.mark.parametrize("p", [13, 251, 4093, 32749, 2 ** 31 - 1])
def test_small_primes_by_trial_division(p):
    r = certify_prime(p)
    assert r["verdict"] == "PROVEN PRIME" and "trial division" in r["how"] and prime_accepted(r)


@pytest.mark.parametrize("m", [1, 4, 21, 561, 32747 * 3, 2 ** 31 + 1, 1105, 41041])
def test_composites_refused(m):
    r = certify_prime(m)
    assert r["verdict"] == "NOT A PRIME" and not prime_accepted(r)


def test_large_named_primes():
    pytest.importorskip("sympy")   # Pocklington factors p-1 with sympy (installed with the qualtran extra)
    assert certify_prime(2 ** 61 - 1)["verdict"] == "PROVEN PRIME"
    assert certify_prime(2 ** 64 - 59)["verdict"] == "PROVEN PRIME"
    assert lucas_lehmer(61) and lucas_lehmer(127) and not lucas_lehmer(67)
    r = certify_prime(2 ** 127 - 1)
    assert r["verdict"] == "PROVEN PRIME" and r["lucas_lehmer"]
    r = certify_prime(2 ** 128 - 159)
    assert r["verdict"] == "PROVEN PRIME" and r["pocklington"]
    r = certify_prime(P256)
    assert r["verdict"] == "PROVEN PRIME" and r["pocklington"]


def test_named_standard_primes_without_certificate_are_labelled():
    for p in (2 ** 255 - 19, 2 ** 256 - 2 ** 32 - 977):
        r = certify_prime(p)
        assert r["verdict"].startswith("KNOWN PRIME") and "not certified" in r["verdict"] and prime_accepted(r)


def test_unnamed_large_number_is_not_accepted_on_a_probabilistic_test():
    p = 2 ** 255 + 95                                  # a prime above 2^64 without a short certificate
    r = certify_prime(p)
    assert not r["verdict"].startswith("KNOWN PRIME")
    if r["verdict"].startswith("PROBABLE"):
        assert not prime_accepted(r)
