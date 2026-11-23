"""Domain of the proofs: the conditions on (n, p) under which a proof of a Qualtran block is meaningful at all.

A proof client must say NOT PROVED (with the reason) when these conditions fail, instead of proving something about a
modulus the contracts do not cover.  The conditions:

``montgomery`` / ``pair-montgomery``   p odd, p >= 3, 2p <= 2^n (the multiplier's boxes keep the values below 2p) and p < 2^n - 1.
                                       Primality is NOT needed (contract C needs only an odd p).
``kaliski``                            p an odd PRIME (certificate or named standard prime, ``kaliski.primes``), 3 <= p < 2^n - 1.
``pair-kaliski``                       odd p, 3 <= p < 2^n - 1.  The mirror identity does not use primality; the forward contract of
                                       the inverse does (see ``kaliski``), and is a separate result.
``ecadd`` / ``ecwindowaddr``           all of the above together: p an odd prime with 2p <= 2^n and p < 2^n - 1.  ECAdd divides modulo p
                                       (inverse, lemma "MM(u,v) = 0 iff u = 0 or v = 0") and multiplies with the Montgomery box.

``check_domain(kind, n, p)`` -> ``dict(ok, problems, prime)``; it never raises for a bad p.
"""
from __future__ import annotations

__all__ = ['check_domain', 'KINDS']

KINDS = ('montgomery', 'pair-montgomery', 'kaliski', 'pair-kaliski', 'ecadd', 'ecwindowaddr')
_PRIME_CACHE = {}


def _prime_row(p):
    if p not in _PRIME_CACHE:
        from .kaliski.primes import certify_prime
        _PRIME_CACHE[p] = certify_prime(p)
    return _PRIME_CACHE[p]


def check_domain(kind, n, p):
    """Check the conditions of ``kind`` for the bit size ``n`` and the modulus ``p``.  -> dict(ok, problems, prime, kind, n, p)."""
    if kind not in KINDS:
        raise ValueError(f'unknown kind {kind!r}')
    problems = []
    try:
        n, p = int(n), int(p)
    except (TypeError, ValueError):
        return dict(kind=kind, n=n, p=str(p), ok=False, problems=['n and p must be integers'], prime=None)
    if n < 2:
        problems.append(f'n = {n} is too small (n >= 2)')
    if p < 3:
        problems.append(f'p = {p} < 3')
    elif p % 2 == 0:
        problems.append('p is even (p must be odd)')
    if n >= 2 and p >= 3:
        if p >= (1 << n) - 1:
            problems.append(f'p = {p} is not below 2^n - 1 = {(1 << n) - 1}')
        if kind in ('montgomery', 'pair-montgomery', 'ecadd', 'ecwindowaddr') and 2 * p > (1 << n):
            problems.append(f'2p > 2^n (the multiplier needs 2p <= 2^{n})')
    prime = None
    if kind in ('kaliski', 'ecadd', 'ecwindowaddr') and p >= 3 and p % 2 == 1:
        from .kaliski.primes import prime_accepted
        prime = _prime_row(p)
        if not prime_accepted(prime):
            problems.append(f"p is not accepted as a prime: {prime.get('verdict')} ({prime.get('how')})")
    return dict(kind=kind, n=n, p=str(p), ok=not problems, problems=problems,
                prime=None if prime is None else dict(verdict=prime.get('verdict'), how=prime.get('how')))
