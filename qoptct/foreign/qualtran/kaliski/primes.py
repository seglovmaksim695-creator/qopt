"""Primality of the modulus p: a certificate where one can be made here, a clear label where not.

Certificates: exhaustive trial division (small p), deterministic Miller-Rabin (p < 2^64), Lucas-Lehmer (Mersenne numbers),
Pocklington (p-1 partly factored).  A named standard prime without a certificate is labelled KNOWN PRIME; a probabilistic
test alone is only PROBABLE PRIME and is never accepted by the verdict."""
import time
from qoptct.verify.compose import primality


def lucas_lehmer(q):
    """2^q - 1 prime iff s_{q-2} = 0 (q odd prime)."""
    m = (1 << q) - 1
    s = 4
    for _ in range(q - 2):
        s = (s * s - 2) % m
    return s == 0


def pocklington(p, limit_sec=120):
    """Pocklington: p-1 = F*R, F fully factored, F > sqrt(p), and for every prime q | F a base a with a^(p-1) = 1 and
    gcd(a^((p-1)/q) - 1, p) = 1. The prime factors q of F must themselves be proved prime (here: < 2^64 by the
    deterministic test, or recursively). -> (ok, certificate dict)."""
    from math import gcd, isqrt
    from sympy import factorint
    t = time.time()
    fac = factorint(p - 1, limit=10 ** 7)
    F, R = 1, p - 1
    proven = {}
    for q, e in sorted(fac.items()):
        v, how = primality(q)
        if v == 'prime':
            proven[q] = how
            F *= q ** e
            R //= q ** e
        elif q.bit_length() < 200 and time.time() - t < limit_sec:
            ok, cert = pocklington(q, limit_sec - (time.time() - t))
            if ok:
                proven[q] = cert
                F *= q ** e
                R //= q ** e
    if F * F <= p:
        return False, dict(F_bits=F.bit_length(), p_bits=p.bit_length(), reason='factored part too small')
    wit = {}
    for q in proven:
        for a in range(2, 200):
            if pow(a, p - 1, p) == 1 and gcd(pow(a, (p - 1) // q, p) - 1, p) == 1:
                wit[q] = a
                break
        else:
            return False, dict(reason='no witness for %d' % q)
    return True, dict(p=p, F=F, factors={str(q): str(h) if isinstance(h, str) else h for q, h in proven.items()},
                      witnesses={str(q): a for q, a in wit.items()})


NAMED = {2 ** 61 - 1: '2^61-1', 2 ** 64 - 59: '2^64-59', 2 ** 127 - 1: '2^127-1', 2 ** 128 - 159: '2^128-159',
         2 ** 255 - 19: '2^255-19', 2 ** 256 - 2 ** 32 - 977: 'secp256k1',
         2 ** 256 - 2 ** 224 + 2 ** 192 + 2 ** 96 - 1: 'P-256'}
_MERSENNE = {2 ** 61 - 1: 61, 2 ** 127 - 1: 127}
TRIAL_LIMIT = 10 ** 12


def certify_prime(p):
    """-> dict(verdict, how, ...).  verdict is 'PROVEN PRIME', 'KNOWN PRIME (...)' (named standard prime, not certified here),
    'PROBABLE PRIME (...)' or 'COMPOSITE' / 'NOT A PRIME'."""
    from math import isqrt
    if p < 2:
        return dict(p=p, verdict='NOT A PRIME', how='p < 2')
    if p % 2 == 0 and p != 2:
        return dict(p=p, verdict='NOT A PRIME', how='even')
    if p < TRIAL_LIMIT:
        r = isqrt(p)
        q = 2
        while q <= r:
            if p % q == 0:
                return dict(p=p, verdict='NOT A PRIME', how='divisible by %d' % q)
            q += 1 if q == 2 else 2
        return dict(p=p, verdict='PROVEN PRIME', how='trial division by all q<=%d (%d^2=%d>%d)' % (r, r, r * r, p))
    v, how = primality(p)
    row = dict(p=p, qoptct=v, how=how)
    if p in _MERSENNE:
        row['lucas_lehmer'] = lucas_lehmer(_MERSENNE[p])
    if v == 'composite':
        row['verdict'] = 'NOT A PRIME'
        return row
    if v != 'prime' or p in NAMED:
        t = time.time()
        try:
            ok, cert = pocklington(p)
        except Exception as e:                # noqa: BLE001
            ok, cert = False, dict(reason='%s: %s' % (type(e).__name__, e))
        row['pocklington'] = ok
        row['pocklington_cert'] = cert
        row['pocklington_sec'] = round(time.time() - t, 1)
    certified = (v == 'prime' and how.startswith('deterministic')) or row.get('lucas_lehmer') or row.get('pocklington')
    row['verdict'] = ('PROVEN PRIME' if certified else
                      'KNOWN PRIME (named standard prime; not certified here)' if (v == 'prime' and p in NAMED) else
                      'PROBABLE PRIME (probabilistic test only)')
    return row


def prime_accepted(row):
    """The verdict accepts only a certificate or a named standard prime."""
    return bool(row) and (row.get('verdict', '').startswith('PROVEN PRIME') or row.get('verdict', '').startswith('KNOWN PRIME'))
