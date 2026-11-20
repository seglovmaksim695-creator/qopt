# SPDX-License-Identifier: Apache-2.0
"""k4/ec.py — классическая арифметика эллиптических кривых y^2 = x^3 + a x + b над F_p (короткая форма Вейерштрасса),
точка на бесконечности O = None. Малые тестовые кривые (нечётный порядок группы => нет точек с y=0) и стандартные P-256, secp256k1."""
import random

P256 = dict(name='P-256',
            p=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff,
            a=0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc,
            b=0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b,
            gx=0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,
            gy=0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5,
            order=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551)
SECP256K1 = dict(name='secp256k1',
                 p=0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f, a=0, b=7,
                 gx=0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
                 gy=0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8,
                 order=0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141)


def is_prime(x):
    if x < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if x % q == 0: return x == q
    d, s = x - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        y = pow(a, d, x)
        if y in (1, x - 1): continue
        for _ in range(s - 1):
            y = y * y % x
            if y == x - 1: break
        else: return False
    return True


def ec_neg(P, p):
    return None if P is None else (P[0], (-P[1]) % p)


def ec_add(P, Q, p, a):
    if P is None: return Q
    if Q is None: return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        l = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        l = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (l * l - x1 - x2) % p
    return x3, (l * (x1 - x3) - y1) % p


def ec_mul(k, P, p, a):
    R = None
    while k:
        if k & 1: R = ec_add(R, P, p, a)
        P = ec_add(P, P, p, a); k >>= 1
    return R


def on_curve(P, cv):
    if P is None: return True
    x, y = P; p = cv['p']
    return (y * y - x ** 3 - cv['a'] * x - cv['b']) % p == 0


def curve_points(p, a, b):
    sq = {}
    for y in range(p): sq.setdefault(y * y % p, []).append(y)
    return [(x, y) for x in range(p) for y in sq.get((x ** 3 + a * x + b) % p, [])]


def small_curve(n, rng, p=None, min_points=6, prime_order=False, pfilter=None):
    """случайная кривая над простым p в (2^(n-1), 2^n) с нечётным порядком группы (#E = len(points)+1)."""
    while True:
        pp = p
        if pp is None:
            pp = rng.randrange(2 ** (n - 1) + 1, 2 ** n, 2)
            if not is_prime(pp) or (pfilter and not pfilter(pp)): continue
        for _ in range(200):
            a, b = rng.randrange(pp), rng.randrange(pp)
            if (4 * a ** 3 + 27 * b * b) % pp == 0: continue
            pts = curve_points(pp, a, b)
            N = len(pts) + 1
            if N % 2 == 0 or len(pts) < min_points: continue
            if prime_order and not is_prime(N): continue
            return dict(name=f'toy{n}', p=pp, a=a, b=b, points=pts, order=N)
