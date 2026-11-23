"""Своя арифметика кривой y^2 = x^3 + A x + B над F_p (точка O = None), кодирование Монтгомери, сценарии входов."""
import random


def inv(a, p):
    return pow(a % p, p - 2, p)


def sqrt_mod(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    q, s = p - 1, 0
    while q % 2 == 0:
        q, s = q // 2, s + 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2, i = t2 * t2 % p, i + 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def add(P, Q, A, p):
    if P is None:
        return Q
    if Q is None:
        return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        lam = (3 * x1 * x1 + A) * inv(2 * y1, p) % p
    else:
        lam = (y2 - y1) * inv(x2 - x1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def on_curve(P, A, B, p):
    return P is None or (P[1] * P[1] - (P[0] ** 3 + A * P[0] + B)) % p == 0


def mont(v, n, p):
    return v * pow(2, n, p) % p


def enc(P, n, p):
    return (0, 0) if P is None else (mont(P[0], n, p), mont(P[1], n, p))


def lam_r(P, A, n, p):
    """Montgomery-форма наклона касательной в P; для O и y=0 — 0 (как ECWindowAddR)."""
    if P is None or P[1] == 0:
        return 0
    return mont((3 * P[0] * P[0] + A) * inv(2 * P[1], p) % p, n, p)


def curve_through(P, A, p):
    return (P[1] * P[1] - P[0] ** 3 - A * P[0]) % p


def good_curve(A, B, p):
    return B % p != 0 and (4 * A ** 3 + 27 * B * B) % p != 0


def random_point(A, B, p, rnd):
    while True:
        x = rnd.randrange(p)
        y = sqrt_mod(x ** 3 + A * x + B, p)
        if y is not None:
            return (x, y if rnd.random() < 0.5 else (-y) % p)


def scenarios(p, rnd, count):
    """Список (имя, A, B, P, Q) — точки на общей кривой; O = None."""
    out = []
    kinds = ['generic', 'double', 'neg', 'Pzero', 'Qzero', 'both_zero', 'minus2P', 'G1', 'G2', 'G3like']
    while len(out) < count:
        k = kinds[len(out) % len(kinds)]
        a, b = rnd.randrange(1, p), rnd.randrange(1, p)
        if k == 'G1':                       # противоположные y при разных x
            x2 = rnd.randrange(p)
            if x2 == a:
                continue
            A = (-(a * a + a * x2 + x2 * x2)) % p
            B = curve_through((a, b), A, p)
            P, Q = (a, b), (x2, (-b) % p)
        elif k == 'G2':                     # удвоение с нулевым наклоном
            A = (-3 * a * a) % p
            B = curve_through((a, b), A, p)
            P = Q = (a, b)
        else:
            A = rnd.randrange(p)
            B = curve_through((a, b), A, p)
            P = (a, b)
            if not good_curve(A, B, p):
                continue
            Q = {'generic': lambda: random_point(A, B, p, rnd), 'double': lambda: P,
                 'neg': lambda: (a, (-b) % p), 'Pzero': lambda: random_point(A, B, p, rnd),
                 'Qzero': lambda: None, 'both_zero': lambda: None,
                 'minus2P': lambda: (lambda R: None if R is None else (R[0], (-R[1]) % p))(add(P, P, A, p)),
                 'G3like': lambda: random_point(A, B, p, rnd)}[k]()
            if k in ('Pzero', 'both_zero', 'G3like'):
                P = None
        if not good_curve(A, B, p) or not on_curve(P, A, B, p) or not on_curve(Q, A, B, p):
            continue
        out.append((k, A, B, P, Q))
    return out
