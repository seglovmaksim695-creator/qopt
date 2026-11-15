# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""Polynomial-basis multiplication in F2[x]/(f): recognition of the tensor (with arbitrary wire permutations),
change of basis between moduli, and the Karatsuba fallback decomposition.

A decomposition is a list of (alpha, beta, gamma) int masks: c = sum_r (alpha_r . a)(beta_r . b) gamma_r.
"""


def pmod(a, f):
    df = f.bit_length() - 1
    while a and a.bit_length() - 1 >= df:
        a ^= f << (a.bit_length() - 1 - df)
    return a


def pmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def mulmod(a, b, f):
    return pmod(pmul(a, b), f)


def powmod(a, e, f):
    r = 1
    while e:
        if e & 1:
            r = mulmod(r, a, f)
        a = mulmod(a, a, f)
        e >>= 1
    return r


def pgcd(a, b):
    while b:
        a, b = b, pmod(a, b)
    return a


def is_irreducible(f):
    """Ben-Or: f of degree d is irreducible iff gcd(x^(2^i) - x mod f, f) = 1 for i = 1..d/2."""
    d = f.bit_length() - 1
    if d < 1:
        return False
    if d == 1:
        return True
    xp = 2
    for _ in range(d // 2):
        xp = mulmod(xp, xp, f)
        if pgcd(f, xp ^ 2) != 1:
            return False
    return True


# ---- polynomials with coefficients in K = F2[x]/(fk) (lists, low degree first) -------------------------------

def _ktrim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def _kinv(a, fk, m):
    return powmod(a, (1 << m) - 2, fk)


def _kmod(a, g, fk, m):
    a = list(a)
    dg = len(g) - 1
    inv = _kinv(g[-1], fk, m)
    while len(_ktrim(a)) - 1 >= dg:
        c = mulmod(a[-1], inv, fk)
        s = len(a) - 1 - dg
        for i, gi in enumerate(g):
            if gi:
                a[s + i] ^= mulmod(c, gi, fk)
    return a


def _kmul(a, b, fk):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[i + j] ^= mulmod(x, y, fk)
    return r


def _kgcd(a, b, fk, m):
    a, b = _ktrim(list(a)), _ktrim(list(b))
    while b:
        a, b = b, _ktrim(_kmod(a, b, fk, m))
    return a


def find_root(g, fk, m, seed=1):
    """a root in K = F2[x]/(fk) (degree m) of g (F2 polynomial as int, irreducible of degree m)."""
    import random
    rng = random.Random(seed)
    G = [((g >> i) & 1) for i in range(g.bit_length())]
    while len(G) - 1 > 1:
        for _ in range(200):
            delta = rng.randrange(1, 1 << m)
            # T(y) = sum_i (delta y)^(2^i) mod G
            cur = _kmod([0, delta], G, fk, m)
            acc = list(cur)
            for _ in range(m - 1):
                cur = _kmod(_kmul(cur, cur, fk), G, fk, m)
                acc = [x ^ y for x, y in zip(acc + [0] * (len(cur) - len(acc)), cur + [0] * (len(acc) - len(cur)))]
            d = _kgcd(G, acc, fk, m)
            if 0 < len(d) - 1 < len(G) - 1:
                if len(d) - 1 <= (len(G) - 1) // 2:
                    G = d
                else:
                    G = _ktrim(_kdiv(G, d, fk, m))
                break
        else:
            raise ValueError("root finding failed")
    # linear: G1 y + G0 = 0
    return mulmod(G[0], _kinv(G[1], fk, m), fk)


def _kdiv(a, b, fk, m):
    a = list(a)
    db = len(b) - 1
    inv = _kinv(b[-1], fk, m)
    q = [0] * max(1, len(a) - db)
    while len(_ktrim(a)) - 1 >= db:
        c = mulmod(a[-1], inv, fk)
        s = len(a) - 1 - db
        q[s] = c
        for i, bi in enumerate(b):
            if bi:
                a[s + i] ^= mulmod(c, bi, fk)
    return q


def field_tensor(m, f):
    """{(i, j): mask of x^(i+j) mod f} (zero entries dropped)."""
    out = {}
    for i in range(m):
        for j in range(m):
            v = pmod(1 << (i + j), f)
            if v:
                out[(i, j)] = v
    return out


def check_dec(dec, m, f):
    from .anf import tensor_of
    return tensor_of(dec) == field_tensor(m, f)


def identify(A, B, T, terms):
    """Recognise c ^= a*b mod f in a permuted polynomial basis.

    Returns (m, f, pa, pb, pc) with pa[k] = index in A of the coefficient of x^k (same for B, T), or None."""
    m = len(A)
    if len(B) != m or len(T) != m or m < 2:
        return None
    single = lambda v: v != 0 and v & (v - 1) == 0
    bit = lambda v: v.bit_length() - 1
    # the element 1 in A and in B: all products are single targets, bijectively
    ones_a = [i for i in range(m) if all(single(terms.get((i, j), 0)) for j in range(m))
              and len({terms[(i, j)] for j in range(m)}) == m]
    ones_b = [j for j in range(m) if all(single(terms.get((i, j), 0)) for i in range(m))
              and len({terms[(i, j)] for i in range(m)}) == m]
    for a0 in ones_a:
        for b0 in ones_b:
            one = terms[(a0, b0)]
            if not single(one):
                continue
            # label maps: T index of a_i * 1, 1 * b_j
            labA = {bit(terms[(i, b0)]): i for i in range(m)}
            labB = {bit(terms[(a0, j)]): j for j in range(m)}
            e0 = bit(one)
            for a1 in range(m):
                if a1 == a0:
                    continue
                chain = [e0]
                ok = True
                for _ in range(m - 1):
                    v = terms.get((a1, labB[chain[-1]]), 0)
                    if not single(v) or bit(v) in chain:
                        ok = False
                        break
                    chain.append(bit(v))
                if not ok:
                    continue
                xm = terms.get((a1, labB[chain[-1]]), 0)
                pos = {t: k for k, t in enumerate(chain)}
                f = (1 << m) | sum(1 << pos[t] for t in range(m) if xm >> t & 1)
                pc = chain
                pa = [labA[t] for t in chain]
                pb = [labB[t] for t in chain]
                ft = field_tensor(m, f)
                good = True
                for i in range(m):
                    for j in range(m):
                        v = ft.get((i, j), 0)
                        want = sum(1 << pc[k] for k in range(m) if v >> k & 1)
                        if terms.get((pa[i], pb[j]), 0) != want:
                            good = False
                            break
                    if not good:
                        break
                if good:
                    return m, f, pa, pb, pc
    return None


def _apply(cols, v):
    r = 0
    k = 0
    while v:
        if v & 1:
            r ^= cols[k]
        v >>= 1
        k += 1
    return r


def _inverse_cols(cols, m):
    """inverse of the linear map given by images of unit vectors."""
    from .gf2lin import express
    inv = []
    for k in range(m):
        S = express(1 << k, cols)
        if S is None:
            raise ValueError("singular")
        inv.append(sum(1 << s for s in S))
    return inv


def change_basis(dec, m, f_from, f_to):
    """decomposition for F2[x]/(f_from) -> same field in the basis of F2[x]/(f_to) (both irreducible)."""
    if f_from == f_to:
        return list(dec)
    def value(cand):
        v = 0
        pw = 1
        for k in range(m + 1):
            if f_from >> k & 1:
                v ^= pw
            pw = mulmod(pw, cand, f_to)
        return v
    xi = None
    if m <= 10:
        xi = next((c for c in range(2, 1 << m) if value(c) == 0), None)
    else:
        r = find_root(f_from, f_to, m)
        conj = [r]
        for _ in range(m - 1):
            conj.append(mulmod(conj[-1], conj[-1], f_to))
        xi = min(conj)  # canonical choice (= the smallest root, as an exhaustive search would give)
    if xi is None or value(xi) != 0:
        raise ValueError("no root: moduli define different rings")
    psi = []  # psi[j] = image of x^j (old basis) in the new basis
    pw = 1
    for _ in range(m):
        psi.append(pw)
        pw = mulmod(pw, xi, f_to)
    psi_inv = _inverse_cols(psi, m)  # new unit vector -> old coordinates
    out = []
    for a, b, c in dec:
        a2 = sum(1 << i for i in range(m) if bin(a & psi_inv[i]).count("1") & 1)
        b2 = sum(1 << i for i in range(m) if bin(b & psi_inv[i]).count("1") & 1)
        out.append((a2, b2, _apply(psi, c)))
    return out


def isotope(dec, m, f, k, lam, mu):
    """a*b = sigma^-k( nu * D(sigma^k(lam a), sigma^k(mu b)) ), nu = (lam mu)^-1, sigma = Frobenius."""
    def frob(v, times):
        for _ in range(times % m):
            v = mulmod(v, v, f)
        return v
    lm = mulmod(lam, mu, f)
    nu = powmod(lm, (1 << m) - 2, f)  # inverse in the field
    if mulmod(lm, nu, f) != 1:
        raise ValueError("not invertible (modulus not irreducible?)")
    Pa = [frob(mulmod(lam, 1 << i, f), k) for i in range(m)]  # images of unit vectors
    Pb = [frob(mulmod(mu, 1 << i, f), k) for i in range(m)]
    out = []
    for a, b, c in dec:
        a2 = sum(1 << i for i in range(m) if bin(a & Pa[i]).count("1") & 1)
        b2 = sum(1 << i for i in range(m) if bin(b & Pb[i]).count("1") & 1)
        c2 = mulmod(nu, frob(c, m - (k % m)), f)
        out.append((a2, b2, c2))
    return out


def karatsuba_poly(n):
    """bilinear algorithm for the product of two polynomials with n coefficients: gamma over 2n-1 coefficients."""
    if n == 1:
        return [(1, 1, 1)]
    if n == 2:
        return [(1, 1, 0b011), (2, 2, 0b110), (3, 3, 0b010)]
    h = n // 2
    lo = (1 << h) - 1
    out = []
    for a, b, c in karatsuba_poly(h):
        out.append((a, b, c ^ (c << h)))
    for a, b, c in karatsuba_poly(n - h):
        out.append((a << h, b << h, (c << (2 * h)) ^ (c << h)))
    for a, b, c in karatsuba_poly(n - h):
        out.append(((a & lo) | (a << h), (b & lo) | (b << h), c << h))
    # cancel duplicate terms
    acc = {}
    for a, b, c in out:
        acc[(a, b)] = acc.get((a, b), 0) ^ c
    return [(a, b, c) for (a, b), c in acc.items() if c]


def karatsuba_field(m, f):
    red = [pmod(1 << k, f) for k in range(2 * m - 1)]
    out = []
    for a, b, c in karatsuba_poly(m):
        g = _apply(red, c)
        if g:
            out.append((a, b, g))
    return out
