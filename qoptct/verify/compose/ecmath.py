# SPDX-License-Identifier: Apache-2.0
"""compose.ecmath -- the machine-checkable parts of Theorem T (PROOF.md section 10): the word-level program W_ec of
the point addition computes ec_add(P, Q) for every point P of the curve and for O, and satisfies the assumptions
D-input consequences, G1..G4 used as preconditions by compose.ecadd.

1. check_identities(): three polynomial identities over Z (hence over every field), by exact expansion with integer
   coefficients (own 40-line polynomial arithmetic; no dependency):
     I1  chord:    the third point of the chord through P != +-Q and Q lies on the curve
     I2  tangent:  4 yq^2 (x^3 + a x + b) - (M (x - xq) + 2 yq^2)^2 = (x - xq)^2 (4 yq^2 x - M^2 + 8 xq yq^2),
                   M = 3 xq^2 + a, b = yq^2 - xq^3 - a xq   (the tangent at Q meets the curve in Q, Q, -2Q)
     I3  if the chord through P and Q gives x' = xq, then P lies on the tangent at Q
   "modulo the curve" means: after the substitution y^2 -> x^3 + a x + b the polynomial is identically 0.
2. theorem_T_violations(p, a, b, Q, D): an independent field-level model of W_ec (not the contract chain, not the
   circuit) executed on ALL points of a small curve; exhaustive_T(...) runs it over all curves / all Q / all D of
   small primes.  This is the "explicit enumeration of small fields" part.
The general statement is the written proof of PROOF.md section 10, which uses only I1..I3 and the field axioms."""
from .ecadd import ec_add, ec_neg, on_curve, curve_points, _two_torsion

__all__ = ['check_identities', 'theorem_T_violations', 'exhaustive_T', 'Pol']


class Pol:
    """multivariate polynomial with integer coefficients: {exponent tuple: coefficient}."""
    VARS = ('a', 'x', 'y', 'xq', 'yq')

    def __init__(self, d=None):
        self.d = {k: v for k, v in (d or {}).items() if v}

    @staticmethod
    def var(name):
        e = [0] * len(Pol.VARS)
        e[Pol.VARS.index(name)] = 1
        return Pol({tuple(e): 1})

    @staticmethod
    def const(c):
        return Pol({(0,) * len(Pol.VARS): c})

    def _c(self, o):
        return o if isinstance(o, Pol) else Pol.const(o)

    def __add__(self, o):
        o = self._c(o)
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d.get(k, 0) + v
        return Pol(d)
    __radd__ = __add__

    def __neg__(self):
        return Pol({k: -v for k, v in self.d.items()})

    def __sub__(self, o):
        return self + (-self._c(o))

    def __rsub__(self, o):
        return self._c(o) - self

    def __mul__(self, o):
        o = self._c(o)
        d = {}
        for k1, v1 in self.d.items():
            for k2, v2 in o.d.items():
                k = tuple(i + j for i, j in zip(k1, k2))
                d[k] = d.get(k, 0) + v1 * v2
        return Pol(d)
    __rmul__ = __mul__

    def __pow__(self, n):
        r = Pol.const(1)
        for _ in range(n):
            r = r * self
        return r

    def is_zero(self):
        return not self.d

    def reduce_y(self, g):
        """substitute y^2 -> g (a polynomial without y) until the degree in y is <= 1."""
        iy = Pol.VARS.index('y')
        cur = self
        while True:
            hi = {k: v for k, v in cur.d.items() if k[iy] >= 2}
            if not hi:
                return cur
            lo = Pol({k: v for k, v in cur.d.items() if k[iy] < 2})
            hp = Pol({k[:iy] + (k[iy] - 2,) + k[iy + 1:]: v for k, v in hi.items()})
            cur = lo + hp * g


def check_identities():
    """-> dict name -> True; raises AssertionError if an identity fails."""
    a, x, y, xq, yq = (Pol.var(v) for v in Pol.VARS)
    b = yq ** 2 - xq ** 3 - a * xq                       # Q on the curve
    g = x ** 3 + a * x + b                               # P on the curve: y^2 = g
    dx, dy = x - xq, y - yq
    out = {}
    # I1: x' = Xn / dx^2, y' = Yn / dx^3 with lambda = dy / dx;  Yn^2 = Xn^3 + a Xn dx^4 + b dx^6  (mod the curve)
    Xn = dy ** 2 - (x + xq) * dx ** 2
    Yn = dy * (xq * dx ** 2 - Xn) - yq * dx ** 3
    out['I1 chord point on the curve'] = (Yn ** 2 - Xn ** 3 - a * Xn * dx ** 4 - b * dx ** 6).reduce_y(g).is_zero()
    # I2: the tangent at Q
    M = 3 * xq ** 2 + a
    out['I2 tangent at Q'] = (4 * yq ** 2 * (x ** 3 + a * x + b) - (M * dx + 2 * yq ** 2) ** 2
                              - dx ** 2 * (4 * yq ** 2 * x - M ** 2 + 8 * xq * yq ** 2)).is_zero()
    # I3: (dy^2 - (x + 2 xq) dx^2) = (M x + a xq - xq^3 + 2 b) - 2 yq y   (mod the curve):
    #     x' = xq  <=>  dy^2 = (x + 2 xq) dx^2  =>  2 yq y = M x + a xq - xq^3 + 2 b = M (x - xq) + 2 yq^2
    out['I3 x\' = xq puts P on the tangent'] = (dy ** 2 - (x + 2 * xq) * dx ** 2
                                                - (M * x + a * xq - xq ** 3 + 2 * b) + 2 * yq * y).reduce_y(g).is_zero()
    out['I3b tangent line form'] = (M * x + a * xq - xq ** 3 + 2 * b - (M * dx + 2 * yq ** 2)).is_zero()
    for k, v in out.items():
        assert v, 'polynomial identity failed: %s' % k
    return out


def theorem_T_violations(p, a, b, Q, D, points=None):
    """field-level model of W_ec on every point (default: all points of the curve and O). -> list of violations
    (empty: Theorem T holds for this curve, Q, D on these points). Independent of compose.ecadd.ECChain."""
    inv = lambda v: pow(v, -1, p)
    Q2 = ec_add(Q, Q, p, a)
    mQ, m2Q = ec_neg(Q, p), ec_neg(Q2, p)
    E = [None, Q, mQ, m2Q]
    img = [Q, Q2, None, mQ]
    DQ = ec_add(D, Q, p, a)
    code = lambda P: (0, 0) if P is None else P
    xq, yq = Q
    if points is None:
        points = [None] + curve_points(p, a, b)
    viol = []
    for P in points:
        x, y = code(P)
        # the flags as the CIRCUIT computes them: comparisons of x and the parity of y, guards by f_O / f_mQ
        fO = int(y == 0)
        g = int(x == xq and not (xq == 0 and fO))
        fQ = g & int(y & 1 == yq & 1)
        fmQ = g ^ fQ
        g = int(x == Q2[0] and not (Q2[0] == 0 and fO))
        fm2Q = g & int(y & 1 == m2Q[1] & 1)
        F = [fO, fQ, fmQ, fm2Q]
        if F != [int(P == e) for e in E]:
            viol.append((P, 'flags', F))
            continue
        if any(F):
            x, y = D
        # generic addition
        dx, dy = (x - xq) % p, (y - yq) % p
        if dx == 0:
            viol.append((P, 'G1 first denominator', None))
            continue
        lam = (-dy * inv(dx)) % p                        # the register holds -lambda
        y1 = (dy + lam * dx) % p                         # 0
        x1 = (lam * lam - dx - 3 * xq) % p               # x' - xq
        y2 = (y1 + lam * x1) % p                         # lambda (xq - x') = y' + yq
        den = (-x1) % p                                  # xq - x'
        if den == 0:
            viol.append((P, 'G1 second denominator', None))
            continue
        if (lam - (-y2 * inv(den))) % p:
            viol.append((P, 'G2 lambda', None))
        xr, yr = (x1 + xq) % p, (y2 - yq) % p
        if y1 != 0:
            viol.append((P, 'y cancel', None))
        if any(F) and (xr, yr) != DQ:
            viol.append((P, 'G3', (xr, yr)))
        if not any(F):
            if not on_curve((xr, yr), p, a, b) or yr == 0 or (xr, yr) in (Q, mQ, Q2):
                viol.append((P, 'G4', (xr, yr)))
        # fix-up
        for j in range(4):
            if F[j]:
                xr, yr = code(img[j])
        # flag erasure as the circuit computes it
        g = int(xr == xq and not (xq == 0 and fmQ))
        if fm2Q ^ g ^ fO:
            viol.append((P, 'erase fm2Q', None))
        if fO != (g & int(yr & 1 == yq & 1)):
            viol.append((P, 'erase fO', None))
        g = int(xr == Q2[0] and not (Q2[0] == 0 and fmQ))
        if fQ != (g & int(yr & 1 == Q2[1] & 1)):
            viol.append((P, 'erase fQ', None))
        if fmQ != int(yr == 0):
            viol.append((P, 'erase fmQ', None))
        if (xr, yr) != code(ec_add(P, Q, p, a)):
            viol.append((P, 'result', (xr, yr)))
    return viol


def exhaustive_T(primes, all_D=False, max_curves=None, progress=None):
    """Theorem T by enumeration: every curve y^2 = x^3 + a x + b over F_p (p in primes) without a point of order 2,
    every Q with 2Q != +-Q, D = 3Q when it is safe (all_D: EVERY safe D), every point P and O.
    -> dict(curves, cases (curve, Q, D), points, violations)."""
    ncv = ncase = npts = 0
    viol = []
    for p in primes:
        done = 0
        for a in range(p):
            for b in range(p):
                if (4 * a ** 3 + 27 * b * b) % p == 0 or _two_torsion(p, a, b):
                    continue
                pts = curve_points(p, a, b)
                if not pts:
                    continue
                if max_curves is not None and done >= max_curves:
                    break
                done += 1
                ncv += 1
                allp = [None] + pts
                for Q in pts:
                    Q2 = ec_add(Q, Q, p, a)
                    if Q2 is None or Q2[0] == Q[0]:
                        continue
                    m2Q = ec_neg(Q2, p)

                    def safe(Dc):
                        if Dc is None or Dc[0] == Q[0] or Dc == m2Q:
                            return False
                        DQ = ec_add(Dc, Q, p, a)
                        return DQ is not None and DQ[0] != Q[0]
                    Ds = [Dc for Dc in pts if safe(Dc)] if all_D else [Dc for Dc in [ec_add(Q2, Q, p, a)] if safe(Dc)]
                    for Dc in Ds:
                        ncase += 1
                        npts += len(allp)
                        v = theorem_T_violations(p, a, b, Q, Dc, allp)
                        if v and len(viol) < 20:
                            viol.append((p, a, b, Q, Dc, v[:3]))
        if progress:
            progress('p = %d: %d curves so far, %d cases, %d violations' % (p, ncv, ncase, len(viol)))
    return dict(curves=ncv, cases=ncase, points=npts, violations=viol)
