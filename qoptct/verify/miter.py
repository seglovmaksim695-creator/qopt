# SPDX-License-Identifier: Apache-2.0
"""Trace miter: one exact scalar identity that decides correctness for ALL inputs x and ALL outcomes m.

Let B_m = (U_ref^dagger (x) I) K_m be the branch operator after composition with the inverse reference, Pi the
projector on admissible inputs (clean model: clean wires = 0), r = Tr Pi = 2^{n_in}, and C_m = <a(m)|_anc B_m Pi
(a(m): ancillas that are classical functions of the outcomes, all other ancillas 0). Completeness of the measurement
gives sum_m ||B_m Pi||_F^2 = r and Cauchy-Schwarz gives
      S := sum_m |Tr C_m|^2  <=  r * sum_m ||C_m||_F^2  <=  r^2,
with S = r^2 iff for EVERY branch m: B_m Pi = t_m Pi (x) |a(m)>, i.e. K_m = t_m U_ref (x) |a(m)> on the inputs.
S is a doubled path sum in which every variable (x, m and the path variables of both copies) is summed. It is
evaluated EXACTLY: reduction rules, then exact variable elimination in the ring Z[w] with Python integers.
Verdict: S == r^2 -> PROVED, S < r^2 -> REFUTED (S > r^2 or a non-real S -> INTERNAL).
"""
import time
from fractions import Fraction

import numpy as np

from ._parse import bits
from .pathsum import PathSum

# ------------------------------------------------------------------ Z[w], w^4 = -1; element = array [..., 4] (object ints


def zw_unit(k):
    k &= 7
    v = [0, 0, 0, 0]
    v[k & 3] = -1 if k >= 4 else 1
    return v


def zw_mul(A, B):
    """Element-wise product of [..., 4] arrays (broadcast)."""
    sh = np.broadcast_shapes(A.shape, B.shape)
    C = np.zeros(sh, dtype=np.result_type(A.dtype, B.dtype))
    for i in range(4):
        for j in range(4):
            k = i + j
            p = A[..., i] * B[..., j]
            if k < 4:
                C[..., k] = C[..., k] + p
            else:
                C[..., k - 4] = C[..., k - 4] - p
    return C


def contract_scalar(ps, max_width=24, deadline=None):
    """sum over ps.ys of w^F, exactly (no open variables allowed).
    -> (Z[w] 4-vector of Python ints, width) or (None, reason)."""
    # Every table entry is a sum over the already eliminated assignments of products of units, so each Z[w]
    # coefficient is bounded by 2^(number of summed variables); a product in zw_mul adds at most 4 such terms per
    # slot, so int64 is exact for up to 56 summed variables (object ints beyond).
    dt = np.int64 if len(ps.ys) <= 56 else object
    factors = []
    for mo, c in ps.ph.items():
        vs = bits(mo)
        if not vs:
            factors.append(([], np.array(zw_unit(c), dtype=dt)))
            continue
        t = np.empty((2,) * len(vs) + (4,), dtype=dt)
        for idx in np.ndindex(*(2,) * len(vs)):
            t[idx] = zw_unit(c if all(idx) else 0)
        factors.append((vs, t))
    scal = np.array([1, 0, 0, 0], dtype=dt)
    present = set(v for vs, _ in factors for v in vs)
    nfree = len([y for y in ps.ys if y not in present])
    fl = []
    for vs, t in factors:
        if not vs:
            scal = zw_mul(scal, t)
        else:
            fl.append((vs, t))
    todo = set(y for y in ps.ys if y in present)
    if present - set(ps.ys):
        return None, 'open variables present'
    maxw = 0
    while todo:
        if deadline is not None and time.time() > deadline:
            raise TimeoutError('trace miter exceeded the time budget')
        best = None
        for y in todo:
            sc = set()
            nf = 0
            for vs, _ in fl:
                if y in vs:
                    sc.update(vs)
                    nf += 1
            k = (len(sc), nf)
            if best is None or k < best[0]:
                best = (k, y, sc)
            if k[0] <= 2:
                break
        _, y, sc = best
        if len(sc) > max_width:
            return None, 'width %d > %d' % (len(sc), max_width)
        maxw = max(maxw, len(sc))
        scope = sorted(sc)
        pos = {v: i for i, v in enumerate(scope)}
        prod = None
        rest = []
        for vs, t in fl:
            if y not in vs:
                rest.append((vs, t))
                continue
            shape = [1] * len(scope) + [4]
            perm = sorted(range(len(vs)), key=lambda i: pos[vs[i]])
            tt = np.transpose(t, perm + [len(vs)])
            for i in range(len(vs)):
                shape[pos[vs[perm[i]]]] = 2
            tt = tt.reshape(shape)
            prod = tt if prod is None else zw_mul(prod, tt)
        prod = np.broadcast_to(prod, (2,) * len(scope) + (4,))
        r = prod.take(0, axis=pos[y]) + prod.take(1, axis=pos[y])
        ns = [v for v in scope if v != y]
        if ns:
            rest.append((ns, np.array(r, dtype=dt)))
        else:
            scal = zw_mul(scal, np.array(r, dtype=dt))
        fl = rest
        todo.discard(y)
    return [int(v) * 2 ** nfree for v in scal], maxw


def doubled(ps, info):
    """Doubled path sum whose total is S. -> (PathSum, meta) or (None, reason)."""
    n, N, xv = info['n'], info['N'], info['xv']
    kinds = ps.kind
    mmask = sum(1 << v for v in range(ps.nv) if kinds[v] == 'm')
    cons = []
    nclass = 0
    for q in range(N):
        W = ps.w[q]
        if q in xv:
            cons.append(W ^ frozenset([1 << xv[q]]))
        else:
            if all((mo & ~mmask) == 0 for mo in W) and q >= n:
                nclass += 1     # classical ancilla a(m): projected out
                continue
            cons.append(W)      # must be 0
    cons += list(ps.guards)
    D = PathSum()
    D.zpick = ps.zpick
    m1, m2 = {}, {}
    for v in range(ps.nv):
        if kinds[v] == 'm':
            if v in ps.ys:
                return None, 'outcome variable among path variables'
            u = D.var('y', 'm')
            m1[v] = m2[v] = u
        elif kinds[v] == 'x' or v in ps.ys:
            m1[v] = D.var('y', 'a')
            m2[v] = D.var('y', 'b')

    def mp(mo, mm):
        r = 0
        for v in bits(mo):
            r |= 1 << mm[v]
        return r
    for mo, c in ps.ph.items():
        D.addc(mp(mo, m1), c)
        D.addc(mp(mo, m2), -c)
    D.e2 = 2 * ps.e2
    for v in ps.frozen:      # cuts stay outside the rules in the doubled sum as well
        if v in m1:
            D.frozen |= {m1[v], m2[v]}
    for Cn in cons:
        for mm in (m1, m2):
            if not Cn:
                continue
            z = D.var('y', 'z')
            zb = 1 << z
            for mo in Cn:
                D.addc(mp(mo, mm) | zb, 4)
            D.e2 -= 2
    D.w = []
    r = 2 ** len(xv)
    return D, dict(r=r, nclass=nclass, ncons=len(cons))


def zw_times_sqrt2pow(v, e):
    """v * sqrt(2)^e as (Z[w] vector, power of two k): value = vec * 2^k (k may be negative)."""
    v = list(v)
    if e % 2:
        # sqrt2 = w - w^3
        v = [int(t) for t in zw_mul(np.array(v, dtype=object), np.array([0, 1, 0, -1], dtype=object))]
        e -= 1
    return v, e // 2


def trace_check(ps, info, max_width=24, deadline=None):
    D, meta = doubled(ps, info)
    if D is None:
        return 'UNKNOWN', meta
    nv0 = len(D.ys)
    D.deadline = deadline
    D.reduce(deadline=deadline)
    if D.zero:
        val = [0, 0, 0, 0]
        w = 0
    else:
        val, w = contract_scalar(D, max_width=max_width, deadline=deadline)
        if val is None:
            return 'UNKNOWN', dict(meta, why=w, vars=len(D.ys))
    vec, k = zw_times_sqrt2pow(val, D.e2)
    r2 = meta['r'] ** 2
    # S = vec * 2^k, vec = a + b w + c w^2 + d w^3; real iff c = 0 and b = -d; then S = (a + b sqrt2) 2^k
    meta.update(r2=r2, vars0=nv0, vars_after_rules=len(D.ys), width=w)
    if vec[2] or vec[1] != -vec[3]:
        return 'INTERNAL', dict(meta, S=vec, k=k, why='S is not real')
    a = Fraction(vec[0]) * (Fraction(2) ** k)
    b = Fraction(vec[1]) * (Fraction(2) ** k)
    meta.update(S=str(a) if not b else '%s + %s*sqrt2' % (a, b))
    if b == 0 and a == r2:
        return 'PROVED', meta
    Sf = float(a) + float(b) * 2 ** 0.5
    if Sf > r2 * (1 + 1e-12):
        return 'INTERNAL', dict(meta, why='S > r^2')
    return 'REFUTED', meta
