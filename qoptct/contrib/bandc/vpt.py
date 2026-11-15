"""Scalable sampled check of the QASM itself for CLASSICAL reversible references (path-sum engine of c10).

For a basis input x and measurement outcomes m the expected output is d* = (classical reference on the data wires,
ancillas 0). The check asserts |A(d*)| = 1 (the normalized branch is concentrated in one point, so all other
amplitudes are 0) and that the phase A(d*) / A_ref(d*) = c_m is the same for all sampled x (=> branch = c_m U).
A zero/inconsistent branch is an error (for a correct protocol all X-measurement outcomes are equally likely).
Unlike enumerating the image of the path sum (2^rank points), this works with hundreds of measurements.

Requirements: every measurement is an X-measurement (``h q[a]; c[k] = measure q[a];``); data wires stay alive
(use the variant with ``reset`` after measuring a data wire, see compiler.to_qasm(ndata=...)).
"""
import math
import random
import time

import numpy as np

from qoptct.contrib.c10 import pathsum as P
from qoptct.qopt_core import load_qc


class Inconclusive(Exception):
    """the sampled check could not be completed (too wide / too slow / unsupported) -- never a pass"""


def _bits(m):
    r = []; i = 0
    while m:
        if m & 1: r.append(i)
        m >>= 1; i += 1
    return r


def classical_target(hl, n, x, with_sign=False):
    """classical reference output: x/cx/ccx directly, z/cz/ccz with one H-frame = X-type on it.
    with_sign: also the sign (-1)^s of the reference amplitude (exact, no path sum)."""
    b = [x.get(q, 0) for q in range(n)]; fr = set(); sg = 0
    for g, q in hl:
        if g == 'h': fr ^= {q[0]}; continue
        F = [u for u in q if u in fr]
        if g in ('x', 'cx', 'ccx'):
            c, t = q[:-1], q[-1]
            if any(u in fr for u in c): raise Inconclusive('control in H-frame')
            if t in fr:                   # X-type on a frame = Z-type: sign (controls * t)
                if all(b[u] for u in c) and b[t]: sg ^= 1
                continue
            if all(b[u] for u in c): b[t] ^= 1
        elif g in ('z', 'cz', 'ccz'):
            if not F:
                if all(b[u] for u in q): sg ^= 1
                continue
            if len(F) != 1: raise Inconclusive('two H-frames')
            t = F[0]
            if all(b[u] for u in q if u != t): b[t] ^= 1
        else: raise Inconclusive('gate ' + g)
    if fr: raise Inconclusive('open H-frame')
    return (b, sg) if with_sign else b


OM = [complex(math.cos(math.pi * k / 4), math.sin(math.pi * k / 4)) for k in range(8)]
_XOR3 = np.zeros((2, 2, 2), dtype=complex)
for _a in (0, 1):
    for _b in (0, 1): _XOR3[_a, _b, _a ^ _b] = 1.0


def amp_sum(ps, targets, maxscope=26, lg_extra=0.0, deadline=None):
    """Amplitude <targets| branch> of a path sum. Wide parities become chains of XOR tensors with auxiliary variables
    (width 3, shared prefixes); constants are propagated and equal variables merged; the remaining variables are
    eliminated greedily (min scope). The scale is tracked in log2 (2^(-nH/2) with thousands of H under/overflows).
    maxscope bounds the largest intermediate tensor (2^25 complex = 512 MB); wider -> Inconclusive."""
    if not ps.consistent: return 0j
    LG = [float(lg_extra)]

    def _nz(x):
        a = abs(x)
        if a > 0 and (a > 1e30 or a < 1e-30): LG[0] += math.log2(a); return x / a
        return x
    scal = OM[ps.gph % 8]; factors = []; nxt = [ps.k]; node = {}

    def pv(m):
        if m in node: return node[m]
        vs = _bits(m)
        if len(vs) == 1: node[m] = vs[0]; return vs[0]
        pre = m ^ (1 << vs[-1]); a = pv(pre); s = nxt[0]; nxt[0] += 1
        factors.append(([a, vs[-1], s], _XOR3)); node[m] = s
        return s
    for m, c in ps.lin.items():
        factors.append(([pv(m)], np.array([1.0 + 0j, OM[c % 8]])))
    for (ma, mb) in ps.quad:
        a, b = pv(ma), pv(mb)
        factors.append(([a, b], np.array([[1, 1], [1, -1]], dtype=complex)))
    for m, bit in ps.cons:
        t = np.zeros(2, dtype=complex); t[bit] = 1; factors.append(([pv(m)], t))
    for q, bit in targets.items():
        m, c = ps.w[q]
        if m == 0:
            if c != bit: return 0j
            continue
        t = np.zeros(2, dtype=complex); t[bit ^ c] = 1; factors.append(([pv(m)], t))
    var_f = {}; fl = {}
    for i, (vs, t) in enumerate(factors):
        if len(set(vs)) < len(vs):   # repeated variable in one factor: take the diagonal
            u = sorted(set(vs)); lab = {v: j for j, v in enumerate(u)}
            t = np.einsum(t, [lab[v] for v in vs], list(range(len(u)))); vs = u
        fl[i] = (list(vs), t)
        for v in vs: var_f.setdefault(v, set()).add(i)
    used = set(var_f)
    LG[0] += len([v for v in range(ps.k) if v not in used])
    alive = set(var_f); nf = len(factors)
    cand = [i for i, (vs, t) in fl.items() if len(vs) <= 2]

    def fix_all():
        """constant propagation: a one-variable factor [a, 0] / [0, a] fixes the variable; a two-variable
        permutation factor merges two variables (u = v xor c)."""
        nonlocal scal
        while cand:
            i = cand.pop()
            if i not in fl: continue
            vs, t = fl[i]
            if len(vs) == 0:
                scal = _nz(scal * complex(t)); del fl[i]; continue
            if len(vs) == 2:
                a00, a01, a10, a11 = complex(t[0, 0]), complex(t[0, 1]), complex(t[1, 0]), complex(t[1, 1])
                z = lambda x: abs(x) < 1e-12
                if z(a01) and z(a10) and not z(a00) and abs(a00 - a11) < 1e-12: c = 0; amp = a00
                elif z(a00) and z(a11) and not z(a01) and abs(a01 - a10) < 1e-12: c = 1; amp = a01
                else: continue
                u, v = vs
                if u == v: continue
                scal = _nz(scal * amp); del fl[i]; var_f[u].discard(i); var_f[v].discard(i)
                for j in list(var_f[u]):
                    vs2, t2 = fl[j]; k = vs2.index(u)
                    if c: t2 = np.flip(t2, axis=k)
                    if v in vs2:
                        u_ = sorted(set(vs2) - {u}); lab = {x: n_ for n_, x in enumerate(u_)}
                        ins = [lab[x] if x != u else lab[v] for x in vs2]
                        t2 = np.einsum(t2, ins, list(range(len(u_)))); vs2 = u_
                    else:
                        vs2 = list(vs2); vs2[k] = v; var_f.setdefault(v, set()).add(j)
                    fl[j] = (list(vs2), t2)
                    if len(vs2) <= 2: cand.append(j)
                del var_f[u]; alive.discard(u)
                if not var_f.get(v):   # v lost all its factors: the sum over it gives 2
                    LG[0] += 1; var_f.pop(v, None); alive.discard(v)
                continue
            if len(vs) != 1: continue
            v = vs[0]
            a0, a1 = complex(t[0]), complex(t[1])
            if abs(a0) < 1e-12 and abs(a1) < 1e-12:
                scal = 0j; return
            if abs(a0) > 1e-12 and abs(a1) > 1e-12: continue
            val = 0 if abs(a1) < 1e-12 else 1
            scal = _nz(scal * (a0 if val == 0 else a1))
            del fl[i]; var_f[v].discard(i)
            for j in list(var_f[v]):
                vs2, t2 = fl[j]; k = vs2.index(v)
                t2 = np.take(t2, val, axis=k); vs2 = vs2[:k] + vs2[k + 1:]
                fl[j] = (vs2, t2)
                if len(vs2) <= 2: cand.append(j)
            del var_f[v]; alive.discard(v)
    fix_all()
    if scal == 0: return 0j
    alive = {u for u in alive if var_f.get(u)}
    steps = 0
    while alive:
        steps += 1
        if deadline is not None and steps % 32 == 0 and time.time() > deadline: raise Inconclusive('timeout')
        best = None
        for v in alive:
            sc = set()
            for i in var_f[v]: sc |= set(fl[i][0])
            key = (len(sc), len(var_f[v]))
            if best is None or key < best[0]: best = (key, v, sc)
        _, v, sc = best
        if len(sc) > maxscope: raise Inconclusive('elimination width %d' % len(sc))
        fi = sorted(var_f[v]); scope = sorted(sc); lab = {u: j for j, u in enumerate(scope)}
        args = []
        for i in fi: args += [fl[i][1], [lab[u] for u in fl[i][0]]]
        newscope = [u for u in scope if u != v]
        r = np.einsum(*args, [lab[u] for u in newscope], optimize=False)
        mx = float(np.max(np.abs(r))) if r.size else 0.0
        if mx > 0: r /= mx; LG[0] += math.log2(mx)   # in place (large tensors)
        for i in fi:
            for u in fl[i][0]: var_f[u].discard(i)
            del fl[i]
        alive.discard(v); del var_f[v]
        if newscope:
            fl[nf] = (newscope, r)
            for u in newscope: var_f.setdefault(u, set()).add(nf)
            nf += 1
            if len(newscope) <= 2:
                cand.append(nf - 1); fix_all()
                if scal == 0: return 0j
        else: scal = _nz(scal * complex(r))
        alive = {u for u in alive if var_f.get(u)}
    for i, (vs, t) in fl.items(): scal = _nz(scal * complex(t))
    if not np.isfinite(scal.real) or not np.isfinite(scal.imag): raise Inconclusive('non-finite amplitude')
    return scal * 2.0 ** LG[0]


def check_text(qasm_text, qc_path, x_trials=4, m_trials=2, seed=1, clean=False, tol=1e-7, deadline=None):
    """-> {'ok': bool, 'F': float, ...}. Raises Inconclusive if the check cannot be completed."""
    names, hl = load_qc(qc_path); n = len(names)
    try:
        prog = P.parse_text(qasm_text) if hasattr(P, 'parse_text') else _parse_text(qasm_text)
    except Exception as e:  # noqa: BLE001
        raise Inconclusive('parse: %s' % e)
    N, nb, ops = prog
    for i, op in enumerate(ops):   # the |A| = 1 normalization holds only for X-measurements (p = 1/2)
        if op[0] == 'm' and not (i > 0 and ops[i - 1][0] == 'h' and ops[i - 1][1] == op[1]):
            raise Inconclusive('Z-measurement')
    rng = random.Random(seed)
    ins = None
    for line in open(qc_path):
        if line.startswith('.i'): ins = line.split()[1:]
    xs = []
    for _ in range(x_trials):
        x = {q: rng.randint(0, 1) for q in range(n)}
        if clean and ins is not None:
            for q, nm in enumerate(names):
                if nm not in ins: x[q] = 0
        xs.append(x)
    worst = 0.0; refc = {}
    for _ in range(m_trials):
        m = [rng.randint(0, 1) for _ in range(nb)]; cs = []
        for x in xs:
            if deadline is not None and time.time() > deadline: raise Inconclusive('timeout')
            try:
                ps, _c = P.build(prog, x, outcomes=m)
            except Exception as e:  # noqa: BLE001
                raise Inconclusive('path sum: %s' % e)
            al = [q for q in range(N) if ps.alive[q]]
            if al[:n] != list(range(n)): return dict(ok=False, F=0.0, error='data wire not alive')
            if not ps.consistent: return dict(ok=False, F=0.0, error='inconsistent (zero) branch')
            tb = classical_target(hl, n, x)
            tg = {q: (tb[q] if q < n else 0) for q in al}
            amp = amp_sum(ps, tg, lg_extra=(-ps.nH + ps.nmeas) / 2, deadline=deadline)
            xk = tuple(sorted(x.items()))
            if xk not in refc:
                refc[xk] = -1.0 if classical_target(hl, n, x, with_sign=True)[1] else 1.0
            worst = max(worst, abs(abs(amp) - 1.0))
            cs.append(amp / refc[xk])
        for c in cs: worst = max(worst, abs(c - cs[0]))
    return dict(ok=worst < tol, F=max(0.0, 1 - worst), worst_dev=worst, n=n, N=N, measurements=nb)


def _parse_text(txt):
    """pathsum.parse takes a file path; write the text to a temporary file."""
    import os
    import tempfile
    fd, tp = tempfile.mkstemp(suffix='.qasm')
    try:
        os.write(fd, txt.encode()); os.close(fd)
        return P.parse(tp)
    finally:
        try: os.unlink(tp)
        except OSError: pass
