"""c10: independent check of large deterministic Clifford+T circuits with measurements (30-70+ qubits), no state vector.

Method: sum-over-paths (Amy) + exact variable elimination (bucket elimination on tensors, width ~10-15).
For a basis input x and CHOSEN measurement outcomes m (a branch) the data state is
   2^{-nH/2} * sum_y w^{F(y)} |W y + e>       (see pathsum.py)
and it is compared with the reference (the .qc circuit lowered by lower_std, also as a path sum): in every branch
   A_det(d) * 2^{nm/2} = c_m * A_ref(d)  for all d, with c_m independent of x  (=> U_branch = c_m U).
XOR/AND conditions in `if` are supported (pathsum.parse_cond).
Branches are chosen uniformly (not by weight): for correct deterministic protocols all outcomes are equally likely;
a zero branch is counted separately (zero_branches) and another branch is drawn.

v0.4 fixes (lab agents b9, b1, b8):
  * AmpPlan: the elimination plan and factor tables are built once per (path sum, target wires); each point d only
    contracts tensors (x13 faster on qcla_com_7, identical results on 120 mutants);
  * scale factors are kept as powers of two (exponents), and intermediate tensors are renormalised by powers of two:
    2^(-nH/2) in double underflowed to 0 for thousands of H, 0*inf gave nan, and nan comparisons made the check
    vacuous (F=1 with undetected canaries). Non-finite values now fail the check;
  * sampled mode: when the image of the data map has more than `cap` points the norm on a subsample is not 1, which
    gave false rejections; now the ratio test is applied on all sampled points (zeros included) and the subsample
    norm is replaced by the EXACT branch norm (one doubled path sum, branch_norm2): ||psi||^2 must equal |c_m|^2,
    so weight on unsampled points is detected (a ratio test alone accepted a T->Tdg mutant of cuccaro_16).
    The result carries sampled>0;
  * for classical references the expected output point is always added to the sample (otherwise random image points
    miss the single support point of a correct circuit and every branch looked like a zero branch).
"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, math, random, cmath, time
import numpy as np
sys.path[:0] = [_QPKG, (_QPKG + '/portfolio'), (_QPKG + '/contrib/c10')]
import pathsum as P

OM = [cmath.exp(1j * math.pi * k / 4) for k in range(8)]
_RENORM_HI = 2.0 ** 200
_RENORM_LO = 2.0 ** -200


def _bits(m):
    r = []; i = 0
    while m:
        if m & 1: r.append(i)
        m >>= 1; i += 1
    return r


def _par_table(w):
    """parity of w bits -> array of shape (2,)*w"""
    idx = np.indices((2,) * w).reshape(w, -1).sum(0) & 1
    return idx.reshape((2,) * w)


_PT = {}


def par(w):
    if w not in _PT: _PT[w] = _par_table(w)
    return _PT[w]


def _pow2(e):
    """2.0**e for a real exponent (exact for integers and half-integers up to rounding of sqrt2)."""
    f = math.floor(e)
    return math.ldexp(1.0, int(f)) * 2.0 ** (e - f)


class AmpPlan:
    """Plan for S(d) = sum_{y: W y + e = d} w^{F(y)} for a FIXED path sum ps and a set of target wires tkeys.
    The factor structure and the elimination order do not depend on the values of d (only the target-factor tables
    change: parity == want), so they are computed once; eval(targets) only contracts tensors.
    Intermediate tensors are renormalised by powers of two and the exponent is tracked separately."""

    def __init__(self, ps, tkeys, maxw=24):
        self.ps = ps; self.tkeys = list(tkeys)
        self.scal0 = OM[ps.gph % 8]
        self.dead = not ps.consistent
        factors = []   # (scope, table or None for target factors)
        for m, c in ps.lin.items():
            vs = _bits(m)
            assert len(vs) <= maxw, "parity too wide: %d" % len(vs)
            factors.append((vs, np.where(par(len(vs)) == 1, OM[c % 8], 1.0 + 0j)))
        for (ma, mb) in ps.quad:
            va, vb = _bits(ma), _bits(mb); vs = sorted(set(va) | set(vb))
            assert len(vs) <= maxw, "parity too wide: %d" % len(vs)
            pos = {v: i for i, v in enumerate(vs)}
            grid = np.indices((2,) * len(vs)).reshape(len(vs), -1)
            pa = grid[[pos[v] for v in va]].sum(0) & 1; pb = grid[[pos[v] for v in vb]].sum(0) & 1
            factors.append((vs, np.where((pa & pb) == 1, -1.0 + 0j, 1.0 + 0j).reshape((2,) * len(vs))))
        for m, bit in ps.cons:
            vs = _bits(m); assert len(vs) <= maxw, "parity too wide: %d" % len(vs)
            factors.append((vs, (par(len(vs)) == bit).astype(complex)))
        self.tinfo = []   # (wire, factor index or None, const)
        self.ttab = {}    # factor index -> {want: table}
        for q in self.tkeys:
            m, c = ps.w[q]
            if m == 0:
                self.tinfo.append((q, None, c)); continue
            vs = _bits(m); assert len(vs) <= maxw, "parity too wide: %d" % len(vs)
            idx = len(factors); factors.append((vs, None))
            p = par(len(vs)); self.ttab[idx] = {0: (p == 0).astype(complex), 1: (p == 1).astype(complex)}
            self.tinfo.append((q, idx, c))
        self.base = {i: f[1] for i, f in enumerate(factors) if f[1] is not None}
        var_f = {}; alive = set(); fl = {}
        for i, (vs, t) in enumerate(factors):
            fl[i] = vs
            for v in vs: var_f.setdefault(v, set()).add(i); alive.add(v)
        nxt = len(factors)
        free_vars = set(range(ps.k)) - alive
        self.exp2 = len(free_vars)          # factor 2^(number of free variables), kept as an exponent
        self.steps = []
        while alive:
            best = None
            for v in alive:
                sc = set()
                for i in var_f[v]: sc |= set(fl[i])
                key = (len(sc), len(var_f[v]))
                if best is None or key < best[0]: best = (key, v, sc)
            _, v, sc = best
            fi = sorted(var_f[v]); scope = sorted(sc)
            lab = {u: j for j, u in enumerate(scope)}
            out = [lab[u] for u in scope if u != v]
            assert len(scope) <= 50, "elimination width %d" % len(scope)
            labs = [(i, [lab[u] for u in fl[i]]) for i in fi]
            newscope = [u for u in scope if u != v]
            for i in fi:
                for u in fl[i]: var_f[u].discard(i)
                del fl[i]
            alive.discard(v); del var_f[v]
            if newscope:
                fl[nxt] = newscope; newid = nxt
                for u in newscope: var_f.setdefault(u, set()).add(nxt)
                nxt += 1
            else:
                newid = None
            self.steps.append((labs, out, newid))
            alive = {u for u in alive if var_f.get(u)}

    def eval(self, targets, log2_scale=0.0):
        """S(targets) * 2^log2_scale, computed with power-of-two renormalisation (no overflow/underflow)."""
        if self.dead: return 0j
        scal = complex(self.scal0); e2 = float(self.exp2) + log2_scale
        tb = dict(self.base)
        for q, idx, c in self.tinfo:
            bit = targets[q]
            if idx is None:
                if c != bit: return 0j
            else:
                tb[idx] = self.ttab[idx][bit ^ c]
        ein = np.einsum
        for labs, out, newid in self.steps:
            args = []
            for i, lb in labs:
                args.append(tb.pop(i)); args.append(lb)
            r = ein(*args, out, optimize=False)
            if newid is not None:
                mx = float(np.max(np.abs(r))) if r.size else 0.0
                if mx > _RENORM_HI or 0.0 < mx < _RENORM_LO:
                    k = math.frexp(mx)[1]
                    r = r * math.ldexp(1.0, -k); e2 += k
                tb[newid] = r
            else:
                scal *= complex(r)
                mx = abs(scal)
                if mx > _RENORM_HI or 0.0 < mx < _RENORM_LO:
                    k = math.frexp(mx)[1]
                    scal = scal * math.ldexp(1.0, -k); e2 += k
        if scal == 0: return 0j
        return scal * _pow2(e2)


class _Doubled:
    """Doubled path sum for the branch norm (copy y and copy y' of the path variables)."""
    w = {}


def branch_norm2(ps, alive, log2_scale=0.0, maxw=24):
    """||psi||^2 = sum_d |S(d)|^2 * 2^log2_scale, exactly as one doubled sum:
    sum_{y, y'} w^{F(y) - F(y')} [W y = W y' on the alive wires] (constraints on both copies)."""
    k = ps.k
    D = _Doubled(); D.k = 2 * k; D.gph = 0; D.consistent = ps.consistent
    D.lin = {}
    for m, c in ps.lin.items():
        D.lin[m] = c % 8; D.lin[m << k] = (-c) % 8
    D.quad = {}
    for (a, b) in ps.quad:
        D.quad[(a, b)] = 1; D.quad[(a << k, b << k)] = 1
    D.cons = list(ps.cons) + [(m << k, bit) for m, bit in ps.cons]
    for q in alive:
        m = ps.w[q][0]
        if m: D.cons.append((m | (m << k), 0))
    v = AmpPlan(D, [], maxw).eval({}, log2_scale)
    return v.real


def amplitude_sum(ps, targets, maxw=24):
    """S = sum_{y: W y + e = targets} w^{F(y)} (without the 2^{-nH/2} normalisation). targets: {wire: bit}."""
    return AmpPlan(ps, list(targets.keys()), maxw).eval(targets)


# ------------------------------------------------------------ reference
def ref_ps(gates, n, x):
    ps = P.PS(n, x)
    for g, q in gates:
        if g == "h": ps.hadamard(q[0])
        else: ps.gate(g, q)
    return ps


def classical_target(hl, n, x):
    """Expected output bits for classical references (x/cx/ccx, h-z-h frames as Toffoli; diagonal gates outside
    frames only add phases). Raises ValueError if the reference is not of this kind."""
    b = [x.get(q, 0) for q in range(n)]
    frame = set()
    for g, q in hl:
        q = tuple(q)
        if g == "h":
            frame ^= {q[0]}
            continue
        fr = [u for u in q if u in frame]
        if g in ("x", "cx", "ccx"):
            if fr: raise ValueError("%s on an H frame" % g)
            if all(b[u] for u in q[:-1]): b[q[-1]] ^= 1
        elif g in ("z", "cz", "ccz"):
            if not fr: continue
            if len(fr) != 1: raise ValueError("two H frames in one Z gate")
            t = fr[0]; ctl = [u for u in q if u != t]
            if all(b[u] for u in ctl): b[t] ^= 1
        elif g in ("t", "tdg", "s", "sdg"):
            if fr: raise ValueError("phase gate on an H frame")
        else:
            raise ValueError("not a classical gate " + g)
    if frame: raise ValueError("open H frame")
    return b


def _with_low(basis, e, nlow, want_low, rng, extra=16):
    """Image points (e + span(basis)) whose low nlow bits equal want_low (plus a few random kernel shifts)."""
    lowmask = (1 << nlow) - 1
    piv = []; kern = []
    for b in basis:
        v = b
        for pb, pv in piv:
            if v >> pb & 1: v ^= pv
        if v & lowmask:
            piv.append(((v & lowmask).bit_length() - 1, v))
        elif v:
            kern.append(v)
    s = (e ^ want_low) & lowmask; acc = 0
    for pb, pv in piv:
        if s >> pb & 1:
            s ^= pv & lowmask; acc ^= pv
    if s: return []
    p0 = e ^ acc; out = [p0]
    for _ in range(extra if kern else 0):
        p = p0
        for k in kern:
            if rng.random() < 0.5: p ^= k
        out.append(p)
    return out


def image_points(ps, wires, cap, rng, want_low=None, nlow=0):
    """Points of the image of the affine map y -> (W y + e) on the given wires: list of ints (bit j = wire wires[j]).
    If want_low is given and the image is sampled, points whose low nlow bits equal want_low are always included."""
    cols = {}
    for j, q in enumerate(wires):
        for v in _bits(ps.w[q][0]): cols[v] = cols.get(v, 0) | (1 << j)
    basis = []
    for c in cols.values():
        for b in basis: c = min(c, c ^ b)
        if c: basis.append(c)
    e = sum(ps.w[q][1] << j for j, q in enumerate(wires))
    r = len(basis)
    if 2 ** r <= cap:
        pts = [e]
        for b in basis: pts = pts + [p ^ b for p in pts]
        return pts, True
    pts = {e}
    if want_low is not None:
        pts.update(_with_low(basis, e, nlow, want_low, rng))
    while len(pts) < cap:
        p = e
        for b in basis:
            if rng.random() < 0.5: p ^= b
        pts.add(p)
    return sorted(pts), False


def check(qasm, qc_path, x_trials=4, m_trials=4, seed=1, verbose=True, cap=2048, tol=1e-7, max_zero_tries=40):
    """Main check. Returns dict(ok, F, worst_dev, sampled, ...). F = 1 - largest deviation.
    For each pair (x, branch m): amplitudes of the circuit on ALL image points of the data map (2^rank points, or a
    random subsample of size cap -> sampled mode), the reference on the support points; the branch vector must be
    c_m * (reference output), c_m independent of x (also in modulus: the branch probability must not depend on x)."""
    import qoptct, evo_lib
    names, hl = qoptct.load_qc(qc_path); n = len(names); ref = evo_lib.lower_std(hl)
    prog = P.parse(qasm); N, nb, ops = prog
    rng = random.Random(seed)
    worst = 0.0; zero_branches = 0; t0 = time.time(); complete = True; maxsupp = 0; sampled = 0
    xs = [{q: rng.randint(0, 1) for q in range(n)} for _ in range(x_trials)]
    ms = []
    while len(ms) < m_trials: ms.append([rng.randint(0, 1) for _ in range(nb)])
    rps = [ref_ps(ref, n, x) for x in xs]
    refcache = [dict() for _ in xs]
    refplan = [None] * len(xs)
    targets = []
    for x in xs:
        try:
            b = classical_target(hl, n, x); targets.append(sum(v << q for q, v in enumerate(b)))
        except ValueError:
            targets.append(None)

    def aref(xi, d):
        if d not in refcache[xi]:
            rp = rps[xi]
            if refplan[xi] is None: refplan[xi] = AmpPlan(rp, list(range(n)))
            refcache[xi][d] = refplan[xi].eval({q: (d >> q) & 1 for q in range(n)}, -rp.nH / 2)
        return refcache[xi][d]

    def bad(v):
        return not (math.isfinite(v.real) and math.isfinite(v.imag))
    pending = list(ms)
    done_m = []
    while pending:
        m = pending.pop(0)
        vecs = []; branch_zero = False
        for xi, x in enumerate(xs):
            try: ps, c = P.build(prog, x, outcomes=m)
            except ValueError as e:
                return dict(ok=False, F=0.0, error=str(e))
            al = [q for q in range(N) if ps.alive[q]]
            assert al[:n] == list(range(n)), "data wires 0..n-1 must stay alive"
            log2s = (ps.nmeas - ps.nH) / 2
            pts, full = image_points(ps, al, cap, rng, targets[xi], n); complete &= full
            plan = AmpPlan(ps, al)
            amp = {}
            for d in pts:
                tg = {q: (d >> j) & 1 for j, q in enumerate(al)}
                amp[d] = plan.eval(tg, log2s)
            if any(bad(v) for v in amp.values()):
                return dict(ok=False, F=0.0, error="non-finite amplitude", sampled=sampled)
            mx = max(abs(v) for v in amp.values())
            if mx < 1e-9:
                branch_zero = True; break
            supp = [d for d, v in amp.items() if abs(v) > 1e-9 * max(1.0, mx)]
            maxsupp = max(maxsupp, len(supp))
            # strictness: alive ancillas must end in one basis state (not entangled with the data)
            hi = {d >> n for d in supp}
            if len(hi) > 1: worst = max(worst, 1.0)
            if not full:
                # sampled mode: ratio test on ALL sampled points (zeros included) + exact branch norm (below)
                sampled += 1
                d0 = max(supp, key=lambda d: abs(amp[d])); h0 = d0 >> n
                ar = {d: (aref(xi, d & ((1 << n) - 1)) if (d >> n) == h0 else 0j) for d in pts}
                if any(bad(v) for v in ar.values()):
                    return dict(ok=False, F=0.0, error="non-finite reference amplitude", sampled=sampled)
                if abs(ar[d0]) < 1e-9: worst = max(worst, 1.0); vecs.append(complex(1.0)); continue
                cc = amp[d0] / ar[d0]
                dev = max(abs(amp[d] - cc * ar[d]) / abs(cc) for d in pts)
                # exact branch norm: psi = cc * ref (x) |h0>  requires ||psi||^2 = |cc|^2 (no weight outside the sample)
                try:
                    nb2 = branch_norm2(ps, al, ps.nmeas - ps.nH)
                except AssertionError as e:
                    return dict(ok=False, F=None, error="sampled mode: branch norm not computable (%s)" % e,
                                sampled=sampled)
                dn = abs(nb2 - abs(cc) ** 2) / abs(cc) ** 2
                dev = max(dev, dn)
                worst = max(worst, dev if math.isfinite(dev) else 1.0); vecs.append(cc); continue
            h0 = min(hi)
            ar = {d: (aref(xi, d & ((1 << n) - 1)) if (d >> n) == h0 else 0j) for d in supp}
            if any(bad(v) for v in ar.values()):
                return dict(ok=False, F=0.0, error="non-finite reference amplitude", sampled=sampled)
            nrm = sum(abs(v) ** 2 for v in ar.values())
            worst = max(worst, abs(nrm - 1))
            d0 = max(supp, key=lambda d: abs(ar[d]))
            if abs(ar[d0]) < 1e-9:
                worst = max(worst, 1.0); vecs.append(complex(1.0)); continue
            cc = amp[d0] / ar[d0]
            dev = max(abs(amp[d] - cc * ar[d]) / abs(cc) for d in supp)
            worst = max(worst, dev if math.isfinite(dev) else 1.0)
            vecs.append(cc)
        if branch_zero:
            zero_branches += 1
            if zero_branches > max_zero_tries * m_trials:
                return dict(ok=False, F=0.0, error="too many zero branches", sampled=sampled)
            pending.append([rng.randint(0, 1) for _ in range(nb)]); continue
        for cc in vecs:
            worst = max(worst, abs(cc - vecs[0]) / abs(vecs[0]))
        done_m.append(abs(vecs[0]) ** 2 * _pow2(-ps.nmeas))
        if verbose: print("  branch done t=%.1fs worst=%.2e supp<=%d" % (time.time() - t0, worst, maxsupp), flush=True)
    if not math.isfinite(worst): worst = float("inf")
    return dict(ok=bool(worst < tol), F=max(0.0, 1 - worst), worst_dev=worst, zero_branches=zero_branches,
                complete_image=complete, sampled=sampled, branch_prob=done_m, n=n, N=N, measurements=nb,
                time=time.time() - t0)


if __name__ == "__main__":
    qasm, qc = sys.argv[1], sys.argv[2]
    xt = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    mt = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    print(check(qasm, qc, x_trials=xt, m_trials=mt))
