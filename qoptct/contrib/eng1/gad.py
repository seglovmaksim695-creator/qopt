"""eng1: Hadamard-gadgetization + TODD-reduction of Clifford+T circuits.

Whole circuit  ->  sum over paths:  state = sum_y  w^{f(x,y)} |L(x,y)>,  f = phase polynomial over N=n+h variables
(h = number of internal H gates; every H creates a fresh variable y).  Equivalent circuit (H-gadgets):
   h ancillas in |+>  ->  CX/phase network on N wires  ->  last h wires measured in the X basis (outcome 0 = the unitary).
f is reduced (TODD: odd parity columns, exact Clifford remainder), resynthesised (graysynth) onto the same final
linear map.  Equality is proven (a) algebraically: coefficients of f (canonical Z8-multilinear form), final wire masks and
constants are identical, extracted independently from the new circuit; (b) numerically by state-vector simulation of the
postselected circuit vs the original (N<=22).
"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, itertools, random, time, json, math
from pathlib import Path
ROOT = Path(_QPKG)
for p in (str(ROOT), str(ROOT / "portfolio")):
    if p not in sys.path: sys.path.insert(0, p)
import numpy as np
import qoptct, evo_lib, tpoly
pp = tpoly.pp
kernel = tpoly.kernel
coeffs = tpoly.coeffs
lin_ops = tpoly.lin_ops


# ---------------------------------------------------------------- path sum
def prep_gates(hl):
    """Toffoli-level list -> gates {h,x,cx,t,tdg,s,sdg,z,cz,ccz}; ccx -> h ccz h, adjacent H cancelled."""
    g = evo_lib.fold_h(hl)
    out = []
    for nm, q in g:
        q = tuple(q)
        if nm == "ccx":
            out += [("h", (q[2],)), ("ccz", q), ("h", (q[2],))]
        else:
            out.append((nm, q))
    return evo_lib.cancel_h(out)


def build(gates, n):
    """-> N, P (mask->weight mod 8), W (N masks: outputs 0..n-1, then measured forms), K (consts of outputs)."""
    W = [1 << i for i in range(n)]; K = [0] * n; P = {}; nv = n; meas = []
    for nm, qs in gates:
        if nm == "h":
            q = qs[0]; z = W[q]; c = K[q]; y = 1 << nv; nv += 1
            pp.add(P, z, 2); pp.add(P, y, 2); pp.add(P, z ^ y, -2)
            if c: pp.add(P, y, 4)
            meas.append(z); W[q] = y; K[q] = 0
        else:
            pp.block_poly([(nm, qs)], W, K, P)
    return nv, P, W + meas, K


def canon(P):
    return {k: v % 8 for k, v in coeffs(list(P.items())).items() if v % 8}


def rank_ok(W, N):
    basis = []
    for v in W:
        for b in basis:
            v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis) == N


# ---------------------------------------------------------------- TODD
def _rows_for(z, cols, rowsP, act, m, Q):
    rows = [(1 << m) - 1] + [rowsP[i] for i in act]
    zs = [i for i in act if z >> i & 1]
    zset = set(zs)
    seen = set()
    for i in zs:
        for j, k in itertools.combinations([x for x in act if x != i], 2):
            key = tuple(sorted((i, j, k)))
            if key in seen: continue
            seen.add(key)
            r = 0
            if z >> key[0] & 1: r ^= Q[key[1], key[2]]
            if z >> key[1] & 1: r ^= Q[key[0], key[2]]
            if z >> key[2] & 1: r ^= Q[key[0], key[1]]
            if r: rows.append(r)
    return rows


def _norm(cols):
    cnt = {}
    for c in cols:
        if c: cnt[c] = cnt.get(c, 0) ^ 1
    return [c for c, v in cnt.items() if v]


def todd(cols, N, rng, time_limit=20.0, plateau_max=40):
    """TODD with singles (z = col), pairs (z = col^col), kernel-vector choice by best count and plateau moves."""
    t0 = time.time()
    cols = _norm(cols); best = list(cols); plateau = 0
    while time.time() - t0 < time_limit and len(cols) >= 1:
        m = len(cols)
        rowsP = [0] * N
        for j, c in enumerate(cols):
            x = c
            while x:
                b = x & -x; rowsP[b.bit_length() - 1] |= 1 << j; x ^= b
        act = [i for i in range(N) if rowsP[i]]
        Q = {}
        for a, b in itertools.combinations(act, 2): Q[a, b] = rowsP[a] & rowsP[b]
        zs = {}
        for a in range(m):
            zs.setdefault(cols[a], []).append((a, None))
            for b in range(a + 1, m):
                zs.setdefault(cols[a] ^ cols[b], []).append((a, b))
        order = list(zs); rng.shuffle(order); order.sort(key=lambda z: bin(z).count("1"))
        improved = False
        moves = []
        for z in order:
            if z == 0: continue
            if time.time() - t0 > time_limit: break
            rows = _rows_for(z, cols, rowsP, act, m, Q)
            ker = kernel(rows, m)
            if not ker: continue
            # candidate y: basis vectors and a few random combos
            cand = list(ker)
            for _ in range(min(6, 2 * len(ker))):
                y = 0
                for kv in ker:
                    if rng.random() < 0.5: y ^= kv
                if y: cand.append(y)
            bestn, besty = m, None
            for y in cand:
                new = _norm([cols[j] ^ (z if y >> j & 1 else 0) for j in range(m)])
                if len(new) < bestn: bestn, besty = len(new), new
                elif len(new) == m and len(moves) < 30 and set(new) != set(cols): moves.append(new)
            if besty is not None:
                cols = besty; improved = True; plateau = 0; break
        if improved:
            if len(cols) < len(best): best = list(cols)
            continue
        if not moves or plateau >= plateau_max: break
        cols = rng.choice(moves); plateau += 1
    return best


def clifford_rest(old_odd, new_odd, even):
    """Total phase polynomial = new odd columns (weight 1) + even terms + exact Clifford difference."""
    Co, Cn = coeffs([(m, 1) for m in old_odd]), coeffs([(m, 1) for m in new_odd])
    D = {k: (Co.get(k, 0) - Cn.get(k, 0)) % 8 for k in set(Co) | set(Cn)}
    P = {}
    for m in new_odd: pp.add(P, m, 1)
    for m, w in even.items():
        if w: pp.add(P, m, w)
    for S, v in D.items():
        if not v: continue
        if len(S) == 3: assert v % 8 == 0, "cubic part"; continue
        if len(S) == 1: assert v % 2 == 0; pp.add(P, 1 << S[0], v)
        else:
            assert v % 4 == 0, "pair part"; a, b = S; h = v // 2
            pp.add(P, 1 << a, h); pp.add(P, 1 << b, h); pp.add(P, (1 << a) | (1 << b), -h)
    return P


def split_odd(P):
    odd = [m for m, w in P.items() if w % 2]
    even = {m: ((w - 1) % 8 if w % 2 else w) for m, w in P.items()}
    return odd, even


def synth(P, Wend, Kend, N, n, seeds=8, base_seed=0):
    """Phase polynomial P + final wires -> circuit on N wires (body only).
    Variants: graysynth / greedy with end-state penalty (c0005) + best linear repair; best by (T, CX)."""
    best = None
    LAM = (0.0, 0.3, 0.6, 1.0, 1.5)
    for s in range(seeds):
        rng = random.Random(base_seed + s)
        W = [1 << i for i in range(N)]
        try:
            if s % 3 == 0: out = pp.graysynth(dict(P), W, N, rng)
            elif s % 3 == 1: out = pp.greedy(dict(P), W, N, rng, Wend, LAM[(s // 3) % 5])
            else: out = pp.greedy(dict(P), W, N, rng, Wend, rng.uniform(0, 1.5), rng.uniform(0, 1.5), rng.uniform(0.1, 1.0))
        except Exception:
            out = None
        if out is None: continue
        M = [pp.solve(Wend, N, W[r]) for r in range(N)]
        ops = pp.lin_best(M, N, range(N), rng, 6)
        out = out + [("cx", o) for o in ops] + [("x", (q,)) for q in range(n) if Kend[q]]
        out = evo_lib.fold_cancel(out, N)
        k = (evo_lib.tcount(out), evo_lib.cx(out))
        if best is None or k < best[0]: best = (k, out)
    return None if best is None else best[1]


# ---------------------------------------------------------------- verification
def extract(body, N):
    """Independently re-extract (canonical f, wire masks, consts) from an H-free body on N wires."""
    W = [1 << i for i in range(N)]; K = [0] * N; P = {}
    pp.block_poly(body, W, K, P)
    return canon(P), W, K


def check_alg(body, N, n, P0, W0, K0):
    c1, W1, K1 = extract(body, N)
    return c1 == canon(P0) and W1 == list(W0) and K1[:n] == list(K0[:n])


def full_circuit(body, N, n):
    h = N - n
    return [("h", (q,)) for q in range(n, N)] + list(body) + [("h", (q,)) for q in range(n, N)]


def rand_state(k, rng):
    v = rng.normal(size=2 ** k) + 1j * rng.normal(size=2 ** k)
    return v / np.linalg.norm(v)


def check_sim(body, ref, N, n, trials=3, seed=1):
    """Postselected (all ancillas -> |0> after final H) full circuit vs ref: overlap 1 and norm^2 = 2^-h."""
    import gidney
    if N > 22: return None
    h = N - n; rng = np.random.default_rng(seed); circ = full_circuit(body, N, n)
    for _ in range(trials):
        v = rand_state(n, rng)
        psi = np.zeros([2] * N, complex)
        psi[tuple([slice(None)] * n + [0] * h)] = v.reshape([2] * n)
        psi = psi.reshape(-1)
        for g, q in circ: psi = gidney.apply(psi, N, g, tuple(q))
        got = psi.reshape([2] * N)[tuple([slice(None)] * n + [0] * h)].reshape(-1)
        ex = v.copy()
        for g, q in ref: ex = gidney.apply(ex, n, g, tuple(q))
        nrm = np.linalg.norm(got)
        if abs(nrm ** 2 - 2.0 ** (-h)) > 1e-9: return False
        if abs(abs(np.vdot(ex, got)) / nrm - 1) > 1e-9: return False
    return True


# ---------------------------------------------------------------- main entry
def improve_gates(hl, n, time_limit=20.0, restarts=3, seeds=6, log=None, exploit=None):
    t0 = time.time()
    gates = prep_gates(hl)
    N, P, W, K = build(gates, n)
    assert rank_ok(W, N)
    odd, even = split_odd(P)
    res = dict(n=n, h=N - n, N=N, t_poly=len(odd))
    best = None
    for r in range(restarts):
        rng = random.Random(100 + r)
        new = todd(odd, N, rng, time_limit)
        P2 = clifford_rest(odd, new, even)
        body = synth(P2, W, K, N, n, seeds=seeds, base_seed=r * 100)
        if body is None: continue
        k = (evo_lib.tcount(body), evo_lib.cx(body))
        if log: log(f"  restart {r}: cols {len(odd)}->{len(new)} T={k[0]} CX={k[1]} ({time.time()-t0:.0f}s)")
        if best is None or k < best[0]: best = (k, body)
    body = best[1]
    res.update(T=best[0][0], CX=best[0][1], extra_qubits=N - n, measurements=N - n)
    ok_alg = check_alg(body, N, n, P, W, K)
    ref = evo_lib.lower_std(hl)
    ok_sim = check_sim(body, ref, N, n)
    res.update(alg_proof=bool(ok_alg), sim=ok_sim, seconds=round(time.time() - t0, 1))
    return res, body, (N, n)


def improve(path, **kw):
    """Entry point: improve(path) -> dict (T, CX, extra qubits, measurements, proof flags)."""
    names, hl = qoptct.load_qc(path)
    res, body, _ = improve_gates(hl, len(names), **kw)
    res["name"] = Path(path).stem
    return res


if __name__ == "__main__":
    for p in sys.argv[1:]:
        print(json.dumps(improve(p, log=lambda s: print(s, flush=True))), flush=True)
