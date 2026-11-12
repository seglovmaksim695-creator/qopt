"""Fast TODD with randomised constraint rows (numpy).
Condition for A' = A + z y^T to keep the cubic tensor sum_k a_k^{(x)3} (mod 2): for all r,s,t
   sum_k y_k c_k(r,s,t) = 0,  c_k = z_r a_s a_t + a_r z_s a_t + a_r a_s z_t + z_r z_s a_t + z_r a_s z_t + a_r z_s z_t + z_r z_s z_t
(a_r = a_k . r mod 2, z_r = z . r).  The rows c(r,s,t) for random r,s,t span the full constraint space w.h.p. when
there are ~6m of them; kernel found by GF(2) elimination; then re-checked on fresh random rows (cutting plane).
The exact tensor/Clifford remainder is verified later (gad.clifford_rest asserts, final algebraic proof)."""
import time, random
import numpy as np


def cols_to_A(cols, N):
    m = len(cols)
    A = np.zeros((N, m), np.uint8)
    for j, c in enumerate(cols):
        x = c
        while x:
            b = x & -x; A[b.bit_length() - 1, j] = 1; x ^= b
    return A


def A_to_cols(A):
    N, m = A.shape
    return [int(sum(1 << int(i) for i in np.nonzero(A[:, j])[0])) for j in range(m)]


def make_rows(A, Af, z, S, nrng):
    N, m = A.shape
    R = nrng.integers(0, 2, size=(3, S, N)).astype(np.float32)
    zf = z.astype(np.float32)
    Ar, As, At = [((R[i] @ Af) % 2).astype(np.uint8) for i in range(3)]
    zr, zs, zt = [((R[i] @ zf) % 2).astype(np.uint8)[:, None] for i in range(3)]
    c = (zr & As & At) ^ (Ar & zs & At) ^ (Ar & As & zt) ^ (zr & zs & At) ^ (zr & zt & As) ^ (zs & zt & Ar) ^ (zr & zs & zt)
    return c


def kernel_np(M):
    M = M.copy(); S, m = M.shape; piv = []; r = 0
    for c in range(m):
        if r >= S: break
        idx = np.nonzero(M[r:, c])[0]
        if len(idx) == 0: continue
        p = idx[0] + r
        if p != r: M[[r, p]] = M[[p, r]]
        rows = np.nonzero(M[:, c])[0]; rows = rows[rows != r]
        if len(rows): M[rows] ^= M[r]
        piv.append(c); r += 1
    ps = set(piv); basis = []
    for f in range(m):
        if f in ps: continue
        y = np.zeros(m, np.uint8); y[f] = 1
        for i, pc in enumerate(piv):
            if M[i, f]: y[pc] = 1
        basis.append(y)
    return basis


def kernel_for(A, Af, z, nrng, mult=3):
    """Basis of {y : tensor preserved} (up to negligible failure probability), via cutting planes."""
    N, m = A.shape
    M = make_rows(A, Af, z, mult * m + 30, nrng)
    for _ in range(6):
        basis = kernel_np(M)
        if not basis: return basis
        # test a random kernel element on fresh rows
        y = np.zeros(m, np.uint8)
        for b in basis:
            if nrng.random() < 0.5: y ^= b
        T = make_rows(A, Af, z, 2 * m + 30, nrng)
        bad = np.nonzero((T.astype(np.int32) @ y.astype(np.int32)) % 2)[0]
        if len(bad) == 0: return basis
        M = np.vstack([M, T])
    return kernel_np(M)


def tens(cols, N):
    A = cols_to_A(cols, N).astype(np.float32); m = A.shape[1]
    if m == 0: return np.zeros((N * N, N), np.uint8)
    B = (A[:, None, :] * A[None, :, :]).reshape(N * N, m)
    return ((B @ A.T) % 2).astype(np.uint8)


def norm_cols(cols):
    cnt = {}
    for c in cols:
        if c: cnt[c] = cnt.get(c, 0) ^ 1
    return [c for c, v in cnt.items() if v]


def todd_fast(cols, N, rng, time_limit=30.0, plateau_max=300, tries=None, log=None, hook=None):
    t0 = time.time(); nrng = np.random.default_rng(rng.randrange(1 << 30))
    cols = norm_cols(cols); best = list(cols); plateau = 0; stuck_rounds = 0
    while time.time() - t0 < time_limit and cols:
        m = len(cols)
        A = cols_to_A(cols, N); Af = A.astype(np.float32); Tcur = None
        pairs = [(a, b) for a in range(m) for b in range(a + 1, m)] + [(a, None) for a in range(m)]
        rng.shuffle(pairs)
        seenz = set(); improved = False; moves = []
        for k, (a, b) in enumerate(pairs):
            if time.time() - t0 > time_limit: break
            if tries and k >= tries: break
            zc = cols[a] if b is None else cols[a] ^ cols[b]
            if zc == 0 or zc in seenz: continue
            seenz.add(zc)
            z = A[:, a].copy() if b is None else A[:, a] ^ A[:, b]
            ker = kernel_for(A, Af, z, nrng)
            if not ker: continue
            cand = []
            for kv in ker:
                if (b is None and kv[a]) or (b is not None and kv[a] != kv[b]): cand.append(kv)
            ys = list(cand)
            for _ in range(8 if hook is None else 40):
                y = np.zeros(m, np.uint8)
                for kv in ker:
                    if nrng.random() < 0.5: y ^= kv
                if (b is None and y[a]) or (b is not None and y[a] != y[b]): ys.append(y)
            cands = []
            for y in ys:
                new = norm_cols([cols[j] ^ (zc if y[j] else 0) for j in range(m)])
                if len(new) < m: cands.append((len(new), new))
            cands.sort(key=lambda t: t[0]); besty = None
            for ln, new in cands[:(4 if hook is None else 25)]:
                if hook is not None and not hook.can(new): continue
                if Tcur is None: Tcur = tens(cols, N)
                if np.array_equal(tens(new, N), Tcur) and (hook is None or hook.commit(new)): besty = new; break
            if besty is not None:
                cols = besty; improved = True; plateau = 0; break
            if len(moves) < 40:  # plateau candidates: random kernel element (same or smaller count, new columns)
                for _ in range(2):
                    y = np.zeros(m, np.uint8)
                    for kv in ker:
                        if nrng.random() < 0.5: y ^= kv
                    if not y.any(): continue
                    new = norm_cols([cols[j] ^ (zc if y[j] else 0) for j in range(m)])
                    if len(new) <= m and set(new) != set(cols): moves.append(new)
        if improved:
            if len(cols) < len(best): best = list(cols)
            if log: log(f"    todd m={len(cols)} t={time.time()-t0:.1f}s")
            continue
        if not moves or plateau >= plateau_max: break
        if moves:
            mv = rng.choice(moves); Tc = tens(cols, N)
            if np.array_equal(tens(mv, N), Tc) and (hook is None or (hook.can(mv) and hook.commit(mv))): cols = mv; plateau += 1
            else: plateau += 1
    return best
