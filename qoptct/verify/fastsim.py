# SPDX-License-Identifier: Apache-2.0
"""Fast dense state-vector simulator for Clifford+T gate lists.

Operations work in place on reshape views (no full copy per gate); the 1/sqrt2 factors of H are accumulated in one
scalar; runs of gates without H (x/cx and diagonal gates) are fused into one affine permutation plus one phase
function, processed in blocks (memory ~ two vectors).
apply_gates(psi, n, gates) -> psi (flat complex128; MODIFIED IN PLACE or replaced by a new array; the caller must pass
its own copy). Qubit q is axis q of reshape([2]*n) (qubit 0 is the most significant index bit).
"""
import numpy as np

S2 = 1 / np.sqrt(2)
_PHV = {"t": 1, "tdg": 7, "s": 2, "sdg": 6, "z": 4}
_W8 = np.exp(1j * np.pi / 4 * np.arange(8))
_PHC = {k: _W8[v] for k, v in _PHV.items()}
_OK = {"h", "x", "cx", "ccx", "cz", "ccz", "t", "tdg", "s", "sdg", "z"}
FUSE_MIN = 6         # runs of at least this many gates are fused
BLOCK = 1 << 20      # elements per block when fusing

_SCR = [None]


def _scratch(size):
    s = _SCR[0]
    if s is None or s.size < size:
        s = _SCR[0] = np.empty(size, np.complex128)
    return s[:size]


def _h(psi, n, q):
    v = psi.reshape(1 << q, 2, 1 << (n - q - 1))
    a = v[:, 0, :]
    b = v[:, 1, :]
    t = _scratch(a.size).reshape(a.shape)
    np.subtract(a, b, out=t)
    np.add(a, b, out=a)
    b[...] = t


def _phase(psi, n, q, f):
    psi.reshape(1 << q, 2, 1 << (n - q - 1))[:, 1, :] *= f


def _x(psi, n, q):
    v = psi.reshape(1 << q, 2, 1 << (n - q - 1))
    a = v[:, 0, :]
    b = v[:, 1, :]
    t = _scratch(a.size).reshape(a.shape)
    t[...] = a
    a[...] = b
    b[...] = t


def _cx(psi, n, c, t):
    if c < t:
        v = psi.reshape(1 << c, 2, 1 << (t - c - 1), 2, 1 << (n - t - 1))
        a = v[:, 1, :, 0, :]
        b = v[:, 1, :, 1, :]
    else:
        v = psi.reshape(1 << t, 2, 1 << (c - t - 1), 2, 1 << (n - c - 1))
        a = v[:, 0, :, 1, :]
        b = v[:, 1, :, 1, :]
    tmp = _scratch(a.size).reshape(a.shape)
    tmp[...] = a
    a[...] = b
    b[...] = tmp


def _multi(psi, n, g, q):
    p = psi.reshape([2] * n)
    if g in ("ccx", "mcx"):
        *c, t = q
        idx = [slice(None)] * n
        for x in c:
            idx[x] = 1
        sub = p[tuple(idx)]
        ax = t - sum(1 for x in c if x < t)
        sub[...] = np.flip(sub, axis=ax).copy()
    else:  # cz, ccz, mcz
        idx = [slice(None)] * n
        for x in q:
            idx[x] = 1
        p[tuple(idx)] *= -1


def _single(psi, n, g, q):
    if g == "h":
        _h(psi, n, q[0])
    elif g == "x":
        _x(psi, n, q[0])
    elif g == "cx":
        _cx(psi, n, q[0], q[1])
    elif g in _PHC:
        _phase(psi, n, q[0], _PHC[g])
    else:
        _multi(psi, n, g, q)


def _table(m, bits):
    """parity table par(x & m) for x < 2^bits (built by doubling)."""
    t = np.zeros(1, np.uint8)
    for b in range(bits):
        t = np.concatenate([t, t ^ np.uint8((m >> b) & 1)])
    return t


def _fuse(psi, n, seg):
    """seg: gates from {x, cx, t, tdg, s, sdg, z, cz, ccz}. Returns a new vector."""
    W = [1 << q for q in range(n)]
    K = [0] * n
    P = {}
    c0 = 0

    def term(qs, v):
        nonlocal c0
        m = 0
        k = 0
        for q in qs:
            m ^= W[q]
            k ^= K[q]
        if k:
            c0 += v
            P[m] = P.get(m, 0) - v
        else:
            P[m] = P.get(m, 0) + v
    for g, q in seg:
        if g == "cx":
            a, b = q
            W[b] ^= W[a]
            K[b] ^= K[a]
        elif g == "x":
            K[q[0]] ^= 1
        elif g in _PHV:
            term(q, _PHV[g])
        elif g == "cz":
            a, b = q
            term((a,), 2)
            term((b,), 2)
            term((a, b), -2)
        else:  # ccz
            a, b, c = q
            for s, v in (((a,), 1), ((b,), 1), ((c,), 1), ((a, b), -1), ((a, c), -1), ((b, c), -1), ((a, b, c), 1)):
                term(s, v)
    L = n // 2
    Hh = n - L

    def flat(m):  # qubit mask -> bit mask of the flat index (qubit i = bit n-1-i)
        r = 0
        for i in range(n):
            if m >> i & 1:
                r |= 1 << (n - 1 - i)
        return r
    # y_q = par(W[q] & x) ^ K[q]; output index = sum y_q << (n-1-q): an affine map of x
    col = [0] * n
    for q in range(n):
        w = W[q]
        for i in range(n):
            if w >> i & 1:
                col[i] |= 1 << (n - 1 - q)
    const = 0
    for q in range(n):
        if K[q]:
            const |= 1 << (n - 1 - q)

    def lintab(bits, off):
        t = np.zeros(1, np.int64)
        for b in range(bits):
            t = np.concatenate([t, t ^ col[n - 1 - (b + off)]])
        return t
    IL = lintab(L, 0)
    IH = lintab(Hh, L) ^ const
    terms = []
    for m, c in P.items():
        c %= 8
        if c == 0:
            continue
        fm = flat(m)
        if fm == 0:
            c0 += c
            continue
        terms.append((c, _table(fm & ((1 << L) - 1), L), _table(fm >> L, Hh)))
    out = np.empty_like(psi)
    pv = psi.reshape(1 << Hh, 1 << L)
    hb = max(1, BLOCK >> L)
    for h0 in range(0, 1 << Hh, hb):
        h1 = min(1 << Hh, h0 + hb)
        e = np.full((h1 - h0, 1 << L), c0 & 7, np.uint8)
        for c, tl, th in terms:
            e += (th[h0:h1, None] ^ tl[None, :]) * np.uint8(c)
        e &= 7
        val = pv[h0:h1] * _W8[e]
        idx = IH[h0:h1, None] ^ IL[None, :]
        out[idx.ravel()] = val.ravel()
    return out


def _other(psi, n, g, q):
    """Gates outside the fused set: mcx/mcz with more controls, y, swap, id."""
    if g in ("mcx", "mcz"):
        _multi(psi, n, g, q)
        return psi
    if g == "y":
        _x(psi, n, q[0])
        p = psi.reshape(1 << q[0], 2, 1 << (n - q[0] - 1))
        p[:, 0, :] *= -1j
        p[:, 1, :] *= 1j
        return psi
    if g == "swap":
        a, b = q
        return np.ascontiguousarray(np.swapaxes(psi.reshape([2] * n), a, b)).reshape(-1)
    if g == "id":
        return psi
    raise ValueError("unsupported gate: " + g)


def apply_gates(psi, n, gates, fuse_min=None):
    fm = FUSE_MIN if fuse_min is None else fuse_min
    psi = np.ascontiguousarray(psi, dtype=np.complex128).reshape(-1)
    scale = 1.0
    seg = []

    def flush(psi):
        nonlocal seg
        if not seg:
            return psi
        s = seg
        seg = []
        if len(s) >= fm and n <= 30:
            return _fuse(psi, n, s)
        for g, q in s:
            _single(psi, n, g, q)
        return psi
    for g, q in gates:
        q = tuple(q)
        if g == "h":
            psi = flush(psi)
            _h(psi, n, q[0])
            scale *= S2
        elif g in ("ccx",) or g not in _OK:
            psi = flush(psi)
            if g == "ccx":
                _multi(psi, n, g, q)
            else:
                psi = _other(psi, n, g, q)
        else:
            seg.append((g, q))
    psi = flush(psi)
    if scale != 1.0:
        psi *= scale
    return psi
