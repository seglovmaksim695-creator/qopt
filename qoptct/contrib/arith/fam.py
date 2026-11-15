# SPDX-License-Identifier: Apache-2.0
"""Arithmetic families: recognition from the truth table (probing + exhaustive confirmation against the definition
of the function) and synthesis on logical ANDs with measurement-based uncomputation.

Families
  add      B <- A + B (mod 2^m), optionally z ^= carry-out           (Gidney 2018: m-1 or m ANDs)
  oopadd   A, B unchanged; clean wires receive sum / carry / propagate bits of A + B   (one AND per carry)
  cmp      t ^= pol ^ [X > Y]  ^  sum_p p * [X mod 2^m_p > Y mod 2^m_p]   (ripple comparator, one AND per bit)
  modadd   B <- ~((A + B + c) mod 2^n), c = [A + B + 1 >= 2^n], e ^= [B == 2^n - 1]   (addition mod 2^n - 1)
Wire roles are found by probing the table with one- and two-bit inputs; a hypothesis is accepted only if the whole
truth table equals the definition of the function. Names of wires or files are never used.
"""
import numpy as np

from .ir import Builder, BuildError


# ---------------------------------------------------------------- helpers on the table
def flips(T, on):
    """set of wires whose output differs from the identity on the input with variables ``on`` set."""
    return T.outset(T.assign(on)) ^ set(on)


def delta(T, on, base=None):
    if base is None:
        base = flips(T, [])
    return flips(T, on) ^ base


def confirm(T, expect):
    """True iff every wire's output equals expect[w] (default: identity / zero)."""
    z = T.zeros()
    for w in range(T.n):
        e = expect.get(w)
        if e is None:
            e = T.inp.get(w, z)
        if not T.eq(T.out[w], e):
            return False
    return True


def ripple(T, X, Y, cin=None):
    """bit-sliced ripple addition: -> (sums, carries) with carries[i] = carry into bit i, carries[n] = carry-out."""
    c = T.zeros() if cin is None else cin
    S, C = [], [c]
    for x, y in zip(X, Y):
        S.append(x ^ y ^ c)
        c = (x & y) | (x & c) | (y & c)
        C.append(c)
    return S, C


def gt_chain(T, X, Y):
    """g[m] = [X mod 2^m > Y mod 2^m], m = 0..n."""
    g = [T.zeros()]
    for x, y in zip(X, Y):
        g.append(T.mask((x & ~y) | (~(x ^ y) & g[-1])))
    return g


def flipvar(T, a, j):
    """words of f(x xor e_j) given words of f."""
    if j >= 6:
        blk = 1 << (j - 6)
        b = a.reshape(-1, 2, blk)
        return b[:, ::-1, :].reshape(-1).copy()
    s = np.uint64(1 << j)
    m = np.uint64(_LOWMASK[j])
    return T.mask(((a & m) >> s) | ((a << s) & m))


_LOWMASK = [0xAAAAAAAAAAAAAAAA, 0xCCCCCCCCCCCCCCCC, 0xF0F0F0F0F0F0F0F0, 0xFF00FF00FF00FF00,
            0xFFFF0000FFFF0000, 0xFFFFFFFF00000000]


def _chain_from_pairs(T, pairs, base, cyclic_ok=False):
    """pairs: dict a -> b (A_i -> B_i). Probe {a, b}: the carry goes to the next pair's B wire.
    Returns (ordered list of (a, b), top successor wire or None, cyclic flag) or None."""
    bset = {b: a for a, b in pairs.items()}
    succ = {}
    for a, b in pairs.items():
        d = delta(T, [a, b], base)
        d.discard(b)
        d.discard(a)
        if len(d) > 1:
            return None
        succ[a] = next(iter(d)) if d else None
    heads = set(pairs) - {bset[s] for s in succ.values() if s in bset}
    if len(heads) == 1:
        a = heads.pop()
        order = []
        seen = set()
        while True:
            if a in seen:
                return None
            seen.add(a)
            order.append((a, pairs[a]))
            s = succ[a]
            if s in bset:
                a = bset[s]
            else:
                break
        if len(order) != len(pairs):
            return None
        return order, succ[order[-1][0]], False
    if not heads and cyclic_ok:
        a0 = min(pairs)
        order = []
        a = a0
        while True:
            order.append((a, pairs[a]))
            s = succ[a]
            if s not in bset:
                return None
            a = bset[s]
            if a == a0:
                break
        if len(order) != len(pairs):
            return None
        return order, None, True
    return None


# ================================================================ in-place adder
def rec_add(T):
    base = flips(T, [])
    if base:
        return []
    pairs = {}
    for u in T.vars:
        d = delta(T, [u], base)
        if not d:
            continue
        if len(d) != 1:
            return []
        v = d.pop()
        if v not in T.inp or v == u:
            return []
        pairs[u] = v
    if not pairs or len(set(pairs.values())) != len(pairs) or set(pairs) & set(pairs.values()):
        return []
    ch = _chain_from_pairs(T, pairs, base)
    if ch is None:
        return []
    order, top, _ = ch
    if top is not None and top in pairs.values():
        return []
    A = [a for a, b in order]
    B = [b for a, b in order]
    z = top
    if z is not None and z not in T.inp:
        return []   # carry into a clean wire: handled by oopadd
    X = [T.inp[w] for w in A]
    Y = [T.inp[w] for w in B]
    S, C = ripple(T, X, Y)
    exp = {b: s for b, s in zip(B, S)}
    if z is not None:
        exp[z] = T.inp[z] ^ C[-1]
    if not confirm(T, exp):
        return []
    return [dict(kind='add', A=A, B=B, z=z)]


def gen_add(T, h, last='and'):
    """Gidney in-place adder B <- A + B (z ^= carry-out).  last='ccx': the top carry (carry-out into z, or carry
    c_{m-1} straight into B[m-1] when there is no z) is made by a unitary Toffoli: one ancilla less, +3 T."""
    A, B, z = h['A'], h['B'], h['z']
    m = len(A)
    bd = Builder(T)
    for w in h.get('wrap', ()):
        bd.x(w)
    _add_body(bd, A, B, z, m, last)
    for w in h.get('wrap', ()):
        bd.x(w)
    return bd


def _add_body(bd, A, B, z, m, last):
    if m == 1:
        if z is not None:
            if last == 'ccx':
                bd.ccx(A[0], B[0], z)
            else:
                t = bd.AND(A[0], B[0])
                bd.cx(t, z)
                bd.erase_and(t, A[0], B[0])
        bd.cx(A[0], B[0])
        return
    # number of carries held in ancillas: c_1 .. c_top
    top = m - 1
    if last == 'ccx' and z is None:
        top = m - 2
    c = {}
    if top >= 1:
        c[1] = bd.AND(A[0], B[0])
    for i in range(1, top):
        bd.cx(c[i], A[i])
        bd.cx(c[i], B[i])
        c[i + 1] = bd.AND(A[i], B[i])
        bd.cx(c[i], c[i + 1])
    if z is None and last == 'ccx':
        i = m - 2                         # c_{m-1} = c_i ^ A_i' B_i'  goes straight into B[m-1]
        if i == 0:
            bd.ccx(A[0], B[0], B[m - 1])
        else:
            bd.cx(c[i], A[i])
            bd.cx(c[i], B[i])
            bd.ccx(A[i], B[i], B[m - 1])
            bd.cx(c[i], B[m - 1])
        bd.cx(A[m - 1], B[m - 1])
        if i >= 1:
            bd.cx(c[i], A[i])
            bd.cx(A[i], B[i])
        else:
            bd.cx(A[0], B[0])
        lo = i - 1
    else:
        i = m - 1
        if z is not None:
            bd.cx(c[i], A[i])
            bd.cx(c[i], B[i])
            if last == 'ccx':
                bd.ccx(A[i], B[i], z)
                bd.cx(c[i], z)
            else:
                t = bd.AND(A[i], B[i])
                bd.cx(c[i], z)
                bd.cx(t, z)
                bd.erase_and(t, A[i], B[i])
            bd.cx(c[i], A[i])
            bd.cx(A[i], B[i])
        else:
            bd.cx(A[i], B[i])
            bd.cx(c[i], B[i])
        lo = m - 2
    for i in range(lo, 0, -1):
        # A_i' = a_i ^ c_i, B_i' = b_i ^ c_i ; c_{i+1} = c_i ^ A_i' B_i'
        bd.erase_and(c[i + 1], A[i], B[i], 'ab', [c[i]])
        bd.cx(c[i], A[i])
        bd.cx(A[i], B[i])
    if lo >= 0 and 1 in c:
        bd.erase_and(c[1], A[0], B[0])
        bd.cx(A[0], B[0])
    return


def gen_add_hybrid(T, h, k):
    """In-place adder with only k ancillas (1 <= k <= m-2): carries c_1..c_k as logical ANDs in ancillas, bits
    k..m-2 as Cuccaro MAJ/UMA in place (relative-phase Toffoli pairs, 4 T each), the top carry by a unitary
    Toffoli into z.  T = 8(m-1) - 4k (+7 with z), i.e. the pebbling-optimal 8n-4k-1 of lab agent b6 with z."""
    A, B, z = h['A'], h['B'], h['z']
    m = len(A)
    if not 1 <= k <= m - 2:
        raise BuildError('hybrid needs 1 <= k <= m-2')
    bd = Builder(T)
    for w in h.get('wrap', ()):
        bd.x(w)
    c = {1: bd.AND(A[0], B[0])}
    for i in range(1, k):
        bd.cx(c[i], A[i])
        bd.cx(c[i], B[i])
        c[i + 1] = bd.AND(A[i], B[i])
        bd.cx(c[i], c[i + 1])
    cw = {k: c[k]}                      # wire holding c_i (i >= k)
    for i in range(k, m - 1):           # MAJ: A_i <- c_{i+1}
        cwi = cw[i]
        bd.cx(A[i], B[i])
        bd.cx(A[i], cwi)
        bd.rtof(cwi, B[i], A[i])
        cw[i + 1] = A[i]
    top = cw[m - 1]
    i = m - 1
    if z is not None:
        bd.cx(top, A[i])
        bd.cx(top, B[i])
        bd.ccx(A[i], B[i], z)
        bd.cx(top, z)
        bd.cx(top, A[i])
        bd.cx(A[i], B[i])
    else:
        bd.cx(A[i], B[i])
        bd.cx(top, B[i])
    for i in range(m - 2, k - 1, -1):   # UMA: A_i restored, c wire restored, B_i <- s_i
        cwi = cw[i]
        bd.rtof(cwi, B[i], A[i], dagger=True)
        bd.cx(A[i], cwi)
        bd.cx(cwi, B[i])
    for i in range(k - 1, 0, -1):
        bd.erase_and(c[i + 1], A[i], B[i], 'ab', [c[i]])
        bd.cx(c[i], A[i])
        bd.cx(A[i], B[i])
    bd.erase_and(c[1], A[0], B[0])
    bd.cx(A[0], B[0])
    for w in h.get('wrap', ()):
        bd.x(w)
    return bd


def variants_add(T, h):
    out = [('gidney', lambda: gen_add(T, h, 'and'))]
    if len(h['A']) >= 2 or h['z'] is not None:
        out.append(('gidney+toffoli-top', lambda: gen_add(T, h, 'ccx')))
    m = len(h['A'])
    for k in range(1, m - 1):
        out.append(('hybrid-k%d' % k, (lambda k=k: gen_add_hybrid(T, h, k))))
    return out


# ================================================================ subtraction  B <- B - A  (= ~(~B + A))
def rec_sub(T):
    """conjugate the function by X on the register that single-bit probes flip, then recognise an adder."""
    from .tt import Table
    base = flips(T, [])
    if base:
        return []
    reg = set()
    for u in T.vars:
        reg |= delta(T, [u], base)
    reg = sorted(w for w in reg if w in T.inp)
    if not reg:
        return []
    xs = [('x', (w,)) for w in reg]
    T2 = Table(T.names, xs + list(T.gates) + xs, T.vars, T.plus, T.clean)
    out = []
    for h in rec_add(T2):
        if set(h['B']) | ({h['z']} if h['z'] is not None else set()) >= set(reg) and \
                set(reg) >= set(h['B']):
            h = dict(h, kind='sub', wrap=list(h['B']))
            if h['z'] is not None:
                if h['z'] in reg:
                    h['wrap'] = h['wrap'] + [h['z']]
            if confirm_sub(T, h):
                out.append(h)
    return out


def confirm_sub(T, h):
    A, B, z = h['A'], h['B'], h['z']
    X = [T.inp[w] for w in A]
    Y = [T.mask(~T.inp[w]) if w in h['wrap'] else T.inp[w] for w in B]
    S, C = ripple(T, X, Y)
    exp = {b: (T.mask(~s) if b in h['wrap'] else s) for b, s in zip(B, S)}
    if z is not None:
        v = T.inp[z]
        if z in h['wrap']:
            v = T.mask(~v)
        v = v ^ C[-1]
        if z in h['wrap']:
            v = T.mask(~v)
        exp[z] = v
    return confirm(T, exp)


# ================================================================ out-of-place adder
def rec_oopadd(T):
    ch = T.changed()
    if not ch or any(w in T.inp for w in ch):
        return []
    base = flips(T, [])
    if base:
        return []
    groups = {}
    for u in T.vars:
        d = frozenset(delta(T, [u], base))
        if d:
            groups.setdefault(d, []).append(u)
    pairs = []
    for d, us in groups.items():
        if len(us) != 2:
            return []
        pairs.append(tuple(us))
    n = len(pairs)
    if n < 1:
        return []
    allA = [p[0] for p in pairs]
    score = []
    for a, b in pairs:
        on = [x for x in allA if x != a] + [a]
        d = flips(T, on + [b]) ^ flips(T, on) ^ {b}
        score.append(len(d))
    if len(set(score)) != n:
        return []
    order = [p for _, p in sorted(zip(score, pairs), key=lambda t: -t[0])]
    A = [p[0] for p in order]
    B = [p[1] for p in order]
    X = [T.inp[w] for w in A]
    Y = [T.inp[w] for w in B]
    S, C = ripple(T, X, Y)
    sig = {}
    for i in range(n):
        sig[('s', i)] = S[i]
        sig[('p', i)] = X[i] ^ Y[i]
    for i in range(1, n + 1):
        sig[('c', i)] = C[i]
    want = {}
    for w in ch:
        hit = None
        for key, v in sig.items():
            if T.eq(v, T.out[w]):
                hit = key
                break
        if hit is None:
            return []
        want[w] = hit
    return [dict(kind='oopadd', A=A, B=B, want=want)]


def gen_oopadd(T, h):
    A, B, want = h['A'], h['B'], h['want']
    n = len(A)
    bd = Builder(T, reserved=list(want))
    need = 0
    for (k, i) in want.values():
        if k == 'c':
            need = max(need, i)
        elif k == 's' and i >= 1:
            need = max(need, i)
    holder = {}
    conv = {}            # holder wire -> position i to convert c_i -> s_i
    used = set()
    for w, (k, i) in sorted(want.items()):
        if k == 'c' and i not in holder:
            holder[i] = w
            used.add(w)
    for w, (k, i) in sorted(want.items()):
        if k == 's' and i >= 1 and i not in holder and w not in used:
            holder[i] = w
            conv[w] = i
            used.add(w)
    temps = []
    for i in range(1, need + 1):
        if i not in holder:
            holder[i] = bd.alloc()
            temps.append(i)
    # forward: carries c_1..c_need
    if need >= 1:
        bd.AND(A[0], B[0], holder[1])
    for i in range(1, need):
        ci = holder[i]
        bd.cx(ci, A[i])
        bd.cx(ci, B[i])
        bd.AND(A[i], B[i], holder[i + 1])
        bd.cx(ci, holder[i + 1])
        bd.cx(ci, A[i])
        bd.cx(ci, B[i])
    # erase temporary carries, highest first (c_i still in carry form below)
    for i in sorted(temps, reverse=True):
        t = holder[i]
        if i == 1:
            bd.erase_and(t, A[0], B[0])
            continue
        cj = holder[i - 1]
        bd.cx(cj, A[i - 1])
        bd.cx(cj, B[i - 1])
        bd.cx(cj, t)
        bd.erase_and(t, A[i - 1], B[i - 1])
        bd.cx(cj, A[i - 1])
        bd.cx(cj, B[i - 1])
    # convert c_i -> s_i
    for w, i in conv.items():
        bd.cx(A[i], w)
        bd.cx(B[i], w)
    # remaining outputs: copies and p / s_0
    for w, (k, i) in sorted(want.items()):
        if w in used:
            continue
        if k == 'p' or (k == 's' and i == 0):
            bd.cx(A[i], w)
            bd.cx(B[i], w)
        else:
            src = None
            for u in used:
                if want[u] == (k, i):
                    src = u
                    break
            if src is None:
                raise BuildError('no source for %s' % ((k, i),))
            bd.cx(src, w)
    return bd


def variants_oopadd(T, h):
    return [('oop-ripple', lambda: gen_oopadd(T, h))]


# ================================================================ comparator
def rec_cmp(T):
    ch = T.changed()
    if len(ch) != 1:
        return []
    t = ch[0]
    f = T.out[t] ^ T.inp[t] if t in T.inp else T.out[t]
    plus = [p for p in T.plus if p != t]
    terms = []
    f0 = f
    for p in plus:
        j = T.pos[p]
        fl = flipvar(T, f0, j)
        d = f0 ^ fl
        if not T.eq(flipvar(T, d, j), d):
            return []
        P = T.inp[p]
        f0 = T.mask((f0 & ~P) | (fl & P))
        if not T.is_zero(d):
            terms.append((p, d))
    rest = [u for u in T.vars if u != t and u not in plus]

    def val(on):
        x = T.assign(on)
        return int((int(f0[x >> 6]) >> (x & 63)) & 1)
    pol = val([])
    Xs = [u for u in rest if val([u]) != pol]
    Ys = [u for u in rest if val([u]) == pol]
    n = len(Xs)
    if n == 0 or len(Ys) != n:
        return []
    rx = {u: sum(1 for y in Ys if val([u, y]) != pol) for u in Xs}
    ry = {y: n - 1 - sum(1 for u in Xs if val([u, y]) != pol) for y in Ys}
    if sorted(rx.values()) != list(range(n)) or sorted(ry.values()) != list(range(n)):
        return []
    X = sorted(Xs, key=lambda u: rx[u])
    Y = sorted(Ys, key=lambda u: ry[u])
    g = gt_chain(T, [T.inp[w] for w in X], [T.inp[w] for w in Y])
    want0 = g[n] if not pol else T.mask(~g[n])
    if not T.eq(want0, f0):
        return []
    hooks = []
    for p, d in terms:
        m = next((m for m in range(1, n) if T.eq(g[m], d)), None)
        if m is None:
            return []
        hooks.append((p, m))
    return [dict(kind='cmp', X=X, Y=Y, t=t, pol=pol, hooks=hooks, clean_t=t not in T.inp)]


def gen_cmp(T, h):
    X, Y, tgt, pol, hooks = h['X'], h['Y'], h['t'], h['pol'], h['hooks']
    n = len(X)
    clean_t = h['clean_t']
    bd = Builder(T, reserved=[tgt])
    t = [None] + [bd.alloc() for _ in range(n - 1)]
    last = tgt if clean_t else bd.alloc()
    if n == 1:
        bd.AND(X[0], Y[0], last, 'a!b')
    else:
        bd.AND(X[0], Y[0], t[1], 'a!b')
        for i in range(1, n):
            bd.cx(t[i], X[i])
            bd.cx(t[i], Y[i])
            dst = t[i + 1] if i < n - 1 else last
            bd.AND(X[i], Y[i], dst, 'a!b')
            bd.cx(t[i], dst)
    if not clean_t:
        bd.cx(last, tgt)
    for p, m in hooks:
        w = bd.AND(p, t[m])
        bd.cx(w, tgt)
        bd.erase_and(w, p, t[m])
    if not clean_t:
        i = n - 1
        lin = [t[i]] if i >= 1 else []
        bd.erase_and(last, X[i], Y[i], 'a!b', lin)
    for i in range(n - 1, 0, -1):
        bd.cx(t[i], X[i])
        bd.cx(t[i], Y[i])
        j = i - 1
        bd.erase_and(t[i], X[j], Y[j], 'a!b', [t[j]] if j >= 1 else [])
    if pol:
        bd.x(tgt)
    return bd


def variants_cmp(T, h):
    return [('ripple-cmp', lambda: gen_cmp(T, h))]


# ================================================================ addition mod 2^n - 1 (end-around carry)
def rec_modadd(T):
    base = flips(T, [])
    pairs = {}
    for u in T.vars:
        d = delta(T, [u], base)
        if not d:
            continue
        if len(d) != 1:
            return []
        v = d.pop()
        if v not in T.inp or v == u:
            return []
        pairs[u] = v
    if len(pairs) < 2 or set(pairs) & set(pairs.values()):
        return []
    ch = _chain_from_pairs(T, pairs, base, cyclic_ok=True)
    if ch is None or not ch[2]:
        return []
    order = ch[0]
    n = len(order)
    others = [w for w in T.changed() if w not in pairs.values()]
    if len(others) != 1 or others[0] in T.inp:
        return []
    e = others[0]
    for r in range(n):
        o = order[r:] + order[:r]
        A = [a for a, b in o]
        B = [b for a, b in o]
        X = [T.inp[w] for w in A]
        Y = [T.inp[w] for w in B]
        _, C1 = ripple(T, X, Y, T.ones())
        c = C1[-1]
        S, _ = ripple(T, X, Y, c)
        exp = {b: T.mask(~s) for b, s in zip(B, S)}
        allb = T.ones()
        for y in Y:
            allb = allb & y
        exp[e] = allb
        if confirm(T, exp):
            return [dict(kind='modadd', A=A, B=B, e=e)]
    return []


def _fwd_bit(bd, a, b, tin, tout, pol='ab'):
    bd.cx(tin, a)
    bd.cx(tin, b)
    bd.AND(a, b, tout, pol)
    bd.cx(tin, tout)


def _and_tree(bd, xs, out):
    k = len(xs)
    if k == 1:
        bd.cx(xs[0], out)
        return
    if k == 2:
        bd.AND(xs[0], xs[1], out)
        return
    tmp = [bd.alloc() for _ in range(k - 2)]
    bd.AND(xs[0], xs[1], tmp[0])
    for i in range(2, k - 1):
        bd.AND(tmp[i - 2], xs[i], tmp[i - 1])
    bd.AND(tmp[k - 3], xs[k - 1], out)
    for i in range(k - 3, 0, -1):
        bd.erase_and(tmp[i], tmp[i - 1], xs[i + 1])
    bd.erase_and(tmp[0], xs[0], xs[1])


def gen_modadd(T, h):
    """b <- ~((a+b+c) mod 2^n), c = carry(a+b+1); e = [b == 1...1]  (a5 construction, T = 4(4n-3))."""
    a, b, e = h['A'], h['B'], h['e']
    n = len(a)
    bd = Builder(T, reserved=[e])
    C = bd.alloc()
    # S0: e = AND(b)
    _and_tree(bd, list(b), e)
    # S1: C = carry-out(a + b + 1)
    t = [bd.alloc() for _ in range(n - 1)]
    bd.AND(a[0], b[0], t[0], 'or')
    for i in range(1, n):
        dst = t[i] if i < n - 1 else C
        _fwd_bit(bd, a[i], b[i], t[i - 1], dst)
    for i in range(n - 1, 0, -1):
        bd.cx(t[i - 1], a[i])
        bd.cx(t[i - 1], b[i])
        j = i - 1
        if j == 0:
            bd.erase_and(t[0], a[0], b[0], 'or')
        else:
            bd.erase_and(t[j], a[j], b[j], 'ab', [t[j - 1]])
    # S2: b <- a + b + C
    tt = [C] + t
    for i in range(0, n - 1):
        _fwd_bit(bd, a[i], b[i], tt[i], tt[i + 1])
    k = n - 1
    bd.cx(tt[k], a[k])
    bd.cx(tt[k], b[k])
    bd.cx(a[k], b[k])
    bd.cx(tt[k], b[k])
    bd.cx(tt[k], a[k])
    for i in range(n - 2, -1, -1):
        bd.erase_and(tt[i + 1], a[i], b[i], 'ab', [tt[i]])
        bd.cx(a[i], b[i])
        bd.cx(tt[i], b[i])
        bd.cx(tt[i], a[i])
    # S3: erase C = carry-out(a + ~s + e) through a phase chain
    s = b
    for w in s:
        bd.x(w)
    t = [bd.alloc() for _ in range(n - 1)]
    tt = [e] + t
    for i in range(0, n - 1):
        _fwd_bit(bd, a[i], s[i], tt[i], tt[i + 1])
    k = n - 1
    bd.cx(tt[k], a[k])
    bd.cx(tt[k], s[k])
    bd.erase(C, [('z', (tt[k],)), ('cz', (a[k], s[k]))])
    bd.cx(tt[k], a[k])
    bd.cx(tt[k], s[k])
    for i in range(n - 2, -1, -1):
        bd.erase_and(tt[i + 1], a[i], s[i], 'ab', [tt[i]])
        bd.cx(tt[i], a[i])
        bd.cx(tt[i], s[i])
    for w in s:
        bd.x(w)
    for w in b:
        bd.x(w)
    return bd


def variants_modadd(T, h):
    return [('mod-2^n-1', lambda: gen_modadd(T, h))]


FAMILIES = [
    ('add', rec_add, variants_add),
    ('sub', rec_sub, variants_add),
    ('oopadd', rec_oopadd, variants_oopadd),
    ('cmp', rec_cmp, variants_cmp),
    ('modadd', rec_modadd, variants_modadd),
]
