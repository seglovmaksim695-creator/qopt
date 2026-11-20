# SPDX-License-Identifier: Apache-2.0
"""dyn.gen: small arithmetic blocks (4-25 qubits) in the AND / MZ intermediate form of qoptct.crypto.andcirc.

Every block is a classical reversible map on its data wires; ancillas start and end at |0>. An erasure is the op
('MZ', t, corr): "t holds the XOR of the monomials in corr; erase it". The lowering (dyn.lower) decides per site how
the erasure is done: by an X-basis measurement with a classically controlled phase fix-up (dynamic circuit) or by
ordinary reverse computation (static circuit).

Extra ops on top of andcirc: ('ccx', a, b, t) and correction terms ('mcz', w1, ..., wk) (a monomial of any degree).
"""
import numpy as np

from qoptct.crypto.andcirc import Circ, ZERO, ONE, const_bits


class DCirc(Circ):
    def ccx(s, a, b, t):
        s.ops.append(('ccx', a, b, t))


def simulate(C, X):
    """bit-sliced classical simulation with the gadget semantics (see andcirc.simulate), plus ccx and mcz terms.
    X: {wire: uint8 array (B,)}. -> (bits[N, B], ok, why). On a batch with ALL inputs of the block `ok` proves the
    values, the ancilla clean-up and (the MZ rule: correction == erased value) the phases for every outcome."""
    B = len(next(iter(X.values())))
    v = np.zeros((C.N, B), np.uint8)
    for w, a in X.items():
        v[w] = a
    ph = np.zeros(B, np.uint8)

    def mono(c):
        r = np.ones(B, np.uint8)
        for w in c[1:]:
            r &= v[w]
        return r
    for op in C.ops:
        k = op[0]
        if k == 'x':
            v[op[1]] ^= 1
        elif k == 'cx':
            v[op[2]] ^= v[op[1]]
        elif k == 'ccx':
            v[op[3]] ^= v[op[1]] & v[op[2]]
        elif k == 'z':
            ph ^= v[op[1]]
        elif k == 'cz':
            ph ^= v[op[1]] & v[op[2]]
        elif k == 'AND':
            a, b, t = op[1:]
            if t in (a, b):
                return v, False, 'AND target = control'
            if v[t].any():
                return v, False, 'AND target %d not 0' % t
            v[t] = v[a] & v[b]
        elif k == 'MZ':
            t, corr = op[1], op[2]; cp = np.zeros(B, np.uint8)
            for c in corr:
                if t in c[1:]:
                    return v, False, 'correction touches the measured wire'
                cp ^= mono(c)
            d = cp ^ v[t]
            if d.any():
                # a constant difference 1 would only be a global phase for the measurement variant, but the
                # reverse-computation variant needs the exact value, so we require equality.
                return v, False, 'MZ %d: correction != value' % t
            v[t] = 0
    if not (ph == ph[0]).all():
        return v, False, 'phase depends on the input'
    return v, True, ''


class Block:
    """a generated block: C (DCirc), regs {name: [wires]} (data, LSB first), anc (list), fn(dict of ints)->dict,
    domain(dict)->bool (valid inputs), name, variant (generator-level structure)."""

    def __init__(s, name, C, regs, fn, domain=None, variant='', xfree=None):
        s.name = name; s.C = C; s.regs = regs; s.fn = fn; s.domain = domain or (lambda d: True)
        s.variant = variant
        s.data = [w for r in regs.values() for w in r]
        s.anc = [w for w in range(C.N) if w not in set(s.data)]
        # wires that may be put into |+> for the superposition test (default: all data wires)
        s.xfree = list(s.data) if xfree is None else list(xfree)

    @property
    def n(s):
        return s.C.N

    def inputs(s, limit=None, seed=0):
        """list of dicts {reg: int} of valid inputs (all, or a random sample of `limit`)."""
        names = list(s.regs); sizes = [len(s.regs[k]) for k in names]; tot = sum(sizes)
        rng = np.random.default_rng(seed)
        if tot <= 16:
            allv = []
            for z in range(1 << tot):
                d = {}; sh = 0
                for k, m in zip(names, sizes):
                    d[k] = (z >> sh) & ((1 << m) - 1); sh += m
                if s.domain(d):
                    allv.append(d)
            if limit is None or len(allv) <= limit:
                return allv
            idx = rng.choice(len(allv), limit, replace=False)
            return [allv[i] for i in sorted(idx)]
        out = []
        while len(out) < (limit or 64):
            d = {k: int(rng.integers(0, 1 << m)) for k, m in zip(names, sizes)}
            if s.domain(d):
                out.append(d)
        return out

    def check(s, limit=4096, seed=0):
        """classical check on all (or `limit` sampled) valid inputs. -> (ok, why, n_inputs, exhaustive)."""
        ins = s.inputs(limit, seed)
        tot = sum(len(r) for r in s.regs.values())
        nall = len(s.inputs(None)) if tot <= 16 else None
        X = {}
        for k, ws in s.regs.items():
            for j, w in enumerate(ws):
                X[w] = np.fromiter(((d[k] >> j) & 1 for d in ins), np.uint8, len(ins))
        v, ok, why = simulate(s.C, X)
        if not ok:
            return False, why, len(ins), False
        for i, d in enumerate(ins):
            e = s.fn(d)
            for k, ws in s.regs.items():
                got = sum(int(v[w, i]) << j for j, w in enumerate(ws))
                if got != e[k]:
                    return False, 'value %s: in %r got %d want %d' % (k, d, got, e[k]), len(ins), False
        for w in s.anc:
            if v[w].any():
                return False, 'ancilla %d not clean' % w, len(ins), False
        return True, '', len(ins), (nall is not None and len(ins) == nall)


# ------------------------------------------------------------------ building bricks
def cmp_into(C, A, B, out):
    """out ^= (A < B) for equal-length registers; A, B restored. len(A) ANDs, erased by MZ."""
    n = len(A)
    for a in A:
        C.x(a)
    carr = []; c = None
    for i in range(n):
        a, b = A[i], B[i]
        if c is None:
            t = C.AND(a, b)
        else:
            C.cx(c, a); C.cx(c, b); t = C.AND(a, b); C.cx(c, t)
        carr.append(t); c = t
    C.cx(c, out)
    for i in range(n - 1, -1, -1):
        a, b = A[i], B[i]; t = carr[i]; c = carr[i - 1] if i else None
        if c is None:
            C.unAND(a, b, t)
        else:
            C.cx(c, t); C.unAND(a, b, t); C.cx(c, b); C.cx(c, a)
        C.release([t])
    for a in A:
        C.x(a)


def lookup_xor(C, addr, out, table):
    """out ^= table[addr] by unary iteration; 2^c - 2 ANDs (c = len(addr)), each erased by MZ."""
    c = len(addr)

    def leaf(ctrl, base):
        v = table[base]
        for j, o in enumerate(out):
            if (v >> j) & 1:
                C.cx(ctrl, o)

    def rec(ctrl, k, base):
        if k == 0:
            leaf(ctrl, base); return
        a = addr[k - 1]
        C.x(a); t = C.AND(ctrl, a); C.x(a)
        rec(t, k - 1, base)
        C.cx(ctrl, t)
        rec(t, k - 1, base | (1 << (k - 1)))
        C.unAND(ctrl, a, t); C.release([t])
    top = addr[c - 1]
    C.x(top); rec(top, c - 1, 0); C.x(top)
    rec(top, c - 1, 1 << (c - 1))


def anf_terms(table, c, j, addr):
    """ANF monomials (as correction terms on the address wires) of bit j of table[.]; (terms, const)."""
    f = np.array([(table[x] >> j) & 1 for x in range(1 << c)], np.uint8)
    for i in range(c):                                     # Moebius transform
        step = 1 << i
        f = f.reshape(-1, 2, step); f[:, 1, :] ^= f[:, 0, :]; f = f.reshape(-1)
    terms = []
    for m in range(1, 1 << c):
        if f[m]:
            ws = [addr[i] for i in range(c) if (m >> i) & 1]
            terms.append(('z', ws[0]) if len(ws) == 1 else ('cz', ws[0], ws[1]) if len(ws) == 2 else ('mcz',) + tuple(ws))
    return terms, int(f[0])


# ------------------------------------------------------------------ blocks
def adder(n):
    """Gidney adder b += a (mod 2^n): 3n-1 qubits, n-1 ANDs."""
    C = DCirc(); a = C.new(n); b = C.new(n)
    C.add(b, a)
    M = (1 << n) - 1
    return Block('adder%d' % n, C, dict(a=a, b=b), lambda d: dict(a=d['a'], b=(d['a'] + d['b']) & M))


def adder_cuccaro(n):
    """Cuccaro ripple adder b += a (mod 2^n) on Toffolis: 2n+1 qubits, 2(n-1) Toffolis, no erasures (baseline)."""
    C = DCirc(); a = C.new(n); b = C.new(n); z = C.new(1)[0]
    c = z
    for i in range(n - 1):
        C.cx(a[i], b[i]); C.cx(a[i], c); C.ccx(c, b[i], a[i]); c = a[i]
    C.cx(a[n - 1], b[n - 1]); C.cx(c, b[n - 1])
    for i in range(n - 2, -1, -1):
        c = z if i == 0 else a[i - 1]
        C.ccx(c, b[i], a[i]); C.cx(a[i], c); C.cx(c, b[i])
    M = (1 << n) - 1
    return Block('adder%d' % n, C, dict(a=a, b=b), lambda d: dict(a=d['a'], b=(d['a'] + d['b']) & M),
                 variant='toffoli')


def comparator(n):
    """out ^= (a < b): 3n+1 qubits, n ANDs."""
    C = DCirc(); a = C.new(n); b = C.new(n); o = C.new(1)
    cmp_into(C, a, b, o[0])
    return Block('cmp%d' % n, C, dict(a=a, b=b, o=o),
                 lambda d: dict(a=d['a'], b=d['b'], o=d['o'] ^ int(d['a'] < d['b'])))


def default_modulus(n):
    return {2: 3, 3: 7, 4: 13, 5: 29, 6: 61, 7: 127}[n]


def modadd(n, p=None):
    """b := (a + b) mod p for a, b < p (2^(n-1) < p < 2^n): add, subtract p, add p back under a flag, erase the
    flag by a comparison."""
    p = p or default_modulus(n)
    C = DCirc(); a = C.new(n); b = C.new(n); bn = C.new(1)[0]; f = C.new(1)[0]
    acc = b + [bn]
    C.add(acc, a)
    C.sub(acc, const_bits(p, n + 1), bout=f)
    C.add(acc, [f if (p >> i) & 1 else ZERO for i in range(n + 1)])
    C.x(f); cmp_into(C, b, a, f)                           # f ^= 1 ^ (r < a)
    return Block('modadd%d_p%d' % (n, p), C, dict(a=a, b=b),
                 lambda d: dict(a=d['a'], b=(d['a'] + d['b']) % p),
                 domain=lambda d: d['a'] < p and d['b'] < p,
                 xfree=a[:n - 1] + b[:n - 1])


def rand_table(c, w, seed=0):
    rng = np.random.default_rng(1000 * c + 10 * w + seed)
    return [int(x) for x in rng.integers(0, 1 << w, 1 << c)]


def lookup(c, w, table=None, seed=0):
    """out ^= T[addr]: c + w + (c-1) qubits, 2^c - 2 ANDs."""
    T = table or rand_table(c, w, seed)
    C = DCirc(); addr = C.new(c); out = C.new(w)
    lookup_xor(C, addr, out, T)
    return Block('lookup%dx%d' % (c, w), C, dict(addr=addr, out=out),
                 lambda d: dict(addr=d['addr'], out=d['out'] ^ T[d['addr']]))


def lookup_add(c, w, table=None, seed=0, out_erase='meas'):
    """acc += T[addr] (mod 2^w) (one step of windowed arithmetic): look T[addr] up into a temporary register, add it,
    erase the temporary. out_erase='meas': each temporary bit is erased by ONE MZ whose correction is the phase
    polynomial of that table bit (measurement-based unlookup); 'unitary': the lookup is run a second time."""
    T = table or rand_table(c, w, seed)
    C = DCirc(); addr = C.new(c); acc = C.new(w); tmp = C.new(w)
    lookup_xor(C, addr, tmp, T)
    C.add(acc, tmp)
    if out_erase == 'meas':
        for j in range(w):
            terms, const = anf_terms(T, c, j, addr)
            if const:
                C.x(tmp[j])
            C.MZ(tmp[j], terms)
    else:
        lookup_xor(C, addr, tmp, T)
    C.release(tmp)
    M = (1 << w) - 1
    return Block('lookadd%dx%d' % (c, w), C, dict(addr=addr, acc=acc),
                 lambda d: dict(addr=d['addr'], acc=(d['acc'] + T[d['addr']]) & M),
                 variant='out-' + out_erase)


def multiplier(n):
    """acc += x * y (mod 2^(2n)), schoolbook: for every bit of x a row of n ANDs, an addition, the row erased."""
    C = DCirc(); x = C.new(n); y = C.new(n); acc = C.new(2 * n)
    for i in range(n):
        row = [C.AND(x[i], y[j]) for j in range(n)]
        C.add(acc[i:], row[:2 * n - i])
        for j in range(n - 1, -1, -1):
            C.unAND(x[i], y[j], row[j])
        C.release(row)
    M = (1 << (2 * n)) - 1
    return Block('mult%d' % n, C, dict(x=x, y=y, acc=acc),
                 lambda d: dict(x=d['x'], y=d['y'], acc=(d['acc'] + d['x'] * d['y']) & M))


FAMILIES = {
    'adder': [2, 3, 4, 5, 6, 7, 8],
    'cmp': [1, 2, 3, 4, 5, 6, 7, 8],
    'modadd': [2, 3, 4, 5],
    'lookup': [2, 3, 4],
    'lookadd': [2, 3, 4],
    'mult': [2, 3],
}


def build(family, size, **kw):
    """Blocks of a family: the list of generator-level structural variants (most families have one)."""
    if family == 'adder':
        return [adder(size), adder_cuccaro(size)]
    if family == 'cmp':
        return [comparator(size)]
    if family == 'modadd':
        return [modadd(size, kw.get('p'))]
    if family == 'lookadd':
        w = kw.get('w') or size
        return [lookup_add(size, w, out_erase='meas'), lookup_add(size, w, out_erase='unitary')]
    if family == 'lookup':
        return [lookup(size, kw.get('w') or size)]
    if family == 'mult':
        return [multiplier(size)]
    raise ValueError('family %r' % family)
