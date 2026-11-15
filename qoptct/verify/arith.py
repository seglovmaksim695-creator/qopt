# SPDX-License-Identifier: Apache-2.0
"""qoptct.verify.arith: fast checks of arithmetic circuits (Clifford+T, logical AND with measurement-based uncomputation) at 64-2048+ bits.

numpy only. Circuits in OpenQASM 3 (qoptct subset): h x z s sdg t tdg cx cz ccx ccz swap reset, `h q; measure q` (X measurement),
`c[k] = measure q[i]` (Z measurement), `if (c[a] ^ c[b] && c[d]) {...}`.

METHOD. For a basis input x and GIVEN measurement outcomes m, the circuit state is a sum over paths (Amy 2018):
    psi = 2^{-nH/2} sum_y w^{F(y)} |W y + e>,   w = e^{i pi/4},
F is a Z8 polynomial in parities of the variables y (T/S/Z linear terms, CZ/CCZ quadratic/cubic via parities).
A batch of inputs (B of them) uses ONE graph: the structure (W, masks) does not depend on x; only wire constants and
coefficient signs do -> factor tables have a leading batch axis.
The sum over y is exact variable elimination (factor graph, min-degree; the width is small for arithmetic).
Checked (for every branch m and every x of the batch):
  (a) the amplitude A at the expected point (classical function f(x), ancillas = 0);
  (b) the FULL branch norm ||psi||^2 = sum_d |A(d)|^2 (doubled ket-bra network) -> |A|^2 == ||psi||^2  <=>  branch = exactly |f(x)>
      (no assumption that all outcomes are equally likely);
  (c) the phase A/|A| and the modulus |A| do not depend on x (one global phase c_m per branch).
If ALL x are enumerated (n <= 10) this is a proof for the given branch m; otherwise (random x) it is a sampled check over x.
Branches m are random (m_trials of them). check_proof: a miter network over ALL x at once (x are summed variables, the reference is a Boolean network) - a proof for the given branch m.
check_sim / check(method='sim'): fast simulator "classical wires + small quantum register" (speed ~ gate count, 2048 bits in seconds).

Ported from the qoptct lab (agent a8). v0.4 fixes: in the clean model check_sim compared only the input wires
(outputs on clean wires were not checked, so wrong circuits could pass); a measured reference wire is now required to
hold its expected value; single `&` and `!` in conditions are rejected instead of being misread.
Branches m are SAMPLED: use qoptct.verify.prove for a proof over all inputs and all outcomes.
"""
import re, sys, math, time, heapq, random, argparse, json
from collections import defaultdict
import numpy as np

OM = np.exp(1j * math.pi * np.arange(8) / 4)


# ============================================================ parsing
from ._parse import split_if  # noqa: E402
def parse_cond(s):
    s = s.strip()
    if not re.fullmatch(r"[\s()c\[\]\d&^]+", s) or re.search(r"(?<!&)&(?!&)", s):
        raise ValueError("unsupported condition: " + s)
    groups = []
    for g in s.split("&&"):
        g = g.strip()
        while g.startswith("(") and g.endswith(")") and g.count("(") == 1:
            g = g[1:-1].strip()
        g = g.replace("(", "").replace(")", "")
        groups.append([int(x) for x in re.findall(r"c\[(\d+)\]", g)])
    return groups


def _gate(s):
    s = s.rstrip(";").strip()
    parts = s.split(None, 1)
    name = parts[0]
    qs = tuple(int(x) for x in re.findall(r"q\[(\d+)\]", parts[1])) if len(parts) > 1 else ()
    return (name, qs)


def parse_text(txt):
    """-> (N, nbits, ops). ops: ('g',name,qs) | ('h',q) | ('m',q,k) | ('if',groups,[(name,qs)..])"""
    N = int(re.search(r"qubit\[(\d+)\]", txt).group(1))
    mb = re.search(r"(?<!qu)bit\[(\d+)\]", txt)
    nb = int(mb.group(1)) if mb else 0
    ops = []
    for L in txt.splitlines():
        L = L.strip()
        if not L or L.startswith(("OPENQASM", "include", "qubit", "bit", "//")):
            continue
        m = re.match(r"c\[(\d+)\] = measure q\[(\d+)\];", L)
        if m:
            ops.append(("m", int(m.group(2)), int(m.group(1)))); continue
        sp = split_if(L)  # also nested ifs (the Qiskit-readable form of AND conditions)
        if sp is not None:
            for cond, body in sp:
                ops.append(("if", parse_cond(cond), [_gate(b) for b in body]))
            continue
        name, qs = _gate(L)
        if name == "h":
            ops.append(("h", qs[0]))
        else:
            ops.append(("g", name, qs))
    return N, nb, ops


def parse(path):
    return parse_text(open(path).read())


# ============================================================ batched path-sum builder
class Net:
    """wires[q] = (mask:int, const:uint8[B]);  lin[mask] = coefficient mod 8 (int or int16[B]) of the parity P(mask)."""

    def __init__(s, N, init, xvars=None):
        """init: (B,N) uint8 (concrete inputs) OR None + xvars {wire: variable index} (symbolic inputs, B=1)"""
        if init is None:
            init = np.zeros((1, N), np.uint8)
        B = init.shape[0]
        s.B, s.N, s.k = B, N, 0
        s.zero = np.zeros(B, np.uint8)
        s.w = [(0, init[:, q].astype(np.uint8)) for q in range(N)]
        if xvars:
            s.k = len(xvars)
            for q, v in xvars.items(): s.w[q] = (1 << v, s.zero)
        s.lin = {}; s.quad = {}
        s.gph = np.zeros(B, np.int16)
        s.nH = 0          # number of 1/sqrt2 factors
        s.nmeas = 0
        s.cons = []       # (mask, bit[B]) : P(mask)=bit
        s.alive = [True] * N
        s.defs = []       # (u, mask): auxiliary variable u = P(mask) (keeps masks short -> local factors)
        s.defc = {}

    CAP = 3

    def small(s, mask):
        if mask.bit_count() <= s.CAP: return mask
        u = s.defc.get(mask)
        if u is None:
            u = s.k; s.k += 1; s.defs.append((u, mask)); s.defc[mask] = u
        return 1 << u

    def add_lin(s, mask, c, const):
        """w^{ c * (const XOR P(mask)) },  c: int or int16[B]"""
        mask = s.small(mask)
        if isinstance(c, (int, np.integer)):
            c = int(c) % 8
            if c == 0: return
        if mask == 0:
            s.gph = (s.gph + c * const.astype(np.int16)) % 8; return
        if const.any():
            ci = const.astype(np.int16)
            s.gph = (s.gph + c * ci) % 8
            delta = c * (1 - 2 * ci)
        else:
            delta = c
        cur = s.lin.get(mask, 0)
        new = (cur + delta) % 8
        if (new == 0).all() if isinstance(new, np.ndarray) else new == 0:
            s.lin.pop(mask, None)
        else:
            s.lin[mask] = new

    def add_quad(s, wa, wb):
        """(-1)^{va*vb}"""
        (ma, ca), (mb, cb) = wa, wb
        ma = s.small(ma); mb = s.small(mb)
        if ca.any() and cb.any(): s.gph = (s.gph + 4 * (ca & cb).astype(np.int16)) % 8
        if ca.any(): s.add_lin(mb, 4 * ca.astype(np.int16), s.zero)
        if cb.any(): s.add_lin(ma, 4 * cb.astype(np.int16), s.zero)
        if ma == 0 or mb == 0: return
        if ma == mb:
            s.add_lin(ma, 4, s.zero); return
        key = (ma, mb) if ma < mb else (mb, ma)
        if key in s.quad: del s.quad[key]
        else: s.quad[key] = 1

    def gate(s, name, qs):
        w = s.w
        if name == "cx":
            c, t = qs; w[t] = (s.small(w[t][0] ^ w[c][0]), w[t][1] ^ w[c][1])
        elif name == "x":
            q = qs[0]; w[q] = (w[q][0], w[q][1] ^ 1)
        elif name in ("t", "tdg", "s", "sdg", "z"):
            c = {"t": 1, "tdg": 7, "s": 2, "sdg": 6, "z": 4}[name]; m, cn = w[qs[0]]; s.add_lin(m, c, cn)
        elif name == "cz":
            s.add_quad(w[qs[0]], w[qs[1]])
        elif name == "ccz":
            s._ccz(qs)
        elif name == "ccx":
            s.hadamard(qs[2]); s._ccz(qs); s.hadamard(qs[2])
        elif name == "swap":
            a, b = qs; w[a], w[b] = w[b], w[a]
        elif name == "reset":
            q = qs[0]
            if w[q][0] != 0 or w[q][1].any(): raise ValueError("reset of an unmeasured wire is not supported")
            s.alive[q] = True
        elif name == "id":
            pass
        else:
            raise ValueError("unsupported gate: " + name)

    def _ccz(s, qs):
        ws = [s.w[q] for q in qs]
        coef = {1: 1, 2: 1, 4: 1, 3: 7, 5: 7, 6: 7, 7: 1}
        for sub, c in coef.items():
            m = 0; cn = s.zero
            for i in range(3):
                if sub >> i & 1: m ^= ws[i][0]; cn = cn ^ ws[i][1]
            s.add_lin(m, c, cn)

    def hadamard(s, q):
        m, c = s.w[q]
        v = 1 << s.k; s.k += 1; s.nH += 1
        if m == 0 and not c.any():
            s.w[q] = (v, s.zero); return
        s.add_quad((m, c), (v, s.zero))
        s.w[q] = (v, s.zero)

    def measure_after_h(s, q, outcome):
        m, c = s.w[q]
        s.nH += 1; s.nmeas += 1
        if outcome: s.add_lin(m, 4, c)
        s.w[q] = (0, np.full(s.B, outcome, np.uint8)); s.alive[q] = False

    def measure_z(s, q, outcome):
        m, c = s.w[q]; m = s.small(m)
        s.cons.append((m, (c ^ outcome).astype(np.uint8)))   # P(m) = outcome ^ c
        s.w[q] = (0, np.full(s.B, outcome, np.uint8)); s.alive[q] = False
        s.nmeas_z = getattr(s, "nmeas_z", 0) + 1


def build(prog, init, outcomes):
    N, nb, ops = prog
    c = [0] * nb
    net = Net(N, init)
    i = 0
    while i < len(ops):
        op = ops[i]
        if op[0] == "h":
            q = op[1]; nxt = ops[i + 1] if i + 1 < len(ops) else None
            if nxt and nxt[0] == "m" and nxt[1] == q:
                k = nxt[2]; c[k] = outcomes[k]; net.measure_after_h(q, outcomes[k]); i += 2; continue
            net.hadamard(q)
        elif op[0] == "g":
            net.gate(op[1], op[2])
        elif op[0] == "if":
            if all((sum(c[b] for b in g) & 1) for g in op[1]):
                for name, qs in op[2]: net.gate(name, qs)
        elif op[0] == "m":
            k = op[2]; c[k] = outcomes[k]; net.measure_z(op[1], outcomes[k])
        i += 1
    return net, c


# ============================================================ factor graph + elimination
_PT = {}


def _bits(m):
    r = []; i = 0
    while m:
        low = m & -m; i = low.bit_length() - 1; r.append(i); m ^= low
    return r


def _pt(scope, sub):
    """parity table of the subset sub of variables along the axes scope (sorted)"""
    k = len(scope); t = np.zeros((2,) * k, np.int8)
    for v in sub:
        ax = scope.index(v); sh = [1] * k; sh[ax] = 2
        t = t ^ np.arange(2, dtype=np.int8).reshape(sh)
    return t


class TooWide(Exception):
    pass


class FG:
    def __init__(s, B, nv0, maxscope=22, nshared=0):
        s.B = B; s.nv = nv0; s.maxscope = maxscope; s.nshared = nshared
        s.fac = {}; s.nf = 0; s.v2f = defaultdict(set)
        s.scal = np.ones(1, complex); s.exp = 0.0; s.zero = False
        s.maxelems = 1 << 24      # <= 256 MB complex128 per table (limit ~1.5 GB per process)
        s.pcache = {}; s.maxw = 0; s.tag = ''; s.nlast = 0   # variables < nlast (symbolic x) are eliminated last

    def newvar(s):
        s.nv += 1; return s.nv - 1

    def add(s, scope, table):
        """scope: list of variables (any order); table: (Bt,)+(2,)*k in scope order"""
        if len(scope) == 0:
            s.scal = s.scal * table.reshape(-1); return
        order = np.argsort(scope)
        if list(order) != list(range(len(scope))):
            table = np.transpose(table, [0] + [1 + int(i) for i in order])
        sc = tuple(int(scope[i]) for i in order)
        fid = s.nf; s.nf += 1
        s.fac[fid] = (sc, table)
        for v in sc: s.v2f[v].add(fid)

    def pv(s, mask, off=0):
        """variables whose parity = P(mask) (indices shifted by off); heavy masks get an auxiliary variable"""
        vs = [(v + off if v >= s.nshared else v) for v in _bits(mask)]
        if len(vs) <= 10: return tuple(vs)
        key = (mask, off)
        if key in s.pcache: return s.pcache[key]
        z = s.newvar(); acc = [vs[0], vs[1], vs[2]]
        rest = vs[3:]
        s._xor_def(z, acc)
        while rest:
            g = rest[:3]; rest = rest[3:]
            z2 = s.newvar(); s._xor_def(z2, [z] + g); z = z2
        s.pcache[key] = (z,)
        return (z,)

    def _xor_def(s, z, vs):
        scope = [z] + list(vs)
        t = (_pt(sorted(scope), scope) == 0).astype(complex)  # z xor vs = 0
        s.add(sorted(scope), t.reshape((1,) + t.shape))

    def add_lin(s, mask, coef, off=0, conj=False):
        vs = s.pv(mask, off); scope = sorted(vs)
        par = _pt(scope, scope)
        if isinstance(coef, np.ndarray):
            val = OM[coef.astype(np.int64) % 8]
            if conj: val = val.conj()
            if (coef == coef.flat[0]).all(): val = val[:1]
            val = val.reshape((-1,) + (1,) * len(scope))
        else:
            val = OM[int(coef) % 8]
            if conj: val = np.conj(val)
            val = np.array([val]).reshape((1,) + (1,) * len(scope))
        tbl = np.where(par[None] == 1, val, 1.0 + 0j)
        s.add(scope, tbl)

    def add_quad(s, ma, mb, off=0):
        va = s.pv(ma, off); vb = s.pv(mb, off)
        scope = sorted(set(va) | set(vb))
        pa = _pt(scope, va); pb = _pt(scope, vb)
        tbl = np.where((pa & pb) == 1, -1.0 + 0j, 1.0 + 0j)
        s.add(scope, tbl.reshape((1,) + tbl.shape))

    def add_parity_eq(s, parts, bit):
        """P(sum of the masks in parts) = bit(batch). parts: [(mask, off)]"""
        vs = []
        for m, off in parts: vs += list(s.pv(m, off))
        # equal variables cancel modulo 2
        cnt = defaultdict(int)
        for v in vs: cnt[v] ^= 1
        vs = [v for v, c in cnt.items() if c]
        bit = np.asarray(bit, np.uint8)
        if not vs:
            ok = (bit == 0).astype(complex); s.scal = s.scal * ok; return
        scope = sorted(vs); par = _pt(scope, scope)
        tbl = (par[None] == bit.reshape((-1,) + (1,) * len(scope))).astype(complex)
        s.add(scope, tbl)

    # ---------------------------------------------------- elimination
    def _cost(s, v):
        u = set()
        for f in s.v2f[v]: u.update(s.fac[f][0])
        return len(u) - 1

    def eliminate(s):
        used = [v for v, fs in s.v2f.items() if fs]
        pen = lambda v: 1 if v < s.nlast else 0
        heap = [(pen(v), s._cost(v), len(s.v2f[v]), v) for v in used]
        heapq.heapify(heap)
        while heap:
            pn, c, nfac, v = heapq.heappop(heap)
            if not s.v2f.get(v): continue
            c2 = s._cost(v)
            if c2 != c or len(s.v2f[v]) != nfac:
                heapq.heappush(heap, (pn, c2, len(s.v2f[v]), v)); continue
            if c > s.maxscope: raise TooWide("width %d (%s)" % (c, s.tag))
            s.maxw = max(s.maxw, c)
            fids = sorted(s.v2f[v], key=lambda f: -len(s.fac[f][0]))
            U = sorted({u for f in fids for u in s.fac[f][0]})
            bt = max(s.fac[f][1].shape[0] for f in fids)
            if bt * (1 << len(U)) > s.maxelems:
                raise TooWide("memory: batch %d x 2^%d (%s)" % (bt, len(U), s.tag))
            pos = {u: i for i, u in enumerate(U)}
            prod = None
            for f in fids:
                sc, tb = s.fac.pop(f)
                for u in sc: s.v2f[u].discard(f)
                sh = [tb.shape[0]] + [2 if u in sc else 1 for u in U]
                tb = tb.reshape(sh)
                prod = tb if prod is None else prod * tb
            if prod.shape[1:] != (2,) * len(U):
                prod = np.broadcast_to(prod, (prod.shape[0],) + (2,) * len(U))
            res = prod.sum(axis=1 + pos[v])
            newscope = [u for u in U if u != v]
            mx = float(np.abs(res).max())
            if mx == 0.0:
                s.zero = True; return
            e = math.frexp(mx)[1]
            if abs(e) > 16:
                res = res * (2.0 ** -e); s.exp += e
            del s.v2f[v]
            if newscope:
                s.add(newscope, res)
                for u in newscope:
                    heapq.heappush(heap, (pen(u), s._cost(u), len(s.v2f[u]), u))
            else:
                s.scal = s.scal * res.reshape(-1)
        # variables without factors: the sum gives 2
        nfree = s.nv - len(used)
        s.exp += nfree
        # remaining scalar factors
        return

    def result(s):
        """(mant complex[B], exp) : value = mant * 2^exp"""
        if s.zero: return np.zeros(s.B, complex), 0.0
        mx = np.abs(s.scal).max()
        if mx == 0: return np.zeros(s.B, complex), 0.0
        e = math.frexp(mx)[1]
        return s.scal * (2.0 ** -e) * np.ones(s.B), s.exp + e


def net_factors(net, fg, off=0, conj=False):
    """adds the factors of the network to fg (variables shifted by off)"""
    for mask, coef in net.lin.items():
        fg.add_lin(mask, coef, off, conj)
    for (ma, mb) in net.quad:
        fg.add_quad(ma, mb, off)
    for m, bit in net.cons:
        fg.add_parity_eq([(m, off)], bit)
    z0 = np.zeros(net.B, np.uint8)
    for u, m in net.defs:
        fg.add_parity_eq([(m, off), (1 << u, off)], z0)
    g = OM[net.gph % 8]
    fg.scal = fg.scal * (g.conj() if conj else g)
    fg.exp -= 0.5 * net.nH


def amplitude(net, targets):
    """A(target) for the batch: targets {wire: bit[B]}; value = (mant[B], exp)"""
    fg = FG(net.B, net.k); fg.tag = 'amp'
    net_factors(net, fg)
    for q, bit in targets.items():
        m, c = net.w[q]
        fg.add_parity_eq([(m, 0)], bit ^ c)
    fg.eliminate()
    return fg.result()


def norm2(net, maxscope=22, nshared=0, nlast=None):
    """||psi||^2 = sum_d |A(d)|^2 (doubled network); value (mant[B], exp). nshared: the first nshared variables (symbolic x) are shared"""
    off = net.k - nshared
    fg = FG(net.B, net.k + off, maxscope=maxscope, nshared=nshared); fg.tag = 'norm'; fg.nlast = nshared if nlast is None else nlast
    net_factors(net, fg, 0, False)
    net_factors(net, fg, off, True)
    for q in range(net.N):
        m, c = net.w[q]
        if m: fg.add_parity_eq([(m, 0), (m, off)], np.zeros(net.B, np.uint8))
    fg.eliminate()
    return fg.result()


# ============================================================ simulator "classical wires + small quantum register"
class Sim:
    """State = amp * |bits> on classical wires (bits: (N,B) uint8, per x of the batch) and a dense tensor amp (B,2,..,2) over the wires S
    currently in superposition (AND gadgets live there for a few gates). A wire that becomes deterministic returns to bits (collapse).
    Everything is exact (dyadic amplitudes); the global phase of the classical part is ph (mod 8), the normalisation e2 (power of 2)."""
    MAXS = 16
    EPS = 1e-18

    def __init__(s, X, N):
        B = X.shape[0]
        s.B, s.N = B, N
        s.bits = np.zeros((N, B), np.uint8); s.bits[:X.shape[1]] = X.T
        s.S = []; s.amp = np.ones(B, complex); s.ph = np.zeros(B, np.int16); s.e2 = 0.0
        s.maxS = 0

    def bc(s, a):
        return a.reshape((s.B,) + (1,) * len(s.S))

    def promote(s, w):
        if len(s.S) >= s.MAXS: raise TooWide("S=%d wires in superposition" % (len(s.S) + 1))
        b = s.bits[w]; new = np.zeros(s.amp.shape + (2,), complex)
        new[..., 0] = s.amp * s.bc(b == 0); new[..., 1] = s.amp * s.bc(b == 1)
        s.amp = new; s.S.append(w); s.bits[w] = 0
        s.maxS = max(s.maxS, len(s.S))

    def collapse(s, w, force=False):
        ax = s.S.index(w) + 1
        a0 = np.take(s.amp, 0, axis=ax); a1 = np.take(s.amp, 1, axis=ax)
        n0 = (np.abs(a0) ** 2).reshape(s.B, -1).sum(1); n1 = (np.abs(a1) ** 2).reshape(s.B, -1).sum(1)
        tot = n0 + n1; th = s.EPS * np.maximum(tot, 1e-300)
        z0 = n0 <= th; z1 = n1 <= th
        if not (z0 | z1).all(): return False
        one = (z0 & ~z1)
        s.amp = np.where(s.bc_after(one, ax), a1, a0)
        s.S.remove(w); s.bits[w] = one.astype(np.uint8)
        return True

    def bc_after(s, a, ax):   # a (B,) -> shape (B,)+(1,)*(ndim-1) without the axis ax
        return a.reshape((s.B,) + (1,) * (s.amp.ndim - 2))

    def hadamard(s, w):
        if w not in s.S: s.promote(w)
        ax = s.S.index(w) + 1
        a0 = np.take(s.amp, 0, axis=ax); a1 = np.take(s.amp, 1, axis=ax)
        s.amp = np.stack((a0 + a1, a0 - a1), axis=ax); s.e2 -= 0.5
        s.collapse(w); s.renorm()

    def h_measure(s, w, m):
        """h w; measure w -> m : projection on <(-1)^m|, the wire becomes the classical value m"""
        s.e2 -= 0.5
        if w in s.S:
            ax = s.S.index(w) + 1
            a0 = np.take(s.amp, 0, axis=ax); a1 = np.take(s.amp, 1, axis=ax)
            s.amp = a0 + (-1) ** m * a1; s.S.remove(w)
        else:
            if m: s.ph = (s.ph + 4 * s.bits[w].astype(np.int16)) % 8
        s.bits[w] = m; s.renorm()

    def renorm(s):
        mx = float(np.abs(s.amp).max())
        if mx > 0:
            e = math.frexp(mx)[1]
            if abs(e) > 8: s.amp = s.amp * 2.0 ** -e; s.e2 += e

    def z_measure(s, w, m):
        if w in s.S:
            ax = s.S.index(w) + 1
            s.amp = np.take(s.amp, m, axis=ax); s.S.remove(w)
        else:
            s.amp = s.amp * s.bc(s.bits[w] == m)
        s.bits[w] = m

    def phase_n(s, ws, c):
        """w^{c * prod(ws)} (c=4: CZ/CCZ/Z; for a single wire c is the weight)"""
        cl = [w for w in ws if w not in s.S]; sw = [w for w in ws if w in s.S]
        pc = np.ones(s.B, np.int16)
        for w in cl: pc = pc & s.bits[w]
        if not sw:
            s.ph = (s.ph + c * pc) % 8; return
        idx = [slice(None)] * s.amp.ndim
        for w in sw: idx[s.S.index(w) + 1] = 1
        idx = tuple(idx); sub = s.amp[idx]
        fac = OM[(c * pc) % 8].reshape((s.B,) + (1,) * (sub.ndim - 1))
        s.amp[idx] = sub * fac

    def cx(s, c, t):
        if c not in s.S:
            if t not in s.S:
                s.bits[t] ^= s.bits[c]; return
            ax = s.S.index(t) + 1
            s.amp = np.where(s.bc(s.bits[c] == 1), np.flip(s.amp, axis=ax), s.amp); return
        if t not in s.S: s.promote(t)
        axc = s.S.index(c) + 1; axt = s.S.index(t) + 1
        a0 = np.take(s.amp, 0, axis=axc); a1 = np.take(s.amp, 1, axis=axc)
        axt2 = axt - 1 if axt > axc else axt
        s.amp = np.stack((a0, np.flip(a1, axis=axt2)), axis=axc)
        s.collapse(t)

    def gate(s, name, qs):
        if name == "x":
            w = qs[0]
            if w in s.S: s.amp = np.flip(s.amp, axis=s.S.index(w) + 1)
            else: s.bits[w] ^= 1
        elif name in ("t", "tdg", "s", "sdg", "z"):
            c = {"t": 1, "tdg": 7, "s": 2, "sdg": 6, "z": 4}[name]; w = qs[0]
            if w in s.S:
                idx = [slice(None)] * s.amp.ndim; idx[s.S.index(w) + 1] = 1; s.amp[tuple(idx)] *= OM[c]
            else: s.ph = (s.ph + c * s.bits[w].astype(np.int16)) % 8
        elif name == "cx": s.cx(qs[0], qs[1])
        elif name == "cz": s.phase_n(qs, 4)
        elif name == "ccz": s.phase_n(qs, 4)
        elif name == "ccx":
            s.hadamard(qs[2]); s.phase_n(qs, 4); s.hadamard(qs[2])
        elif name == "swap":
            a, b = qs; s.cx(a, b); s.cx(b, a); s.cx(a, b)
        elif name == "reset":
            w = qs[0]
            if w in s.S: raise ValueError("reset of a wire in superposition")
            s.bits[w] = 0
        elif name == "id": pass
        else: raise ValueError("unsupported gate: " + name)

    def finish(s):
        for w in list(s.S): s.collapse(w)
        return not s.S

    def value(s):
        """full amplitude (complex mant[B], exp): amp * w^ph * 2^e2"""
        return s.amp * OM[s.ph % 8], s.e2


def sim_run(prog, X, outcomes):
    N, nb, ops = prog
    sim = Sim(X, N); c = [0] * nb; i = 0
    while i < len(ops):
        op = ops[i]
        if op[0] == "h":
            q = op[1]; nxt = ops[i + 1] if i + 1 < len(ops) else None
            if nxt and nxt[0] == "m" and nxt[1] == q:
                k = nxt[2]; c[k] = outcomes[k]; sim.h_measure(q, outcomes[k]); i += 2; continue
            sim.hadamard(q)
        elif op[0] == "g": sim.gate(op[1], op[2])
        elif op[0] == "if":
            if all((sum(c[b] for b in g) & 1) for g in op[1]):
                for name, qs in op[2]: sim.gate(name, qs)
        elif op[0] == "m":
            k = op[2]; c[k] = outcomes[k]; sim.z_measure(op[1], outcomes[k])
        i += 1
    return sim


def check_sim(qasm, ref, n=None, inw=None, B=64, m_trials=4, seed=1, exhaustive_max=10, verbose=False, anc="zero", tol=1e-7,
              max_zero=60, chunk=1024, X=None, domain_complete=False):
    """Like check, but with the simulator (memory O(B*2^S), speed ~ gate count): for deep circuits (multipliers, 2048 bits)."""
    t0 = time.time(); rgates = None
    if isinstance(ref, str): n, f, inw, rgates = ref_from_qc(ref, clean=(anc == "zero_clean"))
    elif isinstance(ref, tuple):
        n, f, inw = ref[:3]; rgates = ref[3] if len(ref) > 3 else None
    else: f = ref
    if inw is None: inw = list(range(n))
    prog = parse(qasm); N, nb, ops = prog
    rng = random.Random(seed)
    exhaustive = len(inw) <= exhaustive_max
    if X is not None:                      # user input domain (promise circuits: x<p, o=T[a], ...)
        X = np.asarray(X, np.uint8); exhaustive = bool(domain_complete)
    elif exhaustive:
        Bn = 1 << len(inw); X = np.zeros((Bn, n), np.uint8)
        for j, q in enumerate(inw): X[:, q] = (np.arange(Bn) >> j) & 1
    else:
        X = special_inputs(n, inw, B, rng)
    D = f(X); aref = None
    if rgates is not None:
        rs = Sim(X, n)
        for g, q in rgates: rs.hadamard(q[0]) if g == "h" else rs.gate(g, q)
        if not rs.finish() or not (rs.bits[:n].T == D).all(): return dict(ok=False, error="reference: no classical output")
        am, ae = rs.value(); aref = am * 2.0 ** ae
    worst = 0.0; zero_br = 0; res_br = []; done = 0; maxS = 0
    while done < m_trials:
        m = [rng.randint(0, 1) for _ in range(nb)]
        ams = []; E0 = None; allzero = True; bad = 0.0; ancref = {}
        for s0 in range(0, len(X), chunk):
            sl = slice(s0, min(len(X), s0 + chunk))
            try:
                sim = sim_run(prog, X[sl], m)
            except TooWide as e:
                return dict(ok=False, error="TooWide " + str(e))
            except ValueError as e:
                return dict(ok=False, error=str(e))
            maxS = max(maxS, sim.maxS)
            if not sim.finish():
                return dict(ok=False, error="wires %s remain in superposition (the branch is not classical)" % sim.S)
            am, ae = sim.value()
            if aref is not None: am = am / aref[sl]
            nz = np.abs(am) > 1e-150
            if nz.any(): allzero = False
            # bits: data = D on EVERY reference wire (clean wires included), live ancillas = 0
            ok_all = np.ones(sl.stop - sl.start, bool); meas = _measured_wires(prog)
            for q in range(N):
                if q < n:
                    ok_all &= (sim.bits[q] == D[sl, q])   # v0.4 fix: clean output wires were skipped before
                elif q in meas: continue          # measured extra qubit: holds its outcome, not an output
                elif anc == "const":              # strict: the ancilla is one basis state, the same for all x
                    if q not in ancref: ancref[q] = sim.bits[q][0]
                    ok_all &= (sim.bits[q] == ancref[q])
                else: ok_all &= (sim.bits[q] == 0)
            bad = max(bad, float((~ok_all & nz).mean()) if nz.any() else 0.0)
            if E0 is None: E0 = ae
            ams.append(am * 2.0 ** (ae - E0))
        if allzero:
            zero_br += 1
            if zero_br > max_zero: return dict(ok=False, error="too many zero branches")
            continue
        am = np.concatenate(ams); a0 = am[0]
        devph = 1.0 if abs(a0) < 1e-150 else float(np.abs(am / a0 - 1).max())
        logp = None if abs(a0) < 1e-150 else float(2 * (math.log2(abs(a0)) + E0))
        w = max(bad, devph); worst = max(worst, w); done += 1
        res_br.append(dict(outcomes_head=m[:12], dev_bits=bad, dev_phase=devph, log2p=logp))
        if verbose: print("  branch %d: bits_bad=%.2e dev_phase=%.2e log2p=%s maxS=%d t=%.1fs" % (done, bad, devph, logp, maxS, time.time() - t0), flush=True)
    return dict(ok=bool(worst < tol), worst=worst, proved_all_x=exhaustive, n=n, N=N, nb=nb, B=len(X), branches=m_trials,
                zero_branches=zero_br, maxS=maxS, method="sim", time=time.time() - t0, detail=res_br)


_MW = {}


def _measured_wires(prog):
    """wires measured in the circuit (an extra qubit keeps its outcome after measurement and is not an output)"""
    k = id(prog)
    if k not in _MW:
        ws = set()
        for op in prog[2]:
            if op[0] == "m": ws.add(op[1])
        _MW[k] = ws
    return _MW[k]


# ============================================================ references
def load_qc(path):
    names, gates = None, []
    ONE = {"H": "h", "X": "x", "T": "t", "T*": "tdg", "S": "s", "S*": "sdg", "P": "s", "P*": "sdg"}
    ins = None
    for line in open(path):
        p = line.split()
        if not p or p[0].startswith("#"): continue
        if p[0] == ".v":
            names = p[1:]; idx = {q: i for i, q in enumerate(names)}; continue
        if p[0] == ".i":
            ins = p[1:]; continue
        if p[0].startswith(".") or p[0] in ("BEGIN", "END"): continue
        op, qs = p[0], tuple(idx[q] for q in p[1:])
        if op in ("Z", "Zd"):
            qs = tuple(dict.fromkeys(qs))
            gates.append(({1: "z", 2: "cz", 3: "ccz"}[len(qs)], qs))
        elif op == "tof":
            gates.append(({1: "x", 2: "cx", 3: "ccx"}[len(qs)], qs))
        elif op in ONE and len(qs) == 1:
            gates.append((ONE[op], qs))
        else:
            raise ValueError("unsupported: " + line.strip())
    return names, gates, ([idx[q] for q in ins] if ins is not None else None)


def classical_gates(gates, X, lenient=False):
    """Vectorised classical simulation of Toffoli-level gates (with H frames h Z.. h = ccx). X: (B,n) uint8 -> (B,n).
    ValueError if the circuit is not classical (phase gates / open frames)."""
    b = X.astype(np.uint8).copy(); frame = set()
    for g, q in gates:
        if g == "h":
            frame ^= {q[0]}; continue
        fr = [u for u in q if u in frame]
        if g == "x":
            if fr:
                if lenient: continue
                raise ValueError("x inside an H frame = Z phase (not classical)")
            b[:, q[0]] ^= 1
        elif g == "cx":
            if fr: raise ValueError("cx on an H frame")
            b[:, q[1]] ^= b[:, q[0]]
        elif g == "ccx":
            if fr: raise ValueError("ccx on an H frame")
            b[:, q[2]] ^= b[:, q[0]] & b[:, q[1]]
        elif g in ("z", "cz", "ccz"):
            if not fr:
                if lenient: continue
                raise ValueError("Z phase outside an H frame (not classical)")
            if len(fr) != 1: raise ValueError("two H frames in one Z")
            t = fr[0]; ctl = [u for u in q if u != t]
            m = np.ones(len(b), np.uint8)
            for u in ctl: m &= b[:, u]
            b[:, t] ^= m
        else:
            raise ValueError("not a classical gate " + g)
    if frame: raise ValueError("open H frame")
    return b


def ref_from_qc(path, clean=False):
    """-> (n, f, input_wires). f(X(B,n))->(B,n). clean: wires not on the .i line start at 0."""
    names, gates, ins = load_qc(path)
    n = len(names)
    inw = ins if (clean and ins is not None) else list(range(n))
    return n, (lambda X: classical_gates(gates, X, lenient=True)), inw, gates


# ============================================================ sampled check over x
def special_inputs(n, inw, B, rng):
    """Input set: zeros, ones, alternating, then a mix: uniform, sparse, dense, runs (long carry chains)."""
    rows = []
    L = len(inw)

    def mk(bits):
        r = np.zeros(n, np.uint8)
        for j, q in enumerate(inw): r[q] = bits[j]
        return r
    rows += [mk([0] * L), mk([1] * L), mk([j & 1 for j in range(L)]), mk([1 - (j & 1) for j in range(L)])]
    k = 0
    while len(rows) < B:
        kind = k % 4; k += 1
        if kind == 0: bits = [rng.randint(0, 1) for _ in range(L)]
        elif kind == 1: bits = [1 if rng.random() < 1 / 16 else 0 for _ in range(L)]
        elif kind == 2: bits = [0 if rng.random() < 1 / 16 else 1 for _ in range(L)]
        else:
            bits = []; v = rng.randint(0, 1)
            while len(bits) < L:
                bits += [v] * rng.randint(1, max(2, L // 3)); v ^= 1
            bits = bits[:L]
        rows.append(mk(bits))
    return np.array(rows[:B], np.uint8)


def check_pathsum(qasm, ref, n=None, inw=None, B=64, m_trials=4, seed=1, norm=True, exhaustive_max=10, verbose=False,
          anc="zero", tol=1e-7, max_zero=60, maxscope=22, chunk=256, X=None, domain_complete=False):
    """qasm: path; ref: path to .qc OR (n, f, inw[, gates]) OR callable f(X)->D (then n is required).
    Returns dict(ok, proved_all_x (ALL x enumerated), worst, ...). anc: 'zero' (ancillas return to 0), 'const' (sim only: the ancilla is
    one basis state, the same for all x, like strict in verify_qasm), 'zero_clean' (like zero, and the .qc wires not on .i
    are clean ancillas starting at 0; every reference wire must hold its expected value)."""
    t0 = time.time()
    rgates = None
    if isinstance(ref, str):
        n, f, inw, rgates = ref_from_qc(ref, clean=(anc == "zero_clean"))
    elif isinstance(ref, tuple):
        n, f, inw = ref[:3]; rgates = ref[3] if len(ref) > 3 else None
    else:
        f = ref
    if inw is None: inw = list(range(n))
    prog = parse(qasm); N, nb, ops = prog
    rng = random.Random(seed)
    exhaustive = len(inw) <= exhaustive_max
    if X is not None:
        X = np.asarray(X, np.uint8); exhaustive = bool(domain_complete)
    elif exhaustive:
        Bn = 1 << len(inw)
        X = np.zeros((Bn, n), np.uint8)
        for j, q in enumerate(inw): X[:, q] = (np.arange(Bn) >> j) & 1
    else:
        X = special_inputs(n, inw, B, rng)
    D = f(X)
    aref = None
    if rgates is not None:   # reference amplitude (with phases) by the same engine: also checks d*
        rnet = Net(n, X)
        for g, q in rgates:
            if g == "h": rnet.hadamard(q[0])
            else: rnet.gate(g, q)
        ram, rae = amplitude(rnet, {q: D[:, q] for q in range(n)})
        if rae > 100 or rae < -100 or not np.allclose(np.abs(ram) * 2.0 ** rae, 1.0, atol=1e-9):
            return dict(ok=False, error="reference: |A_ref(d*)| != 1 (no classical output)")
        aref = ram * 2.0 ** rae
    init = np.zeros((len(X), N), np.uint8); init[:, :n] = X
    worst = 0.0; zero_br = 0; res_br = []; done = 0; Bc = min(chunk, len(X))
    while done < m_trials:
        m = [rng.randint(0, 1) for _ in range(nb)]
        while True:       # adaptive batch-chunk reduction when memory is short
            try:
                ams = []; devn = 0.0; allzero = True; E0 = None; logp = None
                for s0 in range(0, len(X), Bc):
                    sl = slice(s0, min(len(X), s0 + Bc))
                    try:
                        net, c = build(prog, init[sl], m)
                    except ValueError as e:
                        return dict(ok=False, error=str(e))
                    targets = {q: D[sl, q] for q in range(n)}
                    for q in range(n, N):
                        if net.alive[q]: targets[q] = np.zeros(sl.stop - sl.start, np.uint8)
                    for q in range(n):
                        if not net.alive[q]:   # measured reference wire: must hold its expected constant (v0.4)
                            mq_, cq_ = net.w[q]
                            if mq_ != 0 or (np.asarray(cq_) != D[sl, q]).any(): return dict(ok=False, error="measured reference wire %d does not hold its output" % q)
                            targets.pop(q, None)
                    am, ae = amplitude(net, targets)
                    if norm: nm_, ne = norm2(net, maxscope)
                    if aref is not None: am = am / aref[sl]
                    a2 = np.abs(am) ** 2
                    if norm:
                        if not ((np.abs(nm_) < 1e-300).all() and (a2 < 1e-300).all()): allzero = False
                        dd = 2 * ae - ne
                        lhs = a2 * 2.0 ** dd if abs(dd) < 900 else a2 * 0
                        devn = max(devn, float((np.abs(lhs - np.abs(nm_)) / np.maximum(np.abs(nm_), 1e-300)).max()))
                    else:
                        allzero = allzero and (a2 < 1e-300).all()
                    if E0 is None: E0 = ae
                    ams.append(am * 2.0 ** (ae - E0) if abs(ae - E0) < 900 else am * 0)
                break
            except TooWide as e:
                if Bc > 1 and "memory" in str(e): Bc = max(1, Bc // 4); continue
                return dict(ok=False, error="TooWide " + str(e))
        if allzero and norm:
            zero_br += 1
            if zero_br > max_zero: return dict(ok=False, error="too many zero branches")
            continue
        am = np.concatenate(ams)
        a0 = am[0]
        if abs(a0) < 1e-150: devph = 1.0
        else:
            devph = float(np.abs(am / a0 - 1).max())
            logp = float(2 * (math.log2(abs(a0)) + E0))
        w = max(devn, devph); worst = max(worst, w); done += 1
        res_br.append(dict(outcomes_head=m[:12], dev_norm=devn, dev_phase=devph, log2p=logp))
        if verbose: print("  branch %d: dev_norm=%.2e dev_phase=%.2e log2p=%s t=%.1fs" % (done, devn, devph, logp, time.time() - t0), flush=True)
    return dict(ok=bool(worst < tol), worst=worst, proved_all_x=exhaustive, n=n, N=N, nb=nb, B=len(X),
                branches=m_trials, zero_branches=zero_br, norm_checked=norm, time=time.time() - t0, detail=res_br)


# ============================================================ PROOF over ALL x at once (miter)
def qc_netlist(gates):
    """.qc gates -> classical network (x/cx/ccx), H frames h Z.. h = ccx. Strict: ValueError if not classical."""
    out = []; frame = set()
    for g, q in gates:
        if g == "h": frame ^= {q[0]}; continue
        fr = [u for u in q if u in frame]
        if g in ("x", "cx", "ccx"):
            if fr: raise ValueError(g + " on an H frame")
            out.append((g, tuple(q)))
        elif g in ("z", "cz", "ccz"):
            if len(fr) != 1: raise ValueError("Z outside / in several H frames")
            t = fr[0]; ctl = tuple(u for u in q if u != t)
            out.append(({0: "x", 1: "cx", 2: "ccx"}[len(ctl)], ctl + (t,)))
        else:
            raise ValueError("not a classical gate " + g)
    if frame: raise ValueError("open H frame")
    return out


def spec_from_qc(path, clean=False):
    names, gates, ins = load_qc(path)
    n = len(names)
    return n, qc_netlist(gates), (ins if (clean and ins is not None) else list(range(n)))


def _ref_factors(fg, n, netlist, xv):
    """Boolean network of the reference: one variable per new wire value; xv: {wire: x variable}. -> [(var|None, const)] per wire"""
    W = [(xv.get(q), 0) for q in range(n)]
    def val_table(srcs, fn):
        # srcs: list of (var|None, c); table over the unique variables
        vs = sorted({v for v, c in srcs if v is not None})
        return vs
    for item in netlist:
        g, qs = item[0], item[1]
        if g == "x":
            v, c = W[qs[0]]; W[qs[0]] = (v, c ^ 1); continue
        srcs = [W[q] for q in qs]
        z = fg.newvar()
        vs = sorted({v for v, c in srcs if v is not None})
        scope = sorted(vs + [z])
        k = len(scope); pos = {u: i for i, u in enumerate(scope)}
        tbl = np.zeros((2,) * k, complex)
        for assign in range(1 << (k - 1)):
            val = dict(zip(vs, [(assign >> i) & 1 for i in range(len(vs))]))
            bits = [(val[v] if v is not None else 0) ^ c for v, c in srcs]
            if g == "cx": r = bits[1] ^ bits[0]
            elif g == "lut":     # t ^= tab[ins] ; item[2] is the table (index = bits of ins, ins[0] most significant)
                idx = 0
                for b_ in bits[:-1]: idx = idx * 2 + b_
                r = bits[-1] ^ ((item[2][idx]) & 1)
            else: r = bits[2] ^ (bits[0] & bits[1])
            idx = [0] * k
            for u in vs: idx[pos[u]] = val[u]
            idx[pos[z]] = r
            tbl[tuple(idx)] = 1.0
        fg.add(scope, tbl.reshape((1,) + tbl.shape))
        W[qs[-1]] = (z, 0)
    return W


def check_proof(qasm, spec, m_trials=6, seed=1, tol=1e-6, maxscope=24, verbose=False, clean=False):
    """spec = (n, netlist, inw) or a path to .qc. For each chosen branch m PROVES for ALL x (2^|inw| inputs) at once:
        B_m |x> = t_m |f(x), 0_anc>   (t_m independent of x),   via  ||B_m||_F^2 == |Tr(V^+ B_m)|^2 / 2^{nx}.
    Branches m are sampled (m_trials); no claim is made for all m."""
    t0 = time.time()
    if isinstance(spec, str): spec = spec_from_qc(spec, clean)
    n, netlist, inw = spec
    prog = parse(qasm); N, nb, ops = prog
    rng = random.Random(seed)
    xvars = {q: j for j, q in enumerate(inw)}; nx = len(inw)
    worst = 0.0; details = []; zero_br = 0; done = 0; tries = 0
    mlist = [[0] * nb, [1] * nb]
    while done < m_trials:
        tries += 1
        m = mlist.pop(0) if mlist else [rng.randint(0, 1) for _ in range(nb)]
        net = build_sym(prog, N, xvars, m)
        def buildT(nlast):
            fg = FG(1, net.k, maxscope=maxscope); fg.tag = 'T'; fg.nlast = nlast
            net_factors(net, fg)
            W = _ref_factors(fg, n, netlist, {q: j for q, j in xvars.items() if q < n})
            for q in range(N):
                mq, cq = net.w[q]
                if q < n:
                    # v0.4: a measured reference wire is constrained like a live one (its final value must equal the
                    # reference output); it used to be skipped in the clean model
                    rv, rc = W[q]
                    parts = [(mq, 0)]
                    bit = int(cq[0]) ^ rc
                    vs = list(fg.pv(mq, 0)) + ([rv] if rv is not None else [])
                else:
                    if not net.alive[q]: continue
                    vs = list(fg.pv(mq, 0)); bit = int(cq[0])
                cnt = defaultdict(int)
                for v in vs: cnt[v] ^= 1
                vs = sorted(v for v, cc in cnt.items() if cc)
                if not vs:
                    if bit: fg.zero = True
                    continue
                par = _pt(vs, vs)
                fg.add(vs, (par == bit).astype(complex).reshape((1,) + par.shape))
            return fg
        try:
            fg = buildT(0)
        except ValueError as e:
            return dict(ok=False, error=str(e))
        try:
            try:
                fg.eliminate()
            except TooWide:
                fg = buildT(nx)
                fg.eliminate()
            Tm, Te = fg.result()
            try:
                Nm, Ne = norm2(net, maxscope, nshared=nx, nlast=0)
            except TooWide:
                Nm, Ne = norm2(net, maxscope, nshared=nx, nlast=nx)
        except TooWide as e:
            return dict(ok=False, error="TooWide " + str(e))
        Nv = abs(Nm[0])
        if Nv < 1e-300:
            zero_br += 1
            if zero_br > 40: return dict(ok=False, error="too many zero branches")
            continue
        Tv = abs(Tm[0])
        if Tv < 1e-300:
            dev = 1.0; lr = None
        else:
            lr = 2 * (math.log2(Tv) + Te) - nx - (math.log2(Nv) + Ne)
            dev = abs(2.0 ** lr - 1) if abs(lr) < 50 else 1.0
        worst = max(worst, dev); done += 1
        details.append(dict(dev=dev, log2_norm=math.log2(Nv) + Ne))
        if verbose: print("  branch %d: dev=%.2e log2||B||^2=%.3f t=%.1fs" % (done, dev, math.log2(Nv) + Ne, time.time() - t0), flush=True)
    return dict(ok=bool(worst < tol), worst=worst, proved_all_x=True, proof_of_branches=m_trials, n=n, nx=nx, N=N, nb=nb,
                zero_branches=zero_br, time=time.time() - t0, detail=details)


def build_sym(prog, N, xvars, outcomes):
    net = Net(N, None, xvars)
    N_, nb, ops = prog
    c = [0] * nb; i = 0
    while i < len(ops):
        op = ops[i]
        if op[0] == "h":
            q = op[1]; nxt = ops[i + 1] if i + 1 < len(ops) else None
            if nxt and nxt[0] == "m" and nxt[1] == q:
                k = nxt[2]; c[k] = outcomes[k]; net.measure_after_h(q, outcomes[k]); i += 2; continue
            net.hadamard(q)
        elif op[0] == "g": net.gate(op[1], op[2])
        elif op[0] == "if":
            if all((sum(c[b] for b in g) & 1) for g in op[1]):
                for name, qs in op[2]: net.gate(name, qs)
        elif op[0] == "m":
            k = op[2]; c[k] = outcomes[k]; net.measure_z(op[1], outcomes[k])
        i += 1
    return net


def check(qasm, ref, method="auto", **kw):
    """Main sampled check. method: 'sim' (simulator), 'pathsum' (factor graph + full norm), 'auto' (sim, falling back
    to pathsum when the width is too large). Returns dict(ok, worst, proved_all_x, method, ...)."""
    import inspect
    def pick(fn): return {k: v for k, v in kw.items() if k in inspect.signature(fn).parameters}
    if method == "pathsum": return check_pathsum(qasm, ref, **pick(check_pathsum))
    r = check_sim(qasm, ref, **pick(check_sim))
    if method == "auto" and not r.get("ok") and str(r.get("error", "")).startswith("TooWide"):
        r = check_pathsum(qasm, ref, **pick(check_pathsum)); r["method"] = "pathsum"
    elif "method" not in r: r["method"] = "sim"
    return r


def main(argv=None):
    ap = argparse.ArgumentParser(prog="qoptct-arith", description="qoptct.verify.arith: checks of arithmetic circuits (classical .qc reference) at large sizes")
    ap.add_argument("qasm"); ap.add_argument("qc", help="reference .qc (classical: tof/cx/x and H frames h Z.. h)")
    ap.add_argument("--clean", action="store_true", help=".qc wires not on the .i line are clean ancillas")
    ap.add_argument("--anc", default=None, choices=["zero", "const"], help="final state of live ancillas: 0 (default) or a constant (strict, sim only)")
    ap.add_argument("--method", default="auto", choices=["auto", "sim", "pathsum"])
    ap.add_argument("-B", type=int, default=64, help="number of random inputs (all inputs if there are at most 10 input bits)")
    ap.add_argument("-m", type=int, default=6, help="number of random measurement branches")
    ap.add_argument("--proof", action="store_true", help="also prove over ALL inputs (miter, factor graph) for the sampled branches")
    ap.add_argument("--no-norm", action="store_true"); ap.add_argument("-v", action="store_true")
    a = ap.parse_args(argv)
    anc = "zero_clean" if a.clean else (a.anc or "zero")
    r = check(a.qasm, a.qc if not isinstance(a.qc, str) else a.qc, method=a.method, B=a.B, m_trials=a.m, norm=not a.no_norm, anc=anc, verbose=a.v)
    r.pop("detail", None); print("check:", json.dumps(r, ensure_ascii=False)); ok = bool(r.get("ok"))
    if a.proof:
        p = check_proof(a.qasm, spec_from_qc(a.qc, a.clean), m_trials=a.m, clean=a.clean, verbose=a.v)
        p.pop("detail", None); print("proof:", json.dumps(p, ensure_ascii=False)); ok = ok and bool(p.get("ok"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
