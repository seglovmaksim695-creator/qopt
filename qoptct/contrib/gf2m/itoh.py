# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""itoh: inversion and division in GF(2^m) by the Itoh-Tsujii chain, with FREE erasure of product registers
[qoptct v0.5; lab k6 gf2ec.py].

x^-1 = x^(2^m - 2) = (x^(2^(m-1) - 1))^2.  An addition chain 1 = e_0, ..., e_L = m-1 (e_s = e_i + e_j) gives the
chain registers  b_s = x^(2^(e_s) - 1) = b_i^(2^(e_j)) * b_j  (Frobenius powers are linear: free inside the bilinear
algorithm), one GF(2^m) multiplication each.

Erasure of a product register.  A register t = F(A) * G(B) (F, G linear, A and B live registers, A may equal B) is
erased by measuring every bit in the X basis and applying the phase (-1)^(s.t) -- s.t is a quadratic form in the
bits of A and B, so the correction is Z / CZ gates only (0 T).  Hence chain registers are removed without
re-computation as long as their parents are alive ("pebble game with free removal", `pebble_plan`), and with
p = L chain registers
    division   T ^= y / x :   (L + 1) multiplications       (the classical FLT/IT circuits uncompute: 2L + 1)
    inversion  T ^= x^-1  :    L multiplications
of which each costs R logical ANDs (R = rank of the bilinear algorithm, 908 for m = 163), e.g. m = 163, L = 9:
division 10 * 908 = 9080 ANDs.

Prior art (honest statement).  Measurement-based uncomputation (MBU) of a GF(2^m) product register by X
measurement + phase correction is published: Khattar, Shutty, Gidney et al., arXiv:2510.10967, App. C.3 (Qualtran
GF2MulMBUC).  New here (lab k6) is its application to the Itoh-Tsujii chain (pebbling with free removal) and to
binary-curve point addition; against Taguchi-Takagi (ePrint 2023/553, Table 4: 19 * M) the cost model is the same.

Gate level.  Every product is lowered term by term: parities of the pulled-back masks on fresh wires (CX), one
logical AND into a fresh ancilla, fan-out into the target bits, measurement erasure, parities undone (simple, CX
heavy: ~3.7M CX for the m = 163 division; the in-place CRT compiler of mult.py is not used here).
Register erasure, erase_mode='terms' (default): all bits of t are measured in the X basis (outcomes s), then
(-1)^(s.t) = prod over terms k and bits i of gamma_k of (-1)^(s_i (alpha_k.A)(beta_k.B)) is applied as one
classically controlled CZ between the two parity wires per (k, i): O(rank * m) CZ plus the parity CX of one
product.  erase_mode='monomials' (the k6 lab form): one H / measure / ``if (c) {Z.., CZ..; X}`` per bit with the
monomials of the bit's quadratic form, O(m^3) CZ.  Both 0 T.  For large m the register-level program gives exact
AND counts and a gate-level stream check (gf2m_itoh_gate_check) without storing the circuit.

Why erase='terms' is correct for every m (qoptct 0.5.0).  Let t = F_a(A) * F_b(B) be the register to erase, A, B
live registers, F_a, F_b linear (Frobenius powers).  A bilinear algorithm (alpha_k, beta_k, gamma_k) of the field
product, applied to A' = F_a(A), B' = F_b(B), gives  t = sum_k gamma_k (alpha_k.A')(beta_k.B')  and
alpha.A' = pull(alpha).A (pull_mask), so  t = sum_k gamma_k (a_k.A)(b_k.B)  with pulled-back masks a_k, b_k.
  1. Measuring every bit of t in the X basis (H, measure, outcome s_i, conditional X back to |0>) maps a basis
     state |A, B, t> to 2^(-m/2) (-1)^(s.t) |A, B, 0> in branch s = (s_1..s_m).
  2. s.t = sum_i s_i t_i = sum_k sum_{i in gamma_k} s_i (a_k.A)(b_k.B)  (mod 2), and (-1)^x depends on x mod 2
     only, so (-1)^(s.t) = prod over (k, i in gamma_k, s_i = 1) of (-1)^((a_k.A)(b_k.B)).
  3. Each factor is one CZ between two wires that hold a_k.A and b_k.B at that moment (parity wires: CX of the bits
     into a |0> wire, undone afterwards; CRT lowering: in-place invertible CX transforms of A and B, restored
     exactly), classically controlled by s_i; a term with A = B and a_k = b_k is (a_k.A)^2 = a_k.A -> one Z.
     (The CRT lowering takes its terms from the CRT bilinear algorithm of mult.py, the parity lowering from `dec`:
     any correct bilinear algorithm of the same product gives the same s.t, so the argument is the same.)
  4. Hence after the corrections every branch is 2^(-m/2) |A, B, 0>: the same isometry for every basis input, i.e.
     the erasure equals the coherent uncomputation of t on all inputs (and superpositions), for every m, provided
     the classical precondition (t = the asserted bilinear form, the CX bookkeeping returns A, B and the work wires)
     holds -- which the register program checks on samples, `GF2Circuit.simulate` on given inputs and the
     exhaustive tests on all inputs for small m.
Machine checks of the emitted gates: qcirc.ErasurePhaseChecker recomputes, symbolically and for any m, the phase
polynomial sum_i s_i Q_i(A, B) actually applied by the conditional Z / CZ of every erasure and compares it with the
asserted form (exact for all inputs and outcomes; m = 163 division: all 668 520 conditional phases); kprove proves
whole circuits against an independent truth-table reference (inverse m <= 14, division m <= 11 with a 1.4 GB
budget; tests/test_crypto_gf2m.py), and a removed conditional phase is REFUTED.
"""
import heapq

from .field import pmod, pmul
from .qcirc import GF2Circuit


# ------------------------------------------------------------------------------------------- field helpers

class Field:
    def __init__(self, m, f):
        self.m, self.f = m, f
        self._fc = {}

    def mul(self, a, b):
        return pmod(pmul(a, b), self.f)

    def sq(self, a):
        r = 0
        i = 0
        while a:
            if a & 1:
                r |= 1 << (2 * i)
            a >>= 1
            i += 1
        return pmod(r, self.f)

    def frob(self, a, k):
        for _ in range(k % self.m):
            a = self.sq(a)
        return a

    def frob_cols(self, k):
        k %= self.m
        if k not in self._fc:
            self._fc[k] = [self.frob(1 << j, k) for j in range(self.m)]
        return self._fc[k]

    def inv(self, a):
        r, b, e = 1, a, (1 << self.m) - 2
        while e:
            if e & 1:
                r = self.mul(r, b)
            b = self.sq(b)
            e >>= 1
        return r


def _parity(x):
    return bin(x).count("1") & 1


def pull_mask(alpha, cols):
    """beta with alpha.(L x) = beta.x for the linear map L with columns cols."""
    b = 0
    for j, c in enumerate(cols):
        if _parity(alpha & c):
            b |= 1 << j
    return b


# ------------------------------------------------------------------------------------------- chains, pebbling

def it_chain(n):
    """Itoh-Tsujii chain for n (= m - 1): list of (e, i, j), e = e_i + e_j; node 0 is e = 1."""
    ch = [(1, None, None)]
    cur = 0
    for b in bin(n)[3:]:
        e = ch[cur][0]
        ch.append((2 * e, cur, cur))
        cur = len(ch) - 1
        if b == "1":
            e = ch[cur][0]
            ch.append((e + 1, cur, 0))
            cur = len(ch) - 1
    return ch


def short_chain(n, maxlen=None):
    """a shortest star addition chain for n (iterative deepening; practical for n < ~1000)."""
    import math
    if n == 1:
        return [(1, None, None)]
    lo = math.ceil(math.log2(n))
    best = None
    for L in range(lo, (maxlen or 2 * lo) + 1):
        seq = [1]

        def dfs():
            k = len(seq) - 1
            if seq[-1] == n:
                return True
            if k == L or seq[-1] << (L - k) < n:
                return False
            last = seq[-1]
            tried = set()
            for j in range(k, -1, -1):
                e = last + seq[j]
                if e > n or e in tried or e <= last:
                    continue
                tried.add(e)
                seq.append(e)
                if dfs():
                    return True
                seq.pop()
            return False
        if dfs():
            best = list(seq)
            break
    ch = [(1, None, None)]
    idx = {1: 0}
    for k in range(1, len(best)):
        e = best[k]
        d = e - best[k - 1]
        ch.append((e, k - 1, idx[d]))
        idx[e] = k
    return ch


def pebble_plan(ch, p, use_p=None):
    """fewest placements (= multiplications) to reach a state holding the sink and to return to the empty state with
    at most p chain registers (node 0 = the input, always present).  Removal is free but needs both parents.
    -> (cost_up, cost_down, moves_up, moves_down) with moves ('P', s) / ('E', s); None if impossible."""
    L = len(ch) - 1
    sink = L
    par = [0] * (L + 1)
    for s in range(1, L + 1):
        _, i, j = ch[s]
        mask = 0
        for q in (i, j):
            if q:
                mask |= 1 << q
        par[s] = mask

    def pc(S):
        return bin(S).count("1")

    dist = {0: 0}
    prev = {}
    pq = [(0, 0)]
    while pq:
        d, S = heapq.heappop(pq)
        if d > dist[S]:
            continue
        for s in range(1, L + 1):
            b = 1 << s
            if par[s] & ~S:
                continue
            if S & b:
                T, c, mv = S ^ b, 0, ('E', s)
            elif pc(S) < p:
                T, c, mv = S | b, 1, ('P', s)
            else:
                continue
            if d + c < dist.get(T, 1 << 30):
                dist[T] = d + c
                prev[T] = (S, mv)
                heapq.heappush(pq, (d + c, T))
    distb = {0: 0}
    prevb = {}
    pq = [(0, 0)]
    while pq:
        d, S = heapq.heappop(pq)
        if d > distb[S]:
            continue
        for s in range(1, L + 1):
            b = 1 << s
            if S & b:
                U = S ^ b
                if not (par[s] & ~U) and pc(U) < p and d + 1 < distb.get(U, 1 << 30):
                    distb[U] = d + 1
                    prevb[U] = (S, ('P', s))
                    heapq.heappush(pq, (d + 1, U))
            else:
                U = S | b
                if not (par[s] & ~U) and pc(U) <= p and d < distb.get(U, 1 << 30):
                    distb[U] = d
                    prevb[U] = (S, ('E', s))
                    heapq.heappush(pq, (d, U))
    best = None
    for S, d in dist.items():
        if S >> sink & 1 and S in distb and (use_p is None or pc(S) <= use_p):
            c = (d, distb[S])
            if best is None or sum(c) < sum(best[0]) or (sum(c) == sum(best[0]) and c[0] < best[0][0]):
                best = (c, S)
    if best is None:
        return None
    (cu, cd), S = best
    up = []
    T = S
    while T != 0:
        U, mv = prev[T]
        up.append(mv)
        T = U
    up.reverse()
    down = []
    T = S
    while T != 0:
        U, mv = prevb[T]
        down.append(mv)
        T = U
    return cu, cd, up, down


def make_chain(m, chain="it"):
    if chain in ("it", "itoh-tsujii", None):
        return it_chain(m - 1)
    if chain == "short":
        return short_chain(m - 1)
    ch = [tuple(c) for c in chain]
    if ch[-1][0] != m - 1:
        raise ValueError("chain does not end at m-1")
    for s in range(1, len(ch)):
        e, i, j = ch[s]
        if not (0 <= i < s and 0 <= j < s) or ch[i][0] + ch[j][0] != e:
            raise ValueError("not an addition chain at node %d" % s)
    return ch


# ------------------------------------------------------------------------------------------- register program

class CheckError(AssertionError):
    pass


class StreamCircuit:
    """gate-level executor with the GF2Circuit building interface: every gate is applied at once to a bitsliced
    state (one int per wire, one bit per sample lane) with its precondition enforced, and only counted, not stored.
    Register erasures are checked lazily (field product of the lanes' register values) -- for large m, where the
    explicit Z/CZ bodies (O(m^3) gates) are not built."""

    def __init__(self, lanes):
        self.full = (1 << lanes) - 1
        self.lanes = lanes
        self.st = {}
        self.n = dict(c=0, x=0, a=0, e=0, r=0)
        self.nq = 0

    def g(self, w):
        return self.st.get(w, 0)

    def cx(self, c, t):
        self.st[t] = self.g(t) ^ self.g(c)
        self.n['c'] += 1

    def x(self, t):
        self.st[t] = self.g(t) ^ self.full
        self.n['x'] += 1

    def AND(self, p, q, a):
        if self.g(a):
            raise CheckError("AND target not zero")
        self.st[a] = self.g(p) & self.g(q)
        self.n['a'] += 1

    def erase(self, p, q, a):
        if self.g(a) != self.g(p) & self.g(q):
            raise CheckError("erasure precondition violated")
        self.st[a] = 0
        self.n['e'] += 1

    def read(self, wires):
        out = [0] * self.lanes
        for i, w in enumerate(wires):
            x = self.g(w)
            while x:
                ln = (x & -x).bit_length() - 1
                x &= x - 1
                out[ln] |= 1 << i
        return out

    def write(self, wires, vals):
        for i, w in enumerate(wires):
            self.st[w] = sum(1 << ln for ln, v in enumerate(vals) if v >> i & 1)

    def lazy_erase(self, wt, wa, ka, wb, kb, F):
        A, B, T = self.read(wa), self.read(wb), self.read(wt)
        for a, b, t in zip(A, B, T):
            if F.mul(F.frob(a, ka), F.frob(b, kb)) != t:
                raise CheckError("register-erasure precondition violated (gate level)")
        for w in wt:
            self.st[w] = 0
        self.n['r'] += len(wt)

    def counts(self):
        n = self.n
        return dict(AND=n['a'], T=4 * n['a'], CX_lin=n['c'], measurements=n['e'] + n['r'],
                    CX=n['c'] + 3 * n['a'] + n['e'])     # + the conditional CZ of the register erasures


class Program:
    """Register-level program over GF(2^m).  Values are kept on S sample inputs (field elements) and every erasure
    condition is checked on them; with `gates=True` every operation is also lowered to a GF2Circuit.
    dec: bilinear algorithm (alpha, beta, gamma) of c = a*b mod f."""

    def __init__(self, m, f, dec, samples, gates=False, count_cz=False, erase_mode="terms", crt=None):
        self.F = Field(m, f)
        self.m, self.dec = m, dec
        self.count_cz = count_cz
        self.crt = crt               # mult.CRTData -> CRT lowering of products / erasures (else parity lowering)
        self._fresh = set()
        if erase_mode not in ("terms", "monomials"):
            raise ValueError("erase_mode must be 'terms' or 'monomials'")
        self.erase_mode = erase_mode
        self.S = len(samples[next(iter(samples))]) if samples else 0
        self.val = {}
        self.stream = gates == "stream"
        self.gates = bool(gates)
        self.C = StreamCircuit(self.S) if self.stream else (GF2Circuit(0) if gates else None)
        self.wires = {}
        self._free = []
        self._nw = 0
        self._live = 0
        self.peak_wires = 0
        self.live_regs = set()
        self.peak_regs = 0
        self.n_mul = 0
        self.n_erase = 0
        self.and_count = 0           # exact number of logical ANDs (nonzero bilinear terms after pull-back)
        self.cz_count = 0            # exact number of CZ/Z in register erasures
        self._terms_cache = {}
        for name, vals in samples.items():
            self.reg(name, vals)

    # wires
    def _new(self):
        self._live += 1
        self.peak_wires = max(self.peak_wires, self._live)
        if self._free:
            return self._free.pop()
        w = self._nw
        self._nw += 1
        return w

    def _release(self, w):
        if self.stream and self.C.g(w):
            raise CheckError("released wire q[%d] not zero (gate level)" % w)
        self._live -= 1
        self._free.append(w)

    # registers
    def reg(self, name, vals):
        self.val[name] = list(vals)
        self.live_regs.add(name)
        self.peak_regs = max(self.peak_regs, len(self.live_regs))
        self.wires[name] = [self._new() for _ in range(self.m)]
        if self.stream:
            if any(self.C.g(w) for w in self.wires[name]):
                raise CheckError("new register wires not zero")
            self.C.write(self.wires[name], vals)
        return name

    def alloc(self, name):
        self._fresh.add(name)
        return self.reg(name, [0] * self.S)

    def free(self, name):
        if any(self.val[name]):
            raise CheckError("free of a nonzero register %s" % name)
        del self.val[name]
        self.live_regs.discard(name)
        for w in self.wires.pop(name):
            self._release(w)

    # operations
    def lin(self, T, A, k=0):
        """T ^= A^(2^k)  (CNOT only)."""
        F = self.F
        self.val[T] = [t ^ F.frob(a, k) for t, a in zip(self.val[T], self.val[A])]
        self._fresh.discard(T)
        if self.gates:
            wa, wt = self.wires[A], self.wires[T]
            for j, c in enumerate(F.frob_cols(k)):
                for i in range(self.m):
                    if c >> i & 1:
                        self.C.cx(wa[j], wt[i])

    def _prod(self, A, ka, B, kb):
        F = self.F
        return [F.mul(F.frob(a, ka), F.frob(b, kb)) for a, b in zip(self.val[A], self.val[B])]

    def _terms(self, ka, kb, same):
        key = (ka % self.m, kb % self.m, same)
        if key not in self._terms_cache:
            ca, cb = self.F.frob_cols(ka), self.F.frob_cols(kb)
            lin_terms = {}
            ands = []
            for al, be, ga in self.dec:
                a2, b2 = pull_mask(al, ca), pull_mask(be, cb)
                if not a2 or not b2:
                    continue
                if same and a2 == b2:            # (u.x)^2 = u.x : a linear term
                    lin_terms[a2] = lin_terms.get(a2, 0) ^ ga
                else:
                    ands.append((a2, b2, ga))
            self._terms_cache[key] = (ands, [(u, g) for u, g in lin_terms.items() if g])
        return self._terms_cache[key]

    def mul(self, T, A, ka, B, kb):
        """T ^= A^(2^ka) * B^(2^kb)."""
        if T in (A, B):
            raise ValueError("target must differ from the operands")
        p = self._prod(A, ka, B, kb)
        self.val[T] = [t ^ x for t, x in zip(self.val[T], p)]
        self.n_mul += 1
        ands, lins = self._terms(ka, kb, A == B)
        if self.gates and self.crt is not None:
            self.and_count += self._crt_product(T, A, ka, B, kb, "mul")[0]
        else:
            self.and_count += len(ands)
            if self.gates:
                self._lower_mul(T, A, B, ands, lins)
        self._fresh.discard(T)

    def _crt_product(self, T, A, ka, B, kb, mode, mids=None):
        """CRT lowering (mult.emit_product): A and B transformed in place and restored; for A == B the second
        operand is a temporary copy (m extra wires, 2m CX)."""
        from .mult import emit_product
        C, m = self.C, self.m
        wa, wt = self.wires[A], self.wires[T]
        same = A == B
        if same:
            wb = [self._new() for _ in range(m)]
            for i in range(m):
                C.cx(wa[i], wb[i])
        else:
            wb = self.wires[B]
        anc = self._new() if mode == "mul" else None
        ca = self.F.frob_cols(ka) if ka % m else None
        cb = self.F.frob_cols(kb) if kb % m else None
        r = emit_product(C, self.crt, wa, wb, wt, anc, ca, cb, mode=mode, mids=mids, same=same,
                         zero_target=(mode == "mul" and T in self._fresh))
        if anc is not None:
            self._release(anc)
        if same:
            for i in range(m):
                C.cx(wa[i], wb[i])
            for w in wb:
                self._release(w)
        return r

    def erase(self, T, A, ka, B, kb, free=True):
        """measure T away; condition T == A^(2^ka) * B^(2^kb) (checked on the samples).  0 ANDs."""
        p = self._prod(A, ka, B, kb)
        if p != self.val[T]:
            raise CheckError("erasure of %s: value != product" % T)
        self.val[T] = [0] * self.S
        self.n_erase += 1
        ands, lins = self._terms(ka, kb, A == B)
        if self.erase_mode == "terms" and self.gates and self.crt is not None:
            wt = self.wires[T]
            if self.stream:
                self.C.lazy_erase(wt, self.wires[A], ka, self.wires[B], kb, self.F)
                mids = None
            else:
                self.C.check_product(wt, self.wires[A], self.wires[B], ands, lins)
                mids = [self.C.measure_reset(w) for w in wt]
            self.cz_count += self._crt_product(T, A, ka, B, kb, "erase", mids)[1]
            if free:
                self.free(T)
            return
        if self.erase_mode == "terms":
            self.cz_count += sum(bin(g).count("1") for _, _, g in ands) + sum(bin(g).count("1") for _, g in lins)
        elif self.count_cz and (self.stream or not self.gates):
            self.cz_count += self._quad_count(ka, kb, A == B)
        if self.stream:
            if self.erase_mode == "terms":          # the parity CX of the correction, executed for the counts
                self._lower_erase_terms(T, A, B, ands, lins, dry=True)
            self.C.lazy_erase(self.wires[T], self.wires[A], ka, self.wires[B], kb, self.F)
        elif self.gates and self.erase_mode == "terms":
            self._lower_erase_terms(T, A, B, ands, lins)
        elif self.gates:
            polys = self._quad_forms(A, ka, B, kb)
            self.cz_count += sum(len(pl) for pl in polys)
            wa, wb, wt = self.wires[A], self.wires[B], self.wires[T]
            for i in range(self.m):
                mons = []
                for mon in polys[i]:
                    if len(mon) == 1:
                        mons.append((wa[mon[0]],))
                    else:
                        mons.append((wa[mon[0]], wb[mon[1]]) if A != B else (wa[mon[0]], wa[mon[1]]))
                self.C.rerase(wt[i], mons)
        if free:
            self.free(T)

    def _lower_erase_terms(self, T, A, B, ands, lins, dry=False):
        """register erasure by the bilinear terms: measure every bit t_i in the X basis (outcome s_i, wire reset by
        a conditional X), then apply (-1)^(s.T) = prod_k (-1)^((s.gamma_k)(alpha_k.A)(beta_k.B)) as one
        classically controlled CZ between the two parity wires per (term k, bit i of gamma_k)  (Z for a linear
        term).  O(rank * m) CZ instead of the O(m^3) monomials of the per-bit quadratic forms.
        dry=True (stream mode): only the parity CX are executed (for the counts); the conditional phases do not
        act on basis states."""
        C = self.C
        wa, wb, wt = self.wires[A], self.wires[B], self.wires[T]
        if not dry:
            C.check_product(wt, wa, wb, ands, lins)
            mids = [C.measure_reset(w) for w in wt]
        for u, g in lins:
            pa, ia = self._par(wa, u)
            if not dry:
                for i in range(self.m):
                    if g >> i & 1:
                        C.cphase(mids[i], (pa,))
            self._unpar(pa, ia, wa)
        for al, be, ga in ands:
            pa, ia = self._par(wa, al)
            pb, ib = self._par(wb, be)
            if not dry:
                for i in range(self.m):
                    if ga >> i & 1:
                        C.cphase(mids[i], (pa, pb))
            self._unpar(pb, ib, wb)
            self._unpar(pa, ia, wa)

    def _quad_count(self, ka, kb, same):
        """number of Z/CZ gates of a register erasure (without building the monomials)."""
        key = ("n", ka % self.m, kb % self.m, same)
        if key not in self._terms_cache:
            F, m = self.F, self.m
            ca, cb = F.frob_cols(ka), F.frob_cols(kb)
            pr = [[F.mul(ca[u], cb[v]) for v in range(m)] for u in range(m)]
            if same:
                n = sum(bin(pr[u][u]).count("1") for u in range(m))
                n += sum(bin(pr[u][v] ^ pr[v][u]).count("1") for u in range(m) for v in range(u + 1, m))
            else:
                n = sum(bin(x).count("1") for row in pr for x in row)
            self._terms_cache[key] = n
        return self._terms_cache[key]

    def _quad_forms(self, A, ka, B, kb):
        """bit i of A^(2^ka) * B^(2^kb) as monomials over register bit indices: (u,) linear (A == B, u == v) or
        (u, v) for a_u b_v (A != B) / a_u a_v with u < v (A == B)."""
        key = ("q", ka % self.m, kb % self.m, A == B)
        if key in self._terms_cache:
            return self._terms_cache[key]
        F, m = self.F, self.m
        ca, cb = F.frob_cols(ka), F.frob_cols(kb)
        cnt = [dict() for _ in range(m)]
        for u in range(m):
            for v in range(m):
                pr = F.mul(ca[u], cb[v])
                if not pr:
                    continue
                if A == B:
                    mon = (u,) if u == v else (min(u, v), max(u, v))
                else:
                    mon = (u, v)
                x = pr
                while x:
                    i = (x & -x).bit_length() - 1
                    x &= x - 1
                    cnt[i][mon] = cnt[i].get(mon, 0) ^ 1
        out = [sorted(k for k, c in d.items() if c) for d in cnt]
        self._terms_cache[key] = out
        return out

    # gate lowering of a product
    def _par(self, wires, mask):
        idx = [j for j in range(self.m) if mask >> j & 1]
        if len(idx) == 1:
            return wires[idx[0]], None
        p = self._new()
        for j in idx:
            self.C.cx(wires[j], p)
        return p, idx

    def _unpar(self, p, idx, wires):
        if idx is None:
            return
        for j in idx:
            self.C.cx(wires[j], p)
        self._release(p)

    def _lower_mul(self, T, A, B, ands, lins):
        C = self.C
        wa, wb, wt = self.wires[A], self.wires[B], self.wires[T]
        for u, g in lins:
            pa, ia = self._par(wa, u)
            for i in range(self.m):
                if g >> i & 1:
                    C.cx(pa, wt[i])
            self._unpar(pa, ia, wa)
        for al, be, ga in ands:
            pa, ia = self._par(wa, al)
            pb, ib = self._par(wb, be)
            y = self._new()
            C.AND(pa, pb, y)
            for i in range(self.m):
                if ga >> i & 1:
                    C.cx(y, wt[i])
            C.erase(pa, pb, y)
            self._release(y)
            self._unpar(pb, ib, wb)
            self._unpar(pa, ia, wa)

    def finish(self):
        if self.gates:
            self.C.nq = self._nw
        return self


class Chain:
    def __init__(self, m, chain="it", p=None, use_p=None):
        self.ch = make_chain(m, chain)
        self.L = len(self.ch) - 1
        self.p = p or self.L
        r = pebble_plan(self.ch, self.p, use_p)
        if r is None:
            raise ValueError("no pebbling strategy with p=%d chain registers" % self.p)
        self.cu, self.cd, self.up, self.down = r

    def run(self, P, X, moves, tag):
        for kind, s in moves:
            e, i, j = self.ch[s]
            ri = X if i == 0 else "%sb%d" % (tag, i)
            rj = X if j == 0 else "%sb%d" % (tag, j)
            rs = "%sb%d" % (tag, s)
            ej = self.ch[j][0]
            if kind == 'P':
                P.alloc(rs)
                P.mul(rs, ri, ej, rj, 0)
            else:
                P.erase(rs, ri, ej, rj, 0)

    def sink(self, tag):
        return "%sb%d" % (tag, self.L)


def divide_program(P, ch, T, X, Y, tag="c"):
    """T ^= Y / X  (x = 0 gives T ^= 0).  (cu + 1 + cd) multiplications."""
    ch.run(P, X, ch.up, tag)
    P.mul(T, ch.sink(tag), 1, Y, 0)
    ch.run(P, X, ch.down, tag)


def inverse_program(P, ch, T, X, tag="c"):
    """T ^= X^-1 (0 -> 0).  (cu + cd) multiplications; the last step is the linear Frobenius map."""
    ch.run(P, X, ch.up, tag)
    P.lin(T, ch.sink(tag), 1)
    ch.run(P, X, ch.down, tag)
