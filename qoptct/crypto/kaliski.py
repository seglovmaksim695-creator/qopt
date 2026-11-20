# SPDX-License-Identifier: Apache-2.0
"""qoptct.crypto.kaliski -- reversible Kaliski modular inversion on logical ANDs, 12n^2 - 8n - 3 ANDs (generic odd prime p).

    C = kaliski_inverse(p, n=None, reduction='generic', clean=False, legacy=False)
    C.to_qasm()            OpenQASM 3 in the qoptct subset (logical AND = 4 T, erase = X measurement + Clifford fix)
    C.counts()             {'AND': ..., 'ERASE': ..., 'qubits': ...}
    C.evaluate(xs)         exact bit-sliced classical run on many inputs (checks every AND target and erase condition)
    C.expected(x)          the function: -x^{-1} * 2^{2n} mod p
    kaliski_and_count(n)   12n^2 - 8n - 3 (generic reduction, clean=False; legacy=True: 12n^2 + n)

FUNCTION. Input x in [1, p) on the wires C.x (little-endian), p an odd prime, n = p.bit_length(). The result
K = -x^{-1} * 2^{2n} mod p (Kaliski's "almost inverse" with the 2^{2n} factor and the sign that Litinski 2023
removes with a final negation) appears on the wires C.out.
  * clean=False (default): forward computation only, as in the resource counts of the literature. The input wires
    are consumed; the registers end as u = 1, v = 0 (C.regs['u'], C.regs['v']), the rest (C.garbage: the
    per-step bits w_k, e_k, the quotient bits of the final reduction, the s register) is garbage that a surrounding
    compute-copy-uncompute removes. ANDs: exactly 12n^2 - 8n - 3 for reduction='generic'; qubits 11n - 2
    for reduction='generic' (legacy=True: 12n^2 + n ANDs, 11n qubits).  reduction='special' changes both
    counts: at n = 256 secp256k1 657 503 ANDs / 3 107 qubits, P-256 664 966 ANDs / 5 513 qubits (generic: 784 381
    ANDs / 2 814 qubits); counted by C.counts().
  * clean=True: forward, copy of K to n fresh wires, inverse of the forward part:  |x>|0> -> |x>|K>, every other
    wire returns to |0>. ANDs = forward ANDs + forward erases (the inverse turns erases into ANDs and vice versa);
    qubits 12n - 2 (legacy: 12n). The wires 0..n-1 are x, n..2n-1 the result.
x = 0 and x >= p are outside the specified domain. (Observed, and proved at n = 3..7 by the tests: the clean circuit
then gives K(x mod p), and 0 for x = 0 mod p, with all ancillas clean. Not guaranteed in general.)

COST MODEL. Toffoli count = number of logical ANDs (Gidney 2018, arXiv:1709.06648: an AND into a |0> wire costs
4 T; its uncomputation by X-basis measurement and a CZ correction costs 0 T). This is the model of Litinski 2023
(arXiv:2306.08585, Fig. 8), whose Kaliski inversion costs 26n^2 + 2n Toffoli (2n steps of 13n).

CONSTRUCTION (2n Kaliski steps on (u, v, r, s), u = p, v = x, r = 0, s = 1):
  1. the step's swap condition is merged with the next step's (the registers are kept in the frame of w_{k-1}, one
     conditional swap per step on d = w_{k-1} xor w_k; the last swap is a relabelling since w_{2n} = 1 on [1, p));
  2. widths are truncated by provable bounds: before step k (while v >= 1) u, v < 2^{2n-1-k}, so u, v fit in
     2n-1-k bits (and in n bits always); r, s <= 2^k.  Proof: every step halves exactly one of the two registers
     (u/2, v/2, (u-v)/2 <= u/2 or (v-u)/2 <= v/2), and a subtraction never increases the operated register.  Let
     h_u, h_v be the number of halvings of u and of v in steps 0..k-1, so h_u + h_v = k.  Since u = p < 2^n and
     v = x < 2^n at the start, u < 2^{n-h_u} and v < 2^{n-h_v}.  While v >= 1 both registers are >= 1 (u >= 1
     always), hence n - h_u >= 1 and n - h_v >= 1, i.e. h_u, h_v <= n-1.  Therefore
     u < 2^{n-h_u} = 2^{n-k+h_v} <= 2^{2n-1-k} and likewise v < 2^{n-k+h_u} <= 2^{2n-1-k}.  (Once v = 0 the
     state is u = 1, v = 0, which fits any width >= 1.)  The single product u*v < 2^{2n-k} alone would not suffice
     (for v = 1 it only gives u < 2^{2n-k});
  3. the doubling of s is a relabelling (s grows to 2n bits); the reduction mod p is done once at the end
     (reduction='generic': long division, n+1 trial subtractions; reduction='special': O(n) folding for primes with
     a short signed-binary 2^m mod p such as pseudo-Mersenne primes, P-256, secp256k1);
  4. (new here) when e = 1 the step has just swapped so that u >= v (both odd), hence bitlen(v) <= floor((2n-k)/2):
     the high partial products e*v_i of "u -= e*v" are identically 0 and are omitted (-0.5 n^2 ANDs).  Proof, with
     h_u + h_v = k as in item 2 (counted per register, the swap only renames them): b = bitlen(v) satisfies
     b <= n - h_v because v < 2^{n-h_v}, and b <= n - h_u because v <= u < 2^{n-h_u}; adding, 2b <= 2n - k, so
     b <= floor((2n-k)/2).  (From v^2 <= u*v < 2^{2n-k} alone one would only get ceil((2n-k)/2) for odd 2n-k.)
  5. (default, legacy=False) AND/erase pairs whose product is identically constant on [1, p) are not emitted.
     They were found automatically (bit-sliced run over all x in [1, p) for every prime with n = 3..8, every AND
     whose product is constant); only the classes below are removed, each with a proof for every odd prime p and
     every n (a constant on the tested primes alone is not enough). Facts used: gcd(u, v) = 1 throughout (p prime,
     1 <= x < p); in step k the operated register is the even one if one is even, else the larger odd one (v on a
     tie), i.e. w_k = 1 iff the step operates on v; the registers are kept so that A is the operated one.
     P1 (r += e*s, carry out of bit 0): after every step one of r, s is the freshly doubled register, whose bit 0 is
        a fresh |0> wire (and r = 0 before step 0); the swap exchanges r_0 and s_0 together, so r_0 * s_0 = 0 and
        e*s_0*r_0 = 0: one AND + erase per step (2n - 1).
     P2 (u -= e*v, borrow out of bit 0): e is a copy of u_0, the low bit of e*v is e (v is odd after the swap),
        so bit 0 becomes u_0 ^ e = 0 and the borrow (1 ^ u_0)*e is 0; the borrow chain starts at bit 1 (2n - 2).
     P3 (swap of (u, v), bit 0, steps k >= 1): the register B not operated on in step k-1 is odd (the even one is
        always the operated one, and both even is excluded by gcd = 1; v = 0 forces u = 1). If A is even it is
        operated again, d = 0; so d = 1 implies A odd, A_0 = B_0 = 1 and the bit-0 swap is the identity (2n - 2).
     P4 (step 2n-1 is empty): every step at least halves u or v (u/2, v/2, (u-v)/2 < u/2, (v-u)/2 < v/2); u stays
        >= 1 and v stays >= 1 until the single step (1, 1) -> (1, 0) (v becomes 0 only from u = v, i.e. u = v = 1).
        So u is halved at most floor(log2 p) = n-1 times and v at most floor(log2 x) <= n-1 times before that step:
        v = 0 after at most 2n-1 steps. Hence at step 2n-1 always (u, v) = (1, 0): w_{2n-1} = 1, e = 0, and step
        2n-2 is either (1, 0) or the step (1, 1) -> (1, 0): w_{2n-2} = 1, so d_{2n-1} = 0. Step 2n-1 reduces to the
        relabellings u /= 2, s *= 2 (3n + 4 ANDs fewer before P1-P3; its w and e wires vanish: -2 qubits).
     P5 (comparison of step 2n-2): w_{2n-2} = 1 (P4), so it is not computed; d = not w_{2n-3} (4 ANDs: g and w
        have a constant product, w_{k-1}(u_0 ^ v_0) and the comparison carry do not but only feed w).
     Total -(9n + 3) ANDs. Found but NOT removed (no proof for all p, or p-dependent; O(1) or a few per n): the
     long-division bits j = 0, 1 (r_ext = 0 mod 4 would allow it, 2 ANDs), single partial products and carries at
     specific steps near k = n, 2n - 4 (bounds sharper than those of item 2), the high borrows of u -= e*v at the
     truncation boundary, and the p-dependent constants of the special reduction and of the constant additions.
     legacy=True keeps the old circuit (12n^2 + n, byte-identical OpenQASM to the lab circuits k4).
The construction and the count were developed in the qopt lab (agent k4, on an earlier version a9 with 12.5n^2+0.5n)
and audited independently (lab audits k10, x5) by separate bit-sliced simulators.

VERIFIED (tests/test_crypto_kaliski.py): all x in [1, p) for several primes at n = 3..8 (both legacy and default),
the AND count formulas, exact proofs with qoptct.verify.prove of the clean circuits at p = 5, 7, 13, 31 (slow: 61,
127), rejected mutants. At larger n
(lab runs up to n = 32 on 305 inputs per prime; n = 256 for special primes on a few hundred inputs) the statement is
tested, not proved. The counts 12n^2 + n (legacy) and 12n^2 - 8n - 3 are exact for every n and p checked; they are
properties of the construction, not theorems proved here for all n (the removal P1-P5 itself is proved above).
"""
__all__ = ['kaliski_inverse', 'kaliski_and_count', 'KaliskiCircuit', 'CheckError', 'REDUCTIONS']

REDUCTIONS = ('generic', 'special')


from ._errors import CheckError  # noqa: E402  (shared with qoptct.crypto.mulsq)


def kaliski_and_count(n, legacy=False):
    """ANDs of the forward circuit with reduction='generic': 12n^2 - 8n - 3 (legacy=True: 12n^2 + n)."""
    return 12 * n * n + n if legacy else 12 * n * n - 8 * n - 3


# ================================================================== builder (logical AND + measurement-based erase)
# ops: ('x', q) ('cx', c, t) ('and', a, b, t)  t := a*b, t must be 0
#      ('erase', t, a, b, lin)   t must equal a*b xor XOR(lin); X measurement, on outcome 1: cz a b; z w (w in lin); x t

_ZERO, _ONE = ('c', 0), ('c', 1)


def _W(w):
    return ('w', w)


class _Circ:
    def __init__(self):
        self.ops = []
        self.nw = 0
        self.free = []
        self.marks = []

    def mark(self, label, **info):
        """untrusted segment hint for qoptct.verify.compose: (ops so far, label, info). Does not change the ops."""
        self.marks.append((len(self.ops), label, info))

    def reg(self, n):
        out = list(range(self.nw, self.nw + n))
        self.nw += n
        return out

    def alloc(self):
        if self.free:
            return self.free.pop()
        self.nw += 1
        return self.nw - 1

    def release(self, w):
        assert w not in self.free
        self.free.append(w)

    def x(self, q):
        self.ops.append(('x', q))

    def cx(self, c, t):
        assert c != t
        self.ops.append(('cx', c, t))

    def AND(self, a, b, t=None):
        if t is None:
            t = self.alloc()
        assert len({a, b, t}) == 3, (a, b, t)
        self.ops.append(('and', a, b, t))
        return t

    def erase(self, t, a, b, lin=(), release=True):
        assert t not in (a, b) and t not in lin and a != b
        self.ops.append(('erase', t, a, b, tuple(lin)))
        if release:
            self.release(t)


def _inverse_ops(ops):
    out = []
    for o in reversed(ops):
        k = o[0]
        if k in ('x', 'cx'):
            out.append(o)
        elif k == 'and':
            out.append(('erase', o[3], o[1], o[2], ()))
        elif k == 'erase':
            t, a, b, lin = o[1:]
            out.append(('and', a, b, t))
            for w in lin:
                out.append(('cx', w, t))
        else:
            raise ValueError(k)
    return out


def _simulate(ops, nw, init, S):
    """bit-sliced classical run: init {wire: int mask over S samples}; checks every AND target and erase condition."""
    M = (1 << S) - 1
    v = [0] * nw
    for w, m in init.items():
        v[w] = m
    for i, o in enumerate(ops):
        k = o[0]
        if k == 'x':
            v[o[1]] ^= M
        elif k == 'cx':
            v[o[2]] ^= v[o[1]]
        elif k == 'and':
            if v[o[3]]:
                raise CheckError('op %d: AND target not zero' % i)
            v[o[3]] = v[o[1]] & v[o[2]]
        elif k == 'erase':
            t, a, b, lin = o[1:]
            val = v[a] & v[b]
            for w in lin:
                val ^= v[w]
            if v[t] != val:
                raise CheckError('op %d: erase condition violated' % i)
            v[t] = 0
        else:
            raise ValueError(k)
    return v


def _pack(vals, wires):
    m = [0] * len(wires)
    for s, x in enumerate(vals):
        for i in range(len(wires)):
            if (x >> i) & 1:
                m[i] |= 1 << s
    return dict(zip(wires, m))


def _unpack(v, wires, S):
    out = [0] * S
    for i, w in enumerate(wires):
        m = v[w]
        for s in range(S):
            if (m >> s) & 1:
                out[s] |= 1 << i
    return out


# ------------------------------------------------------------------ carries (Gidney adder) with constant propagation
class _Rec:
    """recipe of a carry bit MAJ(x, y, c): kind op (alias / constant), and, or, gid."""
    def __init__(self, kind, *args):
        self.kind = kind
        self.args = args
        self.t = None
        self.fold = False


def _plan(x, y, c):
    ops = [x, y, c]
    for i in range(3):
        for j in range(i + 1, 3):
            if ops[i] == ops[j]:
                return _Rec('op', ops[i])
    cs = [o for o in ops if o[0] == 'c']
    if len(cs) >= 2:
        return _Rec('op', [o for o in ops if o[0] == 'w'][0] if len(cs) == 2 else _ZERO)
    if len(cs) == 1:
        u, v = [o[1] for o in ops if o[0] == 'w']
        return _Rec('and' if cs[0][1] == 0 else 'or', u, v)
    return _Rec('gid', x[1], y[1], c[1])


def _compute(C, r, fold_ok=False, t=None):
    k = r.kind
    if k == 'op':
        o = r.args[0]
        if t is None:
            return o
        if o == _ONE:
            C.x(t)
        elif o[0] == 'w':
            C.cx(o[1], t)
        r.t = t
        return _W(t)
    if k == 'and':
        u, v = r.args
        r.t = C.AND(u, v, t)
        return _W(r.t)
    if k == 'or':
        u, v = r.args
        C.x(u); C.x(v)
        r.t = C.AND(u, v, t)
        C.x(u); C.x(v); C.x(r.t)
        return _W(r.t)
    x, y, c = r.args
    C.cx(c, x); C.cx(c, y)
    r.t = C.AND(x, y, t)
    C.cx(c, y)
    if fold_ok:
        r.fold = True
    else:
        C.cx(c, x)
    C.cx(c, r.t)
    return _W(r.t)


def _uncompute(C, r, t=None, release=True):
    k = r.kind
    t = r.t if t is None else t
    if k == 'op':
        o = r.args[0]
        if o == _ONE:
            C.x(t)
        elif o[0] == 'w':
            C.cx(o[1], t)
        if release:
            C.release(t)
        return
    if k == 'and':
        u, v = r.args
        C.erase(t, u, v, release=release)
        return
    if k == 'or':
        u, v = r.args
        C.x(t); C.x(u); C.x(v)
        C.erase(t, u, v, release=release)
        C.x(u); C.x(v)
        return
    x, y, c = r.args
    if not r.fold:
        C.cx(c, x)
    C.cx(c, y)
    C.erase(t, x, y, lin=(c,), release=release)
    C.cx(c, y)
    if not r.fold:
        C.cx(c, x)


def _add(C, a, b, cout=None, lsb_nocarry=False):
    """b += a (mod 2^n); cout (a zero wire) receives the carry out. ANDs: n-1 (n with cout).
    lsb_nocarry: the caller guarantees a[0]*b[0] = 0 on every input, so the carry out of bit 0 is 0 and is not
    computed (one AND and its erase fewer; with n = 1 and cout, cout simply stays 0)."""
    n = len(a)
    c = [None] * (n + 1)
    last = n if cout is not None else n - 1
    for i in range(last):
        if i == 0 and lsb_nocarry:
            continue
        t = cout if i == n - 1 else C.alloc()
        if c[i] is None:
            C.AND(a[i], b[i], t)
        else:
            C.cx(c[i], a[i]); C.cx(c[i], b[i]); C.AND(a[i], b[i], t); C.cx(c[i], t)
        c[i + 1] = t
    for i in reversed(range(n)):
        if i < last:
            if i + 1 < n and c[i + 1] is not None:
                C.erase(c[i + 1], a[i], b[i], lin=(() if c[i] is None else (c[i],)))
            if c[i] is not None:
                C.cx(c[i], a[i])
            C.cx(a[i], b[i])
        else:
            C.cx(a[i], b[i])
            if c[i] is not None:
                C.cx(c[i], b[i])


def _add_gen(C, Y, b):
    """b += Y mod 2^len(b); Y[i] an operand ('w', wire) / _ZERO / _ONE; missing high operands are _ZERO."""
    m = len(b)
    Y = list(Y) + [_ZERO] * (m - len(Y))
    X = [_W(w) for w in b]
    c = [_ZERO]
    R = [None]
    for i in range(m - 1):
        r = _plan(X[i], Y[i], c[i])
        c.append(_compute(C, r, fold_ok=True))
        R.append(r)
    for i in reversed(range(m)):
        if i + 1 < m and R[i + 1].kind != 'op':
            _uncompute(C, R[i + 1])
        if Y[i] == _ONE:
            C.x(b[i])
        elif Y[i][0] == 'w':
            C.cx(Y[i][1], b[i])
        folded = (i + 1 < m) and R[i + 1].kind == 'gid' and R[i + 1].fold
        if not folded:
            if c[i] == _ONE:
                C.x(b[i])
            elif c[i][0] == 'w':
                C.cx(c[i][1], b[i])


def _sub_gen(C, Y, b):
    """b -= Y mod 2^len(b)."""
    for w in b:
        C.x(w)
    _add_gen(C, Y, b)
    for w in b:
        C.x(w)


def _add_k(C, k, b, ctrl=None):
    """b += k (or ctrl*k), k a classical constant, mod 2^len(b)."""
    m = len(b)
    Y = [(_W(ctrl) if ctrl is not None else _ONE) if (k >> i) & 1 else _ZERO for i in range(m)]
    _add_gen(C, Y, b)


def _carry_chain(C, X, Y, cin=_ZERO):
    m = len(X)
    c = [cin]
    R = [None]
    for i in range(m - 1):
        r = _plan(X[i], Y[i], c[i])
        c.append(_compute(C, r, fold_ok=False))
        R.append(r)
    return c, R


def _chain_undo(C, R):
    for r in reversed(R[1:]):
        if r.kind != 'op':
            _uncompute(C, r)


def _compute_cout(C, X, Y, f, cin=_ZERO):
    """f (zero wire) := carry out of X + Y + cin."""
    m = len(X)
    c, R = _carry_chain(C, X, Y, cin)
    rf = _plan(X[m - 1], Y[m - 1], c[m - 1])
    _compute(C, rf, t=f)
    _chain_undo(C, R)


def _cswap(C, d, A, B):
    for a, b in zip(A, B):
        C.cx(a, b)
        t = C.AND(d, b)
        C.cx(t, a)
        C.erase(t, d, b)
        C.cx(a, b)


# ------------------------------------------------------------------ reduction mod a special prime (O(n) ANDs)
def _naf(c):
    out = []
    e = 0
    while c:
        if c & 1:
            d = 2 - (c & 3)
            out.append((d, e))
            c -= d
        c >>= 1
        e += 1
    return out


def _cmp_ge_const(C, A, K):
    """f := [A >= K] via the carry of A + (2^L - K); None if K >= 2^L (then f = 0)."""
    L = len(A)
    if K >= (1 << L):
        return None
    f = C.alloc()
    Y = [_ONE if (((1 << L) - K) >> i) & 1 else _ZERO for i in range(L)]
    _compute_cout(C, [_W(w) for w in A], Y, f)
    return f


def _red_special(C, V, p, Vbound):
    """V (wires, value < Vbound) -> (out wires = V mod p, garbage). Rounds V = H*2^m + L == K + L + sum s_k (H << e_k)
    with 2^m mod p = sum s_k 2^{e_k} (NAF) and K a multiple of p that keeps the value non-negative; then a binary
    search of the quotient (f = [A >= 2^j p], A -= f 2^j p). Correct for every odd p; O(n) ANDs when the NAF is short."""
    m = p.bit_length()
    terms = _naf((1 << m) % p)
    gar = []
    cur, B = list(V), Vbound
    while B >= (1 << (m + 2)) and len(cur) > m:
        Lw, Hw = cur[:m], cur[m:]
        Hmax = (B - 1) >> m
        Hw = Hw[:Hmax.bit_length()]
        neg = sum((Hmax << e) for s, e in terms if s < 0)
        pos = sum((Hmax << e) for s, e in terms if s > 0)
        K = -(-neg // p) * p
        Bn = K + (1 << m) + pos
        A = [C.alloc() for _ in range(Bn.bit_length())]
        for i, w in enumerate(A):
            if (K >> i) & 1:
                C.x(w)
        C.mark('sp_const', A=list(A), K=K)
        _add_gen(C, [_W(w) for w in Lw], A)
        C.mark('sp_addL')
        for s, e in sorted(terms, key=lambda t: -t[0]):     # positive terms first
            if e >= len(A):
                continue
            Y = [_W(w) for w in Hw][:len(A) - e]
            (_add_gen if s > 0 else _sub_gen)(C, Y, A[e:])
            C.mark('sp_term', s=s, e=e)
        gar += cur
        cur, B = A, Bn
    A = cur
    j = ((B - 1) // p).bit_length() - 1
    while j >= 0:
        K = p << j
        f = _cmp_ge_const(C, A, K)
        if f is not None:
            _add_k(C, (1 << len(A)) - K, A, ctrl=f)
            gar.append(f)
            C.mark('sp_q', f=f, K=K)
        j -= 1
    return A[:m], gar + A[m:]


# ------------------------------------------------------------------ the Kaliski forward circuit
def _kaliski_forward(C, V, p, reduction, legacy=False):
    """V: n input wires (x in [1, p)). -> (u, v, out, s, garbage) wires; out = -x^{-1} 2^{2n} mod p.
    legacy=True: the 0.5.0rc construction (12n^2 + n ANDs). legacy=False: the same circuit without the AND/erase
    pairs whose product is identically constant on [1, p) for every p (proofs in the module docstring, items P1-P5)."""
    n = len(V)
    new = not legacy
    U = [C.alloc() for _ in range(n)]
    R = [C.alloc() for _ in range(n)]
    S = [C.alloc() for _ in range(n)]
    for i in range(n):
        if (p >> i) & 1:
            C.x(U[i])
    C.x(S[0])
    C.mark('init', U=list(U), R=list(R), S=list(S))
    A, B, Rr, Ss = U, list(V), R, S
    wprev = None
    garbage = []
    for k in range(2 * n):
        if new and k == 2 * n - 1:
            # P4: at step 2n-1 always v = 0, u = 1 and w_{2n-2} = w_{2n-1} = 1, so d = 0 (no swap) and e = 0
            # (nothing is subtracted or added): the step is the two relabellings u /= 2 and s *= 2 only
            A = A[1:] + A[:1]
            Ss = [C.alloc()] + Ss
            C.mark('last', k=k, snew=Ss[0])
            if len(Ss) > 2 * n:
                C.release(Ss.pop())
            continue
        # widths (docstring, item 2): u < 2^(n-h_u), v < 2^(n-h_v), h_u + h_v = k, h_u, h_v <= n-1 while v >= 1
        # => u, v < 2^(2n-1-k) (v = 0 only with u = 1);  max(r, s) <= 2^k
        wa = min(n, max(1, 2 * n - 1 - k))
        wr = min(n, k + 1)
        wr2 = min(n, k + 2)
        if new and k == 2 * n - 2:
            # P5: w_{2n-2} = 1 for every x in [1, p): no comparison; d = w_{2n-3} ^ 1 is the wire w_{2n-3} under X
            d = wprev
            C.x(d)
            _cswap(C, d, A[1:wa], B[1:wa])                                 # P3: bit 0 of (u, v) needs no swap
            _cswap(C, d, Rr[:wr], Ss[:wr])
            C.x(d)
            C.mark('swap', k=k)
            w = None
        else:
            w = _compare_and_swap(C, k, n, A, B, Rr, Ss, wa, wr, wprev, new)
        # (c) e = u'_0
        e = C.alloc(); C.cx(A[0], e)
        # (d) u -= e*v; bit 0 of e*v is e. If e = 1 the step swapped to u >= v, so bitlen(v) <= (2n-k)//2
        #     (docstring, item 4: bitlen(v) <= n - h_v and <= n - h_u, summed 2*bitlen(v) <= 2n - k):
        #     the higher partial products are 0 and omitted
        wb = min(wa, max(1, (2 * n - k) // 2))
        m = [e] + [C.AND(e, B[i]) for i in range(1, wb)]
        if new:
            # P2: e is a copy of u_0, so bit 0 is u_0 - e = 0 with no borrow ((1 ^ u_0) * e = 0): the borrow chain
            # starts at bit 1
            C.cx(e, A[0])
            if wa > 1:
                _sub_gen(C, [_W(t_) for t_ in m[1:]], A[1:wa])
        else:
            _sub_gen(C, [_W(t_) for t_ in m], A[:wa])
        for i in reversed(range(1, wb)):
            C.erase(m[i], e, B[i])
        C.mark('sub', k=k, e=e)
        # (e) u /= 2: relabelling
        A = A[1:] + A[:1]
        # (f) r += e*s  (r, s <= 2^k -> the sum is < 2^(k+2)); P1: r_0 * s_0 = 0, no carry out of bit 0
        m = [C.AND(e, s_) for s_ in Ss[:wr]]
        if wr2 > wr:
            _add(C, m, Rr[:wr], cout=Rr[wr], lsb_nocarry=new)
        else:
            _add(C, m, Rr[:wr], lsb_nocarry=new)
        for i in reversed(range(wr)):
            C.erase(m[i], e, Ss[i])
        # (g) s *= 2: relabelling, the register grows to 2n bits (r_ext < 2^(2n): the top bit is provably 0)
        Ss = [C.alloc()] + Ss
        C.mark('add', k=k, snew=Ss[0])
        if len(Ss) > 2 * n:
            C.release(Ss.pop())
        garbage += [e] if w is None else [w, e]
        wprev = w
    if legacy:
        # w_{2n-1} = 1 for every x in [1, p): the last swap is a relabelling
        C.x(wprev); C.release(wprev); garbage.remove(wprev)
    A, B, Rr, Ss = B, A, Ss, Rr
    if reduction == 'none':
        # internal (qoptct.crypto.ecpoint): no final reduction, the raw 2n-bit register r_ext < 2^(2n) with
        # r_ext = -x^{-1} 2^(2n) (mod p) is returned; the caller reduces it (e.g. with qoptct.crypto.fold)
        return A, B, Rr, Ss, garbage
    if reduction == 'special':
        out, g2 = _red_special(C, Rr, p, 1 << (2 * n))
        return A, B, out, Ss, garbage + g2
    # reduction of r_ext (< 2^(2n)) mod p by long division
    R = Rr + [C.alloc() for _ in range(2 * n + 1 - len(Rr))]
    C.mark('red_init', ext=list(R[len(Rr):]))
    for j in range(n, -1, -1):
        win = R[j:j + n + 1]
        _add_k(C, (1 << (n + 1)) - p, win)
        f = C.alloc(); C.cx(win[n], f); garbage.append(f)
        _add_k(C, p, win[:n], ctrl=f)
        C.cx(f, win[n])
        C.mark('div', j=j, f=f)
    for w_ in R[n:]:
        C.release(w_)
    return A, B, R[:n], Ss, garbage


def _compare_and_swap(C, k, n, A, B, Rr, Ss, wa, wr, wprev, new):
    """steps (a) and (b) of Kaliski step k: w_k into a fresh wire, then the swap of (A, B) and (Rr, Ss) on
    d = w_{k-1} ^ w_k (in place on the wire lists). -> the w_k wire."""
    # (a) w_k = [v even] or ([u odd] and [u > v]) in the standard frame; the registers are in the frame of w_{k-1}
    U0 = C.alloc(); V0 = C.alloc()
    C.cx(A[0], U0); C.cx(B[0], V0)
    if wprev is not None:
        C.cx(A[0], B[0]); t = C.AND(wprev, B[0]); C.cx(A[0], B[0])     # t = w (A0 ^ B0)
        C.cx(t, U0); C.cx(t, V0)                                       # U0, V0: low bits of u, v
        for a, b in zip(A[:wa], B[:wa]):
            C.cx(wprev, a); C.cx(wprev, b)
    for a in A[:wa]:
        C.x(a)
    X_, Y_ = [_W(a) for a in A[:wa]], [_W(b) for b in B[:wa]]
    ch, RR = _carry_chain(C, X_, Y_, _ONE)
    rf = _plan(X_[wa - 1], Y_[wa - 1], ch[wa - 1])
    cmp = _compute(C, rf, t=C.alloc())[1]                              # [u <= v]
    C.x(cmp); g = C.AND(V0, cmp); C.x(cmp)                             # g = V0 & ~cmp
    C.x(g); w = C.AND(U0, g); C.x(g)                                   # w = U0 & ~g
    C.x(cmp); C.erase(g, V0, cmp); C.x(cmp)
    _uncompute(C, rf, t=cmp)
    _chain_undo(C, RR)
    for a in A[:wa]:
        C.x(a)
    if wprev is not None:
        for a, b in zip(A[:wa], B[:wa]):
            C.cx(wprev, a); C.cx(wprev, b)
        C.cx(t, U0); C.cx(t, V0)
        C.cx(A[0], B[0]); C.erase(t, wprev, B[0]); C.cx(A[0], B[0])
    C.cx(A[0], U0); C.cx(B[0], V0); C.release(V0); C.release(U0)
    C.mark('cmp', k=k, w=w)
    # (b) merged conditional swap on d = w_{k-1} ^ w_k
    if wprev is None:
        d = w
    else:
        d = C.alloc(); C.cx(wprev, d); C.cx(w, d)
    _cswap(C, d, A[1:wa] if (new and wprev is not None) else A[:wa],     # P3: bit 0 needs no swap
           B[1:wa] if (new and wprev is not None) else B[:wa])
    _cswap(C, d, Rr[:wr], Ss[:wr])
    if wprev is not None:
        C.cx(wprev, d); C.cx(w, d); C.release(d)
    C.mark('swap', k=k)
    return w


# ------------------------------------------------------------------ output formats
_AND_G = lambda a, b, t: [('h', (t,)), ('t', (t,)), ('cx', (b, t)), ('tdg', (t,)), ('cx', (a, t)), ('t', (t,)),
                          ('cx', (b, t)), ('tdg', (t,)), ('h', (t,)), ('sdg', (t,))]


def _g(name, qs):
    return '%s %s;' % (name, ', '.join('q[%d]' % x for x in qs))


def _qasm_lines(ops):
    nb = 0
    for o in ops:
        k = o[0]
        if k == 'x':
            yield _g('x', (o[1],))
        elif k == 'cx':
            yield _g('cx', o[1:])
        elif k == 'and':
            for a, b in _AND_G(*o[1:]):
                yield _g(a, b)
        elif k == 'erase':
            t, a, b, lin = o[1:]
            yield _g('h', (t,))
            yield 'c[%d] = measure q[%d];' % (nb, t)
            body = [_g('cz', (a, b))] + [_g('z', (w,)) for w in lin] + [_g('x', (t,))]
            yield 'if (c[%d]) {%s}' % (nb, ' '.join(body))
            nb += 1


def _qasm_header(ops, nw):
    nb = sum(1 for o in ops if o[0] == 'erase')
    return 'OPENQASM 3;\ninclude "stdgates.inc";\nqubit[%d] q;\nbit[%d] c;\n' % (nw, max(nb, 1))


def _to_qasm(ops, nw):
    return _qasm_header(ops, nw) + '\n'.join(_qasm_lines(ops)) + '\n'


def _write_qasm(ops, nw, path, chunk=200000):
    """the text of _to_qasm written to `path` in chunks (large circuits: no text of hundreds of MB in memory)."""
    with open(path, 'w') as f:
        f.write(_qasm_header(ops, nw))
        buf = []
        for ln in _qasm_lines(ops):
            buf.append(ln)
            if len(buf) >= chunk:
                f.write('\n'.join(buf) + '\n')
                buf = []
        if buf:
            f.write('\n'.join(buf) + '\n')


def _to_qc(ops, nw, inputs):
    """Toffoli-level .qc of the same operations (and -> tof, erase -> cnots from lin + tof)."""
    nm = ['q%d' % i for i in range(nw)]
    L = ['.v ' + ' '.join(nm), '.i ' + ' '.join(nm[i] for i in inputs), '', 'BEGIN']
    for o in ops:
        k = o[0]
        if k == 'x':
            L.append('tof %s' % nm[o[1]])
        elif k == 'cx':
            L.append('tof %s %s' % (nm[o[1]], nm[o[2]]))
        elif k == 'and':
            L.append('tof %s %s %s' % (nm[o[1]], nm[o[2]], nm[o[3]]))
        elif k == 'erase':
            t, a, b, lin = o[1:]
            for w in lin:
                L.append('tof %s %s' % (nm[w], nm[t]))
            L.append('tof %s %s %s' % (nm[a], nm[b], nm[t]))
    L.append('END')
    return '\n'.join(L) + '\n'


# ------------------------------------------------------------------ public API
class KaliskiCircuit:
    """A built Kaliski inversion circuit (see the module docstring). Attributes:
    p, n, reduction, clean; ops (list of logical operations); nqubits; x (input wires, little-endian);
    out (result wires); regs (dict of named final registers, forward mode: 'u' = 1, 'v' = 0, 's');
    garbage (forward mode: wires left dirty; clean mode: [])."""

    def __init__(self, p, n, reduction, clean, ops, nqubits, x, out, regs, garbage):
        self.p, self.n, self.reduction, self.clean = p, n, reduction, clean
        self.ops, self.nqubits, self.x, self.out = ops, nqubits, x, out
        self.regs, self.garbage = regs, garbage

    def __repr__(self):
        c = self.counts()
        return 'KaliskiCircuit(p=%d, n=%d, reduction=%r, clean=%r, AND=%d, qubits=%d)' % (
            self.p, self.n, self.reduction, self.clean, c['AND'], c['qubits'])

    def expected(self, x):
        """-x^{-1} * 2^{2n} mod p (x in [1, p))."""
        return (-pow(x, -1, self.p) * pow(2, 2 * self.n, self.p)) % self.p

    def counts(self):
        """{'AND': logical ANDs (= Toffoli count), 'ERASE': measurement-based uncomputations, 'X', 'CX', 'qubits'}.
        The T count of to_qasm() is 4*AND."""
        k = {}
        for o in self.ops:
            k[o[0]] = k.get(o[0], 0) + 1
        return dict(AND=k.get('and', 0), ERASE=k.get('erase', 0), X=k.get('x', 0), CX=k.get('cx', 0),
                    qubits=self.nqubits)

    def evaluate(self, xs, return_wires=False):
        """Exact classical run (bit-sliced over all xs at once) with a check of every AND target and every erase
        condition (raises CheckError). -> list of the values on C.out; with return_wires=True also the final
        bit-sliced wire values (list of int masks, one bit per input)."""
        xs = list(xs)
        v = _simulate(self.ops, self.nqubits, _pack(xs, self.x), len(xs))
        got = _unpack(v, self.out, len(xs))
        return (got, v) if return_wires else got

    def read(self, wire_values, wires, count):
        """values of a register from the bit-sliced result of evaluate(..., return_wires=True)."""
        return _unpack(wire_values, wires, count)

    def to_qasm(self, path=None):
        """OpenQASM 3 (qoptct subset): AND = 4 T (h t cx tdg cx t cx tdg h sdg), erase = `h; c[k] = measure;
        if (c[k]) {cz; z..; x}`. Wires 0..n-1 are x (clean mode: n..2n-1 are the result)."""
        txt = _to_qasm(self.ops, self.nqubits)
        if path:
            with open(path, 'w') as f:
                f.write(txt)
        return txt

    def write_qasm(self, path):
        """the same text as to_qasm(), streamed to `path` (for large n: the text is never held in memory)."""
        _write_qasm(self.ops, self.nqubits, path)

    def hints(self):
        """untrusted certificate for qoptct.verify.compose (segment boundaries and wire roles of the forward part)."""
        return dict(p=self.p, n=self.n, reduction=self.reduction, clean=self.clean, nqubits=self.nqubits,
                    x=list(self.x), out=list(self.out), regs={k: list(v) for k, v in self.regs.items()},
                    fwd_out=list(getattr(self, 'fwd_out', self.out)), fwd_len=getattr(self, 'fwd_len', len(self.ops)),
                    marks=[[pos, label, info] for pos, label, info in getattr(self, 'marks', [])])

    def to_qc(self, path=None):
        """Toffoli-level .qc of the same operations (inputs: the x wires). A reference for the measurement-based
        implementation, not an independent specification of the function."""
        txt = _to_qc(self.ops, self.nqubits, self.x)
        if path:
            with open(path, 'w') as f:
                f.write(txt)
        return txt


def kaliski_inverse(p, n=None, reduction='generic', clean=False, legacy=False):
    """Build the Kaliski inversion circuit for an odd prime p (see the module docstring).

    p: odd prime >= 3 (primality is not checked; only primes are tested).
    n: register width; must be p.bit_length() (the default; accepted for explicitness).
    reduction: 'generic' (long division, exactly 12n^2 - 8n - 3 ANDs forward; legacy: 12n^2 + n) or 'special' (O(n)
    folding; cheaper for primes with a short signed-binary 2^n mod p, e.g. pseudo-Mersenne 2^n - c, P-256, secp256k1).
    clean: False -> forward circuit with garbage (the literature's count); True -> |x>|0> -> |x>|K> with all other
    wires restored to |0>.  K = -x^{-1} * 2^{2n} mod p for x in [1, p).
    legacy: True -> the 0.5.0rc circuit (12n^2 + n ANDs, byte-identical OpenQASM to the lab k4 circuits); False
    (default) -> without the AND/erase pairs of identically constant product (P1-P5 in the module docstring).
    -> KaliskiCircuit."""
    if not isinstance(p, int) or p < 3 or p % 2 == 0:
        raise ValueError('p must be an odd integer >= 3')
    if reduction not in REDUCTIONS:
        raise ValueError('reduction must be one of %s' % (REDUCTIONS,))
    nb = p.bit_length()
    if n is None:
        n = nb
    if n != nb:
        # the final long division (n+1 trial subtractions) and the special reduction assume 2^(n-1) < p < 2^n
        raise ValueError('n must be p.bit_length() = %d (got %r)' % (nb, n))
    C = _Circ()
    X = C.reg(n)
    O = C.reg(n) if clean else None
    s = len(C.ops)
    u, v, out, sreg, gar = _kaliski_forward(C, X, p, reduction, legacy=bool(legacy))
    if not clean:
        K = KaliskiCircuit(p, n, reduction, False, C.ops, C.nw, X, list(out),
                           dict(u=list(u), v=list(v), s=list(sreg)), list(gar))
    else:
        e = len(C.ops)
        for a, b in zip(out, O):
            C.cx(a, b)
        C.ops += _inverse_ops(C.ops[s:e])
        K = KaliskiCircuit(p, n, reduction, True, C.ops, C.nw, X, O, {}, [])
        K.fwd_out, K.fwd_len = list(out), e
    K.legacy = bool(legacy)
    K.marks = list(C.marks)
    return K
