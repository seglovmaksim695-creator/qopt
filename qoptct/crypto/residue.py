# SPDX-License-Identifier: Apache-2.0
"""residue: compression of a wide accumulator modulo a classical q on logical ANDs ("loop2" of the residue arithmetic
of Gidney 2025, arXiv:2505.15917, App. A.1), and the width of the accumulator that feeds it (ported from the qoptct
lab, agents k2 / k3; audited by k10 and x5).

    F = fold_accumulator(q, clen, nbits, w=None)     -> ResidueFold
    F.qasm()            OpenQASM 3: wires 0..nbits-1 = S (input), then the garbage registers, then ancillas
    F.n_and, F.counts() AND count (= T/4), cost-model counts
    F.reference(S)      expected values of all registers (independent python-int formulas)
    F.check(values)     exact bit-sliced simulation of the gate list on the given inputs (andcirc.simulate)
    fold_accumulator_ands(q, clen, nbits, w)         AND count only
    accumulator_bits(terms, value_bits)              width of the accumulator for a sum of `terms` values < 2^value_bits
    accumulator_bound(terms, value_bits)             the sum bound terms * (2^value_bits - 1) itself

FUNCTION. S is an nbits-bit register (S < 2^nbits, any S is allowed). Let q be a classical integer with
q <= 2^clen (typically q = p - 1 for a clen-bit prime p, the exponent modulus of Gidney 2025). The circuit leaves on
wires 0..clen-1 a value r with

        r < 2^clen      and      r = S (mod q)            (NOT necessarily r < q: one final subtraction of q is saved),

and the high part H = S >> clen on wires clen..nbits-1 (unchanged) plus two small garbage registers: top (e bits) and
g (1 bit). Every other wire ends at |0>. The garbage is removed by running the block backwards (the paper's unloop2).

CONSTRUCTION ("fold1"). Write S = H 2^clen + L, split H into J = ceil((nbits-clen)/w) windows H_j of w bits.
  1. acc = L (+ e fresh wires): for each window, acc += T_j[H_j], T_j[a] = a 2^(clen+wj) mod q  -- one table lookup
     (ANF / Davio: 2^w - w - 1 ANDs) + a Gidney adder + measurement-based erase of the lookup output. Then
     acc = L + sum_j T_j[H_j] < 2^clen + J (q-1) and acc = S (mod q). e = the least number with q 2^e > that bound;
     the wires clen..clen+e-1 of acc are the garbage register 'top'.
  2. u = (acc mod 2^clen) + T[top] with T[a] = a 2^clen mod q (< q), one more lookup + adder into clen+1 bits;
     u < 2^clen + q, the carry g = [u >= 2^clen] (garbage wire 'g').
  3. if g: add k = 2^clen - q to the low clen bits (controlled addition of a constant, no carry): the value drops by q
     and, because then u - 2^clen < q, stays below 2^clen. If g = 0, u < 2^clen already.
  Result r < 2^clen and r = S (mod q).
  Why not "r < q": the paper's loop3 needs only the exponent modulo p-1 (a representative < 2^clen is enough), so the
  fully reduced value is not needed; this is also why the block is cheaper than binary division (App. A.1).

AND COUNTS (generated gate lists; the tests pin them to the lab claims of agent k2):
    q = 1574500 (a 21-bit p - 1), clen = 21, nbits = 32:  w = 3: 148,  w = 4: 138     (formula of the paper, Table 3: 682)
    same q, clen = 21, nbits = 29 (S < 2^29 promised -- always true, see accumulator_bits): w = 4: 106
  w = None chooses the window with the fewest ANDs among 2..6.

HONEST NOTES.
  * Compatibility with the NEXT step of the paper (the loop3 that uses r as an exponent) is argued in words (it only
    needs r = S mod q and r < 2^clen, which hold exactly); it was NOT verified end to end by a circuit.
  * The comparison "148 / 106 vs 682" is against the FORMULA of the paper (2 additions x (l + len m) x len m); the exact
    AND count of the paper's own scheme (binary division, App. A.1) is 352-418 by the lab count, so against it the
    gain is about -1.6 % of the whole RSA-2048 estimate instead of -3.1 %.
  * Proof status: the emitted circuits are proved by qoptct kprove (all inputs, all measurement outcomes) against a
    word-level specification for small analogues (see the tests) and checked by exact bit-sliced simulation (ALL inputs
    for small analogues, random + edge inputs at the size of the paper). The agent-k2 attempt at clen = 21 with kprove
    gave REFUTED with witnesses that contradict the specification (a limit of the tool, not a counterexample):
    do not use a kprove REFUTED at that size without re-checking the witness.
  * Cost model: AND = 4 T + 3 CX, erasure by measurement = 0 T (Gidney 2018).
"""
import random

from .andcirc import Circ, ZERO, simulate, bits_of, ints_of, text_counts

__all__ = ['ResidueFold', 'fold_accumulator', 'fold_accumulator_ands', 'accumulator_bound', 'accumulator_bits',
           'edge_inputs']


# ------------------------------------------------------------------ accumulator width
def accumulator_bound(terms, value_bits):
    """largest value of a sum of `terms` integers, each < 2^value_bits: terms * (2^value_bits - 1)."""
    if terms < 0 or value_bits < 0:
        raise ValueError('terms and value_bits must be >= 0')
    return terms * ((1 << value_bits) - 1)


def accumulator_bits(terms, value_bits):
    """number of bits that holds every sum of `terms` values < 2^value_bits, and no fewer: the bit length of
    terms * (2^value_bits - 1) (0 for an empty sum).

    PROOF. Every sum s satisfies 0 <= s <= B := terms (2^value_bits - 1), so s < 2^b whenever B < 2^b, i.e. for
    b = B.bit_length(); and s = B is attained (all terms maximal), which needs b bits, so b is the least width.
    Example of Gidney 2025 (RSA-2048): Q_dlog is the sum of 214 table values < 2^21:
    214 (2^21 - 1) = 448 790 314 > 2^28 = 268 435 456 and < 2^29, so 29 bits, not the 32 bits of the paper
    (32 = ell + len m is a bound of the paper for another purpose). A register of 29 bits saves 3 wires and,
    in this module, 138 - 106 = 32 ANDs of loop2 (window 4) against the 32-bit register."""
    return accumulator_bound(terms, value_bits).bit_length()


# ------------------------------------------------------------------ small lookup engine (ANF, measurement erase)
def _anf(table, w, nb):
    out = []
    for j in range(nb):
        f = [(table[a] >> j) & 1 for a in range(1 << w)]
        for i in range(w):
            for a in range(1 << w):
                if a >> i & 1:
                    f[a] ^= f[a ^ (1 << i)]
        out.append([S for S in range(1 << w) if f[S]])
    return out


def _lookup(C, addr, table, out):
    """out ^= table[addr] (out: fresh |0> wires). Monomials of degree >= 2: AND tree (2^w - w - 1 ANDs), kept alive until
    _unlookup. -> state."""
    w = len(addr); mon = {}
    for i, a in enumerate(addr):
        mon[1 << i] = a
    for S in sorted(range(1, 1 << w), key=lambda s: (bin(s).count('1'), s)):
        if bin(S).count('1') < 2:
            continue
        hi = S.bit_length() - 1
        mon[S] = C.AND(mon[S ^ (1 << hi)], addr[hi])
    anf = _anf(table, w, len(out))
    for j, o in enumerate(out):
        for S in anf[j]:
            if S == 0:
                C.x(o)
            else:
                C.cx(mon[S], o)
    return (addr, mon, anf, out)


def _unlookup(C, st):
    """erase out by X-measurements (the phase correction of output j is z on the monomials of its ANF), then erase the
    monomials."""
    addr, mon, anf, out = st
    for j, o in enumerate(out):
        C.MZ(o, [('z', mon[S]) for S in anf[j] if S != 0])
    C.release(out)
    for S in sorted(mon, key=lambda s: (-bin(s).count('1'), -s)):
        if bin(S).count('1') < 2:
            continue
        hi = S.bit_length() - 1
        C.unAND(mon[S ^ (1 << hi)], addr[hi], mon[S]); C.release([mon[S]])


# ------------------------------------------------------------------ parameters
def _check_params(q, clen, nbits, w):
    if not (isinstance(q, int) and isinstance(clen, int) and isinstance(nbits, int)):
        raise ValueError('q, clen, nbits must be integers')
    if q < 1 or clen < 1 or nbits <= clen:
        raise ValueError('need q >= 1, clen >= 1, nbits > clen')
    if q > (1 << clen):
        raise ValueError('need q <= 2^clen (the final step subtracts q once)')
    if w is not None and w < 1:
        raise ValueError('window w must be >= 1')


def _geometry(q, clen, nbits, w):
    hb = nbits - clen
    J = (hb + w - 1) // w
    vmax = (1 << clen) - 1 + J * (q - 1)
    e = 0
    while (q << e) <= vmax:
        e += 1
    return J, e


def _build(q, clen, nbits, w):
    """gate list; -> (C, S, garbage wires in layout order (top, g))."""
    J, e = _geometry(q, clen, nbits, w)
    hb = nbits - clen
    C = Circ(); S = C.new(nbits)
    L, H = S[:clen], S[clen:]
    acc = L + C.new(e)
    for j in range(J):
        win = H[j * w:(j + 1) * w]; ww = len(win)
        tab = [(a << (clen + w * j)) % q for a in range(1 << ww)]
        out = C.new(clen)
        st = _lookup(C, win, tab, out)
        C.add(acc, out)
        _unlookup(C, st)
    top = acc[clen:]
    g = C.new(1)
    acc2 = L + g
    tab = [(a << clen) % q for a in range(1 << e)]
    out = C.new(clen)
    st = _lookup(C, top, tab, out)
    C.add(acc2, out)
    _unlookup(C, st)
    k = (1 << clen) - q
    C.add(L, [g[0] if (k >> i) & 1 else ZERO for i in range(clen)])
    return C, S, list(top) + list(g)


def fold_accumulator_ands(q, clen, nbits, w):
    """AND count of fold_accumulator(q, clen, nbits, w) (gate list only)."""
    _check_params(q, clen, nbits, w)
    if w is None:
        raise ValueError('w must be given here')
    return _build(q, clen, nbits, w)[0].nand


class ResidueFold:
    """A generated compression circuit. Attributes: q, clen, nbits, w, e, J, circ (andcirc.Circ), n_and, layout
    (register name -> wire numbers in the emitted QASM: 'S' (input), 'r' (= 0..clen-1), 'H', 'top', 'g'), nref (wires that
    carry input / outputs; the rest are ancillas that end at 0), order (internal wire ids in QASM order)."""

    def __init__(s, q, clen, nbits, w=None):
        _check_params(q, clen, nbits, w)
        if w is None:
            best = None
            for ww in range(2, 7):
                a = _build(q, clen, nbits, ww)[0].nand
                if best is None or a < best[0]:
                    best = (a, ww)
            w = best[1]
        s.q, s.clen, s.nbits, s.w = q, clen, nbits, w
        s.J, s.e = _geometry(q, clen, nbits, w)
        C, S, garb = _build(q, clen, nbits, w)
        s.circ = C; s.n_and = C.nand
        s.order = list(S) + garb
        s.nref = len(s.order)
        k = nbits
        s.layout = dict(S=list(range(nbits)), r=list(range(clen)), H=list(range(clen, nbits)),
                        top=list(range(k, k + s.e)), g=[k + s.e])
        assert s.nref == k + s.e + 1

    # -- counts / text
    def counts(s):
        return s.circ.counts()

    def qasm(s):
        return s.circ.qasm(s.order)

    def text_counts(s):
        return text_counts(s.qasm())

    # -- reference
    def reference(s, S):
        """-> dict register -> value for the input S (0 <= S < 2^nbits): r, H, top, g (python-int formulas of the
        construction)."""
        if not 0 <= S < (1 << s.nbits):
            raise ValueError('S out of range')
        clen, q, w = s.clen, s.q, s.w
        L, H = S & ((1 << clen) - 1), S >> clen
        v = L
        for j in range(s.J):
            Hj = (H >> (w * j)) & ((1 << w) - 1)
            v += (Hj << (clen + w * j)) % q
        top, lo = v >> clen, v & ((1 << clen) - 1)
        u = lo + ((top << clen) % q)
        g = u >> clen
        r = (u & ((1 << clen) - 1)) + (((1 << clen) - q) if g else 0)
        return dict(r=r, H=H, top=top, g=g)

    def wires_of(s, reg):
        return [s.order[q] for q in s.layout[reg]]

    def check(s, values):
        """exact simulation of the gate list on the inputs `values` (python ints < 2^nbits). -> dict(ok, why, bad):
        r < 2^clen and r = S (mod q) (the END property), r / H / top / g equal reference(), all other wires 0."""
        values = [int(v) for v in values]
        X = bits_of(values, s.nbits)
        v, ok, why = simulate(s.circ, {w: X[j] for j, w in enumerate(s.order[:s.nbits])})
        if not ok:
            return dict(ok=False, why=why, bad=None)
        refs = [s.reference(S) for S in values]
        for reg in ('r', 'H', 'top', 'g'):
            got = ints_of(v, s.wires_of(reg))
            for i, (a, ref) in enumerate(zip(got, refs)):
                if a != ref[reg]:
                    return dict(ok=False, why='register %s, input %d: %d != %d' % (reg, values[i], a, ref[reg]),
                                bad=values[i])
        got = ints_of(v, s.wires_of('r'))
        for S, r in zip(values, got):
            if r >= (1 << s.clen) or (r - S) % s.q:
                return dict(ok=False, why='end property fails for %d: r = %d' % (S, r), bad=S)
        live = sorted(set(range(s.circ.N)) - set(s.order))
        if live and v[live].any():
            return dict(ok=False, why='ancilla not returned to 0', bad=None)
        return dict(ok=True, why='', bad=None)


def fold_accumulator(q, clen, nbits, w=None):
    """Compression (S, 0) -> (r, H, top, g): r < 2^clen, r = S (mod q); see the module docstring.
    w = None: best window among 2..6."""
    return ResidueFold(q, clen, nbits, w)


def edge_inputs(q, clen, nbits, k=64, seed=1):
    """edge and random inputs S < 2^nbits: 0, 1, multiples of q +- 1, 2^nbits - 1, powers of two, random."""
    rng = random.Random(seed); top = 1 << nbits
    out = {0, 1, top - 1, (1 << clen) - 1, 1 << clen}
    step = max(1, ((top - 1) // q) // 8)
    for m in range(0, (top - 1) // q + 1, step):
        for d in (-1, 0, 1):
            out.add(m * q + d)
    for j in range(nbits):
        out.add(1 << j); out.add((1 << j) - 1)
    while len(out) < k:
        out.add(rng.getrandbits(nbits))
    vals = sorted(x for x in out if 0 <= x < top)
    if len(vals) > k:
        keep = vals[:k // 2]; rest = vals[k // 2:]
        vals = sorted(set(keep) | set(rng.sample(rest, k - len(keep))))
    return vals
