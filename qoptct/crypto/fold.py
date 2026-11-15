# SPDX-License-Identifier: Apache-2.0
"""fold: reduction of a 2n-bit integer (a product) modulo a special prime on logical ANDs, and addition of a
classical constant (ported from the qoptct lab, agent k2; audited by agents k10 and x5).

    R = reduce_special(p, n=None, method='auto', step2='auto')     -> FoldReduction
    R.qasm()            OpenQASM 3: wires 0..2n-1 = P (input), then the garbage registers, then ancillas
    R.n_and, R.counts() AND count (= T/4), cost-model counts
    R.reference(P)      expected values of all output registers (independent python-int formulas)
    R.check(values)     exact bit-sliced simulation of the gate list on the given inputs (see andcirc.simulate)
    fold_and_count(p)   AND count only (builds the gate list, no simulation)
    add_const(n, K), controlled_add_const(n, K)  -> ConstAdd

FUNCTION. (P, 0...) -> (r, H, garbage...), r = P mod p FULLY reduced (0 <= r < p), on wires 0..n-1 (the low half of
P is overwritten in place); H = P >> n stays on wires n..2n-1. All other wires end at 0 except the documented
garbage registers (as in the Montgomery multiplier of Litinski 2023: garbage is removed by running the
multiplication backwards). Every P < 2^(2n) is allowed (not only products of reduced inputs).

METHODS.
'fold' (p = 2^n - c, condition c^2 + 2c + 2 < 2^n; p need not be prime):
    1) t = L + c*H: H<<s is added into L (extended by e fresh wires) for each signed digit (s, +-1) of c (NAF; the
       binary form when it is not heavier);
    2) t_hi += 1 (in place, one extra wire), u' = t_lo + c*(t_hi+1) with carry g ('M': M = c*(t_hi+1) is formed out
       of place -- disjoint positive digits are CX copies -- and added once; 'direct': each digit of c is added
       directly). As u = t_lo + c*t_hi < 2p: g = [u >= p], and for g = 1 the low n bits are already u - p;
    3) if g = 0: subtract c (subtraction of the constant controlled by NOT g).
    Garbage: H, thi1 = t_hi + 1 (e+1 bits), g (1 bit), and M = c*thi1 (wM bits, step2 'M' only).
'mersenne' = 'fold' with c = 1 (p = 2^n - 1, P-521): 3n - 1 ANDs.
'solinas' (P-256 only): step 1 is the NIST formula (FIPS 186-4 D.2.3) t = s1 + 2s2 + 2s3 + s4 + s5 - s6 - s7 - s8
    - s9 on 32-bit words; its exact range is [-4p, 5p], so 5p is added -- for free: the bits of 5p = A - B*2^192 are
    placed into constant-zero words of the addends (and B into a zero word of the subtrahend s8); then steps 2-3 with
    c = 2^224 - 2^192 - 2^96 + 1 added digit by digit ('direct'). Garbage: H, thi1 (5 bits), g.
'auto': P-256 -> 'solinas', c = 1 -> 'mersenne', otherwise 'fold'.

AND COUNTS (measured on the generated gate lists; tests pin them):
    p = 2^n - 19 (Curve25519 at n = 255): 5n + 29 for n = 100..521 -> 1304 at n = 255
    secp256k1, p = 2^256 - 2^32 - 977: 7n + 383 -> 2175
    P-256 (solinas): 2744;  P-384 (fold, 'M'): 3333;  P-521 (mersenne): 3n - 1 = 1562
    add_const(n, K): max(0, n - 2 - tz(K));  controlled_add_const(n, K): max(0, n - 1 - tz(K))   (tz = trailing zeros;
    K = 0 costs nothing). The x5 audit: at K = 2^(n-1) the formula n - 2 - tz gives -1, the circuit 0.
Cost model: AND = 4 T + 3 CX (counts()); the emitted AND gadget text has 6 cx lines per AND (text_counts()).

CONTEXT (read before quoting the numbers). There is no published stand-alone "reduction" block these counts can be
compared with: the 0.25n^2 + 9n often quoted as a baseline (18 551 at n = 255) is the lab's estimate of the
reduction share INSIDE the generic (any p) Montgomery multiplier of Litinski 2023 (arXiv:2306.08585, Fig. 8,
2.25n^2 + 9n in total). The fair statement is at the level of a full multiplication for special p only: schoolbook
2n^2 - n plus this fold reduction is about -11 % against 2.25n^2 + 9n. Arithmetic specialised to pseudo-Mersenne
moduli / secp256k1 already appears in 2026 work (Schrottenloher, arXiv:2606.02235 sec. 4; Haener et al.,
arXiv:2609.05625, with an approximate reduction for squaring mod secp256k1); multiplications on conditionally
inverted adders (2609.05625) are cheaper than schoolbook. Not a record.
"""
import random

import numpy as np

from .andcirc import Circ, ZERO, ONE, const_bits, simulate, bits_of, ints_of, text_counts, tz

P25519 = (1 << 255) - 19
SECP256K1 = (1 << 256) - (1 << 32) - 977
P256 = (1 << 256) - (1 << 224) + (1 << 192) + (1 << 96) - 1
P384 = (1 << 384) - (1 << 128) - (1 << 96) + (1 << 32) - 1
P521 = (1 << 521) - 1
CURVES = {'curve25519': P25519, 'secp256k1': SECP256K1, 'p256': P256, 'p384': P384, 'p521': P521}

METHODS = ('auto', 'fold', 'mersenne', 'solinas')


# ------------------------------------------------------------------ helpers
def signed_digits(c):
    """minimal-weight signed representation of c > 0 (NAF); the binary one when it is not heavier.
    -> [(shift, sign)]"""
    b = [(i, 1) for i in range(c.bit_length()) if (c >> i) & 1]
    naf = []; k = c; i = 0
    while k:
        if k & 1:
            d = 2 - (k & 3); naf.append((i, d)); k -= d
        k >>= 1; i += 1
    return b if len(b) <= len(naf) else naf


def _shifted(ws, s, W):
    v = [ZERO] * s + list(ws)
    return (v + [ZERO] * W)[:W]


def _add_digits(C, acc, src, digits):
    """acc += sum sign * (src << s); positive digits first (the intermediate values stay >= 0)."""
    for s, sg in sorted(digits, key=lambda d: -d[1]):
        if s >= len(acc):
            continue
        a = _shifted(src, s, len(acc))
        if sg > 0:
            C.add(acc, a)
        else:
            C.sub(acc, a)


def _mul_const_oop(C, src, c, digits, width):
    """M = c * src into a fresh register; disjoint positive digits into the empty register are CX copies."""
    M = C.new(width); occupied = set()
    pos = [d for d in digits if d[1] > 0]; neg = [d for d in digits if d[1] < 0]
    for s, _ in sorted(pos, key=lambda d: -d[0]):
        rng = set(range(s, s + len(src)))
        if not (rng & occupied) and s + len(src) <= width:
            for j, w in enumerate(src):
                C.cx(w, M[s + j])
            occupied |= rng
        else:
            C.add(M, _shifted(src, s, width)); occupied |= set(range(s, width))
    for s, _ in neg:
        C.sub(M, _shifted(src, s, width))
    return M


def _final(C, lo, thi, n, c, digits, step2, garbage):
    """steps 2-3: u' = lo + c*(t_hi+1) with carry g; if g = 0: lo -= c.  g = 1 <=> u >= p (u = lo + c*t_hi < 2p)."""
    ext = C.new(1)
    C.add(list(thi) + ext, [ONE])                 # t_hi += 1 in place
    thi = list(thi) + ext; garbage += ext
    e = len(thi)
    g = C.new()[0]; garbage.append(g)
    if step2 == 'M':
        wM = min((c * ((1 << e) - 1)).bit_length(), n)
        M = _mul_const_oop(C, thi, c, digits, wM); garbage += M
        C.add(lo, M, cout=g)
    else:
        _add_digits(C, lo + [g], thi, digits)
    C.x(g)
    C.sub(lo, [g if (c >> i) & 1 else ZERO for i in range(n)])
    C.x(g)


def _fold_e(n, c):
    return ((((1 << n) - 1) * (c + 1)) >> n).bit_length()


def _fold_wM(n, c):
    e = _fold_e(n, c)
    return min((c * ((1 << (e + 1)) - 1)).bit_length(), n)


def _reduce_pm(C, P, n, c, step2):
    L, H = P[:n], P[n:2 * n]
    dg = signed_digits(c)
    Sp = sum(1 << s for s, sg in dg if sg > 0)
    W1 = ((Sp + 1) * ((1 << n) - 1)).bit_length()
    acc = L + C.new(W1 - n)
    _add_digits(C, acc, H, dg)
    e = _fold_e(n, c)
    zero_top = acc[n + e:]                         # t < 2^(n+e): these wires are back at 0
    thi = acc[n:n + e]; lo = acc[:n]
    garbage = list(H) + list(thi)
    _final(C, lo, thi, n, c, dg, step2, garbage)
    C.release(zero_top)
    return lo, garbage


def p256_terms():
    """NIST formula (FIPS 186-4 D.2.3): words A0..A15 (32 bits) of each s_i, most significant word first."""
    S = {
        's1': (7, 6, 5, 4, 3, 2, 1, 0), 's2': (15, 14, 13, 12, 11, None, None, None),
        's3': (None, 15, 14, 13, 12, None, None, None), 's4': (15, 14, None, None, None, 10, 9, 8),
        's5': (8, 13, 15, 14, 13, 11, 10, 9), 's6': (10, 8, None, None, None, 13, 12, 11),
        's7': (11, 9, None, None, 15, 14, 13, 12), 's8': (12, None, 10, 9, 8, 15, 14, 13),
        's9': (13, None, 11, 10, 9, None, 15, 14)}
    coef = {'s1': 1, 's2': 2, 's3': 2, 's4': 1, 's5': 1, 's6': -1, 's7': -1, 's8': -1, 's9': -1}
    return S, coef


def p256_formula(P):
    """t = s1 + 2s2 + 2s3 + s4 + s5 - s6 - s7 - s8 - s9 (python int, t == P mod P-256)."""
    S, coef = p256_terms(); w = [(P >> (32 * i)) & 0xffffffff for i in range(16)]; t = 0
    for k, words in S.items():
        v = 0
        for pos, idx in enumerate(reversed(words)):
            if idx is not None:
                v |= w[idx] << (32 * pos)
        t += coef[k] * v
    return t


def _reduce_p256(C, P):
    n = 256; c = (1 << n) - P256
    S, coef = p256_terms()

    def term(name):
        ws = []
        for idx in reversed(S[name]):
            ws += ([ZERO] * 32) if idx is None else P[32 * idx:32 * idx + 32]
        return ws
    # t + K, K = 5p = A - B*2^192: the bits of A go into constant-zero positions of the addends (2s2: bits 0..96,
    # s4: 96..191, 2s3: 225..256, fresh wires 256..259: X), B into word 6 of the subtrahend s8 -> 0 extra ANDs.
    K = 5 * P256
    assert (K >> 224) & 1
    Kmid = (K % (1 << 225)) >> 192; Bw = (1 << 33) - Kmid
    assert 0 <= Bw < (1 << 32)
    A = K + (Bw << 192)
    assert (A >> 192) & ((1 << 33) - 1) == 0 and A < (1 << 260)
    acc = P[:256] + C.new(4)
    t2, t3 = _shifted(term('s2'), 1, 260), _shifted(term('s3'), 1, 260)
    t4, t8 = term('s4') + [ZERO] * 4, term('s8') + [ZERO] * 4
    for i in range(260):
        if not (A >> i) & 1:
            continue
        if i >= 256:
            C.x(acc[i]); continue
        tgt = t2 if i <= 96 else (t4 if i <= 191 else t3)
        assert tgt[i] == ZERO, i
        tgt[i] = ONE
    for j in range(32):
        if (Bw >> j) & 1:
            assert t8[192 + j] == ZERO
            t8[192 + j] = ONE
    C.add(acc, t2); C.add(acc, t3); C.add(acc, t4); C.add(acc, term('s5') + [ZERO] * 4)
    for nm in ('s6', 's7'):
        C.sub(acc, term(nm) + [ZERO] * 4)
    C.sub(acc, t8); C.sub(acc, term('s9') + [ZERO] * 4)
    lo, thi = acc[:256], acc[256:260]
    garbage = list(P[256:]) + list(thi)
    _final(C, lo, thi, n, c, signed_digits(c), 'direct', garbage)
    return lo, garbage


# ------------------------------------------------------------------ public API
def _resolve(p, n, method, step2):
    if not isinstance(p, int) or p < 3:
        raise ValueError('p must be an integer >= 3')
    if method not in METHODS:
        raise ValueError('method must be one of %s' % (METHODS,))
    if n is None:
        n = 256 if p == P256 else p.bit_length()
    c = (1 << n) - p
    if c <= 0:
        raise ValueError('p must be < 2^n')
    if method == 'auto':
        method = 'solinas' if (p == P256 and n == 256) else ('mersenne' if c == 1 else 'fold')
    if method == 'solinas':
        if p != P256 or n != 256:
            raise ValueError("method 'solinas' is implemented for P-256 only")
        if step2 not in ('auto', 'direct'):
            raise ValueError("method 'solinas' uses step2 'direct'")
        return n, c, method, 'direct'
    if method == 'mersenne' and c != 1:
        raise ValueError("method 'mersenne' needs p = 2^n - 1")
    if not (c * c + 2 * c + 2 < (1 << n)):
        raise ValueError('fold: c = 2^n - p too large for one fold (need c^2 + 2c + 2 < 2^n)')
    dg = signed_digits(c)
    if step2 == 'auto':
        step2 = 'M' if len(dg) >= 2 else 'direct'
    if step2 not in ('M', 'direct'):
        raise ValueError("step2 must be 'auto', 'M' or 'direct'")
    return n, c, method, step2


class FoldReduction:
    """A generated reduction circuit. Attributes: p, n, c, method, step2, circ (andcirc.Circ), n_and, layout
    (register name -> list of wire numbers in the emitted QASM: 'P' (input, 2n), 'r' (= wires 0..n-1),
    'H', 'thi1', 'g', 'M' (fold 'M' only)), nref (number of wires carrying input / outputs; the rest are ancillas that
    end at 0), order (internal wire ids in QASM order)."""

    def __init__(s, p, n=None, method='auto', step2='auto'):
        s.p = p
        s.n, s.c, s.method, s.step2 = _resolve(p, n, method, step2)
        n = s.n
        C = Circ(); P = C.new(2 * n)
        if s.method == 'solinas':
            r, garb = _reduce_p256(C, P)
        else:
            r, garb = _reduce_pm(C, P, n, s.c, s.step2)
        assert r == P[:n]
        Pset = set(P)
        s.circ = C; s.n_and = C.nand
        s.order = list(P) + [w for w in garb if w not in Pset]
        s.nref = len(s.order)
        if s.method == 'solinas':
            e1 = 5; wM = 0
        else:
            e1 = _fold_e(n, s.c) + 1; wM = _fold_wM(n, s.c) if s.step2 == 'M' else 0
        k = 2 * n
        s.layout = dict(P=list(range(2 * n)), r=list(range(n)), H=list(range(n, 2 * n)),
                        thi1=list(range(k, k + e1)), g=[k + e1])
        if wM:
            s.layout['M'] = list(range(k + e1 + 1, k + e1 + 1 + wM))
        assert s.nref == k + e1 + 1 + wM, (s.nref, k, e1, wM)

    # -- counts / text
    def counts(s):
        return s.circ.counts()

    def qasm(s):
        return s.circ.qasm(s.order)

    def text_counts(s):
        return text_counts(s.qasm())

    # -- reference (independent of the gate list)
    def reference(s, P):
        """-> dict register -> value for input P (0 <= P < 2^(2n)): r = P mod p, H = P >> n, and the garbage."""
        n, c, p = s.n, s.c, s.p
        if not 0 <= P < (1 << (2 * n)):
            raise ValueError('P out of range')
        L, H = P & ((1 << n) - 1), P >> n
        t = 5 * P256 + p256_formula(P) if s.method == 'solinas' else L + c * H
        thi1 = (t >> n) + 1; lo = t & ((1 << n) - 1)
        g = (lo + c * thi1) >> n
        out = dict(r=P % p, H=H, thi1=thi1, g=g)
        if 'M' in s.layout:
            out['M'] = c * thi1
        return out

    def wires_of(s, reg):
        """internal wire ids of a register (for simulate)."""
        return [s.order[q] for q in s.layout[reg]]

    def check(s, values):
        """exact simulation of the gate list on the inputs `values` (python ints < 2^(2n)), see andcirc.simulate.
        -> dict(ok, why, bad): all registers equal reference(), every other wire 0."""
        if s.method != 'solinas' and 2 * s.n <= 40:
            return s._check_np(np.asarray(list(values), np.int64))
        values = [int(v) for v in values]
        X = bits_of(values, 2 * s.n)
        v, ok, why = simulate(s.circ, {w: X[j] for j, w in enumerate(s.order[:2 * s.n])})
        if not ok:
            return dict(ok=False, why=why, bad=None)
        refs = [s.reference(P) for P in values]
        for reg in s.layout:
            if reg == 'P':
                continue
            got = ints_of(v, s.wires_of(reg))
            want = [ref['r'] if reg == 'r' else ref[reg] for ref in refs]
            for i, (a, b) in enumerate(zip(got, want)):
                if a != b:
                    return dict(ok=False, why='register %s, input %d: %d != %d' % (reg, values[i], a, b), bad=values[i])
        live = sorted(set(range(s.circ.N)) - set(s.order))
        if live and v[live].any():
            return dict(ok=False, why='ancilla not returned to 0', bad=None)
        return dict(ok=True, why='', bad=None)

    def _check_np(s, P):
        """check() for small n (2n <= 40): vectorised int64 reference, same criteria."""
        n, c, p = s.n, s.c, s.p
        if P.size == 0 or P.min() < 0 or P.max() >= (1 << (2 * n)):
            raise ValueError('P out of range')
        X = {w: ((P >> j) & 1).astype(np.uint8) for j, w in enumerate(s.order[:2 * n])}
        v, ok, why = simulate(s.circ, X)
        if not ok:
            return dict(ok=False, why=why, bad=None)
        mask = (1 << n) - 1
        L, H = P & mask, P >> n
        t = L + c * H; thi1 = (t >> n) + 1
        want = dict(r=P % p, H=H, thi1=thi1, g=((t & mask) + c * thi1) >> n)
        if 'M' in s.layout:
            want['M'] = c * thi1
        for reg, val in want.items():
            ws = s.wires_of(reg)
            if (val >> len(ws)).any():
                return dict(ok=False, why='reference %s does not fit %d bits' % (reg, len(ws)), bad=None)
            got = np.zeros(P.shape, np.int64)
            for j, w in enumerate(ws):
                got |= v[w].astype(np.int64) << j
            badi = np.nonzero(got != val)[0]
            if badi.size:
                i = int(badi[0])
                return dict(ok=False, why='register %s, input %d: %d != %d' % (reg, P[i], got[i], val[i]), bad=int(P[i]))
        live = sorted(set(range(s.circ.N)) - set(s.order))
        if live and v[live].any():
            return dict(ok=False, why='ancilla not returned to 0', bad=None)
        return dict(ok=True, why='', bad=None)


def reduce_special(p, n=None, method='auto', step2='auto'):
    """Reduction (P, 0) -> (P mod p, H, garbage) of a 2n-bit P modulo a special p; see the module docstring."""
    return FoldReduction(p, n, method, step2)


def fold_and_count(p, n=None, method='auto', step2='auto'):
    """AND count of reduce_special(p, ...) (gate list only, no simulation)."""
    return FoldReduction(p, n, method, step2).n_and


def edge_inputs(p, n, k=64, seed=1):
    """edge and random inputs P < 2^(2n): 0, 1, multiples of p +- small, (p-1)^2, 2^(2n)-1, powers of two, random."""
    rng = random.Random(seed); top = 1 << (2 * n); out = {0, 1, top - 1, (p - 1) ** 2, (p - 1) * (p - 2)}
    for m in (1, 2, 3, (1 << n) - 1, p - 1, p, rng.getrandbits(n)):
        for d in (-2, -1, 0, 1, 2):
            out.add(m * p + d)
    for j in range(2 * n):
        out.add(1 << j); out.add((1 << j) - 1)
    out |= {(1 << n) * h + l for h in (0, 1, (1 << n) - 1, p - 1) for l in (0, 1, p - 1, p, (1 << n) - 1)}
    while len(out) < k:
        out.add(rng.getrandbits(2 * n))
    vals = sorted(x for x in out if 0 <= x < top)
    if len(vals) > k:
        keep = vals[:k // 2]; rest = vals[k // 2:]
        vals = sorted(set(keep) | set(rng.sample(rest, k - len(keep))))
    return vals


# ------------------------------------------------------------------ constant addition
class ConstAdd:
    """x := x + K (mod 2^n) on wires 0..n-1; controlled version: control on wire n (x unchanged when it is 0).
    No carry-out, ancillas end at 0. Attributes: n, K, controlled, circ, n_and, nref."""

    def __init__(s, n, K, controlled=False):
        if n < 1 or not 0 <= K < (1 << n):
            raise ValueError('need n >= 1 and 0 <= K < 2^n')
        s.n, s.K, s.controlled = n, K, controlled
        C = Circ(); x = C.new(n)
        if controlled:
            ct = C.new()[0]
            C.add(x, [ct if (K >> i) & 1 else ZERO for i in range(n)])
        else:
            C.add(x, const_bits(K, n))
        s.circ = C; s.n_and = C.nand; s.nref = n + (1 if controlled else 0)

    def expected_and(s):
        if s.K == 0:
            return 0
        return max(0, s.n - (1 if s.controlled else 2) - tz(s.K))

    def counts(s):
        return s.circ.counts()

    def qasm(s):
        return s.circ.qasm(list(range(s.nref)))

    def reference(s, x, ctrl=1):
        return (x + s.K) % (1 << s.n) if (ctrl or not s.controlled) else x

    def check(s, xs, ctrls=None):
        """exact simulation on inputs xs (and control bits ctrls for the controlled version)."""
        xs = [int(v) for v in xs]
        if ctrls is None:
            ctrls = [1] * len(xs)
        X = bits_of(xs, s.n)
        inp = {j: X[j] for j in range(s.n)}
        if s.controlled:
            inp[s.n] = np.array(ctrls, np.uint8)
        v, ok, why = simulate(s.circ, inp)
        if not ok:
            return dict(ok=False, why=why)
        got = ints_of(v, range(s.n))
        for x, cb, gv in zip(xs, ctrls, got):
            if gv != s.reference(x, cb):
                return dict(ok=False, why='x=%d ctrl=%d: %d' % (x, cb, gv))
        if s.controlled and not (v[s.n] == np.array(ctrls, np.uint8)).all():
            return dict(ok=False, why='control changed')
        if s.circ.N > s.nref and v[s.nref:].any():
            return dict(ok=False, why='ancilla not returned to 0')
        return dict(ok=True, why='')


def add_const(n, K):
    """x += K (mod 2^n): max(0, n - 2 - tz(K)) ANDs (0 for K = 0)."""
    return ConstAdd(n, K, controlled=False)


def controlled_add_const(n, K):
    """if ctrl: x += K (mod 2^n); control on wire n: max(0, n - 1 - tz(K)) ANDs (0 for K = 0)."""
    return ConstAdd(n, K, controlled=True)
