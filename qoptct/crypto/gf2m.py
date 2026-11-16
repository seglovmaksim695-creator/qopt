# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""qoptct.crypto.gf2m -- GF(2^m) arithmetic circuits (qoptct v0.5): multiplication, squaring, inversion, division.

Model: Clifford+T with logical ANDs (Gidney 2018: 4 T, target |0>) and measurement-based erasure (X measurement +
classically controlled CZ / Z, 0 T).  "AND" below = Toffoli-equivalent of the literature tables.

    gf2m_multiply(m, f=None, mode='strict'|'clean')   c ^= a*b mod f  /  c = a*b mod f (c starts at 0)
    gf2m_square(m, f=None, relabel=False)              a <- a^2 mod f in place, CNOT only
    gf2m_inverse(m, f=None, chain='it', p=None)        t ^= x^-1 (0 -> 0), x kept
    gf2m_divide(m, f=None, chain='it', p=None)         t ^= y/x (x = 0 -> t ^= 0), x, y kept
    gf2m_mult_count / gf2m_inverse_count / gf2m_divide_count   the same costs without building gates
    gf2m_bilinear(m, f=None)                           the bilinear algorithm (alpha, beta, gamma) used
    verify_multiplier(C, exact=False)                  bitsliced check / exact degree-2 check of a multiplier
    default_modulus(m)                                 NIST modulus for 163/233/283/409/571, else low weight

f is the modulus as an int (bit i = coefficient of x^i); default: default_modulus(m).  The circuit functions are
pure: same arguments -> same circuit; they return qoptct.contrib.gf2m.qcirc.GF2Circuit (ops, .qasm(), .counts(),
.simulate(), .regs = named wire lists, .meta).

Numbers (NIST moduli; reproduced by tests/test_crypto_gf2m.py by counting, without simulation):
    multiplication ANDs  m=163: 908, 233: 1343, 283: 1669, 409: 2501, 571: 3569.
    Against Kim et al. 2022 (ePrint 2022/1095, CRT: 992 / 1441 / 1784 / - / 3813 Toffoli) and Garn-Kan 2025
    (999 / 1448 / 1776 / - / 3860) the cost model is the same (number of ANDs = bilinear rank, products uncomputed
    for free by measurement).  NOT a record of bilinear complexity: smaller symmetric bilinear complexities are
    published (Chudnovsky-type algorithms, arXiv:1303.7082, Table 3: 906 / 1340 / 1668 / 2495 / 3566).
    division ANDs = (L+1) * M with p = L chain registers (L = chain length), inversion L * M; m=163 (IT chain,
    L = 9): division 9080 (Taguchi-Takagi ePrint 2023/553 Table 4: 18848 = 19 * 992), inversion 8172.  The free
    erasure of a GF(2^m) product register by X measurement + CZ correction is published (Khattar, Shutty, Gidney
    et al., arXiv:2510.10967, App. C.3, Qualtran GF2MulMBUC); new in the qoptct lab is its use along the
    Itoh-Tsujii chain (pebbling with free removal).
    division ANDs (exact count / (L+1)*M): 163: 9080, 233: 14773, 283: 20028 (IT chain, L = 9/10/11);
    409: 27503 / 27511, 571: 46397 (shortest star chain, L = 10/12; IT chain: 30004 / 30012, 49966) -- a few
    terms of the doubling steps x^(2^k) * x are linear and need no AND.
    squaring CX (exact wire order, seeds=3): 798 / 825 / 1386 / 1418 / 3239; with free wire relabelling
    408 / 153 / 713 / 204 / 1538 (the earlier lab claims 404 / 709 were not reproduced; the published bound
    used for comparison in the lab, n^2 - n from LUP decomposition, is an upper bound, not a built circuit).
    multiplier CX (cx + cz of the text, seed 0): 163 strict 75823, clean 68371.
"""
import random

from qoptct.contrib.gf2m import crt as _crt
from qoptct.contrib.gf2m import mult as _mult
from qoptct.contrib.gf2m import square as _square
from qoptct.contrib.gf2m import itoh as _itoh
from qoptct.contrib.gf2m.field import is_irreducible
from qoptct.contrib.gf2m.qcirc import GF2Circuit, PreconditionError, degree2_exact

__all__ = ["default_modulus", "gf2m_bilinear", "gf2m_multiply", "gf2m_mult_count", "gf2m_square",
           "gf2m_inverse", "gf2m_divide", "gf2m_inverse_count", "gf2m_divide_count", "gf2m_itoh_gate_check",
           "verify_multiplier", "GF2Circuit", "PreconditionError"]

default_modulus = _crt.default_modulus
EXACT_CHECK_MAX_M = 300          # exact tensor check of the bilinear algorithm up to this m (cost ~ rank * m^2)


def _modulus(m, f):
    if m < 2:
        raise ValueError("m must be >= 2")
    f = _crt.default_modulus(m) if f is None else int(f)
    if f.bit_length() - 1 != m or not is_irreducible(f):
        raise ValueError("f must be an irreducible polynomial of degree m")
    return f


def gf2m_bilinear(m, f=None, algorithm="auto"):
    """-> (dec, info).  dec: list of (alpha, beta, gamma) masks with c = sum (alpha.a)(beta.b) gamma = a*b mod f.
    algorithm: 'crt' (CRT plan; optimal over the stored places for m <= 600), 'table' (stored flip-graph / tower
    decompositions of data/decomps.json, m <= 32), 'auto' (the shorter; CRT on ties).
    info: algorithm, rank, places (crt) or source (table), check (exact up to m = EXACT_CHECK_MAX_M)."""
    f = _modulus(m, f)
    key = (m, f, algorithm)
    if key not in _BIL_CACHE:
        if len(_BIL_CACHE) > 16:
            _BIL_CACHE.clear()
        _BIL_CACHE[key] = _bilinear(m, f, algorithm)
    dec, info = _BIL_CACHE[key]
    return list(dec), dict(info)


_BIL_CACHE = {}


def _bilinear(m, f, algorithm):
    cands = []
    if algorithm in ("auto", "crt"):
        dec, places = _crt.field_algorithm(m, f)
        cands.append((len(dec), 0, dec, dict(algorithm="crt", rank=len(dec), places=places)))
    if algorithm in ("auto", "table"):
        for k, (dec, src) in enumerate(_mult.table_decs(m, f)):
            cands.append((len(dec), 1 + k, dec, dict(algorithm="table", rank=len(dec), source=src)))
    if not cands:
        raise ValueError("no bilinear algorithm (algorithm=%r, m=%d)" % (algorithm, m))
    cands.sort(key=lambda c: (c[0], c[1]))
    _, _, dec, info = cands[0]
    exact = m <= EXACT_CHECK_MAX_M
    ok = _crt.check_algorithm(dec, m, f) if exact else _crt.check_random(dec, m, f)
    if not ok:
        raise AssertionError("bilinear algorithm failed its check")
    info["check"] = "exact (all coefficient pairs)" if exact else "construction + 67 random/edge pairs"
    return dec, info


def gf2m_mult_count(m, f=None, algorithm="auto", n_anc=1):
    """cost of gf2m_multiply without building it: AND = T/4 = rank, qubits 3m + n_anc."""
    dec, info = gf2m_bilinear(m, f, algorithm)
    return dict(AND=len(dec), T=4 * len(dec), qubits=3 * m + n_anc, algorithm=info["algorithm"])


def gf2m_multiply(m, f=None, mode="strict", algorithm="auto", seed=0, n_anc=1):
    """GF(2^m) multiplier: wires a = 0..m-1, b = m..2m-1, c = 2m..3m-1, ancillas after (n_anc for 'crt', 1 for
    'table').  mode 'strict': c ^= a*b mod f for any c;  'clean': c = a*b mod f, c must start at 0.  a, b are
    unchanged, ancillas start and end at 0.  #AND = rank of the bilinear algorithm (gf2m_bilinear)."""
    if mode not in ("strict", "clean"):
        raise ValueError("mode must be 'strict' or 'clean'")
    f = _modulus(m, f)
    clean = mode == "clean"
    dec, info = gf2m_bilinear(m, f, algorithm)
    if info["algorithm"] == "crt":
        C, _ = _mult.build_crt(m, f, clean=clean, seed=seed, n_anc=n_anc, places=info["places"])
    else:
        C = _mult.build_table(m, f, dec, clean=clean, seeds=(seed, seed + 1, seed + 2))
    C.meta.update(kind="gf2m_multiply", m=m, f=f, mode=mode, algorithm=info["algorithm"], rank=info["rank"])
    if C.counts()["AND"] > info["rank"]:
        raise AssertionError("more ANDs than the rank")
    return C


def verify_multiplier(C, lanes=64, seed=1, exact=False):
    """check a gf2m_multiply circuit.  exact=True: degree-2 argument (all inputs, all measurement outcomes; the
    lane integers have 1 + n + n(n-1)/2 bits, n = 3m strict / 2m clean); otherwise `lanes` random inputs plus edge
    cases (bitsliced, every AND / erasure precondition enforced).  -> (ok, reason)."""
    m, f, clean = C.meta["m"], C.meta["f"], C.meta["mode"] == "clean"
    a_w, b_w, c_w = C.regs["a"], C.regs["b"], C.regs["c"]
    if exact:
        red = [[_crt.pmod(1 << (i + j), f) for j in range(m)] for i in range(m)]
        ins = a_w + b_w + ([] if clean else c_w)

        def ref(v):
            out = {c_w[t]: (0 if clean else v[c_w[t]]) for t in range(m)}
            for i in range(m):
                for j in range(m):
                    r = red[i][j]
                    if r:
                        pr = v[a_w[i]] & v[b_w[j]]
                        while r:
                            t = (r & -r).bit_length() - 1
                            r &= r - 1
                            out[c_w[t]] ^= pr
            return out
        return degree2_exact(C, ins, ref)
    rng = random.Random(seed)
    A = [0, (1 << m) - 1, 1, 1 << (m - 1)] + [rng.getrandbits(m) for _ in range(lanes)]
    B = [0, (1 << m) - 1, (1 << m) - 1, 1 << (m - 1)] + [rng.getrandbits(m) for _ in range(lanes)]
    Cv = [0] * len(A) if clean else [rng.getrandbits(m) for _ in range(len(A))]
    try:
        out = C.simulate(dict(a=A, b=B, c=Cv), len(A))
    except PreconditionError as e:
        return False, str(e)
    if out["a"] != A or out["b"] != B:
        return False, "inputs not restored"
    for k in range(len(A)):
        if out["c"][k] != Cv[k] ^ _crt.mulmod(A[k], B[k], f):
            return False, "wrong product on input %d" % k
    return True, "%d random + 4 edge inputs (all preconditions held)" % lanes


def gf2m_square(m, f=None, relabel=False, seeds=3):
    """in-place squaring a <- a^2 mod f as a CNOT network.  relabel=False: bit i of a^2 on wire i;
    relabel=True: bit i of a^2 on wire C.regs['out'][i] (fewer CX).  Checked exactly (linear map: all unit
    vectors)."""
    f = _modulus(m, f)
    C = _square.build(m, f, relabel=relabel, seeds=seeds)
    C.meta.update(kind="gf2m_square", m=m, f=f, relabel=relabel)
    return C


# ------------------------------------------------------------------------------------------- Itoh-Tsujii

def _samples(m, n, seed):
    rng = random.Random(seed)
    xs = [0, 1, (1 << m) - 1] + [rng.getrandbits(m) for _ in range(n)]
    ys = [rng.getrandbits(m) for _ in range(len(xs))]
    return xs, ys


_CRT_CACHE = {}


def _crt_data(m, f, places):
    key = (m, f, tuple(places))
    if key not in _CRT_CACHE:
        if len(_CRT_CACHE) > 8:
            _CRT_CACHE.clear()
        _CRT_CACHE[key] = _mult.CRTData(m, f, places)
    return _CRT_CACHE[key]


def _itoh_run(kind, m, f, chain, p, use_p, gates, samples, seed, algorithm, count_cz, erase="terms",
              lowering="auto"):
    f = _modulus(m, f)
    dec, info = gf2m_bilinear(m, f, algorithm)
    if lowering not in ("auto", "crt", "parity"):
        raise ValueError("lowering must be 'auto', 'crt' or 'parity'")
    if lowering == "crt" and info["algorithm"] != "crt":
        raise ValueError("lowering='crt' needs the CRT algorithm (algorithm='crt')")
    use_crt = gates and erase == "terms" and info["algorithm"] == "crt" and lowering in ("auto", "crt")
    ch = _itoh.Chain(m, chain, p, use_p)
    xs, ys = _samples(m, samples, seed)
    regs = dict(x=xs, y=ys) if kind == "divide" else dict(x=xs)
    P = _itoh.Program(m, f, dec, regs, gates=gates, count_cz=count_cz, erase_mode=erase,
                      crt=_crt_data(m, f, info["places"]) if use_crt else None)
    P.alloc("t")
    P._fresh.discard("t")          # t is an input register (t ^= ...), not known to be 0
    if kind == "divide":
        _itoh.divide_program(P, ch, "t", "x", "y")
    else:
        _itoh.inverse_program(P, ch, "t", "x")
    if set(P.live_regs) != set(regs) | {"t"}:
        raise AssertionError("chain registers left alive")
    Fd = _itoh.Field(m, f)
    want = [Fd.mul(y, Fd.inv(x)) for x, y in zip(xs, ys)] if kind == "divide" else [Fd.inv(x) for x in xs]
    if P.val["t"] != want:
        raise AssertionError("register program gives a wrong result")
    P.finish()
    cost = dict(mult=P.n_mul, AND=P.and_count, T=4 * P.and_count, rank=len(dec), chain_length=ch.L, p=ch.p,
                register_erasures=P.n_erase, peak_registers=P.peak_regs, qubits_registers=P.peak_regs * m,
                peak_wires=P.peak_wires, formula_AND=P.n_mul * len(dec), erase=erase)
    if count_cz or gates or erase == "terms":
        cost["CZ_register_erasures"] = P.cz_count
    if gates == "stream":
        S = P.C
        for name, vals in list(regs.items()) + [("t", want)]:
            if S.read(P.wires[name]) != vals:
                raise AssertionError("gate-level stream: register %s wrong" % name)
        regw = set(w for name in list(regs) + ["t"] for w in P.wires[name])
        if any(v for w, v in S.st.items() if w not in regw):
            raise AssertionError("gate-level stream: ancilla not clean")
        cost["gates"] = S.counts()
        cost["gate_check"] = ("%d inputs at gate level (every AND / erasure precondition, register erasures "
                              "checked lazily)" % len(xs))
        return cost
    if not gates:
        return cost
    C = P.C
    C.regs = {name: list(P.wires[name]) for name in ("x", "y", "t") if name in P.wires}
    C.meta.update(kind="gf2m_" + kind, m=m, f=f, chain=ch.ch, p=ch.p, algorithm=info["algorithm"],
                  rank=len(dec), cost=cost)
    return C


_ITOH_DOC = """
    chain: 'it' (Itoh-Tsujii), 'short' (shortest star chain) or an explicit chain [(e, i, j), ...];
    p: chain registers (default L = chain length -> minimal multiplications; fewer -> recomputation);
    erase: 'terms' (default; register erasure = X measurements + one classically controlled CZ per (bilinear term,
    output bit): O(rank * m) CZ and the parity CX of the terms) or 'monomials' (per-bit quadratic forms, O(m^3) CZ;
    the k6 lab form).
    lowering: 'auto' (= 'crt' when the CRT algorithm is used and erase='terms', else 'parity'), 'crt' (every
    product / erasure with the in-place CRT compiler of mult.py: ~5x fewer CX; a product x^(2^k) * x of a
    register with itself works on a temporary copy, +m qubits) or 'parity' (term by term on parity wires: fewer
    qubits, many more CX).  The register program is checked on `samples` random field elements (every erasure
    condition) before lowering."""


def gf2m_inverse(m, f=None, chain="it", p=None, use_p=None, samples=16, seed=1, algorithm="auto", erase="terms",
                 lowering="auto"):
    """out-of-place inversion circuit t ^= x^-1 (x = 0 -> 0); wires x = 0..m-1, t = m..2m-1, then the chain
    registers and temporaries; L multiplications with p = L chain registers."""
    return _itoh_run("inverse", m, f, chain, p, use_p, True, samples, seed, algorithm, False, erase, lowering)


def gf2m_divide(m, f=None, chain="it", p=None, use_p=None, samples=16, seed=1, algorithm="auto", erase="terms",
                lowering="auto"):
    """division circuit t ^= y / x (x = 0 -> t ^= 0); wires x = 0..m-1, y = m..2m-1, t = 2m..3m-1, then chain
    registers / temporaries; (L + 1) multiplications with p = L chain registers."""
    return _itoh_run("divide", m, f, chain, p, use_p, True, samples, seed, algorithm, False, erase, lowering)


gf2m_inverse.__doc__ += _ITOH_DOC
gf2m_divide.__doc__ += _ITOH_DOC


def gf2m_inverse_count(m, f=None, chain="it", p=None, use_p=None, samples=8, seed=1, algorithm="auto",
                       erase="terms", count_cz=False):
    """costs of gf2m_inverse without gates: exact AND count (= nonzero bilinear terms of every product), CZ of the
    register erasures (erase='terms'; for 'monomials' only with count_cz=True: O(m^2) field products per erasure
    type); the register program is checked on samples."""
    return _itoh_run("inverse", m, f, chain, p, use_p, False, samples, seed, algorithm, count_cz, erase)


def gf2m_divide_count(m, f=None, chain="it", p=None, use_p=None, samples=8, seed=1, algorithm="auto",
                      erase="terms", count_cz=False):
    """costs of gf2m_divide without gates (see gf2m_inverse_count)."""
    return _itoh_run("divide", m, f, chain, p, use_p, False, samples, seed, algorithm, count_cz, erase)


def gf2m_itoh_gate_check(kind, m, f=None, chain="it", p=None, use_p=None, lanes=32, seed=1, algorithm="auto",
                         lowering="auto"):
    """gate-level check for large m without storing the circuit: kind 'divide' / 'inverse'.  Every gate of the
    products (and the parity CX of the 'terms' register erasures) is executed bitsliced on `lanes` random inputs
    (plus 0, 1, all-ones) with every AND / erasure precondition enforced; register erasures are checked lazily.
    -> cost dict with 'gates' (gate counts: CX = cx + 3 per AND + 1 CZ per AND erasure; the conditional CZ of
    the register erasures are 'CZ_register_erasures') and 'gate_check'; raises on any failure."""
    if kind not in ("divide", "inverse"):
        raise ValueError("kind must be 'divide' or 'inverse'")
    return _itoh_run(kind, m, f, chain, p, use_p, "stream", lanes, seed, algorithm, False, "terms", lowering)
