# SPDX-License-Identifier: Apache-2.0
"""qoptct.estimate.logical: end-to-end LOGICAL resource models (Toffoli / logical-AND count and logical qubits) of
published quantum algorithms, and the same algorithms with the building blocks of qoptct.crypto.

This is a resource-estimation tool of the kind used in post-quantum-migration planning ("how many logical operations
does the best published algorithm need"). It builds no attack circuits: every function is arithmetic on published
formulas and on the gate counts of the blocks.

    from qoptct.estimate import logical
    logical.gidney2025_rsa(2048).toffoli          # 6.49e9, 1399 qubits   (paper, Table 5 / its code)
    logical.rsa_g25_qoptct(2048).toffoli          # 4.91e9, 1400 qubits   (blocks of qoptct; model, no end-to-end circuit)
    logical.attack_report("ecdlp-p256")           # both rows + the physical estimate (qoptct.estimate.physical)

Unit.  One Toffoli = one logical AND of Gidney (arXiv:1709.06648: 4 T + 3 CX), erased by measurement at 0 T (the
"AND model"); T = 4 x AND.  A source that counts in another model says so in `model` and is compared in ITS model.

Two kinds of functions
----------------------
* published models (``ge2019_rsa``, ``gidney2025_rsa``, ``litinski2023_ecdlp``, ``haner2020_ecdlp``,
  ``roetteler2017_ecdlp``, ``banegas2020_ecdlp``, ``taguchi_takayasu2023_ecdlp``, ``jaques2020_aes``): re-implement the
  accounting of the paper and reproduce its tables within the tolerances given by ``reproduce_published()``.
* qoptct variants (``rsa_g25_qoptct``, ``ecdlp_prime_qoptct``, ``ecdlp_gf2m_qoptct``, ``aes_oracle_qoptct``): the same
  framework with some blocks replaced by checked blocks.  ONLY variants that survived the lab audit are offered.

What is NOT included (withdrawn in the lab review of 16-20 Nov 2026): the Gidney-Ekera-framework variant with
1.979e9 Toffoli (k1 "chain", reaction depth x2.1: physically worse, +91 % space-time volume) and the point
2.43e9 @ 5443 qubits (not reproduced by the independent check); windowed binary-curve numbers are offered only as
``status='estimate'``.

Status of a result (``LogicalEstimate.status``)
  'published'  reproduces a number of the paper itself
  'qoptct'     paper framework with blocks of qoptct that are proved or checked at the build level; the END-TO-END
               circuit is not assembled (it is a model), see ``note``
  'estimate'   not assembled and not checked in the parts that matter (flagged in ``note``)

Block counts come from qoptct.crypto where it has a counter (``library_blocks()`` recomputes them live and the tests
compare them with ``BLOCKS``); the entries that are lab measurements outside the library (the full point addition k4,
the RSA lookup/loop2 constants, the AES oracle) are listed with their provenance.
"""
import math
from collections import Counter
from dataclasses import dataclass, field, asdict

LG = math.log2


# ================================================================================================ result type
@dataclass
class LogicalEstimate:
    attack: str
    source: str            # 'published' | 'qoptct'
    toffoli: float         # AND/Toffoli count (model: see `model`)
    qubits: int
    model: str             # 'AND' (Gidney AND, erased by measurement) | 'Toffoli' | 'T/4'
    status: str            # 'published' | 'qoptct' | 'estimate'
    ref: str               # where the published number / framework is
    note: str = ''
    parts: dict = field(default_factory=dict)
    extra: dict = field(default_factory=dict)

    def to_dict(self):
        return asdict(self)


def _close(x, y, rel):
    return abs(x - y) <= rel * abs(y)


# ================================================================================================ blocks
# Frozen block counts.  'lib' = qoptct.crypto counter (library_blocks() recomputes it, tests compare);
# 'lab' = lab measurement outside the library (provenance in the comment).
BLOCKS = {
    # Kaliski inversion forward, reduction='special', n=256 (qoptct 0.5.x: 12n^2 - 8n - 3 + special reduction).
    'kaliski_special': {'p256': 664966, 'secp256k1': 657503},                 # lib
    # the 0.5.0rc circuit (12n^2 + n), byte-identical to the lab k4 circuits; point-add numbers below were built with it
    'kaliski_special_legacy': {'p256': 667273, 'secp256k1': 659810},          # lib (legacy=True)
    # GF(2^m) multiplier c ^= a*b, ANDs = bilinear rank (CRT algorithm), NIST moduli; lib gf2m_mult_count
    'gf2m_mult': {163: 908, 233: 1343, 283: 1669, 409: 2501, 571: 3569},
    # GF(2^m) in-place squaring with output-wire relabelling, CX count; lib gf2m_square(m, relabel=True)
    'gf2m_square_relabel_cx': {163: 408, 233: 153, 283: 713, 409: 204, 571: 1538},
    # ANDs of one binary-curve point addition (second point classical, default registers, uncontrolled), exact;
    # lib gf2ec.gf2_point_add_count (the controlled form costs +(2m-1)).  m = 409: 55 006, not 22 * 2501 = 55 022,
    # because some summands of the doubling steps are linear and need no AND.
    'gf2m_pointadd_and': {163: 18160, 233: 29546, 283: 40056, 409: 55006, 571: 92794},
    # logical qubits of the lab binary-curve point addition (k6, division with p = l registers + 3)
    'gf2m_pointadd_qubits': {163: 1959, 233: 3032, 283: 3965, 409: 5320, 571: 8568},
    # narrow 6-bit lookup (k1 bilinear split; lib lookup.narrow_lookup): best / typical on random tables, see lookup.py
    'narrow_lookup': {'loop1_c6_w32': 41, 'loop3_c6_w21': 33},
    # k2: loop2/unloop2 of Gidney 2025 compressed mod p-1, ANDs per prime, l=21, promise S < 2^29 (lab, K_CLAIMS)
    'rsa_loop2_per_prime': 106,
    # k4 + k8: AES-128 oracle (r=1) built and simulated in the lab (work/k8): ANDs, CX, logical qubits
    'aes128_oracle': {'AND': 10046, 'CX': 169714, 'qubits': 1696},
}

# Full point addition of the lab (k4, mixed addition, all exceptional cases, ordinary representation, special prime),
# built with the 0.5.0rc Kaliski (legacy).  (curve, Karatsuba levels) -> (ANDs, logical qubits).  logs/v2/*k4s*.log
K4_POINTADD_LAB = {
    ('p256', 'L3'): (3041421, 6742), ('p256', 'L2'): (3082199, 5730), ('p256', 'L0'): (3293519, 4368),
    ('secp256k1', 'L3'): (3033954, 6934), ('secp256k1', 'L2'): (3074732, 5730), ('secp256k1', 'L0'): (3286052, 4560),
}
# Full point addition of the PROVED scheme: qoptct.crypto.ecpoint.ec_point_add in its default construction (library
# blocks: Kaliski core 12n^2 - 8n - 3, mulsq, fold; all exceptional cases; special prime).  Proved at n = 256 for
# P-256 and secp256k1 (L3) by qoptct.verify.compose.ecadd (PROOF.md, part II).  (curve, Karatsuba levels) ->
# (ANDs, logical qubits); tests/test_estimate_logical.py compares them with ec_point_add(...).counts().
POINTADD_PROVED = {
    ('p256', 'L3'): (3031380, 6724), ('p256', 'L2'): (3072707, 5532), ('p256', 'L0'): (3284271, 4365),
    ('secp256k1', 'L3'): (3023913, 6916), ('secp256k1', 'L2'): (3065240, 5724), ('secp256k1', 'L0'): (3276804, 4557),
}
KALISKI_CALLS_PER_POINTADD = 4      # 2 divisions x (compute + uncompute), cf. Litinski Fig. 5 (4 inversions)


def _gf2m_pointadd_and(m, b):
    """ANDs of one uncontrolled binary-curve point addition: exact table entry, else 2(l+1) products."""
    if m in b['gf2m_pointadd_and']:
        return b['gf2m_pointadd_and'][m]
    return 2 * (it_len(m) + 1) * b['gf2m_mult'][m]


def library_blocks(gf2m_ms=(163, 233, 283, 409, 571)):
    """Recompute the library-backed entries of BLOCKS from qoptct.crypto (about 20 s for all five fields).
    -> dict with the same keys (only the 'lib' entries)."""
    from qoptct.crypto import kaliski_inverse, gf2m_mult_count, gf2m_square
    from qoptct.crypto.fold import CURVES
    from qoptct.crypto.gf2ec import gf2_point_add_count
    out = {'kaliski_special': {}, 'kaliski_special_legacy': {}, 'gf2m_mult': {}, 'gf2m_square_relabel_cx': {},
           'gf2m_pointadd_and': {}}
    for name in ('p256', 'secp256k1'):
        p = CURVES[name]
        out['kaliski_special'][name] = kaliski_inverse(p, reduction='special').counts()['AND']
        out['kaliski_special_legacy'][name] = kaliski_inverse(p, reduction='special', legacy=True).counts()['AND']
    for m in gf2m_ms:
        out['gf2m_mult'][m] = gf2m_mult_count(m)['AND']
        out['gf2m_square_relabel_cx'][m] = gf2m_square(m, relabel=True).counts()['CX']
        out['gf2m_pointadd_and'][m] = gf2_point_add_count(m)['AND']
    return out


# ================================================================================================ 1. GE 2019
# Gidney, Ekera 2019 (arXiv:1905.09749): Sec. 2.8 Eq. (1),(2); ancillary estimate_costs.py (tof_count).  Model: Toffoli
# (Cuccaro adder 2L-1, Babbush lookup 2^c - 1, unlookup 2 sqrt(L)).
def _ge19_ne_rsa(n):
    m = math.ceil(n / 2) - 1
    return m + 2 * (m - 20)


def ge2019_rsa(n=2048, ne=None, cexp=5, cmul=5, csep=1024, dev_off=4):
    """GE 2019 Toffoli count exactly as estimate_costs.py (published blocks).  Table 1/3 parameters:
    1024: (dev_off 5, cmul 5, cexp 5) -> 0.41e9; 2048: (4,5,5) -> 2.7e9; 3072: (6,4,5) -> 9.9e9."""
    ne = _ge19_ne_rsa(n) if ne is None else ne
    pad = int(math.ceil(LG(n * n * ne)) + dev_off)
    pieces = math.ceil(n / csep)
    piece_len = csep + pad
    inner = int(math.ceil(ne * 2 * (n + pad + 2) / (cexp * cmul)))
    c = cexp + cmul
    add = (2 * piece_len - 1) * pieces
    look = 2 ** c - 1
    unl = 2 * math.sqrt(2 ** c - 1)
    q = int(round(3 * n + 0.002 * n * LG(n)))      # abstract circuit, Sec. 2.5 (the lab compared at 6330 with runways)
    return LogicalEstimate(
        'rsa%d' % n, 'published', (add + look + unl) * inner, q, 'Toffoli', 'published',
        'arXiv:1905.09749 Table 1/3, estimate_costs.py', 'Cuccaro adder (2L-1), lookup 2^c-1',
        parts={'adder': add * inner, 'lookup': look * inner, 'unlookup': unl * inner}, extra=dict(inner=inner, pad=pad))


# ================================================================================================ 2. Gidney 2025
# arXiv:2505.15917 Table 5: n: (s, l, w1, w3, w4, f, m, P_deviant, E(shots), Toffolis, qubits)
G25_ROWS = {
    # (the accounting below reproduces the 2-digit Toffoli column of all 7 rows; the lab model k3 dropped the loop3
    #  initial lookup through a duplicate key in a dict literal and missed n=3072 by 2.7 %)
    1024: (8, 18, 6, 3, 6, 28, 640, 0.0287, 9.4, 1.1e9, 742),
    1536: (8, 21, 6, 3, 5, 31, 960, 0.0183, 9.3, 3.1e9, 1074),
    2048: (8, 21, 6, 3, 5, 33, 1280, 0.0125, 9.2, 6.5e9, 1399),
    3072: (8, 21, 6, 3, 5, 35, 1920, 0.0091, 9.2, 1.9e10, 2043),
    4096: (8, 24, 6, 3, 5, 36, 2560, 0.0080, 9.2, 4.0e10, 2692),
    6144: (8, 24, 6, 3, 5, 39, 3840, 0.0042, 9.1, 1.2e11, 3978),
    8192: (8, 24, 6, 3, 5, 40, 5120, 0.0040, 9.1, 2.7e11, 5261),
}


def _g25_account(n, lookup=None, loop2_per_prime=None, len_dlog=None):
    """Exact accounting of Gidney 2025 as in the paper's code (Zenodo 10.5281/zenodo.15347487: _execution_config.py,
    _cost_config.py; read, not run): per subroutine ALL additions are counted at the subroutine's maximal addition
    width, all lookups with N > 2 entries at the maximal N: N - bitlen(N) - 1; phase lookups ~ 0; times
    (shots + 1)/keep_rate.  Note (paper vs code): in loop3/unloop3 the lookup address has 2*w3 bits.
    lookup(sub, N, w) replaces the cost of an N-entry, w-bit lookup in subroutine `sub` (None: as in the paper)."""
    s, l, w1, w3, w4, f, m, *_ = G25_ROWS[n]
    x_size = math.ceil(n * (1 / 2 + 1 / (2 * s)))
    y_size = math.ceil(n / (2 * s))
    assert x_size + y_size == m
    P = math.ceil(math.ceil(m / w1) * n / l)
    Ld = len_dlog or (l + m.bit_length())
    nw1, nw3a, nw3b, nw4 = math.ceil(m / w1), math.ceil(l / w3), math.ceil(l / w3), math.ceil(l / w4)
    inner = max(0, nw3a - 2) * nw3b
    C = {'loop1': Counter({('lookup', 1 << w1, Ld): nw1, ('add', Ld): nw1})}
    c2 = Counter()
    for k in range(l, Ld):
        c2[('add', k)] += 1
        c2[('add', k + 1)] += 1
        c2[('lookup', 2, k)] += 1
    C['loop2'] = c2
    C['unloop2'] = Counter(c2)
    C['loop3'] = Counter({('lookup', 1 << (2 * w3), l + 1): 1, ('add', l + 1): inner, ('add', l): inner})
    C['loop3'][('lookup', 1 << (w3 + w3), l + 1)] += inner
    C['loop3'][('lookup', 2, l)] += inner
    C['loop4'] = Counter({('lookup', 1 << w4, f + 1): nw4, ('add', f + 1): nw4, ('lookup', 2, f): nw4, ('add', f): nw4})
    C['loop4'][('lookup', 1 << w4, f)] += nw4 * 0.5
    C['loop4'][('cmp', f)] += nw4 * 0.5
    C['loop4'][('phase', 1 << w4)] += nw4
    C['unloop3'] = Counter({('phase', 1 << (2 * w3)): 1, ('add', l + 1): inner * 2, ('add', l): inner * 2,
                            ('lookup', 1 << (2 * w3), l + 1): inner * 2, ('lookup', 2, l): inner * 2,
                            ('cmp', l): inner, ('lookup', 1 << (2 * w3), l): inner})
    C['unloop3'][('phase', 1 << (2 * w3))] += inner * 2
    mult = {'loop1': P + 1, 'loop2': P, 'loop3': P, 'loop4': P, 'unloop3': P, 'unloop2': P}

    def ph(N):
        b = N.bit_length()
        n1 = b // 2
        n2 = b - n1
        return max(0, 2 ** (n1 // 2 - n1 - 1)) + max(0, 2 ** (n2 // 2 - n2 - 1))
    parts, per_shot = {}, 0
    for sub, cnt in C.items():
        sc = {k: round(v * mult[sub]) for k, v in cnt.items() if v}
        if sub in ('loop2', 'unloop2') and loop2_per_prime is not None:
            parts[sub] = P * loop2_per_prime
            per_shot += parts[sub]
            continue
        add_n = max([k[1] for k in sc if k[0] in ('add', 'cmp')], default=0)
        adds = sum(v for k, v in sc.items() if k[0] in ('add', 'cmp'))
        lk = [(k, v) for k, v in sc.items() if k[0] == 'lookup' and k[1] > 2]
        Nmax = max([k[1] for k, _ in lk], default=0)
        phs = [(k, v) for k, v in sc.items() if k[0] == 'phase' and k[1] > 4]
        Pmax = max([k[1] for k, _ in phs], default=0)
        tot = adds * add_n
        for k, v in lk:
            tot += v * (lookup(sub, Nmax, k[2]) if lookup else max(0, Nmax - Nmax.bit_length() - 1))
        tot += sum(v for _, v in phs) * (ph(Pmax) if Pmax else 0)
        parts[sub] = tot
        per_shot += tot
    # probability of the deviation (mask) failure: _execution_config.probability_of_deviation_failure
    adds_acc = nw4 * P
    eps = adds_acc * 3 * 2.0 ** -f
    mask_bits = round(math.log2(eps) / 2 + f)
    dev = adds_acc * (2 * 2.0 ** -f + 2.0 ** -f)
    mp = 2.0 ** (mask_bits - f)
    pdev = mp + dev / mp
    keep = 1 - pdev
    return dict(tof=math.ceil(per_shot * (s + 1) / keep), per_shot=per_shot, P=P, pdev=pdev,
                eshots=(s + 1) / keep / 0.99, qubits=m + Ld + l + f + max(l, f, Ld),
                parts={k: v * (s + 1) / keep for k, v in parts.items()})


def gidney2025_rsa(n=2048):
    """Gidney 2025 as in the paper's code; reproduces Table 5 (Toffoli, P_deviant, E(shots), qubits) for the 7 rows."""
    x = _g25_account(n)
    return LogicalEstimate(
        'rsa%d' % n, 'published', x['tof'], x['qubits'], 'AND', 'published',
        'arXiv:2505.15917 Table 5 (+ Zenodo code)', 'addition n-1.., lookup N - bitlen(N) - 1 as in the code',
        parts=x['parts'], extra=dict(P=x['P'], p_deviant=x['pdev'], shots=x['eshots'], per_shot=x['per_shot']))


def gidney2025_rsa_table3(n=2048):
    """Cruder Table-3 accounting (addition r-1, lookup 2^w-w-1, phaseup sqrt 2^w), |P| = n m/(l w1) unrounded;
    reproduces Table 5 within 2.5 % (6.57e9 at n=2048)."""
    s, l, w1, w3, w4, f, m, pdev, eshots, _, _ = G25_ROWS[n]
    lenm = m.bit_length()
    P = n * m / (l * w1)
    W1, W3, W4 = math.ceil(m / w1), math.ceil(l / w3), math.ceil(l / w4)
    add = lambda r: r - 1
    lookup = lambda w: 2 ** w - w - 1
    phaseup = lambda w: math.sqrt(2 ** w)
    rows = {
        'loop1': ((P + 1) * W1, l + lenm, w1, 1, 1, 0),
        'loop2': (P * lenm, l + lenm, 0, 2, 0, 0),
        'loop3_start': (P, l, 2 * w3, 0, 1, 0),
        'loop3_body': (P * (W3 - 2) * W3, l, 2 * w3, 2, 1, 0),
        'loop4': (P * W4, f, w4, 1.5, 2.5, 1),
        'unloop3_body': (P * (W3 - 2) * 2 * W3, l, 2 * w3, 2.5, 1.5, 1),
        'unloop3_clean': (P, l, 2 * w3, 0, 0, 1),
        'unloop2': (P * lenm, l + lenm, 0, 2, 0, 0),
    }
    per_shot, parts = 0.0, {}
    for k, (it, r, w, na, nl, nph) in rows.items():
        parts[k] = it * (na * add(r) + (nl * lookup(w) if nl else 0) + (nph * phaseup(w) if nph else 0))
        per_shot += parts[k]
    e = (s + 1) / (1 - pdev) / 0.99
    return LogicalEstimate('rsa%d' % n, 'published', per_shot * e, m + 2 * f + 2 * l + lenm, 'AND', 'published',
                           'arXiv:2505.15917 Table 3 + Table 5', 'rough accounting; qubits = Table 5 formula',
                           parts={k: v * e for k, v in parts.items()}, extra=dict(shots=e))


def rsa_g25_qoptct(n=2048, blocks=None, narrow_loop3=None):
    """Gidney 2025 framework with checked blocks (the audited positive variant; 2048: 4.91e9 Toffoli @ 1400 qubits,
    -24.2 % against 6.49e9 @ 1399):
      * 6-bit table lookups (loop1 word 32 bits, loop3/unloop3 word 21 bits) by the bilinear narrow lookup of
        qoptct.crypto.lookup: 41 resp. 33 ANDs instead of the code's N - bitlen(N) - 1 = 56;  (+4 qubits)
      * loop2/unloop2 compressed modulo p-1 in 106 ANDs per prime instead of 22 additions x 32 (k2, promise S < 2^29);
      * accumulator of the discrete logarithms l + ceil(lg W1) = 29 bits instead of l + len(m) = 32.
    CAVEATS (audit k10): the block constants are checked (prove/kprove at c<=6 incl. c=6 w=21/32, check_sim, canaries)
    but the end-to-end circuit is NOT assembled, and the loop2 -> loop3 hand-over (the output of the compression is
    only == S mod p-1 and < 2^21) is argued in words, not by a circuit.  In Gidney's own metric q^3 t the gain is
    smaller (about -11 % instead of -24 %).  The 33 for w=21 is the best case of the random search (typical 33.3:
    33 in 15/20 tables, 34 in 4/20, 35 in 1/20 -- pass `narrow_loop3=34` to see the effect).  Constants were
    measured at l=21; for other n they are an extrapolation (note).  Not a record: Sharma-Lau 2026 and others
    are better at more qubits; this is a new point of the Pareto front near 1400 qubits."""
    b = dict(BLOCKS, **(blocks or {}))
    n1 = b['narrow_lookup']['loop1_c6_w32']
    n3 = narrow_loop3 if narrow_loop3 is not None else b['narrow_lookup']['loop3_c6_w21']
    l, w1 = G25_ROWS[n][1], G25_ROWS[n][2]

    def lk(sub, N, w):
        src = max(0, N - N.bit_length() - 1)
        if N == 64:
            return min(src, n1 if sub == 'loop1' else (n3 if sub in ('loop3', 'unloop3') else src))
        return src
    m = G25_ROWS[n][6]
    Ld = l + math.ceil(LG(math.ceil(m / w1)))
    x = _g25_account(n, lookup=lk, loop2_per_prime=b['rsa_loop2_per_prime'], len_dlog=Ld)
    base = _g25_account(n)
    note = 'model from checked blocks, no end-to-end circuit; q^3 t gain smaller'
    if l != 21:
        note += '; block constants measured at l=21 (extrapolation for l=%d)' % l
    return LogicalEstimate(
        'rsa%d' % n, 'qoptct', x['tof'], x['qubits'] + 4, 'AND', 'qoptct',
        'arXiv:2505.15917 framework; blocks: qoptct.crypto.lookup (narrow lookup), k2 loop2, 29-bit accumulator', note,
        parts=x['parts'], extra=dict(gain=1 - x['tof'] / base['tof'], baseline=base['tof'], shots=x['eshots'],
                                     narrow=(n1, n3), loop2_per_prime=b['rsa_loop2_per_prime'], acc_bits=Ld))


# ================================================================================================ 3. ECDLP, prime curves
# Litinski 2023 (arXiv:2306.08585): Fig.3 (elementary), Fig.5 (point addition), Fig.8 (modular blocks), Fig.9c.
LIT_PUB = dict(
    nctof=lambda n: n, mod_add=lambda n: 4 * n, cmod_add=lambda n: 5 * n, mod_sub=lambda n: 6 * n,
    cmod_sub=lambda n: 7 * n, mod_neg=lambda n: 2 * n, cmod_neg=lambda n: 3 * n, mod_dbl=lambda n: 2 * n,
    mod_mul=lambda n: 2.25 * n * n + 9 * n, mod_sqr=lambda n: 2.25 * n * n + 9 * n,
    mod_inv=lambda n: 26 * n * n + 2 * n, lookup=lambda K: K, unlookup=lambda K: 3 * math.sqrt(K))
LIT_PA_COUNTS = dict(nctof=18, mod_add=3, cmod_add=2, mod_sub=2, cmod_sub=4, mod_neg=2, cmod_neg=1,
                     mod_dbl=2, mod_mul=8, mod_sqr=2, mod_inv=4)       # Fig.5: 10 multiplications, 2 of them squares


def _lit_pointadd(n, k=1):
    parts = {}
    for b, c in LIT_PA_COUNTS.items():
        if b == 'mod_inv':     # Fig.9b: k parallel inversions = (3k-3) multiplications + 1 inversion
            parts[b] = c * ((3 * k - 3) * LIT_PUB['mod_mul'](n) + LIT_PUB['mod_inv'](n)) / k
        else:
            parts[b] = c * LIT_PUB[b](n)
    return sum(parts.values()), parts


def litinski2023_ecdlp(n=256, k=1, w=16, reps=13):
    """Toffolis per key: reps x (lookup 2^w + point addition + unlookup) (Sec. 2: 13 repetitions of 16-bit windows
    for 256 bits); k = number of parallel instances sharing the inversion (Fig.9c: (44k+65)e6 / k).  Logical qubits:
    the 3000 of the paper (11 registers of 256 = 2816 plus margin)."""
    pa, parts = _lit_pointadd(n, k)
    K = 2 ** w
    parts = {b: reps * v for b, v in parts.items()}
    parts['lookup'] = reps * LIT_PUB['lookup'](K)
    parts['unlookup'] = reps * LIT_PUB['unlookup'](K)
    return LogicalEstimate('ecdlp-p%d' % n if n in (256, 384, 521) else 'ecdlp-%d' % n, 'published',
                           reps * (pa + LIT_PUB['lookup'](K) + LIT_PUB['unlookup'](K)), 3000 if n == 256 else 11 * n,
                           'AND', 'published', 'arXiv:2306.08585 Fig.5,8,9c', 'exact point addition (all cases)',
                           parts=parts, extra=dict(pointadd=pa, k=k))


def pointadd_ands(curve='p256', variant='L3', kaliski='current'):
    """ANDs of the full point addition.
    'current' (default): the real, PROVED scheme (``POINTADD_PROVED``: ec_point_add on the library blocks).
    'legacy': the number of the lab logs (k4 construction, ``legacy=True``, 0.5.0rc Kaliski 12n^2 + n).
    'estimate': the former estimate = lab number - 4 x (legacy - current Kaliski) = 4 x (9n + 3) = 9228 at n=256
    (kept for comparison: the real scheme is within 10 ANDs of it, see ``pointadd_ands_with_mulsq``)."""
    if kaliski == 'current':
        return POINTADD_PROVED[(curve, variant)][0]
    pa, _ = K4_POINTADD_LAB[(curve, variant)]
    if kaliski == 'legacy':
        return pa
    d = BLOCKS['kaliski_special_legacy'][curve] - BLOCKS['kaliski_special'][curve]
    return pa - KALISKI_CALLS_PER_POINTADD * d


def pointadd_ands_with_mulsq(curve='p256', variant='L3', mul_sq_delta=None):
    """Optional: also take the library's lowered multiplier/square counts (8 multiplications and 2 squares per
    point addition, Litinski Fig.5).  mul_sq_delta = (d_mul, d_sq) per call, default for Karatsuba levels 3 (n=256):
    (-54, -187) = (37826 - 37880, 22145 - 22332) [L2: (-18, -61), L0: (-2, -5)].  An ESTIMATE (upper bound on the gain)."""
    dm, ds = mul_sq_delta or {'L3': (-54, -187), 'L2': (-18, -61), 'L0': (-2, -5)}[variant]
    return pointadd_ands(curve, variant, 'estimate') + 8 * dm + 2 * ds


def _windowed_total(pa, w=16, reps=13):
    K = 2 ** w
    return reps * (pa + K + 3 * math.sqrt(K))


def ecdlp_prime_qoptct(curve='p256', variant='L3', blocks=None, kaliski='current'):
    """Litinski framework (13 windows of 16 bits, k=1) with the full exact point addition of qoptct (all exceptional
    cases, special prime: P-256 or secp256k1): the real scheme ``ecpoint.ec_point_add`` (default construction), proved at
    n = 256 by ``verify.compose.ecadd`` (classical Q and the window w = 4; the numbers of ``POINTADD_PROVED``).
    variant: Karatsuba levels 'L3' (P-256: 6724 qubits, 4.03e7 Toffoli), 'L2', 'L0' (4365 qubits, 4.36e7).
    kaliski='legacy': the lab k4 construction (``legacy=True``, not proved).  NOT a record: the approximate point
    additions of 2026 (Haner et al. 2609.05625, Schrottenloher 2606.02235, Initio) are cheaper and use 3-5 times fewer
    qubits; the gain is against Litinski's EXACT scheme.  The point addition is built classically-controlled (Q
    classical); in the windowed scheme Q comes from the lookup (+O(n), about +0.1 %); the whole ECDLP circuit is not
    assembled or proved, only the point addition.  `blocks` may carry the Kaliski counts of library_blocks()."""
    b = dict(BLOCKS, **(blocks or {}))
    pa0, q0 = K4_POINTADD_LAB[(curve, variant)]
    if kaliski == 'legacy':
        pa, q = pa0, q0
    else:
        pa, q = POINTADD_PROVED[(curve, variant)]
    base = litinski2023_ecdlp(256)
    tot = _windowed_total(pa)
    return LogicalEstimate(
        'ecdlp-%s' % curve, 'qoptct', tot, q, 'AND', 'qoptct', 'arXiv:2306.08585 framework; point addition: ec_point_add',
        'full exact point addition %s (%d ANDs = %.1f n^2); %s; end-to-end circuit not assembled; approximate schemes '
        'of 2026 are better' % (variant, pa, pa / 65536, 'lab k4 construction, not proved' if kaliski == 'legacy'
                                else 'proved at n=256 for the classical Q and the window w=4 (compose.ecadd; L3)'),
        parts={'pointadd': 13 * pa, 'lookup': 13 * 2 ** 16, 'unlookup': 13 * 3 * 256},
        extra=dict(pointadd=pa, pointadd_legacy_kaliski=pa0, variant=variant, gain=1 - tot / base.toffoli,
                   kaliski_per_call=b['kaliski_special'][curve]))


# ================================================================================================ 4. ECDLP, other published
# Haner, Jaques, Naehrig, Roetteler, Soeken 2020 (arXiv:2001.09580), Low-T variant, Table 1 and Table 4 (T = 4 AND).
def haner2020_ecdlp(n=256):
    T = 1115 * n ** 3 / LG(n) - 1.08 * 2 ** 24
    return LogicalEstimate('ecdlp-%d' % n, 'published', T / 4, int(10 * n + 7.4 * math.floor(LG(n)) + 1.3), 'T/4',
                           'published', 'arXiv:2001.09580 Table 1, Table 4 (Low-T)', 'AND = 4 T, AND-dagger = 0 T')


# Roetteler, Naehrig, Svore, Lauter 2017 (arXiv:1706.06752): (448 lg n + 4090) n^3, 9n + 2 ceil(lg n) + 10 qubits.
RNSL_T2 = {110: (1014, 9.44e9), 160: (1466, 2.97e10), 192: (1754, 5.30e10), 224: (2042, 8.43e10),
           256: (2330, 1.26e11), 384: (3484, 4.52e11), 521: (4719, 1.14e12)}


def roetteler2017_ecdlp(n=256):
    return LogicalEstimate('ecdlp-%d' % n, 'published', (448 * LG(n) + 4090) * n ** 3, 9 * n + 2 * math.ceil(LG(n)) + 10,
                           'Toffoli', 'published', 'arXiv:1706.06752 Table 2, Sec. 5',
                           'Toffoli with Toffoli erasure, minimal qubits; blocks of qoptct are NOT substitutable')


# ================================================================================================ 5. ECDLP, binary curves
# Banegas, Bernstein, van Hoof, Lange 2020 (TCHES 2021(1)): Table 2 (division: GCD, FLT), Table 5 (total).
BBHL_T2 = {8: (3641, 243), 16: (10403, 1053), 127: (277195, 50255), 163: (442161, 83353),
           233: (827977, 132783), 283: (1202987, 236279), 571: (4461673, 814617)}     # n: (GCD TOF, FLT TOF)
BBHL_T5 = {8: (68, 7360, 132480), 16: (125, 21016, 714544), 127: (904, 559141, 143140096),
           163: (1157, 893585, 293095880), 233: (1647, 1669299, 781231932),
           283: (1998, 2427369, 1378745592), 571: (4015, 8987401, 10281586744)}      # n: (qubits, step TOF, total TOF)


def banegas2020_ecdlp(n=163):
    """(2n+2) steps x (2 GCD divisions + 2 MODMULT + 3n) (Alg. 3, Table 5); MODMULT inferred from Table 5."""
    gcd, _ = BBHL_T2[n]
    q, step, tot = BBHL_T5[n]
    mm = (step - 2 * gcd - 3 * n) / 2
    return LogicalEstimate('ecdlp-b%d' % n, 'published', (2 * n + 2) * (2 * gcd + 2 * mm + 3 * n), q, 'Toffoli',
                           'published', 'TCHES 2021(1) Table 2, Table 5', 'Toffoli model, no windows',
                           extra=dict(modmult=mm))


# Taguchi, Takayasu (ePrint 2023/553, QIP 2024) Table 6 (no windows) and Table 7 (windows), multiplier of Kim et al.
KIM_M = {163: 992, 233: 1441, 283: 1784, 571: 3813}                          # ePrint 2022/1095 CRT multiplier
TT23_T6 = {163: 13175432, 233: 30000204, 283: 49121208, 571: 228787416}
TT23_T7 = {163: 1781025, 233: 3679975, 283: 5765145, 571: 23390601}          # best windows 9, 9, 10, 11
TT23_Q = {163: 1957, 233: 3030, 283: 3963, 571: 8566}
# Itoh-Tsujii: multiplications in the addition chain for m-1 (for 409 and 571 shorter than the binary method)
IT_L = {163: 9, 233: 10, 283: 11, 409: 10, 571: 12}


def it_len(m):
    return IT_L.get(m, (m - 1).bit_length() - 1 + bin(m - 1).count('1') - 1)


def taguchi_takayasu2023_ecdlp(m=163):
    """No windows: (2m+2) controlled point additions, each (4l+4) multiplications (FLT-type division x2 + 2 MODMULT)
    + 3m.  Reproduces Table 6."""
    l = it_len(m)
    pa = (4 * l + 4) * KIM_M[m] + 3 * m
    return LogicalEstimate('ecdlp-b%d' % m, 'published', (2 * m + 2) * pa, TT23_Q[m], 'Toffoli', 'published',
                           'ePrint 2023/553 Table 6 (+ Kim 2022/1095)', 'Toffoli model, no windows, no exceptions',
                           extra=dict(l=l, pointadd=pa, windowed_table7=TT23_T7[m]))


def ecdlp_gf2m_qoptct(m=163, blocks=None):
    """Binary curve, NO windows (positive variant, -53..-54 % against Taguchi-Takayasu Table 6): (2m+2) controlled
    point additions, each 2(l+1) multiplications (two divisions by the Itoh-Tsujii chain with the register erased by
    X-measurement: (l+1) multiplications each) + (2m-1) ANDs for the control.  Multiplier = qoptct.crypto.gf2m
    (CRT bilinear rank 908/1343/1669/2501/3569).  Like TT/BBHL: exceptional cases x1 = x2, P1 + P2 = +-P2 are NOT
    handled (so it is not comparable with the exact Garn-Kan scheme).  The erase trick for GF(2^m) products is
    published (Khattar et al. 2510.10967); new is its use in the chain.  Qubits: lab k6 (p = l registers)."""
    b = dict(BLOCKS, **(blocks or {}))
    M = b['gf2m_mult'][m]
    l = it_len(m)
    pa = _gf2m_pointadd_and(m, b) + (2 * m - 1)
    tot = (2 * m + 2) * pa
    ref = taguchi_takayasu2023_ecdlp(m).toffoli if m in KIM_M else None
    return LogicalEstimate(
        'ecdlp-b%d' % m, 'qoptct', tot, b['gf2m_pointadd_qubits'][m], 'AND', 'qoptct',
        'ePrint 2023/553 framework; blocks: qoptct.crypto.gf2m', 'no windows; exceptional cases not handled; '
        'end-to-end circuit not assembled (blocks: gate-level m=163,233; kprove m<=9)',
        parts={'multiplications': (2 * m + 2) * _gf2m_pointadd_and(m, b), 'control': (2 * m + 2) * (2 * m - 1)},
        extra=dict(l=l, multiplier=M, pointadd=pa, gain=(1 - tot / ref) if ref else None))


def ecdlp_gf2m_windowed_estimate(m=163, blocks=None, conservative=False):
    """ESTIMATE ONLY (windowed scheme not built; the window model is the lab's, not Taguchi-Takayasu's).
    k6 model: N = 2 ceil((m+1)/s) point additions of 2(l+1)M + 3 lookups (2^s - s - 1); conservative=True: as in TT
    (QIP 2024 Sec. 6.2): N + 12 additions, lookup 2(2^s - 1).  Best s is searched."""
    b = dict(BLOCKS, **(blocks or {}))
    M = b['gf2m_mult'][m]
    l = it_len(m)
    best = None
    for s in range(6, 16):
        N = 2 * math.ceil((m + 1) / s) + (12 if conservative else 0)
        lk = 2 * (2 ** s - 1) if conservative else 2 ** s - s - 1
        t = N * (_gf2m_pointadd_and(m, b) + 3 * lk)
        if best is None or t < best[0]:
            best = (t, s, N)
    t, s, N = best
    return LogicalEstimate('ecdlp-b%d' % m, 'qoptct', t, b['gf2m_pointadd_qubits'][m], 'AND', 'estimate',
                           'ePrint 2023/553 Table 7', 'ESTIMATE: windows not built or checked; window s=%d, %d additions' % (s, N),
                           extra=dict(window=s, additions=N, conservative=conservative,
                                      published_table7=TT23_T7.get(m)))


# ================================================================================================ 6. AES (Grover oracle)
# Jaques, Naehrig, Roetteler, Virdia 2020 (arXiv:1910.01700): Table 1 (S-box 34 AND), Table 8 (oracle), Table 10 (Grover)
AES_PARAMS = {128: (10, 10), 192: (12, 8), 256: (14, 13)}   # (rounds, SubWord calls in the key schedule)
JNRV_T8 = {(128, 1): (54908, 13727, 1665), (192, 1): (61436, 15359, 1985), (256, 1): (75580, 18895, 2305),
           (128, 2): (109820, 27455, 3329), (192, 2): (122876, 30719, 3969), (256, 2): (151164, 37791, 4609),
           (256, 3): (226748, 56687, 6913)}                  # (T, #AND, width), in-place MixColumns


def _aes_oracle_ands(k, r, fwd=34, bwd_state=34, bwd_key=34):
    R, sw = AES_PARAMS[k]
    ns, nk = 16 * R, 4 * sw
    return r * (ns * fwd + nk * fwd + ns * bwd_state + nk * bwd_key) + (128 * r - 1)


def jaques2020_aes(k=128, r=1, sbox_and=34):
    """Oracle with r plaintext-ciphertext pairs: r AES evaluations (forward and backward, 34 AND per S-box) + (128 r -
    1) ANDs for the comparison; Grover: pi/4 2^(k/2) iterations.  Reproduces Table 8 (T = 4 AND) and Table 10."""
    orc = _aes_oracle_ands(k, r, sbox_and, sbox_and, sbox_and)
    gi = math.pi / 4 * 2 ** (k / 2)
    return LogicalEstimate('aes%d' % k, 'published', orc * gi, JNRV_T8[(k, r)][2] if (k, r) in JNRV_T8 else 0, 'AND',
                           'published', 'arXiv:1910.01700 Table 8, 10', 'oracle ANDs per iteration: %d' % orc,
                           extra=dict(oracle=orc, iterations=gi, r=r))


def aes_oracle_qoptct(k=128, r=1, blocks=None):
    """AES-128 oracle (r=1) as built in the lab: S-box 32 AND forward (known: Meuli et al. 2019), backward pass of
    the rounds through X-measurement of the round register + a phase oracle of the S-box (14 AND), 10 046 ANDs
    against 13 727 (Jaques) and 12 927 (Jaques with Meuli's S-box); 169 714 CNOT; 1696 qubits.  NOT independently
    verified (checked by the author's logical simulator with phases and check_sim of the full QASM + one canary).
    AES-192/256 and r > 1: formula only (status 'estimate'; r=2 gets the formula of k8 with S-box 32/14/32)."""
    b = dict(BLOCKS, **(blocks or {}))
    gi = math.pi / 4 * 2 ** (k / 2)
    if (k, r) == (128, 1):
        a = b['aes128_oracle']
        return LogicalEstimate(
            'aes128', 'qoptct', a['AND'] * gi, a['qubits'], 'AND', 'qoptct', 'arXiv:1910.01700 framework; oracle: lab k8',
            'NOT independently verified; oracle built, simulated; 13727 -> %d ANDs (-%.1f %%), CNOT %d' %
            (a['AND'], 100 * (1 - a['AND'] / 13727), a['CX']),
            extra=dict(oracle=a['AND'], cx=a['CX'], iterations=gi,
                       vs_jaques_with_meuli_sbox=1 - a['AND'] / _aes_oracle_ands(128, 1, 32, 32, 32)))
    orc = _aes_oracle_ands(k, r, 32, 14, 32)
    return LogicalEstimate('aes%d' % k, 'qoptct', orc * gi, JNRV_T8[(k, r)][2], 'AND', 'estimate',
                           'arXiv:1910.01700 framework', 'ESTIMATE: formula only, oracle not built', extra=dict(oracle=orc))


# ================================================================================================ reproduction of the papers
def reproduce_published():
    """Check every published model against the numbers printed in its paper.
    -> list of (name, got, want, tolerance, ok, source).  Tolerances: table rounding (GE19 2 %; Gidney 2025: Table 5
    prints 2 digits and the model hits them exactly for all 7 rows); Litinski 0.6 % ('109M' rounding); Haner
    2 %; Roetteler 2.5 % (their own regression), qubits exact; Banegas/Taguchi-Takayasu/Jaques exact (Jaques Table
    10: 0.5 %)."""
    rows = []

    def chk(name, got, want, rel, src):
        rows.append((name, got, want, rel, _close(got, want, rel), src))
    for n, do, cm, ce, want in [(1024, 5, 5, 5, 0.41e9), (2048, 4, 5, 5, 2.7e9), (3072, 6, 4, 5, 9.9e9)]:
        chk('GE19 RSA-%d Toffoli' % n, ge2019_rsa(n, cexp=ce, cmul=cm, dev_off=do).toffoli, want, 0.02, 'arXiv:1905.09749 Table 1/3')
    for n, row in G25_ROWS.items():
        r = gidney2025_rsa(n)
        got = float('%.2g' % r.toffoli)
        chk('G25 RSA-%d Toffoli' % n, got, row[9], 0.0, 'arXiv:2505.15917 Table 5 (printed to 2 digits)')
        chk('G25 RSA-%d P_deviant %%' % n, round(r.extra['p_deviant'] * 100, 2), round(row[7] * 100, 2), 0.0, 'Table 5')
        chk('G25 RSA-%d E(shots)' % n, round(r.extra['shots'], 1), row[8], 0.0, 'Table 5')
        chk('G25 RSA-%d qubits' % n, r.qubits, row[10], 0.0, 'Table 5')
        chk('G25 RSA-%d Toffoli (Table 3)' % n, gidney2025_rsa_table3(n).toffoli, row[9], 0.025, 'Table 3 + 5')
    chk('Litinski P-256 point addition', litinski2023_ecdlp(256).extra['pointadd'], 126.5 * 256 ** 2 + 189 * 256, 0.0,
        'arXiv:2306.08585 Fig.5,8,9c')
    for k in (1, 2, 4):
        chk('Litinski P-256 key k=%d' % k, litinski2023_ecdlp(256, k=k).toffoli, (44 * k + 65) * 1e6 / k, 0.006, 'Fig.9c')
    for n, want in [(256, 1.08 * 2 ** 31), (384, 1.74 * 2 ** 32), (521, 1.00 * 2 ** 34)]:
        chk('Haner Low-T n=%d T' % n, haner2020_ecdlp(n).toffoli * 4, want, 0.02, 'arXiv:2001.09580 Table 1')
    for n, (q, t) in RNSL_T2.items():
        chk('RNSL n=%d Toffoli' % n, roetteler2017_ecdlp(n).toffoli, t, 0.025, 'arXiv:1706.06752 Table 2')
        chk('RNSL n=%d qubits' % n, roetteler2017_ecdlp(n).qubits, q, 0.0, 'Table 2')
    for n in BBHL_T5:
        chk('BBHL n=%d total TOF' % n, banegas2020_ecdlp(n).toffoli, BBHL_T5[n][2], 0.0, 'TCHES 2021(1) Table 5')
    for m in TT23_T6:
        chk('TT23 m=%d no windows' % m, taguchi_takayasu2023_ecdlp(m).toffoli, TT23_T6[m], 0.0, 'ePrint 2023/553 Table 6')
    for (k, r), (T, M, w) in JNRV_T8.items():
        chk('JNRV AES-%d r=%d oracle T' % (k, r), jaques2020_aes(k, r).extra['oracle'] * 4, T, 0.0, 'arXiv:1910.01700 Table 8')
    chk('JNRV AES-128 r=1 Grover T', jaques2020_aes(128, 1).toffoli * 4, 1.32 * 2 ** 79, 0.005, 'Table 10')
    chk('JNRV AES-256 r=2 Grover T', jaques2020_aes(256, 2).toffoli * 4, 1.81 * 2 ** 144, 0.005, 'Table 10')
    return rows


# ================================================================================================ attack registry, report
# name -> (title, published function, qoptct function(variant), physical model, shots, qubit note)
def _rsa_pub(_v):
    return gidney2025_rsa(2048)


def _rsa_q(_v):
    return rsa_g25_qoptct(2048)


def _ec_pub(_v):
    return litinski2023_ecdlp(256)


def _ec_q(curve):
    return lambda v: ecdlp_prime_qoptct(curve, v or 'L3')


def _bin_pub(m):
    return lambda _v: taguchi_takayasu2023_ecdlp(m)


def _bin_q(m):
    return lambda _v: ecdlp_gf2m_qoptct(m)


def _aes_pub(_v):
    return jaques2020_aes(128, 1)


def _aes_q(_v):
    return aes_oracle_qoptct(128, 1)


ATTACKS = {
    'rsa2048-g25': dict(title='RSA-2048 factoring, Gidney 2025 framework', pub=_rsa_pub, q=_rsa_q, phys='gidney2025'),
    'ecdlp-p256': dict(title='ECDLP P-256, Litinski 2023 framework (exact point addition)', pub=_ec_pub,
                       q=_ec_q('p256'), phys='litinski2023', variants=('L3', 'L2', 'L0')),
    'ecdlp-k1': dict(title='ECDLP secp256k1, Litinski 2023 framework (exact point addition)', pub=_ec_pub,
                     q=_ec_q('secp256k1'), phys='litinski2023', variants=('L3', 'L2', 'L0')),
    'ecdlp-b163': dict(title='ECDLP B-163, Taguchi-Takayasu 2023, no windows', pub=_bin_pub(163), q=_bin_q(163), phys='litinski2023'),
    'ecdlp-b233': dict(title='ECDLP B-233, Taguchi-Takayasu 2023, no windows', pub=_bin_pub(233), q=_bin_q(233), phys='litinski2023'),
    'ecdlp-b283': dict(title='ECDLP B-283, Taguchi-Takayasu 2023, no windows', pub=_bin_pub(283), q=_bin_q(283), phys='litinski2023'),
    'ecdlp-b571': dict(title='ECDLP B-571, Taguchi-Takayasu 2023, no windows', pub=_bin_pub(571), q=_bin_q(571), phys='litinski2023'),
    'aes128': dict(title='AES-128 Grover oracle, Jaques et al. 2020', pub=_aes_pub, q=_aes_q, phys=None),
}
PHYSICAL_SHOTS = {'rsa2048-g25': 9.2}


# Reaction depth of the elliptic-curve runs: Litinski 2023 bounds the depth by 2 x Toffoli (Sec. 2.3 of the paper;
# reaction-limited run time).  Without it the run time is only limited by the factories and the gain is overstated.
ECDLP_DEPTH_PER_TOFFOLI = 2.0


def _physical(model, toffoli, qubits, shots=None, reaction_limit=True):
    from . import physical as ph
    kw = {'shots': shots} if model == 'gidney2025' and shots else {}
    depth = ECDLP_DEPTH_PER_TOFFOLI * toffoli if (model == 'litinski2023' and reaction_limit) else None
    r = ph.estimate(model, toffoli, int(qubits), depth, **kw)
    return dict(model=r.model, physical_qubits=r.physical_qubits, code_distance=r.code_distance,
                runtime_s=r.runtime_s, expected_runtime_s=r.expected_runtime_s, failure_prob=r.failure_prob,
                megaqubit_days=r.megaqubit_days, notes=list(r.notes))


def attack_report(name, variant=None, source='both', physical=True, reaction_limit=True):
    """-> dict(name, title, rows=[LogicalEstimate.to_dict() + 'physical' ...]).  source: 'published' | 'qoptct' |
    'both'.  The physical estimate uses qoptct.estimate.physical: gidney2025 for RSA (9.2 shots, as in Table 5),
    litinski2023 for elliptic curves (a generic space-time model from the Toffoli and qubit count; for the binary
    curves and for the 'published' B-curve rows it is an extrapolation of Litinski's architecture); none for AES.
    For the elliptic curves the reaction depth is 2 x Toffoli (Litinski's bound), the run is reaction limited;
    reaction_limit=False drops the limit (only the factories bound the run time; the gain is then overstated)."""
    if name not in ATTACKS:
        raise ValueError('unknown attack %r (known: %s)' % (name, ', '.join(sorted(ATTACKS))))
    a = ATTACKS[name]
    if variant is not None and variant not in a.get('variants', ()):
        raise ValueError('attack %r has no variant %r (variants: %s)' % (name, variant, ', '.join(a.get('variants', ())) or 'none'))
    if source not in ('published', 'qoptct', 'both'):
        raise ValueError("source must be 'published', 'qoptct' or 'both'")
    rows = []
    for src, fn in (('published', a['pub']), ('qoptct', a['q'])):
        if source not in (src, 'both'):
            continue
        r = fn(variant).to_dict()
        r['physical'] = None
        if physical and a['phys']:
            r['physical'] = _physical(a['phys'], r['toffoli'], r['qubits'], PHYSICAL_SHOTS.get(name), reaction_limit)
        rows.append(r)
    return dict(name=name, title=a['title'], variant=variant, rows=rows)


def _fmt(x):
    if x >= 1e4:
        e = int(math.floor(math.log10(x)))
        return '%.3fe%d' % (x / 10 ** e, e)
    return '%.0f' % x


def _fmt_time(s):
    if s >= 2 * 86400:
        return '%.2f days' % (s / 86400)
    if s >= 3600:
        return '%.2f h' % (s / 3600)
    if s >= 120:
        return '%.1f min' % (s / 60)
    return '%.1f s' % s


def format_report(rep):
    L = ['%s%s' % (rep['title'], (' [%s]' % rep['variant']) if rep['variant'] else '')]
    for r in rep['rows']:
        L.append('  %-9s %s Toffoli (%s) @ %d logical qubits   [%s; %s]' % (
            r['source'], _fmt(r['toffoli']), r['model'], r['qubits'], r['status'], r['ref']))
        if r['source'] == 'qoptct' and r['extra'].get('gain') is not None:
            L.append('            change against the published number: %+.1f %%' % (-100 * r['extra']['gain']))
        p = r.get('physical')
        if p:
            L.append('            physical (%s): %.3g qubits, d=%d, expected time %s, %.3g megaqubit-days'
                     % (p['model'], p['physical_qubits'], p['code_distance'], _fmt_time(p['expected_runtime_s']), p['megaqubit_days']))
        if r['note']:
            L.append('            note: ' + r['note'])
    return '\n'.join(L)


def attack_table(physical=False):
    """Markdown table: published / with qoptct blocks for every registered attack (default variants)."""
    L = ['| attack | published (Toffoli @ qubits) | with qoptct blocks | change | status |', '|---|---|---|---|---|']
    for name in ATTACKS:
        rep = attack_report(name, physical=physical)
        p, q = rep['rows']
        L.append('| %s | %s @ %d | %s @ %d | %+.1f %% | %s |' % (
            name, _fmt(p['toffoli']), p['qubits'], _fmt(q['toffoli']), q['qubits'],
            100 * (q['toffoli'] / p['toffoli'] - 1), q['status']))
    return '\n'.join(L)
