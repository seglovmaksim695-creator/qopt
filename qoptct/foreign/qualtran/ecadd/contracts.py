"""Contracts of the opaque boxes used by the point-addition proofs, as DATA with the result that proves each of them.

The step provers treat Montgomery multiplication and Kaliski inversion as black boxes with a stated contract; the
glue (Z3) uses the arithmetic meaning.  A contract is only as good as the result that establishes it, so each
entry names that result (``forward_proved_by`` / ``pair_proved_by``: the name in the registry of ``qoptct.foreign.qualtran.registry``)
and maps every (n, p) instance for which that result exists to its registry id (``forward_evidence`` / ``pair_evidence``;
``tests/foreign_qualtran/test_registry.py`` checks that every id exists and is PROVED).  For any other (n, p) the proof of the point addition is
CONDITIONAL on the contract: ``contracts_status`` in every verdict says so, and the verdict never hides it.

Edit these tables only together with the result that backs the change.
"""

P255 = 2 ** 255 - 19
SECP256K1 = 2 ** 256 - 2 ** 32 - 977
P256 = 2 ** 256 - 2 ** 224 + 2 ** 192 + 2 ** 96 - 1
P64M = 2 ** 63 - 25            # the prime of the Montgomery runs at n = 64
P127 = 2 ** 127 - 1

CONTRACTS = {
    'MM': dict(
        box='DirtyOutOfPlaceMontgomeryModMul(bitsize=n, window_size=1, mod=p), forward and uncompute=True',
        statement=('x, y < p  =>  target < p and target = x*y*2^-n (mod p); the garbage registers (qrom_indices, reduced) '
                   'are arbitrary; the uncompute box returns exactly the registers of its forward box to zero '
                   '(pair identity)'),
        used_for='values of the multiplier in the point-addition glue (MM0, C1, C2) and the box pairing in steps 3-5',
        qualtran_kind='c',
        forward_proved_by='montgomery_mul',          # goal C: contract of the multiplier (task I3)
        pair_proved_by='pair_montgomery',            # goal H2: forward + inverse = identity (pairs.check_pair_montgomery)
        forward_evidence={(16, 32749): 'montgomery_mul/montgomery_n16_w1', (256, P255): 'montgomery_mul/montgomery_n256_w1'},
        pair_evidence={(16, 32749): 'pair_montgomery/pairs_montgomery_n16', (64, P64M): 'pair_montgomery/pairs_montgomery_n64',
                       (128, P127): 'pair_montgomery/pairs_montgomery_n128', (256, P255): 'pair_montgomery/pairs_montgomery_n256'},
    ),
    'KINV': dict(
        box='KaliskiModInverse(bitsize=n, mod=p), forward and uncompute=True',
        statement=('0 <= x < p  =>  x_out < p and x_out = x^-1 * 2^(2n) (mod p) for x != 0, x_out = 0 for x = 0; the junk '
                   'register is arbitrary; the uncompute box returns exactly the registers of its forward box to zero '
                   '(pair identity)'),
        used_for='values of the inverse in the point-addition glue (K0, C1, C2) and the box pairing in steps 2 and 5',
        qualtran_kind=None,
        forward_proved_by='kaliski_inverse',         # goal D: contract of the Kaliski inversion (task I2)
        pair_proved_by='pair_kaliski',               # goal H1: forward + inverse = identity (pairs.check_pair_kaliski)
        forward_evidence={(16, 32749): 'kaliski_inverse/kaliski_n16_p32749',
                          (256, P255): 'kaliski_inverse/kaliski_n256_p57896044618658097711785492504343953926634992332820282019728792003956564819949'},
        pair_evidence={(16, 65521): 'pair_kaliski/pairs_kaliski_n16', (16, 32749): 'pair_kaliski/pairs_kaliski_n16_p32749',
                       (32, 1000000009): 'pair_kaliski/pairs_kaliski_n32', (64, 18446744073709551557): 'pair_kaliski/pairs_kaliski_n64',
                       (128, P127): 'pair_kaliski/pairs_kaliski_n128', (256, P255): 'pair_kaliski/pairs_kaliski_n256',
                       (256, SECP256K1): 'pair_kaliski/pairs_kaliski_n256_secp256k1', (256, P256): 'pair_kaliski/pairs_kaliski_n256_p256'},
    ),
    'ECADD': dict(
        box='ECAdd(n, p, window_size=1), fixed version (flag patch E05)',
        statement=('(a, b, x, y, lam_r) with a, b, x, y, lam_r < p  =>  (x, y) becomes the Montgomery-form point sum of '
                   '(a, b) and (x, y) (P + Q, including P = O, Q = O, P = Q, P = -Q); a, b, lam_r unchanged; every '
                   'other wire clean'),
        used_for='ECWindowAddR (the ECAdd box inside the table lookup)',
        qualtran_kind='e',
        forward_proved_by='ecadd',                   # this package: prove_ecadd(n, p, patched=True)
        pair_proved_by=None,
        forward_evidence={(16, 32749): 'ecadd/ecadd_n16', (256, P255): 'ecadd/ecadd_n256'},
        pair_evidence={},
    ),
}


for _c in CONTRACTS.values():            # the plain instance lists (kept for callers that only need the (n, p) pairs)
    _c['forward_instances'] = list(_c['forward_evidence'])
    _c['pair_instances'] = list(_c['pair_evidence'])


def evidence_id(box, n, p, which='forward'):
    """Registry id of the result that covers this exact (n, p), or ``None``."""
    return CONTRACTS[box]['forward_evidence' if which == 'forward' else 'pair_evidence'].get((n, p))


def instance_covered(box, n, p, which='forward'):
    """True when the result named in ``CONTRACTS[box]`` exists for exactly this (n, p)."""
    return evidence_id(box, n, p, which) is not None


def contracts_status(n, p, boxes=('MM', 'KINV')):
    """Per box: whether its contract and its pair identity are PROVED for this (n, p) or only ASSUMED here."""
    out = {}
    for b in boxes:
        c = CONTRACTS[b]
        fe, pe = evidence_id(b, n, p), evidence_id(b, n, p, 'pair')
        out[b] = dict(
            box=c['box'], statement=c['statement'],
            forward=dict(proved_by=c['forward_proved_by'], evidence=fe,
                         status='PROVED_FOR_THIS_INSTANCE' if fe else 'ASSUMED_FOR_THIS_INSTANCE'),
            pair=(None if c['pair_proved_by'] is None else
                  dict(proved_by=c['pair_proved_by'], evidence=pe,
                       status='PROVED_FOR_THIS_INSTANCE' if pe else 'ASSUMED_FOR_THIS_INSTANCE')))
    return out
