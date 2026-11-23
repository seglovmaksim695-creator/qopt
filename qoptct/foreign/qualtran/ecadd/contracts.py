"""Contracts of the opaque boxes used by the point-addition proofs, as DATA with the result that proves each of them.

The step provers treat Montgomery multiplication and Kaliski inversion as black boxes with a stated contract; the
glue (Z3) uses the arithmetic meaning.  A contract is only as good as the result that establishes it, so each
entry names that result (``proved_by``, a key of the qoptct results registry, assembled by the release step) and
lists the (n, p) instances for which that result exists.  For any other (n, p) the proof of the point addition is
CONDITIONAL on the contract: ``contracts_status`` in every verdict says so, and the verdict never hides it.

Edit these tables only together with the result that backs the change.
"""

P255 = 2 ** 255 - 19
SECP256K1 = 2 ** 256 - 2 ** 32 - 977
P256 = 2 ** 256 - 2 ** 224 + 2 ** 192 + 2 ** 96 - 1

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
        forward_instances=[(16, 32749), (256, P255)],
        pair_instances=[(16, 32749), (64, 2 ** 63 - 25), (128, 2 ** 127 - 1), (256, P255)],
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
        forward_instances=[(256, P255)],             # goal D also did n=64, 128 (other primes); only what ECAdd uses here
        pair_instances=[(16, 65521), (32, 1000000009), (64, 18446744073709551557), (128, 2 ** 127 - 1), (256, P255),
                        (256, SECP256K1), (256, P256)],
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
        forward_instances=[(16, 32749), (256, P255)],
        pair_instances=[],
    ),
}


def instance_covered(box, n, p, which='forward'):
    """True when the result named in ``CONTRACTS[box]`` exists for exactly this (n, p)."""
    c = CONTRACTS[box]
    inst = c['forward_instances' if which == 'forward' else 'pair_instances']
    for n0, p0 in inst:
        if n0 == n and (p0 is None or p0 == p):
            return True
    return False


def contracts_status(n, p, boxes=('MM', 'KINV')):
    """Per box: whether its contract and its pair identity are PROVED for this (n, p) or only ASSUMED here."""
    out = {}
    for b in boxes:
        c = CONTRACTS[b]
        out[b] = dict(
            box=c['box'], statement=c['statement'],
            forward=dict(proved_by=c['forward_proved_by'],
                         status='PROVED_FOR_THIS_INSTANCE' if instance_covered(b, n, p) else 'ASSUMED_FOR_THIS_INSTANCE'),
            pair=(None if c['pair_proved_by'] is None else
                  dict(proved_by=c['pair_proved_by'],
                       status='PROVED_FOR_THIS_INSTANCE' if instance_covered(b, n, p, 'pair') else 'ASSUMED_FOR_THIS_INSTANCE')))
    return out
