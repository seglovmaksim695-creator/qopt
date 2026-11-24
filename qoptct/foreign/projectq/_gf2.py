"""GF(2^m) multiplier check: ``result = x * y mod f`` proved for ALL inputs with qoptct.verify.exact2.

The field arithmetic here is independent of Qualtran, ProjectQ and qoptct (bit i of an integer = coefficient of x^i).
The polynomial ``f`` is an explicit argument; irreducibility is decided by Rabin's exact test for any degree.

Contract proved (status PROVED):  x' = x,  y' = y,  result' = x*y mod f,  every wire that is not x, y, result and not a
PINNED garbage wire ends in |0>, and every garbage wire ends in exactly the value of the pin.  A pin
(:class:`GarbagePin`) is taken from a golden program (the author's circuit); a candidate that differs on a pinned wire is
REFUTED, so a mutant can not hide in the garbage (in the first campaign of this bridge the garbage of the circuit under
test was its own reference, and one mutant was PROVED).  The values of the garbage wires are NOT proved to be
meaningful: they are recorded and compared, nothing more.
"""
from __future__ import annotations

import hashlib
import random

from ._record import Unsupported
from ._sim import PreconditionError, masks_w2, masks_full, simulate, w2_monomials

__all__ = ['GarbagePin', 'rabin_irreducible', 'poly_mod', 'mulmod', 'pin_gf2', 'check_gf2_mul', 'reference_text']


def poly_mod(value, poly):
    degree = poly.bit_length() - 1
    r = value
    for top in range(value.bit_length() - 1, degree - 1, -1):
        if (r >> top) & 1:
            r ^= poly << (top - degree)
    return r


def _pmul(a, b):
    p = 0
    for d in range(b.bit_length()):
        if (b >> d) & 1:
            p ^= a << d
    return p


def mulmod(a, b, poly):
    return poly_mod(_pmul(a, b), poly)


def _gcd(a, b):
    while b:
        if b == 1:
            return 1
        a, b = b, poly_mod(a, b)
    return a


def _primes(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def rabin_irreducible(poly):
    """Exact Rabin test over GF(2): x^(2^n) = x mod f and gcd(f, x^(2^(n/r)) - x) = 1 for every prime r | n."""
    n = poly.bit_length() - 1
    if n < 1 or not poly & 1:
        return False
    x = poly_mod(2, poly)
    at = {n // r for r in _primes(n)}
    fr = x
    for e in range(1, n + 1):
        sq = 0
        for i in range(fr.bit_length()):
            if (fr >> i) & 1:
                sq |= 1 << (2 * i)
        fr = poly_mod(sq, poly)
        if e in at and _gcd(poly, fr ^ x) != 1:
            return False
    return fr == x


class GarbagePin:
    """Final values of the garbage wires of a golden program: wire -> lane integer on the lane set ``kind``.

    ``kind='w2'`` (GF2 multipliers): weight <= 2 lanes of the free bits ordered x then y.  ``kind='full'``: all 2^n inputs of
    the free bits in the order of ``Logical.free_wires``.  ``digest`` identifies the pin."""

    def __init__(self, kind, nbits, wires):
        self.kind, self.nbits, self.wires = kind, nbits, dict(wires)
        h = hashlib.sha256(f'{kind}:{nbits}'.encode())
        for w in sorted(self.wires):
            v = self.wires[w]
            h.update(f'{w}:{v.bit_length()}:'.encode() + v.to_bytes((v.bit_length() + 7) // 8, 'little'))
        self.digest = h.hexdigest()

    def __len__(self):
        return len(self.wires)


def _layout(lg, m):
    ins, outs = lg.inputs, lg.outputs
    if set(ins) != {'x', 'y'} or len(ins['x']) != m or len(ins['y']) != m:
        return 'the program must have exactly the input registers x and y of m bits'
    if 'result' not in outs or len(outs['result']) != m or outs.get('x') != ins['x'] or outs.get('y') != ins['y']:
        return 'the outputs must be x, y (the same wires as the inputs) and result of m bits'
    if set(outs['result']) & (set(ins['x']) | set(ins['y'])):
        return 'result must be made of fresh wires'
    if lg.garbage:
        return 'declared garbage registers are not used by the GF2 multiplier check (use a pin)'
    return None


def _free_order(lg):
    return list(lg.inputs['x']) + list(lg.inputs['y'])


def pin_gf2(lg, m):
    """Pin of the golden program ``lg`` (all wires that are not x, y, result and end non-zero)."""
    err = _layout(lg, m)
    if err:
        raise ValueError(err)
    free = _free_order(lg)
    masks, full, _ = masks_w2(len(free))
    st = simulate(lg.ops, lg.nq, dict(zip(free, masks)), full)
    keep = set(free) | set(lg.outputs['result'])
    return GarbagePin('w2', len(free), {w: st[w] for w in range(lg.nq) if w not in keep and st[w]})


def reference_text(lg, m, poly, pin=None):
    """The schoolbook reference as .qc text (tof only), plus the pinned garbage wires as their own ANF."""
    free = _free_order(lg)
    xs, ys, res = lg.inputs['x'], lg.inputs['y'], lg.outputs['result']
    lines = ['.v ' + ' '.join(f'q{i}' for i in range(lg.nq)), '.i ' + ' '.join(f'q{i}' for i in sorted(free)), 'BEGIN']
    resid = [poly_mod(1 << d, poly) for d in range(2 * m - 1)]
    for i, xw in enumerate(xs):
        for j, yw in enumerate(ys):
            r = resid[i + j]
            for k, t in enumerate(res):
                if (r >> k) & 1:
                    lines.append(f'tof q{xw} q{yw} q{t}')
    if pin is not None:
        if pin.kind != 'w2' or pin.nbits != len(free):
            raise ValueError('the pin does not belong to this multiplier')
        for w in sorted(pin.wires):
            c0, linear, quad = w2_monomials(pin.wires[w], len(free))
            if c0:
                lines.append(f'tof q{w}')
            for u in linear:
                lines.append(f'tof q{free[u]} q{w}')
            for u, v in quad:
                lines.append(f'tof q{free[u]} q{free[v]} q{w}')
    lines.append('END')
    return '\n'.join(lines) + '\n'


def _samples(m, poly, count, seed):
    rng = random.Random(seed)
    top = (1 << m) - 1
    vals = [0, 1, 2, top, top ^ 1, 1 << (m - 1), (1 << (m - 1)) | 1, poly & top, 3]
    pairs = [(a, b) for a in vals for b in vals[:6]]
    while len(pairs) < count:
        pairs.append((rng.getrandbits(m), rng.getrandbits(m)))
    return pairs[:max(count, 1)]


def sample_check(lg, m, poly, pin=None, count=96, seed=1):
    """Independent run on ``count`` explicit input pairs (corner values and random ones) against the schoolbook product.
    Every wire that is neither x, y, result nor pinned must be 0 at the end.  -> (ok, detail)."""
    pairs = _samples(m, poly, count, seed)
    L = len(pairs)
    full = (1 << L) - 1
    init = {}
    for i, w in enumerate(lg.inputs['x']):
        init[w] = sum(((a >> i) & 1) << j for j, (a, _) in enumerate(pairs))
    for i, w in enumerate(lg.inputs['y']):
        init[w] = sum(((b >> i) & 1) << j for j, (_, b) in enumerate(pairs))
    try:
        st = simulate(lg.parsed_ops(), lg.nq, init, full)
    except Unsupported as e:
        return False, str(e)
    except PreconditionError as e:
        return False, f'precondition on the sample lanes: {e}'
    for j, (a, b) in enumerate(pairs):
        got = sum(((st[w] >> j) & 1) << i for i, w in enumerate(lg.outputs['result']))
        if got != mulmod(a, b, poly):
            return False, f'sample {j}: x={a:#x} y={b:#x} gave {got:#x}, expected {mulmod(a, b, poly):#x}'
        for name in ('x', 'y'):
            if any(((st[w] >> j) & 1) != ((a if name == 'x' else b) >> i) & 1 for i, w in enumerate(lg.inputs[name])):
                return False, f'sample {j}: register {name} was changed'
    keep = set(_free_order(lg)) | set(lg.outputs['result']) | (set(pin.wires) if pin else set())
    bad = [w for w in range(lg.nq) if w not in keep and st[w]]
    if bad:
        return False, f'{len(bad)} wires that must be clean are not (first: {bad[0]})'
    return True, f'{L} sample lanes agree with the schoolbook product, x and y unchanged, clean wires clean'


def check_gf2_mul(lg, m, poly, pin=None):
    """-> dict(status, detail, ...) for the multiplier ``lg`` (a classical :class:`Logical`).  Status PROVED needs the exact2
    proof AND the independent sample check; a disagreement is UNKNOWN."""
    out = dict(claim=f'x*y mod f in GF(2^{m}), f={hex(poly)}', m=m, poly=poly)
    err = _layout(lg, m)
    if err:
        return dict(out, status='UNKNOWN', detail=err)
    if poly.bit_length() - 1 != m or not rabin_irreducible(poly):
        return dict(out, status='UNKNOWN', detail='the polynomial is not irreducible of degree m (Rabin test)')
    if not lg.is_classical():
        return dict(out, status='UNKNOWN', detail='the program is not an x/cx/and/erase program')
    from ...verify.exact2 import prove_degree2
    qasm = lg.qasm()
    ref = reference_text(lg, m, poly, pin)
    try:
        r = prove_degree2(qasm, ref, clean=True)
    except Exception as e:  # noqa: BLE001 - never an exception to the caller
        return dict(out, status='UNKNOWN', detail=f'exact2 failed: {type(e).__name__}: {e}')
    st = r.get('status')
    ok, why = sample_check(lg, m, poly, pin)
    out.update(method='exact2', exact2=dict(status=st, detail=str(r.get('detail'))[:300], time=r.get('time'), lanes=r.get('lanes')),
               independent=dict(ok=ok, detail=why), garbage_wires=0 if pin is None else len(pin), pin=None if pin is None else pin.digest)
    if st == 'PROVED' and ok:
        return dict(out, status='PROVED', detail=str(r.get('detail'))[:200])
    if st == 'PROVED':
        return dict(out, status='UNKNOWN', detail='exact2 proved but the independent sample check disagrees: ' + why)
    if st == 'REFUTED':
        return dict(out, status='REFUTED', detail=str(r.get('detail'))[:300], witness=r.get('counterexample'))
    return dict(out, status='UNKNOWN', detail=str(r.get('detail'))[:300])
