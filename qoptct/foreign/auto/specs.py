# SPDX-License-Identifier: Apache-2.0
"""Specifications of composite blocks: what the block computes, as formulas (data), nothing about how to prove it.

Each function returns a ``Contract`` (kind ``spec``) for given parameters and checks the parameter domain first
(``SpecDomainError``: the verdict is NOT_PROVED before any work).  A specification is used only as the statement of
a proof; it is itself checked against the classical simulation of the library (``registry.selfcheck``) on small
sizes.  Register names and packing follow the source signature (``qualtran_adapter``): a scalar register is an
unsigned integer, a one-dimensional bit array is an integer whose element 0 is the most significant bit.

Loops: a block that repeats one child needs a ``Fold`` -- the contract every iteration satisfies, the invariant and the
values the iterations receive.  These are formulas too; without them the mechanism refuses the loop (UNKNOWN).
"""
from __future__ import annotations

from . import expr as X
from .contract import Contract, Fold

__all__ = ['SpecDomainError', 'montgomery_window_body', 'montgomery_spec', 'mod_add_spec', 'cmod_add_spec',
           'mod_sub_spec', 'cmod_sub_spec', 'mod_neg_spec', 'cmod_neg_spec', 'mod_dbl_spec', 'is_probable_odd_prime']

V = X.var


class SpecDomainError(ValueError):
    """The parameters lie outside the domain the specification (and the proof) is about."""


def is_probable_odd_prime(p):
    """deterministic Miller-Rabin for p < 3.3e24, else 40 fixed bases (used only for messages; the composition proofs
    of this package do not depend on primality)"""
    if p < 3 or p % 2 == 0:
        return False
    for q in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if p % q == 0:
            return p == q
    d, s = p - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97):
        x = pow(a, d, p)
        if x in (1, p - 1):
            continue
        for _ in range(s - 1):
            x = x * x % p
            if x == p - 1:
                break
        else:
            return False
    return True


def _need(cond, msg):
    if not cond:
        raise SpecDomainError(msg)


def _odd_modulus(n, p, *, two_p=False):
    _need(type(n) is int and n >= 2, 'n must be an integer >= 2')
    _need(type(p) is int and p >= 3 and p % 2 == 1, 'p must be an odd integer >= 3')
    _need(p < (1 << n), 'p must be < 2^n')
    if two_p:
        _need(2 * p <= (1 << n), 'this block needs 2p <= 2^n')


# ------------------------------------------------------------------------------------------------ Montgomery
def montgomery_table(w, p):
    """the table of SingleWindowModMul: c(m) * p with c(m) = -m * p^-1 mod 2^w (as in Qualtran)"""
    inv = pow(p, 2 ** (w - 1) - 1, 2 ** w)
    N = 2 ** w
    return tuple(((-m * inv) % N) * p for m in range(N))


def montgomery_window_body(n, w, p):
    """Contract of one ``SingleWindowModMul(w, n, p)`` inside the multiplier, for the fold: on the domain
    T < y + p, y < p, qrom_index = 0:  S = T + sum_k x[k] * y * 2^k,  m = S mod 2^w,
    target' = (S + c(m) p) / 2^w,  qrom_index' = m;  x and y unchanged.
    (x is the bit array of the window, element k has weight 2^k; packed with element 0 as the most significant bit.)"""
    table = montgomery_table(w, p)
    terms = [X.ite(X.bit(V('x'), w - 1 - k), X.shl(V('y'), k) if k else V('y'), 0) for k in range(w)]
    S = X.add(V('target'), X.sum_(terms))
    defs = [('S', S), ('m', X.mod(X.ref('S'), 1 << w)),
            ('Tn', X.div(X.add(X.ref('S'), X.lut(X.ref('m'), table)), 1 << w))]
    return Contract(
        inputs={'x': w, 'y': n, 'target': n + w, 'qrom_index': w},
        outputs={'x': w, 'y': n, 'target': n + w, 'qrom_index': w},
        defs=defs,
        pre=X.all_([X.lt(V('target'), X.add(V('y'), p)), X.lt(V('y'), p), X.eq(V('qrom_index'), 0)]),
        out={'x': V('x'), 'y': V('y'), 'target': X.ref('Tn'), 'qrom_index': X.ref('m')},
        name=f'SingleWindowModMul body (n={n}, w={w}, p={p})', kind='body')


def montgomery_spec(n, w, p):
    """``DirtyOutOfPlaceMontgomeryModMul(n, w, p)`` (fixed Qualtran, PR 1997): for every n-bit x and every y < p,
    with the Montgomery recurrence T_0 = 0, T_{i+1} = window(T_i, x_i, y) over the n/w windows of x (window 0 = the
    lowest bits): target = T_N - p [T_N >= p], reduced = [T_N >= p], qrom_indices = the digits m_i (window 0 most
    significant), x and y unchanged.  The same statement as the 0.7 result ``montgomery_mul``.
    Fold 'win': the windows; invariant T < y + p (and y < p)."""
    _odd_modulus(n, p, two_p=True)
    _need(type(w) is int and w >= 1 and n % w == 0, 'the window size must divide n')
    N = n // w
    body = montgomery_window_body(n, w, p)
    inputs = []
    for i in range(N):
        xi = X.sum_([X.shl(X.bit(V('x'), w * i + k), w - 1 - k) if w - 1 - k else X.bit(V('x'), w * i + k)
                     for k in range(w)])
        inputs.append({'x': xi, 'y': V('y'), 'qrom_index': 0})
    fold = Fold(name='win', body='SingleWindowModMul', count=N, body_contract=body, state=['target'],
                init={'target': 0}, inputs=inputs, emit=['qrom_index'],
                invariant=X.and_(X.lt(V('target'), X.add(V('y'), p)), X.lt(V('y'), p)))
    TN = X.fold('win', 'target', N)
    red = X.not_(X.lt(TN, p))
    qrom = X.sum_([X.shl(X.fold('win', 'emit:qrom_index', i), w * (N - 1 - i)) if N - 1 - i else
                   X.fold('win', 'emit:qrom_index', i) for i in reversed(range(N))])
    return Contract(
        inputs={'x': n, 'y': n},
        outputs={'x': n, 'y': n, 'target': n, 'qrom_indices': n, 'reduced': 1},
        pre=X.lt(V('y'), p),
        defs=[('TN', TN), ('red', red)],
        out={'x': V('x'), 'y': V('y'), 'target': X.sub(X.ref('TN'), X.mul(X.ref('red'), p)),
             'qrom_indices': qrom, 'reduced': X.ref('red')},
        folds=[fold], name=f'DirtyOutOfPlaceMontgomeryModMul (n={n}, w={w}, p={p})', kind='spec',
        note='Montgomery recurrence with the final subtraction of p (PR 1997); invariant T < y + p given as data')


# ------------------------------------------------------------------------------------------------ modular add/sub
def _xy_lt_p(p):
    return X.and_(X.lt(V('x'), p), X.lt(V('y'), p))


def _x_le_p_y_lt_p(p):
    # Qualtran: "the construction still works when at most one of inputs equals mod"; the circuit does so for x = p
    # (y < p) -- for y = p it leaves its ancilla dirty, so y = p is NOT part of this domain
    return X.and_(X.le(V('x'), p), X.lt(V('y'), p))


def mod_add_spec(n, p):
    """ModAdd(n, p): y' = (x + y) mod p, x unchanged, for x <= p, y < p"""
    _odd_modulus(n, p)
    return Contract(inputs={'x': n, 'y': n}, outputs={'x': n, 'y': n}, pre=_x_le_p_y_lt_p(p),
                    out={'x': V('x'), 'y': X.mod(X.add(V('x'), V('y')), p)}, name=f'ModAdd(n={n}, p={p})', kind='spec')


def cmod_add_spec(n, p, cv=1):
    """CModAdd(QUInt(n), p, cv): y' = (x + y) mod p when ctrl == cv, else unchanged; x, ctrl unchanged;
    x <= p, y < p"""
    _odd_modulus(n, p)
    _need(cv in (0, 1), 'cv must be 0 or 1')
    act = X.eq(V('ctrl'), cv)
    return Contract(inputs={'ctrl': 1, 'x': n, 'y': n}, outputs={'ctrl': 1, 'x': n, 'y': n}, pre=_x_le_p_y_lt_p(p),
                    out={'ctrl': V('ctrl'), 'x': V('x'), 'y': X.ite(act, X.mod(X.add(V('x'), V('y')), p), V('y'))},
                    name=f'CModAdd(n={n}, p={p}, cv={cv})', kind='spec')


def mod_sub_spec(n, p):
    """ModSub(QUInt(n), p): y' = (y - x) mod p, x unchanged, for x, y < p"""
    _odd_modulus(n, p)
    return Contract(inputs={'x': n, 'y': n}, outputs={'x': n, 'y': n}, pre=_xy_lt_p(p),
                    out={'x': V('x'), 'y': X.mod(X.sub(V('y'), V('x')), p)}, name=f'ModSub(n={n}, p={p})', kind='spec')


def cmod_sub_spec(n, p, cv=1):
    """CModSub(QUInt(n), p, cv): y' = (y - x) mod p when ctrl == cv; x, y < p"""
    _odd_modulus(n, p)
    act = X.eq(V('ctrl'), cv)
    return Contract(inputs={'ctrl': 1, 'x': n, 'y': n}, outputs={'ctrl': 1, 'x': n, 'y': n}, pre=_xy_lt_p(p),
                    out={'ctrl': V('ctrl'), 'x': V('x'), 'y': X.ite(act, X.mod(X.sub(V('y'), V('x')), p), V('y'))},
                    name=f'CModSub(n={n}, p={p}, cv={cv})', kind='spec')


def mod_neg_spec(n, p):
    """ModNeg(QUInt(n), p): x' = (-x) mod p, for x < p"""
    _odd_modulus(n, p)
    return Contract(inputs={'x': n}, outputs={'x': n}, pre=X.lt(V('x'), p),
                    out={'x': X.mod(X.sub(0, V('x')), p)}, name=f'ModNeg(n={n}, p={p})', kind='spec')


def cmod_neg_spec(n, p, cv=1):
    """CModNeg(QUInt(n), p, cv): x' = (-x) mod p when ctrl == cv; x < p"""
    _odd_modulus(n, p)
    act = X.eq(V('ctrl'), cv)
    return Contract(inputs={'ctrl': 1, 'x': n}, outputs={'ctrl': 1, 'x': n}, pre=X.lt(V('x'), p),
                    out={'ctrl': V('ctrl'), 'x': X.ite(act, X.mod(X.sub(0, V('x')), p), V('x'))},
                    name=f'CModNeg(n={n}, p={p}, cv={cv})', kind='spec')


def mod_dbl_spec(n, p):
    """ModDbl(QUInt(n), p): x' = 2x mod p, for x < p"""
    _odd_modulus(n, p)
    return Contract(inputs={'x': n}, outputs={'x': n}, pre=X.lt(V('x'), p),
                    out={'x': X.mod(X.shl(V('x'), 1), p)}, name=f'ModDbl(n={n}, p={p})', kind='spec')


# ------------------------------------------------------------------------------------------------ more blocks
def _width(n):
    _need(type(n) is int and n >= 1, 'the width must be a positive integer')


def subtract_spec(n):
    """Subtract(QUInt(n), QUInt(n)): b' = (a - b) mod 2^n, a unchanged (all inputs)"""
    _width(n)
    return Contract(inputs={'a': n, 'b': n}, outputs={'a': n, 'b': n},
                    out={'a': V('a'), 'b': X.mod(X.sub(V('a'), V('b')), 1 << n)}, name=f'Subtract(n={n})', kind='spec')


def subtract_from_spec(n):
    """SubtractFrom(QUInt(n)): b' = (b - a) mod 2^n, a unchanged (all inputs)"""
    _width(n)
    return Contract(inputs={'a': n, 'b': n}, outputs={'a': n, 'b': n},
                    out={'a': V('a'), 'b': X.mod(X.sub(V('b'), V('a')), 1 << n)}, name=f'SubtractFrom(n={n})',
                    kind='spec')


def negate_spec(n):
    """Negate(QUInt(n)): x' = (-x) mod 2^n (all inputs)"""
    _width(n)
    return Contract(inputs={'x': n}, outputs={'x': n}, out={'x': X.mod(X.sub(0, V('x')), 1 << n)},
                    name=f'Negate(n={n})', kind='spec')


def cmod_add_k_spec(n, k, p):
    """CModAddK(k, p, n): x' = (x + k) mod p if ctrl = 1, else x; ctrl unchanged; for x < p (0 <= k < p)"""
    _odd_modulus(n, p)
    _need(type(k) is int and 0 <= k < p, 'k must satisfy 0 <= k < p')
    act = X.eq(V('ctrl'), 1)
    return Contract(inputs={'ctrl': 1, 'x': n}, outputs={'ctrl': 1, 'x': n}, pre=X.lt(V('x'), p),
                    out={'ctrl': V('ctrl'), 'x': X.ite(act, X.mod(X.add(V('x'), k), p), V('x'))},
                    name=f'CModAddK(n={n}, k={k}, p={p})', kind='spec')


def ctrl_scale_mod_add_spec(n, k, p):
    """CtrlScaleModAdd(k, p, n): y' = (y + x k) mod p if ctrl = 1, else y; ctrl, x unchanged; for y < p (x any;
    k any integer, written as k mod p)"""
    _odd_modulus(n, p)
    _need(type(k) is int, 'k must be an integer')
    km = k % p
    act = X.eq(V('ctrl'), 1)
    return Contract(inputs={'ctrl': 1, 'x': n, 'y': n}, outputs={'ctrl': 1, 'x': n, 'y': n}, pre=X.lt(V('y'), p),
                    out={'ctrl': V('ctrl'), 'x': V('x'),
                         'y': X.ite(act, X.mod(X.add(V('y'), X.mul(V('x'), km)), p), V('y'))},
                    name=f'CtrlScaleModAdd(n={n}, k={k}, p={p})', kind='spec')


def cmod_mul_k_spec(n, k, p):
    """CModMulK(QUInt(n), k, p): x' = x k mod p if ctrl = 1, else x; for x < p (p prime, 0 < k < p)"""
    _odd_modulus(n, p)
    _need(type(k) is int and 0 < k < p, 'k must satisfy 0 < k < p')
    _need(is_probable_odd_prime(p), 'p must be prime (k must be invertible)')
    act = X.eq(V('ctrl'), 1)
    return Contract(inputs={'ctrl': 1, 'x': n}, outputs={'ctrl': 1, 'x': n}, pre=X.lt(V('x'), p),
                    out={'ctrl': V('ctrl'), 'x': X.ite(act, X.mod(X.mul(V('x'), k), p), V('x'))},
                    name=f'CModMulK(n={n}, k={k}, p={p})', kind='spec')


# ------------------------------------------------------------------------------------------------ ECAdd, step 1
def ecadd_step1_spec(n, p):
    """``_ECAddStepOne(n, p)`` of the fixed ECAdd (patch E05, Qualtran kind 'e'): for a, b, x, y < p it sets the
    flags f1 = [a = x], f3 = [(a, b) = (0, 0)], f4 = [(x, y) = (0, 0)], ctrl = [f3 = f4 = 0 and not (f1 and
    b = -y)], f2 = [b = -y] and not ctrl; a, b, x, y unchanged (the same statement as step 1 of the 0.7 proof,
    ``ecadd.ispecs.S1``)."""
    _odd_modulus(n, p, two_p=True)
    e1 = X.eq(V('a'), V('x'))
    e2 = X.eq(V('b'), X.mod(X.sub(0, V('y')), p))
    e3 = X.and_(X.eq(V('a'), 0), X.eq(V('b'), 0))
    e4 = X.and_(X.eq(V('x'), 0), X.eq(V('y'), 0))
    ctrl = X.and_(X.not_(X.ref('e3')), X.and_(X.not_(X.ref('e4')), X.not_(X.and_(X.ref('e1'), X.ref('e2')))))
    return Contract(
        inputs={'a': n, 'b': n, 'x': n, 'y': n},
        outputs={'f1': 1, 'f2': 1, 'f3': 1, 'f4': 1, 'ctrl': 1, 'a': n, 'b': n, 'x': n, 'y': n},
        pre=X.all_([X.lt(V(k), p) for k in 'abxy']),
        defs=[('e1', e1), ('e2', e2), ('e3', e3), ('e4', e4), ('ctrl', ctrl)],
        out={'f1': X.ref('e1'), 'f2': X.and_(X.ref('e2'), X.not_(X.ref('ctrl'))), 'f3': X.ref('e3'),
             'f4': X.ref('e4'), 'ctrl': X.ref('ctrl'), 'a': V('a'), 'b': V('b'), 'x': V('x'), 'y': V('y')},
        name=f'ECAdd step 1 (n={n}, p={p}, fixed)', kind='spec')
