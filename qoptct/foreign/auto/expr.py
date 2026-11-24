# SPDX-License-Identifier: Apache-2.0
"""Contract expressions of ``qoptct.foreign.auto``: formulas as data.

An expression is a Python ``int`` (a constant) or a tuple ``(op, *args)``:

==========================  ===================================================================================
``('var', name)``           value of the input register ``name`` (unsigned, ``0 <= v < 2**width``)
``('ref', name)``           a let-bound intermediate (``Contract.defs``)
``('val', vid)``            composition only: a value of the symbolic state
``('sym', name)``           composition only: a fold symbol (the carried state of a loop)
``('fold', f, key, i)``     specification only: state ``key`` of fold ``f`` before iteration ``i`` (``i = count``:
                            the final state), or ``'emit:<reg>'`` = what iteration ``i`` writes to ``reg``
``('add'|'sub', a, b)``     integer sum / difference (exact)
``('mul', a, b)``           product; one factor must be a constant or 0/1-valued (linear arithmetic only)
``('div'|'mod', a, k)``     floor division / remainder by a positive constant ``k``
``('lt'|'le'|'eq'|'ne', a, b)``   comparisons, value 0 or 1
``('and'|'or', a, b)``, ``('not', a)``   on 0/1 values
``('ite', c, a, b)``        ``a`` if ``c == 1`` else ``b`` (``c`` is 0/1)
``('bxor'|'band'|'bor', a, b)``   bitwise, non-negative operands
``('bit', a, j)``           ``(a >> j) & 1``
``('slice', a, j, l)``      ``(a >> j) mod 2**l``
``('shl'|'shr', a, k)``     ``a * 2**k``, ``floor(a / 2**k)``
``('lut', i, (t0, ...))``   a table of integer constants; ``i`` must lie in ``range(len(table))``
==========================  ===================================================================================

Semantics are those of Python integers (exact, unbounded); booleans are 0 and 1.  Back ends: ``evaluate`` (Python
int), ``Z3Terms`` (z3 ``Int`` terms: linear integer arithmetic with div/mod by constants), ``bdd_word`` (two's
complement ROBDD words for the kprove specifications of the leaves) and ``bounds`` (static intervals).  The module
imports neither Qualtran nor Z3 at import time.
"""
from __future__ import annotations

__all__ = ['ExprError', 'ARITY', 'validate', 'evaluate', 'bounds', 'Z3Terms', 'bdd_word', 'BddWord', 'subst',
           'refs_of', 'var', 'ref', 'add', 'sub', 'mul', 'div', 'mod', 'lt', 'le', 'eq', 'ne', 'and_', 'or_',
           'not_', 'ite', 'bxor', 'band', 'bor', 'bit', 'slc', 'shl', 'shr', 'lut', 'sum_', 'all_', 'fold',
           'to_json', 'from_json', 'digest']

import hashlib
import json


class ExprError(ValueError):
    """An expression outside the language, or an operation the chosen back end cannot express exactly."""


ARITY = {'var': 1, 'ref': 1, 'val': 1, 'sym': 1, 'fold': 3, 'add': 2, 'sub': 2, 'mul': 2, 'div': 2, 'mod': 2,
         'lt': 2, 'le': 2, 'eq': 2, 'ne': 2, 'and': 2, 'or': 2, 'not': 1, 'ite': 3, 'bxor': 2, 'band': 2,
         'bor': 2, 'bit': 2, 'slice': 3, 'shl': 2, 'shr': 2, 'lut': 2}
_LEAF = ('var', 'ref', 'val', 'sym', 'fold')
_CONST_ARGS = {'div': (1,), 'mod': (1,), 'bit': (1,), 'slice': (1, 2), 'shl': (1,), 'shr': (1,)}


def _n_expr(op):
    return ARITY[op] - len(_CONST_ARGS.get(op, ()))


# ----------------------------------------------------------------------------------------------- constructors
def var(name): return ('var', name)
def ref(name): return ('ref', name)
def add(a, b): return ('add', a, b)
def sub(a, b): return ('sub', a, b)
def mul(a, b): return ('mul', a, b)
def div(a, k): return ('div', a, int(k))
def mod(a, k): return ('mod', a, int(k))
def lt(a, b): return ('lt', a, b)
def le(a, b): return ('le', a, b)
def eq(a, b): return ('eq', a, b)
def ne(a, b): return ('ne', a, b)
def and_(a, b): return ('and', a, b)
def or_(a, b): return ('or', a, b)
def not_(a): return ('not', a)
def ite(c, a, b): return ('ite', c, a, b)
def bxor(a, b): return ('bxor', a, b)
def band(a, b): return ('band', a, b)
def bor(a, b): return ('bor', a, b)
def bit(a, j): return ('bit', a, int(j))
def slc(a, j, l): return ('slice', a, int(j), int(l))
def shl(a, k): return ('shl', a, int(k))
def shr(a, k): return ('shr', a, int(k))
def lut(i, table): return ('lut', i, tuple(int(t) for t in table))
def fold(f, key, i): return ('fold', f, key, int(i))


def sum_(terms):
    """sum of a list of expressions (0 for an empty list), as a balanced tree"""
    terms = list(terms)
    if not terms:
        return 0
    while len(terms) > 1:
        terms = [add(terms[i], terms[i + 1]) if i + 1 < len(terms) else terms[i] for i in range(0, len(terms), 2)]
    return terms[0]


def all_(conds):
    """conjunction of 0/1 expressions (1 for an empty list)"""
    conds = list(conds)
    if not conds:
        return 1
    while len(conds) > 1:
        conds = [and_(conds[i], conds[i + 1]) if i + 1 < len(conds) else conds[i] for i in range(0, len(conds), 2)]
    return conds[0]


# ----------------------------------------------------------------------------------------------- structure
def validate(e, allowed=('var', 'ref')):
    """Check the structure; ``allowed``: the permitted leaf kinds."""
    seen = set()
    stack = [e]
    while stack:
        x = stack.pop()
        if type(x) is int:
            continue
        if not isinstance(x, tuple) or not x or x[0] not in ARITY:
            raise ExprError(f'not an expression: {x!r}')
        if id(x) in seen:
            continue
        seen.add(id(x))
        op = x[0]
        if len(x) != ARITY[op] + 1:
            raise ExprError(f'{op}: wrong arity')
        if op in _LEAF:
            if op not in allowed:
                raise ExprError(f'leaf {op} is not allowed here')
            if op == 'fold' and (not isinstance(x[1], str) or not isinstance(x[2], str) or type(x[3]) is not int):
                raise ExprError('fold leaf: (fold, name, key, i)')
            continue
        if op == 'lut':
            if not isinstance(x[2], tuple) or not x[2] or any(type(t) is not int for t in x[2]):
                raise ExprError('lut: the table must be a non-empty tuple of int')
            stack.append(x[1])
            continue
        for pos in _CONST_ARGS.get(op, ()):
            c = x[pos + 1]
            if type(c) is not int or c < 0 or (op in ('div', 'mod') and c <= 0):
                raise ExprError(f'{op}: argument {pos} must be a non-negative constant (positive divisor)')
        stack.extend(x[1:1 + _n_expr(op)])
    return e


def refs_of(e, kind):
    """set of the names of the leaves of one kind in e (for 'fold': tuples (f, key, i))"""
    out, stack, seen = set(), [e], set()
    while stack:
        x = stack.pop()
        if type(x) is int or id(x) in seen:
            continue
        seen.add(id(x))
        if x[0] == kind:
            out.add(x[1] if kind != 'fold' else x[1:])
        elif x[0] in _LEAF:
            continue
        elif x[0] == 'lut':
            stack.append(x[1])
        else:
            stack.extend(x[1:1 + _n_expr(x[0])])
    return out


def subst(e, fn, memo=None):
    """Replace leaves: fn(leaf_tuple) -> expression or None (keep).  Shared sub-expressions stay shared."""
    memo = {} if memo is None else memo

    def go(x):
        if type(x) is int:
            return x
        k = id(x)
        if k in memo:
            return memo[k][1]
        op = x[0]
        if op in _LEAF:
            r = fn(x)
            r = x if r is None else r
        elif op == 'lut':
            a = go(x[1])
            r = x if a is x[1] else ('lut', a, x[2])
        else:
            ne_ = _n_expr(op)
            args = tuple(go(a) for a in x[1:1 + ne_])
            r = x if all(a is b for a, b in zip(args, x[1:1 + ne_])) else (op,) + args + x[1 + ne_:]
        memo[k] = (x, r)          # keep x alive so that its id is not reused
        return r
    return go(e)


def to_json(e):
    """JSON-safe form (tuples -> lists)."""
    if type(e) is int:
        return e
    if e[0] == 'lut':
        return ['lut', to_json(e[1]), list(e[2])]
    if e[0] in _LEAF:
        return list(e)
    ne_ = _n_expr(e[0])
    return [e[0]] + [to_json(a) for a in e[1:1 + ne_]] + list(e[1 + ne_:])


def from_json(j):
    if type(j) is int:
        return j
    op = j[0]
    if op == 'lut':
        return ('lut', from_json(j[1]), tuple(j[2]))
    if op in _LEAF:
        return tuple(j)
    ne_ = _n_expr(op)
    return (op,) + tuple(from_json(a) for a in j[1:1 + ne_]) + tuple(j[1 + ne_:])


def digest(obj):
    """SHA-256 of a canonical JSON form."""
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


# ----------------------------------------------------------------------------------------------- Python
def evaluate(e, env, memo=None):
    """Exact value.  env(leaf_tuple) -> int for the leaves."""
    memo = {} if memo is None else memo

    def go(x):
        if type(x) is int:
            return x
        k = id(x)
        if k in memo:
            return memo[k][1]
        op = x[0]
        if op in _LEAF:
            v = env(x)
            if type(v) is not int:
                raise ExprError(f'leaf {x!r} has no integer value')
        elif op == 'lut':
            i = go(x[1])
            if not 0 <= i < len(x[2]):
                raise ExprError(f'lut index {i} outside the table')
            v = x[2][i]
        elif op == 'ite':
            c = go(x[1])
            if c not in (0, 1):
                raise ExprError('ite: condition is not 0/1')
            v = go(x[2]) if c else go(x[3])
        elif op == 'not':
            a = go(x[1])
            if a not in (0, 1):
                raise ExprError('not: operand must be 0/1')
            v = 1 - a
        elif op in _CONST_ARGS:
            a = go(x[1])
            c = x[2]
            if op == 'div': v = a // c
            elif op == 'mod': v = a % c
            elif op == 'bit': v = (a >> c) & 1
            elif op == 'shl': v = a << c
            elif op == 'shr': v = a >> c
            else: v = (a >> c) % (1 << x[3])
        else:
            a, b = go(x[1]), go(x[2])
            if op == 'add': v = a + b
            elif op == 'sub': v = a - b
            elif op == 'mul': v = a * b
            elif op == 'lt': v = int(a < b)
            elif op == 'le': v = int(a <= b)
            elif op == 'eq': v = int(a == b)
            elif op == 'ne': v = int(a != b)
            elif op in ('and', 'or'):
                if a not in (0, 1) or b not in (0, 1):
                    raise ExprError(f'{op}: operands must be 0/1')
                v = (a & b) if op == 'and' else (a | b)
            elif op in ('bxor', 'band', 'bor'):
                if a < 0 or b < 0:
                    raise ExprError(f'{op}: negative operand')
                v = (a ^ b) if op == 'bxor' else (a & b) if op == 'band' else (a | b)
            else:
                raise ExprError(f'unknown operation {op}')
        memo[k] = (x, v)
        return v
    return go(e)


# ----------------------------------------------------------------------------------------------- intervals
def bounds(e, leaf_bounds, memo=None):
    """Static interval (lo, hi) of e; leaf_bounds(leaf_tuple) -> (lo, hi).  Conservative, never wrong."""
    memo = {} if memo is None else memo

    def go(x):
        if type(x) is int:
            return (x, x)
        k = id(x)
        if k in memo:
            return memo[k][1]
        op = x[0]
        if op in _LEAF:
            r = leaf_bounds(x)
        elif op == 'lut':
            i = go(x[1])
            if i[0] < 0 or i[1] >= len(x[2]):
                raise ExprError(f'lut index range {i} is not inside the table of {len(x[2])} entries')
            vals = x[2][i[0]:i[1] + 1]
            r = (min(vals), max(vals))
        elif op == 'ite':
            c, a, b = go(x[1]), go(x[2]), go(x[3])
            r = a if c == (1, 1) else b if c == (0, 0) else (min(a[0], b[0]), max(a[1], b[1]))
        elif op == 'not':
            a = go(x[1])
            r = (1 - min(a[1], 1), 1 - max(a[0], 0))
        elif op in ('lt', 'le', 'eq', 'ne', 'and', 'or'):
            go(x[1]); go(x[2])
            r = (0, 1)
        elif op == 'bit':
            go(x[1])
            r = (0, 1)
        else:
            a = go(x[1])
            if op == 'add':
                b = go(x[2]); r = (a[0] + b[0], a[1] + b[1])
            elif op == 'sub':
                b = go(x[2]); r = (a[0] - b[1], a[1] - b[0])
            elif op == 'mul':
                b = go(x[2]); ps = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]]; r = (min(ps), max(ps))
            elif op == 'div':
                r = (a[0] // x[2], a[1] // x[2])
            elif op == 'mod':
                r = (a[0], a[1]) if 0 <= a[0] and a[1] < x[2] else (0, x[2] - 1)
            elif op == 'shl':
                r = (a[0] << x[2], a[1] << x[2])
            elif op == 'shr':
                r = (a[0] >> x[2], a[1] >> x[2])
            elif op == 'slice':
                j, l = x[2], x[3]
                r = (a[0] >> j, a[1] >> j) if a[0] >= 0 and a[1] < (1 << (j + l)) else (0, (1 << l) - 1)
            elif op in ('bxor', 'bor', 'band'):
                b = go(x[2])
                if a[0] < 0 or b[0] < 0:
                    raise ExprError(f'{op}: operand may be negative')
                r = (0, min(a[1], b[1])) if op == 'band' else (0, (1 << max(a[1].bit_length(), b[1].bit_length())) - 1)
            else:
                raise ExprError(f'unknown operation {op}')
        memo[k] = (x, r)
        return r
    return go(e)


# ----------------------------------------------------------------------------------------------- Z3 (integers)
class Z3Terms:
    """Translator to z3 Int terms.  leaf(leaf_tuple) -> z3 term.  Bitwise operations are expanded bit by bit and are
    refused above ``max_bitwise`` bits (ExprError: the verdict is then UNKNOWN, never PROVED)."""

    def __init__(self, z3, leaf, leaf_bounds, max_bitwise=64):
        self.z3, self.leaf, self.leaf_bounds, self.max_bitwise = z3, leaf, leaf_bounds, max_bitwise
        self.memo, self.bmemo = {}, {}

    def bounds(self, e):
        return bounds(e, self.leaf_bounds, self.bmemo)

    def _bits_of(self, t, j):
        return (t / (1 << j)) % 2 if j else t % 2

    def term(self, e):
        z3 = self.z3
        if type(e) is int:
            return z3.IntVal(e)
        k = id(e)
        if k in self.memo:
            return self.memo[k][1]
        op = e[0]
        one, zero = z3.IntVal(1), z3.IntVal(0)
        if op in _LEAF:
            t = self.leaf(e)
        elif op == 'lut':
            i = self.term(e[1])
            lo, hi = self.bounds(e[1])
            if lo < 0 or hi >= len(e[2]):
                raise ExprError('lut index range outside the table')
            t = z3.IntVal(e[2][hi])
            for v in range(hi - 1, lo - 1, -1):
                t = z3.If(i == v, z3.IntVal(e[2][v]), t)
        elif op == 'ite':
            t = z3.If(self.term(e[1]) == 1, self.term(e[2]), self.term(e[3]))
        elif op == 'not':
            t = 1 - self.term(e[1])
        elif op in _CONST_ARGS:
            a = self.term(e[1])
            c = e[2]
            if op == 'div': t = a / c
            elif op == 'mod': t = a % c
            elif op == 'bit': t = self._bits_of(a, c)
            elif op == 'shl': t = a * (1 << c)
            elif op == 'shr': t = a / (1 << c) if c else a
            else: t = ((a / (1 << c)) if c else a) % (1 << e[3])
        else:
            a, b = self.term(e[1]), self.term(e[2])
            if op == 'add': t = a + b
            elif op == 'sub': t = a - b
            elif op == 'mul':
                if type(e[1]) is int or type(e[2]) is int:
                    t = a * b
                else:
                    la, lb = self.bounds(e[1]), self.bounds(e[2])
                    if la[0] >= 0 and la[1] <= 1:
                        t = z3.If(a == 1, b, zero)
                    elif lb[0] >= 0 and lb[1] <= 1:
                        t = z3.If(b == 1, a, zero)
                    else:
                        raise ExprError('mul of two non-constant, non-boolean factors (non-linear)')
            elif op == 'lt': t = z3.If(a < b, one, zero)
            elif op == 'le': t = z3.If(a <= b, one, zero)
            elif op == 'eq': t = z3.If(a == b, one, zero)
            elif op == 'ne': t = z3.If(a == b, zero, one)
            elif op == 'and': t = z3.If(z3.And(a == 1, b == 1), one, zero)
            elif op == 'or': t = z3.If(z3.Or(a == 1, b == 1), one, zero)
            elif op in ('bxor', 'band', 'bor'):
                if type(e[1]) is int and e[1] == 0:
                    t = b if op != 'band' else zero
                elif type(e[2]) is int and e[2] == 0:
                    t = a if op != 'band' else zero
                else:
                    la, lb = self.bounds(e[1]), self.bounds(e[2])
                    if la[0] < 0 or lb[0] < 0:
                        raise ExprError(f'{op}: operand may be negative')
                    width = max(la[1].bit_length(), lb[1].bit_length(), 1)
                    if width > self.max_bitwise:
                        raise ExprError(f'{op} of {width}-bit operands is not expanded (limit {self.max_bitwise})')
                    acc = zero
                    for j in range(width):
                        x, y = self._bits_of(a, j), self._bits_of(b, j)
                        if op == 'bxor':
                            bj = (x + y) % 2
                        elif op == 'band':
                            bj = z3.If(z3.And(x == 1, y == 1), one, zero)
                        else:
                            bj = z3.If(z3.Or(x == 1, y == 1), one, zero)
                        acc = acc + bj * (1 << j)
                    t = acc
            else:
                raise ExprError(f'unknown operation {op}')
        self.memo[k] = (e, t)
        return t

    def cond(self, e):
        """z3 Bool for a 0/1 expression"""
        return self.term(e) == 1


# ----------------------------------------------------------------------------------------------- ROBDD words
def _need(lo, hi):
    def bl(x):
        return x.bit_length() if x >= 0 else (-x - 1).bit_length()
    return max(bl(lo), bl(hi)) + 1


class BddWord:
    """A two's complement word of ROBDD edges (LSB first, the last bit is the sign) with a static interval."""
    __slots__ = ('bits', 'lo', 'hi')

    def __init__(self, bits, lo, hi):
        self.bits, self.lo, self.hi = list(bits), lo, hi

    def ext(self, w):
        b = self.bits
        if len(b) >= w:
            return b[:w]
        return b + [b[-1]] * (w - len(b))


def bdd_word(e, B, leaf, max_const_mul_bits=600):
    """ROBDD word of e.  leaf(leaf_tuple) -> BddWord (exact, with its static interval).  Raises ExprError for an
    operation that has no exact, small ROBDD form here (division by a non-power of two, a product of two variables,
    a remainder that needs more than 8 conditional subtractions)."""
    from qoptct.verify.bdd import ONE, ZERO
    AND, XOR, ITE, MAJ = B.AND, B.XOR, B.ITE, B.MAJ
    memo = {}

    def const(c):
        return BddWord([ONE if (c >> i) & 1 else ZERO for i in range(_need(c, c))], c, c)

    def addw(a, b, sub_=False):
        lo, hi = (a.lo - b.hi, a.hi - b.lo) if sub_ else (a.lo + b.lo, a.hi + b.hi)
        w = _need(lo, hi)
        x, y = a.ext(w), b.ext(w)
        if sub_:
            y = [t ^ 1 for t in y]
        c = ONE if sub_ else ZERO
        out = []
        for i in range(w):
            out.append(XOR(XOR(x[i], y[i]), c))
            if i + 1 < w:
                c = MAJ(x[i], y[i], c)
        return BddWord(out, lo, hi)

    def shlw(a, k):
        return BddWord([ZERO] * k + a.bits, a.lo << k, a.hi << k)

    def shrw(a, k):
        b = a.bits[k:] if len(a.bits) > k else [a.bits[-1]]
        return BddWord(b, a.lo >> k, a.hi >> k)

    def mux(c, a, b):
        w = max(len(a.bits), len(b.bits))
        x, y = a.ext(w), b.ext(w)
        return BddWord([ITE(c, s, t) for s, t in zip(x, y)], min(a.lo, b.lo), max(a.hi, b.hi))

    def boolw(f):
        return BddWord([f, ZERO], 0, 1)

    def flag(a, what):
        if a.lo < 0 or a.hi > 1:
            raise ExprError(f'{what}: operand is not 0/1 (static range {a.lo}..{a.hi})')
        return a.bits[0]

    def ltw(a, b):
        return addw(a, b, True).bits[-1]

    def eqw(a, b):
        w = max(len(a.bits), len(b.bits))
        x, y = a.ext(w), b.ext(w)
        r = ONE
        for s, t in zip(x, y):
            r = AND(r, XOR(s, t) ^ 1)
        return r

    def low_bits(a, k):          # a mod 2^k (non-negative)
        lo, hi = (a.lo, a.hi) if 0 <= a.lo and a.hi < (1 << k) else (0, (1 << k) - 1)
        return BddWord(a.ext(k) + [ZERO], lo, hi)

    def go(x):
        if type(x) is int:
            return const(x)
        k = id(x)
        if k in memo:
            return memo[k][1]
        op = x[0]
        if op in _LEAF:
            r = leaf(x)
        elif op == 'lut':
            i = go(x[1])
            table = x[2]
            if i.lo < 0 or i.hi >= len(table):
                raise ExprError('lut index range outside the table')
            level = [const(t) for t in table[:i.hi + 1]]
            nbits = max(1, i.hi.bit_length())
            addr = i.ext(nbits)
            for j in range(nbits):
                level = [mux(addr[j], level[q + 1], level[q]) if q + 1 < len(level) else level[q]
                         for q in range(0, len(level), 2)]
            vals = table[i.lo:i.hi + 1]
            r = BddWord(level[0].bits, min(vals), max(vals))
        elif op == 'ite':
            r = mux(flag(go(x[1]), 'ite'), go(x[2]), go(x[3]))
        elif op == 'not':
            r = boolw(flag(go(x[1]), 'not') ^ 1)
        elif op in _CONST_ARGS:
            a = go(x[1])
            c = x[2]
            if op in ('div', 'shr'):
                kk = c if op == 'div' else (1 << c)
                if kk & (kk - 1):
                    raise ExprError('division by a constant that is not a power of two')
                r = shrw(a, kk.bit_length() - 1)
            elif op == 'shl':
                r = shlw(a, c)
            elif op == 'bit':
                r = boolw(a.bits[c] if c < len(a.bits) else a.bits[-1])
            elif op == 'slice':
                s = shrw(a, c)
                r = low_bits(s, x[3])
            else:   # mod
                if not c & (c - 1):
                    r = low_bits(a, c.bit_length() - 1)
                else:
                    if a.lo < 0:
                        raise ExprError('remainder of a possibly negative value by a non-power of two')
                    steps = a.hi // c
                    if steps > 8:
                        raise ExprError('remainder by a non-power of two needs more than 8 subtractions')
                    r = a
                    ck = const(c)
                    for _ in range(steps):
                        hi_next = max(c - 1, r.hi - c)
                        m = mux(ltw(r, ck) ^ 1, addw(r, ck, True), r)
                        r = BddWord(m.bits, 0, hi_next)
                    if a.hi < c:
                        r = a
        else:
            a, b = go(x[1]), go(x[2])
            if op == 'add':
                r = addw(a, b)
            elif op == 'sub':
                r = addw(a, b, True)
            elif op == 'mul':
                if type(x[1]) is int or type(x[2]) is int:
                    c, w_ = (x[1], b) if type(x[1]) is int else (x[2], a)
                    if abs(c).bit_length() > max_const_mul_bits:
                        raise ExprError('constant factor too large')
                    acc, cc, j = const(0), abs(c), 0
                    while cc >> j:
                        if (cc >> j) & 1:
                            acc = addw(acc, shlw(w_, j))
                        j += 1
                    r = addw(const(0), acc, True) if c < 0 else acc
                elif 0 <= a.lo and a.hi <= 1:
                    f = a.bits[0]
                    r = BddWord([AND(f, t) for t in b.bits], min(0, b.lo), max(0, b.hi))
                elif 0 <= b.lo and b.hi <= 1:
                    f = b.bits[0]
                    r = BddWord([AND(f, t) for t in a.bits], min(0, a.lo), max(0, a.hi))
                else:
                    raise ExprError('mul of two non-constant, non-boolean factors has no small ROBDD here')
            elif op in ('lt', 'le', 'eq', 'ne'):
                f = (ltw(a, b) if op == 'lt' else ltw(b, a) ^ 1 if op == 'le' else
                     eqw(a, b) if op == 'eq' else eqw(a, b) ^ 1)
                r = boolw(f)
            elif op in ('and', 'or'):
                fa, fb = flag(a, op), flag(b, op)
                r = boolw(AND(fa, fb) if op == 'and' else B.OR(fa, fb))
            elif op in ('bxor', 'band', 'bor'):
                if a.lo < 0 or b.lo < 0:
                    raise ExprError(f'{op}: operand may be negative')
                w = max(len(a.bits), len(b.bits))
                xa, xb = a.ext(w), b.ext(w)
                fn = XOR if op == 'bxor' else AND if op == 'band' else B.OR
                hi = min(a.hi, b.hi) if op == 'band' else (1 << max(a.hi.bit_length(), b.hi.bit_length())) - 1
                r = BddWord([fn(s, t) for s, t in zip(xa, xb)] + [ZERO], 0, hi)
            else:
                raise ExprError(f'unknown operation {op}')
        memo[k] = (x, r)
        return r
    return go(e)
