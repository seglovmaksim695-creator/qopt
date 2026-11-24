"""Binary_ECC (Jang, Kim, Seo et al., ProjectQ circuits for binary elliptic curves) as ProjectQ programs.

The source files are NOT part of this repository (Eclipse Public License 2.0).  Point ``source`` (or the environment
variable ``QOPT_BINARY_ECC``) at a checkout of the ``Binary_ECC`` folder.  The source text is read, the part before the
first ``global resource_check`` line (imports and function definitions only) is executed unchanged, and the circuit
functions are called with the flags the authors use for resource counting (``NCT=1``: X / CNOT / Toffoli only, no
measurements).  Nothing in the source is edited.  The SHA-256 of the file is part of every result.
"""
from __future__ import annotations

import hashlib
import os
import re
from pathlib import Path

from ._record import Unsupported

__all__ = ['POLYNOMIALS', 'source_dir', 'multiplier_program', 'point_addition_program', 'point_addition_spec',
           'multiplier_spec', 'file_sha256']

POLYNOMIALS = {8: (8, 4, 3, 1, 0), 16: (16, 5, 3, 1, 0), 127: (127, 1, 0), 163: (163, 7, 6, 3, 0), 233: (233, 74, 0),
               283: (283, 12, 7, 5, 0), 571: (571, 10, 5, 2, 0)}


def source_dir(source=None):
    p = source or os.environ.get('QOPT_BINARY_ECC')
    if not p:
        raise FileNotFoundError('Binary_ECC sources not given: pass source=... or set QOPT_BINARY_ECC to the Binary_ECC folder')
    p = Path(p).expanduser()
    if not (p / 'Multiplication').is_dir():
        raise FileNotFoundError(f'{p} is not a Binary_ECC folder (no Multiplication/)')
    return p


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _namespace(path):
    text = Path(path).read_text()
    cut = text.index('global resource_check')
    ns = {'__name__': 'binary_ecc_source'}
    exec(compile(text[:cut], str(path), 'exec'), ns)  # noqa: S102 - the user's own checkout of the sources
    ns['resource_check'] = 1
    ns['NCT'] = 1
    ns['AND_check'] = 0
    return text, ns


def multiplier_path(m, source=None):
    return source_dir(source) / 'Multiplication' / f'Mul_x{m}.py'


def multiplier_program(m, source=None):
    """Program of ``Multiplication/Mul_x{m}.py``: inputs x, y (m qubits), outputs x, y and the low m bits of the reduced
    product as ``result``; the remaining wires are the clean ancillas and the (pinned) garbage of the Karatsuba scheme."""
    path = multiplier_path(m, source)
    text, ns = _namespace(path)
    fn = text[text.index('def Multiplication'):text.index('def Reduction')]
    n = int(re.search(r'\n\s+n = (\d+)', fn).group(1))
    anc = int(re.search(r'ancilla = eng.allocate_qureg\((\d+)\)', fn).group(1))
    if n != m:
        raise Unsupported(f'{path.name} is a multiplier of width {n}, not {m}')

    def program(eng):
        x1 = eng.allocate_qureg(n)
        y1 = eng.allocate_qureg(n)
        ancilla = eng.allocate_qureg(anc)
        result, _count, ancilla = ns['recursive_karatsuba'](eng, x1, y1, n, 0, ancilla)
        result = ns['Reduction'](eng, result)
        return dict(inputs=dict(x=x1, y=y1), outputs=dict(x=x1, y=y1, result=result[:n]))
    return program


def multiplier_spec(m):
    from .check import GF2MulSpec
    return GF2MulSpec(m, sum(1 << e for e in POLYNOMIALS[m]))


def point_addition_path(m, variant='FLT_in', source=None):
    return source_dir(source) / 'Point_additions' / variant / f'x^{m}' / f'Point_addition_x{m}.py'


def _constants(text):
    g = lambda name: int(re.search(rf'\n\s+{name} = (0x[0-9a-fA-F]+|\d+)', text).group(1), 0)  # noqa: E731
    return g('a'), g('x2'), g('y2'), g('n')


def point_addition_program(m, variant='FLT_in', source=None):
    """Program of ``Point_additions/{variant}/x^{m}/Point_addition_x{m}.py`` with the registers q (1 qubit: 0 = copy, 1 = add),
    x1, y1 (m qubits each): inputs and outputs.  The registers are the first allocations of the source, in this order."""
    path = point_addition_path(m, variant, source)
    text, ns = _namespace(path)
    a, x2, y2, n = _constants(text[text.index('def Point_addition'):])
    if n != m:
        raise Unsupported(f'{path.name} has width {n}, not {m}')

    def program(eng):
        flat = []
        real_qubit, real_qureg = eng.allocate_qubit, eng.allocate_qureg

        def allocate_qubit():
            q = real_qubit()
            flat.extend(q)
            return q

        def allocate_qureg(k):
            q = real_qureg(k)
            flat.extend(q)
            return q
        eng.allocate_qubit, eng.allocate_qureg = allocate_qubit, allocate_qureg
        try:
            ns['Point_addition'](eng)
        finally:
            del eng.allocate_qubit, eng.allocate_qureg
        seen, order = set(), []                  # allocate_qureg calls allocate_qubit: keep the first record of each qubit
        for q in flat:
            if q.id not in seen:
                seen.add(q.id)
                order.append(q)
        q0, x1, y1 = order[:1], order[1:1 + m], order[1 + m:1 + 2 * m]
        regs = dict(q=q0, x1=x1, y1=y1)
        return dict(inputs=regs, outputs=dict(regs), _keep=order)
    program.constants = dict(a=a, x2=x2, y2=y2, n=n)
    return program


def point_addition_spec(m, program, poly=None):
    """Specification of the point addition of (x1, y1) and the constant point (x2, y2) of ``program`` (see
    ``point_addition_program``): q = 0: unchanged;  q = 1 and x1 != x2: (x3, y3) with lambda = (y1+y2)/(x1+x2),
    x3 = lambda^2 + lambda + x1 + x2 + a, y3 = lambda (x1 + x3) + x3 + y1.  The claim is restricted to the domain
    'q = 0 or x1 != x2' (the circuit inverts by Fermat's little theorem; at x1 = x2 the value is not a point addition)."""
    from ._classical import ClassicalSpec
    from ._gf2 import mulmod
    poly = poly or sum(1 << e for e in POLYNOMIALS[m])
    a, x2, y2 = program.constants['a'], program.constants['x2'], program.constants['y2']

    def inv(v):
        r, b, e = 1, v, (1 << m) - 2
        while e:
            if e & 1:
                r = mulmod(r, b, poly)
            b = mulmod(b, b, poly)
            e >>= 1
        return r

    def fn(v):
        q, x1, y1 = v['q'], v['x1'], v['y1']
        if q == 0:
            return dict(q=0, x1=x1, y1=y1)
        lam = mulmod(y1 ^ y2, inv(x1 ^ x2), poly)
        x3 = mulmod(lam, lam, poly) ^ lam ^ x1 ^ x2 ^ a
        y3 = mulmod(lam, x1 ^ x3, poly) ^ x3 ^ y1
        return dict(q=1, x1=x3, y1=y3)
    return ClassicalSpec(fn, promise=lambda v: v['q'] == 0 or v['x1'] != x2, name=f'binary_point_add_x{m}')
