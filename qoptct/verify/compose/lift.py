# SPDX-License-Identifier: Apache-2.0
"""compose.lift -- OpenQASM text -> sequence of logical operations (x, cx, AND, ERASE), strictly.

The text must be a concatenation of exactly these line groups (the dialect of qoptct.crypto.kaliski.to_qasm):

    x q[a];                                                   ('x', a)
    cx q[a], q[b];                                            ('cx', a, b)
    h q[t]; t q[t]; cx q[b], q[t]; tdg q[t]; cx q[a], q[t];
    t q[t]; cx q[b], q[t]; tdg q[t]; h q[t]; sdg q[t];        ('and', a, b, t)      (10 lines, a, b, t distinct)
    h q[t]; c[k] = measure q[t];
    if (c[k]) {cz q[a], q[b]; [z q[l];] x q[t];}              ('erase', a, b, t, l) (3 lines, wires distinct)

Anything else raises LiftError (the circuit is then not a composition of the gadgets; nothing is claimed about it).

GADGET LEMMAS (checked by exact arithmetic in Q(zeta_8) by check_gadget_lemmas(), on the very text templates the
lifter recognises):
  AND:    the 10 gates map |a, b, 0> to |a, b, a*b> exactly (amplitude 1), for all a, b.
  ERASE:  for every basis state |a, b, l, t> with t = a*b xor l and for BOTH measurement outcomes r, the three lines
          give (1/sqrt 2) |a, b, l, 0> exactly (the same vector for r = 0 and r = 1).
Consequence (linearity; reasoning, three lines): on an input basis state on which the classical run of the lifted
program has every AND target equal to 0 and every erase condition t = a*b xor l true, the circuit yields, for every
outcome string, 2^(-m/2) |classical result>; hence on any superposition of such basis states the circuit acts as the
classical reversible map, with no phase and no dependence on the outcomes.

Compact form of the result: array('i') with 5 ints per operation (kind, a, b, t, l):
    kind 0: x a        1: cx a -> b       2: and (a, b) -> t       3: erase t against a*b xor l   (l = -1: none)
"""
import io
import os
import re
from array import array
from fractions import Fraction

K_X, K_CX, K_AND, K_ERASE = 0, 1, 2, 3
__all__ = ['lift', 'Lifted', 'LiftError', 'check_gadget_lemmas', 'and_template', 'erase_template',
           'K_X', 'K_CX', 'K_AND', 'K_ERASE']


class LiftError(Exception):
    pass


class Lifted:
    """nq: number of qubits; ops: array('i'), 5 ints per operation; n: number of operations."""
    __slots__ = ('nq', 'ops', 'n', 'nmeas')

    def __init__(self, nq, ops, nmeas):
        self.nq, self.ops, self.n, self.nmeas = nq, ops, len(ops) // 5, nmeas

    def op(self, i):
        return tuple(self.ops[5 * i:5 * i + 5])


def and_template(a, b, t):
    q = 'q[%d];' % t
    cb, ca = 'cx q[%d], q[%d];' % (b, t), 'cx q[%d], q[%d];' % (a, t)
    return ['h ' + q, 't ' + q, cb, 'tdg ' + q, ca, 't ' + q, cb, 'tdg ' + q, 'h ' + q, 'sdg ' + q]


def erase_template(a, b, t, lin=(), k=0):
    body = ['cz q[%d], q[%d];' % (a, b)] + ['z q[%d];' % w for w in lin] + ['x q[%d];' % t]
    return ['h q[%d];' % t, 'c[%d] = measure q[%d];' % (k, t), 'if (c[%d]) {%s}' % (k, ' '.join(body))]


_RE_IF = re.compile(r'if \(c\[(\d+)\]\) \{cz q\[(\d+)\], q\[(\d+)\];((?: z q\[\d+\];)*) x q\[(\d+)\];\}\n\Z')
_RE_Z = re.compile(r' z q\[(\d+)\];')
_RE_Q = re.compile(r'qubit\[(\d+)\] q;\n\Z')
_RE_B = re.compile(r'bit\[(\d+)\] c;\n\Z')


def lift(src):
    """src: path of a file, or the OpenQASM text itself (if it contains a newline). -> Lifted. Raises LiftError."""
    if isinstance(src, str) and '\n' not in src:
        if not os.path.exists(src):
            raise LiftError('no such file: %s' % src)
        f = open(src)
    else:
        f = io.StringIO(src)
    with f:
        return _lift(f)


def _lift(f):
    lineno = [4]

    def bad(msg):
        raise LiftError('line %d: %s' % (lineno[0], msg))

    def nxt():
        ln = f.readline()
        lineno[0] += 1
        if not ln:
            bad('unexpected end of file inside a gadget')
        return ln if ln[-1] == '\n' else ln + '\n'
    hdr = [f.readline() for _ in range(4)]
    if hdr[0].strip() not in ('OPENQASM 3;', 'OPENQASM 3.0;') or hdr[1].strip() != 'include "stdgates.inc";':
        raise LiftError('header: expected OPENQASM 3 with stdgates.inc')
    m, mb = _RE_Q.match(hdr[2]), _RE_B.match(hdr[3])
    if not m or not mb:
        raise LiftError('header: expected "qubit[N] q;" and "bit[M] c;"')
    nq = int(m.group(1))
    ops = array('i')
    app = ops.extend
    nmeas = 0
    nline = 4
    while True:
        ln = f.readline()
        if not ln:
            break
        nline += 1
        if ln[-1] != '\n':
            ln += '\n'
        c0 = ln[0]
        try:
            if c0 == 'c':
                a, b = ln[5:-3].split('], q[')
                a = int(a); b = int(b)
                if ln != 'cx q[%d], q[%d];\n' % (a, b) or a == b or a >= nq or b >= nq:
                    lineno[0] = nline; bad('not a recognised line: %r' % ln[:60])
                app((K_CX, a, b, 0, -1))
            elif c0 == 'x':
                a = int(ln[4:-3])
                if ln != 'x q[%d];\n' % a or a >= nq:
                    lineno[0] = nline; bad('not a recognised line: %r' % ln[:60])
                app((K_X, a, 0, 0, -1))
            elif c0 == 'h':
                lineno[0] = nline
                t = int(ln[4:-3])
                q = 'q[%d];\n' % t
                if ln != 'h ' + q or t >= nq:
                    bad('not a recognised line: %r' % ln[:60])
                l2 = nxt()
                if l2 == 't ' + q:
                    cb = nxt()
                    b = int(cb[5:-3].split('], q[')[0])
                    if cb != 'cx q[%d], q[%d];\n' % (b, t):
                        bad('AND gadget: bad line 3')
                    if nxt() != 'tdg ' + q:
                        bad('AND gadget: bad line 4')
                    ca = nxt()
                    a = int(ca[5:-3].split('], q[')[0])
                    if ca != 'cx q[%d], q[%d];\n' % (a, t):
                        bad('AND gadget: bad line 5')
                    for j, exp in enumerate(('t ' + q, cb, 'tdg ' + q, 'h ' + q, 'sdg ' + q)):
                        if nxt() != exp:
                            bad('AND gadget: bad line %d' % (6 + j))
                    if len({a, b, t}) != 3 or a >= nq or b >= nq:
                        bad('AND gadget: wires not distinct')
                    app((K_AND, a, b, t, -1))
                    nline += 9
                elif l2.startswith('c['):
                    k = int(l2[2:l2.index(']')])
                    if l2 != 'c[%d] = measure q[%d];\n' % (k, t):
                        bad('ERASE gadget: bad measurement line')
                    l3 = nxt()
                    mm = _RE_IF.match(l3)
                    if not mm or int(mm.group(1)) != k or int(mm.group(5)) != t:
                        bad('ERASE gadget: bad correction line')
                    a, b = int(mm.group(2)), int(mm.group(3))
                    lin = [int(x) for x in _RE_Z.findall(mm.group(4))]
                    if len(lin) > 1:
                        bad('ERASE gadget: more than one linear correction (unsupported)')
                    if len({a, b, t, *lin}) != 3 + len(lin) or max([a, b] + lin) >= nq:
                        bad('ERASE gadget: wires not distinct')
                    app((K_ERASE, a, b, t, lin[0] if lin else -1))
                    nmeas += 1
                    nline += 2
                else:
                    bad('h is not the start of an AND or ERASE gadget')
            else:
                lineno[0] = nline; bad('not a recognised line: %r' % ln[:60])
        except (ValueError, IndexError):
            lineno[0] = max(lineno[0], nline)
            raise LiftError('line %d: not a recognised line' % lineno[0])
    return Lifted(nq, ops, nmeas)


# ------------------------------------------------------------------ exact gadget lemmas, arithmetic in Q(zeta_8)
def _mul(x, y):
    out = [Fraction(0)] * 4
    for i in range(4):
        if x[i]:
            for j in range(4):
                if y[j]:
                    k = i + j
                    out[k & 3] += (-1 if k >= 4 else 1) * x[i] * y[j]
    return tuple(out)


_F0, _F1 = Fraction(0), Fraction(1)
_ONE = (_F1, _F0, _F0, _F0)
_R2 = (_F0, Fraction(1, 2), _F0, Fraction(-1, 2))          # 1/sqrt(2) = (w - w^3)/2, w = exp(i pi/4)
_PH = {'t': (_F0, _F1, _F0, _F0), 'tdg': (_F0, _F0, _F0, -_F1), 's': (_F0, _F0, _F1, _F0),
       'sdg': (_F0, _F0, -_F1, _F0), 'z': (-_F1, _F0, _F0, _F0)}


def _add_amp(st, key, amp):
    cur = st.get(key)
    amp = amp if cur is None else tuple(a + b for a, b in zip(cur, amp))
    if any(amp):
        st[key] = amp
    else:
        st.pop(key, None)


def _apply(st, name, qs):
    out = {}
    for key, amp in st.items():
        if name == 'h':
            q = qs[0]
            a = _mul(amp, _R2)
            _add_amp(out, key[:q] + (0,) + key[q + 1:], a)
            _add_amp(out, key[:q] + (1,) + key[q + 1:], tuple(-v for v in a) if key[q] else a)
        elif name in _PH:
            _add_amp(out, key, _mul(amp, _PH[name]) if key[qs[0]] else amp)
        elif name == 'x':
            q = qs[0]
            _add_amp(out, key[:q] + (1 - key[q],) + key[q + 1:], amp)
        elif name == 'cx':
            c, q = qs
            _add_amp(out, key[:q] + (key[q] ^ key[c],) + key[q + 1:], amp)
        elif name == 'cz':
            _add_amp(out, key, tuple(-v for v in amp) if key[qs[0]] and key[qs[1]] else amp)
        else:
            raise ValueError(name)
    return out


_RE_G = re.compile(r'(\w+) q\[(\d+)\](?:, q\[(\d+)\])?;')


def _run_unitary(st, lines):
    for ln in lines:
        m = _RE_G.fullmatch(ln)
        st = _apply(st, m.group(1), [int(m.group(2))] + ([int(m.group(3))] if m.group(3) else []))
    return st


def check_gadget_lemmas():
    """exact check of the two gadget lemmas of the module docstring. -> True, or raises AssertionError."""
    for a in (0, 1):
        for b in (0, 1):
            st = _run_unitary({(a, b, 0): _ONE}, and_template(0, 1, 2))
            assert st == {(a, b, a & b): _ONE}, ('AND lemma', a, b, st)
    for nl in (0, 1):
        lin = (2,) if nl else ()
        t = 2 + nl
        tpl = erase_template(0, 1, t, lin)
        body = re.match(r'if \(c\[0\]\) \{(.*)\}', tpl[2]).group(1)
        blines = [x.strip() + ';' for x in body.split(';') if x.strip()]
        for a in (0, 1):
            for b in (0, 1):
                for l in ((0, 1) if nl else (0,)):
                    key = (a, b) + ((l,) if nl else ()) + ((a & b) ^ l,)
                    st = _run_unitary({key: _ONE}, tpl[:1])
                    for r in (0, 1):
                        br = {k: v for k, v in st.items() if k[t] == r}
                        if r:
                            br = _run_unitary(br, blines)
                        want = (a, b) + ((l,) if nl else ()) + (0,)
                        assert br == {want: _R2}, ('ERASE lemma', a, b, l, r, br)
    return True
