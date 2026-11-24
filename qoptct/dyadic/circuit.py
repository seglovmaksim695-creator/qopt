# SPDX-License-Identifier: Apache-2.0
"""Circuits with dyadic angles for qoptct.dyadic: exact angles, gate set, expansion to primitives.

Angles.  An angle is kept as an exact multiple of pi: theta = pi * r with r a Fraction (reduced modulo 2).  Inputs:
  * Pi(p, q)            exact: pi * p / q;
  * str                 an expression in radians such as "pi/4", "-3*pi/8", "pi*2**-5", "0.5*pi" (own parser, no eval);
  Angles are kept modulo 4*pi: half-angle gates (RZ, RX, RY, U, CRZ, RZZ) have period 4*pi.
  * float / int         radians.  A float is accepted only if theta/pi, reduced to (-1, 1], equals p/2^k with
                        0 < |p| <= 2^20 within 32 units in the last place of theta/pi (FLOAT_MAX_NUM, FLOAT_ULPS);
                        the count and the largest deviation of recognised float angles are recorded.  Every float is
                        itself a dyadic rational, so without the bound on p any angle would be "recognised"; the bound
                        keeps the reading unambiguous.  The proof is about the circuit with the recognised angles.
A circuit whose angle is not a dyadic multiple of pi (or not recognised) is rejected (NotDyadic): the API answers
UNKNOWN, never PROVED.

Gates (Qiskit names and conventions; qubit order as listed, controls first):
  id x y z h s sdg t tdg sx sxdg                      p(l) u1(l) phase(l) rz(l) rx(l) ry(l) u(t,f,l) u3 u2(f,l)
  gphase(l)  cx cnot cy cz cs csdg csx cp(l) cu1 cphase crz(l) swap dcx rzz(l)
  ccx toffoli ccz cswap fredkin  mcx mct c3x c4x (k controls, target)  mcz  mcp(l) mcphase mcu1 (controls, target)
Every gate is expanded into the exact primitives of qoptct.dyadic.core (Hadamard, XOR, swap, cswap, phase of a
product of wires, phase of an xor of wires); the identities are listed in `expand`.  The inverse of a circuit is the
reversed primitive list with every phase negated (H, XOR, swap, cswap are self-inverse), so no gate-level inverse
table is needed.
"""
import math
import re
from fractions import Fraction

__all__ = ['Pi', 'NotDyadic', 'Unsupported', 'DCircuit', 'to_angle', 'expand', 'invert_prims', 'needed_K',
           'parse_qasm', 'from_qc', 'GATES', 'FLOAT_MAX_NUM', 'FLOAT_ULPS', 'recognise_float']

FLOAT_MAX_NUM = 1 << 20       # largest numerator p of a float angle pi*p/2^k (after reduction to (-pi, pi])
FLOAT_ULPS = 32               # tolerance of a float angle, in units in the last place of theta/pi
AMBIGUITY_RATIO = 64          # with an absolute tolerance (rounding source), nonzero angles must exceed 64 x tolerance
MAX_K = 4096                  # largest 2-power denominator accepted for an exact angle


class NotDyadic(ValueError):
    """An angle is not (recognised as) a dyadic multiple of pi."""


class Unsupported(ValueError):
    """A gate or statement outside the supported set (measurement, reset, classical control, ...)."""


class Pi:
    """Exact angle pi * p / q."""
    __slots__ = ('r',)

    def __init__(self, p, q=1):
        self.r = Fraction(p) / Fraction(q)

    def __repr__(self):
        return 'Pi(%s)' % self.r


def _dyadic(r):
    d = r.denominator
    if d & (d - 1):
        raise NotDyadic('angle pi*%s is not a dyadic multiple of pi' % r)
    if d.bit_length() - 1 > MAX_K:
        raise NotDyadic('angle denominator 2^%d exceeds the limit 2^%d' % (d.bit_length() - 1, MAX_K))
    return r


def _mod2(r):
    """r modulo 2 (exp(i pi r) has period 2 in r), in [0, 2): for phases (primitives, global phase)."""
    return r - 2 * (r.numerator // (2 * r.denominator))


def _mod4(r):
    """r modulo 4, in [0, 4): gate angles.  Half-angle gates (RZ, RX, RY, U, CRZ, RZZ) have period 4*pi, not 2*pi:
    RZ(t + 2 pi) = -RZ(t) and CRZ(t + 2 pi) differs from CRZ(t) by a relative phase, so angles are never reduced
    modulo 2*pi."""
    return r - 4 * (r.numerator // (4 * r.denominator))


def recognise_float(theta, rec=None, abs_tol=None):
    """theta (radians, float) -> Fraction r in [0, 4) with theta = pi*r up to the float precision; NotDyadic if none.
    x = theta/pi is reduced to (-2, 2]; the tolerance is FLOAT_ULPS units in the last place of x (a float carries no
    more information than that), or abs_tol radians if larger (sources that round angles, e.g. ProjectQ keeps 12
    decimals); the simplest candidate p/2^k (smallest k) with 0 < |p| <= FLOAT_MAX_NUM within the tolerance is
    taken; x within the tolerance of a multiple of 4 reads 0.  With abs_tol, a nonzero angle closer than
    AMBIGUITY_RATIO * abs_tol to 0 or +-2 pi is refused (it could be several dyadic angles, or 0)."""
    theta = float(theta)
    if theta != theta or theta in (float('inf'), float('-inf')):
        raise NotDyadic('angle %r is not finite' % theta)
    if theta == 0.0:
        return Fraction(0)
    x = theta / math.pi
    tol = FLOAT_ULPS * math.ulp(x)
    if abs_tol:
        tol = max(tol, abs_tol / math.pi)
    m = round(x / 4.0)
    xr = x - 4.0 * m
    r = dev = None
    amb = AMBIGUITY_RATIO * tol
    for o in (0, 2, -2):
        # x = o + p/2^k with o in {0, 2, -2}: near +-2 (angles near +-2 pi of the half-angle gates) the simple part
        # is the distance to +-2
        y = xr - o
        if o and abs(y) > 1.0:
            continue
        if abs_tol and 0 < abs(y) < amb:
            # with a coarse tolerance a small angle could be several dyadic values (or 0): refuse
            raise NotDyadic('float angle %r is below %d x the tolerance %g of its source' % (theta, AMBIGUITY_RATIO,
                                                                                       abs_tol))
        if abs(y) <= tol:
            r, dev = Fraction(o), abs(y)
            break
        for k in range(0, 1100):
            p = round(math.ldexp(y, k))
            if abs(p) > FLOAT_MAX_NUM:
                break
            if p and abs(y - math.ldexp(p, -k)) <= tol:
                r, dev = o + Fraction(p, 1 << k), abs(y - math.ldexp(p, -k))
                break
        if r is not None:
            break
    if r is None:
        raise NotDyadic('float angle %r is not recognised as pi*p/2^k with |p| <= 2^20' % theta)
    if rec is not None:
        rec['n_float'] = rec.get('n_float', 0) + 1
        rec['max_abs_dev'] = max(rec.get('max_abs_dev', 0.0), dev * math.pi)
    return _mod4(r)


def to_angle(v, rec=None, abs_tol=None):
    """Angle -> Fraction r (multiple of pi), reduced modulo 4 (period 4*pi of the half-angle gates)."""
    if isinstance(v, Pi):
        return _mod4(_dyadic(v.r))
    if isinstance(v, str):
        return expr_angle(v, rec=rec, abs_tol=abs_tol)
    if isinstance(v, (Fraction, bool)):
        raise TypeError('bad angle %r: use Pi(p, q) for pi*p/q, a float for radians or a string' % (v,))
    if isinstance(v, int):
        return Fraction(0) if v == 0 else recognise_float(float(v), rec, abs_tol)
    try:
        return recognise_float(float(v), rec, abs_tol)        # floats, numpy floats
    except (TypeError, ValueError) as e:
        if isinstance(e, NotDyadic):
            raise
        raise TypeError('bad angle %r (%s)' % (v, e))


# ----------------------------------------------------------------------------- angle expressions (no eval)
class _Val:
    """a * pi + b (a rational; b rational or, after a float literal, a float)."""
    __slots__ = ('a', 'b')

    def __init__(self, a=Fraction(0), b=Fraction(0)):
        self.a, self.b = a, b

    def add(self, o, s=1):
        return _Val(self.a + s * o.a, self.b + s * o.b)

    def mul(self, o):
        if self.a and o.a:
            raise NotDyadic('pi*pi in an angle expression')
        return _Val(self.a * o.b + o.a * self.b, self.b * o.b)

    def div(self, o):
        if o.a:
            if o.b or self.b:
                raise NotDyadic('division by an expression with pi')
            return _Val(Fraction(0), Fraction(self.a) / Fraction(o.a))
        if not o.b:
            raise NotDyadic('division by zero in an angle')
        return _Val(self.a / o.b, self.b / o.b)


_TOK = re.compile(r'\s*(?:(\d+\.\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?|\d+[eE][-+]?\d+)|(\d+)|(pi|π|tau|τ)'
                  r'|([A-Za-z_][A-Za-z_0-9]*)|(\*\*|[-+*/^()]))')


def _tokens(s):
    pos, out = 0, []
    s = s.strip()
    while pos < len(s):
        m = _TOK.match(s, pos)
        if not m or m.end() == pos:
            raise NotDyadic('cannot read angle expression %r' % s)
        pos = m.end()
        fl, it, pi, name, op = m.groups()
        if fl is not None:
            out.append(('num', float(fl)))
        elif it is not None:
            out.append(('num', Fraction(int(it))))
        elif pi is not None:
            out.append(('pi', 2 if pi in ('tau', 'τ') else 1))
        elif name is not None:
            out.append(('name', name))
        else:
            out.append(('op', '**' if op == '^' else op))
    return out


def eval_expr(s, env=None):
    """Angle expression -> _Val.  env: name -> _Val (formal parameters of a gate definition)."""
    toks = _tokens(s)
    pos = [0]

    def peek():
        return toks[pos[0]] if pos[0] < len(toks) else (None, None)

    def eat():
        t = toks[pos[0]]
        pos[0] += 1
        return t

    def expr():
        v = term()
        while peek() in (('op', '+'), ('op', '-')):
            _, o = eat()
            v = v.add(term(), 1 if o == '+' else -1)
        return v

    def term():
        v = unary()
        while peek() in (('op', '*'), ('op', '/')):
            _, o = eat()
            w = unary()
            v = v.mul(w) if o == '*' else v.div(w)
        return v

    def unary():
        if peek() == ('op', '-'):
            eat()
            return _Val().add(unary(), -1)
        if peek() == ('op', '+'):
            eat()
            return unary()
        return power()

    def power():
        v = atom()
        if peek() == ('op', '**'):
            eat()
            e = unary()
            if v.a or e.a or not isinstance(e.b, Fraction) or e.b.denominator != 1 or not isinstance(v.b, Fraction):
                raise NotDyadic('only rational ** integer is supported in angles')
            if abs(e.b) > 4096:
                raise NotDyadic('exponent too large in an angle')
            if v.b == 0 and e.b < 0:
                raise NotDyadic('0 ** negative in an angle')
            return _Val(Fraction(0), v.b ** int(e.b))
        return v

    def atom():
        kind, val = peek()
        if kind is None:
            raise NotDyadic('truncated angle expression %r' % s)
        eat()
        if kind == 'num':
            return _Val(Fraction(0), val)
        if kind == 'pi':
            return _Val(Fraction(val), Fraction(0))
        if kind == 'name':
            if env is not None and val in env:
                return env[val]
            raise NotDyadic('unknown name %r in angle %r' % (val, s))
        if val == '(':
            v = expr()
            if pos[0] >= len(toks) or eat() != ('op', ')'):
                raise NotDyadic('missing ) in %r' % s)
            return v
        raise NotDyadic('unexpected %r in angle %r' % (val, s))
    v = expr()
    if pos[0] != len(toks):
        raise NotDyadic('trailing tokens in angle %r' % s)
    return v


def val_angle(v, rec=None, abs_tol=None):
    if isinstance(v.b, Fraction) and v.b == 0:
        return _mod4(_dyadic(Fraction(v.a)))
    return recognise_float(float(v.a) * math.pi + float(v.b), rec, abs_tol)


def expr_angle(s, env=None, rec=None, abs_tol=None):
    return val_angle(eval_expr(s, env), rec, abs_tol)


# ----------------------------------------------------------------------------- gate set
# name -> (number of angle parameters, number of qubits or None = any number >= 1)
GATES = {
    'id': (0, 1), 'i': (0, 1), 'x': (0, 1), 'y': (0, 1), 'z': (0, 1), 'h': (0, 1), 's': (0, 1), 'sdg': (0, 1),
    't': (0, 1), 'tdg': (0, 1), 'sx': (0, 1), 'sxdg': (0, 1),
    'p': (1, 1), 'u1': (1, 1), 'phase': (1, 1), 'rz': (1, 1), 'rx': (1, 1), 'ry': (1, 1),
    'u': (3, 1), 'u3': (3, 1), 'u2': (2, 1), 'gphase': (1, 0),
    'cx': (0, 2), 'cnot': (0, 2), 'cy': (0, 2), 'cz': (0, 2), 'cs': (0, 2), 'csdg': (0, 2), 'csx': (0, 2),
    'cp': (1, 2), 'cu1': (1, 2), 'cphase': (1, 2), 'crz': (1, 2), 'swap': (0, 2), 'dcx': (0, 2), 'rzz': (1, 2),
    'ccx': (0, 3), 'toffoli': (0, 3), 'ccz': (0, 3), 'cswap': (0, 3), 'fredkin': (0, 3), 'c3x': (0, 4),
    'c4x': (0, 5), 'mcx': (0, None), 'mct': (0, None), 'mcz': (0, None), 'mcp': (1, None), 'mcphase': (1, None),
    'mcu1': (1, None),
}
_H = Fraction(1, 2)


def expand(name, qs, rs):
    """Exact expansion of one gate into primitives.  Identities used (time order left to right):
      Y = i X Z: phase i, phase (-1)^q, then X          SX = H S H, SXdg = H Sdg H
      RZ(t) = exp(-i t/2) P(t)                          RX(t) = H RZ(t) H      RY(t) = S RX(t) Sdg
      U(t, f, l) = P(f) RY(t) P(l)  (Qiskit U)         U2(f, l) = U(pi/2, f, l)
      CY = S_t CX Sdg_t      CSX = H_t CS H_t           CRZ(t) = P_c(-t/2) CP(t)    DCX = CX(a,b) then CX(b,a)
      RZZ(t) = exp(-i t/2) exp(i t [a xor b])           CP, MCP, CZ, CCZ, MCZ: phase of the AND of all wires."""
    def ph(qq, r):
        return ('ph', tuple(qq), r)

    def rz(q, r):
        return [ph((), -r / 2), ph((q,), r)]

    def ry(q, r):
        return [ph((q,), -_H), ('h', q)] + rz(q, r) + [('h', q), ph((q,), _H)]
    q = qs[0] if qs else None
    if name in ('id', 'i'):
        return []
    if name == 'x':
        return [('xor', q, ())]
    if name == 'y':
        return [ph((), _H), ph((q,), Fraction(1)), ('xor', q, ())]
    if name in ('z', 's', 'sdg', 't', 'tdg'):
        return [ph((q,), {'z': Fraction(1), 's': _H, 'sdg': -_H, 't': Fraction(1, 4), 'tdg': Fraction(-1, 4)}[name])]
    if name == 'h':
        return [('h', q)]
    if name in ('sx', 'sxdg'):
        return [('h', q), ph((q,), _H if name == 'sx' else -_H), ('h', q)]
    if name in ('p', 'u1', 'phase'):
        return [ph((q,), rs[0])]
    if name == 'rz':
        return rz(q, rs[0])
    if name == 'rx':
        return [('h', q)] + rz(q, rs[0]) + [('h', q)]
    if name == 'ry':
        return ry(q, rs[0])
    if name in ('u', 'u3'):
        return [ph((q,), rs[2])] + ry(q, rs[0]) + [ph((q,), rs[1])]
    if name == 'u2':
        return [ph((q,), rs[1])] + ry(q, _H) + [ph((q,), rs[0])]
    if name == 'gphase':
        return [ph((), rs[0])]
    if name in ('cx', 'cnot'):
        return [('xor', qs[1], (qs[0],))]
    if name == 'cy':
        return [ph((qs[1],), -_H), ('xor', qs[1], (qs[0],)), ph((qs[1],), _H)]
    if name in ('cz', 'cs', 'csdg'):
        return [ph(qs, {'cz': Fraction(1), 'cs': _H, 'csdg': -_H}[name])]
    if name == 'csx':
        return [('h', qs[1]), ph(qs, _H), ('h', qs[1])]
    if name in ('cp', 'cu1', 'cphase', 'mcp', 'mcphase', 'mcu1'):
        return [ph(qs, rs[0])]
    if name == 'crz':
        return [ph((qs[0],), -rs[0] / 2), ph(qs, rs[0])]
    if name == 'swap':
        return [('swap', qs[0], qs[1])]
    if name == 'dcx':
        return [('xor', qs[1], (qs[0],)), ('xor', qs[0], (qs[1],))]
    if name == 'rzz':
        return [ph((), -rs[0] / 2), ('phx', tuple(qs), rs[0])]
    if name in ('ccx', 'toffoli', 'c3x', 'c4x', 'mcx', 'mct'):
        return [('xor', qs[-1], tuple(qs[:-1]))]
    if name in ('ccz', 'mcz'):
        return [ph(qs, Fraction(1))]
    if name in ('cswap', 'fredkin'):
        return [('cswap', qs[0], qs[1], qs[2])]
    raise Unsupported('unsupported gate: %s' % name)


def invert_prims(prims):
    """Inverse of a primitive list: reversed, every phase negated (H, XOR, swap, cswap are self-inverse)."""
    out = []
    for p in reversed(prims):
        if p[0] in ('ph', 'phx'):
            out.append((p[0], p[1], -p[2]))
        else:
            out.append(p)
    return out


def needed_K(prims, at_least=3):
    """Smallest K such that every phase pi*r of the primitives is a multiple of 2*pi/2^K: K-1 >= log2(denominator)."""
    K = at_least
    for p in prims:
        if p[0] in ('ph', 'phx'):
            d = _mod2(p[2]).denominator
            if d & (d - 1):
                raise NotDyadic('non-dyadic phase pi*%s' % p[2])
            K = max(K, d.bit_length())
    return K


class DCircuit:
    """A unitary circuit on n qubits with exact dyadic angles.
        DCircuit(n, ops=[(name, qubits, *angles), ...], global_phase=0)
        c.add('cp', (0, 1), 'pi/4');  c.add('h', (2,));  c.add('mcp', (0, 1, 2), Pi(1, 8))
    Angles: see the module docstring (Pi, str, float radians; kept modulo 4*pi).  global_phase: a circuit-level
    factor exp(i*g).  float_abs_tol: tolerance (radians) for float angles from a source that rounds them."""

    def __init__(self, n, ops=None, global_phase=0, name=None, float_abs_tol=None):
        self.n = int(n)
        if self.n < 0:
            raise ValueError('negative number of qubits')
        self.ops = []
        self.angle_rec = {}
        self.float_abs_tol = float_abs_tol
        self.gphase = Fraction(0) if not global_phase else _mod2(to_angle(global_phase, self.angle_rec,
                                                                           float_abs_tol))
        self.name = name
        for op in ops or ():
            self.add(op[0], op[1], *op[2:])

    def add(self, name, qubits, *angles):
        nm = str(name).lower()
        if nm not in GATES:
            raise Unsupported('unsupported gate: %s' % name)
        npar, arity = GATES[nm]
        qs = tuple(int(q) for q in qubits)
        if len(angles) != npar:
            raise ValueError('%s takes %d angle(s), got %d' % (nm, npar, len(angles)))
        if arity is None:
            if not qs:
                raise ValueError('%s needs at least one qubit' % nm)
        elif len(qs) != arity:
            raise ValueError('%s acts on %d qubit(s), got %d' % (nm, arity, len(qs)))
        if any(q < 0 or q >= self.n for q in qs) or len(set(qs)) != len(qs):
            raise ValueError('bad qubit operands for %s: %s (n=%d)' % (nm, list(qs), self.n))
        rs = tuple(to_angle(a, self.angle_rec, self.float_abs_tol) for a in angles)
        self.ops.append((nm, qs, rs))
        return self

    def add_global_phase(self, angle):
        self.gphase = _mod2(self.gphase + to_angle(angle, self.angle_rec, self.float_abs_tol))
        return self

    def __len__(self):
        return len(self.ops)

    def prims(self):
        out = []
        for nm, qs, rs in self.ops:
            out += expand(nm, qs, rs)
        if self.gphase:
            out.append(('ph', (), self.gphase))
        return out

    def copy(self):
        c = DCircuit(self.n, name=self.name, float_abs_tol=self.float_abs_tol)
        c.ops = list(self.ops)
        c.gphase = self.gphase
        c.angle_rec = dict(self.angle_rec)
        return c

    def to_qasm(self):
        """OpenQASM 3 text with exact angles written as pi*p/q."""
        def ang(r):
            r = _mod4(r)
            return '0' if not r else 'pi*%d/%d' % (r.numerator, r.denominator)
        lines = ['OPENQASM 3.0;', 'include "stdgates.inc";', 'qubit[%d] q;' % self.n]
        if self.gphase:
            lines.append('gphase(%s);' % ang(self.gphase))
        for nm, qs, rs in self.ops:
            par = '(%s)' % ', '.join(ang(r) for r in rs) if rs else ''
            lines.append('%s%s %s;' % (nm, par, ', '.join('q[%d]' % q for q in qs)))
        return '\n'.join(lines) + '\n'


# ----------------------------------------------------------------------------- OpenQASM 2 / 3 subset reader
_IGNORE = ('OPENQASM', 'include', 'creg', 'bit', 'barrier')
_REFUSE = ('measure', 'reset', 'if', 'while', 'for', 'def', 'opaque', 'let', 'box', 'delay', 'switch', 'return',
           'break', 'continue', 'extern', 'defcal', 'cal', 'input', 'output', 'const', 'float', 'int', 'uint',
           'angle', 'bool', 'duration', 'stretch', 'array', 'complex', 'gate_call')


def _strip_comments(s):
    s = re.sub(r'/\*.*?\*/', ' ', s, flags=re.S)
    return re.sub(r'//[^\n]*', '', s)


def _statements(s):
    """Top-level statements: 'gate ... { ... }' blocks and ';'-terminated statements."""
    out, i, n = [], 0, len(s)
    while i < n:
        while i < n and s[i] in ' \t\r\n;':
            i += 1
        if i >= n:
            break
        if re.match(r'gate\s', s[i:]):
            j = s.find('{', i)
            if j < 0:
                raise ValueError('gate definition without a body')
            depth, k = 0, j
            while k < n:
                if s[k] == '{':
                    depth += 1
                elif s[k] == '}':
                    depth -= 1
                    if depth == 0:
                        break
                k += 1
            if k >= n:
                raise ValueError('unbalanced { } in a gate definition')
            out.append(s[i:k + 1])
            i = k + 1
            continue
        if s[i] == '{' or re.match(r'(if|for|while|def|box|switch)\b', s[i:]):
            raise Unsupported('unsupported statement: %s' % s[i:i + 60])
        j = s.find(';', i)
        if j < 0:
            raise ValueError('missing ; after: %s' % s[i:i + 60])
        out.append(s[i:j].strip())
        i = j + 1
    return out


def _split_args(s):
    return [a.strip() for a in s.split(',') if a.strip()]


def _split_params(s):
    out, depth, cur = [], 0, ''
    for ch in s:
        if ch == ',' and depth == 0:
            out.append(cur.strip())
            cur = ''
            continue
        depth += ch == '('
        depth -= ch == ')'
        cur += ch
    if cur.strip():
        out.append(cur.strip())
    return out


_CALL = re.compile(r'^(?P<mods>(?:\s*(?:ctrl|negctrl|inv|pow)\s*(?:\(\s*[^)]*\))?\s*@)*)\s*'
                   r'(?P<name>[A-Za-z_][A-Za-z_0-9]*)\s*(?:\((?P<par>.*)\))?\s*(?P<args>[^()]*)$', re.S)


def parse_qasm(text):
    """OpenQASM 2.0 / 3 subset (text or a path ending in .qasm) -> DCircuit.  Supported: qreg/qubit declarations
    (registers are concatenated in order), gate calls with angle expressions, `gate` definitions (expanded), gphase,
    the modifiers ctrl / ctrl(k) / negctrl on x, z and phase gates; barrier, creg/bit and include are ignored.
    measure, reset, if, loops, inv/pow modifiers, opaque gates and other statements raise Unsupported."""
    if '\n' not in text and text.strip().endswith(('.qasm', '.qasm3', '.txt')):
        with open(text.strip()) as f:
            text = f.read()
    regs, N = {}, 0
    defs = {}
    ops, gph = [], []
    rec = {}
    for st in _statements(_strip_comments(text)):
        hw = re.match(r'[A-Za-z_]+', st)
        hw = hw.group(0) if hw else st[:1]
        if hw in _IGNORE:
            continue
        if hw in _REFUSE or re.search(r'\bmeasure\b|\breset\b', st) or ('=' in st and not st.startswith('gate')):
            raise Unsupported('unsupported statement: %s' % st[:80])
        m = (re.match(r'qreg\s+([A-Za-z_]\w*)\s*\[\s*(\d+)\s*\]$', st) or
             re.match(r'qubit\s*\[\s*(\d+)\s*\]\s*([A-Za-z_]\w*)$', st) or re.match(r'qubit\s+([A-Za-z_]\w*)$', st))
        if m:
            if st.startswith('qreg'):
                nm, sz = m.group(1), int(m.group(2))
            elif '[' in st:
                nm, sz = m.group(2), int(m.group(1))
            else:
                nm, sz = m.group(1), 1
            if nm in regs:
                raise ValueError('register %s declared twice' % nm)
            regs[nm] = (N, sz)
            N += sz
            continue
        if hw == 'qreg' or hw == 'qubit':
            raise ValueError('cannot read declaration: %s' % st[:80])
        if st.startswith('gate'):
            g = re.match(r'gate\s+([A-Za-z_]\w*)\s*(?:\(([^)]*)\))?\s*([^{]*)\{(.*)\}$', st, re.S)
            if not g:
                raise ValueError('cannot read gate definition: %s' % st[:80])
            defs[g.group(1)] = (_split_args(g.group(2) or ''), _split_args(g.group(3)),
                                [b for b in _statements(g.group(4)) if b])
            continue
        _call(st, None, None, regs, defs, ops, gph, rec, 0)
    c = DCircuit(N)
    c.angle_rec = rec
    for nm, qs, r in ops:
        c.add(nm, qs, *[Pi(x.numerator, x.denominator) for x in r])
    for r in gph:
        c.gphase = _mod2(c.gphase + r)
    return c


def _qarg(a, regs, qmap):
    a = a.strip()
    if qmap is not None:
        if a not in qmap:
            raise ValueError('unknown qubit argument %s in a gate body' % a)
        return [qmap[a]]
    m = re.match(r'([A-Za-z_]\w*)\s*\[\s*(\d+)\s*\]$', a)
    if m:
        if m.group(1) not in regs:
            raise ValueError('unknown register in %s' % a)
        base, sz = regs[m.group(1)]
        k = int(m.group(2))
        if k >= sz:
            raise ValueError('index out of range: %s' % a)
        return [base + k]
    if a in regs:
        base, sz = regs[a]
        return list(range(base, base + sz))
    raise ValueError('unknown qubit argument: %s' % a)


def _call(st, env, qmap, regs, defs, ops, gph, rec, depth):
    if depth > 64:
        raise ValueError('gate definitions nested too deep')
    st = st.strip()
    if not st or re.match(r'barrier\b', st):
        return
    m = _CALL.match(st)
    if not m:
        raise ValueError('cannot read statement: %s' % st[:80])
    mods, name = m.group('mods'), m.group('name')
    pars = [eval_expr(p, env) for p in _split_params(m.group('par'))] if m.group('par') is not None else []
    qargs = [_qarg(a, regs, qmap) for a in _split_args(m.group('args'))]
    if any(len(a) != 1 for a in qargs):
        raise Unsupported('register broadcast is not supported: %s' % st[:80])
    qs = [a[0] for a in qargs]
    if name == 'gphase' and not mods.strip():
        if len(pars) != 1 or qs:
            raise ValueError('gphase takes one angle and no qubits')
        gph.append(val_angle(pars[0], rec))
        return
    if mods.strip():
        neg, pos = [], 0
        for mm in re.finditer(r'(ctrl|negctrl|inv|pow)\s*(?:\(\s*([^)]*)\))?\s*@', mods):
            kind, arg = mm.group(1), mm.group(2)
            if kind in ('inv', 'pow'):
                raise Unsupported('modifier %s is not supported' % kind)
            k = int(arg) if arg else 1
            if kind == 'negctrl':
                neg += list(range(pos, pos + k))
            pos += k
        base = {'x': 'mcx', 'z': 'mcz', 'p': 'mcp', 'phase': 'mcp', 'u1': 'mcp'}.get(name)
        if base is None or len(qs) != pos + 1:
            raise Unsupported('ctrl modifier on %s is not supported' % name)
        for i in neg:
            ops.append(('x', (qs[i],), ()))
        ops.append((base, tuple(qs), tuple(val_angle(p, rec) for p in pars)))
        for i in neg:
            ops.append(('x', (qs[i],), ()))
        return
    if name in defs:
        fpar, fq, body = defs[name]
        if len(fpar) != len(pars) or len(fq) != len(qs):
            raise ValueError('wrong arity in a call of %s' % name)
        env2 = dict(zip(fpar, pars))
        qmap2 = dict(zip(fq, qs))
        for b in body:
            _call(b, env2, qmap2, regs, defs, ops, gph, rec, depth + 1)
        return
    nm = {'U': 'u', 'CX': 'cx'}.get(name, name.lower())
    if nm not in GATES:
        raise Unsupported('unsupported gate: %s' % name)
    ops.append((nm, tuple(qs), tuple(val_angle(p, rec) for p in pars)))


def from_qc(path):
    """qoptct .qc reference (read with qoptct.verify._parse.parse_qc) -> (DCircuit, names, inputs)."""
    from ..verify._parse import parse_qc
    names, ins, gates = parse_qc(path)
    c = DCircuit(len(names))
    for nm, qs in gates:
        c.add(nm, qs)
    return c, names, ins
