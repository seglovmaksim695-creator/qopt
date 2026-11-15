# SPDX-License-Identifier: Apache-2.0
"""Shared parsing for qoptct.verify: the OpenQASM 3 subset written by qoptct and the .qc reference format.

Boolean polynomials over F2 are kept in algebraic normal form (ANF): a monomial is an int bit mask of
variables (0 is the constant 1) and a polynomial is a frozenset of monomials (XOR of its monomials).
"""
import os
import re
import tempfile

ONE = frozenset([0])
ZERO = frozenset()


def bits(m):
    """Indices of the set bits of m (low to high)."""
    r = []
    while m:
        lo = m & -m
        r.append(lo.bit_length() - 1)
        m ^= lo
    return r


def pmul(a, b):
    """Product of two ANF polynomials."""
    r = set()
    for i in a:
        for j in b:
            k = i | j
            if k in r:
                r.remove(k)
            else:
                r.add(k)
    return frozenset(r)


def peval(P, val):
    """Evaluate ANF polynomial P at the assignment val (bit v of val = value of variable v)."""
    s = 0
    for mo in P:
        if mo & val == mo:
            s ^= 1
    return s


# ----------------------------------------------------------------------------- OpenQASM 3 subset

def gate_of(s):
    s = s.rstrip(';').strip()
    name, args = s.split(None, 1)
    return name, tuple(int(x) for x in re.findall(r'q\[(\d+)\]', args))


def _cond_tokens(s):
    return re.findall(r'c\[\d+\]|&&|\^|\(|\)|!|&|[01]|\|\|', s)


def cond_poly(s, cvar):
    """Condition over classical bits c[k] (operators && & ^ ! and parentheses) -> ANF polynomial.
    cvar(k) returns the variable mask of bit c[k]."""
    toks = _cond_tokens(s)
    pos = [0]
    if ''.join(toks) != re.sub(r'\s+', '', s):
        raise ValueError('unsupported condition: ' + s)

    def peek():
        return toks[pos[0]] if pos[0] < len(toks) else None

    def eat(t=None):
        tok = toks[pos[0]]
        if t and tok != t:
            raise ValueError('condition: expected ' + t)
        pos[0] += 1
        return tok

    def E():
        P = A()
        while peek() == '&&':
            eat()
            P = pmul(P, A())
        return P

    def A():
        P = X()
        while peek() == '&':
            eat()
            P = pmul(P, X())
        return P

    def X():
        P = B()
        while peek() == '^':
            eat()
            P = P ^ B()
        return P

    def B():
        t = peek()
        if t == '!':
            eat()
            return ONE ^ B()
        if t == '(':
            eat()
            P = E()
            eat(')')
            return P
        if t in ('0', '1'):
            eat()
            return ONE if t == '1' else ZERO
        if t and t.startswith('c['):
            eat()
            return frozenset([cvar(int(t[2:-1]))])
        raise ValueError('unsupported condition: ' + s)

    P = E()
    if pos[0] != len(toks):
        raise ValueError('condition: trailing tokens in ' + s)
    return P


def eval_cond(s, cbits):
    """Evaluate a condition string for concrete classical bits (dict k -> 0/1; missing bits read 0)."""
    P = cond_poly(s, lambda k: 1 << k)
    return peval(P, sum(1 << k for k, b in cbits.items() if b))


def _balanced(s, i, o, c):
    """Index just after the bracket that closes s[i] == o."""
    depth = 0
    for j in range(i, len(s)):
        if s[j] == o:
            depth += 1
        elif s[j] == c:
            depth -= 1
            if depth == 0:
                return j + 1
    raise ValueError('unbalanced %s%s in: %s' % (o, c, s))


def _if_parts(s):
    """s starts with 'if (' -> (condition, inner text of the block, text after the block)."""
    i = s.index('(')
    j = _balanced(s, i, '(', ')')
    t = s[j:].lstrip()
    if not t.startswith('{'):
        raise ValueError('if without a { } block: ' + s)
    k = _balanced(t, 0, '{', '}')
    return s[i + 1:j - 1].strip(), t[1:k - 1], t[k:]


_IF_START = re.compile(r'if\s*\(')


def _if_block(block, cond):
    out, gates, s = [], [], block
    while True:
        s = s.lstrip()
        if not s:
            break
        if _IF_START.match(s):
            if gates:
                out.append((cond, gates))
                gates = []
            c2, b2, s = _if_parts(s)
            out += _if_block(b2, '(%s) && (%s)' % (cond, c2))
            continue
        k = s.find(';')
        if k < 0:
            raise ValueError('missing ; in: ' + s)
        g, s = s[:k].strip(), s[k + 1:]
        if g:
            gates.append(g)
    if gates:
        out.append((cond, gates))
    return out


def split_if(L):
    """'if (C) { ... }' (also with nested ifs, which Qiskit's OpenQASM 3 importer reads, unlike && and ^) ->
    [(condition, [gate statements])] in program order; a nested if gets the condition (outer) && (inner).
    None if L is not an if statement; ValueError if it is malformed (or has an else branch)."""
    if not _IF_START.match(L):
        return None
    cond, block, rest = _if_parts(L)
    if rest.strip():
        raise ValueError('unsupported statement after an if block: ' + L)
    return _if_block(block, cond)


_STMT_GATE = re.compile(r'[a-z][a-z0-9_]*\s+q\[\d+\](\s*,\s*q\[\d+\])*')


def parse_qasm_text(txt):
    """-> (N qubits, number of classical bits, ops).
    ops: ('G', name, qubits) | ('M', qubit, bit) | ('IF', condition string, [(name, qubits), ...])."""
    m = re.search(r'qubit\[(\d+)\]', txt)
    if not m:
        raise ValueError('no qubit register')
    N = int(m.group(1))
    mb = re.search(r'(?<!qu)bit\[(\d+)\]', txt)
    nb = int(mb.group(1)) if mb else 0
    ops = []
    for L in txt.splitlines():
        L = L.split('//')[0].strip()
        if not L or L.startswith(('OPENQASM', 'include', 'qubit', 'bit', 'barrier')):
            continue
        sp = split_if(L)
        if sp is not None:
            for cond, body in sp:
                for g in body:
                    if not _STMT_GATE.fullmatch(g):
                        raise ValueError('malformed statement in an if block: ' + g)
                ops.append(('IF', cond, [gate_of(g) for g in body]))
            continue
        for g in L.split(';'):
            g = g.strip()
            if not g:
                continue
            mm = re.fullmatch(r'c\[(\d+)\]\s*=\s*measure\s+q\[(\d+)\]', g)
            if mm:
                ops.append(('M', int(mm.group(2)), int(mm.group(1))))
                continue
            if '=' in g or g.startswith(('measure', 'if')):
                raise ValueError('unsupported statement: ' + g)
            if not _STMT_GATE.fullmatch(g):
                raise ValueError('malformed statement: ' + g)
            ops.append(('G',) + gate_of(g))
    _check_ops(N, nb, ops)
    return N, nb, ops


def _check_ops(N, nb, ops):
    """Static checks of a parsed circuit: qubit and bit indices in range, conditions only on measured bits."""
    measured = set()
    for op in ops:
        if op[0] == 'M':
            if op[1] >= N or op[2] >= nb:
                raise ValueError('measure out of range: c[%d] = measure q[%d]' % (op[2], op[1]))
            measured.add(op[2])
            continue
        gs = [op[1:]] if op[0] == 'G' else op[2]
        for name, qs in gs:
            if not qs or max(qs) >= N or len(set(qs)) != len(qs):
                raise ValueError('bad qubit operands: %s %s' % (name, list(qs)))
        if op[0] == 'IF':
            used = [int(k) for k in re.findall(r'c\[(\d+)\]', op[1])]
            cond_poly(op[1], lambda k: 1 << k)  # syntax of the condition
            if any(k not in measured for k in used):
                raise ValueError('condition on an unmeasured bit: ' + op[1])


def parse_qasm(path):
    """-> (N, ops) (kept for the path-sum engine)."""
    N, _, ops = parse_qasm_text(open(path).read())
    return N, ops


def nbits_of(path):
    mb = re.search(r'(?<!qu)bit\[(\d+)\]', open(path).read())
    return int(mb.group(1)) if mb else 0


def is_unitary_program(ops):
    return all(op[0] == 'G' and op[1] != 'reset' for op in ops)


# ----------------------------------------------------------------------------- .qc reference

QC1 = {'H': 'h', 'X': 'x', 'Y': 'y', 'T': 't', 'T*': 'tdg', 'S': 's', 'S*': 'sdg', 'P': 's', 'P*': 'sdg', 'Z': 'z'}


def parse_qc(path):
    """.qc reader (tof with any number of controls, Z on 1..k wires).
    -> (names, inputs (list from the .i line or None), gates [(name, qubits)])."""
    names = None
    ins = None
    gates = []
    idx = {}
    for line in open(path):
        p = line.split()
        if not p or p[0].startswith('#'):
            continue
        if p[0] == '.v':
            names = p[1:]
            idx = {q: i for i, q in enumerate(names)}
            continue
        if p[0] == '.i':
            ins = p[1:]
            continue
        if p[0].startswith('.') or p[0] in ('BEGIN', 'END'):
            continue
        if names is None:
            raise ValueError('.qc without a .v line')
        op, qs = p[0], [idx[q] for q in p[1:]]
        if op in ('Z', 'Zd'):
            qs = list(dict.fromkeys(qs))
            gates.append(('z' if len(qs) == 1 else 'cz' if len(qs) == 2 else 'mcz', tuple(qs)))
        elif op == 'tof':
            if len(set(qs)) != len(qs):
                raise ValueError('tof with a repeated wire: ' + line.strip())
            gates.append(('x' if len(qs) == 1 else 'cx' if len(qs) == 2 else 'mcx', tuple(qs)))
        elif op in QC1 and len(qs) == 1:
            gates.append((QC1[op], tuple(qs)))
        else:
            raise ValueError('unsupported in .qc: ' + line.strip())
    if names is None:
        raise ValueError('.qc without a .v line')
    return names, ins, gates


def clean_wires(names, ins, clean):
    """Reference wires that start in |0> and must return to |0> (clean model): wires not listed on the .i line."""
    if not clean or ins is None:
        return set()
    return {i for i, nm in enumerate(names) if nm not in ins}


# ----------------------------------------------------------------------------- exact H-frame folding

DIAG = {'z', 'cz', 'ccz', 'mcz'}
TOX = {'z': 'x', 'cz': 'cx', 'ccz': 'ccx', 'mcz': 'mcx'}


def fold_h(items):
    """Exact identity H_t D H_t = D with target t in the X basis: `h t; {z/cz/ccz on t}; h t` -> x/cx/ccx on t.
    items: list of (name, qubits) or ('#', qubits) for opaque operations (measurements, if). A frame is folded only
    if every operation on t between the two h is diagonal and its other wires are not inside a candidate frame."""
    cand = []
    open_ = {}
    for i, it in enumerate(items):
        nm, qs = it[0], it[1]
        for q in qs:
            if nm == 'h' and len(qs) == 1:
                if q in open_:
                    cand.append((q, open_.pop(q), i))
                else:
                    open_[q] = i
            elif nm in DIAG:
                pass
            else:
                open_.pop(q, None)
    iv = {}
    for t, a, b in cand:
        iv.setdefault(t, []).append((a, b))

    def inside(q, i):
        return any(a < i < b for a, b in iv.get(q, ()))

    ok = []
    for t, a, b in cand:
        good = True
        for i in range(a + 1, b):
            nm, qs = items[i][0], items[i][1]
            if t in qs:
                if nm not in DIAG or any(inside(o, i) for o in qs if o != t):
                    good = False
                    break
        if good:
            ok.append((t, a, b))
    drop = set()
    conv = {}
    for t, a, b in ok:
        drop |= {a, b}
        for i in range(a + 1, b):
            if t in items[i][1]:
                conv[i] = t
    out = []
    for i, it in enumerate(items):
        nm, qs = it[0], it[1]
        if i in drop:
            continue
        if i in conv:
            t = conv[i]
            out.append((TOX[nm], tuple([q for q in qs if q != t] + [t])))
        else:
            out.append(it)
    return out


def fold_ops(ops):
    items = []
    for op in ops:
        if op[0] == 'G' and op[1] != 'reset':
            items.append((op[1], op[2], op))
        elif op[0] == 'G':
            items.append(('#', op[2], op))
        elif op[0] == 'M':
            items.append(('#', (op[1],), op))
        else:
            items.append(('#', tuple(sorted({q for _, qs in op[2] for q in qs})), op))
    out = []
    for it in fold_h(items):
        out.append(it[2] if len(it) == 3 else ('G', it[0], it[1]))
    return out


# ----------------------------------------------------------------------------- input helper

class QasmFile:
    """Context manager: accepts a path or OpenQASM text and yields a file path."""

    def __init__(self, qasm_or_path):
        self.src = qasm_or_path
        self.tmp = None

    def __enter__(self):
        s = self.src if isinstance(self.src, str) else os.fspath(self.src)
        if '\n' in s or s.lstrip().startswith('OPENQASM'):
            fd, self.tmp = tempfile.mkstemp(suffix='.qasm', prefix='qoptct_verify_')
            with os.fdopen(fd, 'w') as f:
                f.write(s)
            return self.tmp
        if not os.path.exists(s):
            raise FileNotFoundError(s)
        return s

    def __exit__(self, *exc):
        if self.tmp:
            try:
                os.unlink(self.tmp)
            except OSError:
                pass
        return False
