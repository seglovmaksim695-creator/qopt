"""qoptct.foreign.auto: the library-independent core (expressions, contracts, composition) without Qualtran."""
import random

import pytest

from qoptct.foreign.auto import expr as X
from qoptct.foreign.auto.contract import Contract, ContractError
from qoptct.foreign.auto.tree import Node, Program

z3 = pytest.importorskip('z3')
V = X.var
W = {'a': 5, 'b': 4, 'c': 1}


def _rand_expr(rng, d):
    if d == 0 or rng.random() < 0.2:
        return ('var', rng.choice(list(W))) if rng.random() < 0.6 else rng.randrange(-3, 40)
    op = rng.choice(['add', 'sub', 'mulc', 'mulb', 'div', 'mod', 'lt', 'le', 'eq', 'ne', 'ite', 'bxor', 'band', 'bit',
                     'slice', 'shl', 'shr', 'lut', 'and', 'not'])
    a = _rand_expr(rng, d - 1)
    if op in ('add', 'sub', 'lt', 'le', 'eq', 'ne'):
        return (op, a, _rand_expr(rng, d - 1))
    if op == 'mulc':
        return ('mul', a, rng.randrange(-5, 9))
    if op == 'mulb':
        return ('mul', ('var', 'c'), a)
    if op in ('div', 'mod'):
        return (op, a, rng.choice([1, 2, 3, 4, 7, 8, 16]))
    if op == 'ite':
        return ('ite', ('lt', a, _rand_expr(rng, d - 1)), _rand_expr(rng, d - 1), _rand_expr(rng, d - 1))
    if op in ('bxor', 'band'):
        return (op, ('mod', a, 32), ('var', 'b'))
    if op == 'bit':
        return ('bit', a, rng.randrange(0, 6))
    if op == 'slice':
        return ('slice', a, rng.randrange(0, 4), rng.randrange(1, 4))
    if op in ('shl', 'shr'):
        return (op, a, rng.randrange(0, 3))
    if op == 'lut':
        return ('lut', ('mod', a, 4), (3, 1, 4, 1))
    if op == 'and':
        return ('and', ('lt', a, 3), ('var', 'c'))
    return ('not', ('eq', a, 2))


def test_expression_back_ends_agree():
    """Python, Z3 integers, static bounds and ROBDD words give the same value on random expressions"""
    from qoptct.verify.bdd import BDD, ZERO
    rng = random.Random(5)
    lb = lambda leaf: (0, (1 << W[leaf[1]]) - 1)
    checked = 0
    for _ in range(250):
        e = _rand_expr(rng, 4)
        X.validate(e)
        vs = {k: z3.Int(k) for k in W}
        T = X.Z3Terms(z3, lambda leaf: vs[leaf[1]], lb)
        B = BDD()
        bits = {k: [B.new_var()[1] for _ in range(w)] for k, w in W.items()}
        order = [(k, j) for k in W for j in range(W[k])]
        try:
            te = T.term(e)
            word = X.bdd_word(e, B, lambda leaf: X.BddWord(bits[leaf[1]] + [ZERO], 0, (1 << W[leaf[1]]) - 1))
        except X.ExprError:
            continue
        lo, hi = X.bounds(e, lb)
        for _ in range(4):
            env = {k: rng.randrange(1 << w) for k, w in W.items()}
            v = X.evaluate(e, lambda leaf: env[leaf[1]])
            assert lo <= v <= hi
            s = z3.Solver()
            s.add(*[vs[k] == env[k] for k in W])
            s.add(te != v)
            assert s.check() == z3.unsat
            asg = {B.level(bits[k][j]): (env[k] >> j) & 1 for k, j in order}
            got = 0
            for j, bt in enumerate(word.bits):
                if B.eval(bt, asg):
                    got |= 1 << j
            if got >> (len(word.bits) - 1):
                got -= 1 << len(word.bits)
            assert got == v, (e, env)
            checked += 1
    assert checked > 300


def test_contract_needs_every_output_and_checks_ranges():
    with pytest.raises(ContractError):
        Contract({'x': 3}, {'x': 3, 'junk': 1}, {'x': V('x')}).validate()
    c = Contract({'x': 3}, {'x': 3}, {'x': X.add(V('x'), 1)}).validate()
    with pytest.raises(ContractError):
        c.evaluate({'x': 7})           # 8 does not fit 3 bits: no silent truncation
    assert c.evaluate({'x': 6}) == {'x': 7}


# ------------------------------------------------------------------------------------------- synthetic programs
class _Adapter:
    def __init__(self):
        self.kids, self.contracts = {}, {}

    def children(self, node):
        return self.kids.get(node.path, [])

    def contract(self, node):
        return self.contracts.get(node.path)


def _fake_leaf(ops, ins, outs, contract, **kw):
    return dict(status='PROVED', method='fake', ops=0, wall=0)


def _prim(ins, outs, out, name='p'):
    return Contract({k: len(v) for k, v in ins.items()}, {k: len(v) for k, v in outs.items()}, out, name=name,
                    kind='primitive')


def _modadd(n, p, mutate=None):
    """the structure of Qualtran's ModAdd with fake leaves (only the composition is tested)"""
    A = _Adapter()
    cnt = [0]

    def wires(k):
        r = list(range(cnt[0], cnt[0] + k))
        cnt[0] += k
        return r
    x, y = wires(n), wires(n)
    junk, sign = wires(1)[0], wires(1)[0]
    root = Node('M', 'ModAdd', 'op', {'x': x, 'y': y}, {'x': x, 'y': y}, None)
    kids = []

    def add(typ, kind, ins, outs, c=None):
        nd = Node(f'M/{len(kids)}:{typ}', typ, kind, ins, outs, None)
        kids.append(nd)
        if c is not None:
            A.contracts[nd.path] = c
    add('Allocate', 'alloc', {}, {'reg': [junk]})
    add('Allocate', 'alloc', {}, {'reg': [sign]})
    a1, b1 = x + [junk], y + [sign]
    if mutate == 'swap':
        a1 = [x[1], x[0]] + x[2:] + [junk]
    add('Add', 'op', {'a': a1, 'b': b1}, {'a': a1, 'b': b1},
        _prim({'a': a1, 'b': b1}, {'a': a1, 'b': b1}, {'a': V('a'), 'b': X.mod(X.add(V('a'), V('b')), 1 << (n + 1))}))
    k = -p if mutate != 'wrongk' else -(p - 1)
    add('AddK', 'op', {'x': b1}, {'x': b1}, _prim({'x': b1}, {'x': b1}, {'x': X.mod(X.add(V('x'), k), 1 << (n + 1))}))
    add('Controlled', 'op', {'ctrl': [sign], 'x': y}, {'ctrl': [sign], 'x': y},
        _prim({'ctrl': [sign], 'x': y}, {'ctrl': [sign], 'x': y},
              {'ctrl': V('ctrl'), 'x': X.ite(X.eq(V('ctrl'), 1), X.mod(X.add(V('x'), p), 1 << n), V('x'))}))
    add('GT', 'op', {'a': x, 'b': y, 'target': [sign]}, {'a': x, 'b': y, 'target': [sign]},
        _prim({'a': x, 'b': y, 'target': [sign]}, {'a': x, 'b': y, 'target': [sign]},
              {'a': V('a'), 'b': V('b'), 'target': X.bxor(V('target'), X.lt(V('b'), V('a')))}))
    if mutate != 'nox':
        add('XGate', 'op', {'q': [sign]}, {'q': [sign]}, _prim({'q': [sign]}, {'q': [sign]}, {'q': X.sub(1, V('q'))}))
    add('Free', 'free', {'reg': [junk]}, {})
    add('Free', 'free', {'reg': [sign]}, {})
    A.kids[root.path] = kids
    pre = X.and_(X.lt(V('x'), p), X.lt(V('y'), p)) if mutate != 'weakpre' else X.lt(V('x'), p)
    spec = Contract({'x': n, 'y': n}, {'x': n, 'y': n}, {'x': V('x'), 'y': X.mod(X.add(V('x'), V('y')), p)},
                    pre=pre, name='ModAdd spec', kind='spec')
    return Program([], cnt[0], root, A), spec


@pytest.mark.parametrize('n,p', [(8, 251), (256, 2 ** 255 - 19)])
@pytest.mark.parametrize('mutate', [None, 'swap', 'wrongk', 'nox', 'weakpre'])
def test_composition_of_a_modular_adder(n, p, mutate):
    from qoptct.foreign.auto.compose import Engine
    P, spec = _modadd(n, p, mutate)
    E = Engine(P, identity=lambda P_, node: ('id:' + node.path, ''), leaf_prover=_fake_leaf)
    r = E.prove(P.root, spec)
    assert (r['status'] == 'PROVED') == (mutate is None), r.get('detail')
