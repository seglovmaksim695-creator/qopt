"""Спецификации шагов в третьей записи: один код для двух «алгебр» —
IntOps (целые Python: прогон и сверка с specs/ispecs) и SmtOps (термы SMT-LIB для склейки в Z3).
Структура формул — построчно как в specs.py. XOR в SmtOps заменён свежей переменной с верными свойствами XOR
(u=0 -> v=w; w=0 -> v=u; u=w -> v=0; 0<=v<2^n) — это огрубление: доказанное с ним верно и для настоящего XOR."""


class IntOps:
    def __init__(self, n, p):
        self.n, self.p = n, p
    def lt(self, x, c): return x < c
    def eq(self, x, y): return x == y
    def eqz(self, x): return x == 0
    def and_(self, *bs): return all(bs)
    def or_(self, *bs): return any(bs)
    def not_(self, b): return not b
    def xorb(self, a, b): return bool(a) != bool(b)
    def ite(self, c, t, e): return t if c else e
    def modadd(self, x, y): return (x + y) % self.p
    def modsub(self, x, y): return (x - y) % self.p
    def modneg(self, y): return 0 if y == 0 else self.p - y
    def moddbl(self, x): return (2 * x) % self.p
    def xor(self, u, w): return u ^ w
    def zero(self): return 0


class SmtOps:
    def __init__(self, n, p):
        self.n, self.p, self.lines, self.k, self.sorts = n, p, [], 0, {}
    def _new(self, sort, term=None):
        self.k += 1
        v = f'v{self.k}'
        self.sorts[v] = sort
        self.lines.append(f'(declare-const {v} {sort})')
        if term is not None:
            self.lines.append(f'(assert (= {v} {term}))')
        return v
    def word(self, name):
        self.sorts[name] = 'Int'
        self.lines.append(f'(declare-const {name} Int)')
        self.lines.append(f'(assert (and (<= 0 {name}) (< {name} {self.p})))')
        return name
    def lt(self, x, c): return self._new('Bool', f'(< {x} {c})')
    def eq(self, x, y): return self._new('Bool', f'(= {x} {y})')
    def eqz(self, x): return self._new('Bool', f'(= {x} 0)')
    def and_(self, *bs): return self._new('Bool', '(and ' + ' '.join(map(str, bs)) + ')') if len(bs) > 1 else bs[0]
    def or_(self, *bs): return self._new('Bool', '(or ' + ' '.join(map(str, bs)) + ')') if len(bs) > 1 else bs[0]
    def not_(self, b): return self._new('Bool', f'(not {b})')
    def xorb(self, a, b): return self._new('Bool', f'(xor {a} {b})')
    def ite(self, c, t, e):
        st = self.sorts.get(t, 'Int') if isinstance(t, str) else 'Int'
        se = self.sorts.get(e, 'Int') if isinstance(e, str) else 'Int'
        assert st == se == 'Int', ('ite только над словами', t, e)
        return self._new('Int', f'(ite {c} {t} {e})')
    def modadd(self, x, y):
        s = self._new('Int', f'(+ {x} {y})')
        return self._new('Int', f'(ite (>= {s} {self.p}) (- {s} {self.p}) {s})')
    def modsub(self, x, y):
        return self._new('Int', f'(ite (>= {x} {y}) (- {x} {y}) (+ (- {x} {y}) {self.p}))')
    def modneg(self, y): return self._new('Int', f'(ite (= {y} 0) 0 (- {self.p} {y}))')
    def moddbl(self, x): return self.modadd(x, x)
    def xor(self, u, w):
        v = self._new('Int')
        self.lines.append(f'(assert (and (<= 0 {v}) (< {v} {1 << self.n})))')
        self.lines.append(f'(assert (=> (= {u} 0) (= {v} {w})))')
        self.lines.append(f'(assert (=> (= {w} 0) (= {v} {u})))')
        self.lines.append(f'(assert (=> (= {u} {w}) (= {v} 0)))')
        return v
    def zero(self): return 0


def rng(o, *ws):
    return [o.lt(w, o.p) for w in ws]


def step1(o, W, O, patched):
    a, b, x, y = W['a'], W['b'], W['x'], W['y']
    e1 = o.eq(a, x); e2 = o.eq(b, o.modneg(y))
    e3 = o.and_(o.eqz(a), o.eqz(b)); e4 = o.and_(o.eqz(x), o.eqz(y))
    if patched:
        ctrl = o.and_(o.not_(e3), o.not_(e4), o.not_(o.and_(e1, e2))); f2 = o.and_(e2, o.not_(ctrl))
    else:
        ctrl = o.and_(o.not_(e2), o.not_(e3), o.not_(e4)); f2 = e2
    return rng(o, a, b, x, y), [], dict(f1=e1, f2=f2, f3=e3, f4=e4, ctrl=ctrl, a=a, b=b, x=x, y=y)


def step2(o, W, O, patched):
    f1, ctrl, lr = W['f1'], W['ctrl'], W['lam_r']
    x2 = o.modsub(W['x'], W['a'])
    y2 = o.ite(ctrl, o.modsub(W['y'], W['b']), W['y'])
    calls = [('KINV', {'x': x2}), ('MM', {'x': O[0], 'y': y2})]
    Z = O[1]
    lam = o.ite(ctrl, o.ite(f1, lr, Z), o.zero())
    if patched:
        f1o = o.xorb(f1, o.and_(ctrl, o.eqz(x2)))
    else:
        anc = o.eq(Z, lr)
        f1o = o.xorb(f1, o.and_(o.not_(anc), o.eq(lam, lr)))
    return (rng(o, W['a'], W['b'], W['x'], W['y'], lr), calls,
            dict(f1=f1o, ctrl=ctrl, a=W['a'], b=W['b'], x=x2, y=y2, lam=lam, lam_r=lr))


def step3(o, W, O, patched):
    ctrl, a = W['ctrl'], W['a']
    calls = [('MM', {'x': W['x'], 'y': W['lam']})]
    Z1 = O[0]
    y3 = o.ite(ctrl, o.modsub(W['y'], Z1), W['y'])
    a3 = o.modadd(a, o.moddbl(a))
    x3 = o.ite(ctrl, o.modadd(W['x'], a3), W['x'])
    return rng(o, a, W['b'], W['x'], W['y'], W['lam']), calls, dict(ctrl=ctrl, a=a, b=W['b'], x=x3, y=y3, lam=W['lam'])


def step4(o, W, O, patched):
    x4 = o.modsub(W['x'], O[0])
    calls = [('MM', {'x': W['lam'], 'y': W['lam']}), ('MM', {'x': x4, 'y': W['lam']})]
    return rng(o, W['x'], W['lam']), calls, dict(x=x4, y=o.xor(W['y'], O[1]), lam=W['lam'])


def step5(o, W, O, patched):
    ctrl, lr = W['ctrl'], W['lam_r']
    K, Z = O[0], O[1]
    p5 = o.eq(W['lam'], o.ite(ctrl, o.xor(Z, o.ite(o.eqz(K), lr, o.zero())), o.zero()))
    calls = [('KINV', {'x': W['x']}), ('MM', {'x': K, 'y': W['y']})]
    xn = o.ite(ctrl, o.modneg(W['x']), W['x'])
    x5 = o.modadd(W['a'], xn)
    y5 = o.ite(ctrl, o.modsub(W['y'], W['b']), W['y'])
    return (rng(o, W['a'], W['b'], W['x'], W['y'], lr) + [p5], calls,
            dict(ctrl=ctrl, a=W['a'], b=W['b'], x=x5, y=y5, lam_r=lr))


def step6(o, W, O, patched):
    f1, f2, f3, f4, ctrl = (W[k] for k in ('f1', 'f2', 'f3', 'f4', 'ctrl'))
    a, b, x, y = W['a'], W['b'], W['x'], W['y']
    za, zb, zx, zy = o.eqz(a), o.eqz(b), o.eqz(x), o.eqz(y)
    if patched:
        nn = o.and_(o.not_(f3), o.not_(f4))
        ctrl1 = o.xorb(o.xorb(ctrl, nn), o.and_(o.and_(f1, f2), nn))
    else:
        ctrl1 = o.xorb(ctrl, o.and_(o.not_(f2), o.not_(f3), o.not_(f4)))
    f1a = o.xorb(o.xorb(f1, o.and_(za, zx, zy)), o.and_(zx, za, zb))
    f2a = o.xorb(o.xorb(f2, o.and_(za, zb, zy)), o.and_(zx, zy, zb))
    x1 = o.xor(x, o.ite(f4, a, o.zero()))
    y1 = o.xor(y, o.ite(f4, b, o.zero()))
    anc = o.and_(f1a, f2a)
    x2 = o.ite(anc, o.modsub(x1, a), x1)
    y2 = o.ite(anc, o.modadd(y1, b), y1)
    f4a = o.xorb(f4, o.and_(o.eq(a, x2), o.eq(b, y2)))
    f3a = o.xorb(f3, o.and_(za, zb))
    m = o.and_(o.eqz(x2), o.eqz(y2))
    f1b, f2b = o.xorb(f1a, m), o.xorb(f2a, m)
    pre = rng(o, a, b, x1, y1) + [o.not_(f) for f in (ctrl1, f1b, f2b, f3a, f4a)]
    return pre, [], dict(a=a, b=b, x=x2, y=y2)


ZSTEPS = {1: step1, 2: step2, 3: step3, 4: step4, 5: step5, 6: step6}
