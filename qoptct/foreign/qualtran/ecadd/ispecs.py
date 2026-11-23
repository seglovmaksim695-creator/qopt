"""The same step specifications on integers (a second, independent record; for simulation and for checking the glue).
v: dict name -> int; o: list of outputs of the forward boxes (dict).  They return the same as specs.py."""


def _neg(y, p):
    return 0 if y == 0 else p - y


class S1:
    def __init__(s, n, p, patched=True): s.n, s.p, s.patched = n, p, patched
    def pre(s, v, o): return all(v[k] < s.p for k in 'abxy')
    def calls(s, v, o): return []
    def post(s, v, o):
        a, b, x, y, p = v['a'], v['b'], v['x'], v['y'], s.p
        e1, e2, e3, e4 = int(a == x), int(b == _neg(y, p)), int(a == 0 and b == 0), int(x == 0 and y == 0)
        if s.patched:
            ctrl = int(not e3 and not e4 and not (e1 and e2)); f2 = e2 & (1 - ctrl)
        else:
            ctrl = int(not e2 and not e3 and not e4); f2 = e2
        return dict(f1=e1, f2=f2, f3=e3, f4=e4, ctrl=ctrl, a=a, b=b, x=x, y=y)


class S2:
    def __init__(s, n, p, patched=True): s.n, s.p, s.patched = n, p, patched
    def pre(s, v, o): return all(v[k] < s.p for k in ('a', 'b', 'x', 'y', 'lam_r'))
    def _xy(s, v):
        p = s.p
        return (v['x'] - v['a']) % p, ((v['y'] - v['b']) % p if v['ctrl'] else v['y'])
    def calls(s, v, o):
        x2, y2 = s._xy(v)
        return [('KINV', {'x': x2}), ('MM', {'x': o[0]['x'], 'y': y2})]
    def post(s, v, o):
        x2, y2 = s._xy(v)
        Z = o[1]['target']; f1, ctrl, lr = v['f1'], v['ctrl'], v['lam_r']
        lam = (lr if f1 else Z) if ctrl else 0
        if s.patched:
            f1o = f1 ^ (ctrl & int(x2 == 0))
        else:
            f1o = f1 ^ (int(Z != lr) & int(lam == lr))
        return dict(f1=f1o, ctrl=ctrl, a=v['a'], b=v['b'], x=x2, y=y2, lam=lam, lam_r=lr)


class S3:
    def __init__(s, n, p, patched=True): s.n, s.p, s.patched = n, p, patched
    def pre(s, v, o): return all(v[k] < s.p for k in ('a', 'b', 'x', 'y', 'lam'))
    def calls(s, v, o): return [('MM', {'x': v['x'], 'y': v['lam']})]
    def post(s, v, o):
        p, c = s.p, v['ctrl']; Z1 = o[0]['target']
        return dict(ctrl=c, a=v['a'], b=v['b'], x=(v['x'] + 3 * v['a']) % p if c else v['x'],
                    y=(v['y'] - Z1) % p if c else v['y'], lam=v['lam'])


class S4:
    def __init__(s, n, p, patched=True): s.n, s.p, s.patched = n, p, patched
    def pre(s, v, o): return v['x'] < s.p and v['lam'] < s.p
    def calls(s, v, o):
        x4 = (v['x'] - o[0]['target']) % s.p
        return [('MM', {'x': v['lam'], 'y': v['lam']}), ('MM', {'x': x4, 'y': v['lam']})]
    def post(s, v, o):
        return dict(x=(v['x'] - o[0]['target']) % s.p, y=v['y'] ^ o[1]['target'], lam=v['lam'])


class S5:
    def __init__(s, n, p, patched=True): s.n, s.p, s.patched = n, p, patched
    def pre(s, v, o):
        K, Z = o[0]['x'], o[1]['target']
        p5 = v['lam'] == (v['ctrl'] * (Z ^ (v['lam_r'] if K == 0 else 0)))
        return all(v[k] < s.p for k in ('a', 'b', 'x', 'y', 'lam_r')) and p5
    def calls(s, v, o): return [('KINV', {'x': v['x']}), ('MM', {'x': o[0]['x'], 'y': v['y']})]
    def post(s, v, o):
        p, c = s.p, v['ctrl']
        xn = _neg(v['x'], p) if c else v['x']
        return dict(ctrl=c, a=v['a'], b=v['b'], x=(v['a'] + xn) % p, y=(v['y'] - v['b']) % p if c else v['y'],
                    lam_r=v['lam_r'])


class S6:
    def __init__(s, n, p, patched=True): s.n, s.p, s.patched = n, p, patched
    def literal(s, v):
        p = s.p
        f1, f2, f3, f4, ctrl = (v[k] for k in ('f1', 'f2', 'f3', 'f4', 'ctrl'))
        a, b, x, y = v['a'], v['b'], v['x'], v['y']
        za, zb, zx, zy = (int(t == 0) for t in (a, b, x, y))
        if s.patched:
            nn = (1 - f3) & (1 - f4); ctrl1 = ctrl ^ nn ^ (f1 & f2 & nn)
        else:
            ctrl1 = ctrl ^ ((1 - f2) & (1 - f3) & (1 - f4))
        f1a = f1 ^ (za & zx & zy) ^ (zx & za & zb)
        f2a = f2 ^ (za & zb & zy) ^ (zx & zy & zb)
        x1, y1 = x ^ (a if f4 else 0), y ^ (b if f4 else 0)
        anc = f1a & f2a
        x2 = (x1 - a) % p if anc else x1
        y2 = (y1 + b) % p if anc else y1
        f4a = f4 ^ int(a == x2 and b == y2)
        f3a = f3 ^ (za & zb)
        m = int(x2 == 0 and y2 == 0)
        return dict(x1=x1, y1=y1, x2=x2, y2=y2, flags=(ctrl1, f1a ^ m, f2a ^ m, f3a, f4a))
    def pre(s, v, o):
        L = s.literal(v)
        return (v['a'] < s.p and v['b'] < s.p and L['x1'] < s.p and L['y1'] < s.p and not any(L['flags']))
    def calls(s, v, o): return []
    def post(s, v, o):
        L = s.literal(v)
        return dict(a=v['a'], b=v['b'], x=L['x2'], y=L['y2'])


ISTEPS = {1: S1, 2: S2, 3: S3, 4: S4, 5: S5, 6: S6}
