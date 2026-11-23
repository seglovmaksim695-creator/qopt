"""Canaries of the E2 glue: corrupt the formula of ONE step contract (as if the step had been proved with another
formula).  Each one must give NOT proved.  Applied by replacing the step function in zspecs.ZSTEPS for the run."""
from . import zspecs as Z


def _wrap(k, post_fix=None, pre_fix=None, calls_fix=None):
    orig = Z.ZSTEPS[k]
    def f(o, W, O, patched):
        pre, calls, post = orig(o, W, O, patched)
        if post_fix:
            post = dict(post, **post_fix(o, W, O, post))
        if pre_fix:
            pre = pre_fix(o, W, O, pre)
        if calls_fix:
            calls = calls_fix(o, W, O, calls)
        return pre, calls, post
    return f


MUTANTS = {
    # шаг 1: ctrl без условия на f4 (Q = O)
    's1_ctrl_no_f4': (1, lambda o, W, O, post: dict(ctrl=o.and_(o.not_(post['f3']), o.not_(o.and_(post['f1'], o.eq(W['b'], o.modneg(W['y']))))))),
    # шаг 1: f1 = [b == y] вместо [a == x]
    's1_f1_wrong': (1, lambda o, W, O, post: dict(f1=o.eq(W['b'], W['y']))),
    # шаг 2: x2 = x + a вместо x - a
    's2_x_plus': (2, lambda o, W, O, post: dict(x=o.modadd(W['x'], W['a']))),
    # шаг 2: lam = Z даже при f1 (без lam_r для удвоения)
    's2_lam_no_lr': (2, lambda o, W, O, post: dict(lam=o.ite(W['ctrl'], O[1], o.zero()))),
    # шаг 2: f1 не снимается
    's2_f1_keep': (2, lambda o, W, O, post: dict(f1=W['f1'])),
    # шаг 3: x3 = x + 2a вместо 3a
    's3_2a': (3, lambda o, W, O, post: dict(x=o.ite(W['ctrl'], o.modadd(W['x'], o.moddbl(W['a'])), W['x']))),
    # шаг 3: y3 не обнуляется (вычитание пропущено)
    's3_no_sub': (3, lambda o, W, O, post: dict(y=W['y'])),
    # шаг 4: x4 = x + Z3 вместо x - Z3
    's4_plus': (4, lambda o, W, O, post: dict(x=o.modadd(W['x'], O[0]))),
    # шаг 5: x5 = a - x без условия ctrl
    's5_neg_always': (5, lambda o, W, O, post: dict(x=o.modadd(W['a'], o.modneg(W['x'])))),
    # шаг 5: y5 не вычитает b
    's5_no_b': (5, lambda o, W, O, post: dict(y=W['y'])),
    # шаг 6: результат без вычитания a для P = -Q (x2 = x1)
    's6_no_cmodsub': (6, lambda o, W, O, post: dict(x=W['x'])),
    # шаг 6: y не меняется вовсе
    's6_y_keep': (6, lambda o, W, O, post: dict(y=W['y'])),
}


def apply(name):
    k, fix = MUTANTS[name]
    saved = Z.ZSTEPS[k]
    Z.ZSTEPS[k] = _wrap(k, post_fix=fix)
    return lambda: Z.ZSTEPS.__setitem__(k, saved)
