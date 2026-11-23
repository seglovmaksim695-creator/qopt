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
    # step 1: ctrl without the condition on f4 (Q = O)
    's1_ctrl_no_f4': (1, lambda o, W, O, post: dict(ctrl=o.and_(o.not_(post['f3']), o.not_(o.and_(post['f1'], o.eq(W['b'], o.modneg(W['y']))))))),
    # step 1: f1 = [b == y] instead of [a == x]
    's1_f1_wrong': (1, lambda o, W, O, post: dict(f1=o.eq(W['b'], W['y']))),
    # step 2: x2 = x + a instead of x - a
    's2_x_plus': (2, lambda o, W, O, post: dict(x=o.modadd(W['x'], W['a']))),
    # step 2: lam = Z even with f1 (no lam_r for doubling)
    's2_lam_no_lr': (2, lambda o, W, O, post: dict(lam=o.ite(W['ctrl'], O[1], o.zero()))),
    # step 2: f1 is not cleared
    's2_f1_keep': (2, lambda o, W, O, post: dict(f1=W['f1'])),
    # step 3: x3 = x + 2a instead of 3a
    's3_2a': (3, lambda o, W, O, post: dict(x=o.ite(W['ctrl'], o.modadd(W['x'], o.moddbl(W['a'])), W['x']))),
    # step 3: y3 is not zeroed (the subtraction is skipped)
    's3_no_sub': (3, lambda o, W, O, post: dict(y=W['y'])),
    # step 4: x4 = x + Z3 instead of x - Z3
    's4_plus': (4, lambda o, W, O, post: dict(x=o.modadd(W['x'], O[0]))),
    # step 5: x5 = a - x without the condition ctrl
    's5_neg_always': (5, lambda o, W, O, post: dict(x=o.modadd(W['a'], o.modneg(W['x'])))),
    # step 5: y5 does not subtract b
    's5_no_b': (5, lambda o, W, O, post: dict(y=W['y'])),
    # step 6: the result without subtracting a for P = -Q (x2 = x1)
    's6_no_cmodsub': (6, lambda o, W, O, post: dict(x=W['x'])),
    # step 6: y is not changed at all
    's6_y_keep': (6, lambda o, W, O, post: dict(y=W['y'])),
}


def apply(name):
    k, fix = MUTANTS[name]
    saved = Z.ZSTEPS[k]
    Z.ZSTEPS[k] = _wrap(k, post_fix=fix)
    return lambda: Z.ZSTEPS.__setitem__(k, saved)
