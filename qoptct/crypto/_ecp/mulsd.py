# SPDX-License-Identifier: Apache-2.0
"""k4/mulsd.py — умножение и квадрат на «условно-инвертируемых» сложениях (знаковые цифры s_i = 2y_i − 1 = ±1, знак — маска
CX, без частичных произведений И). Окно каждой строки фиксированной ширины, знак расширяется копией CX (растущий регистр).
  sd_mul(C, x, y, out):  2xy = Σ s_i 2^i x + (2^n − 1)x;          И = n(n+1) + n + 2n  ≈ n² + 4n
  sd_sq(C, x, out):      x² = Σ_{i≥1} s_i 2^i L_i − Σ x_i 4^i + 2^n x, L_i = x mod 2^i;  И ≈ n²/2 + 3.5n
(идея условно-инвертируемых сложений для квадрата — arXiv:2609.05625 §VII.A: ½n(n+3)−1; здесь — своя реализация + умножение)."""
from .engine import *
from .mulk import add_gen, sub_gen


def add_cin(C, Y, b, cin):
    """b += Y + cin (mod 2^len(b)); cin — провод."""
    m = len(b); Y = list(Y) + [ZERO] * (m - len(Y))
    X = [W(w) for w in b]
    c = [W(cin)]; R = [None]
    for i in range(m - 1):
        r = plan(X[i], Y[i], c[i]); c.append(compute(C, r, fold_ok=True)); R.append(r)
    for i in reversed(range(m)):
        if i + 1 < m and R[i + 1].kind != 'op': uncompute(C, R[i + 1])
        if Y[i] == ONE: C.x(b[i])
        elif Y[i][0] == 'w': C.cx(Y[i][1], b[i])
        folded = (i + 1 < m) and R[i + 1].kind == 'gid' and R[i + 1].fold
        if not folded:
            if c[i] == ONE: C.x(b[i])
            elif c[i][0] == 'w': C.cx(c[i][1], b[i])


def _signed_row(C, A, off, src, ctrl):
    """A[off : off+len(src)+2] += (ctrl ? +src : −src) (знаковое окно ширины len(src)+2). ctrl — провод y_i."""
    C.x(ctrl)                                   # m = ¬y_i: вычитать при y_i = 0
    for w in src: C.cx(ctrl, w)
    Y = [W(w) for w in src] + [W(ctrl), W(ctrl)]
    add_cin(C, Y, A[off:off + len(src) + 2], ctrl)
    for w in src: C.cx(ctrl, w)
    C.x(ctrl)


def sd_mul(C, x, y):
    """возвращает 2n проводов = x*y; |x| == |y| == n; x и y — разные провода."""
    n = len(x); assert len(y) == n and not (set(x) & set(y))
    A = [C.alloc() for _ in range(n + 1)]       # биты [0, n]
    for i in range(n):
        top = C.alloc(); C.cx(A[-1], top); A.append(top)   # расширение знака: теперь биты [0, n+i+1]
        _signed_row(C, A, i, x, y[i])
    # A: 2n+1 бит (знаковое Σ s_i 2^i x); + 2^n x − x = 2xy >= 0
    add_gen(C, [W(w) for w in x], A[n:2 * n + 1])
    sub_gen(C, [W(w) for w in x], A[:2 * n + 1])
    C.release(A[0])                              # 2xy чётно: младший бит 0
    return A[1:2 * n + 1]


def sd_sq(C, x):
    n = len(x)
    A = [C.alloc() for _ in range(2)]           # биты [0,1]
    for i in range(1, n):
        while len(A) < 2 * i + 2:
            t = C.alloc(); C.cx(A[-1], t); A.append(t)
        _signed_row(C, A, i, x[:i], x[i])
    while len(A) < 2 * n:
        t = C.alloc(); C.cx(A[-1], t); A.append(t)
    A = A[:2 * n] if len(A) == 2 * n else A
    add_gen(C, [W(w) for w in x], A[n:2 * n])
    D = []
    for i in range(n): D += [W(x[i]), ZERO]
    sub_gen(C, D[:len(A)], A[:2 * n])
    return A[:2 * n]


def _half_sum(C, x0, x1):
    L = max(len(x0), len(x1)) + 1
    s = [C.alloc() for _ in range(L)]
    lo, hi = (x0, x1) if len(x0) >= len(x1) else (x1, x0)
    for w, t in zip(lo, s): C.cx(w, t)
    add_gen(C, [W(w) for w in hi], s)
    return s


def kara(C, a, b, L, base=8, sq=False, gar=None):
    """Карацуба с мусором над базой sd_mul/sd_sq. Возвращает (out[2m], gar). |a|==|b|==m (при sq b игнорируется)."""
    if gar is None: gar = []
    m = len(a)
    if L == 0 or m <= base:
        return (sd_sq(C, a) if sq else sd_mul(C, a, b)), gar
    h = m // 2
    z0, _ = kara(C, a[:h], None if sq else b[:h], L - 1, base, sq, gar)
    z2, _ = kara(C, a[h:], None if sq else b[h:], L - 1, base, sq, gar)
    out = list(z0) + list(z2)
    sa = _half_sum(C, a[:h], a[h:])
    if sq:
        T, _ = kara(C, sa, None, L - 1, base, True, gar); gar += sa
    else:
        sb = _half_sum(C, b[:h], b[h:])
        T, _ = kara(C, sa, sb, L - 1, base, False, gar); gar += sa + sb
    sub_gen(C, [W(w) for w in z0], T)
    sub_gen(C, [W(w) for w in z2], T)
    add_gen(C, [W(w) for w in T[:m + 1]], out[h:2 * m])
    gar += list(T)
    return out, gar
