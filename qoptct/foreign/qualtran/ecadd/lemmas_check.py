"""Check of the glue lemmas, separately from the schemes.
(1) sympy: the chord/tangent identity (the basis of EC_TAN) over the field of rational functions;
(2) full enumeration for small p: all curves (B!=0, non-singular), all pairs of points (with O), all lam_r for P=O or y=0:
    with the semantic MM/KINV (contracts of goals C/D) the instances of the lemmas MM0, K0, C1, C2, EC_DBL, EC_TAN,
    EC_RNOTO, EC_RNOTP are true; the chain of the six step contracts (zspecs, IntOps) gives P+Q and clean flags (fixed
    scheme); the original one outside G1-G3."""
import itertools, json, sys, time
from .zspecs import IntOps, ZSTEPS
from . import ec


def sympy_check():
    import sympy as sp
    xP, yP, xQ, A, X = sp.symbols('xP yP xQ A X')
    # B through P; yQ is free but lies on the curve: yQ^2 = xQ^3 + A xQ + B
    B = yP**2 - xP**3 - A * xP
    yQ = sp.symbols('yQ')
    lam = (yQ - yP) / (xQ - xP)
    xR = lam**2 - xP - xQ
    c = X**3 + A * X + B - (yP + lam * (X - xP))**2
    target = (X - xP) * (X - xQ) * (X - xR)
    diff = sp.together(sp.expand(c - target))
    num = sp.numer(diff)
    # reduction by the curve equation for Q
    num_red = sp.rem(sp.Poly(sp.expand(num), yQ), sp.Poly(yQ**2 - (xQ**3 + A * xQ + B), yQ))
    ok1 = sp.simplify(num_red.as_expr()) == 0
    # derivative at xP: c'(xP) = 3xP^2 + A - 2 lam yP ; and (from the factorisation) c'(xP) = (xP-xQ)(xP-xR)
    d1 = sp.diff(c, X).subs(X, xP)
    ok2 = sp.simplify(d1 - (3 * xP**2 + A - 2 * lam * yP)) == 0
    return dict(cubic_factorization=bool(ok1), derivative_formula=bool(ok2))


def semantic(n, p):
    R = pow(2, n, p)
    Ri = pow(R, p - 2, p)
    MM = lambda u, v: u * v * Ri % p
    KINV = lambda u: 0 if u == 0 else pow(u, p - 2, p) * R * R % p
    return MM, KINV


def chain(n, p, patched, a, b, x, y, lr, MM, KINV):
    """The six contracts (zspecs.IntOps); the boxes are the meaningful functions. -> (result, flag obligations, lemma instances)."""
    o = IntOps(n, p)
    v = {}
    pre_ok = True

    def step(k, W, O_fn):
        nonlocal pre_ok
        O = []
        for _ in range(3):
            pre, cl, post = ZSTEPS[k](o, W, O + [0] * 3, patched)
            if len(O) == len(cl):
                break
            kind, args = cl[len(O)]
            O.append(KINV(args['x']) if kind == 'KINV' else MM(args['x'], args['y']))
        pre, cl, post = ZSTEPS[k](o, W, O, patched)
        if not all(pre):
            pre_ok = False
        for kind, args in cl:
            if any(w >= p for w in args.values()):
                pre_ok = False
        return post, cl, O
    r1, _, _ = step(1, dict(a=a, b=b, x=x, y=y), None)
    r2, c2, O2 = step(2, dict(f1=r1['f1'], ctrl=r1['ctrl'], a=a, b=b, x=x, y=y, lam_r=lr), None)
    r3, c3, O3 = step(3, dict(ctrl=r2['ctrl'], a=a, b=b, x=r2['x'], y=r2['y'], lam=r2['lam']), None)
    r4, c4, O4 = step(4, dict(x=r3['x'], y=r3['y'], lam=r3['lam']), None)
    r5, c5, O5 = step(5, dict(ctrl=r3['ctrl'], a=a, b=b, x=r4['x'], y=r4['y'], lam_r=lr, lam=r4['lam']), None)
    r6, _, _ = step(6, dict(f1=r2['f1'], f2=r1['f2'], f3=r1['f3'], f4=r1['f4'], ctrl=r5['ctrl'], a=a, b=b,
                            x=r5['x'], y=r5['y']), None)
    return (r6['x'], r6['y']), pre_ok, dict(Z2=O2[1], Z3=O4[0], Z3p=O4[1], c4=c4)


def exhaustive(primes, max_seconds=600):
    t0 = time.time()
    out = dict(points=0, patched_bad=0, orig_bad_outside=0, orig_ok_inside=0, lemma_bad=0, groups=[0, 0, 0],
               primes=[])
    for p in primes:
        n = (2 * p - 1).bit_length()
        MM, KINV = semantic(n, p)
        R = pow(2, n, p)
        mont = lambda t: t * R % p
        # lemmas MM0, K0, C1, C2: exhaustive search over the arguments
        for u in range(p):
            if (KINV(u) == 0) != (u == 0):
                out['lemma_bad'] += 1
            for v_ in range(p):
                if (MM(u, v_) == 0) != (u == 0 or v_ == 0):
                    out['lemma_bad'] += 1
                if u and MM(u, MM(KINV(u), v_)) != v_:
                    out['lemma_bad'] += 1
                if u and MM(KINV(u), MM(u, v_)) != v_:
                    out['lemma_bad'] += 1
        for A in range(p):
            for B in range(1, p):
                if not ec.good_curve(A, B, p):
                    continue
                pts = [None] + [(X, Y) for X in range(p) for Y in range(p) if (Y * Y - X**3 - A * X - B) % p == 0]
                for P in pts:
                    lrs = ([mont((3 * P[0]**2 + A) * pow(2 * P[1], p - 2, p) % p)] if P is not None and P[1]
                           else sorted({0, 1, p - 1}))
                    for Q in pts:
                        for lr in lrs:
                            a, b = ec.enc(P, n, p); x, y = ec.enc(Q, n, p)
                            want = ec.enc(ec.add(P, Q, A, p), n, p)
                            out['points'] += 1
                            res, ok, info = chain(n, p, True, a, b, x, y, lr, MM, KINV)
                            if res != want or not ok:
                                out['patched_bad'] += 1
                            e1, e2 = a == x, b == (p - y) % p
                            e3, e4 = a == 0 and b == 0, x == 0 and y == 0
                            g1 = (not e3) and (not e4) and (not e1) and e2
                            g2 = (not e3) and e1 and b == y and b != 0 and lr == 0
                            g3 = e3 and x != 0 and y != 0 and lr == 0
                            res0, ok0, _ = chain(n, p, False, a, b, x, y, lr, MM, KINV)
                            good0 = res0 == want and ok0
                            if g1 or g2 or g3:
                                out['groups'][0 if g1 else 1 if g2 else 2] += 1
                                if ok0:          # on the group the release must fail
                                    out['orig_ok_inside'] += 1
                            elif not good0:
                                out['orig_bad_outside'] += 1
                            # instances of the EC lemmas (ordinary branch, binding done in the fixed chain)
                            gen = not e3 and not e4 and not (e1 and e2)
                            if gen:
                                Lam = lr if e1 else info['Z2']
                                XR = ((info['Z3'] - a) - x) % p
                                YR = (info['Z3p'] - b) % p
                                if e1 and not (y == b and b != 0):
                                    out['lemma_bad'] += 1           # EC_DBL
                                if XR == a and Lam != lr:
                                    out['lemma_bad'] += 1           # EC_TAN
                                if (XR, YR) == (0, 0) or (XR, YR) == (a, b):
                                    out['lemma_bad'] += 1           # EC_RNOTO / EC_RNOTP
                            if time.time() - t0 > max_seconds:
                                out['stopped_early'] = True
                                return out
        out['primes'].append(p)
    out['seconds'] = round(time.time() - t0, 1)
    return out


if __name__ == '__main__':
    res = dict(sympy=sympy_check())
    print(res, flush=True)
    res['exhaustive'] = exhaustive([5, 7, 11, 13, 17, 19, 23], max_seconds=float(sys.argv[2]) if len(sys.argv) > 2 else 900)
    print(json.dumps(res))
    if len(sys.argv) > 1:
        json.dump(res, open(sys.argv[1], 'w'), indent=1)
