"""pq -- P_Q: the integer recursion given by the contracts of ``chain`` (one iteration = S1..S6).

It is the SAME text as the SMT text of ``lemmas_z3`` (the function ``step`` below was copied from there line by line), and it
is compared with the execution of the contracts (mode ``run`` of the client) after every iteration.  Plus the ghost
quantities of the invariant (ga, gb -- the exponents of the bounds u < 2^ga, v < 2^gb; cR, cS -- the comparison quotients mod p).
"""


def step(N, p, u, v, r, s, f):
    """one Qualtran iteration (contracts S1..S6), integers. -> (u, v, r, s, f, m, t, checks)
    checks: values of the assumptions and cross checks (all must be True)."""
    z1 = (v == 0)
    m1 = f and z1
    f1 = f and not z1                       # f xor m1
    t1 = m1
    a2 = f1 and u % 2 == 0
    m2 = m1 ^ (f1 and v % 2 == 0 and not a2)
    b2 = a2 ^ m2
    c3 = (u > v) and f1 and not b2
    a3 = a2 ^ c3
    m3 = m2 ^ c3
    u4, v4, r4, s4 = (v, u, s, r) if a3 else (u, v, r, s)
    c5 = f1 and not b2
    v5 = v4 - u4 if c5 else v4
    s5 = s4 + r4 if c5 else s4
    b6 = b2 ^ (m3 ^ a3)
    v6 = v5 // 2 if f1 else v5
    r6 = 2 * r4 if 2 * r4 < p else 2 * r4 - p
    u7, v7, r7, s7 = (v6, u4, s5, r6) if a3 else (u4, v6, r6, s5)
    a7 = a3 ^ (s7 % 2 == 0)
    checks = {'M-dbl': 0 <= r4 < p, 'M-a': (not a3) or 2 * r4 < p, 'M-s': a3 or s5 % 2 == 1,
              'M-carry': (not c5) or s4 + r4 < 2 ** N, 'no borrow': (not c5) or u4 <= v4,
              'v5 even': (not f1) or v5 % 2 == 0, 'b cleared': not b6, 'a cleared': not a7,
              'fits': all(0 <= w < 2 ** N for w in (u7, v7, r7, s7, v5, s5))}
    return u7, v7, r7, s7, f1, m3, t1, checks


def run(N, p, x):
    """-> (out, junk, states, bad): states[k] = (u, v, r, s, f) before iteration k (k = 0..2N); bad = failed checks."""
    u, v, r, s, f = p, x, 0, 1, True
    ms, ts, states, bad = [], [], [(u, v, r, s, f)], []
    for k in range(2 * N):
        u, v, r, s, f, m, t, ch = step(N, p, u, v, r, s, f)
        bad += [(k, name) for name, ok in ch.items() if not ok]
        ms.append(int(m))
        ts.append(int(t))
        states.append((u, v, r, s, f))
    t0 = ts[0]
    if not (v == 0 and not f and ((t0 and (u, s, r) == (p, 1, 0)) or (not t0 and (u, s) == (1, p)))):
        bad.append((2 * N, 'M-end'))
    out = (((p - r) % 2 ** N) ^ (p if t0 else 0))
    junk = 0
    for b in ms + ts:                        # MSB first: m_0 ... m_{2N-1}, t_0 ... t_{2N-1}
        junk = 2 * junk + b
    return out, junk, states, bad


def ghosts(N, p, x):
    """the states with the ghost values of the invariant of qkal_lemmas_z3: (k, f, u, v, r, s, ga, gb, cR, cS)."""
    nb = p.bit_length()
    out = []
    _, _, states, _ = run(N, p, x)
    ga = gb = nb
    for k, (u, v, r, s, f) in enumerate(states):
        if k > 0 and x and v >= 1 and f:
            pu, pv = states[k - 1][0], states[k - 1][1]
            if u != pu:
                ga -= 1
            elif v != pv:
                gb -= 1
        T = 2 ** k
        if x == 0:
            cR = cS = 0
        else:
            cR = (x * r + (u * T if v >= 1 and f else T)) // p
            cS = (x * s - v * T) // p
        out.append((k, f, u, v, r, s, ga, gb, cR, cS))
    return out
