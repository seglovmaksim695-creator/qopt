# SPDX-License-Identifier: Apache-2.0
"""compose.kal -- the contract chain of the Kaliski inversion circuit (qoptct.crypto.kaliski, legacy=False).

One class, Chain, walks the algorithm P_std (Kaliski's almost-inverse in the "frame of w", exact integer semantics)
segment by segment. Every segment has a CONTRACT written once over a boolean algebra (compose.bv):

    spec(al, val) -> (posts, asserts, prefacts, postfacts, proved)
        val(wire)   value of a wire at the entry of the segment
        posts       {wire: value at the exit}; every other wire keeps its value
        proved      [(name, bit)]: properties of the entry state that are PROVED under the precondition (they make
                    the literal contract an exact integer statement, e.g. "no borrow")
        asserts     [(name, bit)]: MATHEMATICAL ASSUMPTIONS on the entry state (properties of the algorithm that are
                    proved on paper, compose/PROOF.md: M1..M6), used as the precondition
        prefacts    facts (see fact_bit) that must have been ESTABLISHED by earlier segments (machine-proved
                    postconditions carried over wires that were not touched in between), used as the precondition
        postfacts   facts proved at the exit of this segment

Two modes of the same walk:
  * prove (lifted circuit given): for every segment the slice of the REAL lifted circuit between two certificate
    marks is executed symbolically on ROBDDs over all entry values at full width, and it is proved that under the
    precondition  (i) every AND target is 0,  (ii) every erase condition holds,  (iii) every wire ends with the value
    of the contract (ancillas: 0),  (iv) the postfacts hold;  in clean mode also  (v) the mirror slice of the inverse
    half, started from the exit state, satisfies (i), (ii) and returns every wire to its entry value.
    The certificate (marks: positions and wire names, from KaliskiCircuit.hints()) is NOT trusted: a wrong mark
    makes a proof fail, it cannot make a wrong circuit pass.
  * run (inputs xs given): P_std itself is executed bit-sliced on the inputs; all asserts, facts and the final
    result are checked. This is the TEST of the mathematical assumptions (exhaustive for small p, sampled at n=256).
"""
from . import bv
from .walker import Walker, ProofFailure, register_fact, fact_bit, simulate_flat, inter as _inter

__all__ = ['Chain', 'ProofFailure', 'wspec', 'fact_bit', 'simulate_flat']


def wspec(al, a, b, wp):
    """w_k of Kaliski's step in the frame of wp = w_{k-1}: (u, v) = (b, a) if wp else (a, b);
    w = [the step operates on v] = [u odd] and ([v even] or [u <= v])."""
    u0 = al.ite(wp, b[0], a[0])
    v0 = al.ite(wp, a[0], b[0])
    lee = al.ite(wp, bv.le(al, b, a), bv.le(al, a, b))
    return al.and_(u0, al.or_(al.not_(v0), lee))


def _geodd(al, val, key):
    a = [val(w) for w in key[1]]
    b = [val(w) for w in key[2]]
    return al.or_(al.not_(a[0]), bv.le(al, b, a))


def _wdef(al, val, key):
    a = [val(w) for w in key[3]]
    b = [val(w) for w in key[4]]
    return al.not_(al.xor(val(key[1]), wspec(al, a, b, val(key[2]))))


# ('geodd', A, B): A_0 implies value(B) <= value(A);   ('wdef', w, wprev, A, B): w = wspec(A, B, wprev)
register_fact('geodd', lambda key: key[1] + key[2], _geodd)
register_fact('wdef', lambda key: (key[1],) + ((key[2],) if isinstance(key[2], int) else ()) + key[3] + key[4], _wdef)


class Chain(Walker):
    def __init__(self, p, reduction, hints, lifted=None, xs=None, clean=False, progress=None):
        self.p = p
        self.n = n = p.bit_length()
        self.red = reduction
        init = None
        if lifted is None:
            xs = list(xs)
            init = {}
            for i in range(n):
                m = 0
                for j, x in enumerate(xs):
                    if (x >> i) & 1:
                        m |= 1 << j
                init[i] = m
        Walker.__init__(self, hints, range(n), lifted=lifted, init=init, labels=xs, clean=clean,
                        reserved=range(n, 2 * n) if clean else (), progress=progress)

    # ------------------------------------------------------------------ the chain: P_std
    def run(self):
        n, p = self.n, self.p
        need = self.need
        X = list(range(n))
        pos, info = self.mark('init')
        U, R, S = self.fresh_list(info['U']), self.fresh_list(info['R']), self.fresh_list(info['S'])
        need(len(U) == n and len(R) == n and len(S) == n, 'certificate: bad register lengths')

        def spec(al, val):
            posts = {U[i]: (al.T if (p >> i) & 1 else al.F) for i in range(n)}
            posts[S[0]] = al.T
            return posts, [], [], [], []
        self.seg('init', None, pos, U + [S[0]], spec)                  # u = p, v = x, r = 0, s = 1
        A, B, wprev = U, X, None
        for k in range(2 * n):
            if self.progress and k % 16 == 0:
                self.progress('step %d/%d' % (k, 2 * n))
            wa = min(n, max(1, 2 * n - 1 - k))
            wb = min(wa, max(1, (2 * n - k) // 2))
            self.assume_zero(A[wa:] + B[wa:], 'M1-width')              # u, v < 2^wa
            if k == 2 * n - 1:
                # M5: the state is (u, v) = (1, 0) with w_{2n-2} = 1: A = v is even, w = 1, d = 0, e = 0: the step
                # is A /= 2 (relabelling) and s *= 2 (relabelling) only
                pos, info = self.mark('last')
                self.assume_zero([A[0]], 'M5-last-even')

                def spec(al, val, A=A, B=B, wprev=wprev):
                    wv = wspec(al, [val(x) for x in A], [val(x) for x in B], val(wprev))
                    return {}, [('M5-last-w', wv)], [], [], []
                self.seg('last', k, pos, [wprev] + A + B, spec)
                A = A[1:] + A[:1]
                S = self._double(S, info.get('snew'))
                continue
            if k == 2 * n - 2:
                w, wkey = None, None
            else:
                pos, info = self.mark('cmp')
                w = self.fresh(info.get('w'))
                wkey = ('wdef', w, wprev, tuple(A), tuple(B))

                def spec(al, val, A=A, B=B, wprev=wprev, w=w, wkey=wkey):
                    wv = wspec(al, [val(x) for x in A], [val(x) for x in B], val(wprev))
                    return {w: wv}, [], [], [wkey], []
                self.seg('cmp', k, pos, [wprev] + _inter(A[:wa], B[:wa]) + [w], spec)
            pos, info = self.mark('swap')

            def spec(al, val, A=A, B=B, R=R, S=S, wprev=wprev, w=w, wkey=wkey, k=k):
                a = [val(x) for x in A]
                b = [val(x) for x in B]
                r = [val(x) for x in R]
                s = [val(x) for x in S]
                wp = val(wprev)
                asserts, prefacts = [], []
                if w is None:
                    wv = al.T
                    asserts.append(('M5-w', wspec(al, a, b, wp)))      # w_{2n-2} = 1
                else:
                    wv = val(w)
                    prefacts.append(wkey)
                d = al.xor(wp, wv)
                nd = al.not_(d)
                if k >= 1:
                    asserts.append(('M2-P3', al.or_(nd, al.not_(al.xor(a[0], b[0])))))
                asserts.append(('M4-rs', al.or_(nd, bv.allzero(al, s[n:]))))
                posts = {}
                for i in range(n):
                    posts[A[i]] = al.ite(d, b[i], a[i])
                    posts[B[i]] = al.ite(d, a[i], b[i])
                    posts[R[i]] = al.ite(d, s[i], r[i])
                    posts[S[i]] = al.ite(d, r[i], s[i])
                return posts, asserts, prefacts, [('imp', (A[0], B[0])), ('geodd', tuple(A), tuple(B)),
                                                  ('nand', (R[0], S[0]))], []
            self.seg('swap', k, pos, [wprev, w] + _inter(A, B, R, S[:n]) + S[n:], spec)
            pos, info = self.mark('sub')
            e = self.fresh(info.get('e'))

            def spec(al, val, A=A, B=B, e=e, wa=wa):
                # A and B are < 2^wa here (M1): the subtraction is written on wa bits, and "no borrow" (proved
                # from the fact geodd: e -> B <= A) makes it the exact integer subtraction
                a = [val(x) for x in A[:wa]]
                b = [val(x) for x in B[:wa]]
                ev = a[0]
                diff, bor = bv.sub(al, a, [al.and_(ev, x) for x in b])
                posts = {A[i]: diff[i] for i in range(wa)}
                posts[e] = ev
                return posts, [], [('geodd', tuple(A), tuple(B))], [], [('no borrow in u - e*v', al.not_(bor))]
            # parametrised preconditions: fact imp (e -> B_0: P2) and M3 (e -> B < 2^wb)
            sb = [(B[0], lambda al, val, A=A, B=B: al.or_(val(B[0]), val(A[0])), None, ('imp', (A[0], B[0])),
                   lambda al, val, A=A, B=B: fact_bit(al, val, ('imp', (A[0], B[0]))))]
            for i in range(wb, wa):
                sb.append((B[i], lambda al, val, A=A, B=B, i=i: al.and_(val(B[i]), al.not_(val(A[0]))), 'M3-vshort',
                           None, lambda al, val, A=A, B=B, i=i: al.or_(al.not_(val(A[0])), al.not_(val(B[i])))))
            self.seg('sub', k, pos, _inter(A, B) + [e], spec, sb)      # e = u_0; u -= e*v (exact: u >= v if e)
            self.require_zero(A[0], 'u is even after the subtraction')
            A = A[1:] + A[:1]                                          # u /= 2
            pos, info = self.mark('add')

            def spec(al, val, R=R, S=S, e=e):
                ev = val(e)
                r = [val(x) for x in R]
                s = [val(x) for x in S]
                asserts = [('M4-sfit', al.or_(al.not_(ev), bv.allzero(al, s[n:])))]
                sm, c = bv.add(al, r, [al.and_(ev, x) for x in s[:n]])
                asserts.append(('M4-rfit', al.not_(c)))
                return {R[i]: sm[i] for i in range(n)}, asserts, [], [], []
            # parametrised precondition: fact nand (not (R_0 and S_0): P1)
            sb = [(S[0], lambda al, val, R=R, S=S: al.and_(val(S[0]), al.not_(val(R[0]))), None, ('nand', (R[0], S[0])),
                   lambda al, val, R=R, S=S: fact_bit(al, val, ('nand', (R[0], S[0]))))]
            self.seg('add', k, pos, [e] + _inter(R, S[:n]) + S[n:], spec, sb)   # r += e*s (exact, fits n bits)
            S = self._double(S, info.get('snew'))                      # s *= 2
            wprev = w if w is not None else 'ONE'
        need(len(S) == 2 * n, 'internal: s register length')
        # frame of w_{2n-1} = 1: u = B, v = A, r = S, s = R
        self.u, self.v, self.sreg = B, A, R
        if self.red == 'generic':
            out = self._reduce_generic(S)
        else:
            out = self._reduce_special(S)
        self.out = out
        return self.finish(out)

    def _double(self, S, snew):
        n = self.n
        S = [self.fresh(snew)] + S
        if len(S) > 2 * n:
            top = S.pop()
            self.assume_zero([top], 'M6-rtop')                         # r < 2^k: bit 2n of the doubled register is 0
            self.owned.discard(top)
        return S

    def _reduce_generic(self, S):
        """long division of r (2n bits) by p: n+1 steps  W = 2*rem + next bit < 2p,  rem' = W mod p."""
        n, p = self.n, self.p
        pos, info = self.mark('red_init')
        ext = self.fresh_list(info.get('ext'))
        R2 = S + ext
        self.need(len(R2) == 2 * n + 1, 'certificate: bad extension of r')
        self.need((1 << (n - 1)) < p < (1 << n), 'p must have exactly n bits')
        self.establish(('lt', tuple(R2[n + 1:2 * n + 1]), p))
        for j in range(n, -1, -1):
            win = R2[j:j + n + 1]
            pos, info = self.mark('div')
            f = self.fresh(info.get('f'))

            def spec(al, val, win=win, f=f):
                wv = [val(x) for x in win]
                d, bor = bv.sub(al, wv, bv.const(al, p, n + 1))
                # literally: f = bit n of (W - p) mod 2^(n+1); rem = (W - p) + f*p on n bits = (f ? W : W - p).
                # Proved under W < 2p: f = [W < p], i.e. rem = W mod p; and (postfact) rem < p
                fl = d[n]
                posts = {win[i]: al.ite(fl, wv[i], d[i]) for i in range(n)}
                posts[win[n]] = al.F
                posts[f] = fl
                return (posts, [], [('lt', tuple(win[1:]), p)], [('lt', tuple(win[:n]), p)],
                        [('f = [W < p]', al.not_(al.xor(fl, bor)))])
            self.seg('div', j, pos, win + [f], spec)
        return R2[:n]

    def _reduce_special(self, S):
        """folding rounds (literal modular additions; exactness by the fold lemma, PROOF.md, with its numeric side
        conditions checked here), then the binary search of the quotient."""
        n, p = self.n, self.p
        need = self.need
        m = n
        cur, Bd = S, 1 << (2 * n)
        self.establish(('lt', tuple(cur), Bd))
        while self.peek() == 'sp_const':
            pos, info = self.mark('sp_const')
            K = info.get('K')
            Aw = self.fresh_list(info.get('A'))
            L = len(Aw)
            need(isinstance(K, int) and 0 <= K < (1 << L) and L >= m, 'certificate: bad fold constant')
            need(('lt', tuple(cur), Bd) in self.facts, 'internal: bound fact lost')
            Hmax = (Bd - 1) >> m
            Lw, Hw = cur[:m], cur[m:][:Hmax.bit_length()]

            def spec(al, val, Aw=Aw, K=K):
                return {Aw[i]: (al.T if (K >> i) & 1 else al.F) for i in range(len(Aw))}, [], [], [], []
            self.seg('sp_const', None, pos, Aw, spec)
            pos, info = self.mark('sp_addL')

            def spec(al, val, Aw=Aw, Lw=Lw):
                a = [val(x) for x in Aw]
                l = [val(x) for x in Lw] + [al.F] * (len(Aw) - len(Lw))
                sm, _ = bv.add(al, a, l)
                return {Aw[i]: sm[i] for i in range(len(Aw))}, [], [], [], []
            self.seg('sp_addL', None, pos, _inter(Aw, Lw), spec)
            terms = []
            while self.peek() == 'sp_term':
                pos, info = self.mark('sp_term')
                s_, e_ = info.get('s'), info.get('e')
                need(s_ in (1, -1) and isinstance(e_, int) and 0 <= e_ < L, 'certificate: bad fold term')
                sub_, h = Aw[e_:], Hw[:L - e_]

                def spec(al, val, sub_=sub_, h=h, s_=s_):
                    a = [val(x) for x in sub_]
                    hv = [val(x) for x in h] + [al.F] * (len(sub_) - len(h))
                    res = bv.add(al, a, hv)[0] if s_ > 0 else bv.sub(al, a, hv)[0]
                    return {sub_[i]: res[i] for i in range(len(sub_))}, [], [], [], []
                self.seg('sp_term', e_, pos, _inter(sub_, h), spec)
                terms.append((s_, e_))
            neg = sum(Hmax << e_ for s_, e_ in terms if s_ < 0)
            pos_ = sum(Hmax << e_ for s_, e_ in terms if s_ > 0)
            Bn = K + (1 << m) + pos_
            need((sum(s_ << e_ for s_, e_ in terms) - (1 << m)) % p == 0, 'fold lemma: terms are not 2^m mod p')
            need(K % p == 0 and K >= neg and Bn <= (1 << L), 'fold lemma: side conditions fail')
            self.admit(('lt', tuple(Aw), Bn), 'fold lemma')
            cur, Bd = Aw, Bn
        while self.peek() == 'sp_q':
            pos, info = self.mark('sp_q')
            f = self.fresh(info.get('f'))
            Kq = info.get('K')
            L = len(cur)
            need(isinstance(Kq, int) and 0 < Kq < (1 << L) and Kq % p == 0, 'certificate: bad quotient constant')
            Bn = max(Kq, Bd - Kq) if Bd > Kq else Bd

            def spec(al, val, cur=cur, f=f, Kq=Kq, Bd=Bd, Bn=Bn):
                a = [val(x) for x in cur]
                ge = al.not_(bv.lt_const(al, a, Kq))
                d, _ = bv.sub(al, a, bv.const(al, Kq, len(a)))
                posts = {cur[i]: al.ite(ge, d[i], a[i]) for i in range(len(a))}
                posts[f] = ge
                return posts, [], [('lt', tuple(cur), Bd)], [('lt', tuple(cur), Bn)], []
            self.seg('sp_q', None, pos, cur + [f], spec)
            Bd = Bn
        need(Bd <= p, 'special reduction: the final bound %d bits is not <= p' % Bd.bit_length())
        return cur[:m]


