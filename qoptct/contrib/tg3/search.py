"""tg3/search.py: hill-climb / annealing over gauge freedoms (F1 initial ancilla wires, F2 garbage measured rows)."""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, time, random, math
sys.path.insert(0, (_QPKG + "/contrib/tg3"))
from core import *
from core import _canon_to_masks


class Search:
    def __init__(self, ctx, seed=0, body_seeds=5, move_garbage=True, move_init=True):
        self.ctx = ctx; self.rng = random.Random(seed); self.bs = body_seeds
        self.corr_cache = {}; self.body_cache = {}; self.net_cache = {}
        self.mg, self.mi = move_garbage, move_init
        self.n, self.N, self.h = ctx["n"], ctx["N"], ctx["h"]

    def corr_cost(self, Wend, body_for_prog=None):
        key = tuple(Wend)
        if key in self.corr_cache: return self.corr_cache[key]
        ctx = self.ctx
        Cs, order = corr_info(ctx, Wend)
        if order is None: r = (None, None, None)
        else:
            # cost of the conditional part only: use a dummy empty body
            try:
                if self.N > 20:   # fast proxy: number of CZ pairs of the fragments (frag_net alternatives only help at the end)
                    if any(v % 4 for C in Cs for k, v in C.items() if len(k) == 2): raise AssertionError
                    cc = sum(1 for C in Cs for k, v in C.items() if len(k) == 2 and v == 4)
                else:
                    prog = post.free_cz_pass(make_prog(ctx, [], Wend, Cs, order), ctx["n"])
                    cc = sum(1 for o in prog if o[0] == "g" and o[1] in ("cx", "cz"))
                r = (cc, Cs, order)
            except AssertionError:
                r = (None, None, None)
        self.corr_cache[key] = r
        return r

    def _prep_split(self):
        if hasattr(self, "PT"): return
        ctx = self.ctx; P2 = ctx["P2"]
        new = [m for m, w in P2.items() if w % 2]; self.PT = {m: 1 for m in new}
        Cq = dict(canon(P2))
        for k, v in canon(self.PT).items(): Cq[k] = (Cq.get(k, 0) - v) % 8
        self.Cq_masks = _canon_to_masks({k: v for k, v in Cq.items() if v})

    def _pairs(self, basis):
        N = self.N; C2 = {}
        for m, w in self.Cq_masks.items():
            for k2, v2 in coeffs([(pp.solve(basis, N, m), w)]).items(): C2[k2] = (C2.get(k2, 0) + v2) % 8
        cz = 0
        for k, v in C2.items():
            v %= 8
            if not v: continue
            if len(k) == 1 and v % 2: return None
            if len(k) == 2:
                if v % 4: return None
                if v == 4: cz += 1
            if len(k) > 2 and v: return None
        return cz

    def body_cost(self, Winit, Wend):
        key = (tuple(Winit), tuple(Wend))
        if key in self.body_cache: return self.body_cache[key]
        if self.N > 20: return self._fast_body(Winit, Wend, key)
        r = synth3(self.ctx, Winit, Wend, seeds=self.bs, base_seed=7, variants=("start", "end"))
        r = (cost_body(r[1]), r[1]) if r is not None else (None, None)
        self.body_cache[key] = r
        return r

    def _fast_body(self, Winit, Wend, key):
        ctx = self.ctx; N = self.N; best = None; self._prep_split()
        nets = self.net_cache.get(tuple(Winit))
        if nets is None:
            nets = []
            for sd in range(self.bs):
                W = list(Winit); rng = random.Random(11 + sd)
                try: net = pp.graysynth(dict(self.PT), W, N, rng)
                except Exception: net = None
                if net is not None: nets.append((sum(g == "cx" for g, _ in net), W))
            self.net_cache[tuple(Winit)] = nets
            if len(self.net_cache) > 300: self.net_cache.pop(next(iter(self.net_cache)))
            self.pairs_start = getattr(self, "pairs_start", {}); self.pairs_start[tuple(Winit)] = self._pairs(Winit)
        ps = self.pairs_start.get(tuple(Winit)); pe = self._pairs(Wend)
        rng = random.Random(3)
        for cxn, W in nets:
            M = [pp.solve(Wend, N, W[r]) for r in range(N)]
            ops = pp.lin_best(M, N, range(N), rng, 3)
            for extra in (ps, pe):
                if extra is None: continue
                c = cxn + len(ops) + extra
                if best is None or c < best: best = c
        r = (best, None) if best is not None else (None, None)
        self.body_cache[key] = r
        return r

    def cost(self, Winit, Wend):
        cc, Cs, order = self.corr_cost(Wend)
        if cc is None: return None
        bc, body = self.body_cost(Winit, Wend)
        if bc is None: return None
        return bc + cc

    def neighbors(self, Winit, Wend):
        n, N = self.n, self.N; rng = self.rng
        if self.mi and (not self.mg or rng.random() < 0.5):
            j = rng.randrange(n, N); r = rng.choice([k for k in range(N) if k != j and (k < n or True)])
            W2 = list(Winit)
            # keep a-part invertible: allowed unit ops (add any row, ancilla rows only target)
            W2[j] ^= W2[r]
            if not rank_ok(W2, N): return None
            return W2, Wend
        elif self.mg:
            i = rng.randrange(n, N); r = rng.choice([k for k in range(N) if k != i])
            E2 = list(Wend); E2[i] ^= E2[r]
            if not rank_ok(E2, N): return None
            return Winit, E2
        return None

    def run(self, seconds, T0=1.0, log=None):
        ctx = self.ctx; t_end = time.process_time() + seconds
        Winit = identity_init(ctx); Wend = list(ctx["W"])
        cur = self.cost(Winit, Wend); best = (cur, Winit, Wend); it = 0
        T = T0
        while time.process_time() < t_end:
            it += 1
            nb = self.neighbors(Winit, Wend)
            if nb is None: continue
            c = self.cost(*nb)
            if c is None: continue
            frac = 1 - (t_end - time.process_time()) / seconds; T = T0 * (1 - frac) + 0.05
            if c <= cur or self.rng.random() < math.exp(-(c - cur) / T):
                Winit, Wend = nb; cur = c
                if c < best[0]:
                    best = (c, Winit, Wend)
                    if log: log(f"it {it} best {c} ({time.process_time()-(t_end-seconds):.0f}s)")
        return best


def finalize(ctx, Winit, Wend, seeds=40, base_seed=1000):
    Cs, order = corr_info(ctx, Wend)
    body = synth2(ctx, Winit, Wend, seeds=seeds, base_seed=base_seed)
    ok = check_body(ctx, body, Winit, Wend)
    prog = make_prog(ctx, body, Wend, Cs, order)
    return body, prog, ok
