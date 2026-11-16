# SPDX-License-Identifier: Apache-2.0
"""kprove — exact (proof-level) verification of large Clifford+T circuits with mid-circuit measurements and feed-forward,
for classical (reversible) references: arithmetic at 256-2048 bits, GF(2^m) multipliers, AES S-boxes/rounds.

    from qoptct.verify.kprove import kprove
    r = kprove("circ.qasm", "ref.qc", clean=True, budget_sec=600)      # -> dict(status, method, detail, ...)
    qoptct kprove circ.qasm ref.qc [--clean] [--budget 600] [--engine auto|bdd|anf|sample] [--json]

status: PROVED   - for EVERY basis input x and EVERY combination of measurement outcomes m:
                   branch K_m = c_m * (U_ref (x) |a(m)>)  (strict; clean: on inputs whose clean wires are 0)
        REFUTED  - a canonical violation (exact), with a witness (x, m) when one is extracted
        SAMPLED  - not proved; a sampled check passed (method/detail say which)
        UNKNOWN  - not proved and no sampled verdict
        ERROR    - parse / unsupported input

ENGINE "bdd" (main, new). Sum-over-paths (Amy, arXiv:1805.06908) with SYMBOLIC measurement outcomes, as in qoptct.verify,
but wire functions are ROBDDs (complement edges) instead of ANF, so ripple carries stay linear in size; the phase is a
sum of Z8-weighted Boolean atoms c*[f], f a BDD (exact: c*[a]*[b] = c*[a AND b] for Boolean a, b).
  h q       : new path variable y, phase += 4*[y AND w], w := y, amplitude sqrt2^-1
  measure   : new outcome variable m, new path variable z, phase += 4*[z AND (w XOR m)], w := m, sqrt2^-2
  if (cond) : gate controlled by the BDD of the condition over outcome variables
  reset     : symbolic measurement + X controlled by its outcome
Reduction (all exact identities), for a path variable y occurring in no wire and no guard:
  D(x) := F|y=1 - F|y=0 as a Z8 function (3 bit-BDDs, only the terms containing y)
  D == 4*[Q] (bits 0,1 vanish): sum_y = 2*delta(Q).  Q == 0 -> factor 2; Q == 1 -> zero;
        Q = z XOR R for a path variable z -> substitute z := R (sum over z of a delta);  Q free of path vars -> guard Q = 0
  D == 2 + 4*[Q] (bit0 = 0, bit1 = 1): sum_y = sqrt2 * w^(1 - 2[Q])
Normal form (no path variables left): data wires = reference outputs (canonical BDD equality), extra wires depend only on
outcomes, phase and guards depend only on outcomes  =>  PROOF for all x and all m. Any violation is canonical (REFUTED) and
a witness (x, m) is read off the BDD; it is confirmed by the independent sparse simulator of qoptct.verify when possible.
Variable order: inputs in order of first use, NEWEST VARIABLE ON TOP (ripple carries reuse the previous carry's BDD).

ENGINE "exact2": qoptct.verify.exact2 (circuits of cx / x / logical AND / AND erasure whose wires stay of GF(2) degree
<= 2, against a bilinear reference: all inputs of weight <= 2 decide everything; GF(2^163..571) multipliers in seconds).
Tried first by "auto" when the reference is bilinear.

MEMORY: kprove(max_mb=...) / --max-mb / $KPROVE_MAX_MB (default 900): budget of ONE call, counted as the growth of the
resident set over its value at the start of the call (MemGuard).  The budget applies to the BDD engine (and to the
counterexample search of this process): over budget the BDD engine stops with UNKNOWN, with the share of circuit ops
processed.  With engine="auto" the run then goes on: the "anf" engine runs in a separate worker process
(qoptct.verify.prove) WITHOUT this memory limit, and then the sampled check, so a small budget usually ends in SAMPLED
(exit code 2), not UNKNOWN.  To keep the whole call inside the budget use engine="bdd" (--engine bdd).

ENGINE "anf": qoptct.verify.prove (ANF path sum + trace miter; good for GF(2^m), degree-2 outputs).
ENGINE "sample": qoptct.verify.arith (a8 simulator, random inputs and branches; clean model compares all reference wires).
"""
import json
import os
import random
import re
import sys
import threading
import time

from .bdd import BDD, Z8, ONE, ZERO

from ._parse import parse_qasm_text, parse_qc, fold_h, fold_ops, cond_poly, bits


class Stuck(Exception):
    pass


class Budget(Exception):
    pass


# ============================================================================ memory guard
DEFAULT_MAX_MB = 900
NODE_BYTES = 200          # BDD storage per node (arrays + unique table; ~285 B/node incl. caches measured on GF(2^163))


def current_rss_mb():
    """Resident set size of this process NOW, in MB (10^6 bytes); None if the platform gives no way to read it.
    psutil if installed, else /proc/self/statm (Linux), else mach task_info (macOS, ctypes)."""
    try:
        import psutil
        return psutil.Process().memory_info().rss / 1e6
    except Exception:
        pass
    try:
        with open('/proc/self/statm') as fh:
            return int(fh.read().split()[1]) * os.sysconf('SC_PAGE_SIZE') / 1e6
    except Exception:
        pass
    if sys.platform == 'darwin':
        try:
            import ctypes
            import ctypes.util

            class _Info(ctypes.Structure):      # struct mach_task_basic_info
                _fields_ = [('virtual_size', ctypes.c_uint64), ('resident_size', ctypes.c_uint64),
                            ('resident_size_max', ctypes.c_uint64), ('user_time', ctypes.c_int32 * 2),
                            ('system_time', ctypes.c_int32 * 2), ('policy', ctypes.c_int32),
                            ('suspend_count', ctypes.c_int32)]
            lib = ctypes.CDLL(ctypes.util.find_library('c'))
            lib.mach_task_self.restype = ctypes.c_uint32
            info = _Info()
            cnt = ctypes.c_uint32(ctypes.sizeof(_Info) // 4)
            if lib.task_info(lib.mach_task_self(), 20, ctypes.byref(info), ctypes.byref(cnt)) == 0:   # 20 = MACH_TASK_BASIC_INFO
                return info.resident_size / 1e6
        except Exception:
            pass
    return None


def _peak_rss_mb():
    import resource
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return (peak if sys.platform == 'darwin' else peak * 1024) / 1e6


def resolve_max_mb(max_mb=None):
    """explicit max_mb > $KPROVE_MAX_MB > DEFAULT_MAX_MB (900)."""
    if max_mb is not None:
        return float(max_mb)
    return float(os.environ.get('KPROVE_MAX_MB', DEFAULT_MAX_MB))


class MemGuard:
    """Memory budget of ONE kprove call, counted as GROWTH of the resident set over its value when the call started
    (so a test process that already used 1 GB in earlier tests is not charged for it).  Without a current-RSS reading
    the growth of the peak RSS is used (never charges earlier peaks either).  Freed memory of earlier work can be
    reused without any RSS growth, so the growth is taken as max(RSS growth, NODE_BYTES * live BDD nodes) once the
    BDD engine has attached its node counter (nodes_fn)."""

    def __init__(self, max_mb=None, base_mb=None):
        self.max_mb = resolve_max_mb(max_mb)
        self.nodes_fn = None
        cur = current_rss_mb()
        self.use_peak = cur is None
        if base_mb is not None:
            self.base = float(base_mb)
        else:
            self.base = _peak_rss_mb() if self.use_peak else cur
        self.peak_growth = 0.0

    def growth(self):
        cur = _peak_rss_mb() if self.use_peak else current_rss_mb()
        if cur is None:
            cur = _peak_rss_mb()
        g = cur - self.base
        if self.nodes_fn is not None:
            g = max(g, self.nodes_fn() * NODE_BYTES / 1e6)
        self.peak_growth = max(self.peak_growth, g)
        return g

    def check(self):
        g = self.growth()
        if g > self.max_mb:
            raise Budget('memory limit: RSS grew by %.0f MB since the start of kprove (base %.0f MB) > %d MB '
                         '(max_mb / KPROVE_MAX_MB)' % (g, self.base, self.max_mb))
        return g


# ============================================================================ symbolic path sum on BDDs
class PS:
    def __init__(self, max_nodes=1_500_000, deadline=None, mem=None):
        self.B = BDD()
        self.kind = []           # per level: 'x' | 'y' | 'm'
        self.label = []
        self.wires = []
        self.wpv = []            # cached live-path-variable sets of the wires
        self.busy = {}           # path var -> number of wires/guards containing it
        self.terms = {}          # regular edge -> coefficient mod 8
        self.tpv = {}            # term edge -> path vars
        self.occ = {}            # path var -> set of term edges
        self.c0 = 0
        self.guards = []         # edges G with constraint G = 0
        self.e2 = 0
        self.zero = False
        self.live = set()
        self.thr = 0
        self.pvc = {}            # node -> frozenset of path vars (levels >= thr at computation time)
        self.stats = dict(h=0, meas=0, elim=0, hh=0, omega=0, guard=0, gc=0)
        self.max_nodes = max_nodes
        self.deadline = deadline
        self.dirty = set()       # wires changed since the last cut check
        self.cut_K = None        # cut wires whose BDD exceeds this many nodes (None = off)
        self.defs = {}           # function edge -> cut variable edge (for the reference run)
        self.cutvars = set()
        self.ncut = 0
        self.mset = set()
        self.max_cache = int(os.environ.get('KPROVE_MAX_CACHE', '1000000'))
        self.force_gc = int(os.environ.get('KPROVE_FORCE_GC', '0'))   # regression tests: collect every k checks
        self.gc_tick = 0
        self.mem = mem if mem is not None else MemGuard()
        self.max_mb = self.mem.max_mb
        self.mem.nodes_fn = self.B.nodes
        self.B.hook = self.mem_guard

    def mem_guard(self):
        self.mem.check()
        if self.deadline is not None and time.time() > self.deadline:
            raise Budget('time budget exceeded')

    # -------------------------------------------------------------- variables
    def newvar(self, kind, label=None):
        v, e = self.B.new_var()
        self.kind.append(kind)
        self.label.append(label)
        if kind == 'y':
            if not self.live:
                self.thr = v
            self.live.add(v)
        return v, e

    def kill(self, y):
        self.live.discard(y)
        if self.live:
            self.thr = min(self.live)
        else:
            self.thr = self.B.nvars

    # -------------------------------------------------------------- path-variable support
    def pv(self, e):
        n = e >> 1
        B = self.B
        if n == 0 or B.var[n] < self.thr:
            return frozenset()
        r = self.pvc.get(n)
        if r is not None:
            return r & self.live if r else r
        # iterative post-order
        thr = self.thr
        kind = self.kind
        stack = [n]
        while stack:
            m = stack[-1]
            if m in self.pvc:
                stack.pop()
                continue
            kids = []
            for ch in (B.hi[m] >> 1, B.lo[m] >> 1):
                if ch != 0 and B.var[ch] >= thr and ch not in self.pvc:
                    kids.append(ch)
            if kids:
                stack.extend(kids)
                continue
            stack.pop()
            s = set()
            for ch in (B.hi[m] >> 1, B.lo[m] >> 1):
                if ch != 0 and B.var[ch] >= thr:
                    s |= self.pvc[ch]
            if kind[B.var[m]] == 'y':
                s.add(B.var[m])
            self.pvc[m] = frozenset(s)
        r = self.pvc[n]
        return r & self.live if r else r

    # -------------------------------------------------------------- wires
    def set_wire(self, q, e):
        for y in self.wpv[q]:
            self.busy[y] -= 1
        s = self.pv(e)
        self.wires[q] = e
        self.wpv[q] = s
        self.dirty.add(q)
        for y in s:
            self.busy[y] = self.busy.get(y, 0) + 1

    def add_guard(self, g):
        if g == ZERO:
            return
        if g == ONE:
            self.zero = True
            return
        self.guards.append(g)

    # -------------------------------------------------------------- phase
    def add_phase(self, c, f):
        c &= 7
        if not c or f == ZERO:
            return
        if f == ONE:
            self.c0 = (self.c0 + c) & 7
            return
        if f & 1:
            self.c0 = (self.c0 + c) & 7
            f ^= 1
            c = (-c) & 7
        old = self.terms.get(f)
        if old is None:
            self.terms[f] = c
            s = self.pv(f)
            self.tpv[f] = s
            for y in s:
                self.occ.setdefault(y, set()).add(f)
        else:
            nv = (old + c) & 7
            if nv:
                self.terms[f] = nv
            else:
                self.del_term(f)

    def del_term(self, f):
        del self.terms[f]
        for y in self.tpv.pop(f, ()):
            s = self.occ.get(y)
            if s is not None:
                s.discard(f)

    # -------------------------------------------------------------- gates
    def gate(self, name, qs, cond=ONE):
        B = self.B
        w = self.wires
        if name == 'x':
            self.set_wire(qs[0], B.XOR(w[qs[0]], cond))
        elif name in ('cx', 'cnot', 'ccx', 'toffoli', 'mcx'):
            *cs, t = qs
            P = cond
            for c in cs:
                P = B.AND(P, w[c])
            self.set_wire(t, B.XOR(w[t], P))
        elif name in ('z', 's', 'sdg', 't', 'tdg'):
            k = {'z': 4, 's': 2, 'sdg': 6, 't': 1, 'tdg': 7}[name]
            self.add_phase(k, B.AND(cond, w[qs[0]]))
        elif name in ('cz', 'ccz', 'mcz'):
            P = cond
            for q in qs:
                P = B.AND(P, w[q])
            self.add_phase(4, P)
        elif name == 'swap':
            a, b = qs
            if cond == ONE:
                wa, wb = w[a], w[b]
                self.set_wire(a, wb)
                self.set_wire(b, wa)
            else:
                for (c, t) in ((a, b), (b, a), (a, b)):
                    self.gate('cx', (c, t), cond)
        elif name == 'h':
            if cond != ONE:
                raise ValueError('conditional H is not supported')
            self.hadamard(qs[0])
        elif name == 'y':
            self.add_phase(2, cond)
            self.gate('z', qs, cond)
            self.gate('x', qs, cond)
        elif name == 'id':
            pass
        else:
            raise ValueError('unsupported gate: ' + name)

    def hadamard(self, q):
        y, ye = self.newvar('y', 'h')
        self.add_phase(4, self.B.AND(ye, self.wires[q]))
        self.set_wire(q, ye)
        self.e2 -= 1
        self.stats['h'] += 1

    def measure(self, q, label=None):
        m, me = self.newvar('m', label)
        self.mset.add(m)
        z, ze = self.newvar('y', 'meas')
        self.add_phase(4, self.B.AND(ze, self.B.XOR(self.wires[q], me)))
        self.set_wire(q, me)
        self.e2 -= 2
        self.stats['meas'] += 1
        return m, me

    # -------------------------------------------------------------- substitution
    def subst(self, z, R):
        B = self.B
        if os.environ.get('KPROVE_DEBUG'):
            _t = time.time(); _n = B.nodes()
            ws = [q for q in range(len(self.wires)) if z in self.wpv[q]]
            _sup = B.support(R)
            print('supp x=%d m=%d y=%d' % tuple(sum(1 for v in _sup if self.kind[v] == k) for k in 'xmy'), end=' ')
            print('subst z=%d %s |R|=%d Rpv=%d wires=%s sizes=%s terms=%d' % (z, self.label[z], B.size(R), len(self.pv(R)),
                  ws[:6], [B.size(self.wires[q]) for q in ws[:6]], len(self.occ.get(z, ()))), flush=True)
        for q in range(len(self.wires)):
            if z in self.wpv[q]:
                self.set_wire(q, B.compose(self.wires[q], z, R))
        for f in list(self.occ.get(z, ())):
            if f not in self.terms:
                continue
            c = self.terms[f]
            self.del_term(f)
            self.add_phase(c, B.compose(f, z, R))
        ng = []
        for g in self.guards:
            if z in self.pv(g):
                for y in self.pv(g):
                    self.busy[y] -= 1
                g = B.compose(g, z, R)
                for y in self.pv(g):
                    self.busy[y] = self.busy.get(y, 0) + 1
            ng.append(g)
        self.guards = []
        for g in ng:
            if g == ONE:
                self.zero = True
            elif g != ZERO:
                self.guards.append(g)

    # -------------------------------------------------------------- reduction
    def check_budget(self):
        if self.deadline is not None and time.time() > self.deadline:
            raise Budget('time budget exceeded')
        B = self.B
        self.mem.check()
        if len(B.c_and) + len(B.c_xor) > self.max_cache:
            B.c_and = {}
            B.c_xor = {}
        if len(self.pvc) > self.max_cache:
            self.pvc = {}
        if self.force_gc:
            self.gc_tick += 1
            if self.gc_tick % self.force_gc == 0:
                self.collect()
        if self.B.nodes() > self.max_nodes:
            self.collect()
            if self.B.nodes() > self.max_nodes * 0.7:
                raise Budget('BDD too large (%d nodes)' % self.B.nodes())

    def try_eliminate(self, y, strict=False):
        B = self.B
        fs = [f for f in self.occ.get(y, ()) if f in self.terms]
        if not fs:
            self.kill(y)
            self.e2 += 2
            self.stats['elim'] += 1
            return True
        D = Z8.zero()
        F0 = Z8.zero()
        for f in fs:
            c = self.terms[f]
            f1 = B.restrict(f, y, 1)
            f0 = B.restrict(f, y, 0)
            D = Z8.add_term(B, D, c, f1)
            D = Z8.add_term(B, D, (-c) & 7, f0)
            F0 = Z8.add_term(B, F0, c, f0)
        # canonical remainder: value = b0 + 2*b1 + 4*b2 (mod 8), so local cancellations become atom merges
        rest = [(1, F0[0]), (2, F0[1]), (4, F0[2])]
        d0, d1, d2 = D
        if d0 != ZERO:
            return False
        if d1 == ZERO:
            Q = d2
            if Q not in (ZERO, ONE):
                qpv = self.pv(Q) - {y}
                z = R = None
                best = None
                for cand in sorted(qpv, reverse=True):
                    Q1 = B.restrict(Q, cand, 1)
                    Q0 = B.restrict(Q, cand, 0)
                    if Q1 == Q0 ^ 1:
                        key = (len(self.pv(Q0)), B.size(Q0))
                        if best is None or key < best[0]:
                            best = (key, cand, Q0)
                if best is not None:
                    if strict and best[0][0] > 0:
                        return False          # deferred: R still contains path variables
                    z, R = best[1], best[2]
                if z is None and qpv:
                    return False
            for f in fs:
                self.del_term(f)
            for c, f0 in rest:
                self.add_phase(c, f0)
            self.kill(y)
            self.e2 += 2
            if Q == ZERO:
                self.stats['elim'] += 1
            elif Q == ONE:
                self.zero = True
            elif z is not None:
                self.subst(z, R)
                self.kill(z)
                self.stats['hh'] += 1
            else:
                for v in self.pv(Q):
                    self.busy[v] = self.busy.get(v, 0) + 1
                self.add_guard(Q)
                self.stats['guard'] += 1
            return True
        if d1 == ONE:
            for f in fs:
                self.del_term(f)
            for c, f0 in rest:
                self.add_phase(c, f0)
            self.c0 = (self.c0 + 1) & 7
            self.add_phase(6, d2)
            self.kill(y)
            self.e2 += 1
            self.stats['omega'] += 1
            return True
        return False

    def reduce(self, final=False):
        """Phase 1 (strict): only eliminations whose substitution R is free of path variables. If nothing more is
        possible, ONE non-strict elimination (smallest R) is done and phase 1 resumes. Non-strict steps are taken only
        at the end of the circuit (final=True) or when many path variables pile up."""
        while not self.zero:
            changed = True
            while changed and not self.zero:
                changed = False
                for y in sorted(self.live, reverse=True):
                    if y not in self.live or self.busy.get(y, 0) > 0:
                        continue
                    if self.try_eliminate(y, strict=True):
                        changed = True
                        if self.zero:
                            return
                self.check_budget()
            if not (final or len(self.live) > 64):
                return
            done = False
            for y in sorted(self.live, reverse=True):
                if y in self.live and self.busy.get(y, 0) == 0 and self.try_eliminate(y, strict=False):
                    done = True
                    break
            if not done:
                return

    # -------------------------------------------------------------- cut points
    def cut_check(self, mlevels_fn):
        """Replace large classical wire functions F (no path variables, depending on inputs) by fresh variables v;
        defs[F] = v lets the reference run substitute the same v when it computes the same F (canonical BDDs)."""
        if not self.cut_K:
            self.dirty.clear()
            return
        B = self.B
        for q in sorted(self.dirty):
            e = self.wires[q]
            if self.wpv[q] or (e >> 1) == 0:
                continue
            if e in self.defs:
                self.set_wire(q, self.defs[e])
                continue
            if (e ^ 1) in self.defs:
                self.set_wire(q, self.defs[e ^ 1] ^ 1)
                continue
            if B.size(e) <= self.cut_K:
                continue
            if not (B.support(e) - mlevels_fn()):
                continue
            v, ve = self.newvar('x', 'cut')
            self.cutvars.add(v)
            self.defs[e] = ve
            self.ncut += 1
            self.set_wire(q, ve)
        self.dirty.clear()

    # -------------------------------------------------------------- GC
    def collect(self, extra=()):
        roots = list(self.wires) + list(self.terms) + list(self.guards) + list(extra) + list(self.defs) + \
            list(self.defs.values()) + list(getattr(self, 'extra_roots', []))
        mp = self.B.gc(roots)
        self.defs = {mp(k): mp(v) for k, v in self.defs.items()}
        if hasattr(self, 'extra_roots'):
            self.extra_roots = [mp(e) for e in self.extra_roots]
        self.wires = [mp(e) for e in self.wires]
        self.terms = {mp(f): c for f, c in self.terms.items()}
        self.guards = [mp(g) for g in self.guards]
        self.pvc = {}
        self.tpv = {f: self.pv(f) for f in self.terms}
        self.occ = {}
        for f, s in self.tpv.items():
            for y in s:
                self.occ.setdefault(y, set()).add(f)
        self.wpv = [self.pv(e) for e in self.wires]
        self.busy = {}
        for s in self.wpv + [self.pv(g) for g in self.guards]:
            for y in s:
                self.busy[y] = self.busy.get(y, 0) + 1
        self.stats['gc'] += 1
        return mp


# ============================================================================ circuit + reference
def count_gates(txt):
    t = len(re.findall(r'(?<![a-z])(t|tdg)\s+q\[', txt))
    cx = len(re.findall(r'(?<![a-z])cx\s+q\[', txt))
    ccx = len(re.findall(r'(?<![a-z])ccx\s+q\[', txt))
    meas = len(re.findall(r'measure', txt))
    N = int(re.search(r'qubit\[(\d+)\]', txt).group(1))
    return dict(T=t, CX=cx, CCX=ccx, measurements=meas, qubits=N)


def run_bdd(qasm_path, qc_path, clean=False, pre=None, budget_sec=600, max_nodes=1_500_000, reduce_every=64,
            order=None, verbose=False, spec=None, n=None, clean_wires=(), free_wires=(), cut_K=None, max_mb=None,
            mem=None):
    """BDD engine. pre: optional callable(ps, xedge) -> BDD edge of the input precondition (promise), where xedge(q)
    returns the BDD literal of input wire q. max_mb / mem: memory budget (growth of the RSS, see MemGuard).
    Returns a result dict."""
    t0 = time.time()
    if mem is None:
        mem = MemGuard(max_mb)
    txt = open(qasm_path).read()
    N, nb, ops = parse_qasm_text(txt)
    if spec is not None:
        n = n if n is not None else N
        names, ins, ref = ['w%d' % q for q in range(n)], None, []
    elif isinstance(qc_path, tuple):
        names, ins, ref = qc_path
        n = len(names)
    else:
        names, ins, ref = parse_qc(qc_path)
        n = len(names)
    if N < n:
        return dict(status='ERROR', detail='circuit has fewer qubits than the reference')
    ops = fold_ops(ops)
    if not isinstance(qc_path, tuple):
        ref = fold_h(ref)
    if any(g == 'h' for g, _ in ref):
        return dict(status='ERROR', detail='reference is not classical (H outside a foldable frame)')
    cleanw = {i for i, nm in enumerate(names) if nm not in ins} if (clean and ins is not None) else set()
    if spec is not None:
        cleanw = set(clean_wires)
    ps = PS(max_nodes=max_nodes, deadline=t0 + budget_sec if budget_sec else None, mem=mem)
    ps.cut_K = cut_K
    B = ps.B
    xv = {}

    def xedge(q):
        if q not in xv:
            v, e = ps.newvar('x', names[q])
            xv[q] = (v, None)
        return B.lit(xv[q][0])      # never cache edges: node numbers change at garbage collection
    if order:
        for q in order:
            if q < n and q not in cleanw:
                xedge(q)
    ps.wires = [ZERO] * N
    ps.wpv = [frozenset()] * N
    lazy = set(q for q in range(n) if q not in cleanw)
    cvar = {}

    def touch(qs):
        for q in qs:
            if q in lazy:
                lazy.discard(q)
                ps.wires[q] = xedge(q)

    def cond_bdd(s):
        def cv(k):
            if k not in cvar:
                raise ValueError('condition on unmeasured bit c[%d]' % k)
            return 1 << k
        P = cond_poly(s, cv)
        r = ZERO
        for mo in P:
            t = ONE
            for k in bits(mo):
                t = B.AND(t, B.lit(cvar[k][0]))
            r = B.XOR(r, t)
        return r
    since = 0
    try:
        for i, op in enumerate(ops):
            if op[0] == 'G':
                name, qs = op[1], op[2]
                touch(qs)
                if name == 'reset':
                    m, me = ps.measure(qs[0], 'reset')
                    ps.gate('x', qs, me)
                else:
                    ps.gate(name, qs)
            elif op[0] == 'M':
                _, q, k = op
                touch((q,))
                cvar[k] = ps.measure(q, 'c[%d]' % k)
                since = reduce_every
            else:
                c = cond_bdd(op[1])
                for name, qs in op[2]:
                    touch(qs)
                    ps.gate(name, qs, c)
            since += 1
            if (i & 15) == 0:
                ps.check_budget()
            if since >= reduce_every:
                ps.reduce()
                ps.cut_check(lambda: ps.mset)
                since = 0
                if ps.zero:
                    break
            if verbose and i % (verbose if verbose > 1 else 5000) == 0:
                print('  op %d/%d live=%d terms=%d nodes=%d t=%.1fs' % (i, len(ops), len(ps.live), len(ps.terms),
                                                                         B.nodes(), time.time() - t0), flush=True)
        touch(range(n))
        i = len(ops)
        ps.reduce(final=True)
        ps.mem_guard()
    except Budget as e:
        done = locals().get('i', 0)
        return dict(status='UNKNOWN', method='bdd', detail='%s [circuit ops processed: %d/%d, %d BDD nodes]'
                    % (e, done, len(ops), B.nodes()), time=round(time.time() - t0, 2), stats=ps.stats,
                    progress=(done, len(ops)))
    except (ValueError, KeyError, IndexError) as e:
        return dict(status='ERROR', method='bdd', detail=str(e), time=round(time.time() - t0, 2))
    res = dict(method='bdd', stats=ps.stats, nodes=B.nodes(), time=None)
    # ---------------- reference (classical netlist, same input variables) + its phase
    try:
        rw = [ZERO] * n
        for q in range(n):
            if q not in cleanw:
                rw[q] = xedge(q)
        rterms = {}
        rc0 = 0

        def rphase(c, f):
            nonlocal rc0
            if f == ZERO:
                return
            if f & 1:
                rc0 += c
                f ^= 1
                c = -c
            rterms[f] = (rterms.get(f, 0) + c) & 7
        defs = ps.defs

        def sub_def(q):
            e = rw[q]
            if e in defs:
                rw[q] = defs[e]
            elif (e ^ 1) in defs:
                rw[q] = defs[e ^ 1] ^ 1
        rcount = 0
        for g, qs in ref:
            rcount += 1
            if defs and (rcount & 255) == 0 and B.nodes() > max_nodes:
                raise Budget('reference netlist BDD too large (cut points not matched?)')
            if g == 'x':
                rw[qs[0]] ^= 1
            elif g in ('cx', 'ccx', 'mcx'):
                P = ONE
                for c in qs[:-1]:
                    P = B.AND(P, rw[c])
                rw[qs[-1]] = B.XOR(rw[qs[-1]], P)
                if defs:
                    sub_def(qs[-1])
            elif g in ('z', 'cz', 'ccz', 'mcz', 's', 'sdg', 't', 'tdg'):
                P = ONE
                for c in qs:
                    P = B.AND(P, rw[c])
                rphase({'s': 2, 'sdg': 6, 't': 1, 'tdg': 7}.get(g, 4), P)
            elif g == 'swap':
                rw[qs[0]], rw[qs[1]] = rw[qs[1]], rw[qs[0]]
            elif g == 'zero':           # reference program: wire is known to be erased here
                rw[qs[0]] = ZERO
            else:
                raise ValueError('reference gate not classical: ' + g)
        if spec is not None:
            from .ksym import SpecContext
            exp = spec(SpecContext(B, xedge, cleanw))
            for q, e in exp.items():
                if not 0 <= q < n:
                    raise ValueError('spec names wire %d outside 0..%d' % (q, n - 1))
                rw[q] = e
    except Budget as e:
        res.update(status='UNKNOWN', detail=str(e), time=round(time.time() - t0, 2))
        return res
    except ValueError as e:
        res.update(status='ERROR', detail=str(e))
        return res

    def expected_fn(x):
        xa = {xv[q][0]: x.get(q, 0) for q in xv}
        return sum(int(B.eval(rw[q], xa)) << q for q in range(n))
    res['expected_fn'] = expected_fn
    if ps.zero:
        res.update(status='REFUTED', detail='the path sum is identically zero (no branch is realised)')
        res['time'] = round(time.time() - t0, 2)
        return res
    if ps.live:
        res.update(status='UNKNOWN', detail='reduction stuck: %d path variables left (T-type phases not cancelled)'
                   % len(ps.live))
        res['time'] = round(time.time() - t0, 2)
        return res
    P = pre(ps, xedge) if pre is not None else ONE
    xlevels = {v for v, _ in xv.values()} | ps.cutvars
    mlevels = {v for v in range(B.nvars) if ps.kind[v] == 'm'}
    viol = []
    allowed = ONE
    for g in ps.guards:
        allowed = B.AND(allowed, g ^ 1)
    # data wires
    free = set(free_wires)
    for q in range(n):
        if q in free:
            continue
        d = B.AND(B.XOR(ps.wires[q], rw[q]), P)
        if d != ZERO:
            viol.append(('data', q, B.AND(d, allowed)))
    # extra wires: functions of the outcomes only
    for q in range(n, N):
        if q in free:
            continue
        sup = B.support(ps.wires[q])
        if sup & xlevels:
            viol.append(('ancilla', q, None))
    # guards: outcomes only
    for g in ps.guards:
        if B.support(g) & xlevels:
            viol.append(('guard', None, None))
    # phase: circuit - reference, must not depend on x
    terms = dict(ps.terms)
    for f, c in rterms.items():
        terms[f] = (terms.get(f, 0) - c) & 7
    terms = {f: c for f, c in terms.items() if c}
    # Phase S(x, m) = sum of c*[f] must not depend on x. Without forming the global sum (its BDD explodes when many
    # outcome variables are mixed): S(x, m) - S(x, 0) is a sum of one-variable differences along a path from 0 to m,
    # so it suffices that S(x, 0) and every d_k S = S|m_k=1 - S|m_k=0 are free of x. Each d_k S involves only the
    # terms containing m_k.
    xt = {}
    for f, c in terms.items():
        sup = B.support(f)
        if sup & xlevels:
            xt[f] = (c, sup & mlevels)
    phase_bad = False
    if xt:
        def xfree(Z):
            return not any(B.support(b) & xlevels for b in Z)
        S0 = Z8.zero()
        bym = {}
        for f, (c, ms) in xt.items():
            g = f
            for mv in ms:
                g = B.restrict(g, mv, 0)
            S0 = Z8.add_term(B, S0, c, g)
            for mv in ms:
                bym.setdefault(mv, []).append((c, f))
        if not xfree(S0):
            phase_bad = True
        for mv, items in bym.items():
            if phase_bad:
                break
            Dk = Z8.zero()
            for c, f in items:
                Dk = Z8.add_term(B, Dk, c, B.restrict(f, mv, 1))
                Dk = Z8.add_term(B, Dk, (-c) & 7, B.restrict(f, mv, 0))
            if not xfree(Dk):
                phase_bad = True
        if phase_bad and pre is not None:
            # check on the domain only: for every outcome assignment, each phase bit is constant over {x : P(x)}
            S = Z8.zero()
            for f, (c, ms) in xt.items():
                S = Z8.add_term(B, S, c, f)
            for b in S:
                E = B.exists(B.AND(b, P), xlevels)              # some x in the domain with bit 1
                A = B.exists(B.AND(b ^ 1, P), xlevels)          # some x in the domain with bit 0
                if B.AND(E, A) != ZERO:
                    break
            else:
                phase_bad = False
                res['note'] = 'phase constant on the domain (checked with the precondition)'
        if phase_bad:
            viol.append(('phase', None, None))
    res['time'] = round(time.time() - t0, 2)
    res['cuts'] = ps.ncut
    if viol and ps.ncut:
        res.update(status='UNKNOWN', detail='with %d cut points: %s (a mismatch may come from an unmatched cut; '
                   'not a refutation)' % (ps.ncut, [(v[0], v[1]) for v in viol[:6]]))
        return res
    if viol:
        res.update(status='REFUTED', detail='normal-form violations: %s' % [(v[0], v[1]) for v in viol[:6]])
        for v in viol:
            if v[2] not in (None, ZERO):
                a = B.pick_sat(v[2])
                x = [int(a.get(xv[q][0], 0)) if q in xv else 0 for q in range(n)]
                m = [0] * max(nb, 1)
                for k, (lv, _) in cvar.items():
                    m[k] = int(a.get(lv, 0))
                xa = {xv[q][0]: x[q] for q in xv}
                expd = sum(int(B.eval(rw[q], xa)) << q for q in range(n))
                res['witness'] = dict(kind=v[0], wire=v[1], x=x, m=m[:nb], expected=expd)
                break
        return res
    # normalisation self-check: (#allowed outcome assignments) * 2^e2 == 1
    try:
        cnt = B.satcount(allowed, mlevels) if mlevels else (1 if allowed == ONE else 0)
        ok = cnt > 0 and (cnt * 2 ** ps.e2 == 1 if ps.e2 <= 0 else False)
    except Exception:
        ok, cnt = None, None
    if ok is False:
        res.update(status='UNKNOWN', detail='INTERNAL: normalisation check failed (count %s, e2 %d)' % (cnt, ps.e2))
        return res
    res.update(status='PROVED', detail=('[%d cut points] ' % ps.ncut if ps.ncut else '') +
               'BDD path sum reached the normal form for all inputs and outcomes '
               '(%d measurements, %d H, %d HH, %d omega; %d BDD nodes)' % (ps.stats['meas'], ps.stats['h'],
                                                                            ps.stats['hh'], ps.stats['omega'], B.nodes()))
    return res


def _in_thread(fn, *a, **k):
    out = {}

    def run():
        try:
            out['r'] = fn(*a, **k)
        except Budget as e:
            out['r'] = dict(status='UNKNOWN', method='bdd', detail=str(e))
        except BaseException as e:  # noqa
            out['r'] = dict(status='ERROR', detail='%s: %s' % (type(e).__name__, e))
    threading.stack_size(512 * 1024 * 1024)
    th = threading.Thread(target=run)
    th.start()
    th.join()
    return out['r']


def confirm_spec_witness(qasm, w):
    """Spec mode: simulate branch m from input x independently; weight of outputs whose data bits differ from
    the specified value (1.0 also if the branch vanishes)."""
    N, nb, ops = parse_qasm_text(open(qasm).read())
    n = len(w['x'])
    k0 = sum(b << q for q, b in enumerate(w['x']))
    m = list(w['m']) + [0] * max(0, nb - len(w['m']))
    out, lg = branch_run(ops, k0, m)
    if out is None:
        return 1.0
    mask = (1 << n) - 1
    tot = sum(abs(a) ** 2 for a in out.values())
    good = sum(abs(a) ** 2 for k, a in out.items() if (k & mask) == w['expected'])
    return 1.0 - good / tot


def confirm_witness(qasm, qc, w):
    """Independent confirmation with the sparse simulator of qoptct.verify (classical-ish circuits scale to 1000s
    of qubits because the superposition stays small)."""
    n = len(w['x'])
    x = {q: b for q, b in enumerate(w['x'])}
    x2 = {q: 0 for q in range(n)}
    try:
        return confirm_witness2(qasm, qc, [x, x2], w['m'])
    except Exception as e:
        return 'error: %r' % e


def run_anf(qasm, qc, clean, budget_sec):
    from qoptct.verify import prove
    r = prove(qasm, qc, clean=clean, budget_sec=budget_sec)
    st = {'PROVED': 'PROVED', 'REFUTED': 'REFUTED'}.get(r['status'], 'UNKNOWN')
    return dict(status=st, method='anf:' + r.get('method', ''), detail=r.get('detail', ''), time=r.get('time'))


def run_sample(qasm, qc, clean, B=64, m_trials=6):
    from qoptct.verify import arith as VA
    t0 = time.time()
    try:
        ref = VA.ref_from_qc(qc, clean=clean)
        r = VA.check(qasm, ref, anc='zero_clean' if clean else 'zero', B=B, m_trials=m_trials)
    except Exception as e:
        return dict(status='UNKNOWN', method='sample', detail='sampler error: %r' % e)
    r.pop('detail', None)
    st = 'SAMPLED' if r.get('ok') else ('REFUTED' if r.get('worst', 0) and r.get('worst', 0) > 1e-6 and
                                         not r.get('error') else 'UNKNOWN')
    return dict(status=st, method='sample:' + str(r.get('method')),
                detail='a8 simulator: ok=%s worst=%s all_x=%s branches=%s err=%s' % (
                    r.get('ok'), r.get('worst'), r.get('proved_all_x'), r.get('branches'), r.get('error')),
                time=round(time.time() - t0, 2))


def branch_run(ops, x_key, outcomes, cap=1 << 16):
    """Independent sparse simulation of ONE branch (fixed outcomes) from basis input x_key. Renormalises after every
    measurement and keeps the scale as log2, so thousands of measurements do not underflow.
    -> (state dict normalised to max |amp| = 1, log2 of the true scale) or (None, None) if the branch has amplitude 0."""
    from .branchsim import run_branch_scaled   # renormalising branch simulator of the package (0.5.0)
    st, lg = run_branch_scaled(ops, {x_key: 1.0 + 0j}, outcomes, cap)
    if lg is None or not st:
        return None, None
    return st, lg


def confirm_witness2(qasm, qc, xs, m, expected_fn=None, n=None):
    """Branch m must equal t * (U_ref|x> (x) |a>) with the same t (modulus and phase) for all x in xs.
    Returns the largest deviation (> 1e-6: counterexample confirmed). Independent of the BDD engine."""
    import math
    from qoptct.verify.branchsim import apply_ref
    N, nb, ops = parse_qasm_text(open(qasm).read())
    if expected_fn is None:
        names, ins, ref = parse_qc(qc)
        n = len(names)
    mask = (1 << n) - 1
    m = list(m) + [0] * max(0, nb - len(m))
    ts = []
    worst = 0.0
    for x in xs:
        k0 = sum((x.get(q, 0) & 1) << q for q in range(n))
        out, lg = branch_run(ops, k0, m)
        want = {expected_fn(x): 1.0 + 0j} if expected_fn else apply_ref({k0: 1.0 + 0j}, ref, 1 << 16)
        if out is None:
            ts.append(None)
            continue
        wa = {}
        for k, a in out.items():
            wa[k >> n] = wa.get(k >> n, 0) + abs(a) ** 2
        anc = max(wa, key=wa.get)
        tot = sum(wa.values())
        worst = max(worst, (tot - wa[anc]) / tot)
        col = {k & mask: a for k, a in out.items() if k >> n == anc}
        i = max(want, key=lambda k: abs(want[k]))
        t = col.get(i, 0j) / want[i]
        if abs(t) < 1e-9:
            worst = max(worst, 1.0)
            ts.append(None)
            continue
        keys = set(col) | set(want)
        worst = max(worst, max(abs(col.get(k, 0j) - t * want.get(k, 0j)) for k in keys) / abs(t))
        ts.append((math.log2(abs(t)) + lg, t / abs(t)))
    nz = [t for t in ts if t is not None]
    if nz and len(nz) < len(ts):
        worst = max(worst, 1.0)
    for lgt, ph in nz:
        worst = max(worst, abs(2 ** (lgt - nz[0][0]) - 1), abs(ph - nz[0][1]))
    return float(min(worst, 1e6))


def cex_search(qasm, qc, clean=False, tries=40, seconds=60, seed=1, pre_fn=None, expected_fn=None, n=None,
               clean_wires=(), mem=None):
    """Random (x, m) pairs checked by the independent sparse simulator; dev > 1e-6 = confirmed counterexample.
    mem: MemGuard; the search stops (no counterexample) once the memory budget is exceeded."""
    from qoptct.verify._parse import nbits_of
    if expected_fn is None:
        names, ins, _ = parse_qc(qc)
        n = len(names)
        cleanw = {i for i, nm in enumerate(names) if nm not in ins} if (clean and ins is not None) else set()
    else:
        cleanw = set(clean_wires)
    nb = nbits_of(qasm)
    rng = random.Random(seed)
    t0 = time.time()
    for t in range(tries):
        if time.time() - t0 > seconds:
            break
        if mem is not None and mem.growth() > mem.max_mb:
            return None, 'memory budget exceeded'
        xs = []
        for _ in range(2):
            x = {q: (0 if q in cleanw else rng.randint(0, 1)) for q in range(n)}
            if pre_fn is not None:
                x = pre_fn(x, rng)
            xs.append(x)
        m = [rng.randint(0, 1) for _ in range(nb)] if t % 3 else [0] * nb
        try:
            dev = confirm_witness2(qasm, qc, xs, m, expected_fn=expected_fn, n=n)
        except Exception as e:
            return None, 'simulator error: %r' % e
        if dev > 1e-6:
            return dict(x=[xs[0][q] for q in range(n)], x2=[xs[1][q] for q in range(n)], m=m, dev=dev), None
    return None, None


def reg_value_lt(names, prefix, N):
    """pre builder: value of register `prefix`0.. (LSB = index 0) < N (BDD comparator)."""
    idx = {}
    for q, nm in enumerate(names):
        m = re.fullmatch(re.escape(prefix) + r'(\d+)', nm)
        if m:
            idx[int(m.group(1))] = q
    if not idx:
        raise ValueError('no wires named %s<i> in the reference' % prefix)
    width = max(idx) + 1

    def build(ps, xedge):
        B = ps.B
        lt = ZERO
        for i in range(width):
            xi = xedge(idx[i]) if i in idx else ZERO
            if (N >> i) & 1:
                lt = B.OR(xi ^ 1, B.AND(xi, lt))
            else:
                lt = B.AND(xi ^ 1, lt)
        if N >> width:
            lt = ONE
        return lt

    def sample(x, rng):
        v = rng.randrange(min(N, 1 << width))
        for i, q in idx.items():
            x[q] = (v >> i) & 1
        return x
    return build, sample


def parse_pre(spec, qc):
    """'lt:PREFIX:N[,lt:PREFIX2:N2]' -> (bdd builder, sampler) for kprove."""
    names, _, _ = parse_qc(qc)
    parts = []
    for item in spec.split(','):
        kind, prefix, N = item.split(':')
        if kind != 'lt':
            raise ValueError('unknown precondition ' + kind)
        parts.append(reg_value_lt(names, prefix, int(N, 0)))

    def build(ps, xedge):
        r = ONE
        for b, _ in parts:
            r = ps.B.AND(r, b(ps, xedge))
        return r

    def sample(x, rng):
        for _, s in parts:
            x = s(x, rng)
        return x
    return build, sample


def is_bilinear_ref(qc):
    """Reference of the form z ^= B(x, y) (+ linear): no Toffoli control is ever a target."""
    try:
        names, ins, ref = parse_qc(qc)
    except Exception:
        return False
    targets = set()
    ctrls = set()
    nt = 0
    for g, qs in ref:
        if g in ('mcx', 'ccx'):
            nt += 1
            ctrls |= set(qs[:-1])
            targets.add(qs[-1])
        elif g in ('cx', 'x'):
            targets.add(qs[-1])
        elif g == 'h':
            return False
    return nt > 0 and not (ctrls & targets)


def kprove(qasm, qc=None, clean=False, budget_sec=600, engine='auto', pre=None, confirm=True, sample=True, spec=None,
           max_mb=None, **kw):
    """Main entry. engine: auto | bdd | anf | sample | exact2.  auto: for a bilinear reference first exact2
    (qoptct.verify.exact2: exact for circuits of GF(2) degree <= 2, any size) and anf; then bdd; anf if still open
    and <= 4000 qubits; then sample.
    max_mb: memory budget of this call in MB, counted as the GROWTH of the process RSS over its value at the start of
    the call (default: $KPROVE_MAX_MB, else 900).  Over budget, the BDD engine stops with UNKNOWN (never a wrong
    verdict) and the counterexample search stops early.  The limit does not apply to the "anf" engine (a separate
    worker process): with engine='auto' a BDD UNKNOWN is followed by anf (unlimited memory) and by sampling, so the
    result is then usually SAMPLED; engine='bdd' keeps the call within max_mb.
    result['mem'] = dict(max_mb, base_mb, peak_growth_mb)."""
    mem = MemGuard(max_mb)
    kw['mem'] = mem
    qasm = os.path.abspath(qasm)
    qc = os.path.abspath(qc) if isinstance(qc, str) else qc
    if spec is not None:
        engine = 'bdd'
        kw['spec'] = spec
    pre_s = None
    if isinstance(pre, str):
        pre, pre_s = parse_pre(pre, qc)
    elif isinstance(pre, tuple):
        pre, pre_s = pre
    txt = open(qasm).read()
    cnt = count_gates(txt)
    t0 = time.time()
    tried = []
    res = None
    if engine == 'exact2' or (engine == 'auto' and isinstance(qc, str) and pre is None and is_bilinear_ref(qc)):
        from .exact2 import prove_degree2
        r0 = prove_degree2(qasm, qc, clean=clean)
        tried.append('exact2 (degree <= 2): %s %s' % (r0['status'], r0.get('detail', '')[:120]))
        if r0['status'] in ('PROVED', 'REFUTED') or engine == 'exact2':
            res = r0
            engine = 'done2'
    if engine == 'auto' and isinstance(qc, str) and pre is None and is_bilinear_ref(qc) and cnt['qubits'] <= 4000:
        r0 = run_anf(qasm, qc, clean, budget_sec / 2)
        tried.append('anf (bilinear reference): %s %s' % (r0['status'], r0.get('detail', '')[:120]))
        if r0['status'] in ('PROVED', 'REFUTED'):
            res = r0
            engine = 'done'
    if engine in ('auto', 'bdd'):
        res = _in_thread(run_bdd, qasm, qc, clean=clean, pre=pre, budget_sec=budget_sec, **kw)
        mem.growth()
        mem.nodes_fn = None          # release the BDD of this run (the counterexample search runs without it)
        tried.append('bdd: %s %s' % (res['status'], res.get('detail', '')[:120]))
        if res['status'] == 'UNKNOWN' and (qc or res.get('expected_fn')):
            cex, err = cex_search(qasm, qc, clean, seconds=min(120, budget_sec / 4), pre_fn=pre_s,
                                  expected_fn=None if qc else res.get('expected_fn'), n=kw.get('n'),
                                  clean_wires=kw.get('clean_wires', ()), mem=mem)
            tried.append('counterexample search: %s' % ('found' if cex else (err or 'none')))
            if cex:
                res = dict(status='REFUTED', method='counterexample',
                           detail='confirmed by the independent sparse simulator (dev %.3g)' % cex['dev'],
                           counterexample=cex)
    nonclassical = res is not None and res['status'] == 'ERROR' and 'not classical' in str(res.get('detail'))
    if engine == 'anf' or (engine == 'auto' and (res['status'] in ('UNKNOWN',) or nonclassical) and cnt['qubits'] <= 4000
                           and pre is None and not any(t.startswith('anf') for t in tried)):
        left = max(30.0, budget_sec - (time.time() - t0))
        r2 = run_anf(qasm, qc, clean, left)
        tried.append('anf: %s %s' % (r2['status'], r2.get('detail', '')[:120]))
        if res is None or r2['status'] in ('PROVED', 'REFUTED'):
            res = r2
    if qc and (res is None or (res['status'] in ('UNKNOWN',) and sample and engine in ('auto', 'sample', 'done')
                               and pre is None)):
        r3 = run_sample(qasm, qc, clean)
        tried.append('sample: %s %s' % (r3['status'], r3.get('detail', '')[:120]))
        if res is None or r3['status'] in ('SAMPLED', 'REFUTED'):
            res = dict(r3, status='SAMPLED' if r3['status'] == 'SAMPLED' else r3['status'])
    if (res['status'] == 'REFUTED' and confirm and not res.get('witness') and not res.get('counterexample')
            and res.get('method') == 'bdd' and (qc or res.get('expected_fn'))):
        cex, err = cex_search(qasm, qc, clean, seconds=min(120, budget_sec / 4), pre_fn=pre_s,
                              expected_fn=None if qc else res.get('expected_fn'), n=kw.get('n'), clean_wires=kw.get('clean_wires', ()),
                              mem=mem)
        res['counterexample'] = cex or ('not found (%s)' % (err or 'random search'))
    if res['status'] == 'REFUTED' and res.get('witness') and confirm:
        res['witness_confirm_dev'] = (confirm_witness(qasm, qc, res['witness']) if qc
                                      else confirm_spec_witness(qasm, res['witness']))
    res.pop('expected_fn', None)
    res['counts'] = cnt
    res['tried'] = tried
    res['total_time'] = round(time.time() - t0, 2)
    mem.growth()
    res['mem'] = dict(max_mb=mem.max_mb, base_mb=round(mem.base, 1), peak_growth_mb=round(mem.peak_growth, 1),
                      measure='peak RSS' if mem.use_peak else 'current RSS')
    return res


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(prog='qoptct kprove', description=__doc__.split('\n')[0],
                                 epilog='Exit code: 0 PROVED, 1 REFUTED, 2 SAMPLED/UNKNOWN (not proved), '
                                        '3 usage or input error. Memory guard: --max-mb, else $KPROVE_MAX_MB, else '
                                        '900 (MB of RSS growth during the proof). The guard limits the BDD engine; '
                                        'with --engine auto the anf engine (separate process, not limited) and '
                                        'sampling follow, so the result may be SAMPLED. Use --engine bdd to stay '
                                        'within the budget.')
    ap.add_argument('qasm')
    ap.add_argument('qc')
    ap.add_argument('--clean', action='store_true')
    ap.add_argument('--budget', type=float, default=600)
    ap.add_argument('--engine', default='auto', choices=['auto', 'bdd', 'anf', 'sample', 'exact2'])
    ap.add_argument('--pre', default=None, help="input promise, e.g. 'lt:a:0xFFFF,lt:b:0xFFFF' (register a<N, b<N)")
    ap.add_argument('--max-mb', type=float, default=None,
                    help='memory budget in MB of the BDD engine (growth of the RSS during the proof; default '
                         '$KPROVE_MAX_MB or 900); not applied to the anf engine')
    ap.add_argument('--json', action='store_true')
    ap.add_argument('-v', action='store_true')
    try:
        a = ap.parse_args(argv)
    except SystemExit as e:
        return 3 if e.code else 0
    for p in (a.qasm, a.qc):
        if not os.path.isfile(p):
            print('error: no such file: %s' % p, file=sys.stderr)
            return 3
    try:
        if a.max_mb is not None and a.max_mb <= 0:
            raise ValueError('--max-mb must be positive')
        r = kprove(a.qasm, a.qc, clean=a.clean, budget_sec=a.budget, engine=a.engine, pre=a.pre, verbose=a.v,
                   max_mb=a.max_mb)
    except (ValueError, KeyError, OSError) as e:
        print('error: %s' % e, file=sys.stderr)
        return 3
    if a.json:
        print(json.dumps(r, default=str))
    else:
        print('%s  (%s, %.1f s)  %s' % (r['status'], r.get('method'), r['total_time'], r.get('counts')))
        print('  ' + str(r.get('detail')))
        for t in r['tried']:
            print('  tried ' + t)
        if r.get('witness'):
            print('  witness', r['witness'], 'confirm dev', r.get('witness_confirm_dev'))
    return {'PROVED': 0, 'REFUTED': 1, 'ERROR': 3}.get(r['status'], 2)


if __name__ == '__main__':
    sys.exit(main())
