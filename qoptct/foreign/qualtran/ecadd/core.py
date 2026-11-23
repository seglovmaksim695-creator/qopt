"""Symbolic execution of an ECAdd step on ROBDD (qoptct.verify.bdd, unchanged).

The step program = the real logical operations (Qualtran -> bridge lower_bloqgraph -> translate rules -> QASM ->
qoptct lift; the gadget lemmas for AND/ERASE are checked by qoptct.check_gadget_lemmas) + black-box markers.
For ALL values of the variables that satisfy the assumption A (step precondition + ranges of the box outputs) it checks:
the target of every and is 0; the target of every erase equals a*b (xor l) EXACTLY; the inputs of every forward box equal
the inputs stated by the specification and lie in its domain (<p); a backward box receives exactly the outputs of its
paired forward call; the outputs equal the specification; every other wire is 0.  A is checked to be non-empty.
"""
import time
from qoptct.verify.bdd import BDD, ONE, ZERO
from qoptct.verify.compose.lift import lift, check_gadget_lemmas
from qoptct.foreign.qualtran.bridge.emit import emit_qasm
from . import bv

KIND = {'x': 0, 'cx': 1, 'and': 2, 'erase': 3}
import os as _os
MEMLOG = bool(_os.environ.get('EW_MEMLOG'))
MSB_ROOT = _os.environ.get('EW_LSB_ROOT') is None      # порядок переменных (по умолчанию старшие биты у корня)


class Fail(Exception):
    def __init__(self, msg, **detail):
        super().__init__(msg)
        self.detail = detail


def lsb(ws):
    return list(reversed(ws))


def disjoint(B, f, g):
    """True <=> f AND g == ZERO; без создания узлов (обход пар с ранним выходом)."""
    if f == ZERO or g == ZERO:
        return True
    var, hi, lo = B.var, B.hi, B.lo
    seen = set()
    stack = [(f, g)]
    while stack:
        f, g = stack.pop()
        if f == ZERO or g == ZERO or f == (g ^ 1):
            continue
        if f == ONE or g == ONE or f == g:
            return False
        key = (f << 40) | g if f < g else (g << 40) | f
        if key in seen:
            continue
        seen.add(key)
        nf, ng = f >> 1, g >> 1
        vf, vg = var[nf], var[ng]
        if vf >= vg:
            cf = f & 1
            f1, f0 = hi[nf] ^ cf, lo[nf] ^ cf
        else:
            f1 = f0 = f
        if vg >= vf:
            cg = g & 1
            g1, g0 = hi[ng] ^ cg, lo[ng] ^ cg
        else:
            g1 = g0 = g
        stack.append((f0, g0))
        stack.append((f1, g1))
    return True


def lower_step(step_bloq):
    from qoptct.foreign.qualtran.bridge.qualtran_layout import lower_bloqgraph
    from .opaque import opaque_atom_rule
    from .translate_ew import translate_ew
    low = lower_bloqgraph(step_bloq, extra_atom_rule=opaque_atom_rule, max_leaves=5_000_000, max_nodes=20_000_000,
                          record_events='top')
    tr = translate_ew(low)
    low['top_events'] = [dict(key=ev['path'].split('/', 1)[1], inputs={k: list(reversed(v)) for k, v in ev['inputs_msb'].items()},
                              outputs={k: list(reversed(v)) for k, v in ev['outputs_msb'].items()})
                         for ev in low.get('graph_events', []) if ev['path'].count('/') == 1]
    for key in ('wire_names', 'graph_events', 'source_wire_names'):
        low.pop(key, None)
        tr['layout'].pop(key, None)
    return low, tr


def build_program(step_bloq, use_lift=True, mutate=None, lowered=None, keep_ops=True):
    """mutate(tr, low) -> (ops, leaves): только для канареек (порча настоящей программы)."""
    low, tr = lowered if lowered is not None else lower_step(step_bloq)
    N = tr['nq']
    ops_all, leaves = tr['ops'], low['leaves']
    if mutate is not None:
        ops_all, leaves = mutate(tr, low)
    items, run = [], []

    def flush():
        if not run:
            return
        if use_lift:
            L = lift(emit_qasm(N, run))
            ops = L.ops
            if L.n != len(run):
                raise Fail('lift: число операций не совпало')
            for j, op in enumerate(run):
                if op[0] == 'erase':
                    (cz, a, b), = op[2]
                    want = (3, a, b, op[1], -1)
                elif op[0] == 'and':
                    want = (2, op[1], op[2], op[3], -1)
                elif op[0] == 'cx':
                    want = (1, op[1], op[2], 0, -1)
                else:
                    want = (0, op[1], 0, 0, -1)
                if tuple(ops[5 * j:5 * j + 5]) != want:
                    raise Fail('lift: операция не совпала с переводом', j=j, op=op)
            items.append(('lifted', ops, gstart[0]))
        else:
            items.append(('ops', list(run), gstart[0]))
        run.clear()

    gstart = [0]
    for g, op in enumerate(ops_all):
        if op[0] == 'opaque':
            flush()
            items.append(('opaque', leaves[op[1]], g))
            gstart[0] = g + 1
        else:
            if not run:
                gstart[0] = g
            run.append(op)
    flush()
    from .regions import regions
    reg = [(key, lo, hi) for key, (lo, hi, cnt) in regions(tr).items()]
    top_cx = [(key, tr['ops'][lo][1], tr['ops'][lo][2]) for key, lo, hi in reg
              if key.endswith(':CNOT') and hi == lo + 1 and tr['ops'][lo][0] == 'cx']
    if not keep_ops:
        ops_all_len = len(ops_all)
        return dict(N=N, items=items, nops=ops_all_len, inputs=low['inputs'], outputs=low['outputs'],
                    regions=reg, top_cx=top_cx, top_events=low.get('top_events', []))
    return dict(N=N, items=items, nops=len(ops_all), ops=ops_all, leaves=leaves, regions=reg, top_cx=top_cx,
                top_events=low.get('top_events', []),
                inputs=low['inputs'], outputs=low['outputs'])


def region_end(prog, type_name, occurrence=0):
    """Индекс (в переведённых операциях) сразу после occurrence-го верхнеуровневого блока данного типа."""
    hits = [hi for key, lo, hi in prog['regions'] if key.split(':', 1)[1] == type_name]
    if len(hits) <= occurrence:
        raise Fail(f'нет блока {type_name}#{occurrence}')
    return hits[occurrence] - 1


def region_outputs(prog, type_name, occurrence, reg):
    """Провода (LSB первым) выходного регистра reg у occurrence-го верхнеуровневого блока типа type_name
    (по настоящему графу соединений Qualtran: учитывает перестановки проводов внутри блока, напр. ModDbl)."""
    hits = [ev for ev in prog['top_events'] if ev['key'].split(':', 1)[1] == type_name]
    if len(hits) <= occurrence or reg not in hits[occurrence]['outputs']:
        raise Fail(f'нет выхода {reg} у блока {type_name}#{occurrence}')
    return hits[occurrence]['outputs'][reg]


def copy_targets(prog, src_wires):
    """Цели первой верхнеуровневой копии CNOT(src_i -> t_i) (в порядке src, LSB первым)."""
    m = {}
    for key, c, t in prog['top_cx']:
        if c in src_wires and c not in m:
            m[c] = t
    if any(w not in m for w in src_wires):
        raise Fail('копия регистра не найдена')
    return [m[w] for w in src_wires]


class Run:
    """Одно доказательство шага. spec: объект с методами words(), pre(), calls(), post() (см. specs.py)."""

    def __init__(self, prog, spec, n, p, gc_nodes=900_000, log=print, max_nodes=60_000_000):
        self.prog, self.spec, self.n, self.p = prog, spec, n, p
        self.B = BDD()
        self.log, self.gc_nodes, self.max_nodes = log, gc_nodes, max_nodes
        self.stats = dict(ops=0, ands=0, erases=0, opaque=0, gc=0, checks_nontrivial=0)

    # ---------------------------------------------------------------- variables
    def _alloc(self):
        B, n = self.B, self.n
        ins = self.prog['inputs']
        fwd = [it[1] for it in self.prog['items'] if it[0] == 'opaque' and it[1]['kind'] in ('MM', 'KINV')]
        self.fwd_leaves = fwd
        arith, garbage, flags = [], [], []
        for name, ws in ins.items():
            (arith if len(ws) == n else flags if len(ws) == 1 else garbage).append(('in', name, len(ws)))
        for nm in getattr(self.spec, 'named', ()):
            arith.append(('nm', nm, n))
        for c, leaf in enumerate(fwd):
            if leaf['kind'] == 'MM':
                arith.append(('op', (c, 'target'), n))
                garbage.append(('op', (c, 'qrom_indices'), len(leaf['regs_out']['qrom_indices'])))
                garbage.append(('op', (c, 'reduced'), 1))
            else:
                arith.append(('op', (c, 'x'), n))
                garbage.append(('op', (c, 'junk'), len(leaf['regs_out']['junk'])))
        words = {}
        for kind, key, w in garbage:            # самый низ
            words[(kind, key)] = [B.new_var()[1] for _ in range(w)]
        cols = {(kind, key): [None] * n for kind, key, w in arith}
        order = range(n) if MSB_ROOT else reversed(range(n))
        for i in order:                          # MSB_ROOT: младшие биты внизу, старшие — у корня (общие цепочки переносов)
            for kind, key, w in arith:
                cols[(kind, key)][i] = B.new_var()[1]
        words.update(cols)
        for kind, key, w in flags:              # флаги у корня
            words[(kind, key)] = [B.new_var()[1]]
        self.W = {key: words[('in', key)] for kind, key, w in arith + garbage + flags if kind == 'in'}
        self.O = [dict() for _ in fwd]
        for kind, key, w in arith + garbage:
            if kind == 'op':
                self.O[key[0]][key[1]] = words[('op', key)]
        self.N = {key: words[('nm', key)] for kind, key, w in arith if kind == 'nm'}
        self.spec.N = self.N
        self.nvars = B.nvars

    # ---------------------------------------------------------------- helpers
    def holds(self, d):
        """d == ZERO при допущении A (конъюнкция частей). Части добавляются «слоями» по общим переменным с d:
        если d несовместно уже с частью условий — несовместно и со всем A (достаточное условие, корректно);
        если не хватило всей связной компоненты — это настоящее нарушение (остальные части от d не зависят)."""
        if d == ZERO:
            return True
        self.stats['checks_nontrivial' if 'checks_nontrivial' in self.stats else 'checks'] += 1
        B = self.B
        used, frontier = set(), B.support(d)
        while True:
            new = [i for i, (e, s_) in enumerate(self.conj) if i not in used and not s_.isdisjoint(frontier)]
            if not new:
                return False
            used.update(new)
            key = frozenset(used)
            Ak = self.conj_cache.get(key)
            if Ak is None:
                Ak = ONE
                for i in sorted(used):
                    Ak = B.AND(Ak, self.conj[i][0])
                self.conj_cache[key] = Ak
            if disjoint(B, Ak, d):
                return True
            frontier = set()
            for i in new:
                frontier |= self.conj[i][1]

    def word_eq_holds(self, X, Y):
        if len(X) != len(Y):
            return False
        return all(self.holds(self.B.XOR(x, y)) for x, y in zip(X, Y))

    def full_A(self):
        A = ONE
        for e, _ in self.conj:
            A = self.B.AND(A, e)
        return A

    def components(self, sup):
        """Разбить части A на связную компоненту носителя sup и остальные компоненты (по общим переменным)."""
        rest = list(range(len(self.conj)))
        comps = []
        seeds = [set(sup)]
        while rest:
            cur = seeds.pop() if seeds else set(self.conj[rest[0]][1])
            members, grow = [], True
            while grow:
                grow = False
                for i in list(rest):
                    if not self.conj[i][1].isdisjoint(cur):
                        members.append(i); rest.remove(i); cur |= self.conj[i][1]; grow = True
            comps.append(members)
        return comps

    def component_and(self, members):
        A = ONE
        for i in members:
            A = self.B.AND(A, self.conj[i][0])
        return A

    def counterexample(self, d):
        """Выполняющий набор для (A и d) без построения полного A: компонента d + по набору на каждую прочую."""
        B = self.B
        comps = self.components(B.support(d))
        Ad = self.component_and(comps[0]) if comps else ONE
        self.A_d = Ad
        asg = dict(B.pick_sat(B.AND(Ad, d)) or {})
        for members in comps[1:]:
            asg.update(B.pick_sat(self.component_and(members)) or {})
        level = {}
        for name, ws in self.W.items():
            level[name] = ws
        out = {}
        for name, ws in list(self.W.items()) + [(f'o{c}.{r}', w) for c, o in enumerate(self.O) for r, w in o.items()]:
            v = 0
            for i, e in enumerate(ws):
                if e in (ONE, ZERO):
                    bit = int(e == ONE)
                else:
                    bit = asg.get(self.B.var[e >> 1], 0)
                if bit:
                    v |= 1 << i
            out[name] = v
        return out

    def fail(self, msg, d=None, **kw):
        self.viol = d
        cex = self.counterexample(d) if d is not None and d != ZERO else None
        raise Fail(msg, counterexample=cex, **kw)

    def maybe_gc(self, force=False):
        B = self.B
        if len(B.c_and) + len(B.c_xor) > 2_000_000:      # кэши операций растут быстрее узлов
            B.c_and, B.c_xor = {}, {}
        limit = max(self.gc_nodes, getattr(self, 'gc_limit', 0))
        if not force and B.nodes() < limit:
            return
        roots = list(self.val) + ([self.A] if self.A is not None else []) + [e for e, _ in self.conj]
        roots += [e for w in getattr(self, 'stack_saved', []) for e in w]
        roots += [e for _, e in getattr(self, 'preds', [])]
        stores = [self.W, *self.O, getattr(self, 'N', {})]
        for st in stores:
            for w in st.values():
                roots.extend(w)
        for rec in self.calls_rec:
            for w in rec['in'].values():
                roots.extend(w)
        for kind, exp in self.expected:
            for w in exp.values():
                roots.extend(w)
        nb = B.nodes()
        f = B.gc(roots)
        if MEMLOG:
            import resource
            self.log(f'    [gc] узлов {nb} -> {B.nodes()}, кэш conj {len(self.conj_cache)}, '
                     f'rss_max {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1000000} МБ, оп {self.stats["ops"]}',
                     flush=True)
        self.val = [f(e) for e in self.val]
        self.A = f(self.A) if self.A is not None else None
        self.conj = [(f(e), s_) for e, s_ in self.conj]
        self.conj_cache = {}
        for st in stores:
            for k in list(st):
                st[k] = [f(e) for e in st[k]]
        for rec in self.calls_rec:
            for k in list(rec['in']):
                rec['in'][k] = [f(e) for e in rec['in'][k]]
        self.expected = [(kind, {k: [f(e) for e in w] for k, w in exp.items()}) for kind, exp in self.expected]
        self.preds = [(nm, f(e)) for nm, e in getattr(self, 'preds', [])]
        self.stack_saved = [[f(e) for e in w] for w in getattr(self, 'stack_saved', [])]
        self.stats['gc'] += 1
        self.gc_limit = 2 * B.nodes()            # следующая сборка — когда узлов станет вдвое больше живых
        if B.nodes() > self.max_nodes:
            raise MemoryError(f'BDD nodes {B.nodes()} > {self.max_nodes}')

    # ---------------------------------------------------------------- main
    def run(self, split_flags=True, only_flags=None):
        """Разбор по значениям однобитных входов (флагов): для каждого набора — отдельное доказательство
        с константами вместо флагов (свой менеджер BDD). Наборы покрывают все входы; пустые — «vacuous»."""
        t0 = time.time()
        if check_gadget_lemmas() is not True:
            raise Fail('gadget lemmas failed')
        flags = [name for name, ws in self.prog['inputs'].items() if len(ws) == 1] if split_flags else []
        self.cases = []
        tot = dict(ops=0, ands=0, erases=0, checks_nontrivial=0, gc=0, max_bdd_nodes=0)
        masks = range(1 << len(flags))
        if _os.environ.get('EW_REVERSE_CASES'):      # для канареек: сначала ветви с ctrl=1 (порча чаще видна там)
            masks = reversed(masks)
        for mask in masks:
            asg = {f: (mask >> i) & 1 for i, f in enumerate(flags)}
            if only_flags is not None and asg not in only_flags:
                continue                       # ограниченный прогон: только перечисленные наборы флагов (это сужает утверждение)
            self.B = BDD()
            self.stats = dict(ops=0, ands=0, erases=0, opaque=0, gc=0, checks_nontrivial=0, simp=0)
            try:
                for st, sp in self._run_case(asg):
                    self.cases.append(dict(flags=asg, splits=sp, status=st, seconds=round(time.time() - t0, 1),
                                           bdd_nodes=self.B.nodes()))
                    self.log(f'  случай {asg} {sp}: {st} ({round(time.time() - t0, 1)} с, узлов {self.B.nodes()})',
                             flush=True)
            finally:
                bv.SIMP = None
            for k in ('ops', 'ands', 'erases', 'checks_nontrivial', 'gc'):
                tot[k] += self.stats[k]
            tot['max_bdd_nodes'] = max(tot['max_bdd_nodes'], self.B.nodes())
        if all(c['status'] == 'vacuous' for c in self.cases):
            raise Fail('все случаи пусты (доказательство было бы пустым)')
        self.stats = dict(tot, cases=len(self.cases), vacuous=sum(c['status'] == 'vacuous' for c in self.cases),
                          seconds=round(time.time() - t0, 2), nvars=self.nvars, nq=self.prog['N'])
        return 'PROVED'

    def simp(self, f):
        """f постоянна при текущем A? -> константа, иначе f (замена значения, равного на A)."""
        if f == ONE or f == ZERO:
            return f
        self.stats['simp'] += 1
        if self.holds(f):
            return ZERO
        if self.holds(f ^ 1):
            return ONE
        return f

    def _run_case(self, asg):
        B, n, p = self.B, self.n, self.p
        self._alloc()
        for f, v in asg.items():
            self.W[f] = [ONE if v else ZERO]
        # сначала только диапазоны входных слов и выходов ящиков; полное предусловие — в листе разбиения (с SIMP)
        conj = [bv.lt_const(B, w, p) for name, w in self.W.items() if len(w) == n and name in self.spec_range_words()]
        for c, leaf in enumerate(self.fwd_leaves):     # диапазоны выходов ящиков (контракт C/D)
            key = 'target' if leaf['kind'] == 'MM' else 'x'
            conj.append(bv.lt_const(B, self.O[c][key], p))
        fp = getattr(self.spec, 'flag_prune', None)
        if fp is not None and fp(B, self.W) == ZERO:
            yield 'vacuous', {}
            return
        self.A = None                      # полное A не строится (дорого); при отказе — лениво (full_A)
        self.conj = [(e, frozenset(B.support(e))) for e in conj if e != ONE]
        self.conj_cache = {}
        splits = getattr(self.spec, 'splits', None)
        self.pred_thunks = list(splits(B, self.W, self.O)) if splits else []
        self.preds = [(t[0], ONE) for t in self.pred_thunks]
        self.stack_saved = []
        self.calls_rec, self.expected, self.val = [], [], []
        yield from self._split(0, {})

    def _split(self, i, chosen):
        """Разбиение по предикатам спецификации: f=0 и f=1 (если f не постоянна при A). Всё состояние — в self,
        чтобы сборка мусора BDD переименовывала его (A, части, предикаты, стек)."""
        B = self.B
        if i == len(self.preds):
            yield self._exec_case(), dict(chosen)
            return
        bv.SIMP = self.simp
        entry = self.pred_thunks[i]
        nm, th = entry[0], entry[1]
        subst = entry[2] if len(entry) > 2 else None    # слово, которое в ветви «=1» заменяется нулями (оно = 0)
        self.preds[i] = (nm, th(B, self.W, self.O) if callable(th) else th)   # вычисляется при уже выбранных ветвях
        for cval in (0, 1):
            f = self.simp(self.preds[i][1])
            if f in (ONE, ZERO) and (f == ONE) != bool(cval):
                continue
            pushed = False
            if subst and cval == 1:
                self.stack_saved.append(self.W[subst])       # переименовывается сборкой мусора
                self.W[subst] = [ZERO] * len(self.W[subst])
                pushed = True
            c = ONE if cval else ZERO
            mark = len(self.conj)
            e = B.XOR(f, c) ^ 1
            label = dict(chosen, **{self.preds[i][0]: cval})
            if e == ZERO or (e != ONE and self.holds(e)):      # ветвь несовместна с уже принятыми условиями
                if pushed:
                    self.W[subst] = self.stack_saved.pop()
                yield 'vacuous', label
                continue
            if e != ONE:
                self.conj.append((e, frozenset(B.support(e))))
                self.conj_cache = {}
            yield from self._split(i + 1, label)
            del self.conj[mark:]
            self.conj_cache = {}
            if pushed:
                self.W[subst] = self.stack_saved.pop()

    def spec_range_words(self):
        return getattr(self.spec, 'range_words', ())

    def _exec_case(self):
        B, n, p = self.B, self.n, self.p
        bv.SIMP = self.simp
        self.val, self.calls_rec, self.expected = [], [], []
        self.maybe_gc()
        pre = self.spec.pre(B, self.W, self.O)            # полное предусловие (при выбранных ветвях разбиения)
        mark0 = len(self.conj)
        try:
            for e in (list(pre) if isinstance(pre, (list, tuple)) else [pre]):
                e = self.simp(e)
                if e == ZERO or (e != ONE and self.holds(e)):
                    return 'vacuous'
                if e != ONE:
                    self.conj.append((e, frozenset(B.support(e))))
                    self.conj_cache = {}
            self.A = None
            return self._exec_body()
        finally:
            del self.conj[mark0:]
            self.conj_cache = {}

    def _exec_body(self):
        B, n, p = self.B, self.n, self.p
        N = self.prog['N']
        self.val = [ZERO] * N
        for name, ws in self.prog['inputs'].items():
            for q, e in zip(ws, self.W[name]):
                self.val[q] = e
        self.internal = bytearray(N)    # провода-кандидаты на упрощение: внутренние и однобитные регистры
        io = {q for ws in list(self.prog['inputs'].values()) + list(self.prog['outputs'].values())
              if len(ws) > 1 for q in ws}
        for q in range(N):
            if q not in io:
                self.internal[q] = 1
        self.expected = self.spec.calls(B, self.W, self.O)
        if [k for k, _ in self.expected] != [l['kind'] for l in self.fwd_leaves]:
            raise Fail('список прямых ящиков не совпал со спецификацией',
                       got=[l['kind'] for l in self.fwd_leaves], want=[k for k, _ in self.expected])
        self.calls_rec, fwd_i = [], 0
        val = self.val
        internal = self.internal
        cps = {}
        for pos, wires, fn, name in (self.spec.checkpoints(self.prog) if hasattr(self.spec, 'checkpoints') else []):
            cps.setdefault(pos, []).append((wires, fn, name))
        for item in self.prog['items']:
            if item[0] == 'opaque' and item[2] in cps:
                raise Fail('контрольная точка на ящике не поддержана')
            if item[0] == 'lifted':
                ops, g0 = item[1], item[2]
                for j in range(0, len(ops), 5):
                    if cps and (g0 + j // 5 - 1) in cps:
                        self._checkpoint(cps.pop(g0 + j // 5 - 1))
                        val = self.val
                    k, a, b, t, l = ops[j], ops[j + 1], ops[j + 2], ops[j + 3], ops[j + 4]
                    if k == 1:
                        val[b] = B.XOR(val[b], val[a])
                    elif k == 0:
                        val[a] ^= 1
                    elif k == 2:
                        self.stats['ands'] += 1
                        if internal[a] and val[a] > 1:
                            val[a] = self.simp(val[a])
                        if internal[b] and val[b] > 1:
                            val[b] = self.simp(val[b])
                        if not self.holds(val[t]):
                            self.fail('and: цель не 0', val[t], op=self.stats['ops'], wire=t)
                        val[t] = B.AND(val[a], val[b])
                    else:
                        self.stats['erases'] += 1
                        ab = B.AND(val[a], val[b])
                        if l >= 0:
                            ab = B.XOR(ab, val[l])
                        d = B.XOR(val[t], ab)
                        if not self.holds(d):
                            self.fail('erase: цель не равна a*b', d, op=self.stats['ops'], wire=t)
                        val[t] = ZERO
                    self.stats['ops'] += 1
                    if (self.stats['ops'] & 63) == 0:
                        self.maybe_gc()
                        val = self.val
                last = g0 + len(ops) // 5 - 1
                if last in cps:
                    self._checkpoint(cps.pop(last))
                    val = self.val
            else:
                leaf = item[1]
                self.stats['opaque'] += 1
                fwd_i = self._opaque(leaf, fwd_i)
                val = self.val
            if (self.stats['ops'] & 0xFFFF) == 0:
                pass
        if cps:
            raise Fail('контрольные точки не достигнуты', left=sorted(cps))
        # выходы и нули
        outs = self.prog['outputs']
        post = self.spec.post(B, self.W, self.O)
        if set(post) != set(outs):
            raise Fail('выходы спецификации не совпали с сигнатурой', got=sorted(outs), want=sorted(post))
        for name, ws in outs.items():
            if len(post[name]) != len(ws):
                raise Fail('ширина выхода', name=name)
            for i, (q, e) in enumerate(zip(ws, post[name])):
                d = B.XOR(self.val[q], e)
                if not self.holds(d):
                    self.fail(f'выход {name}[{i}] не равен спецификации', d, wire=q)
        outw = {q for ws in outs.values() for q in ws}
        obs = self.spec.observe(B, self.W, self.O) if getattr(self.spec, 'observe_on', False) else {}
        for name, exp in obs.items():
            ws = self.prog['inputs'][name]
            if set(ws) & outw or len(ws) != len(exp):
                raise Fail('observe: провода не освобождаемый вход', name=name)
            for i, (q, e) in enumerate(zip(ws, exp)):
                d = B.XOR(self.val[q], e)
                if not self.holds(d):
                    self.fail(f'observe {name}[{i}] не равен формуле', d, wire=q)
            outw = outw | set(ws)
        for q in range(N):
            if q not in outw and not self.holds(self.val[q]):
                self.fail('провод вне выходов не 0', self.val[q], wire=q)
        if any(not r['consumed'] for r in self.calls_rec):
            raise Fail('прямой ящик без парного обратного')
        return 'PROVED'

    def _checkpoint(self, lst):
        """Контрольная точка спецификации: провода регистра = формуле (доказывается под A); затем значение заменяется
        именованным словом (определение T == формула добавляется к A: свежие переменные, допущение не сужается)
        или самой формулой. Держит BDD «мелкими» (без вложенных приведений)."""
        B = self.B
        for wires, fn, name in lst:
            expr = fn(B, self.W, self.O, self.N)
            if len(expr) != len(wires):
                raise Fail('ширина контрольной точки')
            for i, (q, e) in enumerate(zip(wires, expr)):
                d = B.XOR(self.val[q], e)
                if not self.holds(d):
                    self.fail(f'контрольная точка {name or "expect"}[{i}]: провод не равен формуле', d, wire=q)
            if name:
                c = bv.eq(B, self.N[name], expr)         # одно определение T == формула (одна часть A)
                self.conj.append((c, frozenset(B.support(c))))
                for q, t in zip(wires, self.N[name]):
                    self.val[q] = t
                self.conj_cache = {}
            else:
                for q, e in zip(wires, expr):
                    self.val[q] = e
        self.maybe_gc()

    def _opaque(self, leaf, fwd_i):
        B, p, val = self.B, self.p, self.val
        kind = leaf['kind']
        rin = {r: lsb(ws) for r, ws in leaf['regs_in'].items()}
        rout = {r: lsb(ws) for r, ws in leaf['regs_out'].items()}
        cur = {r: [val[q] for q in ws] for r, ws in rin.items()}
        if kind in ('MM', 'KINV'):
            c = fwd_i
            ekind, exp = self.expected[c]
            if ekind != kind or self.fwd_leaves[c] is not leaf:
                raise Fail('порядок ящиков', c=c)
            argregs = ('x', 'y') if kind == 'MM' else ('x',)
            if set(exp) != set(argregs):
                raise Fail('спецификация ящика: неверные аргументы', c=c)
            for r in argregs:                       # входы = заявленным и в области
                if not self.word_eq_holds(cur[r], exp[r]):
                    self.fail(f'ящик {c} ({kind}): вход {r} не равен заявленному', bv.eq(B, cur[r], exp[r]) ^ 1, c=c)
                lt = bv.lt_const(B, cur[r], p)
                if not self.holds(lt ^ 1):
                    self.fail(f'ящик {c} ({kind}): вход {r} вне области <p', lt ^ 1, c=c)
            O = self.O[c]
            if kind == 'MM':
                if rout['x'] != rin['x'] or rout['y'] != rin['y']:
                    raise Fail('MM: x,y должны быть THRU')
                for r in ('target', 'qrom_indices', 'reduced'):
                    if any(val[q] != ZERO for q in rout[r]):
                        raise Fail('MM: новый провод не ноль')
                    for q, e in zip(rout[r], O[r]):
                        val[q] = e
            else:
                if rout['x'] != rin['x']:
                    raise Fail('KINV: x должен быть THRU')
                for q, e in zip(rout['x'], O['x']):
                    val[q] = e
                for q, e in zip(rout['junk'], O['junk']):
                    if val[q] != ZERO:
                        raise Fail('KINV: новый провод не ноль')
                    val[q] = e
            self.calls_rec.append(dict(kind=kind, c=c, consumed=False, **{'in': {r: cur[r] for r in argregs}}))
            return fwd_i + 1
        # обратные ящики: ровно выходы парного прямого вызова
        fk = kind[:-4]
        cands = [r for r in self.calls_rec if r['kind'] == fk and not r['consumed']]
        for rec in reversed(cands):
            O = self.O[rec['c']]
            if fk == 'MM':
                ok = (self.word_eq_holds(cur['x'], rec['in']['x']) and self.word_eq_holds(cur['y'], rec['in']['y'])
                      and all(self.word_eq_holds(cur[r], O[r]) for r in ('target', 'qrom_indices', 'reduced')))
            else:
                ok = self.word_eq_holds(cur['x'], O['x']) and self.word_eq_holds(cur['junk'], O['junk'])
            if ok:
                break
        else:
            raise Fail(f'{kind}: нет парного прямого вызова с теми же значениями', n_candidates=len(cands))
        rec['consumed'] = True
        if fk == 'MM':
            if rout['x'] != rin['x'] or rout['y'] != rin['y'] or set(rout) != {'x', 'y'}:
                raise Fail('MM_adj: сигнатура')
            for r in ('target', 'qrom_indices', 'reduced'):
                for q in rin[r]:
                    val[q] = ZERO
        else:
            if rout['x'] != rin['x'] or set(rout) != {'x'}:
                raise Fail('KINV_adj: сигнатура')
            for q, e in zip(rin['x'], rec['in']['x']):
                val[q] = e
            for q in rin['junk']:
                val[q] = ZERO
        return fwd_i
