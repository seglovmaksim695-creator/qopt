"""mont_compose -- пара DirtyOutOfPlaceMontgomeryModMul(w=1) прямой + обратный в макро-виде: n/w прямых окон SingleWindowModMul и
n/w обратных (Adjoint окна) -- макро-листья; одно прямое и одно обратное окно раскрываются мостом C (правило QROAMCleanAdjoint).
Особенность Монтгомери: в обратном окне QROAMClean и QROAMCleanAdjoint ПОМЕНЯЛИСЬ ролями (зеркало загрузки -- снятие, и наоборот),
а снятие -- примитив по правилу C1, не побитовое зеркало загрузки. Поэтому оба блока (w=1, sel = 1 бит) заменяются на ПОБИТНЫЕ
абстрактные операции с проверяемым смыслом:
  LD(sel, t, код) : t ^= b0 ^ sel*(b0^b1)         (блок QROAMClean: только x/cx, управление только sel)
  CLR(sel, t, код): предусловие t == b0 ^ sel*(b0^b1), результат t = 0   (блок QROAMCleanAdjoint по правилу C1: исполнением при sel=0,1)
Зеркало: LD <-> CLR с тем же кодом. Смысл блоков проверяется ПОЛНЫМ перебором sel in {0,1} (исполнение операций моста), цели LD --
свежие нулевые провода. Всё остальное -- обычное зеркало (mirror.py).
Среда: env_h2.sh (исправленная Qualtran H/qlocal, мост C)."""
import sys, json, time
import numpy as np
from qoptct.foreign.qualtran.bridge.qroam_clean_adjoint import chained_atom_rule
from qoptct.foreign.qualtran.bridge.qualtran_layout import lower_bloqgraph, AtomExpansion
from qoptct.foreign.qualtran.bridge.rules import translate
from qoptct.verify.compose.lift import lift
from . import mirror as M
from .mirror import K_X, K_CX, K_AND, K_ERASE, K_LD, K_CLR
from . import pair_compose as PC
from .pair_compose import CErr, need
from .pair_kal import ops_list

BIG = 10 ** 9


def classify(path):
    names = [q.split(':', 1)[1] if ':' in q else q for q in path.split('/')[1:]]
    if any(nm == 'QROAMClean' for nm in names):
        return 'QC'
    if any(nm.startswith('QROAMCleanAdjoint') for nm in names):
        return 'QA'
    return None


def token_ranges(prog):
    """непрерывные диапазоны операций блоков QC/QA по трассе: [(lo, hi, kind, ключ)]"""
    out = []
    for row in prog['trace']:
        kind = classify(row['leaf']['path'])
        if kind is None or row['logical_end'] == row['logical_start']:
            continue
        key = row['leaf']['path'].split('/')[1]
        if out and out[-1][3] == key and out[-1][1] == row['logical_start']:
            out[-1][1] = row['logical_end']
        else:
            need(all(o[3] != key for o in out), 'блок %s встречается дважды' % key)
            out.append([row['logical_start'], row['logical_end'], kind, key])
    return out


def extract_token(kind, tops):
    """-> (sel, {t: (b0, b1)}, внутренние провода). Проверяет смысл блока перебором sel in {0,1}."""
    wset = sorted({w for o in tops for w in M.wires(o)})
    ctrl = {o[1] for o in tops if o[0] == K_CX}
    tgt_written = {o[2] for o in tops if o[0] == K_CX} | {o[3] for o in tops if o[0] in (K_AND, K_ERASE)}
    sel = sorted(ctrl - tgt_written)
    need(len(sel) == 1, 'блок %s: ожидался один sel-провод (w=1), найдено %r' % (kind, sel))
    sel = sel[0]
    if kind == 'QC':
        need(all(o[0] in (K_X, K_CX) for o in tops), 'QC: только x/cx')
        need(all(o[1] == sel for o in tops if o[0] == K_CX), 'QC: cx управляются только sel')
        need(all(o[1] == sel for o in tops if o[0] == K_X), 'QC: x только на sel')
        need(all(o[2] != sel for o in tops if o[0] == K_CX), 'QC: sel не мишень')
        code = {}
        for s in (0, 1):
            st = {w: 0 for w in wset}
            st[sel] = s
            for o in tops:
                if o[0] == K_X:
                    st[o[1]] ^= 1
                else:
                    st[o[2]] ^= st[o[1]]
            need(st[sel] == s, 'QC: sel изменён')
            for w, v in st.items():
                if w != sel:
                    code.setdefault(w, [0, 0])[s] = v
        tab = {w: tuple(c) for w, c in code.items() if tuple(c) != (0, 0)}
        return sel, tab, set()
    # QA: исполнение с записью требуемых значений мишеней erase
    need(all(o[0] in (K_X, K_CX, K_AND, K_ERASE) for o in tops), 'QA: виды вентилей')
    erased = [o[3] for o in tops if o[0] == K_ERASE]
    need(len(erased) == len(set(erased)), 'QA: мишень erase стирается дважды')
    need(all(o[0] != K_ERASE or o[4] < 0 for o in tops), 'QA: поправка l не поддержана')
    read_as_ctrl = {o[1] for o in tops if o[0] == K_CX} | {w for o in tops if o[0] in (K_AND, K_ERASE) for w in (o[1], o[2])}
    need(not (set(erased) & read_as_ctrl), 'QA: стёртый провод читается как управляющий')
    code = {}
    internal = set(wset) - {sel} - set(erased)
    for s in (0, 1):
        st = {w: 0 for w in wset}
        st[sel] = s
        for o in tops:
            k = o[0]
            if k == K_X:
                st[o[1]] ^= 1
            elif k == K_CX:
                st[o[2]] ^= st[o[1]]
            elif k == K_AND:
                need(st[o[3]] == 0, 'QA: цель AND не 0')
                st[o[3]] = st[o[1]] & st[o[2]]
            else:
                code.setdefault(o[3], [0, 0])[s] = st[o[1]] & st[o[2]]
                st[o[3]] = 0
        need(st[sel] == s and all(v == 0 for w, v in st.items() if w != sel), 'QA: провода не вернулись в 0 при sel=%d' % s)
    tab = {w: tuple(c) for w, c in code.items() if tuple(c) != (0, 0)}
    return sel, tab, internal


def abstract_ops(ops, prog, interface_in):
    """заменяет блоки QC/QA абстрактными LD/CLR; проверяет: цели LD -- свежие провода (не интерфейс, первое использование в блоке);
    таблица CLR парного LD в окне совпадает. -> (новые ops, сведения)"""
    toks = token_ranges(prog)
    need(len(toks) == 2 and [t[2] for t in toks] in (['QC', 'QA'], ['QA', 'QC']), 'окно: ожидались блоки QC и QA, есть %r' % [t[2] for t in toks])
    info = []
    new = []
    pos = 0
    tabs = {}
    for lo, hi, kind, key in toks:
        new.extend(ops[pos:lo])
        pos = hi
        sel, tab, internal = extract_token(kind, ops[lo:hi])
        tabs[kind] = (sel, tab)
        absops = []
        for w in sorted(tab):
            code = tab[w][0] * 2 + tab[w][1]
            absops.append((K_LD if kind == 'QC' else K_CLR, sel, w, code, -1))
        new.extend(absops)
        info.append(dict(kind=kind, ops=hi - lo, sel=sel, table=tab, internal=sorted(internal)))
        if kind == 'QC':
            first_use = {}
            for i, o in enumerate(ops):
                for w in M.wires(o):
                    first_use.setdefault(w, i)
            for w in tab:
                need(w not in interface_in, 'LD: цель %d -- провод интерфейса' % w)
                need(first_use[w] >= lo, 'LD: цель %d использована до блока QC' % w)
    new.extend(ops[pos:])
    need(tabs['QC'] == tabs['QA'], 'таблицы QC и QA окна различаются: %r %r' % (tabs['QC'], tabs['QA']))
    return new, info


class PieceW:
    def __init__(self, n, w, p, adjoint):
        from qualtran.bloqs.mod_arithmetic.mod_multiplication import SingleWindowModMul
        self.n, self.w, self.p, self.adjoint = n, w, p, adjoint
        b = SingleWindowModMul(w, n, p)
        self.bloq = b.adjoint() if adjoint else b
        self.low = low = lower_bloqgraph(self.bloq, extra_atom_rule=chained_atom_rule(), max_leaves=BIG, max_nodes=BIG)
        self.prog = prog = translate(low, low['inputs'], low['outputs'])
        self.L = L = lift(prog['qasm'])
        need(L.n == len(prog['ops']), 'lift != translate')
        self.ops_real = ops_list(L)
        self.nsrc = low['nq']
        self.scratch = prog['layout']['bridge_scratch']
        need(self.scratch in ([], [self.nsrc]), 'scratch')
        regs_in, regs_out = low['registers']['inputs'], low['registers']['outputs']
        self.in_flat = {d['name']: [q for e in d['elements'] for q in e['wires_msb']] for d in regs_in}
        self.out_flat = {d['name']: [q for e in d['elements'] for q in e['wires_msb']] for d in regs_out}
        self.ins, self.outs = self.in_flat, self.out_flat
        self.nin = sum(len(v) for v in self.in_flat.values())
        need(sorted(q for v in self.in_flat.values() for q in v) == list(range(self.nin)), 'входы не 0..nin-1')
        self.ops, self.tokens = abstract_ops(self.ops_real, prog, set(range(self.nin)))
        self.n_meas = L.nmeas


class MontPair(PC.PairComposed):
    def __init__(self, n, w, p, pieces=None):
        from qualtran import Register, QAny, Adjoint
        from qualtran.bloqs.mod_arithmetic.mod_multiplication import SingleWindowModMul
        from .pair_mont import pair_bloq
        self.n, self.w, self.p = n, w, p
        self.Pf, self.Pb = pieces or (PieceW(n, w, p, False), PieceW(n, w, p, True))
        Pf, Pb = self.Pf, self.Pb
        win = SingleWindowModMul(w, n, p)
        base = chained_atom_rule()
        macros = []

        def rule(bloq, inputs, allocate, path):
            kind, P = None, None
            if type(bloq) is SingleWindowModMul:
                kind, P = 'f', Pf
                need(bloq == Pf.bloq, 'прямое окно отличается')
            elif type(bloq) is Adjoint and bloq.subbloq == win:
                kind, P = 'b', Pb
                need(bloq == Pb.bloq, 'обратное окно отличается')
            if kind is None:
                return base(bloq, inputs, allocate, path)
            nanc = P.nsrc - P.nin
            anc = allocate([Register('anc', QAny(nanc))], path + '/macro')['anc'].ravel() if nanc else []
            g = {}
            for d in P.low['registers']['inputs']:
                loc = [q for e in d['elements'] for q in e['wires_msb']]
                glob = [int(q) for q in inputs[d['name']].flat]
                need(len(loc) == len(glob), 'макро: ширина регистра %s' % d['name'])
                g.update(zip(loc, glob))
            for j in range(nanc):
                g[P.nin + j] = int(anc[j])
            need(len(g) == P.nsrc and len(set(g.values())) == P.nsrc, 'макро: отображение не инъективно')
            outs = {}
            for d in P.low['registers']['outputs']:
                loc = [q for e in d['elements'] for q in e['wires_msb']]
                outs[d['name']] = np.asarray([g[q] for q in loc], dtype=object).reshape(inputs[d['name']].shape)
            tw = {g[q] for o in P.ops_real for q in M.wires(o) if q in g}
            macros.append(dict(path=path, g=g, kind=kind, tw=tw))
            leaf = dict(name='macro', qubits=sorted(g.values()), path=path, macro=len(macros) - 1)
            return AtomExpansion(outs, [leaf], tuple(g[q] for q in P.low['retired_wires']))
        self.bloq = pair_bloq(n, w, p)
        self.top = top = lower_bloqgraph(self.bloq, extra_atom_rule=rule, max_leaves=BIG, max_nodes=BIG)
        nwin = n // w
        need(len(macros) == 2 * nwin, 'ожидалось %d оконных узлов, есть %d' % (2 * nwin, len(macros)))
        self.nwin = nwin
        leaves, mpos = [], []
        for lf in top['leaves']:
            if lf['name'] == 'macro':
                mpos.append((len(leaves), lf['macro']))
            else:
                leaves.append(lf)
        need([m for _, m in mpos] == list(range(2 * nwin)), 'порядок макро')
        freew = sorted({q for v in top['inputs'].values() for q in v} |
                       {q for lf in top['leaves'] if lf['name'] == 'macro' for q in lf['qubits']})
        tprog = translate(dict(nq=top['nq'], leaves=leaves), {'all': freew}, {'all': freew})
        TL = lift(tprog['qasm'])
        need(TL.n == len(tprog['ops']), 'top lift')
        starts = [row['logical_start'] for row in tprog['trace']] + [len(tprog['ops'])]
        self.nsrc = top['nq']
        uses_scr = bool(tprog['layout']['bridge_scratch']) or bool(Pf.scratch) or bool(Pb.scratch)
        need(tprog['layout']['bridge_scratch'] in ([], [self.nsrc]), 'top scratch')
        self.nq = self.nsrc + (1 if uses_scr else 0)
        self.scratch = [self.nsrc] if uses_scr else []
        self.top_ops = ops_list(TL)
        self.top_cut = [starts[i] for i, _ in mpos]
        self.macros = macros
        halves = [row['leaf']['path'].split('/')[1].split(':')[0] for row in tprog['trace']]
        need(set(halves) <= {'0', '1'}, 'halves')
        first_b = next((i for i, h in enumerate(halves) if h == '1'), len(halves))
        need(all(h == '0' for h in halves[:first_b]) and all(h == '1' for h in halves[first_b:]), 'halves not contiguous')
        self.F_end = starts[first_b]
        self.top_low, self.tprog = top, tprog
        for m in macros:
            m['half'] = 'f' if m['path'].split('/')[1].split(':')[0] == '0' else 'b'
            need((m['half'] == 'f') == (m['kind'] == 'f'), 'half/kind')
        self.nmac = len(macros)

    def gfull(self, m):
        P = self.Pf if m['kind'] == 'f' else self.Pb
        g = np.zeros(P.nsrc + (1 if P.scratch else 0), dtype=np.int64)
        for a, b in m['g'].items():
            g[a] = b
        if P.scratch:
            g[P.nsrc] = self.scratch[0]
        return g

    def materialize(self):
        out, pos = [], 0
        for k, cut in enumerate(self.top_cut):
            out.extend(self.top_ops[pos:cut])
            pos = cut
            m = self.macros[k]
            P = self.Pf if m['kind'] == 'f' else self.Pb
            g = self.gfull(m)
            for (kk, a, b, t, l) in P.ops_real:
                if kk == 0:
                    out.append((0, int(g[a]), 0, 0, -1))
                elif kk == 1:
                    out.append((1, int(g[a]), int(g[b]), 0, -1))
                else:
                    out.append((kk, int(g[a]), int(g[b]), int(g[t]), int(g[l]) if l >= 0 else -1))
        out.extend(self.top_ops[pos:])
        return out


def check_pair_mont(n, w, p, compare_full=False):
    t0 = time.time()
    C = MontPair(n, w, p)
    r0 = dict(n=n, w=w, p=p, piece_fwd_ops=len(C.Pf.ops_real), piece_bwd_ops=len(C.Pb.ops_real), top_ops=len(C.top_ops), macros=C.nmac,
              tokens_fwd=C.Pf.tokens, tokens_bwd=C.Pb.tokens)
    if compare_full:
        from .pair_mont import lower_pair
        low, tr = lower_pair(n, w, p)
        L = lift(tr['qasm'])
        r0['equals_full_lowering'] = (ops_list(L) == C.materialize()) and L.nq == C.nq
    r = PC.check_composed(C)
    r.update(r0)
    r['sec'] = round(time.time() - t0, 1)
    return r


if __name__ == '__main__':
    n, w, p = int(sys.argv[1]), int(sys.argv[2]), int(eval(sys.argv[3], {'__builtins__': {}}))
    cmp_full = len(sys.argv) > 4 and sys.argv[4] == 'full'
    try:
        r = check_pair_mont(n, w, p, compare_full=cmp_full)
    except CErr as e:
        r = dict(n=n, w=w, p=p, mirror_ok=False, error=str(e))
    print(json.dumps(r, ensure_ascii=False, default=str))
