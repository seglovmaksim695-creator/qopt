"""Canaries: single-operation mutants of the REAL lowered Montgomery program (window 0, the sub-nodes of the generic
window, LessThanConstant, the controlled AddK).  Port of c2_canaries_split of round 3.

A mutant is kept only if the independent simulator sees that it differs from the contract (otherwise it may be
equivalent).  Then the real prover (``kprove`` with a negative control, i.e. it must not be fooled by the witness) and
the exact domain guard run on it.  A canary is CAUGHT if NOT (kprove PROVED and domain PASS).  A canary that is not
caught is a false PROVED and is a failure of the checker.

Run in a process whose Qualtran is the patched copy::

    python -m qoptct.foreign.qualtran.montgomery.regression.canaries --n 16 --w 2 --count 12 --out canaries.json
"""
import argparse
import json
import random
import time
from hashlib import sha256

from .. import _core as C
from .. import _split as SP
from .._stages import Sl
from ...bridge.emit import emit_qasm
from ...bridge.strict_proof import prove_program, digest_json
from ...bridge.domain_guard import check_domains


def run_canaries(n, w, count=12, seed=11, budget=600, max_mb=2800, p=None, quiet=False):
    p = p or C.prime_for(n)
    P = C.Params(n, w, p)
    rng = random.Random(seed)
    B = C.build_segs(n, w, p)
    low, tr, top, order, wins, segs = B['low'], B['tr'], B['top'], B['order'], B['wins'], B['segs']
    info = SP.children_info(low, tr, top, wins, order)
    segmain = {s[0]: s for s in segs}
    win = info[1]
    off = win['span'][0]
    wops = segs[1][2]
    subs = SP.build_sub(P, win['children'], Sl(wops, off))
    real = [c for c in win['children'] if c['name'] not in ('Join', 'Split')]

    def child_of(gidx):
        for c in real:
            if c['a'] <= gidx < c['b']: return c['name']
        return ''
    items = {}
    for (label, kind, o, swi, swo) in subs:
        if label == 'qaa':
            a0 = real[w + 1]['a']; names = [child_of(a0 + j) for j in range(len(o))]
        elif label == 'xor': names = ['Xor'] * len(o)
        else: names = ['CAdd'] * len(o)
        items[label] = (kind, list(o), swi, swo, names)
    for lab in ('win0', 'lt', 'add'):
        s = segmain[lab]; items[lab] = (s[1], list(s[2]), s[3], s[4], [lab] * len(s[2]))
    TARGETS = [
        ('targeted:xor_qrom_index', ['xor'], 'Xor', ('cx',), ['rewire_cx', 'del_cx', 'swap_target']),
        ('targeted:qrom_compute_table', ['qaa'], 'QROAMClean', ('cx',), ['del_cx', 'rewire_cx']),
        ('targeted:adjoint_rule_correction', ['qaa'], 'QROAMCleanAdjointWrapper', ('erase',), ['drop_corr']),
        ('targeted:adjoint_rule_table', ['qaa'], 'QROAMCleanAdjointWrapper', ('cx',), ['del_cx', 'rewire_cx']),
        ('targeted:adder_carry_chain', ['qaa'], 'Add', ('and', 'cx'), ['del_and', 'del_cx', 'rewire_cx']),
        ('targeted:adder_erase_correction', ['qaa'], 'Add', ('erase',), ['drop_corr']),
        ('targeted:controlled_add_window_ctrl', ['cadd%d' % j for j in range(w)], 'CAdd', ('cx', 'and'), ['rewire_cx', 'del_and', 'del_cx']),
        ('targeted:lt_constant', ['lt'], 'lt', ('x',), ['del_x']),
        ('targeted:lt_carry', ['lt'], 'lt', ('and', 'cx'), ['del_and', 'del_cx']),
        ('targeted:final_addk', ['add'], 'add', ('cx', 'and', 'x'), ['del_cx', 'del_and', 'del_x', 'rewire_cx']),
        ('targeted:final_addk_correction', ['add'], 'add', ('erase',), ['drop_corr']),
        ('targeted:win0_whole_window', ['win0'], 'win0', ('cx', 'and', 'erase'), ['del_cx', 'del_and', 'rewire_cx', 'drop_corr']),
    ]
    rows = []; tries = 0; tgt = 0; t_all = time.time()
    want_t = count // 2 + 2
    while len(rows) < count and tries < count * 80:
        tries += 1
        targeted = len([r for r in rows if r['class'].startswith('targeted')]) < want_t
        if targeted:
            cls, labs, sub, kinds, muts = TARGETS[tgt % len(TARGETS)]; tgt += 1
        else:
            cls, labs, sub, kinds, muts = 'random', list(items), '', ('x', 'cx', 'and', 'erase'), ['del_cx', 'rewire_cx', 'drop_corr', 'extra_x', 'del_and', 'swap_target', 'del_x']
        lab = rng.choice(labs); kind, ops, wi_g, wo_g, names = items[lab]
        cand = [i for i, o in enumerate(ops) if o[0] in kinds and (sub == '' or sub == names[i] or sub == lab)]
        if not cand: continue
        i = rng.choice(cand); o = ops[i]; mk = rng.choice(muts)
        new = list(ops); pool = sorted({q for oo in ops if oo[0] == 'cx' for q in oo[1:3]} - {o[1], o[2] if len(o) > 2 else -1})
        if mk == 'del_cx' and o[0] == 'cx': del new[i]
        elif mk == 'del_and' and o[0] == 'and': del new[i]
        elif mk == 'del_x' and o[0] == 'x': del new[i]
        elif mk == 'drop_corr' and o[0] == 'erase': new[i] = ('erase', o[1], ())
        elif mk == 'extra_x': new.insert(i, ('x', rng.choice(sorted({q for oo in ops if oo[0] == 'cx' for q in oo[1:3]}))))
        elif mk == 'rewire_cx' and o[0] == 'cx' and pool: new[i] = ('cx', rng.choice(pool), o[2])
        elif mk == 'swap_target' and o[0] == 'cx' and pool: new[i] = ('cx', o[1], rng.choice(pool))
        else: continue
        nqc, cops, wi, wo, mp = C.compact(new, wi_g, wo_g)
        layout = C.make_layout(nqc, wi, wo, cops); qasm = emit_qasm(nqc, cops)
        first = (lab == 'win0')
        if kind == 'win': pre, spec = C.make_win(P, first, wi, wo)
        elif kind == 'lt': pre, spec = C.make_lt(P, wi, wo)
        elif kind == 'add': pre, spec = C.make_add(P, wi, wo)
        else:
            pre = None
            spec = (lambda wi_, wo_, kind_: (lambda S: SP.apply_sub(P, kind_, S.B, S.reg, wi_, wo_)))(wi, wo, kind)
        contract = dict(source='C2 segment MUTANT', kind=kind, n=n, w=w, p=p, mutation=[cls, mk, lab, i])
        try:
            if kind in ('win', 'lt', 'add'): cnt, bad, _fb = C.sim_check(kind, P, nqc, cops, wi, wo, first, rng, 12)
            else: cnt, bad = SP.sim_sub(P, kind, nqc, cops, wi, wo, rng, 24)
        except Exception:
            cnt, bad = 1, 1   # the simulator rejected the program: a violation found independently
        if not bad: continue
        t0 = time.time()
        b0 = sha256(emit_qasm(*C.compact(ops, wi_g, wo_g)[:2]).encode()).hexdigest()
        rec = dict(status='WITNESS', qasm_sha256=sha256(qasm.encode()).hexdigest(), contract_sha256=digest_json(contract), layout_sha256=digest_json(layout),
                   independent_inputs=max(cnt, 1), all_wires_checked=True, phase_checked=True, baseline_qasm_sha256=b0, witness={'sim_bad': bad})
        v = prove_program(qasm, layout, spec, pre=pre, validation=rec, contract=contract, budget=budget, max_mb=max_mb, negative_control=True)
        try: dg = check_domains(qasm, layout, pre=pre, contract=contract, budget=budget, max_mb=max_mb)['status']
        except Exception as e: dg = 'ERROR:' + str(e)[:50]
        caught = not (v['status'] == 'PROVED' and dg == 'PASS')
        row = dict(n=len(rows), **{'class': cls}, mutation=mk, segment=lab, sub=names[i], op_index=i, op=str(o), sim_cases=cnt, sim_bad=bad,
                   kprove=v['status'], domain=dg, caught=caught, invoked=v.get('prover_invoked'), wall=round(time.time() - t0, 1))
        rows.append(row)
        if not quiet:
            print({k_: row[k_] for k_ in ('n', 'class', 'mutation', 'segment', 'kprove', 'domain', 'caught', 'wall')}, flush=True)
    return dict(n=n, w=w, p=p, seed=seed, canaries=rows, false_proved=[r for r in rows if not r['caught']],
                targeted=len([r for r in rows if r['class'].startswith('targeted')]), seconds=round(time.time() - t_all, 1))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--w', type=int, required=True)
    ap.add_argument('--count', type=int, default=12); ap.add_argument('--seed', type=int, default=11)
    ap.add_argument('--out', required=True)
    a = ap.parse_args(argv)
    res = run_canaries(a.n, a.w, a.count, a.seed)
    json.dump(res, open(a.out, 'w'), indent=1)
    print('DONE', len(res['canaries']), 'false PROVED:', len(res['false_proved']), 'targeted', res['targeted'], flush=True)


if __name__ == '__main__':
    main()
