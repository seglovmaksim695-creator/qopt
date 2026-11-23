"""mutants -- wiring mutants for the corrected client (the wire-flow check).

    python -m qoptct.foreign.qualtran.kaliski.mutants N P OUT [--noprepass] [--aud]
Every mutant is built on the REAL composed circuit ``Composed(N, P)``: either the wire map of an iteration is changed (gmaps: the
operations AND the certificate describe the new wires, like the auditor's mutants F1, F2, F4), or only the certificate (marks) or
the list outputs['junk'] is damaged.  For each: the client ``QChain(...).run()`` in mode prove; CHAIN_FAILED with the reason
"wire flow" is expected.
--noprepass: switch the preliminary static pass off (only the check built into _begin/_fin is tested).
--aud: add the auditor's mutants F1-F4 (at n = 256 with the parameters of the audit; at other n their analogues, see ``aud_params``).
Result: OUT (json): per mutant the status and message.
"""
import copy, json, os, sys, time

from . import chain as qkal_chain
from .lower import Composed, OpsView
from .chain import QChain, ProofFailure, build_marks
from qoptct.verify.compose.lift import K_X, K_CX


def aud_params(n):
    """Parameters of the auditor's mutants F1-F4.  At n = 256 exactly those of the audit (AUD results: F1 swaps bits 1 and 255 of v
    in iteration 507; F4 writes t of iteration 269 onto the wire t_262; F2 swaps m and t in iteration 100).  At other n the
    analogues: the last-but-four iteration and the bits 1 and n-1; t of iteration j+7 onto t_j with j = 2n/3 (same relation,
    j+7 < 2n); the swap of m and t in iteration 2n/4."""
    if n == 256:
        return dict(vswap_k=507, vswap_bit=1, t_src=262, t_dst=269, mt_k=100)
    N2 = 2 * n
    j = N2 // 3
    return dict(vswap_k=N2 - 5, vswap_bit=1, t_src=j, t_dst=j + 7, mt_k=N2 // 4)


def with_g(C, k, newmap):
    M = copy.copy(C)
    M.gmaps = list(C.gmaps)
    g = C.gmaps[k].copy()
    for loc, glob_ in newmap.items():
        g[loc] = glob_
    M.gmaps[k] = g
    M.ops = OpsView(M)
    return M


def swap_in_g(C, k, a, b):
    g = C.gmaps[k]
    return with_g(C, k, {a: int(g[b]), b: int(g[a])})


def idx(k, label):
    return 1 + 6 * k + ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'].index(label)


def edit_marks(marks, i, fn):
    mk = [(pos, lab, copy.deepcopy(inf)) for pos, lab, inf in marks]
    fn(mk[i][2])
    return mk


def build(C, aud):
    n, p = C.n, C.p
    P = C.piece
    lsb = lambda reg: list(reversed([int(q) for q in P.ins[reg]]))
    lr, ls, lu, lv = lsb('r'), lsb('s'), lsb('u'), lsb('v')
    lm, lt = int(P.ins['m'][0]), int(P.ins['is_terminal'][0])
    N2 = 2 * n
    ka, kb, kc = N2 // 2, N2 // 5, N2 - 2
    out = []

    def gm(name, what, M):
        out.append((name, what, M, None, M.outputs))

    def cm(name, what, mk_fn):
        out.append((name, what, C, lambda: mk_fn(build_marks(C)), C.outputs))
    # --- own damages (>= 6) ---
    gm('O1 r between iterations', 'gmaps[%d]: bits 1 and 3 of r swapped (input and output of the iteration)' % ka, swap_in_g(C, ka, lr[1], lr[3]))
    cm('O2 r in S4->S5', 'certificate: wires 0 and 1 of r swapped at S5 of iteration %d' % kb,
       lambda mk: [e if i != idx(kb, 'S5') else (e[0], e[1], dict(e[2], ins=dict(e[2]['ins'], r=[e[2]['ins']['r'][1], e[2]['ins']['r'][0]] + e[2]['ins']['r'][2:]))) for i, e in enumerate(mk)])
    cm('O3 s in S5->S6', 'certificate: wires 2 and 5 of s swapped at S6 of iteration %d' % kb,
       lambda mk: [e if i != idx(kb, 'S6') else (e[0], e[1], dict(e[2], ins=dict(e[2]['ins'], s=[(e[2]['ins']['s'][5] if j == 2 else e[2]['ins']['s'][2] if j == 5 else w) for j, w in enumerate(e[2]['ins']['s'])]))) for i, e in enumerate(mk)])
    gm('O4 u between iterations', 'gmaps[%d]: bits 2 and 4 of u swapped' % kb, swap_in_g(C, kb, lu[2], lu[4]))
    gm('O5 v between iterations', 'gmaps[%d]: bits 2 and 3 of v swapped' % kc, swap_in_g(C, kc, lv[2], lv[3]))
    gm('O6 s between iterations', 'gmaps[%d]: bits 0 and 1 of s swapped' % ka, swap_in_g(C, ka, ls[0], ls[1]))
    jk = [int(w) for w in C.outputs['junk']]
    jk2 = list(jk)
    jk2[3], jk2[2 * n + 5] = jk2[2 * n + 5], jk2[3]
    M = copy.copy(C)
    M.outputs = dict(C.outputs, junk=jk2)
    out.append(('O7 order of junk', 'outputs[junk]: positions 3 and 2N+5 swapped', M, None, M.outputs))
    jk3 = list(jk)
    jk3[0], jk3[1] = jk3[1], jk3[0]
    M = copy.copy(C)
    M.outputs = dict(C.outputs, junk=jk3)
    out.append(('O8 order of junk (neighbours)', 'outputs[junk]: positions 0 and 1 swapped', M, None, M.outputs))
    cm('O9 t0 of the end', 'certificate: t0 = t of iteration 1 instead of iteration 0',
       lambda mk: edit_marks(mk, len(mk) - 1, lambda inf: inf.update(t0=mk[idx(1, 'S1')][2]['outs']['is_terminal'][0])))
    gm('O10 "old" m instead of a new one', 'gmaps[%d]: m is written onto the m wire of iteration %d' % (kc, kc - 1),
       with_g(C, kc, {lm: int(C.gmaps[kc - 1][lm])}))
    gm('O11 "old" t instead of a new one', 'gmaps[%d]: t is written onto the t wire of iteration %d' % (kb, kb - 3),
       with_g(C, kb, {lt: int(C.gmaps[kb - 3][lt])}))
    gm('O12 r of iteration 0 occupies a u wire', 'gmaps[0]: r enters on the wire of bit 0 of register u',
       with_g(C, 0, {lr[0]: int(C.gmaps[0][lu[0]])}))
    cm('O13 f of the end', 'certificate: f of the end = the m wire of the last iteration',
       lambda mk: edit_marks(mk, len(mk) - 1, lambda inf: inf.update(f=mk[idx(N2 - 1, 'S6')][2]['outs']['m'][0])))
    cm('O14 u between iterations (certificate only)', 'certificate: wires 0 and 1 of u swapped at S2 of iteration %d' % kc,
       lambda mk: [e if i != idx(kc, 'S2') else (e[0], e[1], dict(e[2], ins=dict(e[2]['ins'], u=[e[2]['ins']['u'][1], e[2]['ins']['u'][0]] + e[2]['ins']['u'][2:]))) for i, e in enumerate(mk)])
    # --- the auditor's mutants ---
    if aud:
        ap = aud_params(n)
        k, i = ap['vswap_k'], ap['vswap_bit']
        g = C.gmaps[k]
        gm('F1 (audit) iteration %d, bits %d and %d of v' % (k, i, n - 1), 'v[%d] <-> v[%d]' % (i, n - 1), with_g(C, k, {lv[i]: int(g[lv[n - 1]]), lv[n - 1]: int(g[lv[i]])}))
        j_, kt = ap['t_src'], ap['t_dst']
        gm('F4 (audit) iteration %d writes t onto the wire t_%d' % (kt, j_), 'gmaps[%d]: t -> wire t_%d' % (kt, j_),
           with_g(C, kt, {lt: int(C.gmaps[j_][lt])}))
        k2 = ap['mt_k']
        g = C.gmaps[k2]
        gm('F2 (audit) iteration %d, m and t swapped' % k2, 'gmaps[%d]: m <-> t' % k2, with_g(C, k2, {lm: int(g[lt]), lt: int(g[lm])}))
        # F3: the branch x=0 of the end is controlled by t_1 instead of t_0 (operations and certificate)
        from .canary import mutant, op_at, seg_ranges
        marks = build_marks(C)
        R = seg_ranges(marks)
        flo, fhi = R[(None, 'fin')]
        t0w = marks[-1][2]['t0']
        t1w = marks[idx(1, 'S1')][2]['outs']['is_terminal'][0]
        pt = {}
        for j in range(flo, fhi):
            op = list(op_at(C, j))
            if t0w in op[1:3] and op[0] != K_X:
                op[1:3] = [t1w if q == t0w else q for q in op[1:3]]
                pt[j] = tuple(op)
        mk = list(marks)
        mk[-1] = (mk[-1][0], mk[-1][1], dict(mk[-1][2], t0=t1w))
        del marks
        out.append(('F3 (audit) the end x=0 is controlled by t_1', '%d operations of fin: t0 -> t1; the certificate too' % len(pt), mutant(C, pt), lambda: mk, C.outputs))
    return out


def run_mutants(n, p, *, noprepass=False, aud=False, outp=None, only=None, verbose=True):
    """Run every mutant through the client; returns the result dictionary (rows with status, reason, flow_reason)."""
    saved = qkal_chain.check_flow
    if noprepass:
        qkal_chain.check_flow = lambda *a, **k: None
    try:
        C = Composed(n, p)
        res = dict(n=n, p=p, prepass=not noprepass, started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), rows=[])
        for name, what, M, marks, outputs in build(C, aud):
            if only and not name.startswith(only):
                continue
            t = time.time()
            M2 = copy.copy(M)
            M2.outputs = outputs
            marks = marks() if callable(marks) else build_marks(M2)
            ch = QChain(M2, lifted=M2.lifted, marks=marks)
            try:
                ch.run()
                st, msg = 'CHAIN_PROVED', ''
            except ProofFailure as e:
                st, msg = 'NOT_PROVED', str(e)[:300]
            row = dict(name=name, what=what, status=st, reason=msg, flow_reason=msg.startswith('wire flow'), sec=round(time.time() - t, 2))
            res['rows'].append(row)
            del ch, marks, M2
            if verbose:
                print(json.dumps(row, ensure_ascii=False), flush=True)
    finally:
        qkal_chain.check_flow = saved
    res['all_rejected'] = all(r['status'] == 'NOT_PROVED' for r in res['rows'])
    res['all_flow_reason'] = all(r['flow_reason'] for r in res['rows'])
    if outp:
        json.dump(res, open(outp, 'w'), indent=1, ensure_ascii=False)
    return res


def main():
    n, p, outp = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    res = run_mutants(n, p, noprepass='--noprepass' in sys.argv, aud='--aud' in sys.argv, outp=outp)
    print('rejected %d of %d; reason "wire flow" for all: %s' % (sum(r['status'] == 'NOT_PROVED' for r in res['rows']), len(res['rows']), res['all_flow_reason']))


if __name__ == '__main__':
    main()
