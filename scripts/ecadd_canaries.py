# SPDX-License-Identifier: Apache-2.0
"""Canaries of the compositional proof of the point addition (qoptct.verify.compose.ecadd).

    python scripts/ecadd_canaries.py small [--count N] [--seed S] [--curves n4,n5,...] [-o out.json]
        single-op mutants of small circuits (all segment kinds, both halves of the compute/uncompute blocks, segment
        boundaries, doctored certificates); the verdict of the CHAIN ALONE (no sampled simulation) is compared with
        the truth = exhaustive classical run over all points of the curve and O.  A wrong circuit must never be
        PROVED.
    python scripts/ecadd_canaries.py big QASM HINTS --curve secp256k1|p256 [--only LABELS] [-o out.jsonl]
        targeted in-memory mutants of a full-size (n = 256) circuit lifted from its OpenQASM text, one per run of the
        verifier; none may be PROVED.
"""
import argparse
import json
import random
import sys
import time
from array import array

from qoptct.crypto import ecpoint as E
from qoptct.crypto._ecp.ec import small_curve
from qoptct.verify.compose import lift, lift_circuit, verify_ecadd
from qoptct.verify.compose.ecadd import ECChain, _circuit_run, curve_points, check_parameters
from qoptct.verify.compose.lift import Lifted, K_X, K_CX, K_AND, K_ERASE
from qoptct.verify.compose.walker import ProofFailure

KINDS = ('cx_flip', 'cx_target', 'cx_control', 'cx_to_x', 'x_move', 'and_ctrl', 'erase_ctrl', 'erase_lin', 'drop',
         'insert_x', 'insert_cx', 'swap_adjacent', 'and_to_erase', 'erase_to_and')


APPLICABLE = {K_X: ('x_move', 'drop', 'insert_x', 'insert_cx', 'swap_adjacent'),
              K_CX: ('cx_flip', 'cx_target', 'cx_control', 'cx_to_x', 'drop', 'insert_x', 'insert_cx', 'swap_adjacent'),
              K_AND: ('and_ctrl', 'and_ctrl', 'and_to_erase', 'insert_x', 'swap_adjacent'),
              K_ERASE: ('erase_ctrl', 'erase_lin', 'erase_lin', 'erase_to_and', 'insert_cx', 'swap_adjacent')}


def mutate(L, kind, i, rnd):
    """-> (new op array, delta of the op count) or None if the mutation does not apply at op i."""
    ops = L.ops
    o = list(ops[5 * i:5 * i + 5])
    kd, a, b, t, l = o
    nq = L.nq

    def other(*avoid):
        for _ in range(50):
            w = rnd.choice([a, b, t, a + 1, a - 1, b + 1, b - 1, t + 1, t - 1, rnd.randrange(nq)])
            if 0 <= w < nq and w not in avoid:
                return w
        return None
    new = None
    ins = None
    if kind == 'cx_flip' and kd == K_CX:
        new = [K_CX, b, a, 0, -1]
    elif kind == 'cx_target' and kd == K_CX:
        w = other(a, b)
        new = None if w is None else [K_CX, a, w, 0, -1]
    elif kind == 'cx_control' and kd == K_CX:
        w = other(a, b)
        new = None if w is None else [K_CX, w, b, 0, -1]
    elif kind == 'cx_to_x' and kd == K_CX:
        new = [K_X, b, 0, 0, -1]
    elif kind == 'x_move' and kd == K_X:
        w = other(a)
        new = None if w is None else [K_X, w, 0, 0, -1]
    elif kind == 'and_ctrl' and kd == K_AND:
        w = other(a, b, t)
        new = None if w is None else [K_AND, w, b, t, -1]
    elif kind == 'erase_ctrl' and kd == K_ERASE:
        w = other(a, b, t, l)
        new = None if w is None else [K_ERASE, a, w, t, l]
    elif kind == 'erase_lin' and kd == K_ERASE:
        if l >= 0:
            new = [K_ERASE, a, b, t, -1]
        else:
            w = other(a, b, t)
            new = None if w is None else [K_ERASE, a, b, t, w]
    elif kind == 'and_to_erase' and kd == K_AND:
        new = [K_ERASE, a, b, t, -1]
    elif kind == 'erase_to_and' and kd == K_ERASE:
        new = [K_AND, a, b, t, -1]
    elif kind == 'drop' and kd in (K_X, K_CX):
        out = array('i', ops[:5 * i])
        out.extend(ops[5 * i + 5:])
        return out, -1
    elif kind == 'insert_x':
        ins = [K_X, rnd.choice([a, b if kd != K_X else a, rnd.randrange(nq)]), 0, 0, -1]
    elif kind == 'insert_cx':
        w = other(a)
        ins = None if w is None else [K_CX, w, a, 0, -1]
        if ins is None:
            return None
    elif kind == 'swap_adjacent' and 5 * i + 10 <= len(ops):
        o2 = list(ops[5 * i + 5:5 * i + 10])
        if o2 == o:
            return None
        out = array('i', ops)
        out[5 * i:5 * i + 5] = array('i', o2)
        out[5 * i + 5:5 * i + 10] = array('i', o)
        return out, 0
    if ins is not None:
        out = array('i', ops[:5 * i])
        out.extend(ins)
        out.extend(ops[5 * i:])
        return out, 1
    if new is None or new == o:
        return None
    out = array('i', ops)
    out[5 * i:5 * i + 5] = array('i', new)
    return out, 0


def commute(o1, o2):
    """True if the two lifted ops act on disjoint wires or only share wires that both only read."""
    def rw(o):
        kd, a, b, t, l = o
        if kd == K_X:
            return {a}, {a}
        if kd == K_CX:
            return {b}, {a, b}
        return {t}, {a, b, t} | ({l} if l >= 0 else set())
    w1, a1 = rw(o1)
    w2, a2 = rw(o2)
    return not (w1 & a2) and not (w2 & a1)


def shift_marks(hints, i, delta):
    if not delta:
        return hints
    h = dict(hints)
    h['marks'] = [[pos + delta if pos > i else pos, lab, info] for pos, lab, info in hints['marks']]
    return h


def doctor(hints, rnd):
    """a doctored certificate: one mark moved, or two wire names of a mark exchanged / changed."""
    h = dict(hints)
    marks = [[pos, lab, json.loads(json.dumps(info))] for pos, lab, info in hints['marks']]
    j = rnd.randrange(len(marks))
    what = rnd.randrange(3)
    if what == 0:
        marks[j][0] += rnd.choice([-2, -1, 1, 2])
    else:
        info = marks[j][2]
        keys = [k for k, v in info.items() if isinstance(v, list) and len(v) >= 2 and all(isinstance(x, int) for x in v)]
        ints = [k for k, v in info.items() if isinstance(v, int) and not isinstance(v, bool)]
        if keys and what == 1:
            k = rnd.choice(keys)
            u, v = rnd.sample(range(len(info[k])), 2)
            info[k][u], info[k][v] = info[k][v], info[k][u]
        elif ints:
            k = rnd.choice(ints)
            info[k] += rnd.choice([-1, 1])
        else:
            marks[j][0] += 1
    h['marks'] = marks
    return h


def chain_verdict(p, a, b, Q, hints, L):
    """the verdict of the chain alone: True = PROVED (no sampled simulation, no refuter)."""
    try:
        ECChain(p, a, b, Q, tuple(hints['D']), hints, lifted=L, budget=120.0).run()
        return True, ''
    except ProofFailure as e:
        return False, str(e)[:160]
    except (KeyError, TypeError, IndexError, AttributeError, ValueError) as e:
        return False, 'certificate: %s' % type(e).__name__


def positions(hints, nops, rnd, count):
    """op positions stratified over the segments: boundaries (first / last op of every segment) and random ops."""
    marks = hints['marks']
    segs = []
    prev = 0
    for pos, lab, info in marks:
        if pos > prev:
            segs.append((lab, prev, pos))
        prev = pos
    if prev < nops:
        segs.append(('tail', prev, nops))
    bylab = {}
    for s in segs:
        bylab.setdefault(s[0], []).append(s)
    out = []
    labs = sorted(bylab)
    while len(out) < count:
        lab = rnd.choice(labs)
        _, lo, hi = rnd.choice(bylab[lab])
        r = rnd.random()
        i = lo if r < 0.15 else (hi - 1 if r < 0.3 else rnd.randrange(lo, hi))
        out.append((lab, i))
    return out


def small_curves(ns, seed=1000):
    out = []
    for n in ns:
        rng = random.Random(seed + n)
        got = []
        while len(got) < 2:
            cv = small_curve(n, rng, prime_order=True, min_points=6,
                             pfilter=lambda q: (1 << n) - q < (1 << (n // 2 + 1)))
            if all((cv['p'], cv['a'], cv['b']) != (o['p'], o['a'], o['b']) for o in got):
                got.append(cv)
        out += got
    return out


def run_small(count, seed, ns, levels=1, out=None, progress=True):
    rnd = random.Random(seed)
    stats = dict(mutants=0, wrong=0, wrong_proved=0, correct=0, correct_proved=0, by_label={}, by_kind={},
                 doctored=0, doctored_wrong_proved=0, base_proved=0, base=0)
    false_proved = []
    t0 = time.time()
    for cv in small_curves(ns):
        p, a, b = cv['p'], cv['a'], cv['b']
        n = p.bit_length()
        allp = [None] + curve_points(p, a, b)
        for Q in rnd.sample(cv['points'], 2):
            C = E.ec_point_add(p, a, b, Q, karatsuba_levels=levels)
            hints = json.loads(json.dumps(C.hints()))
            assert check_parameters(p, a, b, Q, tuple(hints['D'])) is None
            L0 = lift_circuit(C.circ)
            ok, why = chain_verdict(p, a, b, Q, hints, L0)
            stats['base'] += 1
            stats['base_proved'] += ok
            assert _circuit_run(L0, n, allp, p, a, Q) is None
            if not ok:
                print('BASE NOT PROVED', p, a, b, Q, why)
            for lab, i in positions(hints, L0.n, rnd, count):
                if rnd.random() < 0.5:                       # half of the mutants on AND / ERASE gadgets
                    j = i
                    while j < L0.n and L0.ops[5 * j] not in (K_AND, K_ERASE):
                        j += 1
                    if j < L0.n and j - i < 12:
                        i = j
                kind = rnd.choice(APPLICABLE[L0.ops[5 * i]])
                m = mutate(L0, kind, i, rnd)
                if m is None:
                    continue
                L = Lifted(L0.nq, m[0], 0)
                h = shift_marks(hints, i, m[1])
                doc = rnd.random() < 0.15
                if doc:
                    h = doctor(h, rnd)
                truth = _circuit_run(L, n, allp, p, a, Q) is None
                ok, why = chain_verdict(p, a, b, Q, h, L)
                stats['mutants'] += 1
                stats['doctored'] += doc
                key = 'correct' if truth else 'wrong'
                stats[key] += 1
                stats[key + '_proved'] += ok
                for d, k in ((stats['by_label'], lab), (stats['by_kind'], kind)):
                    e = d.setdefault(k, [0, 0, 0, 0])            # wrong, wrong proved, correct, correct proved
                    e[0 if not truth else 2] += 1
                    e[1 if not truth else 3] += ok
                if ok and not truth:
                    stats['doctored_wrong_proved'] += doc
                    false_proved.append(dict(p=p, a=a, b=b, Q=list(Q), op=i, label=lab, kind=kind, doctored=doc))
                    print('FALSE PROVED', false_proved[-1], flush=True)
            if progress:
                print('p=%d a=%d b=%d Q=%s: %d mutants so far, wrong %d (proved %d), correct %d (proved %d)  %.0f s'
                      % (p, a, b, Q, stats['mutants'], stats['wrong'], stats['wrong_proved'], stats['correct'],
                         stats['correct_proved'], time.time() - t0), flush=True)
    stats['false_proved'] = false_proved
    stats['time'] = round(time.time() - t0, 1)
    if out:
        with open(out, 'w') as f:
            json.dump(stats, f, indent=1)
    return stats


# ---------------------------------------------------------------------------------------------- n = 256
def pick_big(hints, L, rnd, only=None):
    """targeted mutants: (label of the segment, op index, kind) spread over all segment kinds, first / last ops of
    segments (joints) and both halves of the compute / uncompute blocks."""
    marks = hints['marks']
    segs = []
    prev = 0
    for pos, lab, info in marks:
        if pos > prev:
            segs.append((lab, prev, pos))
        prev = pos
    bylab = {}
    for s in segs:
        bylab.setdefault(s[0], []).append(s)
    plan = []
    want = {'w_pre': 4, 'addc': 3, 'cm_const': 1, 'cm_term': 2, 'fold': 4, 'init': 1, 'cmp': 3, 'swap': 3, 'sub': 3,
            'add': 3, 'last': 1, 'mul': 4, 'sq': 2, 'apply': 5, 'neg': 2, 'mirror': 6, 'w_post': 5, 'sp_const': 1,
            'sp_term': 1, 'sp_q': 1}
    for lab, cnt in want.items():
        if lab not in bylab or (only and lab not in only):
            continue
        for c in range(cnt):
            _, lo, hi = rnd.choice(bylab[lab])
            where = c % 3
            i = lo if where == 1 else (hi - 1 if where == 2 else rnd.randrange(lo, hi))
            plan.append((lab, i))
    return plan


def run_big(qasm, hints_path, curve, only=None, out=None, seed=5):
    cv = E.CURVES[{'secp256k1': 'secp256k1', 'p256': 'P-256'}[curve]]
    p, a, b = cv['p'], cv['a'], cv['b']
    hints = json.load(open(hints_path))
    Q = tuple(hints['Q'])
    t0 = time.time()
    L = lift(qasm)
    print('lifted %d ops in %.0f s' % (L.n, time.time() - t0), flush=True)
    rnd = random.Random(seed)
    res = []
    kinds = ('cx_flip', 'cx_target', 'cx_to_x', 'x_move', 'and_ctrl', 'erase_ctrl', 'erase_lin', 'swap_adjacent',
             'and_to_erase', 'cx_control')
    for lab, i in pick_big(hints, L, rnd, only):
        m = None
        for _ in range(40):
            kind = rnd.choice(kinds)
            j = min(L.n - 2, i)
            saved = array('i', L.ops[5 * j:5 * j + 10])
            small = Lifted(L.nq, array('i', saved), 0)
            if kind == 'swap_adjacent' and commute(list(saved[:5]), list(saved[5:])):
                continue                                     # an equivalent circuit, not a canary
            m = mutate(small, kind, 0, rnd)
            if m is not None and m[1] == 0:
                break
            m = None
        if m is None:
            continue
        L.ops[5 * j:5 * j + 10] = m[0]
        t1 = time.time()
        r = verify_ecadd(None, p, a, b, Q, hints=hints, lifted=L, refute_samples=64, test_samples=64)
        L.ops[5 * j:5 * j + 10] = saved
        rec = dict(curve=curve, label=lab, op=j, kind=kind, status=r['status'], stage=r.get('stage'),
                   detail=str(r.get('detail'))[:200], counterexample=str(r.get('counterexample'))[:120],
                   time=round(time.time() - t1, 1))
        res.append(rec)
        print(json.dumps(rec), flush=True)
        if out:
            with open(out, 'a') as f:
                f.write(json.dumps(rec) + '\n')
    cnt = {}
    for r in res:
        cnt[r['status']] = cnt.get(r['status'], 0) + 1
    print('SUMMARY', curve, len(res), cnt, 'PROVED mutants: %d' % cnt.get('PROVED', 0), flush=True)
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('small')
    s.add_argument('--count', type=int, default=150)
    s.add_argument('--seed', type=int, default=1)
    s.add_argument('--curves', default='4,5,6')
    s.add_argument('--levels', type=int, default=1)
    s.add_argument('-o', '--out')
    g = sub.add_parser('big')
    g.add_argument('qasm')
    g.add_argument('hints')
    g.add_argument('--curve', required=True)
    g.add_argument('--only')
    g.add_argument('--seed', type=int, default=5)
    g.add_argument('-o', '--out')
    a = ap.parse_args(argv)
    if a.cmd == 'small':
        st = run_small(a.count, a.seed, [int(x) for x in a.curves.split(',')], a.levels, a.out)
        print(json.dumps({k: v for k, v in st.items() if k not in ('by_label', 'by_kind')}))
        print('by label [wrong, wrong PROVED, correct, correct PROVED]:', json.dumps(st['by_label']))
        print('by kind:', json.dumps(st['by_kind']))
        return 1 if st['wrong_proved'] else 0
    res = run_big(a.qasm, a.hints, a.curve, a.only.split(',') if a.only else None, a.out, a.seed)
    return 1 if any(r['status'] == 'PROVED' for r in res) else 0


if __name__ == '__main__':
    sys.exit(main())
