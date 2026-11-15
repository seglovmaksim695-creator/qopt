"""bandc worker: search + verification in a separate process (run by ``candidates``; killed on timeout).

    python -m qoptct.contrib.bandc.worker CIRCUIT.qc OUT.jsonl [--clean] [--budget SEC] [--mem-mb MB]

Every VERIFIED candidate is appended to OUT.jsonl as soon as it is found (one JSON object per line), so a partial
result survives a kill. A watchdog thread terminates the process when the time or memory budget is exceeded.
"""
import json
import os
import random
import re
import resource
import signal
import sys
import tempfile
import threading
import time

from . import compiler as M
from . import lcheck, vpt

MAX_QUBITS = 1600          # data wires
MAX_GATES = 40000          # logical gates
SV_MAX = 16                # total qubits for the statevector check (qoptct.verify_qasm)


# ---------------------------------------------------------------- search space
def _supp(g):
    nm, q = g
    if nm in ('cx', 'ccx'): return set(q[:-1]), {q[-1]}
    if nm == 'x': return set(), set(q)
    return set(q), set()          # phase gates "read" all their wires


def _commute(g, h):
    rg, wg = _supp(g); rh, wh = _supp(h)
    if g[0] == 'x' and h[0] == 'x': return True
    if wg & (rh | wh) - (wg & wh) or wh & (rg | wg) - (wg & wh): return False
    if wg & wh and not (g[0] in ('cx', 'ccx', 'x') and h[0] in ('cx', 'ccx', 'x')): return False
    return not (wg & rh) and not (wh & rg)


def shuffle(lg, rng, steps=None):
    """random permutation of adjacent commuting logical gates (same circuit)"""
    lg = list(lg); steps = steps or 4 * len(lg)
    for _ in range(steps):
        i = rng.randrange(len(lg) - 1)
        if _commute(lg[i], lg[i + 1]): lg[i], lg[i + 1] = lg[i + 1], lg[i]
    return lg


T_, F_ = True, False
# (K, lazy, max_anc, lctl, qmax, lin_force); comments: lab (agent b1) Pareto points found with the setting
_CLEAN = [
    (None, T_, 16, T_, 48, T_),     # cycle_17_3 612/754
    (None, 0.4, 4, T_, 48, F_),     # mod_adder_1024 276/355
    (512, T_, 4, T_, 48, F_),       # mod_adder_1024 288/380
    (None, 0.4, 4, T_, 48, T_),     # mod_adder_1024 308/388
    (512, T_, 4, T_, 48, T_),       # qcla_adder_10 96/161
    (None, T_, 2, F_, 48, T_),      # csum_mux_9 56/100
    (None, T_, 64, T_, 48, T_),     # ham15-high 720/1142
    (None, T_, 16, T_, 48, F_),     # ham15-high 720/1308
    (None, F_, 16, F_, 48, T_),     # adder_8 108/196
    (None, 0.4, 16, F_, 48, T_),    # adder_8 112/220, ham15-high 816/1049
    (None, F_, 4, F_, 48, T_),      # ham15-high 848/983, qcla_mod_7 172/225, qcla_adder_10 116/149
    (None, T_, 4, F_, 48, T_),      # qcla_mod_7 144/269, ham15-med 164/387
    (None, 0.4, 4, F_, 48, T_),     # qcla_mod_7 164/238
    (None, 0.4, 4, T_, 24, T_),     # cycle_17_3 631/727 (variant 1)
    (None, 0.4, 16, T_, 24, T_),    # cycle_17_3 628/726 (variant 1)
    (None, 0.4, 0, T_, 24, T_),     # cycle_17_3 758/777
    (None, 0.4, 16, T_, 48, T_),    # ham15-med 136/346 (variant 2)
    (None, T_, 4, T_, 48, T_),      # mod_adder_1024 320/413, ham15-med 140/400
    (64, F_, 16, F_, 48, T_), (64, T_, 16, F_, 48, T_), (64, F_, 4, F_, 48, T_), (64, T_, 16, T_, 48, T_),
    (None, F_, 0, F_, 48, T_), (None, T_, 0, F_, 48, T_), (512, T_, 0, T_, 48, T_), (None, T_, 16, F_, 48, T_),
    (None, T_, 16, F_, 24, T_),
]
_STRICT = [
    (64, F_, 16, F_, 48, T_),       # gf2^n_mult, ham15-med 148/240, adder_8 116/213, qcla_mod_7 184/260
    (64, 0.4, 16, F_, 48, T_),      # mod_adder_1024 608/860, cycle_17_3 1789/2359
    (None, T_, 16, T_, 48, T_),     # csum_mux_9 56/112 (no cut: can be slow in the strict model)
    (512, F_, 16, F_, 48, T_),      # ham15-high 752/1014
    (64, T_, 16, T_, 48, T_),       # qcla_adder_10 100/194
    (None, T_, 2, F_, 48, T_),      # csum_mux_9 56/112 (first lab run)
    (64, T_, 16, F_, 48, T_),       # ham15-med 180/425
    (512, T_, 16, T_, 24, T_),      # qcla_mod_7 192/392
    (512, 0.4, 16, F_, 24, T_), (64, F_, 4, F_, 48, T_), (512, T_, 4, F_, 48, T_), (512, T_, 0, F_, 48, T_),
    (None, F_, 16, F_, 48, T_), (None, 0.4, 16, F_, 48, T_), (64, T_, 64, T_, 24, T_), (64, F_, 0, F_, 48, T_),
]
_VARIANTS = [(None, 0.4, 16, T_, 24, T_), (None, 0.4, 4, T_, 24, T_), (None, 0.4, 16, T_, 48, T_),
             (None, T_, 4, T_, 48, T_), (None, 0.4, 4, T_, 48, T_), (None, T_, 16, F_, 48, T_),
             (None, T_, 16, F_, 24, T_), (None, T_, 4, F_, 48, T_), (512, 0.4, 4, F_, 48, T_)]


def configs(n, G, clean):
    """(variant, K, lazy, max_anc, lctl, qmax, lin_force, strict_compile) in priority order (lab agent b1 Pareto
    settings first). variant > 0: random permutation of commuting logical gates; lin_force=False: compiler version
    b1/v7; strict_compile: in the clean model compile without the knowledge of clean wires."""
    half = n // 2 + 8       # max_anc = -1 in the tables: about one ancilla per Toffoli of a ripple adder
    if n > 100:   # large: cut K=64 keeps the polynomials small (lab scaling runs: adders up to 1536 wires)
        big = [(64, T_, 16, F_, 48, T_), (64, F_, -1, F_, 48, T_)]
        if not clean: big.reverse()
        C = big + [(512, T_, 16, T_, 48, T_), (64, 0.4, 16, F_, 48, T_)]
    else:
        C = list(_CLEAN if clean else _STRICT)
        C.insert(1 if not clean else len(C), (64, F_, -1, F_, 48, T_))
    out = [(0,) + c + (F_,) for c in C]
    if n <= 100:
        if clean:   # strict compilations are valid in the clean model too (ham15-high: 752/1014 came from strict)
            out += [(0,) + c + (T_,) for c in _STRICT[:5]]
        for v in (1, 2):
            out += [(v,) + c + (F_,) for c in (_VARIANTS if clean else C[:6])]
    return [c[:3] + ((half if c[3] == -1 else c[3]),) + c[4:] for c in out]


# ---------------------------------------------------------------- verification
def strip_resets(txt):
    """Remove ``reset q[a];`` lines emitted after the X-measurement of a data wire. Sound only because the wire is
    provably |0> there: the line before is ``if (c[k]) {... x q[a]; ...}`` and before it ``c[k] = measure q[a];``."""
    L = txt.split('\n'); out = []
    for i, line in enumerate(L):
        m = re.fullmatch(r'reset q\[(\d+)\];', line.strip())
        if not m: out.append(line); continue
        a = m.group(1)
        mi = re.fullmatch(r'if \(c\[(\d+)\]\) \{(.*)\}', L[i - 1].strip()) if i >= 2 else None
        mm = re.fullmatch(r'c\[(\d+)\] = measure q\[(\d+)\];', L[i - 2].strip()) if i >= 2 else None
        if not (mi and mm and mm.group(1) == mi.group(1) and mm.group(2) == a
                and re.search(r'(^|;)\s*x q\[%s\];' % a, mi.group(2))):
            raise ValueError('reset not preceded by measure + x correction')
    return '\n'.join(out)


def mutants(txt, seed=7, ncx=4):
    """qasm-level canaries: replace one t by tdg; delete one cx (up to ncx different ones)."""
    L = txt.split('\n'); rng = random.Random(seed); out = []
    ts = [i for i, l in enumerate(L) if l.startswith('t ')]
    if ts:
        i = rng.choice(ts); Mx = list(L); Mx[i] = 'tdg' + Mx[i][1:]; out.append(('t_to_tdg', '\n'.join(Mx)))
    cxs = [i for i, l in enumerate(L) if l.startswith('cx ')]
    for i in rng.sample(cxs, min(ncx, len(cxs))):
        Mx = list(L); del Mx[i]; out.append(('drop_cx', '\n'.join(Mx)))
    return out


def _ops_canaries(ops, rng, k=4):
    """op-level canaries: drop one CZ of a measurement correction; drop one cx (up to k random ones of each kind)."""
    out = []
    mzc = [i for i, (g, q) in enumerate(ops) if g == 'MZ' and any(gg == 'cz' for gg, _ in q[1])]
    for i in rng.sample(mzc, min(k, len(mzc))):
        g, q = ops[i]; c = list(q[1]); j = rng.choice([u for u, x in enumerate(c) if x[0] == 'cz']); del c[j]
        m = list(ops); m[i] = (g, (q[0], tuple(c))); out.append(('drop_corr_cz', m))
    cxs = [i for i, (g, q) in enumerate(ops) if g == 'cx']
    for i in rng.sample(cxs, min(k, len(cxs))):
        m = list(ops); del m[i]; out.append(('drop_cx', m))
    return out


def check_ops(ops, N, qc, clean, deadline):
    """exact logical-level check + op canaries (for each kind at least one mutant must be rejected; a single mutant
    can be neutral, e.g. a CZ on a wire that is known to be 0 in the clean model). -> (ok, note)"""
    try:
        r = lcheck.check_ops(ops, N, qc, clean=clean, deadline=deadline)
        if not r['ok']: return False, 'lcheck: ' + r.get('why', '')
        res = {}
        for nm, m in _ops_canaries(ops, random.Random(3)):
            if res.get(nm) is False: continue          # this kind is already caught
            res[nm] = lcheck.check_ops(m, N, qc, clean=clean, deadline=deadline)['ok']
        bad = [nm for nm, v in res.items() if v]
        if bad: return False, 'lcheck canary not caught: %s' % ','.join(bad)
        if not res: return False, 'no op canary available'
    except (lcheck.TooBig, ValueError) as e:
        return False, 'lcheck inconclusive: %s' % e
    return True, 'lcheck(exact; canaries %s caught)' % ','.join(sorted(res))


def _sv(txt, qc, clean, trials):
    from qoptct import verify_qasm
    fd, tp = tempfile.mkstemp(suffix='.qasm')
    try:
        os.write(fd, txt.encode()); os.close(fd)
        return float(verify_qasm.check(tp, qc, trials=trials, clean=clean))
    finally:
        try: os.unlink(tp)
        except OSError: pass


def check_qasm(txt_r, txt_out, N, qc, clean, deadline, big):
    """sampled check of the emitted QASM + qasm canaries. -> (ok, note)"""
    muts = mutants(txt_out if N <= SV_MAX else txt_r)
    if not muts: return False, 'no canary available'
    try:
        if N <= SV_MAX:
            F = _sv(txt_out, qc, clean, 8)
            if F < 1 - 1e-6: return False, 'statevector F=%.6f' % F
            res = {}
            for nm, mt in muts:
                try: res.setdefault(nm, []).append(_sv(mt, qc, clean, 8) >= 1 - 1e-3)
                except Exception:  # noqa: BLE001  (a mutant that does not even run counts as rejected)
                    res.setdefault(nm, []).append(False)
            meth = 'statevector(verify_qasm, 8 random states)'
        else:
            xt, mt_ = (3, 2) if big else (4, 2)
            r = vpt.check_text(txt_r, qc, x_trials=xt, m_trials=mt_, clean=clean, deadline=deadline)
            if not r['ok']: return False, 'vpt: %s' % (r.get('error') or 'F=%.6f' % r['F'])
            res = {}
            for nm, mt in muts:
                # same samples as the main check: the canary shows that THIS check has power
                try: res.setdefault(nm, []).append(vpt.check_text(mt, qc, x_trials=xt, m_trials=mt_, clean=clean,
                                                                  deadline=deadline)['ok'])
                except vpt.Inconclusive: raise
                except Exception:  # noqa: BLE001
                    res.setdefault(nm, []).append(False)
            meth = 'vpt(path-sum at classical output, %dx%d samples)' % (xt, mt_)
    except vpt.Inconclusive as e:
        return _prove_fallback(txt_out, qc, clean, deadline, str(e))
    # canary policy: t->tdg must be caught; at least one of the dropped-cx mutants must be caught
    if 't_to_tdg' in res and any(res['t_to_tdg']): return False, 'canary t_to_tdg not caught'
    if 'drop_cx' in res and all(res['drop_cx']): return False, 'canary drop_cx not caught'
    return True, meth + ', canaries caught: ' + ','.join(sorted(res))


def _prove_fallback(txt, qc, clean, deadline, why):
    """vpt inconclusive (too wide): exact proof by qoptct.verify.prove if available; canary t->tdg must not prove."""
    try:
        from qoptct.verify import prove
    except Exception:  # noqa: BLE001
        return False, 'qasm check inconclusive: %s' % why
    b = min(deadline - time.time() - 1, 120.0)
    if b < 5: return False, 'qasm check inconclusive: %s (no time for exact proof)' % why
    try:
        r = prove(txt, qc, clean=clean, budget_sec=b)
        if r.get('status') != 'PROVED': return False, 'vpt inconclusive (%s); prove: %s' % (why, r.get('status'))
        muts = [t for nm, t in mutants(txt) if nm == 't_to_tdg'] or [t for nm, t in mutants(txt) if nm == 'drop_cx'][:1]
        if not muts: return False, 'no canary available'
        b = min(deadline - time.time() - 1, 120.0)
        if b < 2 or prove(muts[0], qc, clean=clean, budget_sec=b).get('status') == 'PROVED':
            return False, 'exact-proof canary not rejected'
        return True, 'qoptct.verify.prove PROVED (%s), canary rejected' % r.get('method', '')
    except Exception as e:  # noqa: BLE001
        return False, 'exact proof failed: %s' % e


# ---------------------------------------------------------------- main
class _Timeout(BaseException):
    pass


def _alarm(*_):
    raise _Timeout()


def _watchdog(deadline, mem_mb):
    scale = 1 if sys.platform == 'darwin' else 1024      # ru_maxrss: bytes on macOS, KiB on Linux
    while True:
        if time.time() > deadline + 3: os._exit(124)
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale > mem_mb * 2 ** 20: os._exit(125)
        time.sleep(0.25)


def run(qc, out, clean, budget, mem_mb=3000, log=None):
    t0 = time.time(); deadline = t0 + budget
    threading.Thread(target=_watchdog, args=(deadline, mem_mb), daemon=True).start()
    names, hl = M.load_qc(qc); n = len(names)
    if n > MAX_QUBITS or len(hl) > MAX_GATES: return
    lg0 = M.logical(hl, n)                                   # raises Fail for non-classical circuits
    lcheck.load_qc(qc)                                        # independent parser must accept it too
    if not any(g in ('ccx', 'ccz') for g, _ in lg0): return   # nothing non-linear to gain
    big = n > 100 or len(lg0) > 1500
    ins = M.inputs_idx(qc, names) if clean else None
    signal.signal(signal.SIGALRM, _alarm)
    variants = {0: lg0}; rng_v = random.Random(1)
    done = {}            # key (T, CX, A) -> verified?
    verified = []
    model = 'clean' if clean else 'strict'
    seen = set()
    for cfg in configs(n, len(lg0), clean):
        left = deadline - time.time()
        if left < 2: break
        v, K, lazy, anc, lctl, qm, lf, sc = cfg
        if (lctl and lazy is False) or cfg in seen: continue
        seen.add(cfg)
        while v not in variants: variants[len(variants)] = shuffle(lg0, rng_v)
        M.reset_monos()
        # time cap per configuration: the first one may use most of the budget (large circuits), later ones less
        if n > 100 and not verified: cap = left * 0.6
        elif not verified: cap = min(left * 0.5, max(20.0, 0.25 * budget))
        else: cap = min(left * 0.5, max(15.0, 0.2 * budget))
        signal.setitimer(signal.ITIMER_REAL, max(1.0, cap))
        try:
            C = M.Compiler(n, variants[v], anc, random.Random(v), lazy=lazy, linmax=3, inputs=None if sc else ins, K=K,
                           lctl=lctl, qmax=qm, lin_force=lf)
            ops = C.run(); ops, N = M.compact(ops, n)
        except (_Timeout, Exception) as e:  # noqa: BLE001  (a failed configuration is skipped)
            signal.setitimer(signal.ITIMER_REAL, 0)
            if log: log('cfg %s: %s %s' % (cfg, type(e).__name__, str(e)[:60]))
            continue
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        texts = []
        for fold in (False, True):
            try:
                tr = M.to_qasm(ops, N, fold=fold, ndata=n); to = strip_resets(tr)
            except Exception as e:  # noqa: BLE001
                if log: log('emit failed: %s' % e)
                continue
            T, CX = M.count(to); key = (T, CX, N - n)
            if key in done or any(key == t[0] for t in texts): continue
            if any(o[0] <= key[0] and o[1] <= key[1] and o[2] <= key[2] for o in verified): continue
            texts.append((key, fold, tr, to))
        if log: log('cfg %s -> %s (%.1fs)' % (cfg, [k for k, *_ in texts], time.time() - t0))
        if not texts: continue
        ok, note1 = check_ops(ops, N, qc, clean, deadline)
        if log and not ok: log('  ops check: %s' % note1)
        for key, fold, tr, to in texts:
            done[key] = False
            if not ok: continue
            if time.time() > deadline - 1: break
            ok2, note2 = check_qasm(tr, to, N, qc, clean, deadline, big)
            if log: log('  verify %s: %s %s' % (key, ok2, note2))
            if not ok2: continue
            done[key] = True; verified.append(key)
            note = ('bandc %s: T=%d CX=%d extra=%d; K=%s lazy=%s anc=%d lctl=%d qmax=%d linforce=%d var=%d fold=%d%s; '
                    'checked: %s + %s' % (model, key[0], key[1], key[2], K, lazy, anc, lctl, qm, lf, v, fold,
                                          ' strict-compiled' if sc else '', note1, note2))
            with open(out, 'a') as f:
                f.write(json.dumps(dict(qasm=to, note=note, T=key[0], CX=key[1], extra=key[2])) + '\n')


def main(argv):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('qc'); ap.add_argument('out')
    ap.add_argument('--clean', action='store_true'); ap.add_argument('--budget', type=float, default=60.0)
    ap.add_argument('--mem-mb', type=int, default=3000); ap.add_argument('-v', action='store_true')
    a = ap.parse_args(argv)
    lg = (lambda s: print(s, file=sys.stderr, flush=True)) if a.v else None
    try:
        run(a.qc, a.out, a.clean, a.budget, a.mem_mb, log=lg)
    except M.Fail:
        pass
    except Exception as e:  # noqa: BLE001
        if lg: lg('worker error: %s: %s' % (type(e).__name__, e))
    sys.stdout.flush()
    os._exit(0)


if __name__ == '__main__':
    main(sys.argv[1:])
