"""run -- the evidence run for one pair (n, p): composition, Walker proof, execution mode, independent simulation.

    python -m qoptct.foreign.qualtran.kaliski.run --n N --p P [--prove] [--run all|sample] [--sim K] [--out FILE]
The final status PROVED is set only by verdict.py (it also needs the Z3 lemmas, the primality of p and the canaries).
"""
import argparse, hashlib, json, os, random, resource, sys, time


def samples(p, n, k, seed=20261123):
    """k distinct inputs in [0, p): the rare ones first (0, 1, 2, 3, p-1, p-2, p-3, (p+-1)/2, p//3, 2^j, 2^j - 1,
    2^j + 1, and the inverses of 2^j), then fixed pseudo-random ones."""
    k = min(k, p)
    rnd = random.Random(seed)
    must = [0, 1, 2, 3, p - 1, p - 2, p - 3, (p - 1) // 2, (p + 1) // 2, p // 3]
    extra = []
    for j in range(p.bit_length()):
        extra += [(1 << j) % p, ((1 << j) - 1) % p, ((1 << j) + 1) % p, pow(2, -j, p)]
    rnd.shuffle(extra)
    xs = []
    for x in must + extra:
        if 0 <= x < p and x not in xs and len(xs) < max(k // 2, min(k, len(must))):
            xs.append(x)
    seen = set(xs)
    while len(xs) < k:
        x = rnd.randrange(p)
        if x not in seen:
            seen.add(x)
            xs.append(x)
    return sorted(xs)


def expected(n, p, x):
    out = (pow(x, -1, p) * pow(2, 2 * n, p)) % p if x else 0
    return out


def qualtran_classical(n, p, xs):
    from qualtran.bloqs.mod_arithmetic.mod_division import KaliskiModInverse
    b = KaliskiModInverse(n, p)
    res = {}
    for x in xs:
        try:
            r = b.call_classically(x=x)
            res[x] = (int(r[0]), int(r[1]))
        except Exception as e:                  # noqa: BLE001
            res[x] = ('ERROR', '%s: %s' % (type(e).__name__, str(e)[:200]))
    return res


def word(store_or_vals, wires, j):
    return sum(((store_or_vals[w] >> j) & 1) << i for i, w in enumerate(wires))


def manifest():
    """SHA-256 of every source file the result depends on: this client, the bridge, the qoptct provers, Qualtran."""
    import glob, qualtran
    pkg = os.path.dirname(os.path.abspath(__file__))                 # .../qoptct/foreign/qualtran/kaliski
    qoptct_root = os.path.dirname(os.path.dirname(os.path.dirname(pkg)))   # .../qoptct
    qroot = os.path.dirname(qualtran.__file__)
    used = ('lower.py', 'chain.py', 'pq.py', 'run.py', 'lemmas_z3.py', 'canary.py')
    groups = dict(src=(os.path.dirname(pkg), sorted(os.path.join(pkg, f) for f in used)),
                  bridge=(os.path.dirname(pkg), sorted(glob.glob(os.path.join(os.path.dirname(pkg), 'bridge', '*.py')))),
                  qoptct=(os.path.dirname(qoptct_root),
                          sorted(glob.glob(os.path.join(qoptct_root, 'verify', '**', '*.py'), recursive=True))),
                  qualtran=(qroot, sorted(glob.glob(os.path.join(qroot, '**', '*.py'), recursive=True))))
    out = {}
    for g, (root, files) in groups.items():
        h = {os.path.relpath(f, root): hashlib.sha256(open(f, 'rb').read()).hexdigest() for f in files}
        out[g] = dict(files=len(h), sha256=hashlib.sha256(json.dumps(h, sort_keys=True).encode()).hexdigest())
        if g in ('src', 'qoptct'):
            out[g]['each'] = h
    out['qualtran_version'] = qualtran_version()
    out['mod_division_sha256'] = hashlib.sha256(open(os.path.join(qroot, 'bloqs', 'mod_arithmetic', 'mod_division.py'), 'rb').read()).hexdigest()
    return out


def qualtran_version():
    from importlib.metadata import version
    return version('qualtran')


def run_pair(n, p, *, prove=True, run='sample', sim=127, out=None, verbose=True, run_inputs=256):
    """Compose KaliskiModInverse(n, p) from the real Qualtran pieces and run the Walker proof, the execution mode and the
    independent simulation.  Returns the evidence dictionary (the final verdict is made by ``verdict.judge``)."""
    a = argparse.Namespace(n=n, p=p, prove=prove, run=run, sim=sim, out=out or '')
    from .lower import Composed
    from .chain import QChain, ProofFailure, build_marks
    t0 = time.time()
    rep = dict(n=n, p=p, qualtran=qualtran_version(), started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
               manifest_before=manifest())
    C = Composed(n, p)
    rep['compose'] = dict(sec=round(time.time() - t0, 2), nq=C.nq, nops=C.nops, piece_ops=C.piece.L.n,
                          piece_qasm_sha256=C.piece.qasm_sha256, nmeas=C.nmeas,
                          top_ops=C.TL.n, top_ops_sha256=hashlib.sha256(C.TL.ops.tobytes()).hexdigest())
    if verbose:
        print('composed', rep['compose'], flush=True)
    marks = build_marks(C)
    log = (lambda m: print('%7.1fs %s' % (time.time() - t0, m), flush=True)) if verbose else (lambda m: None)
    if a.prove:
        t = time.time()
        ch = QChain(C, lifted=C.lifted, marks=marks, progress=log)
        try:
            ch.run()
            st = 'CHAIN_PROVED'
            err = None
        except ProofFailure as e:
            st, err = 'CHAIN_FAILED', dict(msg=str(e), detail={k: (v if not isinstance(v, dict) else {str(a): b for a, b in list(v.items())[:40]}) for k, v in e.detail.items()})
        rep['prove'] = dict(status=st, error=err, sec=round(time.time() - t, 2), segments=ch.nseg, maxnodes=ch.maxnodes,
                            nodes=ch.nodes, assumed=dict(ch.assumed), admitted=list(ch.admitted), pos=ch.pos, F=ch.F,
                            flow=ch.flow.summary())
        log('prove: %s %s' % (st, err['msg'] if err else ''))
    if a.run:
        try:
            t = time.time()
            xs = list(range(p)) if a.run == 'all' else samples(p, n, run_inputs)
            ch = QChain(C, xs=xs, marks=marks)
            ch.run()
            Xo, J = C.outputs['x'], C.outputs['junk']
            ref = qualtran_classical(n, p, xs)
            bad_out, bad_junk = [], []
            for j, x in enumerate(xs):
                o, jk = word(ch.store, Xo, j), word(ch.store, J, j)
                if o != expected(n, p, x):
                    bad_out.append(x)
                if ref[x] != (o, jk):
                    bad_junk.append(x)
            nviol = len(ch.viol)
            # the integer recursion P_Q of qkal_pq (= the SMT text of qkal_lemmas_z3) against the executed contracts,
            # after EVERY iteration, register by register; and the assumptions of P_Q evaluated in python
            from . import pq as qkal_pq
            bad_pq, bad_pq_checks = [], []
            dec = lambda masks, j: sum(((mk >> j) & 1) << i for i, mk in enumerate(masks))
            for j, x in enumerate(xs):
                out_pq, junk_pq, states, badc = qkal_pq.run(n, p, x)
                if badc:
                    bad_pq_checks.append((x, badc[:3]))
                ms, ts = [], []
                for k, tr in enumerate(ch.trace):
                    got = (dec(tr['u'], j), dec(tr['v'], j), dec(tr['r'], j), dec(tr['s'], j), bool(dec(tr['f'], j)))
                    ms.append(dec(tr['m'], j))
                    ts.append(dec(tr['t'], j))
                    if got != states[k + 1]:
                        bad_pq.append((x, k))
                        break
                jk = 0
                for b in ms + ts:
                    jk = 2 * jk + b
                if (out_pq, junk_pq) != (word(ch.store, Xo, j), jk) or jk != word(ch.store, J, j):
                    bad_pq.append((x, 'out/junk'))
            nviol = len(ch.viol)
            rep['run'] = dict(inputs=len(xs), violations=nviol, first_violations=ch.viol[:10], wrong_output=bad_out[:20],
                              n_wrong_output=len(bad_out), differs_from_qualtran_classical=bad_junk[:20],
                              n_differs=len(bad_junk), pq_vs_contracts_mismatch=bad_pq[:20], n_pq_mismatch=len(bad_pq),
                              pq_assumption_failures=bad_pq_checks[:10], iterations_traced=len(ch.trace),
                              xs_sha256=hashlib.sha256(json.dumps(xs).encode()).hexdigest(), sec=round(time.time() - t, 2),
                              flow=ch.flow.summary())
            log('run: %s' % {k: v for k, v in rep['run'].items() if k != 'first_violations'})
        except ProofFailure as e:  # wire flow or another client error in mode run
            rep['run'] = dict(error=str(e))
    if a.sim:
        from qoptct.verify.compose.walker import simulate_flat
        t = time.time()
        xs = samples(p, n, a.sim)
        X = C.inputs['x']
        init = {}
        for i, w in enumerate(X):
            init[w] = sum(1 << j for j, x in enumerate(xs) if (x >> i) & 1)
        vals, err = simulate_flat(C.lifted, init, len(xs))
        Xo, J = C.outputs['x'], C.outputs['junk']
        keep = set(Xo) | set(J)
        dirty = [w for w in range(C.nq) if w not in keep and vals[w]]
        ref = qualtran_classical(n, p, xs)
        bad = [x for j, x in enumerate(xs) if (word(vals, Xo, j), word(vals, J, j)) != ref[x]
               or word(vals, Xo, j) != expected(n, p, x)]
        rep['sim'] = dict(inputs=len(xs), gadget_error=None if err is None else [err[0], err[1]], dirty_wires=len(dirty),
                          mismatches=bad[:20], n_mismatch=len(bad), sec=round(time.time() - t, 2),
                          xs_sha256=hashlib.sha256(json.dumps(xs).encode()).hexdigest())
        log('sim: %s' % rep['sim'])
    rep['sec'] = round(time.time() - t0, 2)
    rep['peak_rss_mb'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 1)
    after = manifest()
    rep['sources_unchanged_during_run'] = all(after[g]['sha256'] == rep['manifest_before'][g]['sha256']
                                              for g in ('src', 'bridge', 'qoptct', 'qualtran'))
    if a.out:
        with open(a.out, 'w') as f:
            json.dump(rep, f, indent=1, default=str)
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description='Evidence run for one pair (n, p); the verdict is made by verdict.py')
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', type=int, required=True)
    ap.add_argument('--prove', action='store_true')
    ap.add_argument('--run', default='')
    ap.add_argument('--sim', type=int, default=0)
    ap.add_argument('--out', default='')
    a = ap.parse_args(argv)
    rep = run_pair(a.n, a.p, prove=a.prove, run=a.run, sim=a.sim, out=a.out)
    print(json.dumps({k: rep[k] for k in rep if k in ('prove', 'sec', 'peak_rss_mb')}, default=str)[:3000])


if __name__ == '__main__':
    main()
