"""Prove one ECAdd step with the ew method.
Example: python -m qoptct.foreign.qualtran.ecadd.run_step --step 2 --n 16 --p 32749 --out step2.json"""
import argparse, hashlib, json, os, sys, time, resource


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def code_hashes():
    """SHA-256 (16 hex) of the code that decides the verdict: prover, specs, bridge translation."""
    here = os.path.dirname(os.path.abspath(__file__))
    br = os.path.join(os.path.dirname(here), 'bridge')
    files = {'ecadd/core.py': os.path.join(here, 'core.py'), 'ecadd/specs.py': os.path.join(here, 'specs.py'),
             'ecadd/bv.py': os.path.join(here, 'bv.py'), 'ecadd/translate_ew.py': os.path.join(here, 'translate_ew.py'),
             'ecadd/opaque.py': os.path.join(here, 'opaque.py'),
             'bridge/qualtran_layout.py': os.path.join(br, 'qualtran_layout.py'),
             'bridge/emit.py': os.path.join(br, 'emit.py'), 'bridge/rules.py': os.path.join(br, 'rules.py')}
    return {k: sha(v)[:16] for k, v in files.items()}


def version_info():
    """Qualtran that is imported: kind ('original'/'e'/'e4'/'c'/'other'), version label, hashes of the touched files."""
    from qoptct.foreign.qualtran._variant import detect_variant
    v = detect_variant()
    v['version'] = 'patched' if v['kind'] in ('e', 'e4') else ('original' if v['kind'] == 'original' else v['kind'])
    return v


def step_bloq(k, n, p):
    from qualtran.bloqs.cryptography.ecc import ec_add as EC
    cls = [EC._ECAddStepOne, EC._ECAddStepTwo, EC._ECAddStepThree, EC._ECAddStepFour, EC._ECAddStepFive,
           EC._ECAddStepSix][k - 1]
    return cls(n=n, mod=p)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--step', type=int, required=True)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', type=str, required=True)
    ap.add_argument('--out', default=None)
    ap.add_argument('--nolift', action='store_true')
    ap.add_argument('--observe', action='store_true', help='steps 5,6: without P5/P6, final values of the released wires = the formulas')
    ap.add_argument('--only-flags', default=None,
                    help='restrict the flag sets, e.g. "f1=0,f2=1,f3=0,f4=0,ctrl=0;f1=1,..." (the claim holds only for them)')
    a = ap.parse_args(argv)
    p = eval(a.p, {}, {})  # allows 2**255-19
    from qoptct.foreign.qualtran.ecadd.core import build_program, Run, Fail
    from qoptct.foreign.qualtran.ecadd.specs import STEPS
    vi = version_info()
    spec = STEPS[a.step](a.n, p, patched=(vi['version'] == 'patched'), observe=a.observe)
    t0 = time.time()
    res = dict(step=a.step, n=a.n, p=str(p), observe=a.observe, code_sha256_16=code_hashes(), **vi, method='ew (ROBDD, qoptct.bdd) + opaque MM/KINV contracts',
               assumptions=['MM: x,y<p -> target<p; garbage arbitrary; MM_adj exact inverse on paired outputs',
                            'KINV: x<p -> x_out<p; junk arbitrary; KINV_adj exact inverse on paired outputs'])
    try:
        prog = build_program(step_bloq(a.step, a.n, p), use_lift=not a.nolift, keep_ops=False)
        res.update(build_seconds=round(time.time() - t0, 2), logical_ops=prog['nops'],
                   opaque=[it[1]['kind'] for it in prog['items'] if it[0] == 'opaque'])
        R = Run(prog, spec, a.n, p)
        only = None
        if a.only_flags:
            only = [{kv.split('=')[0]: int(kv.split('=')[1]) for kv in grp.split(',')} for grp in a.only_flags.split(';')]
            res['only_flags'] = only
        res['status'] = R.run(only_flags=only)
        res['stats'] = R.stats
    except Fail as e:
        res.update(status='REFUTED', detail=str(e), info={k: (v if k != 'counterexample' else
                   {kk: str(vv) for kk, vv in (v or {}).items()}) for k, v in e.detail.items()})
    except MemoryError as e:
        res.update(status='UNKNOWN', detail=f'memory: {e}')
    except Exception as e:
        import traceback
        res.update(status='ERROR', detail=f'{type(e).__name__}: {e}', tb=traceback.format_exc()[-2000:])
    res['wall_seconds'] = round(time.time() - t0, 2)
    res['peak_rss_mb'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 1)
    line = json.dumps(res, ensure_ascii=False)
    print(line[:3000])
    if a.out:
        with open(a.out, 'w') as f:
            f.write(json.dumps(res, ensure_ascii=False, indent=1))
    return res


if __name__ == '__main__':
    main()
