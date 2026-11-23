"""Proof of Qualtran's Montgomery modular multiplication (``DirtyOutOfPlaceMontgomeryModMul``).

``prove_montgomery(n, w, p, patched=True)`` lowers the multiplier gate by gate with the strict bridge, cuts it at the
boundaries of the Qualtran nodes and proves every piece with the unchanged qoptct prover.  The verdict is a dictionary
(``verdict['verdict']`` is ``'PROVED'`` or ``'NOT_PROVED'``) with the segments, the sub-nodes of the generic window,
the wire glue, the hashes of the programs and the Qualtran version (installed, or with the patches PR 1997 + exact
integer tables).  Anything that is not PROVED (UNKNOWN, a skipped check, a failed structural check) gives NOT_PROVED.

What is proved (assumptions: ``verdict['assumptions']``): for odd ``p`` with ``2p <= 2**n``, every ``y < p`` and every
n-bit ``x``, the program is the exact Montgomery recurrence with window ``w``: outputs ``x``, ``y``, ``qrom_indices``
(the digits ``m_i``), ``reduced`` and ``target = T_N - p * reduced`` with ``reduced = [T_N >= p]``; all other wires
are returned to zero and no measurement outcome leaves a phase.  By the telescoping identity checked here,
``2**n * T_N = x*y + p * sum(c_i 2**(w i))`` and the invariant ``T_i < y + p``, this is ``x * y * 2**-n mod p`` with
the result below ``p``.

The ORIGINAL Qualtran (``patched=False``) adds ``p`` instead of subtracting it in the final reduction; it is NOT_PROVED
against this contract for small ``n`` (witness: any input with ``T_N >= p``, about a quarter of all inputs).
"""
import json
import os
import subprocess
import sys
import time
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

from .._cache import cache_dir

__all__ = ['prove_montgomery', 'prime_for', 'exhaustive_small', 'SCHEMA']

SCHEMA = 'qoptct.montgomery/1'
_HERE = Path(__file__).resolve().parent
_PKG = _HERE.parent
_ROOT = _PKG.parents[1]          # the directory that contains the ``qoptct`` package


def prime_for(n):
    """p used by the proofs: ``2**255 - 19`` for n = 256, else the largest prime <= 2**(n-1) (then 2p-1 < 2**n)"""
    from ._core import prime_for as _pf
    return _pf(n)


def _source_digest():
    h = sha256()
    files = sorted(_HERE.glob('*.py')) + sorted((_PKG / 'bridge').glob('*.py')) + [_PKG / 'patches' / 'manifest.json']
    for f in files:
        h.update(f.name.encode()); h.update(f.read_bytes())
    return h.hexdigest()


def _versions():
    from importlib.metadata import version, PackageNotFoundError

    def v(name):
        try:
            return version(name)
        except PackageNotFoundError:
            return None
    z3v = v('z3-solver')
    if z3v is None:
        try:
            import z3
            z3v = z3.get_version_string()
        except Exception:
            z3v = None
    return dict(python=sys.version.split()[0], qualtran=v('qualtran'), z3=z3v, numpy=v('numpy'))


def _env(patched):
    from ..patches import patched_env
    env = patched_env('c') if patched else dict(os.environ)
    old = env.get('PYTHONPATH')
    env['PYTHONPATH'] = str(_ROOT) + (os.pathsep + old if old else '')
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env.setdefault('BRIDGE_SKIP_QASM', '1')
    env.setdefault('BRIDGE_DOMAIN_MAX_NODES', '12000000')   # resource knob of the exact domain guard (default 700000)
    env['QOPT_CACHE'] = str(cache_dir())
    return env


def _run_stage(args, env, logpath):
    cmd = [sys.executable, '-B', '-m', 'qoptct.foreign.qualtran.montgomery._worker'] + [str(a) for a in args]
    t0 = time.time()
    with open(logpath, 'ab') as lg:
        lg.write(('\n=== %s\n' % ' '.join(cmd[4:])).encode())
        lg.flush()
        r = subprocess.run(cmd, env=env, stdout=lg, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
    tail = ''
    if r.returncode:
        tail = Path(logpath).read_text(errors='replace')[-600:]
    return r.returncode, tail, round(time.time() - t0, 1)


def telescoping(n, w):
    """exact symbolic check of 2^n T_N = x y + p sum c_i 2^(w i) for T_{i+1} = (T_i + x_i y + c_i p) / 2^w, T_0 = 0"""
    from ._core import telescoping as _t
    return _t(n, w)


def prove_montgomery(n, w, p=None, *, patched=True, mutant=None, budget=3600, max_mb=2800, seed=7, out=None, fresh=False,
                     work_dir=None, e2e_max_n=64):
    """Prove the Montgomery multiplier of Qualtran for bit size ``n``, window ``w`` and modulus ``p``.

    ``patched=True`` uses the copy of Qualtran with the patches the proofs were made with (``patches.patched_qualtran('c')``);
    the installed Qualtran is not modified.  ``patched=False`` uses the installed (original) Qualtran.
    ``mutant`` (a tag of ``montgomery.mutants``) runs the proof on a deliberately damaged program; it must not be PROVED.
    Heavy runs (n = 256) should be started under ``lab/tools/one`` (3 GB).  Work files go to ``$QOPT_CACHE/montgomery``.
    Returns the verdict dictionary; with ``out`` it is also written there as JSON.
    """
    if w < 1 or n % w:
        raise ValueError('window size w must divide the bit size n')
    p = p or prime_for(n)
    if p % 2 == 0 or p < 3 or 2 * p > (1 << n):
        raise ValueError('p must be odd with 2p <= 2**n')
    t_all = time.time()
    kind = 'c' if patched else ''
    env = _env(patched)
    src = _source_digest()
    key = 'n%d_w%d_p%s_%s_%s%s' % (n, w, ('%x' % p)[:16], kind or 'orig', src[:12], ('_mut' + mutant) if mutant else '')
    wd = Path(work_dir) if work_dir else cache_dir('montgomery', key)
    wd.mkdir(parents=True, exist_ok=True)
    log = wd / 'run.log'
    reasons = []
    stages = {}
    lower_p = wd / 'lower.pkl'
    inv_p = {'first': wd / 'inv_first.json', 'generic': wd / 'inv_generic.json'}
    prove_p = wd / 'prove.json'

    def cached_ok():
        sh = Path(str(lower_p) + '.sha256')
        return (not fresh) and lower_p.is_file() and sh.is_file() and sh.read_text().strip() == sha256(lower_p.read_bytes()).hexdigest()
    if not cached_ok():
        for f in [lower_p, *inv_p.values(), prove_p]:
            if f.exists(): f.unlink()
        rc, tail, secs = _run_stage(['lower', '--n', n, '--w', w, '--p', p, '--out', lower_p, '--kind', kind, '--mutant', mutant or '',
                                     '--e2e-max-n', e2e_max_n], env, log)
        stages['lower'] = dict(rc=rc, seconds=secs, cached=False)
        if rc:
            return _fail(n, w, p, patched, mutant, 'lowering/structural stage failed: ' + tail[-300:], stages, src, out, t_all)
    else:
        stages['lower'] = dict(rc=0, seconds=0, cached=True)
    import pickle
    D = pickle.load(open(lower_p, 'rb'))
    if (D['n'], D['w'], D['p'], D['kind'], D['mutant']) != (n, w, p, kind or None, mutant):
        return _fail(n, w, p, patched, mutant, 'cache does not belong to this request', stages, src, out, t_all)
    inv = {}
    for part in ('first', 'generic'):
        if fresh or not inv_p[part].is_file():
            rc, tail, secs = _run_stage(['inv', '--n', n, '--w', w, '--p', p, '--part', part, '--out', inv_p[part]], env, log)
            stages['inv_' + part] = dict(rc=rc, seconds=secs)
            if rc:
                return _fail(n, w, p, patched, mutant, 'invariant stage failed: ' + tail[-300:], stages, src, out, t_all)
        inv.update(json.load(open(inv_p[part])))
    rc, tail, secs = _run_stage(['prove', '--n', n, '--w', w, '--p', p, '--inp', lower_p, '--out', prove_p, '--budget', budget,
                                 '--max-mb', max_mb, '--seed', seed], env, log)
    stages['prove'] = dict(rc=rc, seconds=secs)
    if rc:
        return _fail(n, w, p, patched, mutant, 'proof stage failed: ' + tail[-300:], stages, src, out, t_all)
    R = json.load(open(prove_p))
    ck = D['checks']
    # ---- verdict: every piece must hold
    if not ck['glue']['ok']: reasons.append('wire glue: %s' % ck['glue']['notes'][:3])
    if not ck['bookkeeping']['ok']: reasons.append('live-wire bookkeeping: %s' % ck['bookkeeping']['notes'][:3])
    if not ck['windows_ok']: reasons.append('coverage/freshness of sub-node wires failed in windows %s' % [x['key'] for x in ck['windows'] if not x['ok']][:4]
                                           + ' first notes: %s' % next((x['notes'] for x in ck['windows'] if not x['ok']), []))
    if ck['tiles_bad']: reasons.append('children of a window do not tile its span: %s' % list(ck['tiles_bad'])[:3])
    if not ck['prefix']['ok']: reasons.append('X(reduced) before window 0: %s' % ck['prefix']['notes'])
    if not ck['same_hash']: reasons.append('windows 1..N-1 are not identical programs (QASM/layout hashes differ)')
    if not ck['same_structure']: reasons.append('windows 1..N-1 differ in the structure of their sub-nodes')
    inv_ok = bool(inv) and all(all(v.values()) for v in inv.values())
    if not inv_ok: reasons.append('interval invariants failed: %s' % inv)
    tel = telescoping(n, w)
    if not tel: reasons.append('telescoping identity failed')
    reasons += R['reasons']
    e2e = D.get('e2e')
    if e2e and e2e.get('mismatches'):
        reasons.append('end-to-end simulation differs from the Montgomery contract on %d of %d inputs' % (e2e['mismatches'], e2e['cases']))
    ver = _versions()
    if not ver['z3']:
        reasons.append('z3 is not available (the verdict is never PROVED without it)')
    if ver['qualtran'] != '0.7.0':
        reasons.append('Qualtran %s, the proofs are for 0.7.0' % ver['qualtran'])
    witness = None
    if e2e and e2e.get('first_witness'):
        witness = dict(kind='end-to-end', **e2e['first_witness'])
    for r in R['whole']:
        if r.get('witness') is not None and witness is None:
            witness = dict(kind='segment ' + r['label'], input=r['witness'], status=r['status'])
    verdict = dict(
        schema=SCHEMA, verdict='PROVED' if not reasons else 'NOT_PROVED', reasons=reasons, n=n, w=w, p=p, p_hex=hex(p),
        patched=bool(patched), mutant=mutant,
        contract='x,y,qrom_indices(digits m_i),reduced=[T_N>=p],target=T_N-p*reduced (PR1997 contract); other wires zero; no phase',
        qualtran=D['qualtran'], versions=ver,
        program=dict(nq=D['meta']['nq'], ops=D['meta']['ops'], ops_sha256=D['meta']['ops_sha256'], windows=D['nw'],
                     generic_window_qasm_sha256=(ck['generic_window_hash'] or [None, None])[0],
                     generic_window_layout_sha256=(ck['generic_window_hash'] or [None, None])[1]),
        glue=ck['glue'], bookkeeping=ck['bookkeeping'], prefix_X_reduced=ck['prefix'],
        window_checks=dict(all_windows_ok=ck['windows_ok'], windows=ck['n_windows'], same_hash=ck['same_hash'], same_structure=ck['same_structure'],
                           tiles_bad=ck['tiles_bad']),
        invariants=inv, telescoping=tel, end_to_end=e2e,
        segments=R['whole'], window0_by_subnodes=R['win0_split'], generic_window=R['generic'], witness=witness,
        assumptions=[
            'Qualtran %s with the patches %s' % (ver['qualtran'], 'PR1997 + exact integer tables (kind c)' if patched else 'NONE (installed original)'),
            'p odd, 2p <= 2^n, y < p, x any n-bit number',
            'the rules of the bridge (A1 controlled-And, C1 QROAMCleanAdjoint) and qoptct kprove/BDD with the exact partial-domain guard',
            'windows 1..N-1 are identical programs (equal QASM and layout hash, checked) with their own wires (glue and bookkeeping checked)',
            'invariant T_i < y + p of the generic window: BDD for first window / small n, interval arithmetic where BDD does not fit in memory',
            'w >= 4: composition of the sub-nodes by wire positions + interval arithmetic + exact integer runs (not BDD)',
        ],
        stages=stages, source_sha256=src, lower_cache_sha256=Path(str(lower_p) + '.sha256').read_text().strip(), work_dir=('$QOPT_CACHE/montgomery/' + wd.name) if not work_dir else wd.name,
        seconds=round(time.time() - t_all, 1))
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text(json.dumps(verdict, indent=1, default=str))
    return verdict


def _fail(n, w, p, patched, mutant, why, stages, src, out, t0):
    v = dict(schema=SCHEMA, verdict='NOT_PROVED', reasons=[why], n=n, w=w, p=p, patched=bool(patched), mutant=mutant, stages=stages,
             source_sha256=src, versions=_versions(), seconds=round(time.time() - t0, 1))
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text(json.dumps(v, indent=1, default=str))
    return v


def exhaustive_small(n, w, p, *, patched=True):
    """All inputs (x < 2^n for n <= 8, else x < p; y < p) through the lowered program against the exact recurrence.

    Runs in a fresh process with the right Qualtran.  Returns dict(total, mismatch, bad_phase, bad_rest, first).
    ``patched=False`` compares the installed Qualtran with the FIXED contract: the mismatches are exactly the inputs with T_N >= p.
    """
    env = _env(patched)
    code = ('import json,sys\n'
            'from qoptct.foreign.qualtran.montgomery import _core as C\n'
            'n,w,p=%d,%d,%d\n'
            'low,tr=C.lower_translate(C.build(n,w,p))\n'
            'tot=bad=bph=brest=bad_red=0;first=None\n'
            'for x in range(1<<n if n<=8 else p):\n'
            '  for y in range(p):\n'
            '    res,ph,rest=C.run_basis(low,tr,n,x,y); e=C.reference(n,w,p,x,y); tot+=1\n'
            '    if any(ph): bph+=1\n'
            '    if rest: brest+=1\n'
            '    if any(res[k]!=e[k] for k in ("x","y","target","qrom_indices","reduced")):\n'
            '      bad+=1; bad_red+=int(e["T"]>=p)\n'
            '      if first is None: first=dict(x=x,y=y,got=res["target"],expected=e["target"],T_N_ge_p=bool(e["T"]>=p))\n'
            'print("RESULT"+json.dumps(dict(total=tot,mismatch=bad,mismatch_with_T_N_ge_p=bad_red,bad_phase=bph,bad_rest=brest,first=first)))\n') % (n, w, p)
    r = subprocess.run([sys.executable, '-B', '-c', code], env=env, capture_output=True, text=True, stdin=subprocess.DEVNULL)
    for line in r.stdout.splitlines():
        if line.startswith('RESULT'):
            return json.loads(line[6:])
    raise RuntimeError('exhaustive_small failed: ' + (r.stderr or r.stdout)[-400:])
