# SPDX-License-Identifier: Apache-2.0
"""Pluggable proposers. A proposer enumerates parameter dicts for a target (`params(target)`, deterministic order) and
builds a candidate circuit text from one dict (`build(target, params, workdir)` -> qasm text or None).

Proposers shipped:
  lookup    parameter search for table lookups (xor / add): seeds and splits of the bilinear plan search, Gray-code
            fan-out batch sizes, QROAM with one select bit (qoptct.crypto.lookup builders; resynthesis of the lookup)
  optimize  the optimizer portfolio: `qoptct optimize` on the reference .qc with several modes (cx / balanced / t,
            measurement-based AND compiler, post-pass) under a time limit

A proposer's output is never trusted: the loop proves every candidate that would improve the front.
"""
import itertools
import os
import subprocess
import sys

from .costs import pad_qubits


def _wires(ref_text):
    for ln in ref_text.splitlines():
        if ln.startswith('.v'):
            return len(ln.split()) - 1
    return 0


# ------------------------------------------------------------------------------------------------------ lookup
def lookup_params(target):
    g = target.get('gen')
    if not g or g.get('family') != 'lookup':
        return
    c, kind = g['c'], g['kind']
    split_opts = [None] + [[h] for h in (1, 2, 3) if h <= c - h] + ([[1, 2]] if c >= 4 else [])
    # round 0: the cheap variants first
    if kind == 'xor':
        yield {'p': 'lookup', 'm': 'qroam1'}
        for b in range(2, 9):
            yield {'p': 'lookup', 'm': 'gray', 'b': b, 'l': None}
        for b, l in itertools.product(range(2, 9), range(1, c)):
            yield {'p': 'lookup', 'm': 'gray', 'b': b, 'l': l}
    else:
        for b in range(2, 9):
            yield {'p': 'lookup', 'm': 'add-gray', 'b': b, 'l': None}
    for seed, tries in itertools.product(range(0, 60), (None, 4)):
        for sp in split_opts:
            yield {'p': 'lookup', 'm': 'narrow' if kind == 'xor' else 'add-narrow', 'seed': seed, 'splits': sp,
                   'tries': tries}


def lookup_build(target, params, workdir):
    from ..crypto import lookup as L
    g = target['gen']
    T, c, w = g['table'], g['c'], g['w']
    m = params['m']
    if m == 'qroam1':
        lc = L.qroam1_lookup(T, c, w)
    elif m == 'gray':
        lc = L.davio_lookup(T, c, w, fanout='gray', gray_b=params['b'], gray_l=params['l'], gray_force=True)
    elif m == 'narrow':
        lc = L.narrow_lookup(T, c, w, seed=params['seed'], splits=params['splits'], tries=params['tries'])
    elif m == 'add-gray':
        lc = L.lookup_add(T, c, w, fanout='gray', gray_b=params['b'], gray_l=params['l'], gray_force=True)
    elif m == 'add-narrow':
        lc = L.lookup_add(T, c, w, seed=params['seed'], splits=params['splits'])
    else:
        raise ValueError(m)
    return pad_qubits(lc.qasm, _wires(target['ref_text']))


# ------------------------------------------------------------------------------------------------------ optimizer
def optimize_params(target):
    for mode, meas, post, secs in [('cx', False, False, 60), ('balanced', True, False, 90), ('t', True, False, 90),
                                   ('cx', False, True, 120), ('balanced', True, True, 150), ('t', True, True, 150)]:
        yield {'p': 'optimize', 'mode': mode, 'meas': meas, 'post': post, 'secs': secs}


def optimize_build(target, params, workdir):
    ref = os.path.join(workdir, 'ref.qc')
    out = os.path.join(workdir, 'opt_out.qasm')
    open(ref, 'w').write(target['ref_text'])
    if os.path.exists(out):
        os.remove(out)
    cmd = [sys.executable, '-c', 'import sys; from qoptct.qopt_core import main; sys.argv[0]="qoptct"; main()',
           'optimize', ref, '-o', out, '--mode', params['mode'], '--time-limit', str(params['secs']),
           '--budget', str(params['secs'])]
    if target['clean']:
        cmd.append('--clean')
    if params['meas']:
        cmd.append('--meas')
    if params['post']:
        cmd.append('--post')
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(
        [os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))] + [os.environ.get('PYTHONPATH', '')]))
    try:
        subprocess.run(cmd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=params['secs'] * 3 + 300,
                       cwd=workdir)
    except subprocess.TimeoutExpired:
        return None
    if not os.path.exists(out):
        return None
    return pad_qubits(open(out).read(), _wires(target['ref_text']))


PROPOSERS = {'lookup': (lookup_params, lookup_build), 'optimize': (optimize_params, optimize_build)}


def key_of(target_id, params):
    return target_id + '|' + '|'.join('%s=%s' % (k, params[k]) for k in sorted(params))


def stream(target, names):
    """Round-robin over the proposers (lookup candidates first, a few optimizer runs interleaved)."""
    gens = [(n, PROPOSERS[n][0](target)) for n in names if n in PROPOSERS]
    live = list(gens)
    while live:
        nxt = []
        for n, g in live:
            for _ in range(8 if n == 'lookup' else 1):
                try:
                    yield next(g)
                except StopIteration:
                    break
            else:
                nxt.append((n, g))
        live = nxt


def build(target, params, workdir):
    return PROPOSERS[params['p']][1](target, params, workdir)
