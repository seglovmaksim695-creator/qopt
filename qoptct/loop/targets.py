# SPDX-License-Identifier: Apache-2.0
"""Targets of the loop: a reference circuit (.qc, the specification), a model (strict or clean) and, when known, the
generator parameters that the proposers use.

Target specs (strings):
    crypto:<name>           a circuit of results/crypto (baseline = its .qasm; lookup_* get generator parameters)
    lookup:c<c>w<w>s<seed>  random table (xor lookup, z ^= T[a]); baseline = the Davio lookup
    lookupadd:c<c>w<w>s<seed>   random table, acc += T[a] mod 2^w; baseline = the plain narrow lookup-add
    qc:<path>[:clean]       any .qc file; no baseline, the first proved candidate becomes the baseline
"""
import json
import os
import random
import re

from .costs import pad_qubits
from .journal import sha256_bytes


class TargetError(Exception):
    pass


def repo_root():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(here, '..', '..'))


def _ref_wires(qc_text):
    for ln in qc_text.splitlines():
        if ln.startswith('.v'):
            return len(ln.split()) - 1
    return 0


def table_from_ref(qc_text, c, w):
    """z ^= T[a] reference (a0.., z0.., clean wires): read T by running the Toffoli network on every address."""
    names, ins, gates = [], [], []
    for ln in qc_text.splitlines():
        p = ln.split()
        if not p:
            continue
        if p[0] == '.v':
            names = p[1:]
        elif p[0] in ('tof', 'cnot', 'X', 'not', 'Tof'):
            gates.append(p[1:])
    idx = {n: i for i, n in enumerate(names)}
    T = []
    for a in range(1 << c):
        st = [0] * len(names)
        for i in range(c):
            st[idx['a%d' % i]] = (a >> i) & 1
        for g in gates:
            t = idx[g[-1]]
            if all(st[idx[x]] for x in g[:-1]):
                st[t] ^= 1
        T.append(sum(st[idx['z%d' % i]] << i for i in range(w)))
        if any(st[idx['a%d' % i]] != (a >> i) & 1 for i in range(c)):
            raise TargetError('reference changes the address wires')
    return T


def _random_table(c, w, seed):
    rng = random.Random(seed)
    return [rng.getrandbits(w) for _ in range(1 << c)]


def resolve(spec, root=None):
    """spec -> dict(id, spec, clean, ref_text, baseline_qasm (text or None), gen (dict or None))."""
    root = root or repo_root()
    kind, _, rest = spec.partition(':')
    if kind == 'crypto':
        base = os.path.join(root, 'results', 'crypto', rest)
        if not os.path.exists(base + '.qc'):
            raise TargetError('no such circuit in results/crypto: ' + rest)
        idx = json.load(open(os.path.join(root, 'results', 'crypto', 'index.json')))
        ent = [e for e in idx['circuits'] if e['name'] == rest]
        clean = bool(ent) and ent[0].get('model') == 'clean'
        ref = open(base + '.qc').read()
        gen = None
        m = re.match(r'lookup_c(\d+)_w(\d+)_', rest)
        if m:
            c, w = int(m.group(1)), int(m.group(2))
            gen = {'family': 'lookup', 'kind': 'xor', 'c': c, 'w': w, 'table': table_from_ref(ref, c, w)}
        return dict(id='crypto_' + rest, spec=spec, clean=clean, ref_text=ref, baseline_qasm=open(base + '.qasm').read(),
                    gen=gen)
    if kind in ('lookup', 'lookupadd'):
        m = re.fullmatch(r'c(\d+)w(\d+)s(\d+)', rest)
        if not m:
            raise TargetError('bad spec ' + spec)
        c, w, s = map(int, m.groups())
        from ..crypto import lookup as L
        T = _random_table(c, w, s)
        if kind == 'lookup':
            lc = L.davio_lookup(T, c, w)
            gen = {'family': 'lookup', 'kind': 'xor', 'c': c, 'w': w, 'table': T}
        else:
            lc = L.lookup_add(T, c, w, seed=0)
            gen = {'family': 'lookup', 'kind': 'add', 'c': c, 'w': w, 'table': T}
        ref = lc.reference_qc()
        return dict(id='%s_c%dw%ds%d' % (kind, c, w, s), spec=spec, clean=True, ref_text=ref,
                    baseline_qasm=pad_qubits(lc.qasm, _ref_wires(ref)), gen=gen)
    if kind == 'qc':
        path, _, flag = rest.partition(':')
        p = path if os.path.isabs(path) else os.path.join(root, path)
        if not os.path.exists(p):
            raise TargetError('no such file: ' + path)
        return dict(id='qc_' + re.sub(r'[^A-Za-z0-9_]+', '_', os.path.basename(path)[:-3] if path.endswith('.qc') else path),
                    spec=spec, clean=(flag == 'clean'), ref_text=open(p).read(), baseline_qasm=None, gen=None)
    raise TargetError('unknown target spec: ' + spec)


def register(state, spec, root=None):
    """Store the target definition (reference copy by hash). Returns the definition dict (with ref_sha, ref path)."""
    t = resolve(spec, root)
    ref_sha = state.put(t['ref_text'].encode(), '.qc')
    d = dict(id=t['id'], spec=spec, clean=t['clean'], ref_sha=ref_sha, gen=t['gen'], has_baseline=t['baseline_qasm'] is not None)
    try:
        old = state.load_target(t['id'])
        if old['ref_sha'] != ref_sha:
            raise TargetError('target %s already exists with another reference (hash differs)' % t['id'])
    except FileNotFoundError:
        pass
    state.save_target(t['id'], d)
    return d, t['baseline_qasm']
