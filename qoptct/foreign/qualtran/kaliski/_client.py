"""Hashes of the client that makes (and later judges) the Kaliski evidence.

The client is: the six source files of this sub-package that the evidence depends on, the whole strict bridge
(``qoptct/foreign/qualtran/bridge``) and the unchanged qoptct provers (``qoptct/verify``).  ``run.manifest`` records these
hashes next to every piece of evidence; ``verdict`` accepts canary evidence only when ALL of them equal the hashes of the
code that is judging.  A change in the bridge or in a prover therefore makes old canaries invalid until they are regenerated.
"""
import glob
import hashlib
import json
import os

USED = ('lower.py', 'chain.py', 'pq.py', 'run.py', 'lemmas_z3.py', 'canary.py', '_client.py')


def client_groups():
    """-> {'src': {files, sha256, each}, 'bridge': {files, sha256}, 'qoptct': {files, sha256, each}} of the code on disk now."""
    pkg = os.path.dirname(os.path.abspath(__file__))                         # .../qoptct/foreign/qualtran/kaliski
    qoptct_root = os.path.dirname(os.path.dirname(os.path.dirname(pkg)))     # .../qoptct
    groups = dict(src=(os.path.dirname(pkg), sorted(os.path.join(pkg, f) for f in USED)),
                  bridge=(os.path.dirname(pkg), sorted(glob.glob(os.path.join(os.path.dirname(pkg), 'bridge', '*.py')))),
                  qoptct=(os.path.dirname(qoptct_root),
                          sorted(glob.glob(os.path.join(qoptct_root, 'verify', '**', '*.py'), recursive=True))))
    out = {}
    for g, (root, files) in groups.items():
        h = {os.path.relpath(f, root): hashlib.sha256(open(f, 'rb').read()).hexdigest() for f in files}
        out[g] = dict(files=len(h), sha256=hashlib.sha256(json.dumps(h, sort_keys=True).encode()).hexdigest())
        if g in ('src', 'qoptct'):
            out[g]['each'] = h
    return out


def same_client(manifest, current=None):
    """True when ``manifest`` (the ``manifest`` / ``manifest_before`` of a piece of evidence) was made by exactly the code on
    disk now: every file of ``src`` and ``qoptct`` individually, and the group hash of the bridge (with its file count)."""
    cur = current or client_groups()
    m = manifest or {}
    for g in ('src', 'bridge', 'qoptct'):
        a, b = m.get(g) or {}, cur[g]
        if not a or a.get('sha256') != b['sha256'] or a.get('files') != b['files']:
            return False
        if g != 'bridge' and a.get('each') != b['each']:
            return False
    return True
