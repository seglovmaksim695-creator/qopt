"""Montgomery, window 0 by sub-nodes (the fallback for w <= 2 when the proof of the whole window 0 returns UNKNOWN; audit AUD7, remark 2).

The whole-window proof is forced to UNKNOWN (as the auditor did).  The fallback must then PROVE the control (n = 16, w = 1 and 2) and
must NOT prove any of three damaged programs (the control X of window 0 swapped with its target at the start, the middle and the end
of the window).  Structural checks do not see these damages: only the proof of the sub-node does."""
import json
import os
import pickle
import subprocess
import sys

import pytest

pytest.importorskip("qualtran")
pytest.importorskip("qoptct.verify.bdd")

from qoptct.foreign.qualtran._variant import child_env
from qoptct.foreign.qualtran.montgomery import prime_for

LOWER = r'''
import sys
from qoptct.foreign.qualtran.montgomery import _stages as S, mutants as MUT
tag, out, n, w, p = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
def my_apply(t, P, segs, info, base0, nw):
    segs = [list(s) for s in segs]
    a_ = info[0]['span'][0] - base0
    ops = list(segs[0][2])
    idx = [i for i, o in enumerate(ops) if i >= a_ and o[0] == 'cx']
    i = idx[0] if t == 'W0a' else idx[-1] if t == 'W0b' else idx[len(idx) // 2]
    o = ops[i]; ops[i] = (o[0], o[2], o[1]) + tuple(o[3:]); segs[0][2] = ops
    return [tuple(s) for s in segs], info
MUT.apply = my_apply
S.stage_lower(n, w, p, 'c', out, mutant=(None if tag == 'ctrl' else tag))
'''
PROVE = r'''
import sys
from qoptct.foreign.qualtran.montgomery import _stages as S, _core as C
real = C.run_segment
def fake(label, kind, *a, **k):
    if label == 'win0':
        return dict(label=label, kind=kind, status='UNKNOWN', detail='forced: the whole-window proof gives UNKNOWN')
    return real(label, kind, *a, **k)
C.run_segment = fake
S.stage_prove(sys.argv[1], sys.argv[2], 600, 2800, 7)
'''


def _run(tmp_path, tag, w, n=16):
    p = prime_for(n)
    lp, pp = str(tmp_path / f"lower_{tag}_{w}.pkl"), str(tmp_path / f"prove_{tag}_{w}.json")
    env = child_env("c")
    env["QOPT_CACHE"] = str(tmp_path)
    r1 = subprocess.run([sys.executable, "-B", "-c", LOWER, tag, lp, str(n), str(w), str(p)], env=env, capture_output=True, text=True)
    assert r1.returncode == 0, r1.stderr[-600:]
    r2 = subprocess.run([sys.executable, "-B", "-c", PROVE, lp, pp], env=env, capture_output=True, text=True)
    assert os.path.exists(pp), r2.stderr[-600:]
    return json.load(open(pp)), pickle.load(open(lp, "rb"))


def _proved(out, D):
    """The part of the verdict of prove_montgomery that depends on window 0 and the structural checks."""
    ck = D["checks"]
    structural = (ck["glue"]["ok"] and ck["bookkeeping"]["ok"] and ck["windows_ok"] and ck["prefix"]["ok"] and ck["same_hash"]
                  and not ck["tiles_bad"] and not (D.get("e2e") or {}).get("mismatches"))
    return bool(structural and out.get("win0_ok") and not out["reasons"])


@pytest.mark.parametrize("w", [1, 2])
def test_fallback_proves_control_and_rejects_three_damages(tmp_path, w):
    out, D = _run(tmp_path, "ctrl", w)
    assert out["win0_split"] is not None, "the fallback path (window 0 by sub-nodes) did not run"
    assert out["win0_split"]["ok"] and _proved(out, D), out["reasons"]
    for tag in ("W0a", "W0m", "W0b"):
        out, D = _run(tmp_path, tag, w)
        assert out["win0_split"] is not None
        assert not _proved(out, D), f"damage {tag} (w={w}) was PROVED by the fallback"
