"""Patched Qualtran copies: exactly the proven files, installed Qualtran untouched, hashes as in the proof rounds."""
import json
import os
import subprocess
import sys
from hashlib import sha256
from pathlib import Path

import pytest

pytest.importorskip("qualtran")

from qoptct.foreign.qualtran import patches
from qoptct.foreign.qualtran.patches import KINDS, PatchError, patched_env, patched_qualtran, verify_copy

# SHA-256 of the patched files as they were proved (copied from the round directories C/qlocal, E/qlocal, E/qlocal4).
EXPECTED = {
    "bloqs/data_loading/qrom_base.py": "fa02bc765c0aea68538fc2bb383d363bdc4dbfd111b9ff46fea8d0999f173f9e",
    "bloqs/data_loading/select_swap_qrom.py": "6ff73893944bdfd58c0ae150d5a65ebede6ee4ce1430ae0d52133de4ada1c116",
    "bloqs/cryptography/ecc/ec_add.py": "73bdcb891a4bae73a3b978d713e66a1d0d8825add323febff0af6e0d02c89cc2",
}
EXPECTED_MM = {"c": "c648332760a6f7c114ecedfa9415cfecfc6acae28a901fdb86fe5343461b7432",
               "e": "bc489ab85232a20fe83c78650b8bed36c7dc0ff1208baddbcab03b25d986296b"}
MM = "bloqs/mod_arithmetic/mod_multiplication.py"


@pytest.fixture(autouse=True)
def _cache(tmp_path, monkeypatch):
    monkeypatch.setenv("QOPT_CACHE", str(tmp_path))


def _file_hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def test_kinds():
    assert set(KINDS) == {"c", "e", "e4"}


@pytest.mark.parametrize("kind,files", [("c", ["bloqs/data_loading/qrom_base.py", "bloqs/data_loading/select_swap_qrom.py", MM]),
                                        ("e", ["bloqs/cryptography/ecc/ec_add.py", MM]),
                                        ("e4", ["bloqs/cryptography/ecc/ec_add.py", MM, "bloqs/data_loading/qrom_base.py",
                                                "bloqs/data_loading/select_swap_qrom.py"])])
def test_copy_has_proven_hashes_and_other_files_are_untouched(kind, files):
    import importlib.util
    installed = Path(list(importlib.util.find_spec("qualtran").submodule_search_locations)[0])
    base = patched_qualtran(kind)
    assert (base / "qualtran").is_dir() and base != installed.parent
    got = verify_copy(kind, base)
    assert set(got) == set(files)
    for rel in files:
        want = EXPECTED_MM["c" if kind == "c" else "e"] if rel == MM else EXPECTED[rel]
        assert _file_hash(base / "qualtran" / rel) == want
    # the installed Qualtran is not modified, and all non-patched files of the copy are identical to it
    for rel in installed.rglob("*.py"):
        r = rel.relative_to(installed).as_posix()
        if "__pycache__" in r or r in files:
            continue
        assert _file_hash(base / "qualtran" / r) == _file_hash(rel), r
    for rel in files:
        assert _file_hash(installed / rel) != _file_hash(base / "qualtran" / rel)


def test_copy_is_reused_and_repaired(tmp_path):
    first = patched_qualtran("c")
    stamp = (first / "qualtran" / MM).stat().st_mtime_ns
    assert patched_qualtran("c") == first and (first / "qualtran" / MM).stat().st_mtime_ns == stamp
    (first / "qualtran" / MM).write_text("# tampered\n")
    again = patched_qualtran("c")
    assert _file_hash(again / "qualtran" / MM) == EXPECTED_MM["c"]


def test_tampered_copy_is_refused():
    base = patched_qualtran("e")
    (base / "qualtran" / "bloqs/cryptography/ecc/ec_add.py").write_text("x = 1\n")
    with pytest.raises(PatchError):
        verify_copy("e", base)


def test_unknown_kind():
    with pytest.raises(PatchError):
        patched_qualtran("nope")


def test_subprocess_sees_the_patched_qualtran_and_the_env_switch():
    code = ("import qualtran, qualtran.bloqs.mod_arithmetic.mod_multiplication as m;"
            "print(qualtran.__file__); print(m._FINAL_K(7))")
    out = subprocess.run([sys.executable, "-c", code], env=patched_env("c"), capture_output=True, text=True, check=True)
    path, k = out.stdout.split("\n")[:2]
    assert "qualtran_patched" in path and k == "-7"
    env = dict(patched_env("c"), QLOCAL_ORIG_SIGN="1")
    out = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True)
    assert out.stdout.split("\n")[1] == "7"
    plain = subprocess.run([sys.executable, "-c", "import qualtran; print(qualtran.__file__)"], capture_output=True,
                           text=True, check=True)
    assert "qualtran_patched" not in plain.stdout


def test_pr1997_fix_changes_the_montgomery_final_reduction():
    """With the fix the final reduction adds -p; the env switch QLOCAL_ORIG_SIGN=1 restores the original +p."""
    code = ("from qualtran.bloqs.mod_arithmetic.mod_multiplication import DirtyOutOfPlaceMontgomeryModMul as M;"
            "b=M(bitsize=5,window_size=1,mod=13).decompose_bloq();"
            "ks=sorted({getattr(getattr(i.bloq,'subbloq',None),'k',None) for i in b.bloq_instances}-{None});"
            "print(ks)")
    fixed = subprocess.run([sys.executable, "-c", code], env=patched_env("c"), capture_output=True, text=True)
    assert fixed.returncode == 0, fixed.stderr[-500:]
    assert fixed.stdout.strip() == "[-13]", fixed.stdout
    orig = subprocess.run([sys.executable, "-c", code], env=dict(patched_env("c"), QLOCAL_ORIG_SIGN="1"),
                          capture_output=True, text=True)
    assert orig.stdout.strip() == "[13]", orig.stdout
    plain = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert plain.stdout.strip() == "[13]", plain.stdout
