"""ProjectQ bridge tests: skipped without ProjectQ (optional extra ``projectq``), except the import-guard test."""
import os
import pathlib

import pytest

REPO = pathlib.Path(__file__).resolve().parents[2]
RESULTS = REPO / 'results' / 'projectq'


def have_projectq():
    try:
        import projectq  # noqa: F401
        return True
    except ImportError:
        return False


needs_projectq = pytest.mark.skipif(not have_projectq(), reason='ProjectQ is not installed (extra "projectq")')


def binary_ecc_source():
    p = os.environ.get('QOPT_BINARY_ECC')
    return p if p and (pathlib.Path(p) / 'Multiplication').is_dir() else None


needs_binary_ecc = pytest.mark.skipif(binary_ecc_source() is None, reason='QOPT_BINARY_ECC does not point at a Binary_ECC folder')
