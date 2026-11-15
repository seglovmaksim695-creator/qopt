"""Benchmark location for the module tests.

The tests read benchmark circuits from $QOPT_BENCH/<set>/<name>.qc (sets: std, large, crypto, extra). Without
$QOPT_BENCH, a temporary folder is used whose `std` set is the redistributable `bench_qc/` of this repository;
tests that need circuits of the other sets are then skipped."""
import os
import pathlib
import tempfile

if not os.environ.get("QOPT_BENCH"):
    _root = pathlib.Path(tempfile.mkdtemp(prefix="qopt_bench_"))
    _qc = pathlib.Path(__file__).resolve().parents[1] / "bench_qc"
    if _qc.is_dir():
        (_root / "std").symlink_to(_qc, target_is_directory=True)
    os.environ["QOPT_BENCH"] = str(_root)
