"""F1: check the (wrong) hwb6 51-T measurement-based circuit with the sum-over-paths checker.
Neither hwb6.qc nor circuits derived from it are redistributed with this release (see bench_qc/SOURCES.md):
  sh scripts/fetch_benchmarks.sh bench_extra           # downloads hwb6.qc from the Feynman repository
  python docs/findings/f1_hwb6_51T.py CIRCUIT.qasm bench_extra/hwb6.qc
Expected for the 51 T / 535 CX circuit: ok=False (rejected in every sampled measurement branch)."""
import sys

import qopt  # noqa: F401  (sets up the package paths)
from qopt.contrib.c10 import verify_big

if len(sys.argv) != 3:
    sys.exit(__doc__)
r = verify_big.check(sys.argv[1], sys.argv[2])
print({k: r[k] for k in r if k in ("ok", "F")})
