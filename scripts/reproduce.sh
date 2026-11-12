#!/bin/sh
# Example: run qopt on all shipped benchmarks.
#   sh scripts/reproduce.sh [MODE] [OUTDIR]       MODE: cx (default) | net | balanced | t | t-clean | all
# Writes OUTDIR/<mode>/<name>.(qc|qasm) and one JSON report per circuit, then prints a summary.
# Time: 'cx' takes minutes for small circuits and up to ~1 h for the largest ones (time limits per rewriting
# program: QOPT_PROG_SEC=120 s by default); 'balanced'/'t' take minutes per circuit (<= 22 qubits only).
# The published tables in results/ were obtained with this pipeline plus longer development searches
# (see results/TABLE.md, column "source"); a single run of this script is not expected to match every entry.
set -e
MODE=${1:-cx}
OUT=${2:-reproduce_out}
HERE=$(cd "$(dirname "$0")/.." && pwd)
BENCH="$HERE/bench_qc"
QOPT=${QOPT:-qopt}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

run_mode() {
  m=$1; mkdir -p "$OUT/$m"
  for f in "$BENCH"/*.qc; do
    n=$(basename "$f" .qc)
    case $m in
      cx)       args="";                         ext=qc ;;
      net)      args="--net";                    ext=qc ;;
      balanced) args="--mode balanced";          ext=qasm ;;
      t)        args="--mode t";                 ext=qasm ;;
      t-clean)  args="--mode t --clean";         ext=qasm ;;
    esac
    echo "== $m $n"
    $QOPT optimize "$f" -o "$OUT/$m/$n.$ext" --report "$OUT/$m/$n.json" $args > "$OUT/$m/$n.log" 2>&1 || echo "   failed (see $OUT/$m/$n.log)"
  done
  python3 - "$OUT/$m" <<'EOF'
import json, sys, glob, os
rows = []
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.json"))):
    r = json.load(open(f)); n = os.path.basename(f)[:-5]
    best = (r["t_out"], r["cx_out"], 0)
    for k in ("balanced", "tgad", "and"):
        x = r.get(k) or {}
        if x.get("used") or x.get("chosen"):
            best = (x.get("t"), x.get("cx"), x.get("extra_qubits"))
    rows.append((n, r["t_in"], r["cx_in"], r["t_out"], r["cx_out"], r["verified"]) + best)
print(f"{'circuit':16s} {'T_in':>5} {'CX_in':>6} {'T_unit':>6} {'CX_unit':>7} proven  T_out CX_out extra")
for r in rows:
    print(f"{r[0]:16s} {r[1]:5d} {r[2]:6d} {r[3]:6d} {r[4]:7d} {str(r[5]):6s} {r[6]!s:>6} {r[7]!s:>6} {r[8]!s:>5}")
EOF
}

if [ "$MODE" = all ]; then
  for m in cx net balanced t t-clean; do run_mode $m; done
else
  run_mode "$MODE"
fi
