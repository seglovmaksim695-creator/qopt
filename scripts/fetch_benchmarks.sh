#!/bin/sh
# Download benchmark circuits that are not redistributed with qopt (their headers forbid redistribution)
# from the Feynman repository, for local use only. Usage: sh scripts/fetch_benchmarks.sh [DEST_DIR]
set -e
DEST=${1:-bench_extra}
BASE=https://raw.githubusercontent.com/meamy/feynman/master/benchmarks/qc
mkdir -p "$DEST"
for f in hwb6; do
  curl -fsSL "$BASE/$f.qc" -o "$DEST/$f.qc"
  echo "downloaded $DEST/$f.qc (read its copyright header: redistribution is not permitted)"
done
