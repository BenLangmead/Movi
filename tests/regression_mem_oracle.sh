#!/bin/bash
#
# Author: Ben Langmead (ben.langmead@gmail.com)
# Date: Sep 22, 2026
#
# Distributed under the GPL3 license.
# See accompanying LICENSE or https://opensource.org/license/gpl-3-0
#
# Regression test: `movi query --mem` must report exactly the MEMs a brute-force oracle
# finds (tests/mem_oracle.py), for every minimum length L and ftab depth Movi accepts,
# and every execution path must print the same records.
#
# tests/mem_differential.py builds small random references (one to five records, with
# planted repeats) and random queries, including queries copied from the reference with
# substitutions and Ns and queries that straddle record boundaries. For each L in
# {1, 2, 3, 4, 6, 10} and ftab-k in {2, 3, 4, 6} (ftab-k above L is rejected by Movi
# for L >= 2, so for L >= 2 this covers ftab-k below and equal to L, and L = 1 covers
# ftab-k above it) it runs the sequential path
# with -t 1 and -t 4, the coroutine path with -t 1 and -t 4, and a strand-scheduler
# request, and requires:
#   * all paths print the same records, and seq1 and co1 the same bytes;
#   * the records equal the oracle's (MEMs, counts and covered/length);
#   * every MEM query on an index built without --separators is refused, since such
#     an index lets matches run across record boundaries;
#   * --mem --kmer-out membership views equal direct --kmer presence queries.
# It also checks that -L is the short form of --min-mem-length.
#
# Usage:
#   regression_mem_oracle.sh [BUILD_DIR]
# Env overrides:
#   MOVI_MEM_SEEDS="1 2 3"   seeds to run (one set of scenarios per seed)
#   MOVI_MEM_SCENARIOS=6     random references per separator setting and seed
#
# Requires: the 'movi' launcher and the threshold build pipeline (pfp) to build the
# toy indexes, movi-regular-thresholds, and python3.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
BUILD="${1:-$ROOT/build-release}"
MV="$BUILD/bin/movi-regular-thresholds"
[ -x "$MV" ] || { echo "ERROR: missing $MV (build first)"; exit 2; }
[ -x "$BUILD/movi" ] || { echo "ERROR: missing launcher $BUILD/movi"; exit 2; }
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

SEEDS="${MOVI_MEM_SEEDS:-1 2}"
SCEN="${MOVI_MEM_SCENARIOS:-6}"
status=0
for seed in $SEEDS; do
  echo "[regression] [mem-oracle] seed $seed ..."
  if ! python3 "$HERE/mem_differential.py" --build-dir "$BUILD" --seed "$seed" \
        --scenarios "$SCEN" --queries 60 --lengths 1,2,3,4,6,10 --ftab-ks 2,3,4,6 \
        --require-refusal --kmer-out-check --work "$WORK/seed$seed"; then
    status=1
  fi
done

echo "[regression] [mem-oracle] -L is the short form of --min-mem-length ..."
IDX=$(ls -d "$WORK"/seed*/sep_0/idx | head -1)
READS="$(dirname "$IDX")/reads.fa"
"$MV" query --index "$IDX" --read "$READS" --mem --ftab-k 3 -L 3 -t 1 -o "$WORK/short" >/dev/null 2>&1
"$MV" query --index "$IDX" --read "$READS" --mem --ftab-k 3 --min-mem-length 3 -t 1 -o "$WORK/long" >/dev/null 2>&1
if cmp -s "$WORK/short.mems" "$WORK/long.mems"; then
  echo "PASS: -L 3 and --min-mem-length 3 agree"
else
  echo "FAIL: -L 3 and --min-mem-length 3 disagree"; status=1
fi

[ $status -eq 0 ] && echo "[regression] all MEM oracle checks passed" || echo "[regression] MEM oracle checks FAILED"
exit $status
