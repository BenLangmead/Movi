#!/bin/bash
# Regression test: a query's result must not depend on how many threads or strands
# ran it. For every query the strand/sequential paths serve -- PML, ZML and whole-read
# exact count -- the per-read results at `-t 1` must equal those from several threaded
# configurations.
#
# Why this exists: the byte-gate and the coroutine regressions both run at `-t 1`, so
# neither can see a race between worker threads. A concrete example this test catches:
# std::cerr is tied to std::cout, so emitting a progress line flushes std::cout; doing
# that outside the output lock lets one thread flush while another is mid-write, which
# corrupts stdout only when more than one thread is running. That defect passes every
# other check in the suite.
#
# Results are compared per read id rather than as raw files, because emit order across
# threads is not defined -- only the content per read is. Trailing run-summary stats
# that some queries append to stdout are ignored for the same reason (ff_count is
# documented as not meaningful on the prefetch path).
#
# Usage:
#   regression_thread_determinism.sh [BUILD_DIR]
# Env overrides:
#   MOVI_TEST_SEP_INDEX=<index dir, mode 6>   (required unless tests_data has one)
#   MOVI_TEST_READS=<reads fasta/fastq>
#   MOVI_TEST_CONFIGS="t:s t:s ..."           (defaults below)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
BUILD="${1:-$ROOT/build-release}"
BIN="$BUILD/bin"
DATA="$ROOT/tests_data"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

MV="$BIN/movi-regular-thresholds"
[ -x "$MV" ] || { echo "ERROR: missing $MV (build first)"; exit 2; }

READS="${MOVI_TEST_READS:-$DATA/reads.fasta}"
IDX="${MOVI_TEST_SEP_INDEX:-}"
[ -n "$IDX" ] || { echo "ERROR: set MOVI_TEST_SEP_INDEX to a mode-6 index directory"; exit 2; }
[ -d "$IDX" ] || { echo "ERROR: index directory not found: $IDX"; exit 2; }

# Thread:strand pairs. The mix matters: a single strand with many threads and many
# strands with many threads stress different parts of the scheduler.
CONFIGS="${MOVI_TEST_CONFIGS:-8:8 8:1 16:4 48:8 48:1}"

# Compare two result files by read id. Handles both shapes the strand path emits:
# '>id' followed by a values line (PML, ZML), and one tab-separated line per read
# beginning with the id (count).
compare_by_id () {  # $1=reference file  $2=candidate file  $3=label
  python3 - "$1" "$2" "$3" <<'PY'
import sys
ref_path, cand_path, label = sys.argv[1], sys.argv[2], sys.argv[3]

def parse(path):
    d = {}
    rid = None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            rid = line[1:]
        elif rid is not None:
            d[rid] = line          # values line for the id just seen
            rid = None
        elif "\t" in line:
            d[line.split("\t", 1)[0]] = line   # count-style row
    return d

a, b = parse(ref_path), parse(cand_path)
missing = sorted(set(a) - set(b))
extra   = sorted(set(b) - set(a))
differ  = sorted(r for r in set(a) & set(b) if a[r] != b[r])
if not a:
    print("  FAIL [%s]: reference produced no per-read records" % label); sys.exit(1)
if missing or extra or differ:
    print("  FAIL [%s]: reads=%d missing=%d extra=%d value_mismatches=%d"
          % (label, len(a), len(missing), len(extra), len(differ)))
    for r in (missing[:3] + extra[:3] + differ[:3]):
        print("    id=%s ref=%s cand=%s" % (r, a.get(r, "<none>")[:60], b.get(r, "<none>")[:60]))
    sys.exit(1)
print("  PASS [%s]: %d reads identical" % (label, len(a)))
PY
}

status=0
for Q in --pml --zml --count; do
  tag=${Q#--}
  echo "[thread-determinism] [$tag] reference run at -t 1 ..."
  "$MV" query --index "$IDX" --read "$READS" $Q --stdout -t 1 -s 8 > "$WORK/$tag.ref" 2>/dev/null
  for cfg in $CONFIGS; do
    t=${cfg%%:*}; s=${cfg##*:}
    "$MV" query --index "$IDX" --read "$READS" $Q --stdout -t "$t" -s "$s" > "$WORK/$tag.$t.$s" 2>/dev/null
    compare_by_id "$WORK/$tag.ref" "$WORK/$tag.$t.$s" "$tag t=$t s=$s" || status=1
  done
done

if [ $status -ne 0 ]; then
  echo "[thread-determinism] FAILED: results depend on thread/strand count"
  exit 1
fi
echo "[thread-determinism] all queries identical across $(echo $CONFIGS | wc -w) threaded configurations"
