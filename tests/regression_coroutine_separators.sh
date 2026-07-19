#!/bin/bash
# Regression test: coroutine PML (`movi query --pml --coroutine`) must produce
# EXACTLY the same pseudo-matching lengths as the sequential path on a
# SEPARATOR-using index.
#
# Why this exists: the coroutine PML body re-implements the PML inner loop inline
# (suspend points can't cross function calls), and that copy has historically
# diverged from the sequential path in subtle ways. One bug was a separator-run
# threshold off-by-one (values[r_ch_id] vs values[r_ch_id_adj]) that produced
# *silent* single-position errors, only when repositioning off a separator run.
# A single-sequence reference does NOT create separator runs (only the terminal
# $), so a meaningful test needs a MULTI-sequence reference built with --separators.
#
# Usage:
#   regression_coroutine_separators.sh [BUILD_DIR]
# Env overrides (skip the build pipeline; run against a prebuilt separator index):
#   MOVI_TEST_SEP_INDEX=<dir with index.movi built --separators, mode 6>
#   MOVI_TEST_READS=<reads fasta/fastq>
#
# Requires: movi-regular-thresholds (always); the 'movi' launcher and the
# threshold build pipeline (pfp/r-permute) only when building the index from
# scratch (same requirement as the existing build/pml Catch2 tests).
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

if [ -n "${MOVI_TEST_SEP_INDEX:-}" ]; then
  IDX="$MOVI_TEST_SEP_INDEX"
  echo "[regression] using prebuilt separator index: $IDX"
else
  # Locate the 'movi' launcher (drives the threshold build pipeline).
  LAUNCH=""
  for c in "$BUILD/movi" "$BIN/movi"; do [ -x "$c" ] && LAUNCH="$c" && break; done
  [ -n "$LAUNCH" ] || { echo "ERROR: 'movi' launcher not found under $BUILD (needed to build the index)"; exit 2; }
  # Build a multi-sequence reference so --separators yields real separator runs.
  python3 - "$DATA/ref.fasta" "$WORK/ref_multi.fasta" <<'PY'
import sys
seq="".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
n=10; L=len(seq)//n
with open(sys.argv[2],"w") as o:
    for i in range(n):
        o.write(">chunk%d\n%s\n"%(i, seq[i*L:(i+1)*L] if i<n-1 else seq[i*L:]))
PY
  IDX="$WORK/idx_sep"
  echo "[regression] building separator index via launcher..."
  "$LAUNCH" build --type regular-thresholds --index "$IDX" --fasta "$WORK/ref_multi.fasta" --separators --verify >/dev/null 2>&1 \
    || { echo "ERROR: separator index build failed (threshold pipeline available?)"; exit 2; }
fi

echo "[regression] running coroutine PML and sequential PML ..."
# Both emit `>id\n<space-separated PMLs>` to stdout via the shared output_base_stats.
# The coroutine computes the same forward-read PMLs as the sequential path, so
# the per-read value arrays must match exactly (compared by read id, so the emit
# order -- reorder-buffer input order vs sequential batch order -- does not matter).
"$MV" query --index "$IDX" --read "$READS" --pml --coroutine --stdout -s 8 -t 1 > "$WORK/co.out"  2>/dev/null
"$MV" query --index "$IDX" --read "$READS" --pml             --stdout       -t 1 > "$WORK/seq.out" 2>/dev/null

python3 - "$WORK/co.out" "$WORK/seq.out" <<'PY'
import sys
def parse(path):
    d={}; rid=None
    for line in open(path):
        if line.startswith(">"): rid=line[1:].strip()
        elif rid is not None: d[rid]=list(map(int,line.split())); rid=None
    return d
co=parse(sys.argv[1]); seq=parse(sys.argv[2])
shared=[r for r in co if r in seq]
unmatched=[r for r in seq if r not in co]+[r for r in co if r not in seq]
bad=[r for r in shared if co[r]!=seq[r]]
print("reads: coroutine=%d sequential=%d shared=%d unmatched_ids=%d value_mismatches=%d"
      % (len(co),len(seq),len(shared),len(unmatched),len(bad)))
if not shared or unmatched or bad:
    print("FAIL: coroutine PML diverges from sequential on a separator index")
    for r in bad[:5]:
        a,b=co[r],seq[r]
        i=next((j for j in range(min(len(a),len(b))) if a[j]!=b[j]),-1)
        print("  %s len=%d first_diff_idx=%d co=%s seq=%s"%(r,len(a),i,a[i:i+1],b[i:i+1]))
    sys.exit(1)
print("PASS: coroutine PML matches sequential exactly on the separator index (%d reads)"%len(shared))
PY
