#!/bin/bash
# Regression test: the coroutine PML binary (movi-co) must produce EXACTLY the
# same pseudo-matching lengths as production Movi on a SEPARATOR-using index.
#
# Why this exists: movi-co re-implements the PML inner loop inline, and that
# copy has diverged from production in subtle ways. One bug was a separator-run
# threshold off-by-one (values[r_ch_id] vs values[r_ch_id_adj]) that produced
# *silent* single-position errors, only when repositioning off a separator run.
# A single-sequence reference does NOT create separator runs (only the
# terminal $), so a meaningful test needs a MULTI-sequence reference built with
# --separators.
#
# Usage:
#   regression_movi_co_separators.sh [BUILD_DIR]
# Env overrides (skip the build pipeline; useful to run against a prebuilt
# separator index, e.g. a large one on scratch):
#   MOVI_TEST_SEP_INDEX=<dir with index.movi built --separators, mode 6>
#   MOVI_TEST_READS=<reads fasta/fastq>
#
# Requires: movi-co and movi-regular-thresholds (always); the 'movi' launcher
# and the threshold build pipeline (pfp/r-permute) only when building the index
# from scratch (same requirement as the existing build/pml Catch2 tests).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
BUILD="${1:-$ROOT/build-release}"
BIN="$BUILD/bin"
DATA="$ROOT/tests_data"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

CO="$BIN/movi-co"
MV="$BIN/movi-regular-thresholds"
[ -x "$CO" ] || { echo "ERROR: missing $CO (build first)"; exit 2; }
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
  # Build a MULTI-sequence reference so --separators yields real separator runs.
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

echo "[regression] running movi-co and production (--reverse) ..."
"$CO" "$READS" "$IDX" 8 > "$WORK/co.out" 2>/dev/null
"$MV" query --index "$IDX" --read "$READS" --pml --reverse -s 8 -t 1 -o "$WORK/mv" >/dev/null 2>&1
"$MV" view --bpf "$WORK/mv.pml.bpf" > "$WORK/mv.view" 2>/dev/null

python3 - "$WORK/co.out" "$WORK/mv.view" <<'PY'
import sys
co={}
for line in open(sys.argv[1]):
    if line.startswith("#") or not line.strip(): continue
    p=line.split(); co[p[0]]=list(map(int,p[1:]))
mv={}; rid=None
for line in open(sys.argv[2]):
    if line.startswith(">"): rid=line[1:].strip()
    elif rid is not None: mv[rid]=list(map(int,line.split())); rid=None
shared=[r for r in co if r in mv]
unmatched=[r for r in mv if r not in co]+[r for r in co if r not in mv]
# movi-co reverses the read; production --reverse + view prints read order:
# accept exact OR reversed (orientation convention), reject anything else.
bad=[r for r in shared if not (co[r]==mv[r] or co[r]==mv[r][::-1])]
print("reads: movi-co=%d production=%d shared=%d unmatched_ids=%d value_mismatches=%d"
      % (len(co),len(mv),len(shared),len(unmatched),len(bad)))
if not shared or unmatched or bad:
    print("FAIL: movi-co diverges from production on a separator index")
    for r in bad[:5]:
        a,b=co[r],mv[r][::-1]
        i=next((j for j in range(min(len(a),len(b))) if a[j]!=b[j]),-1)
        print("  %s len=%d first_diff_idx=%d co=%s prod=%s"%(r,len(a),i,a[i:i+1],b[i:i+1]))
    sys.exit(1)
print("PASS: movi-co matches production exactly on the separator index (%d reads)"%len(shared))
PY
