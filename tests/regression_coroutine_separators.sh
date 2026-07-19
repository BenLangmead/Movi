#!/bin/bash
# Regression test: the coroutine latency-hiding path (`movi query --coroutine`)
# must produce exactly the same results as the sequential path on a separator-using
# index, for every query it routes: PML, k-mer presence, and MEM. PML is checked
# per-read (exact value arrays by id); k-mer and MEM are checked as sorted-line
# multisets (order-independent, since the coroutine reorder buffer emits in input
# order while the sequential path emits in batch order).
#
# Why a separator index: the coroutine PML body inlines the PML inner loop (suspend
# points can't cross function calls), so it is easy for that inline copy to diverge
# from the sequential path in ways that only surface when repositioning off a
# separator run -- e.g. a separator-run threshold off-by-one (values[r_ch_id] vs
# values[r_ch_id_adj]) causes silent single-position errors. A single-sequence
# reference does not create separator runs (only the terminal $), so a meaningful
# test needs a multi-sequence reference built with --separators.
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

# k for the k-mer sub-test; reads must be at least this long. Override via env.
K="${MOVI_TEST_K:-21}"

echo "[regression] [PML] running coroutine PML and sequential PML ..."
# Both emit `>id\n<space-separated PMLs>` to stdout via the shared output_base_stats.
# The coroutine computes the same forward-read PMLs as the sequential path, so
# the per-read value arrays must match exactly (compared by read id, so the emit
# order -- reorder-buffer input order vs sequential batch order -- does not matter).
"$MV" query --index "$IDX" --read "$READS" --pml --coroutine --stdout -s 8 -t 1 > "$WORK/pml.co"  2>/dev/null
"$MV" query --index "$IDX" --read "$READS" --pml             --stdout       -t 1 > "$WORK/pml.seq" 2>/dev/null

python3 - "$WORK/pml.co" "$WORK/pml.seq" <<'PY'
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
print("  reads: coroutine=%d sequential=%d shared=%d unmatched_ids=%d value_mismatches=%d"
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

# MEM and plain k-mer presence are also routed through the coroutine (movi.cpp
# dispatch), both via the shared output_mems / output_kmers emitters. Their per-read
# content must be byte-identical to the sequential path. These are compared in file
# mode (-o), not --stdout: the sequential --stdout path appends a trailing run-summary
# stats block (total_kmers/backward_search_*/...) that the coroutine does not emit, and
# which would otherwise show up as spurious "extra lines"; the -o output files hold only
# the per-read records. The emit order can still differ (reorder-buffer input order vs
# sequential batch order), so compare the sorted line multiset -- the project's
# order-independent "byte-identical content" standard.
compare_file_mode () {  # $1=label  $2=output-file suffix  $3...=query flags (no --coroutine/-o/-t/-s)
  local label="$1" suf="$2"; shift 2
  echo "[regression] [$label] running coroutine vs sequential (file mode) ..."
  "$MV" query --index "$IDX" --read "$READS" "$@" --coroutine -s 8 -t 1 -o "$WORK/$label.co"  >/dev/null 2>&1
  "$MV" query --index "$IDX" --read "$READS" "$@"             -t 1 -o "$WORK/$label.seq" >/dev/null 2>&1
  local cof="$WORK/$label.co$suf" sef="$WORK/$label.seq$suf"
  [ -f "$cof" ] || { echo "FAIL: coroutine $label produced no output file ($cof)"; exit 1; }
  [ -f "$sef" ] || { echo "FAIL: sequential $label produced no output file ($sef)"; exit 1; }
  sort "$cof" > "$cof.sorted"; sort "$sef" > "$sef.sorted"
  local nco nseq
  nco=$(wc -l < "$cof.sorted"); nseq=$(wc -l < "$sef.sorted")
  echo "  records: coroutine=$nco sequential=$nseq"
  if ! diff -q "$cof.sorted" "$sef.sorted" >/dev/null; then
    echo "FAIL: coroutine $label diverges from sequential on a separator index"
    diff "$cof.sorted" "$sef.sorted" | head -10
    exit 1
  fi
  [ "$nco" -gt 0 ] || { echo "FAIL: coroutine $label produced no output"; exit 1; }
  echo "PASS: coroutine $label matches sequential exactly (file content, $nco records)"
}

compare_file_mode kmer ".kmers.$K" --kmer -k "$K"
compare_file_mode mem  ".mems"     --mem --ftab-k 10 --min-mem-length 25

# Scope note: this test covers the coroutine queries that run on a plain
# regular-thresholds --separators index (PML, k-mer presence, MEM). The remaining
# routed coroutine query, bitvector count (--kmer-count --kmer-bv), needs a kmer-bv
# index (movi build-kmerbv, with the per-k B_k/C_k structures) that this lightweight
# test does not build; its coroutine-vs-sequential equivalence is validated by the
# cluster byte-gate baseline.

echo "[regression] all coroutine-vs-sequential checks passed (PML, k-mer, MEM)"
