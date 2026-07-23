#!/bin/bash
# Regression test: coroutine PML (`movi query --pml --coroutine`) must produce exactly
# the same pseudo-matching lengths as the sequential path in every threshold index mode
# -- regular-thresholds (6), sampled-thresholds (7), and blocked-thresholds (8).
#
# Why this exists: the coroutine PML body uses the mode-portable move-structure
# primitives (get_c / reposition_thresholds / get_id) rather than a mode-6-specific
# inline, so it must stay byte-identical to the sequential path in all three modes.
#
# No pfp / grlBWT is needed: the test ships a checked-in, mode-independent PRECURSOR set
# (the run-length BWT heads/lengths + thresholds for a small reference) in
# tests_data/pml_modes/. Each mode's index is serialized from those precursors by the
# corresponding movi-*-thresholds binary (build --preprocessed), so the only build tools
# required are the Movi binaries themselves.
#
# Usage: regression_coroutine_pml_modes.sh [BUILD_DIR]   (BUILD_DIR default: <root>/build-release)
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
BUILD="${1:-$ROOT/build-release}"
BIN="$BUILD/bin"
DATA="$ROOT/tests_data/pml_modes"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

REF="$DATA/ref.fa"          # cleaned reference; its .bwt.heads/.bwt.len/.thr_pos sit beside it
READS="$DATA/reads.fa"
for f in "$REF" "$REF.bwt.heads" "$REF.bwt.len" "$REF.thr_pos" "$READS"; do
  [ -f "$f" ] || { echo "ERROR: missing fixture $f"; exit 2; }
done

# mode -> binary (all three are the same MODE=6/7/8 build of the query engine)
declare -A BINOF=( [6]=movi-regular-thresholds [7]=movi-sampled-thresholds [8]=movi-blocked-thresholds )

ran=0
fail=0
declare -A SEQOUT   # mode -> sequential PML output file (for the cross-mode check)

for m in 6 7 8; do
  MV="$BIN/${BINOF[$m]}"
  if [ ! -x "$MV" ]; then
    echo "[pml-modes] mode $m: ${BINOF[$m]} not built, skipping"
    continue
  fi
  echo "[pml-modes] mode $m (${BINOF[$m]}): serializing index from precursors ..."
  idx="$WORK/idx$m"; mkdir -p "$idx"
  if ! "$MV" build --fasta "$REF" --index "$idx" --preprocessed >"$WORK/build$m.log" 2>&1; then
    echo "FAIL: mode $m index build from precursors failed (see below)"; tail -5 "$WORK/build$m.log"; fail=1; continue
  fi
  [ -f "$idx/index.movi" ] || { echo "FAIL: mode $m produced no index.movi"; fail=1; continue; }

  "$MV" query -i "$idx" -r "$READS" --pml --coroutine -s 8 -t 1 --stdout > "$WORK/cor$m.txt" 2>/dev/null
  "$MV" query -i "$idx" -r "$READS" --pml             -t 1 --stdout > "$WORK/seq$m.txt" 2>/dev/null
  SEQOUT[$m]="$WORK/seq$m.txt"
  ran=$((ran+1))

  # Compare per-read PML arrays by read id (order-independent: the coroutine reorder
  # buffer emits in input order, the sequential path in batch order).
  if ! python3 - "$WORK/cor$m.txt" "$WORK/seq$m.txt" "$m" <<'PY'
import sys
def parse(p):
    d={}; rid=None
    for line in open(p):
        if line.startswith(">"): rid=line[1:].strip()
        elif rid is not None: d[rid]=line.split(); rid=None
    return d
co=parse(sys.argv[1]); sq=parse(sys.argv[2]); m=sys.argv[3]
shared=[r for r in co if r in sq]
unmatched=[r for r in sq if r not in co]+[r for r in co if r not in sq]
bad=[r for r in shared if co[r]!=sq[r]]
nz=sum(1 for r in sq for v in sq[r] if v!='0')
print("  mode %s: reads=%d shared=%d unmatched=%d value_mismatches=%d nonzero_PMLs=%d"
      % (m,len(sq),len(shared),len(unmatched),len(bad),nz))
if not shared or unmatched or bad:
    print("FAIL: coroutine PML diverges from sequential in mode %s" % m); sys.exit(1)
print("PASS: coroutine PML matches sequential in mode %s (%d reads)" % (m,len(shared)))
PY
  then fail=1; fi

  # k-mer count, coroutine vs sequential (file mode, so the sequential run-stats block on
  # stdout does not confound the diff). This fixture's mostly-present reads form long
  # positive-skip runs, which exercise the count branch's reported position -- a divergence
  # that ecoli100 (this check's usual index) happens not to trigger. k is small because the
  # fixture reference is small.
  KC=10
  "$MV" query -i "$idx" -r "$READS" --kmer-count -k "$KC" --coroutine -s 8 -t 1 -o "$WORK/kccor$m" >/dev/null 2>&1
  "$MV" query -i "$idx" -r "$READS" --kmer-count -k "$KC" --no-prefetch      -o "$WORK/kcseq$m" >/dev/null 2>&1
  if diff -q <(sort "$WORK/kccor$m.kmers.$KC") <(sort "$WORK/kcseq$m.kmers.$KC") >/dev/null 2>&1; then
    echo "PASS: coroutine k-mer count matches sequential in mode $m"
  else
    echo "FAIL: coroutine k-mer count diverges from sequential in mode $m"; fail=1
  fi
done

# Cross-mode sanity: PML values are mode-independent, so the sequential outputs must
# agree across the modes that ran. (Also catches a broken index in any one mode.)
mode_list=("${!SEQOUT[@]}")
if [ "${#mode_list[@]}" -ge 2 ]; then
  base=${mode_list[0]}
  for m in "${mode_list[@]}"; do
    [ "$m" = "$base" ] && continue
    if diff -q <(sort "${SEQOUT[$base]}") <(sort "${SEQOUT[$m]}") >/dev/null; then
      echo "[pml-modes] cross-check: sequential mode $base == mode $m  (PASS)"
    else
      echo "FAIL: sequential PML differs between mode $base and mode $m"; fail=1
    fi
  done
fi

[ "$ran" -ge 1 ] || { echo "ERROR: no threshold-mode binaries were available to test"; exit 2; }
if [ "$fail" -ne 0 ]; then echo "[pml-modes] FAILED"; exit 1; fi
echo "[pml-modes] all coroutine-vs-sequential PML checks passed across $ran threshold mode(s)"
