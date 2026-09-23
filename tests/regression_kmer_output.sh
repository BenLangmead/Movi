#!/bin/bash
#
# Author: Ben Langmead (ben.langmead@gmail.com)
# Date: Sep 22, 2026
#
# Distributed under the GPL3 license.
# See accompanying LICENSE or https://opensource.org/license/gpl-3-0
#
# Regression test for the per-read records `movi query --kmer` writes, covering two
# properties that every execution path has to share.
#
# 1. Reads shorter than k. A read of length L holds max(L - k + 1, 0) length-k windows,
#    so a read shorter than k holds none and its record reads "0/0". The denominator was
#    once computed as L - k + 1 on unsigned values, which wrapped for L < k and printed a
#    number near 2^64 (a read three bases short of k printed 18446744073709551614). The
#    reads here bracket k at L = k-3, k-1, k and k+1, and include an all-N read, so both
#    the wrap and any off-by-one at L = k are visible.
#
# 2. Run tokens are maximal. A token is "<start>:<length>", a start position and the
#    number of consecutive present k-mers from there. The presence search walks a read
#    from its end and restarts wherever its backward search runs out, so one run of
#    present k-mers can be found in several pieces that abut each other. Tokens must
#    still be the maximal runs, both so a token means a property of the read rather than
#    of the search, and so the file matches the MEM-derived view `--mem --kmer-out <k>`
#    writes, which docs/mem_output.md promises is byte-identical. The reference here is
#    built as rec0 = P+Q and rec1 = Q+R with a read P+Q+R, which is exactly the shape
#    that makes the search restart in the middle of a run of present k-mers.
#
# Both are checked for --kmer and --kmer-count, over every path that serves k-mer
# queries: sequential at -t 1 and -t 4, --coroutine at -t 1 and -t 4, and a strand
# scheduler request. Paths are compared as sorted line multisets, since a multi-threaded
# path may emit in batch order.
#
# Usage:
#   regression_kmer_output.sh [BUILD_DIR]
# Env overrides:
#   MOVI_KMER_KS="8 21"   the k values to test
#
# Requires: the 'movi' launcher and the threshold build pipeline (pfp) to build the toy
# index, movi-regular-thresholds, and python3.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
BUILD="${1:-$ROOT/build-release}"
MV="$BUILD/bin/movi-regular-thresholds"
[ -x "$MV" ] || { echo "ERROR: missing $MV (build first)"; exit 2; }
LAUNCH=""
for c in "$BUILD/movi" "$BUILD/bin/movi"; do [ -x "$c" ] && LAUNCH="$c" && break; done
[ -n "$LAUNCH" ] || { echo "ERROR: 'movi' launcher not found under $BUILD"; exit 2; }

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

KS="${MOVI_KMER_KS:-8 21}"
# The ftab depth used for the index and for the MEM run that produces the --kmer-out
# view. A MEM search reports membership for every k at or above --min-mem-length, so it
# has to stay at or below the smallest k tested.
FTAB_K=4

echo "[regression] [kmer-output] building the toy reference and reads ..."
python3 - "$WORK" "$KS" <<'PY'
import os, random, sys
work, ks = sys.argv[1], [int(x) for x in sys.argv[2].split()]
rng = random.Random(20260922)
def seq(n): return "".join(rng.choice("ACGT") for _ in range(n))

# rec0 and rec1 share the middle block Q, so the read P+Q+R below is present only in
# pieces: its suffix Q+R sits in rec1 and its prefix P+Q sits in rec0, which is what
# makes the presence search restart inside a run of present k-mers.
P, Q, R = seq(60), seq(60), seq(60)
filler = seq(600)
with open(os.path.join(work, "ref.fa"), "w") as o:
    o.write(">rec0\n%s\n>rec1\n%s\n>rec2\n%s\n" % (P + Q, Q + R, filler))

reads = [("spanning", P + Q + R)]
for k in ks:
    # Bracket k. k-3 and k-1 are short (no window at all); k gives exactly one window.
    for label, n in (("short3", k - 3), ("short1", k - 1), ("exact", k), ("plus1", k + 1)):
        st = rng.randrange(0, len(filler) - n)
        reads.append(("k%d_%s" % (k, label), filler[st:st + n]))
    reads.append(("k%d_allN" % k, "N" * (k + 4)))
    # A read that is not in the reference at all, so found is 0 with a nonzero total.
    reads.append(("k%d_absent" % k, seq(k + 10)))
with open(os.path.join(work, "reads.fa"), "w") as o:
    for name, s in reads:
        o.write(">%s\n%s\n" % (name, s))
PY

echo "[regression] [kmer-output] building the index and ftab-$FTAB_K ..."
"$LAUNCH" build --type regular-thresholds --index "$WORK/idx" --fasta "$WORK/ref.fa" --separators >/dev/null 2>&1 \
  || { echo "ERROR: index build failed (threshold pipeline available?)"; exit 2; }
"$LAUNCH" ftab --index "$WORK/idx" --ftab-k "$FTAB_K" >/dev/null 2>&1 \
  || { echo "ERROR: ftab-$FTAB_K build failed"; exit 2; }

status=0
for K in $KS; do
  for MODE in kmer kmer-count; do
    echo "[regression] [kmer-output] k=$K --$MODE over every path ..."
    files=""
    for path in "seq1:-t 1" "seq4:-t 4" "co1:--coroutine -t 1" "co4:--coroutine -t 4" "str:-s 8 -t 1"; do
      name="${path%%:*}"; flags="${path#*:}"
      # shellcheck disable=SC2086
      "$MV" query --index "$WORK/idx" --read "$WORK/reads.fa" "--$MODE" -k "$K" \
            $flags -o "$WORK/$MODE.$K.$name" >/dev/null 2>&1 \
        || { echo "  FAIL: --$MODE k=$K path $name exited nonzero"; status=1; continue; }
      files="$files $name=$WORK/$MODE.$K.$name.kmers.$K"
    done

    # The MEM-derived view of the same reads, for the byte comparison in check 2.
    ko=""
    if "$MV" query --index "$WORK/idx" --read "$WORK/reads.fa" --mem --ftab-k "$FTAB_K" \
          --min-mem-length "$FTAB_K" --kmer-out "$K" -t 1 -o "$WORK/ko.$K" >/dev/null 2>&1; then
      ko="$WORK/ko.$K.kmers.$K"
    else
      echo "  FAIL: --mem --kmer-out $K exited nonzero"; status=1
    fi

    # shellcheck disable=SC2086
    python3 "$HERE/check_kmer_records.py" --reads "$WORK/reads.fa" --k "$K" --mode "$MODE" \
        ${ko:+--kmer-out "$ko"} $files || status=1
  done
done

[ $status -eq 0 ] && echo "[regression] all k-mer output checks passed" \
                  || echo "[regression] k-mer output checks FAILED"
exit $status
