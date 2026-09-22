#!/usr/bin/env python3
#
# Author: Ben Langmead (ben.langmead@gmail.com)
# Date: Sep 22, 2026
#
# Distributed under the GPL3 license.
# See accompanying LICENSE or https://opensource.org/license/gpl-3-0
#
"""Randomized differential test of `movi query --mem` against tests/mem_oracle.py.

For each scenario (a random reference with one or more records, built with or without
--separators) the script builds a Movi index and a set of ftabs, generates random
queries and queries derived from the reference with substitutions and Ns, and runs
the MEM query for every combination of minimum length L and ftab depth that Movi
accepts, over every execution path:

  seq1  sequential, -t 1          co1  --coroutine, -t 1
  seq4  sequential, -t 4          co4  --coroutine, -t 4
  str   strand scheduler requested (-s 8, prefetching left on; MEM turns it off)

Checks, per (scenario, L, ftab-k):
  * the per-read records of all paths are identical (sorted-line comparison, since
    multi-threaded paths may emit in batch order; seq1 and co1 are also compared
    byte for byte, both emitting in input order);
  * the seq1 records equal the oracle's.  On an index built with --separators the
    oracle respects record boundaries.  On an index built without them the oracle is
    run in --junctions mode, which models the concatenated text such an index holds.
  * with --require-refusal, every MEM query on an index built without --separators
    must instead fail with Movi's message asking for --separators.

For indexes built without --separators the script also reports how many reads get a
different answer from the junction oracle than from the record-respecting oracle,
that is, reads whose MEMs such an index gets wrong even when Movi computes exactly
what its text holds.

With --kmer-out-check, the --mem --kmer-out membership views are also compared with a
direct --kmer presence query on the same reads.

Exit status is 0 when every check passes, 1 otherwise.
"""
import argparse
import os
import random
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mem_oracle  # noqa: E402

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
PATHS = {
    "seq1": ["-t", "1"],
    "seq4": ["-t", "4"],
    "co1": ["--coroutine", "-t", "1"],
    "co4": ["--coroutine", "-t", "4"],
    "str": ["-s", "8", "-t", "1"],
}


def rand_seq(rng, n, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(n))


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def mutate(rng, s, nsub, nN):
    s = list(s)
    for _ in range(nsub):
        s[rng.randrange(len(s))] = rng.choice("ACGT")
    for _ in range(nN):
        s[rng.randrange(len(s))] = "N"
    return "".join(s)


def make_reference(rng, nrec, lo, hi, repeat_frac):
    recs = []
    for i in range(nrec):
        n = rng.randint(lo, hi)
        seq = rand_seq(rng, n)
        # Plant copies of earlier segments so matches have multiplicity > 1 and the
        # BWT has longer runs.
        pool = "".join(recs) + seq
        for _ in range(int(repeat_frac * n / 20)):
            a = rng.randrange(0, max(1, len(pool) - 20))
            seg = pool[a:a + rng.randint(5, 20)]
            if rng.random() < 0.5:
                seg = revcomp(seg)
            b = rng.randrange(0, max(1, n - len(seg)))
            seq = seq[:b] + seg + seq[b + len(seg):]
        recs.append(seq[:n])
    return recs


def make_queries(rng, recs, nq):
    qs = []
    full = recs + [revcomp(r) for r in recs]
    for i in range(nq):
        kind = rng.random()
        n = rng.randint(4, 80)
        if kind < 0.15:
            q = rand_seq(rng, n)
        elif kind < 0.30:
            # Straddle a record boundary (the junction of the concatenated text).
            a, b = rng.choice(full), rng.choice(full)
            h = rng.randint(1, min(n - 1, len(a), len(b)))
            q = a[len(a) - h:] + b[:n - h]
            q = mutate(rng, q, rng.randint(0, 2), 0)
        else:
            src = rng.choice(full)
            st = rng.randrange(0, max(1, len(src) - n))
            q = src[st:st + n]
            q = mutate(rng, q, rng.randint(0, max(1, len(q) // 8)),
                       1 if rng.random() < 0.15 else 0)
        if not q:
            q = rand_seq(rng, 5)
        qs.append(("q%d" % i, q))
    return qs


def write_fasta(path, recs):
    with open(path, "w") as fh:
        for name, s in recs:
            fh.write(">%s\n%s\n" % (name, s))


def run(cmd, timeout=20):
    # A runaway query (for example one that stops advancing along the read) would
    # otherwise hang the test, so every command gets a time limit.
    try:
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                              timeout=timeout)
    except subprocess.TimeoutExpired:
        return subprocess.CompletedProcess(cmd, 124, "", "timed out after %ds" % timeout)


def build_index(launcher, ref, idx, separators, ftab_ks):
    cmd = [launcher, "build", "--index", idx, "--fasta", ref]
    if separators:
        cmd.append("--separators")
    r = run(cmd)
    if r.returncode != 0 or not os.path.exists(os.path.join(idx, "index.movi")):
        raise RuntimeError("index build failed for %s:\n%s" % (ref, r.stderr[-2000:]))
    for k in ftab_ks:
        r = run([launcher, "ftab", "--index", idx, "--ftab-k", str(k)])
        if r.returncode != 0 or not os.path.exists(os.path.join(idx, "ftab.%d.bin" % k)):
            raise RuntimeError("ftab-%d build failed:\n%s" % (k, r.stderr[-2000:]))


def movi_mem(binary, idx, reads, L, ftab_k, path_flags, out_prefix, extra=()):
    cmd = [binary, "query", "--index", idx, "--read", reads, "--mem",
           "--min-mem-length", str(L), "--ftab-k", str(ftab_k), "-o", out_prefix]
    cmd += list(path_flags) + list(extra)
    r = run(cmd)
    if r.returncode != 0:
        return None, r.stderr[-1500:]
    with open(out_prefix + ".mems") as fh:
        return fh.read(), None


def lines_by_id(text):
    d = {}
    for line in text.splitlines():
        d[line.split("\t", 1)[0]] = line
    return d


def classify(movi_line, oracle_line):
    """Names the kinds of disagreement between one read's Movi and oracle records."""
    def parse(line):
        parts = line.split("\t")
        return [tuple(map(int, t.split(":"))) for t in parts[2].split()] if len(parts) > 2 else []
    m, o = parse(movi_line), parse(oracle_line)
    ospans = {(s, e) for s, e, _ in o}
    mspans = {(s, e) for s, e, _ in m}
    kinds = set()
    for s, e in mspans - ospans:
        if any(os_ <= s and e <= oe for os_, oe in ospans):
            kinds.add("under")  # a reported span strictly inside a true MEM
        elif any(s <= os_ and oe <= e for os_, oe in ospans):
            kinds.add("over")   # a reported span strictly containing a true MEM
        else:
            kinds.add("other")
    if ospans - mspans and not kinds:
        kinds.add("missing")
    if len([1 for s, e, _ in m]) != len(mspans):
        kinds.add("dup")
    starts = [s for s, e, _ in m]
    if len(starts) != len(set(starts)):
        kinds.add("same-start")
    if mspans == ospans and {x for x in m} != {x for x in o}:
        kinds.add("count")
    return kinds


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--build-dir", default=os.path.join(HERE, "..", "build-release"))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scenarios", type=int, default=6,
                    help="random references per separator setting")
    ap.add_argument("--queries", type=int, default=60)
    ap.add_argument("--lengths", default="1,2,3,4,6,10")
    ap.add_argument("--ftab-ks", default="2,3,4,6")
    ap.add_argument("--paths", default=",".join(PATHS))
    ap.add_argument("--separators", default="yes,no", help="which index kinds to test")
    ap.add_argument("--kmer-out-check", action="store_true")
    ap.add_argument("--require-refusal", action="store_true",
                    help="MEM queries on indexes built without --separators must be refused")
    ap.add_argument("--keep", action="store_true", help="keep the work directory")
    ap.add_argument("--work")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    build = os.path.abspath(args.build_dir)
    launcher = os.path.join(build, "movi")
    binary = os.path.join(build, "bin", "movi-regular-thresholds")
    for p in (launcher, binary):
        if not os.access(p, os.X_OK):
            print("ERROR: missing %s" % p)
            return 2
    lengths = [int(x) for x in args.lengths.split(",")]
    ftab_ks = [int(x) for x in args.ftab_ks.split(",")]
    paths = args.paths.split(",")
    work = args.work or tempfile.mkdtemp(prefix="movi_mem_diff_")
    os.makedirs(work, exist_ok=True)
    rng = random.Random(args.seed)

    # agree[(sep, L)] = [reads agreeing with oracle, reads compared]
    agree = defaultdict(lambda: [0, 0])
    kinds_seen = defaultdict(lambda: defaultdict(int))
    failures = []
    path_mismatch = 0
    configs = 0
    examples = []
    refused = 0
    # junction_effect[nrec] = [reads whose junction and record oracles differ, reads]
    junction_effect = defaultdict(lambda: [0, 0])

    sep_kinds = [s == "yes" for s in args.separators.split(",")]
    for sep in sep_kinds:
        for sc in range(args.scenarios):
            nrec = [1, 1, 2, 3, 4, 5][sc % 6]
            recs = make_reference(rng, nrec, 60, 400, 0.6)
            tag = "%s_%d" % ("sep" if sep else "nosep", sc)
            sdir = os.path.join(work, tag)
            os.makedirs(sdir, exist_ok=True)
            ref = os.path.join(sdir, "ref.fa")
            write_fasta(ref, [("rec%d" % i, s) for i, s in enumerate(recs)])
            reads = os.path.join(sdir, "reads.fa")
            queries = make_queries(rng, recs, args.queries)
            write_fasta(reads, queries)
            idx = os.path.join(sdir, "idx")
            build_index(launcher, ref, idx, sep, ftab_ks)
            texts = mem_oracle.build_texts(ref, junctions=not sep)
            rec_texts = mem_oracle.build_texts(ref, junctions=False)
            scenario_refused = False
            for L in lengths:
                oracle = {qid: mem_oracle.format_mems(qid, q, mem_oracle.mems(texts, q, L)).rstrip("\n")
                          for qid, q in queries}
                if not sep:
                    for qid, q in queries:
                        rec = mem_oracle.format_mems(qid, q, mem_oracle.mems(rec_texts, q, L)).rstrip("\n")
                        junction_effect[nrec][0] += rec != oracle[qid]
                        junction_effect[nrec][1] += 1
                for fk in ftab_ks:
                    if L > 1 and fk > L:
                        continue  # rejected by movi for length-thresholded search
                    configs += 1
                    outs = {}
                    for p in paths:
                        out, err = movi_mem(binary, idx, reads, L, fk, PATHS[p],
                                            os.path.join(sdir, "o_%d_%d_%s" % (L, fk, p)))
                        if not sep and args.require_refusal:
                            if out is None and "--separators" in err:
                                refused += 1
                                scenario_refused = True
                            else:
                                failures.append("%s L=%d k=%d %s: not refused on an index without separators"
                                                % (tag, L, fk, p))
                            out = None
                        elif out is None:
                            failures.append("%s L=%d k=%d %s: movi failed: %s" % (tag, L, fk, p, err.strip()[-300:]))
                        outs[p] = out
                    ok = [p for p in paths if outs[p] is not None]
                    if not ok:
                        continue
                    base = sorted(outs[ok[0]].splitlines())
                    for p in ok[1:]:
                        if sorted(outs[p].splitlines()) != base:
                            path_mismatch += 1
                            failures.append("%s L=%d k=%d: path %s differs from %s" % (tag, L, fk, p, ok[0]))
                    if "seq1" in ok and "co1" in ok and outs["seq1"] != outs["co1"]:
                        if sorted(outs["seq1"].splitlines()) == sorted(outs["co1"].splitlines()):
                            failures.append("%s L=%d k=%d: seq1 and co1 differ in order only" % (tag, L, fk))
                    ref_path = "seq1" if "seq1" in ok else ok[0]
                    got = lines_by_id(outs[ref_path])
                    for qid, q in queries:
                        agree[(sep, L)][1] += 1
                        if got.get(qid) == oracle[qid]:
                            agree[(sep, L)][0] += 1
                        else:
                            for kd in classify(got.get(qid, qid + "\t0/0\t"), oracle[qid]):
                                kinds_seen[(sep, L)][kd] += 1
                            if len(examples) < 12:
                                examples.append((tag, L, fk, qid, q, got.get(qid), oracle[qid]))
            if args.kmer_out_check and not scenario_refused:
                failures += kmer_out_check(binary, idx, reads, sdir, tag, ftab_ks)

    print("configs run: %d, path mismatches: %d, refused (no separators): %d"
          % (configs, path_mismatch, refused))
    for nrec, (d, n) in sorted(junction_effect.items()):
        print("  without separators, %d record(s): junction oracle differs from record oracle on %d / %d reads"
              % (nrec, d, n))
    print("oracle agreement (reads matching / reads compared, summed over ftab-k):")
    for (sep, L), (a, n) in sorted(agree.items()):
        kinds = ", ".join("%s=%d" % kv for kv in sorted(kinds_seen[(sep, L)].items()))
        print("  %-12s L=%-3d %6d / %-6d %s" % ("separators" if sep else "no-separators", L, a, n,
                                               ("[" + kinds + "]") if kinds else ""))
    for ex in examples[: (12 if args.verbose else 4)]:
        tag, L, fk, qid, q, got, want = ex
        print("  example %s L=%d ftab-k=%d %s %s\n    movi:   %s\n    oracle: %s" % (tag, L, fk, qid, q, got, want))
    total_bad = sum(n - a for a, n in agree.values())
    for f in failures[:20]:
        print("  FAIL: " + f)
    if not args.keep and not args.work:
        shutil.rmtree(work, ignore_errors=True)
    else:
        print("work directory: " + work)
    if total_bad or failures:
        print("FAIL: %d read records disagree with the oracle, %d other failures" % (total_bad, len(failures)))
        return 1
    print("PASS")
    return 0


def membership(path, k):
    """Reads a k-mer membership file into {id: (found, total, present k-mer starts)}.

    Each line is "<id>\t<found>/<total>\t<start>:<run> ...". The two producers may split
    one run of consecutive present k-mers differently, so runs are expanded to the set
    of present starts. A read shorter than k has no k-mers; its total is taken as 0
    whatever the file says.
    """
    out = {}
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        found, total = map(int, parts[1].split("/"))
        if total > 1 << 40:
            total = 0
        starts = set()
        for tok in (parts[2].split() if len(parts) > 2 else []):
            st, run = map(int, tok.split(":"))
            starts.update(range(st, st + run))
        out[parts[0]] = (found, total, frozenset(starts))
    return out


def kmer_out_check(binary, idx, reads, sdir, tag, ftab_ks):
    """Compares --mem --kmer-out membership views with direct --kmer presence queries."""
    fails = []
    fk = max(ftab_ks)
    for k in (fk, fk + 3, 12):
        if k < fk:
            continue
        pre = os.path.join(sdir, "ko_%d" % k)
        r = run([binary, "query", "--index", idx, "--read", reads, "--mem", "--ftab-k", str(fk),
                 "--min-mem-length", str(fk), "--kmer-out", str(k), "-t", "1", "-o", pre])
        if r.returncode != 0:
            fails.append("%s kmer-out k=%d failed: %s" % (tag, k, r.stderr.strip()[-300:]))
            continue
        pre2 = os.path.join(sdir, "kd_%d" % k)
        r2 = run([binary, "query", "--index", idx, "--read", reads, "--kmer", "-k", str(k),
                  "-t", "1", "-o", pre2])
        if r2.returncode != 0:
            fails.append("%s kmer k=%d failed: %s" % (tag, k, r2.stderr.strip()[-300:]))
            continue
        a_path, b_path = pre + ".kmers.%d" % k, pre2 + ".kmers.%d" % k
        if not (os.path.exists(a_path) and os.path.exists(b_path)):
            fails.append("%s kmer-out k=%d: missing %s or %s" % (tag, k, a_path, b_path))
            continue
        a, b = membership(a_path, k), membership(b_path, k)
        bad = [q for q in set(a) | set(b) if a.get(q) != b.get(q)]
        if bad:
            q = sorted(bad)[0]
            fails.append("%s kmer-out k=%d disagrees with --kmer on %d reads, e.g. %s: %r vs %r"
                         % (tag, k, len(bad), q, a.get(q), b.get(q)))
    return fails


if __name__ == "__main__":
    sys.exit(main())
