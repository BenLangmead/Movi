#!/usr/bin/env python3
#
# Author: Ben Langmead (ben.langmead@gmail.com)
# Date: Sep 22, 2026
#
# Distributed under the GPL3 license.
# See accompanying LICENSE or https://opensource.org/license/gpl-3-0
#
"""Checker for regression_kmer_output.sh: validates the records `movi query --kmer` writes.

Given the reads and one k-mer output file per execution path, it requires that

  * every read gets exactly one record, and no record names a read that was not asked
    for;
  * the denominator of the found/total field is the number of length-k windows the read
    holds, max(len - k + 1, 0), so a read shorter than k reads "0/0" rather than a value
    that wrapped around zero;
  * found never exceeds total, and the tokens account for exactly `found` k-mers: in
    presence mode a token is a run and the run lengths sum to found, while in count mode
    a token is one k-mer whose value is its occurrence multiplicity, so the tokens are
    counted rather than summed;
  * in presence mode the run tokens descend by start and no two of them abut, that is,
    each run is maximal;
  * every path writes the same records (compared as a sorted line multiset, since a
    multi-threaded path may emit in batch order);
  * with --kmer-out, the presence records equal the MEM-derived membership view byte for
    byte, which is what docs/mem_output.md promises.

Exit status is 0 when every check passes, 1 otherwise.
"""
import argparse
import sys


def read_fasta(path):
    """Yields (id, sequence) pairs, taking the id as the header up to the first space."""
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def parse_records(path):
    """Reads a k-mer output file into {id: (found, total, [(start, value), ...])}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError("malformed record: %r" % line[:120])
            found, total = (int(x) for x in parts[1].split("/"))
            toks = []
            for tok in (parts[2].split() if len(parts) > 2 else []):
                fields = tok.split(":")
                toks.append((int(fields[0]), int(fields[1])))
            out[parts[0]] = (found, total, toks)
    return out


def check_one(name, recs, lengths, k, mode, fails):
    missing = sorted(set(lengths) - set(recs))
    extra = sorted(set(recs) - set(lengths))
    if missing:
        fails.append("%s: %d read(s) have no record, e.g. %s" % (name, len(missing), missing[0]))
    if extra:
        fails.append("%s: %d record(s) name an unknown read, e.g. %s" % (name, len(extra), extra[0]))

    for rid in sorted(set(recs) & set(lengths)):
        found, total, toks = recs[rid]
        want_total = max(lengths[rid] - k + 1, 0)
        if total != want_total:
            fails.append("%s: read %s of length %d at k=%d has total %d, expected %d"
                         % (name, rid, lengths[rid], k, total, want_total))
        if found > total:
            fails.append("%s: read %s reports found %d > total %d" % (name, rid, found, total))
        if mode == "kmer":
            if sum(length for _, length in toks) != found:
                fails.append("%s: read %s run lengths sum to %d, but found is %d"
                             % (name, rid, sum(length for _, length in toks), found))
            starts = [start for start, _ in toks]
            if starts != sorted(starts, reverse=True):
                fails.append("%s: read %s run tokens do not descend by start: %s" % (name, rid, starts))
            # Tokens descend, so the token after one starting at s covers positions
            # ending just below s; it abuts when its start plus its length reaches s.
            for (s1, l1), (s2, l2) in zip(toks, toks[1:]):
                if s2 + l2 == s1:
                    fails.append("%s: read %s has abutting run tokens %d:%d and %d:%d, "
                                 "which are one maximal run" % (name, rid, s1, l1, s2, l2))
                elif s2 + l2 > s1:
                    fails.append("%s: read %s has overlapping run tokens %d:%d and %d:%d"
                                 % (name, rid, s1, l1, s2, l2))
        else:
            if len(toks) != found:
                fails.append("%s: read %s has %d count tokens but found is %d"
                             % (name, rid, len(toks), found))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reads", required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--mode", choices=("kmer", "kmer-count"), required=True)
    ap.add_argument("--kmer-out", help="the --mem --kmer-out view of the same reads")
    ap.add_argument("paths", nargs="+", metavar="NAME=FILE")
    args = ap.parse_args()

    lengths = {rid: len(seq) for rid, seq in read_fasta(args.reads)}
    if not lengths:
        print("  FAIL: no reads in %s" % args.reads)
        return 1

    fails = []
    baseline_name, baseline_lines = None, None
    for spec in args.paths:
        name, _, path = spec.partition("=")
        try:
            recs = parse_records(path)
        except (OSError, ValueError) as exc:
            fails.append("%s: could not read %s (%s)" % (name, path, exc))
            continue
        check_one(name, recs, lengths, args.k, args.mode, fails)
        lines = sorted(open(path).read().splitlines())
        if baseline_lines is None:
            baseline_name, baseline_lines = name, lines
        elif lines != baseline_lines:
            differ = set(lines) ^ set(baseline_lines)
            fails.append("%s differs from %s on %d line(s), e.g. %r"
                         % (name, baseline_name, len(differ), sorted(differ)[0][:120]))

    if args.kmer_out and args.mode == "kmer":
        # The MEM-derived view and a direct presence query are documented as byte
        # identical, so they are compared as text rather than as a set of positions.
        want = open(args.kmer_out).read().splitlines()
        got = open(args.paths[0].partition("=")[2]).read().splitlines()
        if want != got:
            differ = [(a, b) for a, b in zip(got, want) if a != b]
            fails.append("--kmer and --mem --kmer-out %d disagree on %d line(s)%s"
                         % (args.k, len(differ) + abs(len(want) - len(got)),
                            (", e.g. %r vs %r" % differ[0]) if differ else ""))

    for f in fails[:15]:
        print("  FAIL: " + f)
    if fails:
        print("  %d check(s) failed for k=%d --%s" % (len(fails), args.k, args.mode))
        return 1
    print("  PASS: k=%d --%s, %d reads, %d path(s) agree%s"
          % (args.k, args.mode, len(lengths), len(args.paths),
             ", and match --kmer-out" if (args.kmer_out and args.mode == "kmer") else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
