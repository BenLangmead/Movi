#!/usr/bin/env python3
#
# Author: Ben Langmead (ben.langmead@gmail.com)
# Date: Sep 22, 2026
#
# Distributed under the GPL3 license.
# See accompanying LICENSE or https://opensource.org/license/gpl-3-0
#
"""Brute-force oracle for `movi query --mem`.

Definitions used here (they match what the Movi MEM search is meant to report):

* The searchable text is every reference record and its reverse complement, taken
  one record at a time.  A match never spans two records or a record and its reverse
  complement, which is what an index built with --separators provides.  Reference
  characters other than ACGT are replaced by A, as movi-prepare-ref does.  With
  --junctions the records are instead concatenated in movi-prepare-ref order without
  separators (R1 rc(R1) R2 rc(R2) ...), which models an index built without
  --separators; that mode exists to diagnose index-configuration effects.
* For a query Q, rmax(s) is the largest e such that Q[s:e] occurs in the text (e == s
  when Q[s] is not A, C, G or T, or does not occur).
* A MEM (in the query-centric sense, also called an SMEM) is a half-open interval
  [s, e) with e = rmax(s) > s and either s == 0 or rmax(s - 1) < e.  Such an interval
  occurs in the text and cannot be extended left or right, so no MEM contains another.
  MEMs of length >= L are reported, in increasing order of s.  L <= 1 means all MEMs.
* The count of a MEM is the number of occurrences of Q[s:e] in the text (overlapping
  occurrences counted separately, both strands of each record included).

Output mirrors the default `movi query --mem --stdout` format:
    <id>\t<covered>/<len>\t<s>:<e>:<count> <s>:<e>:<count> ... \n
where <covered> is the number of query positions inside at least one reported MEM.

With --ms the oracle instead prints the matching statistics of each query, one line
per query: the id and then, for every position i, the length of the longest suffix of
Q[0:i+1] that occurs in the text.
"""
import argparse
import sys

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0] if len(line) > 1 else "", []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def clean_ref(seq):
    seq = seq.upper()
    return "".join(c if c in COMP else "A" for c in seq)


def revcomp(seq):
    return "".join(COMP[c] for c in reversed(seq))


class Texts:
    """The searchable strands, with a k-mer index so substring queries stay fast on
    references of tens of kilobases.  Strands are joined with '#', which never occurs
    in a query, so no match spans two strands."""

    K = 12

    def __init__(self, strands):
        self.joined = "#".join(strands)
        self.short = set()
        self.index = {}
        t, K = self.joined, self.K
        for i in range(len(t)):
            for j in range(i + 1, min(i + K, len(t)) + 1):
                if t[j - 1] == "#":
                    break
                self.short.add(t[i:j])
            if i + K <= len(t) and "#" not in t[i:i + K]:
                self.index.setdefault(t[i:i + K], []).append(i)

    def positions(self, pat):
        if len(pat) < self.K:
            return None
        return [p for p in self.index.get(pat[:self.K], ()) if self.joined.startswith(pat, p)]

    def occurs(self, pat):
        if len(pat) <= self.K:
            return pat in self.short
        return bool(self.positions(pat))

    def count(self, pat):
        pos = self.positions(pat)
        if pos is not None:
            return len(pos)
        t, total = self.joined, 0
        start = t.find(pat)
        while start != -1:
            total += 1
            start = t.find(pat, start + 1)
        return total


def build_texts(ref_path, junctions):
    strands = []
    for _, seq in read_fasta(ref_path):
        seq = clean_ref(seq)
        strands.append(seq)
        strands.append(revcomp(seq))
    if junctions:
        return Texts(["".join(strands)])
    return Texts(strands)


def occurs(texts, pat):
    return texts.occurs(pat)


def count_occ(texts, pat):
    return texts.count(pat)


def rmax_array(texts, q):
    m = len(q)
    rmax = [0] * m
    e = 0
    for s in range(m):
        # rmax is non-decreasing in s, so the scan resumes from the previous end.
        if e < s:
            e = s
        while e < m and q[e] in COMP and occurs(texts, q[s:e + 1]):
            e += 1
        rmax[s] = e
    return rmax


def mems(texts, q, min_len):
    rmax = rmax_array(texts, q)
    out = []
    for s in range(len(q)):
        e = rmax[s]
        if e <= s:
            continue
        if s > 0 and rmax[s - 1] >= e:
            continue
        if e - s < max(min_len, 1):
            continue
        out.append((s, e, count_occ(texts, q[s:e])))
    return out


def matching_statistics(texts, q):
    ms, cur = [], 0
    for i in range(len(q)):
        if q[i] not in COMP:
            cur = 0
        else:
            cur += 1
            while cur > 0 and not occurs(texts, q[i - cur + 1:i + 1]):
                cur -= 1
        ms.append(cur)
    return ms


def format_mems(qid, q, ms):
    covered, run_s, run_e, in_run = 0, 0, 0, False
    for s, e, _ in ms:
        if not in_run:
            run_s, run_e, in_run = s, e, True
        elif s <= run_e:
            run_e = max(run_e, e)
        else:
            covered += run_e - run_s
            run_s, run_e = s, e
    if in_run:
        covered += run_e - run_s
    body = "".join("%d:%d:%d " % m for m in ms)
    return "%s\t%d/%d\t%s\n" % (qid, covered, len(q), body)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ref")
    ap.add_argument("reads")
    ap.add_argument("-l", "--min-mem-length", type=int, default=25)
    ap.add_argument("--junctions", action="store_true",
                    help="concatenate records without separators, as an index built without --separators does")
    ap.add_argument("--ms", action="store_true", help="print matching statistics instead of MEMs")
    args = ap.parse_args()
    texts = build_texts(args.ref, args.junctions)
    out = sys.stdout
    for qid, q in read_fasta(args.reads):
        if args.ms:
            out.write(qid + "\t" + " ".join(map(str, matching_statistics(texts, q))) + "\n")
        else:
            out.write(format_mems(qid, q, mems(texts, q, args.min_mem_length)))


if __name__ == "__main__":
    main()
