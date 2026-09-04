# The `--mem` output format

`movi query --mem` writes **one line per read**: the read id, how much of the read
matched, then a space-separated token per maximal exact match found in it.

```
NC_096936.1	150/150	0:150:1
m54238_180903_015530/30540296/ccs	11890/13007	0:29:8 14:31:100 15:33:97
m54238_180925_225123/53412816/ccs	0/11402
```

This mirrors the k-mer output line for line:

```
--kmer   NC_096936.1	120/120	0:120
--mem    NC_096936.1	150/150	0:150:1
```

**Fields** are tab separated: the read id (the header up to the first whitespace), then
`covered/read_len`, then the MEM list. `covered` is the number of bases of the read lying
inside at least one MEM, with overlapping matches counted once, so the field reads as
"how much of this read matched" exactly as the k-mer output's `found/total` reads as "how
many of its k-mers are present". The read length is the denominator. The third field is
empty when the search found no MEM, as in the third line above.

**MEM tokens** are `start:end:count`:

- `start` is the 0-based offset of the match in the read.
- `end` is **exclusive**, so the match is `read[start .. end)` and its length is
  `end - start`.
- `count` is how many places the MEM occurs in the collection, on **either strand**,
  because the index holds the reverse complement alongside the forward sequence. On a
  100-haplotype pangenome a MEM present in every haplotype reads `100`, and one specific
  to a single haplotype reads `1`. It is a 64-bit value: a MEM conserved across a large
  bacterial collection occurs far more often than a 16- or 32-bit field would hold.

Lines are as long as a read is rich in matches. On HiFi reads against a 100-haplotype
pangenome at `--min-mem-length 12`, the median read gives about 2 KB and the longest
seen is 18 KB. All-MEM mode (`--min-mem-length 1`) is far denser and will produce
correspondingly longer lines.

## Guarantees

- **Every read gets a record**, including reads where the search found no MEM. Those
  produce a header and nothing else, so a consumer sees exactly the reads it submitted.
- **Reads appear in input order.** With `--coroutine` a reorder buffer restores that
  order, so the output does not depend on which worker finished first.
- **Within a read, MEMs have strictly increasing `start`.** A maximal match beginning at
  a given offset is unique, so no two lines of one record share a `start`.
- **One line per read**, so `wc -l` counts reads and the file joins line for line
  against the k-mer views or any other per-read output.

## Where this still differs from the k-mer output

Two differences are deliberate, and both matter when reading a `kmers.<k>` file beside
the `.mems` file it came from:

- **Token order.** MEM tokens ascend by `start`. The k-mer output's run tokens *descend*,
  because that query scans each read from its end.
- **What the token holds.** A MEM token is `start:end:count`, an interval plus an
  occurrence multiplicity. A k-mer run token is `start:length`, a position plus how many
  consecutive k-mers are present. The last field of each therefore means different
  things: a multiplicity in one, a length in the other.

The k-mer format is left as released rather than changed to match.

## Which MEMs are reported

`--min-mem-length L` reports the MEMs of length `L` or more. It acts purely as a filter
on the output: running at a smaller `L` and discarding the shorter matches afterwards
gives byte-identical results for the longer ones. Smaller is also much faster, because
the search cost is set by how far it must extend past the ftab seed. Running at
`L = ftab-k` (the deepest ftab in the index) skips that extension entirely, and is the
setting to prefer unless there is a reason not to.

## Deriving k-mer membership

Because a MEM `[start, end)` with `end - start >= k` means every k-mer starting in
`[start, end - k]` is present, the record set answers membership for **every k at or
above `--min-mem-length`**. Two ways to get it:

- `movi query --mem --kmer-out 31,63 ...` writes the membership views during the query,
  one file per k, byte-identical to what `movi query --kmer -k <k>` would write.
- `movi kmers-from-mems --mems out.mems --kmer-out 31,63` rebuilds the same views from a
  saved MEM stream, with no index. Pass `--mems-of` as well only if you want the
  aggregate report's invalid-window count, which needs the bases: whether a window holds
  a non-ACGT base depends on k, so it cannot be recorded in the MEM stream.

## Compatibility

Movi 2.0.0 wrote one line per MEM, `read_id start end count`, with no read length and no
line at all for reads without a MEM. `--legacy-mems` still writes exactly that, byte for
byte, for consumers built against it.

That format is write-only. It cannot say how long a read was, and it cannot represent a
read with no MEMs, so `kmers-from-mems` refuses it with an explicit message rather than
misreading it. Anything that needs to be read back should use the default format.
