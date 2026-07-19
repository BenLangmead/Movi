# Query implementations and latency hiding

Movi runs each hot query in up to three styles. This doc is the index: it names
the styles, says where each lives, and maps every query to its implementation(s).

## The three styles

- **(a) Sequential** — no latency hiding; straight-line LF walks. The reference
  implementation that the other two must match exactly.
- **(b) Manual "strands"** — prefetch the next move-structure row, then hand-switch
  among many in-flight reads ("strands") at each prefetch point, so other reads
  make progress while a cache line fills. A round-robin state machine.
- **(c) Coroutines** — prefetch, then `co_yield` at each prefetch point; a scheduler
  resumes another coroutine while the line fills. Same idea as (b), expressed with
  C++20 coroutines instead of a hand-written state machine.

Only long-range, non-sequential move-structure accesses (LF destinations, ftab
entries) are worth hiding; sequential run scans are left to the hardware prefetcher.

## Where each style lives

| Style | Files |
|-------|-------|
| (a) sequential | `query_pml.cpp` (`query_pml`), `query_zml.cpp` (`query_zml`), `query_kmer.cpp` (`query_all_kmers`, `query_kmers_from`), `query_kmer_bv.cpp` (`query_kmers_count_bv`, `query_kmers_id_bv`), `query_mem.cpp` (`query_mems`, `query_mem_bml`); shared repositioning primitives in `move_structure_query.cpp` (`reposition_thresholds`, `reposition_randomly`, `handle_reposition_*`, `reposition_up/down`), search primitives in `move_structure_search.cpp` (`backward_search`, `extend_bidirectional`, `*_search_step`), and `move_structure.cpp` (`LF_move`, `fast_forward`) |
| (b) manual strands | `read_processor.cpp` (`process_latency_hiding`); `struct Strand` in `read_processor.hpp`. Serves PML, ZML and whole-read exact-count only |
| (c) coroutines | `coroutine_processor.cpp` (the `query_*_coroutine` bodies + scheduler + `run_coroutine_query`; guarded to the threshold index modes, `(MODE==6 \|\| MODE==7 \|\| MODE==8) && COLOR_MODE==0`, i.e. regular-, sampled-, and blocked-thresholds). The bodies use the mode-portable move-structure primitives (`get_c`/`get_id`/`reposition_thresholds`/`get_n`/`get_offset`) throughout. Reached via **`movi query --coroutine`** (dispatch in `movi.cpp` `query()`), which routes PML, MEM, plain k-mer presence, plain k-mer count, and count-bv through the coroutine — all byte-identical to the sequential path. |

The k-mer-bitvector add-ons (SSHash-style MPHF IDs and the run-local exact count)
live in `query_kmer_bv.cpp`: the build/load (`build_kmerbv`, `load_kmerbv`) and the
bv queries (`query_kmers_count_bv`, `query_kmers_id_bv`). The small
`kmer_count_from_bv` predecessor/successor helper is an inline in
`include/move_structure.hpp`.

## Coverage grid

"N/A" = not implemented in that style.

| Query | (a) sequential | (b) strands | (c) coroutines |
|-------|----------------|-------------|----------------|
| Exact match (whole-read count) | `handle_count` → `backward_search(...).count()` | `process_latency_hiding` | N/A |
| ZML | `query_zml` | `process_latency_hiding` | N/A |
| PML | `query_pml` | `process_latency_hiding` | `query_pml_coroutine` |
| k-mer presence | `query_all_kmers`/`query_kmers_from` | N/A | `query_kmer_coroutine` |
| k-mer count (no bitvectors) | `query_all_kmers` count branch → `MoveInterval::count()` | N/A | `query_kmer_coroutine` (count mode) |
| k-mer count (bitvector) | `query_kmers_count_bv` → `kmer_count_from_bv` (`--kmer-count --kmer-bv`) | N/A | `query_kmer_coroutine` (`--kmer-bv`) |
| MEMs | `query_mems`/`query_mem_bml` → `extend_bidirectional` | N/A | `query_mem_coroutine` |

`movi query --coroutine` routes PML, MEM, and every k-mer query except the MPHF-id lookup
(`--kmer --kmer-bv`): k-mer presence, plain count (`--kmer-count`), and bitvector count
(`--kmer-count --kmer-bv`) all run through `query_kmer_coroutine`, byte-identical to the
sequential path. The MPHF-id query has no coroutine variant and falls through to the
sequential path (as do ZML and whole-read exact-count).

The k-mer-bitvector index also gives each k-mer a dense, collision-free MPHF id via
`movi query --kmer --kmer-bv` (a drop-in for SSHash's lookup), built by
`movi build-kmerbv`.

## Notes

- MEM has no manual-strand path; its coroutine version is correct but compute-bound
  (the per-run scans dominate, so there is little latency to hide).
- **k-mer has no manual-strand path either.** Latency hiding for k-mer queries is the
  coroutine path, which is byte-identical to the sequential one. The strand scheduler
  serves PML, ZML and exact-count; a k-mer query always takes the sequential path (the
  k-mer flags turn prefetching off in the parser) unless `--coroutine` is given.
- **ftab auto-selection (MEM and k-mer):** the ftab accelerates both MEM and k-mer
  queries (it only accelerates -- results are ftab-independent), and the deepest
  applicable ftab is fastest. When `--ftab-k` is not given, `movi query` auto-selects the
  deepest `ftab.<k>.bin` in the index that fits the query:
  - `--mem`: the deepest ftab **no longer than `--min-mem-length`**. `--ftab-k ==
    min-mem-length` is ideal (the ftab seeds the whole minimum window in one lookup, zero
    bidirectional extend steps); a longer ftab over-extends the BML seed left of the
    intended MEM start, so an explicit `--ftab-k > --min-mem-length` is rejected. The
    all-MEM path (`--min-mem-length <= 1`, `query_all_mems`) is uncapped.
  - `--kmer`: the deepest ftab **no longer than `k`**. `--ftab-k == k` resolves the whole
    k-mer in one lookup; a longer ftab would look up a k-mer longer than the query, so an
    explicit `--ftab-k > k` is rejected. The look-ahead ("fishing") heuristic shrinks its
    step as the ftab deepens and skips the look-ahead entirely at `--ftab-k == k`.

  Build an ftab-12 (`movi ftab --ftab-k 12`) for the best MEM throughput; ecoli100 k-mer
  is about 1.76x faster with ftab-10 than with no ftab. With ftab-12 the
  length-thresholded (BML, `--min-mem-length` above ~1) MEM query is within ~5.5x of k-mer
  at HPRC scale and faster than k-mer at E. coli scale. The one slow case is all-MEM
  enumeration (`--min-mem-length 1`, `query_all_mems`), where the `O(#runs)` skip-count +
  rc-walk scan explodes on large, highly-repeated intervals.
- **Counting k-mers: `--kmer-count` vs `--kmer-count --kmer-bv`.** Both emit
  byte-identical counts. The bitvector path is faster, and its lead widens with k
  (ecoli100, 8000 x 150 bp reads, load-excluded, one thread):

  | k | `--kmer-count` | `--kmer-count --kmer-bv` | speedup |
  |---|---|---|---|
  | 15 | 0.742 s | 0.632 s | 1.17x |
  | 23 | 0.418 s | 0.227 s | 1.84x |
  | 31 | 0.498 s | 0.243 s | 2.05x |
  | 39 | 0.519 s | 0.217 s | 2.39x |
  | 47 | 0.585 s | 0.233 s | 2.51x |
  | 63 | 0.597 s | 0.202 s | 2.96x |

  Plain count pays a per-k-mer backward search that gets dearer as k grows, while the
  bitvector count resolves each present k-mer with a predecessor/successor lookup whose
  cost is essentially independent of k.

  The tradeoff is that `--kmer-bv` needs the per-k `B_k`/`C_k` structures
  (`kmerbv.<k>.*`, built by `movi build-kmerbv` for that one k). On ecoli100 they are
  44 MB at k=15, 68 MB at k=31 and 98 MB at k=63 -- 13-29% of the 341 MB base index
  apiece -- and keeping all twelve k values costs 875 MB, 2.6x the base index. So use
  **`--kmer-count --kmer-bv` when k is known and queried repeatedly**, and **plain
  `--kmer-count` when k varies or is not known ahead of time**: the latter is the
  k-independent path, one index answering every k with no per-k build or storage.
- **Which PML style to use: the coroutine, at every read length.** Measured on one
  exclusive node (hprc94, 35 GB index; the same 25 M bases re-cut into fixed-length reads
  so only read length varies; load-excluded, one thread, 5 reps, medians, CV <= 0.01 on
  warm cells):

  | read length | (a) sequential | (b) strand | (c) coroutine | winner |
  |---|---|---|---|---|
  | 150 bp | 14.18 | 21.15 | **25.13** | coroutine |
  | 500 bp | 15.80 | 25.28 | **30.12** | coroutine |
  | 2500 bp | 17.05 | 27.56 | **32.36** | coroutine |
  | 13000 bp (HiFi) | 17.18 | 28.60 | **32.47** | coroutine |

  *(M bases/s.)* The coroutine leads everywhere, by 1.13x to 1.19x over the strand
  scheduler and 1.77x to 1.90x over sequential, with byte-identical output (verified
  against the sequential path on the 13 kb reads at this scale). Both latency-hiding
  styles pull away from sequential as reads lengthen, which is what the memory-latency
  argument predicts: longer reads mean longer LF chains and more latency to hide.

  Earlier profiling reported a read-length **crossover**, with the strand scheduler
  overtaking the coroutine on 13 kb HiFi. That no longer holds. The strand and sequential
  figures here reproduce the earlier ones closely (28.60 vs ~28.3, 17.18 vs ~17.8), but
  the coroutine roughly doubled -- ~16.5 then, 32.47 now -- so it now leads at 13 kb
  rather than trailing. Treat any "strand wins on HiFi" guidance as superseded.
- The coroutine PML body is an inline reimplementation of the PML loop; the MEM and
  k-mer coroutines instead reuse the shared search primitives.
- The manual-strand path (b) is deliberately **one shared scheduler**, not per-query
  files: `process_latency_hiding` is a single round-robin over `struct Strand` state
  machines that serves PML, ZML and whole-read exact-count together via per-query
  branches. The scheduler itself is query-agnostic, so splitting it per query would
  duplicate it; the coverage grid above maps which queries it serves. (k-mer and MEM
  are not served by the strand scheduler — see the grid.)
- **Output ordering:** all coroutine query output is **content byte-identical** to the
  sequential path (per-read values), and is emitted in deterministic **input order** via
  the `OrderedEmitter` reorder buffer (`output_base_stats`/`output_kmers`/`output_mems`
  write into an `ostringstream`, not straight to `std::cout`). The sequential `--pml`
  stdout stream is itself emitted in batch-loader order, which varies with read length,
  so exact stream order can differ; downstream consumers key on the read name (`>id`).
