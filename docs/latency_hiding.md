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
| (b) manual strands | `read_processor.cpp` (`process_latency_hiding`, `kmer_search_latency_hiding`); `struct Strand` in `read_processor.hpp` |
| (c) coroutines | `coroutine_processor.cpp` (the `query_*_coroutine` bodies + scheduler + `run_coroutine_query`; guarded to `MODE==6 && COLOR_MODE==0`). Reached via **`movi query --coroutine`** (dispatch in `movi.cpp` `query()`), which routes PML, MEM, plain k-mer presence, and count-bv through the coroutine — all byte-identical to the sequential path. |

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
| k-mer presence | `query_all_kmers`/`query_kmers_from` | `kmer_search_latency_hiding` (disabled by default) | `query_kmer_coroutine` |
| k-mer count (no bitvectors) | `query_all_kmers` count branch → `MoveInterval::count()` | N/A | not routed¹ |
| k-mer count (bitvector) | `query_kmers_count_bv` → `kmer_count_from_bv` (`--kmer-count --kmer-bv`) | N/A | `query_kmer_coroutine` (`--kmer-bv`) |
| MEMs | `query_mems`/`query_mem_bml` → `extend_bidirectional` | N/A | `query_mem_coroutine` |

¹ `query_kmer_coroutine` contains a plain-count branch, but `movi query --coroutine
--kmer-count` (no `--kmer-bv`) is **not** dispatched to it (`movi.cpp` routes it to the
sequential path): the coroutine's `found` tally diverges from the sequential count. Fixing
that divergence is a tracked coroutine-correctness task; until then only presence and
bitvector-count are routed to the coroutine. The MPHF-id query (`--kmer --kmer-bv`) has no
coroutine variant at all and is likewise routed to the sequential path.

The k-mer-bitvector index also gives each k-mer a dense, collision-free MPHF id via
`movi query --kmer --kmer-bv` (a drop-in for SSHash's lookup), built by
`movi build-kmerbv`.

## Notes

- MEM has no manual-strand path; its coroutine version is correct but compute-bound
  (the per-run scans dominate, so there is little latency to hide).
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
