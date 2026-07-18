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
| (c) coroutines | `coroutine_processor.cpp` (the `query_*_coroutine` bodies + scheduler + `run_coroutine_query`; guarded to `MODE==6 && COLOR_MODE==0`). Reached two ways: **`movi query --coroutine`** (dispatch in `movi.cpp` `query()`, for the queries whose output is byte-identical to sequential — MEM, plain k-mer presence, count-bv) and the standalone `movi-co` binary (thin CLI wrapper in `movi_co.cpp`, all coroutine queries incl. PML) |

The k-mer-bitvector add-ons (SSHash-style MPHF IDs and the run-local exact count)
all live in `query_kmer_bv.cpp`: the build/load (`build_kmerbv`, `load_kmerbv`),
the bv queries (`query_kmers_count_bv`, `query_kmers_id_bv`), and the
`kmer_count_from_bv` helper.

## Coverage grid

"before June" = predates this work. "N/A" = not implemented in that style.

| Query | (a) sequential | (b) strands | (c) coroutines |
|-------|----------------|-------------|----------------|
| Exact match (whole-read count) | `handle_count` → `backward_search(...).count()` | `process_latency_hiding` | N/A |
| ZML | `query_zml` | `process_latency_hiding` | N/A |
| PML | `query_pml` | `process_latency_hiding` | `query_pml_coroutine` |
| k-mer presence | `query_all_kmers`/`query_kmers_from` | `kmer_search_latency_hiding` (disabled by default) | `query_kmer_coroutine` |
| k-mer count (no bitvectors) | `query_all_kmers` count branch → `MoveInterval::count()` | N/A | `query_kmer_coroutine` (count mode) |
| k-mer count (bitvector) | `query_kmers_count_bv` → `kmer_count_from_bv` (`--kmer-count --kmer-bv`) | N/A | `query_kmer_coroutine` (`--kmer-bv`) |
| MEMs | `query_mems`/`query_mem_bml` → `extend_bidirectional` | N/A | `query_mem_coroutine` |

The k-mer-bitvector index also gives each k-mer a dense, collision-free MPHF id via
`movi query --kmer --kmer-bv` (a drop-in for SSHash's lookup), built by
`movi build-kmerbv`.

## Notes

- MEM has no manual-strand path; its coroutine version is correct but compute-bound
  (the per-run scans dominate, so there is little latency to hide).
- The coroutine PML body is an inline reimplementation of the PML loop; the MEM and
  k-mer coroutines instead reuse the shared search primitives.
- **Fold status (Movi 2.5):** the coroutine engine links into the main `movi` binary;
  `movi query --coroutine` runs MEM, plain k-mer presence, and count-bv through it
  (verified byte-identical to the sequential path). Three cases are NOT yet routed via
  `movi query --coroutine` and fall back to the sequential/strand path (the standalone
  `movi-co` still runs them all): plain **k-mer count** (`--kmer-count` without `--bv`
  — the coroutine over-counts the `found` tally vs sequential), the **MPHF-id query**
  (`--kmer --kmer-bv` — no coroutine implementation), and **PML** (coroutine output
  format and read orientation not yet reconciled with the mainline binary format).
