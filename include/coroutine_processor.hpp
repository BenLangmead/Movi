#pragma once

// Coroutine latency-hiding query engine ("style (c)" in docs/latency_hiding.md).
// The engine (scheduler + SharedFastqReader + OrderedEmitter + the
// query_*_coroutine bodies) lives in movi_co.cpp / coroutine_processor.cpp; this
// header is the entry point so the main `movi` binary can drive it, not just the
// standalone movi-co binary. Everything the engine needs to know about a query is
// carried in CoroutineQueryOptions rather than file-scope globals.

#include <string>

struct CoroutineQueryOptions {
    bool debug = false;
    bool bpf_output = false;      // write binary BPF output instead of text (PML only)
    bool ordered_output = false;  // emit reads in input order (reorder buffer); else completion order
    bool mem_query = false;       // --mem
    bool kmer_query = false;      // --kmer / --kmer-count
    bool kmer_count = false;      // --kmer-count
    bool kmer_bv = false;         // use the count bitvector for k-mer counts
    int ftab_k = 0;
    int min_mem_length = 0;
    int k = 0;                    // k-mer length (required for kmer_query)
};

// Run the coroutine latency-hiding query over fastq_file against the index in
// index_dir, with `concurrency` in-flight coroutines. Default query mode (no
// mem_query/kmer_query set) is PML.
void run_coroutine_query(const std::string& fastq_file, const std::string& index_dir,
                         int concurrency, const CoroutineQueryOptions& opts);
