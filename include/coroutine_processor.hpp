#pragma once

// Coroutine latency-hiding query engine ("style (c)" in docs/latency_hiding.md).
// The engine (scheduler + SharedFastqReader + OrderedEmitter + the
// query_*_coroutine bodies) lives in coroutine_processor.cpp; this header is the
// entry point so the main `movi` binary can drive it (via `movi query --coroutine`),
// not just the standalone movi-co binary. Everything the engine needs to know about
// a query is carried in CoroutineQueryOptions rather than file-scope globals.

#include <string>
#include <iosfwd>

class MoveStructure;  // forward declaration; the core entry point takes it by reference

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

// Entry point: drive the coroutine scheduler over an already-loaded MoveStructure
// (index deserialized, ftab/kmerbv set up by the caller). The main movi binary
// calls this from its `movi query --coroutine` dispatch. Default query mode (no
// mem_query/kmer_query set) is PML. `out` is the output destination (a file stream
// from output_files, or std::cout).
void run_coroutine_query(MoveStructure& mv, const std::string& fastq_file,
                         int concurrency, const CoroutineQueryOptions& opts, std::ostream& out);
