#pragma once

// Coroutine latency-hiding query engine ("style (c)" in docs/latency_hiding.md).
// The engine (scheduler + SharedFastqReader + OrderedEmitter + the
// query_*_coroutine bodies) lives in coroutine_processor.cpp; this header is the
// entry point so the main `movi` binary can drive it (via `movi query --coroutine`).
// Everything the engine needs to know about a query is carried in
// CoroutineQueryOptions rather than file-scope globals.

#include <string>
#include <iosfwd>

class MoveStructure;  // forward declaration; the core entry point takes it by reference

struct CoroutineQueryOptions {
    bool debug = false;
    bool ordered_output = false;  // emit reads in input order (reorder buffer); else completion order
    bool mem_query = false;       // run the MEM coroutine (else k-mer if kmer_query, else PML)
    bool kmer_query = false;      // run the k-mer coroutine
    // Per-query parameters (k, ftab_k, min-mem-length, count/bv mode) are read from the
    // MoveStructure's MoviOptions by the query bodies, so they are not duplicated here.
};

// Entry point: drive the coroutine scheduler over an already-loaded MoveStructure
// (index deserialized, ftab/kmerbv set up by the caller). The main movi binary
// calls this from its `movi query --coroutine` dispatch. Default query mode (no
// mem_query/kmer_query set) is PML. `out` is the output destination (a file stream
// from output_files, or std::cout).
//
// `concurrency` is the per-thread in-flight coroutine count (latency hiding within a
// thread); `nthreads` is the number of worker threads (cores). They are independent:
// total in-flight coroutines = nthreads * concurrency. nthreads == 1 runs the original
// single-threaded path with no locking overhead.
void run_coroutine_query(MoveStructure& mv, const std::string& fastq_file,
                         int concurrency, int nthreads,
                         const CoroutineQueryOptions& opts, std::ostream& out);
