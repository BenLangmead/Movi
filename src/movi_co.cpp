// movi-co: thin CLI wrapper around the coroutine latency-hiding query engine.
//
// The engine itself (scheduler + SharedFastqReader + OrderedEmitter + the
// query_*_coroutine bodies + run_coroutine_query) lives in coroutine_processor.cpp
// and also links into the main `movi` binary. This file is just the standalone
// movi-co executable's argument parsing; it builds a CoroutineQueryOptions and
// calls run_coroutine_query. See docs/latency_hiding.md.
//
// Distributed under the GPL3 license.
// See accompanying LICENSE or https://opensource.org/license/gpl-3-0

#include <cstdio>
#include <string>
#include <exception>

#include "coroutine_processor.hpp"

using std::string;
using std::exception;

static void print_usage(const char* program_name) {
    fprintf(stderr, "Usage: %s [--debug] [--bpf] [--reorder] [query opts] <fastq_file> <index_dir> [concurrency]\n", program_name);
    fprintf(stderr, "  --debug:     Enable debug output\n");
    fprintf(stderr, "  --bpf:       Write binary BPF output (like mainline Movi) instead of text\n");
    fprintf(stderr, "  --reorder:   Emit reads in input order (enable the reorder buffer; default is completion order)\n");
    fprintf(stderr, "  --mem --ftab-k N --min-mem-length N:  MEM query mode\n");
    fprintf(stderr, "  --kmer | --kmer-count, -k/--k-length N [--ftab-k N]:  k-mer query mode\n");
    fprintf(stderr, "               (default query mode is PML)\n");
    fprintf(stderr, "  fastq_file: Input FASTQ file with reads\n");
    fprintf(stderr, "  index_dir:  Directory containing the Movi index\n");
    fprintf(stderr, "  concurrency: Number of concurrent coroutines (default: 1)\n");
}

int main(int argc, char* argv[]) {
    int arg_idx = 1;
    CoroutineQueryOptions opts;

    // Parse leading flags in any order. Accept both long (--foo) and short (-k)
    // flags; a positional (fastq path or concurrency) never starts with '-'.
    while (arg_idx < argc && argv[arg_idx][0] == '-' && string(argv[arg_idx]) != "-") {
        string flag = argv[arg_idx];
        if (flag == "--debug") { opts.debug = true; }
        else if (flag == "--bpf") { opts.bpf_output = true; }
        else if (flag == "--reorder") { opts.ordered_output = true; }
        else if (flag == "--mem") { opts.mem_query = true; }
        else if (flag == "--kmer") { opts.kmer_query = true; }
        else if (flag == "--kmer-count") { opts.kmer_query = true; opts.kmer_count = true; }
        else if (flag == "--kmer-bv") { opts.kmer_bv = true; }
        else if (flag == "-k" || flag == "--k-length") { if (arg_idx + 1 >= argc) { print_usage(argv[0]); return 1; } opts.k = std::stoi(argv[++arg_idx]); }
        else if (flag == "--ftab-k") { if (arg_idx + 1 >= argc) { print_usage(argv[0]); return 1; } opts.ftab_k = std::stoi(argv[++arg_idx]); }
        else if (flag == "--min-mem-length") { if (arg_idx + 1 >= argc) { print_usage(argv[0]); return 1; } opts.min_mem_length = std::stoi(argv[++arg_idx]); }
        else { fprintf(stderr, "Unknown flag: %s\n", flag.c_str()); print_usage(argv[0]); return 1; }
        arg_idx++;
    }

    // Validate query-mode flag combinations.
    if (opts.mem_query && opts.kmer_query) {
        fprintf(stderr, "Error: --mem and --kmer/--kmer-count are mutually exclusive\n");
        return 1;
    }
    if (opts.kmer_query && opts.k < 1) {
        fprintf(stderr, "Error: --kmer/--kmer-count requires -k/--k-length N (N >= 1)\n");
        return 1;
    }

    // Check remaining positional arguments.
    if (argc - arg_idx < 2 || argc - arg_idx > 3) {
        fprintf(stderr, "argc = %d\n", argc);
        print_usage(argv[0]);
        return 1;
    }

    string fastq_file{argv[arg_idx]}, index_dir{argv[arg_idx + 1]};
    int concurrency = 1;

    if (argc - arg_idx == 3) {
        try {
            concurrency = std::stoi(argv[arg_idx + 2]);
            if (concurrency < 1) {
                fprintf(stderr, "Error: Concurrency must be at least 1\n");
                return 1;
            }
        } catch (const exception& e) {
            fprintf(stderr, "Error: Invalid concurrency value: %s\n", argv[arg_idx + 2]);
            return 1;
        }
    }

    try {
        run_coroutine_query(fastq_file, index_dir, concurrency, opts);
    } catch (const exception& e) {
        fprintf(stderr, "Error: %s\n", e.what());
        return 1;
    }
    return 0;
}
