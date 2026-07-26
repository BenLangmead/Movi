#include <cstdint>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <cstddef>
#include <unistd.h>
#include <sys/stat.h>

#include <omp.h>
#include "sdsl_wrapper.hpp"
#include "cxxopts.hpp"

#include "utils.hpp"
#include "move_structure.hpp"
#include "move_query.hpp"
#include "read_processor.hpp"
#include "movi_options.hpp"
#include "movi_parser.hpp"
#include "coroutine_processor.hpp"
#include "classifier.hpp"
#include "batch_loader.hpp"
#include "mmap_batch_source.hpp"

// Function to handle PML/ZML processing for a single read
uint64_t handle_pml_zml(MoveQuery& mq, MoviOptions& movi_options,
                        MoveStructure& mv_, OutputFiles& output_files, Classifier& classifier) {
    uint64_t ff_count = 0;

    if (movi_options.is_pml()) {
        ff_count += mv_.query_pml(mq);
    } else if (movi_options.is_zml()) {
        ff_count += mv_.query_zml(mq);
    }

    #pragma omp critical
    {
        if (movi_options.is_classify()) {
            std::vector<uint16_t> matching_lens_16(mq.get_matching_lengths().begin(), mq.get_matching_lengths().end());
            // Classification with 32 bits pmls works if the pml values are less than 2^16.

            bool found = classifier.classify(mq.get_query_id(), matching_lens_16, movi_options);

            if (movi_options.is_filter() && !movi_options.is_no_output()) {
                if (found && !movi_options.is_invert()) {
                    output_read(mq);
                } else if (!found && movi_options.is_invert()) {
                    output_read(mq);
                }
            }
        }

        if (movi_options.write_output_allowed()) {
            std::ostream& mls_dest = movi_options.write_stdout_enabled()
                ? static_cast<std::ostream&>(std::cout)
                : static_cast<std::ostream&>(output_files.mls_file);
            output_base_stats(DataType::match_length, movi_options.write_stdout_enabled(), mls_dest, mq);

            if (movi_options.is_get_sa_entries()) {
                output_base_stats(DataType::sa_entry, movi_options.write_stdout_enabled(), output_files.sa_entries_file, mq);
            }
        }
    }

    return ff_count;
}

// Function to handle count processing for a single read
void handle_count(MoveQuery& mq, MoviOptions& movi_options,
                  MoveStructure& mv_, OutputFiles& output_files) {
    int32_t pos_on_r = mq.query().length() - 1;
    // uint64_t match_count = mv_.backward_search(query_seq, pos_on_r);
    // if (pos_on_r != 0) pos_on_r += 1;
    uint64_t match_count = mv_.query_backward_search(mq, pos_on_r);

    if (movi_options.write_output_allowed()) {
        #pragma omp critical
        {
            output_counts(movi_options.write_stdout_enabled(), output_files.matches_file, mq.query().length(), pos_on_r, match_count, mq);
        }
    }
}

// Function to handle k-mer processing for a single read
void handle_kmer(MoveQuery& mq, MoviOptions& movi_options,
                 MoveStructure& mv_, OutputFiles& output_files) {
    mv_.query_all_kmers(mq, movi_options.is_kmer_count());

    // For the SSHash-style aggregate report, accumulate window-accurate tallies:
    // every length-k window is positive (found), invalid (non-ACGT), or negative.
    if (movi_options.is_output_format_sshash()) {
        size_t k = movi_options.get_k();
        size_t L = mq.query().length();
        size_t total = (L >= k) ? (L - k + 1) : 0;
        size_t invalid = count_invalid_kmer_windows(mq.query(), k);
        #pragma omp atomic
        mv_.kmer_stats.agg_num_kmers += total;
        #pragma omp atomic
        mv_.kmer_stats.agg_invalid += invalid;
        #pragma omp atomic
        mv_.kmer_stats.agg_positive += mq.found_kmer_count;
    }

    // Emit per-k-mer output in both presence and count modes.  (Count mode now
    // populates the same per-k-mer string via add_kmer, so it is no longer
    // suppressed here.)  The sshash format suppresses per-read lines (handled in
    // output_kmers) and prints its aggregate report once at the end.
    if (movi_options.write_output_allowed()) {
        #pragma omp critical
        {
            output_kmers(movi_options.write_stdout_enabled(), output_files.kmer_file,
                         mq.query().length() - movi_options.get_k() + 1, mq, movi_options);
        }
    }
}

// Function to handle mem processing for a single read
void handle_mem(MoveQuery& mq, MoviOptions& movi_options,
                MoveStructure& mv_, OutputFiles& output_files) {
    mv_.query_mems(mq);
    if (movi_options.write_output_allowed()) {
        #pragma omp critical
        {
            output_mems(movi_options.write_stdout_enabled(), output_files.mems_file, mq);
        }
    }
}

// Helper function to setup input file stream
void setup_input_file(std::ifstream& input_file, const std::string& read_file) {
    if (read_file == "-") {
        input_file.copyfmt(std::cin);
        input_file.clear(std::cin.rdstate());
        input_file.basic_ios<char>::rdbuf(std::cin.rdbuf());
    } else {
        // This check is already done in the parser too.
        if (!std::filesystem::exists(read_file)) {
            throw std::runtime_error(ERROR_MSG("The input file " + read_file + " does not exist."));
        }
        input_file.open(read_file.c_str());
    }
}

// Function to load color table/document sets
void load_color_table(MoveStructure& mv_, MoviOptions& movi_options) {
    auto begin = std::chrono::system_clock::now();

    if (movi_options.is_full_color()) {
        mv_.fill_run_offsets();
        std::string fname = movi_options.get_index_dir() + "/doc_pats.bin";
        mv_.deserialize_doc_pats(fname);
    } else {
        if (movi_options.is_doc_sets_vector_of_vectors()) {
            if (!movi_options.is_freq_compressed() and !movi_options.is_tree_compressed()) {
                std::string fname = movi_options.get_index_dir() + "/doc_sets.bin";
                mv_.deserialize_doc_sets(fname);
            } else if (movi_options.is_freq_compressed()) {
                std::string fname = movi_options.get_index_dir() + "/compress_doc_sets.bin";
                mv_.deserialize_doc_sets(fname);
            } else if (movi_options.is_tree_compressed()) {
                std::string fname = movi_options.get_index_dir() + "/tree_doc_sets.bin";
                mv_.deserialize_doc_sets(fname);
            }
        } else {
            mv_.deserialize_doc_sets_flat();
        }
        mv_.load_document_info();
    }
    mv_.initialize_classify_cnts();

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    TIMING_MSG(elapsed, "loading the document sets");
}

void build_ftab(MoveStructure& mv_, MoviOptions& movi_options) {
    if (movi_options.is_multi_ftab() and movi_options.get_ftab_k() > 1) {
        int max_ftab = movi_options.get_ftab_k();
        for (int i = 2; i <= max_ftab; i++) {
            movi_options.set_ftab_k(i);
            mv_.build_ftab();
            mv_.write_ftab();
            SUCCESS_MSG("The ftab table for k = " + std::to_string(i) + " is built and stored in the index directory.");
        }
    } else if (movi_options.get_ftab_k() > 1) {
        mv_.build_ftab();
        mv_.write_ftab();
    }
}

void color(MoveStructure& mv_, MoviOptions& movi_options) {
    mv_.load_document_info();

    auto begin = std::chrono::system_clock::now();
    if (movi_options.is_full_color()) {
        // Build document patterns (full information)
        mv_.fill_run_offsets();
        mv_.build_doc_pats();
        mv_.serialize_doc_pats(movi_options.get_index_dir() + "/doc_pats.bin");

        mv_.build_doc_sets();
        SUCCESS_MSG("Done building document sets.");
        mv_.serialize_doc_sets(movi_options.get_index_dir() + "/doc_sets.bin");
    } else {
        if (!movi_options.is_compressed()) {
            mv_.fill_run_offsets();

            std::string doc_pats_name = movi_options.get_index_dir() + "/doc_pats.bin";
            std::ifstream doc_pats_file(doc_pats_name);
            if (doc_pats_file.good()) {
                mv_.deserialize_doc_pats(doc_pats_name);
                INFO_MSG("Done reading document pattern information");
            } else {
                INFO_MSG("Doc patterns are not available, building...");

                auto begin = std::chrono::system_clock::now();
                mv_.build_doc_pats();
                auto end = std::chrono::system_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
                TIMING_MSG(elapsed, "building the document patterns");
                mv_.serialize_doc_pats(movi_options.get_index_dir() + "/doc_pats.bin");
            }
            mv_.build_doc_sets();
            if (movi_options.is_doc_sets_vector_of_vectors()) {
                mv_.serialize_doc_sets(movi_options.get_index_dir() + "/doc_sets.bin");
            } else {
                mv_.flat_and_serialize_colors_vectors();
            }
        } else {
            mv_.deserialize_doc_sets(movi_options.get_index_dir() + "/doc_sets.bin");
            mv_.compress_doc_sets();
            mv_.serialize_doc_sets(movi_options.get_index_dir() + "/compress_doc_sets.bin");

            //mv_.build_doc_set_similarities();
            //mv_.build_tree_doc_sets();
            //mv_.serialize_doc_sets(movi_options.get_index_dir() + "/tree_doc_sets.bin");
        }
    }

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    TIMING_MSG(elapsed, "building colors");
}

void query(MoveStructure& mv_, MoviOptions& movi_options) {

    // The ftab accelerates MEM and k-mer queries (it only accelerates -- results are
    // ftab-independent), and the deepest applicable ftab is fastest. So when the user
    // runs `movi query --mem` or `--kmer` without an explicit --ftab-k, auto-select the
    // deepest ftab.<k>.bin in the index that fits the query, giving fast queries by
    // default without the user needing to know which depths were built.
    //
    // Cap on the ftab depth (a longer ftab over-extends the seed and is wrong/unusable):
    //   - MEM (length-thresholded, BML): ftab-k <= --min-mem-length; ftab-k ==
    //     min-mem-length is ideal (whole seed in one lookup). all-MEM (min-mem <= 1) is
    //     uncapped.
    //   - k-mer: ftab-k <= k; ftab-k == k resolves the whole k-mer in one lookup.
    // An explicit --ftab-k is validated against the same caps below.
    if (movi_options.get_ftab_k() == 0 && (movi_options.is_mem() || movi_options.is_kmer())) {
        uint32_t cap = std::numeric_limits<uint32_t>::max();
        std::string what = "query";
        if (movi_options.is_mem()) {
            const uint32_t mm = movi_options.get_min_mem_length();
            if (mm > 1) cap = mm;
            what = "MEM";
        } else {
            cap = static_cast<uint32_t>(movi_options.get_k());
            what = "k-mer";
        }
        uint32_t best = 0;
        std::error_code ec;
        for (const auto& entry : std::filesystem::directory_iterator(movi_options.get_index_dir(), ec)) {
            const std::string name = entry.path().filename().string();
            // Match "ftab.<k>.bin" with a numeric k (skip e.g. ftab.chosen.k).
            if (name.rfind("ftab.", 0) == 0 && name.size() > 9 &&
                name.compare(name.size() - 4, 4, ".bin") == 0) {
                const std::string kstr = name.substr(5, name.size() - 9);
                if (!kstr.empty() && std::all_of(kstr.begin(), kstr.end(), ::isdigit)) {
                    uint32_t fk = static_cast<uint32_t>(std::stoul(kstr));
                    if (fk > best && fk <= cap) best = fk;
                }
            }
        }
        if (best > 1) {
            movi_options.set_ftab_k(best);
            INFO_MSG(what + ": auto-selected ftab-k=" + std::to_string(best) +
                     " (deepest ftab that fits the query).");
        }
    }

    // Reject an explicit --ftab-k longer than the k-mer length (a longer ftab looks up a
    // k-mer longer than the query, over-constraining the seed).
    if (movi_options.is_kmer() && movi_options.get_ftab_k() != 0 &&
        movi_options.get_ftab_k() > static_cast<uint32_t>(movi_options.get_k())) {
        throw std::runtime_error(ERROR_MSG("For k-mer queries, --ftab-k (" +
            std::to_string(movi_options.get_ftab_k()) + ") must not exceed k (" +
            std::to_string(movi_options.get_k()) + "). Use --ftab-k <= k (ftab-k == k resolves the "
            "whole k-mer in one lookup)."));
    }

    if (movi_options.get_ftab_k() != 0) {
        mv_.read_ftab();
        INFO_MSG("Ftab was read!");
    }

    // Counting with the k-mer bitvector resolves each present k-mer with a
    // predecessor/successor lookup instead of a per-k-mer backward search, so it is
    // faster and its lead grows with k (about 2x at k=31 and 3x at k=63 on ecoli100),
    // for byte-identical counts. It needs the per-k kmerbv.<k>.* structures, so point
    // the user at it only when this index already has them for the k being queried.
    if (movi_options.is_kmer_count() && !movi_options.is_kmer_bv()) {
        const std::string bv_meta = movi_options.get_index_dir() + "/kmerbv." +
                                    std::to_string(movi_options.get_k()) + ".meta";
        if (std::filesystem::exists(bv_meta)) {
            INFO_MSG("This index has a k-mer bitvector for k=" +
                     std::to_string(movi_options.get_k()) +
                     "; adding --kmer-bv counts faster with identical output.");
        }
    }

#if TALLY_MODES
    if (movi_options.is_zml() or movi_options.is_count()) {
        //TODO: Implement tally modes for zml and count queries.
        throw std::runtime_error(ERROR_MSG("This query type is not suppported by the sampled modes yet."));
    }
#endif

    if (movi_options.is_mem()) {
        if (!movi_options.no_prefetch()) {
            movi_options.set_prefetch(false);
            WARNING_MSG("MEM finding does not support prefetching. Continuing with prefetching disabled.");
        }
        if (movi_options.get_ftab_k() == 0) {
            throw std::runtime_error(ERROR_MSG("MEM finding requires an ftab, but none was found in the index. Build one with `movi ftab --ftab-k 12` (a deeper ftab means faster MEM search; ftab-12 is recommended). It is then auto-selected, or pass --ftab-k <k> explicitly."));
        } else if (movi_options.get_min_mem_length() > 1 &&
                   movi_options.get_ftab_k() > movi_options.get_min_mem_length()) {
            // ftab-k must not exceed min-mem-length for the BML search; a longer ftab
            // over-extends the seed left of the intended MEM start. (Auto-select already
            // caps at min-mem-length; this catches an explicit, out-of-range --ftab-k.)
            throw std::runtime_error(ERROR_MSG("For length-thresholded MEM search, --ftab-k (" +
                std::to_string(movi_options.get_ftab_k()) + ") must not exceed --min-mem-length (" +
                std::to_string(movi_options.get_min_mem_length()) + "). A longer ftab over-extends the seed; "
                "use --ftab-k <= min-mem-length (ftab-k == min-mem-length is fastest)."));
        }
    }

    omp_set_num_threads(movi_options.get_threads());
    omp_set_nested(0);

    // This is for the no-prefetch mode (the prefetch mode has its own classifier)
    Classifier classifier;
    mv_.set_classifier(&classifier);
    if (movi_options.is_classify()) {
        classifier.initialize_report_file(movi_options);
    }

    std::ifstream input_file;
    setup_input_file(input_file, movi_options.get_read_file());

    // Maps the input when it is a raw regular FASTA/FASTQ file, so the batch loops below
    // serialize only a record-boundary claim and build their batches in the worker
    // thread. It reports is_mmap() == false for gzip, stdin and pipes, and those keep
    // reading input_file under the lock.
    MmapBatchSource batch_source(movi_options.get_read_file());

    OutputFiles output_files;
    open_output_files(movi_options, output_files);
    mv_.set_output_files(&output_files);

#if (MODE == 6 || MODE == 7 || MODE == 8) && COLOR_MODE == 0
    // Coroutine latency-hiding dispatch (threshold index modes, non-color only). Route
    // the queries whose coroutine output is byte-identical to the sequential path:
    // PML, MEM, and every k-mer query except the MPHF-id lookup. Among k-mer queries
    // that is presence (F,F), plain count (T,F), and bitvector count (T,T); only the
    // MPHF-id query (--kmer --kmer-bv, i.e. bv without count) is excluded, since it has
    // no coroutine implementation and falls through to the sequential path. ZML and
    // exact-count also fall through, having no coroutine variant.
    if (movi_options.is_coroutine()) {
        const bool coroutine_routable =
            movi_options.is_pml() || movi_options.is_mem() ||
            (movi_options.is_kmer() && !(movi_options.is_kmer_bv() && !movi_options.is_kmer_count()));
        if (coroutine_routable) {
            CoroutineQueryOptions copts;
            copts.mem_query      = movi_options.is_mem();
            copts.kmer_query     = movi_options.is_kmer();
            copts.ordered_output = true;  // emit in input order, matching the sequential path
            // k, ftab_k, min-mem-length and count/bv mode are read from mv_'s MoviOptions
            // by the coroutine bodies, so they need not be copied into CoroutineQueryOptions.
            std::ostream& dest = movi_options.is_mem()  ? output_files.mems_file
                               : movi_options.is_kmer() ? output_files.kmer_file
                                                        : output_files.mls_file;
            run_coroutine_query(mv_, movi_options.get_read_file(),
                                static_cast<int>(movi_options.get_strands()),
                                static_cast<int>(movi_options.get_threads()), copts,
                                movi_options.write_stdout_enabled() ? std::cout : dest);
            close_output_files(movi_options, output_files);
            return;
        }
        // --coroutine was requested for a query the coroutine engine cannot reproduce
        // byte-identically; warn so a benchmark does not silently attribute sequential
        // timings to the coroutine path.
        std::cerr << "[movi] note: --coroutine has no coroutine implementation for this query "
                     "(MPHF-id --kmer --kmer-bv, ZML, or exact count); running the sequential "
                     "path instead." << std::endl;
    }
#else
    if (movi_options.is_coroutine()) {
        std::cerr << "[movi] note: --coroutine latency hiding is only implemented for the "
                     "threshold index modes (regular-, sampled-, or blocked-thresholds, non-color); "
                     "running the sequential path instead." << std::endl;
    }
#endif

    uint64_t total_ff_count = 0;

    auto begin = std::chrono::system_clock::now();

    if (!movi_options.no_prefetch()) {

        ReadProcessor rp(mv_, movi_options.get_strands(), movi_options.is_verbose(), movi_options.is_reverse(), output_files, classifier);

        // Reads claimed per strand batch (the per-thread work-claim granularity for the
        // strand scheduler). Default 4*strands; MOVI_STRAND_BATCH overrides it for
        // benchmarking the batch-size vs load-balance tradeoff, mirroring MOVI_CO_BATCH
        // on the coroutine path.
        size_t strand_batch_reads = 4 * movi_options.get_strands();
        if (const char* e = std::getenv("MOVI_STRAND_BATCH")) {
            long v = std::atol(e);
            if (v > 0) strand_batch_reads = static_cast<size_t>(v);
        }

#pragma omp parallel
        {
            BatchLoader reader;

            // Iterates over batches of data until none left
            while (true) {
                bool valid_batch = true;
                // With a mapped input only the record-boundary claim is serialized and the
                // batch is built in this thread; otherwise the whole parse stays under the
                // lock, as it must for gzip and stdin.
                if (batch_source.is_mmap()) {
                    const char* p = nullptr; const char* e = nullptr;
                    #pragma omp critical // one claimer at a time
                    {
                        valid_batch = batch_source.claim(p, e, 1000, strand_batch_reads);
                    }
                    if (!valid_batch) {
                        break;
                    }
                    reader.loadBatchFromRange(p, e, batch_source.format());
                } else {
                    #pragma omp critical // one reader at a time
                    {
                        valid_batch = reader.loadBatch(input_file, 1000, strand_batch_reads);
                    }
                    if (!valid_batch) {
                        break;
                    }
                }

                // The strand scheduler serves PML, ZML and whole-read exact-count. k-mer
                // and MEM are not routed here: the k-mer query flags turn prefetching off
                // in the parser, so a k-mer query reaches the sequential path below, and
                // its latency-hiding variant is the coroutine one dispatched above.
#if TALLY_MODES
                rp.process_latency_hiding_tally(reader);
#else
                rp.process_latency_hiding(reader);
#endif
            }
        }

        // TODO: total ff_count is not correct in the prefetch mode.
        total_ff_count += rp.get_total_ff_count();
        rp.end_process();

        SUCCESS_MSG(format_number_with_commas(rp.get_read_processed()) + " reads are processed.");

    } else {
        if (!movi_options.is_kmer()) {
            // k-mer queries always arrive here (the parser turns prefetching off for
            // them), so the message would be noise; it is meant for a query that could
            // have used the strand scheduler but was asked not to.
            INFO_MSG("Latency hiding is disabled...");
        }

        uint64_t read_processed = 0;

        #pragma omp parallel
        {
            BatchLoader reader;

            // Iterates over batches of data until none left
            while (true) {
                bool valid_batch = true;
                if (batch_source.is_mmap()) {
                    const char* p = nullptr; const char* e = nullptr;
                    #pragma omp critical // one claimer at a time
                    {
                        valid_batch = batch_source.claim(p, e, 1000, 1);
                    }
                    if (!valid_batch) break;
                    reader.loadBatchFromRange(p, e, batch_source.format());
                } else {
                    #pragma omp critical // one reader at a time
                    {
                        valid_batch = reader.loadBatch(input_file, 1000, 1);
                    }
                    if (!valid_batch) break;
                }
                
                Read read_struct;
                bool valid_read = false;

                // Iterates over reads in a single batch
                while (true) {
                    valid_read = reader.grabNextRead(read_struct);
                    if (!valid_read) break;

                    #pragma omp atomic
                    read_processed += 1 ;

                    #pragma omp critical
                    {
                        if (read_processed % 1000 == 0) {
                            QUERY_PROGRESS_MSG("Number of reads processed: " + format_number_with_commas(read_processed));
                        }
                    }

                    // std::string query_seq = seq->seq.s;
                    std::string query_seq = std::string(read_struct.seq);

                    if (movi_options.is_reverse())
                        std::reverse(query_seq.begin(), query_seq.end());

                    MoveQuery mq = MoveQuery(query_seq);
                    mq.set_query_id(read_struct.id);

                    if (movi_options.is_pml() or movi_options.is_zml()) {
                        total_ff_count += handle_pml_zml(mq, movi_options, mv_, output_files, classifier);

                    } else if (movi_options.is_count()) {

                        handle_count(mq, movi_options, mv_, output_files);

                    } else if (movi_options.is_kmer()) {

                        handle_kmer(mq, movi_options, mv_, output_files);


                    } else if (movi_options.is_mem()) {
                        handle_mem(mq, movi_options, mv_, output_files);
                    }

                    if (movi_options.is_logs()) {
                        if (movi_options.write_output_allowed()) {
                            #pragma omp critical
                            {
                                output_logs(output_files.costs_file, output_files.scans_file, output_files.fastforwards_file, mq);
                            }
                        }
                    }
                }
            }
        }
        SUCCESS_MSG(format_number_with_commas(read_processed) + " reads are processed.");

    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    TIMING_MSG(elapsed, "processing the reads");

    if (movi_options.write_output_allowed()) {
        print_query_stats(movi_options, total_ff_count, mv_);
    }

    // SSHash-style aggregate query report (membership tallies), printed once.
    if (movi_options.is_kmer() && movi_options.is_output_format_sshash()) {
        mv_.kmer_stats.print_sshash_report();
    }

    if (movi_options.is_classify()) {
        classifier.close_report_file();
    }

    close_output_files(movi_options, output_files);
}

void view(MoviOptions& movi_options) {
    std::ifstream mls_file(movi_options.get_mls_file(), std::ios::in | std::ios::binary);
    if (!mls_file.good()) {
        throw std::runtime_error(ERROR_MSG("Failed to open the MLS file: " + movi_options.get_mls_file()));
    }

    mls_file.seekg(0, std::ios::beg);

    Classifier classifier;
    if (movi_options.is_classify()) {
        classifier.initialize_report_file(movi_options);
    }

    uint8_t entry_size = 32;
    if (!movi_options.is_no_header()) {
        // Read BPF header
        BPFHeader header;
        mls_file.read(reinterpret_cast<char*>(&header), sizeof(header));
        if (header.magic != BPF_MAGIC) {
            throw std::runtime_error("Invalid BPF header.");
        }
        if (header.version != BPF_VERSION_MAJOR) {
            throw std::runtime_error("Invalid BPF version.");
        }
        entry_size = header.entry_size;
    } else {
        if (movi_options.is_small_pml_lens()) {
            entry_size = 16;
        } else if (movi_options.is_large_pml_lens()) {
            entry_size = 64;
        }
    }

    while (true) {
        uint16_t st_length = 0;
        mls_file.read(reinterpret_cast<char*>(&st_length), sizeof(st_length));
        if (mls_file.eof()) break;

        std::string read_name;
        read_name.resize(st_length);
        mls_file.read(reinterpret_cast<char*>(&read_name[0]), st_length);
        read_name.erase(std::find(read_name.begin(), read_name.end(), '\0'), read_name.end());
        std::cout << ">" << read_name << "\n";
        uint64_t mq_pml_lens_size = 0;
        mls_file.read(reinterpret_cast<char*>(&mq_pml_lens_size), sizeof(mq_pml_lens_size));

        // TODO: There is a lot of duplicate code here to handle 16 and 32 bits pml variations.
        if (entry_size == 16) {
            std::vector<uint16_t> pml_lens;
            pml_lens.resize(mq_pml_lens_size);
            mls_file.read(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));

            for (int64_t i = mq_pml_lens_size - 1; i >= 0; i--) {
                std::cout << pml_lens[i] << " ";
            }

            std::cout << "\n";

            if (movi_options.is_classify()) {
                classifier.classify(read_name, pml_lens, movi_options);
            }

        } else if (entry_size == 64) {

            std::vector<uint64_t> pml_lens;
            pml_lens.resize(mq_pml_lens_size);
            mls_file.read(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));

            for (int64_t i = mq_pml_lens_size - 1; i >= 0; i--) {
                std::cout << pml_lens[i] << " ";
            }

            std::cout << "\n";

            // TODO: used for color offsets (instad of color ids)

        } else if (entry_size == 32) {
            std::vector<uint32_t> pml_lens;
            pml_lens.resize(mq_pml_lens_size);
            mls_file.read(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));

            for (int64_t i = mq_pml_lens_size - 1; i >= 0; i--) {
                std::cout << pml_lens[i] << " ";
            }

            std::cout << "\n";

            if (movi_options.is_classify()) {
                // Classification with 32 bits pmls works if the pml values are less than 2^16.
                std::vector<uint16_t> matching_lens_16(pml_lens.begin(), pml_lens.end());
                classifier.classify(read_name, matching_lens_16, movi_options);
            }
        } else {
            throw std::runtime_error("Invalid BPF entry size.");
        }

    }

    if (movi_options.is_classify()) {
        classifier.close_report_file();
    }
}

void build_rlbwt(MoviOptions& movi_options) {
    INFO_MSG("The run and len files are being built.");

    std::ifstream bwt_file(movi_options.get_bwt_file());
    if (!bwt_file.good()) {
        throw std::runtime_error(ERROR_MSG("[build_rlbwt] Failed to open the BWT file: " + movi_options.get_bwt_file()));
    }

    bwt_file.clear();
    bwt_file.seekg(0,std::ios_base::end);
    std::streampos end_pos = bwt_file.tellg();
    uint64_t length = static_cast<uint64_t>(end_pos);

    if (movi_options.is_verbose()) {
        INFO_MSG("end_pos: " + std::to_string(end_pos));
    }

    bwt_file.seekg(0);
    char current_char = bwt_file.get();
    char last_char = current_char;
    uint64_t r = 0;
    size_t len = 0;

    std::ofstream len_file(movi_options.get_bwt_file() + ".len", std::ios::out | std::ios::binary);
    if (!len_file.good()) {
        throw std::runtime_error(ERROR_MSG("[build_rlbwt] Failed to open the length file: " + movi_options.get_bwt_file() + ".len"));
    }

    std::ofstream heads_file(movi_options.get_bwt_file() + ".heads");
    if (!heads_file.good()) {

        throw std::runtime_error(ERROR_MSG("[build_rlbwt] Failed to open the heads file: " + movi_options.get_bwt_file() + ".heads"));
    }

    uint64_t i = 0;
    while (current_char != EOF) {
        if (i % 1000000 == 0 or i == length - 1 or i == 1) {
            print_progress_bar(i, length - 1, "Building the rlbwt", 1, 1);
        }

        if (current_char != last_char) {
            r += 1;
            // write output
            heads_file << last_char;
            len_file.write(reinterpret_cast<char*>(&len), 5);
            len = 0;
        }
        len += 1;
        last_char = current_char;
        current_char = bwt_file.get();

        i += 1;
    }
    PROGRESS_MSG("Successfully built the rlbwt.");
    INFO_MSG("n:\t" + std::to_string(length));
    INFO_MSG("r:\t" + std::to_string(r));
    INFO_MSG("n/r:\t" + std::to_string(static_cast<double>(length)/r));

    // write output
    heads_file << last_char;
    len_file.write(reinterpret_cast<char*>(&len), 5);

    heads_file.close();
    len_file.close();

    SUCCESS_MSG("The rlbwt is successfully stored at \n" + movi_options.get_bwt_file() + ".heads and \n" + movi_options.get_bwt_file() + ".len");
}

int main(int argc, char** argv) {

    try {

        MoviOptions movi_options;

        if (!parse_command(argc, argv, movi_options)) {
            return 1;
        }

        // If validate flags is set, just return 0, no execution is needed.
        if (movi_options.is_validate_flags()) {
            return 0;
        }

        if (movi_options.is_stdout()) {
            // Disable sync for faster I/O
            std::ios_base::sync_with_stdio(false);

            // Untie std::cin from std::cout for better performance
            std::cin.tie(nullptr);

            // Define a buffer of 1 MB
            constexpr size_t BUFFER_SIZE = 1024 * 1024; // 1 MB
            char buffer[BUFFER_SIZE];
            // Set the custom buffer for std::cout
            std::cout.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

            // Set custom buffer for stdout (printf uses stdout)
            if (setvbuf(stdout, buffer, _IOFBF, BUFFER_SIZE) != 0) {
                perror("Failed to set buffer for stdout");
                return 1;
            }
        }

        std::string command = movi_options.get_command();

        if (command == "build") {
            if (movi_options.use_separators()) {
                if (!SUPPORTS_SEPARATORS) {
                    // TODO: Fully support separators for large, split, and constant indexes
                    throw std::runtime_error(ERROR_MSG("[build] Separators are not supported for the " + program() + " index."));
                }
            }

            MoveStructure mv_(&movi_options, SPLIT_ARRAY, CONSTANT_INDEX);

            mv_.build();

            if (movi_options.is_verify()) {
                INFO_MSG("Verifying the LF_move results...");
                mv_.verify_lf_loop();
            }
            mv_.serialize();
            build_ftab(mv_, movi_options);
            SUCCESS_MSG("The Movi index is successfully stored at " + movi_options.get_index_dir());
            if (movi_options.is_output_ids()) {
                mv_.output_ids();
            }

            INFO_MSG("Generating the null statistics...");
            Classifier classifier;

            // generate pml null database
            movi_options.set_pml();
            movi_options.set_generate_null_reads(true);
            classifier.generate_null_statistics(mv_, movi_options);
            INFO_MSG("Successfully generated null statistics with PML");

            // generate zml null database
            movi_options.set_zml();
            movi_options.set_generate_null_reads(false); // do not regenerate the null reads
            classifier.generate_null_statistics(mv_, movi_options);
            INFO_MSG("Successfully generated null statistics with ZML");

            if (movi_options.is_color()) {
                color(mv_, movi_options);
            }

        } else if (command == "build-SA") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            mv_.find_sampled_SA_entries();
            mv_.serialize_sampled_SA();
            SUCCESS_MSG("Successfully stored sampled SA entries at " + movi_options.get_index_dir());
        } else if (command == "build-kmerbv") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            mv_.build_kmerbv(movi_options.get_build_kmerbv_ks());
            SUCCESS_MSG("Successfully built k-mer bitvector(s) at " + movi_options.get_index_dir());
        } else if (command == "color") {
            MoveStructure mv_(&movi_options);
            auto begin = std::chrono::system_clock::now();
            mv_.deserialize();
            auto end = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            TIMING_MSG(elapsed, "loading the index");

            color(mv_, movi_options);
        } else if (command == "query") {
            // Check if the input file exists
            if (movi_options.get_read_file() != "-" and !std::filesystem::exists(movi_options.get_read_file())) {
                throw std::runtime_error(ERROR_MSG("The input file " + movi_options.get_read_file() + " does not exist."));
            }

            MoveStructure mv_(&movi_options);

            auto begin = std::chrono::system_clock::now();
            mv_.deserialize();
            auto end = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            TIMING_MSG(elapsed, "loading the index");

            if (movi_options.is_get_sa_entries()) {
                begin = std::chrono::system_clock::now();
                mv_.deserialize_sampled_SA();
                end = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
                TIMING_MSG(elapsed, "loading the sampled SA entries");
            }

            if (movi_options.is_kmer_bv()) {
                begin = std::chrono::system_clock::now();
                mv_.load_kmerbv(movi_options.get_k());
                end = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
                TIMING_MSG(elapsed, "loading the k-mer bitvector");
            }

            if (movi_options.is_multi_classify()) {
                load_color_table(mv_, movi_options);

                if (movi_options.is_report_colors() or movi_options.is_report_color_ids()) {
                    // Copmuter color ids for the output
                    mv_.compute_color_ids_from_flat();
                }
            }

            query(mv_, movi_options);

            // quick_exit() skips static destructors, so the std::cout buffer (custom
            // buffer set above) is NOT flushed automatically -- trailing output such
            // as the end-of-run k-mer stats / SSHash report would be lost. Flush
            // explicitly before exiting.
            std::cout.flush();

            // Avoid taking too long for dealloction of large data structures at the end of the program
            std::quick_exit(0);

        } else if (command == "view") {
            view(movi_options);
        } else if (command == "rlbwt") {
            build_rlbwt(movi_options);
        } else if (command == "color-move-rows") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();

            load_color_table(mv_, movi_options);

            mv_.add_colors_to_rlbwt();

            mv_.serialize();

        } else if (command == "LF") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            INFO_MSG("The Movi index is read from the file successfully.");
            auto begin = std::chrono::system_clock::now();
            if (movi_options.get_LF_type() == "sequential") {
                mv_.sequential_lf();
            } else if (movi_options.get_LF_type() == "random") {
                mv_.random_lf();
            } else if (movi_options.get_LF_type() == "reconstruct") {
                mv_.reconstruct_lf();
            }
            auto end = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            TIMING_MSG(elapsed, "LF-mapping for all the BWT characters in the" + movi_options.get_LF_type() + " order");
        } else if (command == "inspect") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            mv_.print_stats();
            if (movi_options.is_flat_color_vectors()) {
                std::string fname = movi_options.get_index_dir() + "/doc_sets.bin";
                mv_.deserialize_doc_sets(fname);
                mv_.load_document_info();
                INFO_MSG("The color table is read successfully.");
                mv_.flat_and_serialize_colors_vectors();
            }
            // mv_.compute_run_lcs();
            // mv_.analyze_rows();
        } else if (command == "ftab") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            build_ftab(mv_, movi_options);
        } else if (command == "null") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            Classifier classifier;
            classifier.generate_null_statistics(mv_, movi_options);
        } else {
            const std::string message = "Invalid action: \"" + command + "\"";
            throw std::runtime_error(message);
        }

        return 0;

    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
