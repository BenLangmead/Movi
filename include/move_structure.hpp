#ifndef MOVE_STRUCTURE_HPP
#define MOVE_STRUCTURE_HPP

#include <fstream>
#include <cstdint>
#include <stdio.h>
#include <chrono>
#include <cstddef>
#include <unistd.h>
#include <sys/stat.h>
#include <vector>
#include <omp.h>          // per-thread k-mer statistics slots
#include <unordered_map>
#include <map>
#include <sstream>
#include <filesystem>
#include <span>
#include <sys/stat.h>


#include <span>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "fastcluster.h"

#include "sdsl_wrapper.hpp"

#include "movi_options.hpp"
#include "move_row.hpp"
#include "move_row_colored.hpp"
#include "move_query.hpp"
#include "query_kmer.hpp"
#include "query_mem.hpp"
#include "utils.hpp"
#include "move_intervals.hpp"
#include "doc_set.hpp"
#include <coroutine>

class Classifier;

// Forward declaration for SharedFastqReader (defined in coroutine_processor.cpp)
class SharedFastqReader;

struct ThresholdsRow {
    uint16_t values[4];
};

class MoveStructure {
    public:
        MoveStructure(MoviOptions* movi_options_);
        MoveStructure(MoviOptions* movi_options_, uint16_t nt_splitting_, bool constant_);

        std::vector<MoveRow> get_rlbwt();
        char get_char(uint64_t idx);
        uint64_t get_n(uint64_t idx);
        uint64_t get_offset(uint64_t idx);
        uint64_t get_id(uint64_t idx);
        bool use_separator();
        // TODO: The following is useful for mmaping, there is a slowdown though
        // MoveRow& get_move_row(uint64_t idx);

#if USE_THRESHOLDS
        uint64_t get_thresholds(uint64_t idx, uint32_t alphabet_index);
#endif

#if SPLIT_THRESHOLDS_FALSE
        uint16_t get_rlbwt_thresholds(uint64_t idx, uint16_t i);
        void set_rlbwt_thresholds(uint64_t idx, uint16_t i, uint16_t value);
#endif

        uint64_t get_SA_entries(uint64_t idx, uint64_t offset);

        uint64_t LF(uint64_t row_number, uint64_t alphabet_index);
        uint64_t LF_heads(uint64_t run_number, uint64_t alphabet_index);
        // The 3rd argument of LF_move is used in the latency_hiding_tally mode
        uint16_t LF_move(uint64_t& pointer, uint64_t& i, uint64_t id = std::numeric_limits<uint64_t>::max());
        uint64_t fast_forward(uint64_t& offset, uint64_t index, uint64_t x);

        char illegal_char_substitute();
        // Applies the --ignore-illegal-chars substitution across the whole read, once.
        void substitute_illegal_chars(MoveQuery& mq);
        // Walks rc(read) once, leaving for every present k-mer a row inside its reverse
        // complement's BWT interval.
        void prepare_rc_kmer_rows(MoveQuery& mq);
        uint64_t rc_kmer_rows_walk(MoveQuery& mq, int32_t& pos_on_r, size_t read_len);
        // Run-local representation of B_k for the id query. An id is rank_{B_k}(lb), a
        // global prefix count, so unlike the count structure's predecessor and successor
        // a purely local record does not suffice: the decomposition carries a prefix term
        // as well as a within-run term,
        //   id = (marks in runs < run) + (marks in run `run` at offset <= off) - 1.
        // A mark on a run head is a bit in kmerbv_idhd; a mark at a positive offset is an
        // entry in kmerbv_idgioff, delimited per run by kmerbv_idsep. Every access is
        // keyed by `run`, which the search has just touched, so no absolute BWT row and
        // no per-run position table enter the id path.
        bool kmerbv_use_runlocal = false;
        sdsl::bit_vector kmerbv_idhd;   // idhd[run] = 1 iff the run head carries a B_k mark
        sdsl::rank_support_v<> kmerbv_idhd_rank;
        // Run separators over the off-head marks: run i contributes a 1 followed by one 0
        // per off-head mark it holds. select_1(i + 1) lands at i plus the marks in earlier
        // runs, giving that prefix by subtraction, and the zeros following it give this
        // run's own count. Keeping the prefix here leaves the head rank and this select
        // independent, so the two issue in parallel.
        sdsl::bit_vector kmerbv_idsep;
        sdsl::select_support_mcl<1> kmerbv_idsep_sel;
        // Off-head offsets, delta-coded within each run and packed in a dac_vector, which
        // reads any entry in O(1). A run's slice is short, so it is decoded by a scan.
        sdsl::dac_vector<2> kmerbv_idgioff;

        // Number of 0s immediately following position p in kmerbv_idsep, which is the
        // number of off-head marks held by the run whose separator sits at p.
        inline uint64_t kmerbv_idsep_zeros_after(uint64_t p) const {
            uint64_t rest = kmerbv_idsep.size() - (p + 1);
            if (rest == 0) return 0;
            uint8_t w = (rest >= 64) ? 64 : static_cast<uint8_t>(rest);
            uint64_t word = kmerbv_idsep.get_int(p + 1, w);
            if (word != 0) return sdsl::bits::lo(word);
            if (rest < 64) return rest;
            uint64_t z = 64, q = p + 65;
            while (q < kmerbv_idsep.size() and kmerbv_idsep[q] == 0) { ++z; ++q; }
            return z;
        }

        // Number of B_k marks at rows up to and including (run, offset), which is
        // rank_{B_k}(lb + 1) for lb = all_p[run] + offset.
        inline uint64_t kmerbv_rank_runlocal(uint64_t run, uint64_t offset) {
            uint64_t p = kmerbv_idsep_sel(run + 1);
            uint64_t a = p - run;   // off-head marks lying in runs before `run`
            // A head mark sits at offset 0, so it always counts when the run matches.
            uint64_t cnt = kmerbv_idhd_rank(run) + a + (kmerbv_idhd[run] ? 1 : 0);
            uint64_t c = kmerbv_idsep_zeros_after(p);
            if (c > 0) {
                uint64_t cur = kmerbv_idgioff[a];
                if (cur <= offset) {
                    cnt += 1;
                    for (uint64_t i = 1; i < c; ++i) {
                        cur += kmerbv_idgioff[a + i];
                        if (cur > offset) break;
                        cnt += 1;
                    }
                }
            }
            return cnt;
        }
        bool check_alphabet(char& c);
        uint32_t compute_index(char row_char, char lookup_char);

        void analyze_rows();
        void print_stats();
        void print_basic_index_info();

/***************************************************************************/
/****  Beginning of functions implemented in move_structure_search.cpp  ****/
        MoveInterval try_ftab(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len, size_t ftab_k, bool rc = false);
        bool look_ahead_ftab(MoveQuery& mq, uint32_t pos_on_r, int32_t& step);
        bool look_ahead_backward_search(MoveQuery& mq, uint32_t pos_on_r, int32_t step);

        uint64_t query_backward_search(MoveQuery& mq, int32_t& pos_on_r);
        bool backward_search_step(char c, MoveInterval& interval);
        bool forward_search_step(char c, MoveInterval& rc_interval);
        uint64_t backward_search_step(std::string& R, int32_t& pos_on_r, MoveInterval& interval);
        MoveInterval backward_search(std::string& R, int32_t& pos_on_r, MoveInterval interval, int32_t max_length);
        MoveInterval initialize_backward_search(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len, bool rc = false);
        void update_interval(MoveInterval& interval, char next_char);

        bool extend_bidirectional(char c_, MoveInterval& fw_interval, MoveInterval& rc_interval);
        bool extend_left(char c, MoveBiInterval& bi_interval);
        bool extend_right(char c, MoveBiInterval& bi_interval);
        MoveBiInterval backward_search_bidirectional(std::string& R, int32_t& pos_on_r, MoveBiInterval interval, int32_t max_length);
        MoveBiInterval initialize_bidirectional_search(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len);
/*******  End of functions implemented in move_structure_search.cpp  *******/
/***************************************************************************/

/***************************************************************************/
/****  Beginning of functions implemented in move_structure_build.cpp   ****/
        void build();
        uint64_t compute_length_from_bwt();
        void initialize_bits();
        void fill_bits_by_thresholds();
        void detect_move_row_boundaries();
        void add_detected_run(uint64_t scanned_bwt_length, uint64_t run_char, uint16_t& run_length);
        void find_run_heads_information();
        void build_move_rows();
        void find_base_interval_data();
        void build_alphabet(std::vector<uint64_t>& all_possible_chars);

#if USE_THRESHOLDS
        void compute_thresholds();
        // Helper for threshold calculation
        void debug_threshold_calculation(uint64_t i, uint64_t thr_i, uint64_t j, char rlbwt_c,
                                         const std::vector<uint64_t>& alphabet_thresholds);
        void set_threshold_for_one_character(uint64_t i, uint16_t threshold_index, uint64_t value,
                                             uint16_t value_split, std::vector<uint64_t>& current_thresholds);
#endif

#if USE_NEXT_POINTERS
        void compute_nexts();
#endif

#if BLOCKED_MODES
        void compute_blocked_ids(std::vector<uint64_t>& raw_ids);
#endif

        void build_ftab();

        // Finds SA entries of all rows in BWT.
        void find_sampled_SA_entries();
/*******  End of functions implemented in move_structure_build.cpp   *******/
/***************************************************************************/

/***************************************************************************/
/****  Beginning of functions implemented in move_structure_query.cpp   ****/
        void verify_lfs();
        void verify_lf_loop();
        void sequential_lf();
        void random_lf();
        void reconstruct_lf();

        uint64_t query_pml(MoveQuery& mq);
        uint64_t query_zml(MoveQuery& mq);

        // C++20 coroutine versions of the queries (latency-hiding style (c); see coroutine_processor.cpp).
        // All three share one return type and the scheduler scaffolding.
        struct coroutine_task;
        coroutine_task query_pml_coroutine(
            SharedFastqReader& reader,
            std::coroutine_handle<>& my_handle_storage,
            int coroutine_id);

        // Coroutine MEM query (bidirectional); mirrors query_mems / query_mem_bml.
        coroutine_task query_mem_coroutine(
            SharedFastqReader& reader,
            std::coroutine_handle<>& my_handle_storage,
            int coroutine_id);

        // Coroutine k-mer query (--kmer / --kmer-count); mirrors query_all_kmers.
        coroutine_task query_kmer_coroutine(
            SharedFastqReader& reader,
            std::coroutine_handle<>& my_handle_storage,
            int coroutine_id);

        uint64_t reposition_up(uint64_t idx, char c, uint64_t& scan_count);
        uint64_t reposition_down(uint64_t idx, char c, uint64_t& scan_count);
        bool reposition_randomly(uint64_t& idx, uint64_t& offset, char r_char, uint64_t& scan_count);

        uint64_t query_mems(MoveQuery& mq);
        bool query_mem_bml(MoveQuery& mq, int32_t& pos_on_r,int32_t& min_mem_length, std::string& query_seq, size_t& ftab_k);
        uint64_t query_all_mems(MoveQuery& mq);

#if USE_THRESHOLDS
        bool reposition_thresholds(uint64_t& idx, uint64_t offset, char r_char, uint64_t& scan_count);
        // Helper for threshold calculation
        void handle_reposition_up(uint64_t& idx, uint64_t saved_idx, char r_char, uint16_t next_up, uint64_t& scan_count);
        void handle_reposition_down(uint64_t& idx, uint64_t saved_idx, char r_char, uint16_t next_down, uint64_t& scan_count);
#endif
/*******  End of functions implemented in move_structure_query.cpp  *******/
/***************************************************************************/

/***************************************************************************/
/**********  Beginning of functions implemented in query_kmer.cpp  ***********/
        void query_all_kmers(MoveQuery& mq, bool kmer_counts = false);
        uint64_t query_kmers_from(MoveQuery& mq, int32_t& pos_on_r, bool single = false,
                                  MoveInterval* interval_out = nullptr);
        uint64_t query_kmers_count_bv(MoveQuery& mq, int32_t& pos_on_r);
        uint64_t query_kmers_id_bv(MoveQuery& mq, int32_t& pos_on_r);
/*******  End of functions implemented in query_kmer.cpp  ***********/
/***************************************************************************/

/***************************************************************************/
/****  Beginning of functions implemented in move_structure_color.cpp  *****/
        int get_num_docs() { return num_docs; }
        int get_num_species() { return num_species; }

        // Fill the run offsets array (used for building colors among other things)
        void fill_run_offsets();
        // Builds document sets for each run in rlbwt.
        void build_doc_sets();
        uint32_t hash_collapse(std::unordered_map<DocSet, uint32_t> &keep_set, DocSet &bv);
        void build_tree_doc_sets();
        void build_doc_set_similarities();
        void compress_doc_sets();
        void compute_color_ids_from_flat();
        // Finds documents corresponding to rows in BWT.
        void build_doc_pats();
        // Initialize classify counts
        void initialize_classify_cnts();
        void set_classifier(Classifier *cl) { classifier = cl; }
        void set_output_files(OutputFiles *of) { output_files = of; }

        // Document tree functions
        bool is_ancestor(uint16_t x, uint16_t y);
        uint16_t LCA(uint16_t x, uint16_t y);
        void dfs_times(uint16_t cur, uint16_t &t);
        void compute_run_lcs();

        // The following is just used to test how storing colors in the rlbwt affects the performance
        void add_colors_to_rlbwt();
/******** End of functions implemented in move_structure_color.cpp *********/
/***************************************************************************/

/***************************************************************************/
/*****  Beginning of functions implemented in move_structure_io.cpp  *******/

        // Writes frequencies of document sets to file.
        void write_doc_set_freqs(std::string fname);

        void flat_and_serialize_colors_vectors();
        void deserialize_doc_sets_flat();

        void serialize_doc_pats(std::string fname);
        void deserialize_doc_pats(std::string fname);
        void serialize_doc_sets(std::string fname);
        void deserialize_doc_sets(std::string fname);
        void load_document_info();

        // The following two methods are not implemented
        void serialize_doc_rows();
        void deserialize_doc_rows();

        std::ifstream open_index_read();
        std::ofstream open_index_write();
        void read_index_header(std::ifstream& fin);
        void write_index_header(std::ofstream& fout);
        void write_basic_index_data(std::ofstream& fout);
        void read_basic_index_data(std::ifstream& fin);
        void write_overflow_tables(std::ofstream& fout);
        void read_overflow_tables(std::ifstream& fin);
        void write_counts_data(std::ofstream& fout);
        void read_counts_data(std::ifstream& fin);
        void write_main_table(std::ofstream& fout);
        void read_main_table(std::ifstream& fin, std::streamoff rlbwt_offset);
        void write_separators_thresholds(std::ofstream& fout);
        void read_separators_thresholds(std::ifstream& fin);
        void serialize();
        void deserialize();
        void output_ids();

#if BLOCKED_MODES
        void write_id_blocks(std::ofstream& fout);
        void read_id_blocks(std::ifstream& fin);
#endif

#if TALLY_MODES
        void write_tally_table(std::ofstream& fout);
        void read_tally_table(std::ifstream& fin);
#endif

        void serialize_sampled_SA();
        void deserialize_sampled_SA();

        void write_ftab();
        void read_ftab();

        void build_kmerbv(const std::vector<uint32_t>& ks);
        void load_kmerbv(uint32_t k);
        void rebuild_all_p_if_needed();
/*********  End of functions implemented in move_structure_io.cpp  *********/
/***************************************************************************/

	    KmerStatistics kmer_stats;

        // Per-thread k-mer tallies, merged into kmer_stats by merge_kmer_stats() once a
        // query run is over. The k-mer queries update these counters for every k-mer, so
        // a single shared tally puts every thread on the same few cache lines: even as
        // atomics they were what limited how the query scaled, not the query work.
        // Each slot is padded to its own cache line so neighbouring threads stay apart.
        struct alignas(64) PaddedKmerStatistics {
            KmerStatistics stats;
        };
        // Sized at construction so thread_kmer_stats() is always valid; prepare and
        // merge only zero it.
        std::vector<PaddedKmerStatistics> kmer_stats_per_thread =
            std::vector<PaddedKmerStatistics>(omp_get_max_threads());

        // This thread's tally. Callers take the reference once per query call rather
        // than looking it up for each counter update. The vector is always sized to
        // omp_get_max_threads(), so there is no empty case to fall back on: a fallback
        // to the shared tally would reintroduce, silently, the contention and the race
        // this exists to remove.
        KmerStatistics& thread_kmer_stats() {
            return kmer_stats_per_thread[omp_get_thread_num()].stats;
        }

        // Zero the per-thread tallies. Call outside any parallel region, before a query.
        void prepare_kmer_stats() {
            kmer_stats_per_thread.assign(omp_get_max_threads(), PaddedKmerStatistics());
        }

        // Fold the per-thread tallies into kmer_stats and zero them again. Call after
        // the parallel region and before anything reports the statistics.
        void merge_kmer_stats() {
            for (auto& p : kmer_stats_per_thread) {
                kmer_stats.merge(p.stats);
            }
            kmer_stats_per_thread.assign(kmer_stats_per_thread.size(), PaddedKmerStatistics());
        }

        friend class ReadProcessor;
    private:
        // Reference to output files for writing results
        OutputFiles* output_files;

        // The BWT file used for building the index
        std::ifstream bwt_file;

        // Sorted vector of the start offsets of each document.  
        std::vector<uint64_t> doc_offsets;
        // Species ID for each document
        std::vector<uint32_t> doc_ids;
        // Map from taxon id to compressed species index
        std::map<uint32_t, uint32_t> taxon_id_compress;
        // Compressed species index to taxon id
        std::vector<uint32_t> to_taxon_id;
        // log length of each species
        std::vector<double> log_lens;
        uint32_t num_docs;
        uint32_t num_species;
        uint64_t num_colors = 0;

        // Offset of run heads in the rlbwt.
        std::vector<uint64_t> run_offsets;

        // Document sets.
        std::vector<uint16_t> flat_colors;
        std::unordered_map<uint64_t, uint32_t> color_offset_to_id;
        std::vector<std::vector<uint16_t>> unique_doc_sets;
        std::vector<uint32_t> doc_set_inds;
        std::vector<MoveTally> doc_set_flat_inds;
        sdsl::bit_vector compressed;

        // Tree over documents
        std::vector<std::vector<uint16_t>> tree;
        std::vector<std::vector<uint16_t>> tree_doc_sets;
        std::vector<std::vector<uint16_t>> bin_lift;
        std::vector<uint16_t> t_in;
        std::vector<uint16_t> t_out;
        
        // Count of how much each doc set appears (by ID).
        std::vector<uint64_t> doc_set_cnts;

        // Document patterns (species that each row in BWT belongs to).
        std::vector<uint16_t> doc_pats;

        // Counts for multi classification
        std::vector<uint32_t> classify_cnts;
        std::vector<double> doc_scores;
        const double log4 = log(4);

        // Classifier object for binary classification
        Classifier *classifier;
    
        // Basic index configurations
        MoviOptions* movi_options;
        uint16_t nt_splitting;
        bool constant;

        // Vector of sampled SA entries. For experiment purposes.
        std::vector<uint64_t> sampled_SA_entries;

        // Basic index characteristics
        std::string bwt_string;
        uint64_t length;
        uint64_t r;
        uint64_t original_r;

        void compute_number_of_build_steps();
        uint64_t total_build_steps = 0;
        uint64_t current_build_step = 0;

        // The explicit values for the end bwt row
        uint64_t end_bwt_idx;
        uint64_t end_bwt_idx_thresholds[4];
        uint64_t end_bwt_idx_next_up[4];
        uint64_t end_bwt_idx_next_down[4];

        // The explicit thresholds for the separators
        std::vector<ThresholdsRow> separators_thresholds;
        std::unordered_map<uint64_t, uint64_t> separators_thresholds_map;

        // Map from 2bit encoded character to the actual character
        // Example: alphabet[0] -> A, alphabet[1] -> C
        std::vector<unsigned char> alphabet;
        // Number of each character
        std::vector<uint64_t> counts;
        // Map from the character to the index of the character
        // Example: alphamap[A] -> 0, alphamap[C] -> 1
        std::vector<uint64_t> alphamap;

        // The move structure rows
        std::vector<MoveRow> rlbwt;
        std::span<MoveRow> rlbwt_view;
        std::vector<MoveRowColored> rlbwt_colored;
#if TALLY_MODES
        uint32_t tally_checkpoints;
        std::vector<std::vector<MoveTally>> tally_ids;
#endif

#if BLOCKED_MODES
        std::vector<std::vector<uint32_t>> id_blocks;
        uint64_t block_size = BLOCK_SIZE;
#endif

        // auxilary datastructures for the length, offset and thresholds overflow
        std::vector<uint64_t> n_overflow;
        std::vector<uint64_t> offset_overflow;
        std::vector<std::vector<uint64_t> > thresholds_overflow;

        // Values used for starting the beckawrd search
        std::vector<uint64_t> first_runs;
        std::vector<uint64_t> first_offsets;
        std::vector<uint64_t> last_runs;
        std::vector<uint64_t> last_offsets;
        std::vector<std::vector<MoveInterval> > ftabs;
        std::vector<MoveInterval> ftab;

        // Used for gathering statistics
        uint64_t no_ftab;
        uint64_t all_initializations;
        std::unordered_map<uint32_t, uint32_t> repositions;
        std::unordered_map<uint32_t, uint32_t> ff_counts;
        std::unordered_map<uint64_t, uint64_t> run_lengths;

        // The following are only used for construction, not stored in the index
        std::vector<uint64_t> thresholds;
        std::vector<uint64_t> all_p;
        std::vector<char> heads;
        std::vector<char> original_run_heads;
        std::vector<std::unique_ptr<sdsl::bit_vector> > occs;
        std::vector<std::unique_ptr<sdsl::rank_support_v<> > > occs_rank;
        std::vector<uint64_t> heads_rank;
        std::vector<uint64_t> lens;
        std::vector<uint32_t> original_lens;

        sdsl::bit_vector bits;
        sdsl::rank_support_v<> rbits;

        // K-mer boundary bitvector B_k and rank support (for the MPHF ID; rank(lb)).
        // Populated by build_kmerbv() / load_kmerbv(); empty unless those are called.
        sdsl::bit_vector kmerbv;
        sdsl::rank_support_v<> kmerbv_rank;

        // Sparse (Elias-Fano) representation of B_k: ~4x smaller on disk than the
        // dense bit_vector+rank when B_k is sparse (canonical k-mers mark only a few
        // percent of the n BWT positions), making the MPHF-id add-on competitive
        // with SSHash. Built alongside the dense bv; selected at load via MOVI_ID_BV=sd.
        // rank/select are built in. kmerbv_rank1() routes the id rank to the chosen rep.
        bool kmerbv_use_sd = false;
        sdsl::sd_vector<> kmerbv_sd;
        sdsl::sd_vector<>::rank_1_type kmerbv_sd_rank;
        sdsl::sd_vector<>::select_1_type kmerbv_sd_sel;
        // Sampled per-run absolute-position checkpoints for the MPHF id query. Both
        // the sd and dense B_k reps mark ABSOLUTE BWT rows, so id = rank(lb+1)-1 with
        // lb = all_p[run] + offset. Materializing the full per-run all_p costs
        // ~r*log2(n) bits (~ a whole SSHash index); instead we store one checkpoint
        // every kmerbv_allp_S runs and recover an exact run-head position by summing
        // the move structure's run lengths (get_n, already in RAM) from the nearest
        // checkpoint -- a bounded O(kmerbv_allp_S) walk, paid once per present k-mer
        // (not per LF step). Trades a sampled-position RSS factor for that walk.
        sdsl::int_vector<> kmerbv_allp_ckpt;   // all_p at runs 0, S, 2S, ... (sampled)
        uint64_t kmerbv_allp_S = 0;            // checkpoint stride in runs; 0 => not built
        // Exact absolute row of run head `run` = nearest checkpoint + sum of the
        // intervening run lengths.
        inline uint64_t reconstruct_allp(uint64_t run) {
            uint64_t base = run - (run % kmerbv_allp_S);
            uint64_t p = kmerbv_allp_ckpt[run / kmerbv_allp_S];
            for (uint64_t j = base; j < run; ++j) p += get_n(j);
            return p;
        }
        // MPHF id of the k-mer whose interval starts at (run, offset): id =
        // rank(lb+1)-1 with lb = all_p[run]+offset. Holds for an exact OR a subset
        // interval (the k-mer's mark is the only 1 in its group). The run-local
        // representation answers the same rank without forming lb at all.
        inline uint64_t kmerbv_id(uint64_t run, uint64_t offset) {
            if (kmerbv_use_runlocal) return kmerbv_rank_runlocal(run, offset) - 1;
            uint64_t lb = reconstruct_allp(run) + offset;
            return kmerbv_use_sd ? (kmerbv_sd_rank(lb + 1) - 1)
                                 : (kmerbv_rank(lb + 1) - 1);
        }
        // From kmerbv.<k>.meta: whether B_k marks only CANONICAL k-mers (ids in
        // [0,num_kmers), SSHash-compatible) and the num_kmers (= ones in B_k).
        bool kmerbv_is_canonical = false;
        uint64_t kmerbv_num_kmers = 0;
        // Full backward search of a k-mer string; returns its BWT interval (empty if
        // absent). Used by the canonical-id lookup to locate rc(x)'s interval.
        MoveInterval search_kmer_interval(std::string s);  // by value: check_alphabet takes char&

        // Run-local (all_p-free) k-mer count structure. A k-mer's count is
        // resolved from Movi's run lengths plus a per-run record of which run heads
        // are group-starts (ex) and where the interior group-starts fall (hi +
        // gioff/gimark) -- no n-space count bitvector and no absolute-position
        // table. (Earlier revisions kept several interchangeable representations of
        // an n-space count bitvector -- plain bit_vector, sd_vector, rrr_vector, and
        // a run-heads + Gi + E hybrid -- to map the space/speed frontier; this
        // run-local form was both smallest and fastest, so it is the only one kept.)
        sdsl::bit_vector kmerbv_ex;  // ex[run] = 1 iff that run head is not a group-start
        sdsl::bit_vector kmerbv_hi;  // hi[run] = 1 iff that run contains an interior group-start
        // Per-run (all_p-free) interior group-starts. The per-entry offsets are
        // DELTA-TRANSFORMED then variable-length coded: kmerbv_gioff[i] holds the
        // first interior's absolute offset where gimark[i]==1 (run start), and the
        // delta to the previous interior otherwise. Stored as a dac_vector<2>
        // (Elias-style direct-access code, O(1) random read) -- ~5.7 MB vs 14.1 MB
        // for the fixed-width offsets on ecoli10 k=31, since within-run offsets are
        // small and increasing (mean delta ~3, p99 ~27). gimark + hi's rank map a
        // run to its [a,b) slice; the slice is decoded sequentially, accumulating
        // deltas (avg ~1.7 entries/run, so scanning beats binary search).
        sdsl::rank_support_v<> kmerbv_hi_rank;
        uint64_t kmerbv_hi_ones = 0;
        sdsl::dac_vector<2> kmerbv_gioff;   // delta-transformed, variable-length
        sdsl::bit_vector kmerbv_gimark;
        sdsl::select_support_mcl<1> kmerbv_gimark_sel;

        // [a,b) range of run's interior entries within kmerbv_gioff (empty if none).
        inline void gi_range(uint64_t run, uint64_t& a, uint64_t& b) {
            if (!kmerbv_hi[run]) { a = b = 0; return; }
            uint64_t rho = kmerbv_hi_rank(run);                       // 0-indexed hi-run rank
            a = kmerbv_gimark_sel(rho + 1);
            b = (rho + 1 < kmerbv_hi_ones) ? kmerbv_gimark_sel(rho + 2) : kmerbv_gioff.size();
        }
        static const uint64_t GI_NONE = std::numeric_limits<uint64_t>::max();
        // The slice is delta-coded: entry a is the absolute first offset, each
        // later entry adds its delta. Offsets are strictly increasing, so a single
        // forward scan resolves pred/succ/first/last.
        inline uint64_t gi_pred_off(uint64_t run, uint64_t off) {     // largest interior <= off
            uint64_t a, b; gi_range(run, a, b);
            if (a == b) return GI_NONE;
            uint64_t cur = kmerbv_gioff[a];
            if (cur > off) return GI_NONE;
            uint64_t best = cur;
            for (uint64_t i = a + 1; i < b; ++i) { cur += kmerbv_gioff[i]; if (cur > off) break; best = cur; }
            return best;
        }
        inline uint64_t gi_succ_off(uint64_t run, uint64_t off) {     // smallest interior > off
            uint64_t a, b; gi_range(run, a, b);
            if (a == b) return GI_NONE;
            uint64_t cur = kmerbv_gioff[a];
            if (cur > off) return cur;
            for (uint64_t i = a + 1; i < b; ++i) { cur += kmerbv_gioff[i]; if (cur > off) return cur; }
            return GI_NONE;
        }
        inline uint64_t gi_first_off(uint64_t run) {
            uint64_t a, b; gi_range(run, a, b);
            return a == b ? GI_NONE : (uint64_t)kmerbv_gioff[a];      // first entry is absolute
        }
        inline uint64_t gi_last_off(uint64_t run) {
            uint64_t a, b; gi_range(run, a, b);
            if (a == b) return GI_NONE;
            uint64_t cur = kmerbv_gioff[a];
            for (uint64_t i = a + 1; i < b; ++i) cur += kmerbv_gioff[i];
            return cur;
        }

        // Count of the k-mer whose interval is `iv` (any sub-range of its true
        // group works): predecessor of lb and successor of rb bracket the k-mer's
        // single group. pred(p)=select(rank(p+1)); succ(p)=select(rank(p+1)+1).
        inline uint64_t kmer_count_from_bv(const MoveInterval& iv) {
            // Fully run-local count: no all_p, no absolute positions.
            // count = down + iv_size + up, all from run lengths + per-run offsets.
            uint64_t off_s = iv.offset_start, off_e = iv.offset_end;
            uint64_t iv_size;
            if (iv.run_start == iv.run_end) {
                iv_size = off_e - off_s + 1;
            } else {
                iv_size = (get_n(iv.run_start) - off_s) + (off_e + 1);
                for (uint64_t j = iv.run_start + 1; j < iv.run_end; ++j) iv_size += get_n(j);
            }
            // down: expand from lb' (run_start, off_s) to the group start
            uint64_t down;
            uint64_t gd = gi_pred_off(iv.run_start, off_s);
            if (gd != GI_NONE) {
                down = off_s - gd;                       // boundary is an interior in run_start
            } else if (!kmerbv_ex[iv.run_start]) {
                down = off_s;                            // boundary is run_start head
            } else {                                     // excepted: walk into lower runs
                down = off_s;
                uint64_t j = iv.run_start;
                while (true) {
                    --j;
                    uint64_t gl = gi_last_off(j);
                    if (gl != GI_NONE) { down += get_n(j) - gl; break; }
                    down += get_n(j);
                    if (!kmerbv_ex[j]) break;            // run j head is the boundary
                }
            }
            // up: expand from rb' (run_end, off_e) to the group end
            uint64_t up;
            uint64_t gu = gi_succ_off(iv.run_end, off_e);
            if (gu != GI_NONE) {
                up = gu - 1 - off_e;                     // boundary is an interior in run_end
            } else if (iv.run_end + 1 >= r || !kmerbv_ex[iv.run_end + 1]) {
                up = get_n(iv.run_end) - 1 - off_e;      // next run head is the boundary
            } else {                                     // excepted: walk into upper runs
                up = get_n(iv.run_end) - 1 - off_e;
                uint64_t j = iv.run_end + 1;
                while (j < r) {
                    if (!kmerbv_ex[j]) break;            // run j head (offset 0) is the boundary
                    uint64_t gf = gi_first_off(j);       // else nearest is this run's first interior
                    if (gf != GI_NONE) { up += gf; break; }
                    up += get_n(j);
                    ++j;
                }
            }
            return down + iv_size + up;
        }
};

#endif