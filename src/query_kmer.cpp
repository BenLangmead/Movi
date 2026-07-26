#include <sys/stat.h>
#include <cstdlib>

#include "move_structure.hpp"
#include "utils.hpp"

// IMPORTANT: The index has to be built on a reference with separators to get the correct kmer counts
// IMPORTANT: For the bidirectional search, also the index has to be built with the separators


uint64_t MoveStructure::query_kmers_from(MoveQuery& mq, int32_t& pos_on_r, bool single,
                                          MoveInterval* interval_out) {
    KmerStatistics& ks = thread_kmer_stats();
    size_t ftab_k = movi_options->get_ftab_k();
    size_t k = movi_options->get_k();
    auto& query_seq = mq.query();
    int32_t pos_on_r_saved = pos_on_r;

    // An alternative strategy to look ahead for possible skipping
    // int32_t step = 0;
    // if (ftab_k > 1 and !look_ahead_ftab(mq, pos_on_r, step)) {
    //     kmer_stats.look_ahead_skipped += k - ftab_k - step;
    //     pos_on_r = pos_on_r - k + ftab_k + step - 1;
    //     pos_on_r_saved = pos_on_r;
    // }

    uint64_t match_len = 0;
    MoveInterval initial_interval;
    do {
        initial_interval = initialize_backward_search(mq, pos_on_r, match_len);
        if (match_len == 0 and ftab_k > 1) {
                        ks.initialize_skipped += 1;
            pos_on_r -= 1;
            pos_on_r_saved = pos_on_r;
        }
    } while (match_len == 0 and pos_on_r >= k - 1 and ftab_k > 1);

    // I want to check how much slower it gets if turn off the positive skip:
    auto backward_search_result = backward_search(query_seq, pos_on_r, initial_interval, single ? k - match_len - 2 : std::numeric_limits<int32_t>::max());

    if (backward_search_result.is_empty()) {
        // We get here when there is an illegal character at pos_on_r, just skip the current position
        pos_on_r = pos_on_r_saved - 1;
                ks.backward_search_empty += 1;
        return 0;
    } else {
        if (pos_on_r_saved - pos_on_r >= k - 1) {
            // At leat one kmer was found, update the postion and return the count
            uint64_t kmers_found = pos_on_r_saved - pos_on_r - k + 2;
                        ks.positive_skipped += kmers_found - 1;

            if (interval_out) *interval_out = backward_search_result;

            pos_on_r = pos_on_r + k - 2;
	        return kmers_found;
        } else {
            // No kmer was found, update the postion
                        ks.backward_search_failed += 1;
            pos_on_r = pos_on_r_saved - 1;
            return 0;
        }
    }
}


void MoveStructure::query_all_kmers(MoveQuery& mq, bool kmer_counts) {
    KmerStatistics& ks = thread_kmer_stats();
    size_t ftab_k = movi_options->get_ftab_k();
    size_t k = movi_options->get_k();
    auto& query_seq = mq.query();
    int32_t pos_on_r = query_seq.length() - 1;

    // To handle a special case for k equal to 1
    if (k == 1) {
        uint64_t kmers_found = 0;
        while (pos_on_r >= 0) {
            kmers_found += check_alphabet(query_seq[pos_on_r]) ? 1 : 0;
            pos_on_r -= 1;
        }
                ks.positive_kmers += kmers_found;
        return;
    }

    while (!check_alphabet(query_seq[pos_on_r])) {
        pos_on_r -= 1; // Find the first position where the character is legal
    }


    int32_t step = k/3;
    // The look-ahead ("fishing") subproblem has length k - step, which must stay >= ftab_k
    // so the ftab lookup fits inside it; a deep ftab shrinks the step. When ftab_k >= k
    // there is no room for a look-ahead skip, and step = k - ftab_k - 1 would go negative
    // (making the look-ahead peek past the read end). Clamp to 0 and skip the look-ahead
    // in that case -- an ftab_k-mer already resolves the whole k-mer in one lookup.
    if (k - step < static_cast<int32_t>(ftab_k)) {
        step = k - static_cast<int32_t>(ftab_k) - 1;
    }
    if (step < 0) step = 0;

    while (pos_on_r >= k - 1) {
        if (step > 0 && pos_on_r >= k -1 + step and !look_ahead_backward_search(mq, pos_on_r, step)) {
                        ks.look_ahead_skipped += step + 1;
            pos_on_r = pos_on_r - step - 1;
        } else {
            if (kmer_counts and movi_options->is_kmer_bv()) {
                // Fast count path: presence positive-skip walk + bitvector
                // predecessor/successor to resolve each present k-mer's count.
                uint64_t found = query_kmers_count_bv(mq, pos_on_r);
                                ks.positive_kmers += found;
            } else if (kmer_counts) {
                // Count via the single-kmer presence search (same path as --kmer-bv,
                // minus the MPHF id): the BWT interval size is the k-mer's occurrence
                // count on the doubled (fwd+rc) text = occ(x)+occ(rc(x)) = KMC's
                // canonical count.  The old query_kmers_from_bidirectional path assumed
                // an ftab (it set ftab_right = kmer_left + ftab_k - 1 and asserted the
                // init advanced to kmer_left) and threw with the default ftab_k = 0.
                MoveInterval interval;
                uint64_t found_kmer_count = query_kmers_from(mq, pos_on_r, /*single=*/true, &interval);
                if (found_kmer_count > 0 && !interval.is_empty()) {
                    uint64_t occ_count = interval.count(rlbwt);
                    mq.add_kmer(pos_on_r + 2 - k, found_kmer_count,
                                std::numeric_limits<uint64_t>::max(), occ_count);
                } else {
                    mq.add_kmer(pos_on_r + 2 - k, found_kmer_count);
                }
                                ks.positive_kmers += found_kmer_count;
            } else {
                if (movi_options->is_kmer_bv()) {
                    // Fast keep-going MPHF-id walk: positive-skip presence + per-k-mer
                    // canonical id (lazy-rc inside). 2.2x faster than the per-k-mer
                    // single-search path (1.61 vs 0.73 M/s, ecoli100 k31), with
                    // byte-identical ids. The count field shows presence (1); use
                    // --kmer-count --kmer-bv for occurrence counts.
                    uint64_t found = query_kmers_id_bv(mq, pos_on_r);
                                        ks.positive_kmers += found;
                } else {
                    uint64_t found_kmer_count = query_kmers_from(mq, pos_on_r);
                    // Outputing the kmer matches only works for the non-count mode (for now)
                    mq.add_kmer(pos_on_r + 2 - k, found_kmer_count);
                                        ks.positive_kmers += found_kmer_count;
                }
            }
        }

        while (!check_alphabet(query_seq[pos_on_r])) {
            pos_on_r -= 1; // Find the first position where the character is legal
        }
    }
}
