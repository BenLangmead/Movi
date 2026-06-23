#include <sys/stat.h>
#include <cstdlib>

#include "move_structure.hpp"
#include "utils.hpp"

std::ostringstream dbg;

// IMPORTANT: The index has to be built on a reference with separators to get the correct kmer counts
// IMPORTANT: For the bidirectional search, also the index has to be built with the separators

// pos_on_r points to the furthest (to the right) base of the kmer to be searched
// The search is from the right end of the kmer
// All kmers on the right side of the kmer_middle might be found in this call if they exist
uint64_t MoveStructure::query_kmers_from_bidirectional(MoveQuery& mq, int32_t& pos_on_r) {

    uint64_t kmers_found = 0;
    int32_t last_kmer_looked_at = pos_on_r;

    size_t ftab_k = movi_options->get_ftab_k();
    size_t k = movi_options->get_k();
    auto& query_seq = mq.query();

    // To remember the kmer which initiated the search
    int32_t pos_on_r_saved = pos_on_r;
    // The number of bases matched so far
    uint64_t match_len = 0;
    // The position of the left side of the kmer on the read
    int32_t kmer_left = pos_on_r - k + 1;
    // The position on the right side of the match to be found by ftab
    int32_t ftab_right = kmer_left + ftab_k - 1;

    // At the time, the initialization works if ftab with ftab_k exists only
    // And the multi-ftab strategy must be turned off
    movi_options->set_multi_ftab(false);

    MoveBiInterval bi_interval;
    int32_t ftab_right_initialize = ftab_right;

    bi_interval = initialize_bidirectional_search(mq, ftab_right_initialize, match_len);

    if (match_len == 0 and ftab_k > 1) {
        // If the ftab-k-mer at the end of the read is not found,
        // we can skip all the kmers until ftab_right - 1
        #pragma omp atomic
        kmer_stats.initialize_skipped += (pos_on_r - ftab_right + 1);
        pos_on_r = ftab_right - 1;
        return 0;
    }

    if (ftab_right_initialize != kmer_left) {
        // Sanity check: to make sure ftab and initialization is working
        throw std::runtime_error(ERROR_MSG("[Sequitur - query kmers from bidirectional] " +
                                           "ftab_right_initialize is now updated and should have been equal to kmer_left.\n" +
                                           "ftab_right: " + std::to_string(ftab_right) +
                                           "\tftab_right_initialize:" + std::to_string(ftab_right_initialize) +
                                           "\tkmer_left:" + std::to_string(kmer_left) + "\n"));
    }

    // pos_on_r and pos_on_r_saved point to the beginning of the kmer for which the bidirectional initialization was performed above
    // Extend to right until the kmer is found and save the observed intervals beyond k/2
    // The middle point of the kmer
    int32_t partial_count = k/2;
    int32_t non_extension_count = k - partial_count;
    int32_t kmer_middle = pos_on_r_saved - partial_count;
    std::vector<MoveBiInterval> partial_matches;
    if (partial_count > 0) {
        partial_matches.resize(partial_count - 1); // -1 is because the last kmer is the one found by bidirectional
    }
    int32_t partial_saved_count = 0;

    // kmer_right is the last position matched so far
    uint64_t kmer_right = ftab_right;

    while (kmer_right < pos_on_r_saved) {
        // The next k-mer we are looking at, ends at next_pos
        int next_pos = kmer_right + 1;

        if (movi_options->is_debug() and next_pos >= k)
            dbg << "\nkmer at " << next_pos << ": " << query_seq.substr(next_pos - k + 1, k) << " ";

        bool extend_right_res = extend_right(query_seq[next_pos], bi_interval);

        if (!extend_right_res) {
            // The kmer was not found, we can skip kmers until kmer_right
            #pragma omp atomic
            kmer_stats.right_extension_failed += pos_on_r_saved - kmer_right; // b += pos_on_r_saved - kmer_right;
            pos_on_r = kmer_right;
            last_kmer_looked_at = next_pos;

            // Printing all the kmers being skipped, because the extension to the right was not possible
            if (movi_options->is_debug()) {
                int j = next_pos + 1;
                while (j < pos_on_r_saved) {
                    if (j > k)
                        dbg << "\nkmer at " << j << ": " << query_seq.substr(j - k + 1, k) << "-";
                    j += 1;
                }
            }

            // It's important to break here, to avoid false extensions in the following iterations
            break;
        } else {
            match_len += 1;
            kmer_right = next_pos; // increaments kmer_right one position to the right

            // pos_on_r always points to the right end of the right most kmer
            pos_on_r = kmer_right;

            if (match_len != bi_interval.match_len) {
                // Sanity check: to make sure match_len is computed correctly
                throw std::runtime_error(ERROR_MSG("[Sequitur - query kmers from bidirectional] The interval's match_len is not set correctly.\n" +
                                                   "bi_interval.match_len: " + std::to_string(bi_interval.match_len) +
                                                   " match_len: " + std::to_string(match_len) + "\n"));
            }

            // store the intervals for matches beyond half point of the kmer
            if (kmer_right > kmer_middle and kmer_right != pos_on_r_saved) {

                if (bi_interval.match_len - non_extension_count - 1 >= partial_matches.size()) {
                    throw std::runtime_error(ERROR_MSG("[Sequitur - query kmers from bidirectional] The interval's match_len is not set correctly.\n" +
                                                       "bi_interval.match_len: " + std::to_string(bi_interval.match_len) +
                                                       " non_extension_count: " + std::to_string(non_extension_count) +
                                                       "partial_matches.size(): " + std::to_string(partial_matches.size()) + "\n"));
                }

                if (partial_saved_count != bi_interval.match_len - non_extension_count - 1) {
                    // Sanity check
                    throw std::runtime_error(ERROR_MSG("[Sequitur - query kmers from bidirectional] Partial saved count is not set correctly.\n" +
                                                       "partial_saved_count: " + std::to_string(partial_saved_count) +
                                                       " bi_interval.match_len: " + std::to_string(bi_interval.match_len) +
                                                       " non_extension_count: " + std::to_string(non_extension_count) + "\n"));
                }

                partial_matches[bi_interval.match_len - non_extension_count - 1] = bi_interval;
                partial_saved_count += 1;
            }
        }
    }

    if (kmer_right == pos_on_r_saved) {
        // The kmer at pos_on_r was found by k bidirectional backward search
        kmers_found += 1; // a += 1;
        uint64_t kmers_count = bi_interval.fw_interval.count(rlbwt);
        #pragma omp atomic
        kmer_stats.total_counts += kmers_count;
        if (pos_on_r != pos_on_r_saved) {
            // pos_on_r should be equal to pos_on_r_saved at this point
            throw std::runtime_error(ERROR_MSG("[Sequitur - query kmers from bidirectional] Pos_on_r is not equal to pos_on_r_saved.\n" +
                                               "pos_on_r: " + std::to_string(pos_on_r) +
                                               " pos_on_r_saved: " + std::to_string(pos_on_r_saved) + "\n"));
        }

        // To avoid searching this kmer again in the next step
        kmer_right -= 1;
        pos_on_r = kmer_right;

        if (movi_options->is_debug()) {
            if (pos_on_r_saved > k)
                dbg << "\nkmer at " << pos_on_r_saved << ": " << query_seq.substr(pos_on_r_saved - k + 1, k) << "\n";
            dbg << "1";
        }
    } else {
        if (movi_options->is_debug()) {
            if (pos_on_r_saved > k)
                dbg << "\nkmer at " << pos_on_r_saved << ": " << query_seq.substr(pos_on_r_saved - k + 1, k) << "-\n";
            dbg << "0";
        }
    }


    // Use partial matches for finding other overlapping kmers (until half point)
    if (kmer_right > kmer_middle) {
        if (movi_options->is_debug()) {
            for (int i = kmer_right + 1; i < pos_on_r_saved; i++) {
                dbg << "0";
            }
        }

        int skip_kmers = 0;
        int partial_found = 0;

        for (int i = 0; i < partial_saved_count; i += (skip_kmers + 1) ) {
            // Set kmer_left_ext to be the last match position on the left end
            int32_t kmer_left_ext = kmer_left;
            auto& partial_match_interval = partial_matches[i];

            kmer_right = kmer_middle + 1 + i;
            last_kmer_looked_at = kmer_right; // not used

            while (partial_match_interval.match_len < k and kmer_left_ext > 0) {
                // We don't need to do bidirectional left extension here, simple backward search is enough
                bool res = backward_search_step(query_seq[kmer_left_ext - 1], partial_match_interval.fw_interval);
                // bool res = extend_left(query_seq[kmer_left_ext - 1], partial_match_interval);
                if (!res) {
                    // The current kmer is not present, move to the next partial match by breaking from the inner loop
                    break;
                } else {
                    kmer_left_ext -= 1;
                    partial_match_interval.match_len += 1;
                }
            }

            // There is a possibility to a number of next kmers because of their overlap
            skip_kmers = (partial_match_interval.match_len < k - 1) ? k - partial_match_interval.match_len - 1 : 0;

            if (partial_match_interval.match_len == k) {
                // The kmer was found by extending the partial match to left
                kmers_found += 1;
                partial_found += 1;

                uint64_t kmers_count = partial_match_interval.fw_interval.count(rlbwt);
                #pragma omp atomic
                kmer_stats.total_counts += kmers_count;
                #pragma omp atomic
                kmer_stats.positive_skipped += 1; // c += 1;

                if (movi_options->is_debug()) {
                    DEBUG_MSG("kmer at " + std::to_string(kmer_left_ext + k - 1) + " was found.");
                    dbg << "1";
                }
            } else {
                skip_kmers = 0;
                #pragma omp atomic
                kmer_stats.backward_search_failed += skip_kmers + 1; // d += skip_kmers + 1;
                if (movi_options->is_debug()) {
                    dbg << "0";
                }
            }
        }

        // To check the correct count of the total kmers
        // if (a+b+c+d != 15) std::cerr << a << "\t" << b << "\t" << c << "\t" << d << "\t|" << partial_saved_count << "|\n";

        // At this point all the kmers to the right of kmer_middle are checked using the partial matches
        pos_on_r = kmer_middle;

    } else {
        // If we got here, we had to start the first while by breaking because of an unsuccessfull attempt to extend to the right
        // kmers_found = 0;
        if (kmer_right != pos_on_r)
            throw std::runtime_error(ERROR_MSG("[Sequitur - query kmers from bidirectional] This should not happen: "
                                               + std::to_string(kmer_right) + "\t" + std::to_string(pos_on_r)));
        // pos_on_r should have already been assigned to be the last kmer_right
        // So we should never get here to do the assignment in practice
        pos_on_r = kmer_right;

        if (movi_options->is_debug()) {
            for (int i = kmer_middle + 1; i <= pos_on_r_saved - 1; i++) {
                dbg << "0";
            }
        }

    }
    return kmers_found;
}

uint64_t MoveStructure::query_kmers_from(MoveQuery& mq, int32_t& pos_on_r, bool single,
                                          MoveInterval* interval_out) {
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
            #pragma omp atomic
            kmer_stats.initialize_skipped += 1;
            pos_on_r -= 1;
            pos_on_r_saved = pos_on_r;
        }
    } while (match_len == 0 and pos_on_r >= k - 1 and ftab_k > 1);

    // I want to check how much slower it gets if turn off the positive skip:
    auto backward_search_result = backward_search(query_seq, pos_on_r, initial_interval, single ? k - match_len - 2 : std::numeric_limits<int32_t>::max());

    if (backward_search_result.is_empty()) {
        // We get here when there is an illegal character at pos_on_r, just skip the current position
        pos_on_r = pos_on_r_saved - 1;
        #pragma omp atomic
        kmer_stats.backward_search_empty += 1;
        return 0;
    } else {
        if (pos_on_r_saved - pos_on_r >= k - 1) {
            // At leat one kmer was found, update the postion and return the count
            uint64_t kmers_found = pos_on_r_saved - pos_on_r - k + 2;
            #pragma omp atomic
            kmer_stats.positive_skipped += kmers_found - 1;

            if (interval_out) *interval_out = backward_search_result;

            if (movi_options->is_debug()) {
                int32_t pos_on_r_ = pos_on_r_saved;
                auto backward_search_result_one_extra_base = backward_search(query_seq, pos_on_r_, initial_interval, k - match_len - 1);
                dbg << backward_search_result.count(rlbwt) << "----" <<  backward_search_result_one_extra_base.count(rlbwt) << "\n";
                dbg << backward_search_result << "\n";
                if (backward_search_result.count(rlbwt) == backward_search_result_one_extra_base.count(rlbwt)) {
                    // TODO: Use a counter to count the number of such incidents
                } else {
                }
            }

            pos_on_r = pos_on_r + k - 2;
	        return kmers_found;
        } else {
            // No kmer was found, update the postion
            #pragma omp atomic
            kmer_stats.backward_search_failed += 1;
            pos_on_r = pos_on_r_saved - 1;
            return 0;
        }
    }
}

// Count via the presence positive-skip walk + bitvector predecessor/successor.
// This is the fast counting path: it does EXACTLY what the presence query does
// (one continuous unidirectional backward search per run of present k-mers,
// ftab-assisted, no bidirectional/rc bookkeeping), and resolves a count for
// every present k-mer along the way using the count bitvector.
//
// As the backward search extends left, the current interval covers
// read[pos_on_r .. anchor] (anchor = the right end). Whenever the covered
// length reaches k, the LEFTMOST k-mer of the match, read[pos_on_r..pos_on_r+k-1],
// is fully matched and the current (possibly longer-match, hence shrinking)
// interval [lb,rb] is a SUBSET of that k-mer's interval. The k-mer's interval is
// a single group in the count bitvector, so predecessor(lb)=its lb and
// successor(rb)=its rb+1, giving count = succ(rb)-pred(lb) in O(1). When the
// range empties we return; the outer loop re-anchors at the next k-mer.
uint64_t MoveStructure::query_kmers_count_bv(MoveQuery& mq, int32_t& pos_on_r) {
    size_t ftab_k = movi_options->get_ftab_k();
    int32_t k = static_cast<int32_t>(movi_options->get_k());
    auto& query_seq = mq.query();
    int32_t pos_on_r_saved = pos_on_r;

    uint64_t match_len = 0;
    MoveInterval interval;
    do {
        interval = initialize_backward_search(mq, pos_on_r, match_len);
        if (match_len == 0 and ftab_k > 1) {
            pos_on_r -= 1;
            pos_on_r_saved = pos_on_r;
        }
    } while (match_len == 0 and pos_on_r >= k - 1 and ftab_k > 1);

    if (interval.is_empty()) {
        pos_on_r = pos_on_r_saved - 1;
        return 0;
    }

    // anchor = right end of the match; covered length = anchor - pos_on_r + 1.
    const int32_t anchor = pos_on_r_saved;
    uint64_t kmers_found = 0;
    while (true) {
        if (anchor - pos_on_r + 1 >= k) {
            // Resolve the count of the leftmost k-mer (starts at pos_on_r).
            uint64_t cnt = kmer_count_from_bv(interval);
            mq.add_kmer(pos_on_r, /*present=*/1, std::numeric_limits<uint64_t>::max(), cnt);
            kmers_found += 1;
        }
        if (pos_on_r <= 0) break;
        MoveInterval prev = interval;
        if (!backward_search_step(query_seq[pos_on_r - 1], interval)) {
            interval = prev;  // step emptied the interval; stop this run
            break;
        }
        pos_on_r -= 1;
    }

    if (kmers_found > 0) {
        #pragma omp atomic
        kmer_stats.positive_skipped += kmers_found - 1;
        // pos_on_r is the leftmost resolved k-mer's start; the next k-mer to
        // process ends at pos_on_r + k - 2 (mirrors query_kmers_from).
        pos_on_r = pos_on_r + k - 2;
    } else {
        pos_on_r = pos_on_r_saved - 1;
    }
    return kmers_found;
}

// MPHF id via the fast presence positive-skip walk (same keep-going structure as
// query_kmers_count_bv), emitting each present k-mer's canonical id. Needs only
// B_k + all_p -- NO count structure: for a present k-mer the walk's (subset)
// interval [lb_row..] sits inside the k-mer's group, and the k-mer's own B_k mark
// is the only 1 in that group, so id = rank(lb_row + 1) - 1 regardless of how much
// the match overshot k. Canonical mode: a canonical k-mer uses lb_row directly; a
// non-canonical k-mer's id comes from an ftab-assisted rc search (lazy-rc). This
// keeps the positive-skip heuristic that makes presence fast, unlike the per-k-mer
// single search the default --kmer-bv path uses.
uint64_t MoveStructure::query_kmers_id_bv(MoveQuery& mq, int32_t& pos_on_r) {
    size_t ftab_k = movi_options->get_ftab_k();
    int32_t k = static_cast<int32_t>(movi_options->get_k());
    auto& query_seq = mq.query();
    int32_t pos_on_r_saved = pos_on_r;

    uint64_t match_len = 0;
    MoveInterval interval;
    do {
        interval = initialize_backward_search(mq, pos_on_r, match_len);
        if (match_len == 0 and ftab_k > 1) {
            pos_on_r -= 1;
            pos_on_r_saved = pos_on_r;
        }
    } while (match_len == 0 and pos_on_r >= k - 1 and ftab_k > 1);

    if (interval.is_empty()) {
        pos_on_r = pos_on_r_saved - 1;
        return 0;
    }

    const int32_t anchor = pos_on_r_saved;
    uint64_t kmers_found = 0;
    while (true) {
        if (anchor - pos_on_r + 1 >= k) {
            uint64_t lb_for_rank = kmerbv_all_p[interval.run_start] + interval.offset_start;
            int orientation = -1;
            if (kmerbv_is_canonical) {
                std::string xstr = query_seq.substr(pos_on_r, k);
                std::string rcstr = reverse_complement(xstr);
                if (xstr <= rcstr) {
                    orientation = 0;  // canonical: the subset interval's lb_row works
                } else {
                    MoveQuery rcq(rcstr);
                    int32_t rp = k - 1; uint64_t rml = 0;
                    MoveInterval rc_iv = initialize_backward_search(rcq, rp, rml);
                    rc_iv = backward_search(rcq.query(), rp, rc_iv,
                                            std::numeric_limits<int32_t>::max());
                    if (!rc_iv.is_empty()) {
                        lb_for_rank = kmerbv_all_p[rc_iv.run_start] + rc_iv.offset_start;
                        orientation = 1;
                    } else {
                        orientation = 0;
                    }
                }
            }
            uint64_t kmer_id = kmerbv_rank1(lb_for_rank + 1) - 1;
            mq.add_kmer(pos_on_r, /*present=*/1, kmer_id,
                        std::numeric_limits<uint64_t>::max(), orientation);
            kmers_found += 1;
        }
        if (pos_on_r <= 0) break;
        MoveInterval prev = interval;
        if (!backward_search_step(query_seq[pos_on_r - 1], interval)) {
            interval = prev;
            break;
        }
        pos_on_r -= 1;
    }

    if (kmers_found > 0) {
        #pragma omp atomic
        kmer_stats.positive_skipped += kmers_found - 1;
        pos_on_r = pos_on_r + k - 2;
    } else {
        pos_on_r = pos_on_r_saved - 1;
    }
    return kmers_found;
}

void MoveStructure::query_all_kmers(MoveQuery& mq, bool kmer_counts) {
    size_t ftab_k = movi_options->get_ftab_k();
    size_t k = movi_options->get_k();
    auto& query_seq = mq.query();
    int32_t pos_on_r = query_seq.length() - 1;

    // To handle a special case for k equal to 1
    if (k == 1) {
        uint64_t kmers_found = 0;
        while (pos_on_r >= 0) {
            #pragma omp atomic
            kmers_found += check_alphabet(query_seq[pos_on_r]) ? 1 : 0;
            pos_on_r -= 1;
        }
        #pragma omp atomic
        kmer_stats.positive_kmers += kmers_found;
        return;
    }

    while (!check_alphabet(query_seq[pos_on_r])) {
        pos_on_r -= 1; // Find the first position where the character is legal
    }


    int32_t step = k/3;
    // k - step has to be always greater than ftab-k
    if (k - step < ftab_k) {
        step = k - ftab_k - 1;
    }

    while (pos_on_r >= k - 1) {
        if (pos_on_r >= k -1 + step and !look_ahead_backward_search(mq, pos_on_r, step)) {
            #pragma omp atomic
            kmer_stats.look_ahead_skipped += step + 1;
            pos_on_r = pos_on_r - step - 1;
        } else {
            if (kmer_counts and movi_options->is_kmer_bv()) {
                // Fast count path: presence positive-skip walk + bitvector
                // predecessor/successor to resolve each present k-mer's count.
                uint64_t found = query_kmers_count_bv(mq, pos_on_r);
                #pragma omp atomic
                kmer_stats.positive_kmers += found;
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
                #pragma omp atomic
                kmer_stats.positive_kmers += found_kmer_count;
            } else {
                if (movi_options->is_kmer_bv()) {
                    // Fast keep-going MPHF-id walk: positive-skip presence + per-k-mer
                    // canonical id (lazy-rc inside). 2.2x faster than the per-k-mer
                    // single-search path (1.61 vs 0.73 M/s, ecoli100 k31), with
                    // byte-identical ids. The count field shows presence (1); use
                    // --kmer-count --kmer-bv for occurrence counts.
                    uint64_t found = query_kmers_id_bv(mq, pos_on_r);
                    #pragma omp atomic
                    kmer_stats.positive_kmers += found;
                } else {
                    uint64_t found_kmer_count = query_kmers_from(mq, pos_on_r);
                    // Outputing the kmer matches only works for the non-count mode (for now)
                    mq.add_kmer(pos_on_r + 2 - k, found_kmer_count);
                    #pragma omp atomic
                    kmer_stats.positive_kmers += found_kmer_count;
                }
            }
        }

        while (!check_alphabet(query_seq[pos_on_r])) {
            pos_on_r -= 1; // Find the first position where the character is legal
        }
    }
}
