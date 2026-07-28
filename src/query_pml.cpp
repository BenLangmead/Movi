#include "move_structure.hpp"

uint64_t MoveStructure::query_pml(MoveQuery& mq) {
    auto& R = mq.query();
    int32_t pos_on_r = R.length() - 1;
    uint64_t idx = r - 1; // or we can start from a random position in the rlbwt std::rand() % r
    uint64_t offset = get_n(idx) - 1;

    uint64_t match_len = 0;
    uint16_t ff_count = 0;
    uint64_t ff_count_tot = 0;
    uint64_t scan_count = 0;
    auto t1 = std::chrono::high_resolution_clock::now();

    if (movi_options->is_debug()) {
        DEBUG_MSG("beginning of the search \ton query: " + mq.query() + "\tand on BWT, idx(r-1): " + std::to_string(idx)
                  + " offset: " + std::to_string(offset));
    }

    // Multi-class classification
    if (movi_options->is_multi_classify()) {
        for (uint16_t i = 0; i < num_species; i++) {
            if (!movi_options->is_pvalue_scoring()) {
                classify_cnts[i] = 0;
            } else {
                doc_scores[i] = 0;
            }
        }
    }

    uint16_t best_doc = std::numeric_limits<uint16_t>::max(); // for multi-class classification
    uint16_t second_best_doc = std::numeric_limits<uint16_t>::max();
    uint64_t iteration_count = 0;
    uint32_t sum_matching_lengths = 0;
    while (pos_on_r > -1) {
        iteration_count += 1;
        if (movi_options->is_logs() and (iteration_count-1)%200 == 0) {
            t1 = std::chrono::high_resolution_clock::now();
        }

        if (movi_options->is_debug())
            DEBUG_MSG("Searching position " + std::to_string(pos_on_r) + " of the read:");

        auto& row = rlbwt[idx];
        uint64_t row_idx = idx;
        char row_c = alphabet[row.get_c()];
        if (!check_alphabet(R[pos_on_r])) {
            // The character from the read does not exist in the reference
            match_len = 0;
            scan_count = 0;

            if (movi_options->is_debug())
                DEBUG_MSG("\t The character " + std::string(1, R[pos_on_r]) + " does not exist.");
        } else if (row_c == R[pos_on_r]) {
            // Case 1
            match_len += 1;
            scan_count = 0;

            if (movi_options->is_debug()) {
                DEBUG_MSG("\tCase 1: It was a match.\n\t Continue the search...");
                DEBUG_MSG("\tmatch_len: " + std::to_string(match_len));
                DEBUG_MSG("\tcurrent_id: " + std::to_string(idx) + "\trow.id: " + std::to_string(get_id(row_idx)));
                DEBUG_MSG("\trow.get_n: " + std::to_string(get_n(row_idx)) + "\trlbwt[idx].get_n: " + std::to_string(get_n(get_id(row_idx))));
                DEBUG_MSG("\toffset: " + std::to_string(offset) + "\trow.get_offset(): " + std::to_string(get_offset(row_idx)));
            }
        } else {
            // Case 2
            // Repositioning up or down (randomly or with thresholds)
            if (movi_options->is_debug())
                DEBUG_MSG("\t Case 2: Not a match, looking for a match either up or down...");

            uint64_t idx_before_reposition = idx;
#if USE_THRESHOLDS
            bool up = movi_options->is_random_repositioning() ?
                               reposition_randomly(idx, offset, R[pos_on_r], scan_count) :
                               reposition_thresholds(idx, offset, R[pos_on_r], scan_count);
#else
            // When there is no threshold, reposition randomly
            bool up = reposition_randomly(idx, offset, R[pos_on_r], scan_count);
#endif

            match_len = 0;
            // scan_count = (!constant) ? std::abs((int)idx - (int)idx_before_reposition) : 0;

            char c = alphabet[rlbwt[idx].get_c()];

            if (movi_options->is_debug())
                DEBUG_MSG("\tup: " + std::to_string(up) + " idx: " + std::to_string(idx) + " c:" + c);

            // sanity check
            if (c == R[pos_on_r]) {
                // Observing a match after the reposition
                // The right match_len should be:
                // min(new_lcp, match_len + 1)
                // But we cannot compute lcp here
                offset = up ? get_n(idx) - 1 : 0;

                if (movi_options->is_debug())
                    DEBUG_MSG("\tidx: " + std::to_string(idx) + " offset: " + std::to_string(offset));

            } else {
                DEBUG_MSG("\t\tpos: " + std::to_string(pos_on_r) + " r[pos]:" +  R[pos_on_r] + " t[pointer]:" + c);
                DEBUG_MSG("\t\t" + std::to_string(up) + ", " + R[pos_on_r] + ", " + std::to_string(pos_on_r));
                DEBUG_MSG("\t\t");
                for (int k = 10; k > 0; --k)
                    DEBUG_MSG(alphabet[rlbwt[idx - k].get_c()] + "-");
                for (int k = 0; k < 10; k++)
                    DEBUG_MSG(alphabet[rlbwt[idx + k].get_c()] + "-");
                DEBUG_MSG("\n");

#if USE_THRESHOLDS
                movi_options->set_verbose(true);
                movi_options->set_debug(true);
                auto saved_idx = idx;
                reposition_thresholds(saved_idx, offset, R[pos_on_r], scan_count);
#endif
                throw std::runtime_error(ERROR_MSG("[query pml] This should not happen!"));
            }
        }
    
        sum_matching_lengths += match_len;
        mq.add_ml(match_len, movi_options->is_stdout());
        if (movi_options->is_get_sa_entries()) {
            uint64_t sa_entry = get_SA_entries(idx, offset);
            mq.add_sa_entries(sa_entry);
        }
        pos_on_r -= 1;

        // LF step
        ff_count = LF_move(offset, idx);
        ff_count_tot += ff_count;
        if (movi_options->is_logs()) {
            if (iteration_count % 200 == 0) {
                auto t2 = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
                mq.add_cost(elapsed);
            }
            mq.add_fastforward(ff_count);
            mq.add_scan(scan_count);
        }

        if (movi_options->is_multi_classify()) {
            if (match_len >= movi_options->get_min_match_len()) {
                /*uint64_t full_ind = run_offsets[idx] + offset;
                uint16_t cur_doc = doc_pats[full_ind];
                classify_cnts[cur_doc]++;
                if (classify_cnts[cur_doc] >= classify_cnts[best_doc]) {
                    best_doc = cur_doc;
                }*/ 

                uint64_t color_id;
#if COLOR_MODE == 1
                color_id = static_cast<uint64_t>(rlbwt[idx].color_id);
                // Skip doc sets that weren't saved (thrown away by compression).
                if (color_id >= unique_doc_sets.size()) continue;
#else
                if (movi_options->is_doc_sets_vector_of_vectors()) {
                    color_id = static_cast<uint64_t>(doc_set_inds[idx]);
                    // Skip doc sets that weren't saved (thrown away by compression).
                    if (color_id >= unique_doc_sets.size()) continue;
                } else {
                    color_id = doc_set_flat_inds[idx].get();
                    // Skip doc sets that weren't saved (thrown away by compression).
                    if (color_id >= flat_colors.size()) continue;
                }
#endif
                std::span<uint16_t> cur_set;
                if (movi_options->is_doc_sets_vector_of_vectors()) {
                    std::vector<uint16_t> &cur_set_vec = unique_doc_sets[color_id];
                    cur_set = std::span<uint16_t>(cur_set_vec.data(), cur_set_vec.size());
                } else {
                    uint32_t cur_set_size = flat_colors[color_id];
                    cur_set = std::span<uint16_t>(flat_colors.data() + color_id + 1,
                                                  flat_colors.data() + color_id + 1 + cur_set_size);
                }

                for (int doc : cur_set) {
                    if (!movi_options->is_pvalue_scoring()) {
                        classify_cnts[doc]++;
                        if (doc != best_doc) {
                            if (best_doc == std::numeric_limits<uint16_t>::max() || classify_cnts[doc] > classify_cnts[best_doc]) {
                                second_best_doc = best_doc;
                                best_doc = doc;
                            } else if (second_best_doc == std::numeric_limits<uint16_t>::max() || classify_cnts[doc] > classify_cnts[second_best_doc]) {
                                second_best_doc = doc;
                            }
                        }
                    } else {
                        // p value strategy
                        double val = match_len - (log_lens[doc] / log4);
                        if (val >= 0) {
                            doc_scores[doc] += std::min(val, 1.);
                            if (doc != best_doc) {
                                if (best_doc == std::numeric_limits<uint16_t>::max() || doc_scores[doc] > doc_scores[best_doc]) {
                                    second_best_doc = best_doc;
                                    best_doc = doc;
                                } else if (second_best_doc == std::numeric_limits<uint16_t>::max() || doc_scores[doc] > doc_scores[second_best_doc]) {
                                    second_best_doc = doc;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (movi_options->is_multi_classify()) {

        output_files->out_file << mq.get_query_id() << ",";

        float PML_mean = static_cast<float>(sum_matching_lengths) / mq.query().length();
        if (PML_mean < UNCLASSIFIED_THRESHOLD || best_doc == std::numeric_limits<uint16_t>::max()) {
            // Not present
            output_files->out_file << "0,0\n";
        } else {
            if (second_best_doc == std::numeric_limits<uint16_t>::max()) {
                output_files->out_file << to_taxon_id[best_doc] << ",0";
            } else {
                float best_doc_cnt, second_best_doc_cnt, second_best_diff;
                if (movi_options->get_min_match_len() > 0) {
                    best_doc_cnt = classify_cnts[best_doc];
                    second_best_doc_cnt = classify_cnts[second_best_doc];
                    second_best_diff = (best_doc_cnt - second_best_doc_cnt);
                } else {
                    // p-value strategy
                    best_doc_cnt = doc_scores[best_doc];
                    second_best_doc_cnt = doc_scores[second_best_doc];
                    second_best_diff = (best_doc_cnt - second_best_doc_cnt);
                }

                if (second_best_diff < 0.05 * best_doc_cnt) {
                    output_files->out_file << to_taxon_id[best_doc] << "," << to_taxon_id[second_best_doc];
                } else {
                    output_files->out_file << to_taxon_id[best_doc] << ",0";
                }
            }
            output_files->out_file << "\n";
        }
    }

    return ff_count_tot;
}
