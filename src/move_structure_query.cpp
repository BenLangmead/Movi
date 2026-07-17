#include "move_structure.hpp"

void MoveStructure::reconstruct_lf() {
    uint64_t offset = 0;
    uint64_t run_index = 0;
    uint64_t i = 0;
    uint64_t ff_count_tot = 0;

    uint64_t total_elapsed = 0;
    // auto begin = std::chrono::system_clock::now();
    for (; run_index != end_bwt_idx; ) {
        if (i % 1000000 == 0) {
            print_progress_bar(i, length - 1, "Reconstructing the original text");
        }

        ff_count_tot += LF_move(offset, run_index);

        i += 1;
        // orig_string = rlbwt[run_index].get_c() + orig_string;
    }
    // auto end = std::chrono::system_clock::now();
    // auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    // total_elapsed = elapsed.count();
    PROGRESS_MSG("Finished reconstructing the original string.");
    INFO_MSG("Total fast forward: " + std::to_string(ff_count_tot));
    TIMING_MSG(total_elapsed, "reconstructing the original text");
}

void MoveStructure::sequential_lf() {
    uint64_t line_index = 0;
    uint64_t row_index = 0;
    uint64_t ff_count_tot = 0;

    uint64_t total_elapsed = 0;
    for (uint64_t row_index = 0; row_index < r; row_index++) {
        auto& current = rlbwt[row_index];
        for (uint64_t j = 0; j < current.get_n(); j ++) {
            if (line_index % 1000000 == 0) {
                print_progress_bar(line_index, length - 1, "LF-mapping for all the BWT characters");
            }

            uint64_t offset = j;
            uint64_t i = row_index;
            auto begin = std::chrono::system_clock::now();
            ff_count_tot += LF_move(offset, i);
            auto end = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            total_elapsed += elapsed.count();
            line_index += 1;
        }

    }
    PROGRESS_MSG("Finished performing LF-mapping for all the BWT characters.");
    INFO_MSG("Total fast forward: " + std::to_string(ff_count_tot));
    TIMING_MSG(total_elapsed, "LF-mapping for all the BWT characters");
}

void MoveStructure::random_lf() {
    uint64_t ff_count_tot = 0;

    // generate the random order from 1 to length
    std::vector<uint64_t> random_order(length);
    for (uint64_t i = 0; i < length; i++) {
        random_order[i] = i;
    }
    std::srand(time(0));
    std::random_shuffle(random_order.begin(), random_order.end());

    // find the n and id for each random BWT row
    std::vector<uint64_t> n_to_id(length);
    std::vector<uint64_t> id_to_p(r);
    uint64_t current_n = 0;
    for (uint64_t i = 0; i < r; i++) {
        id_to_p[i] = current_n;
        for (uint64_t j = 0; j < get_n(i); j++) {
            n_to_id[current_n] = i;
            current_n += 1;
        }
    }


    uint64_t total_elapsed = 0;
    for (uint64_t i = 0; i < length; i++) {
        if (i % 1000000 == 0) {
            print_progress_bar(i, length - 1, "LF-mapping for all the BWT characters in the random order");
        }

        // uint64_t n = std::rand() % r;
        uint64_t n = random_order[i];
        uint64_t id = n_to_id[n];
        uint64_t offset = n - id_to_p[id];
        auto begin = std::chrono::system_clock::now();
        ff_count_tot += LF_move(offset, id);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        total_elapsed += elapsed.count();
    }
    PROGRESS_MSG("Finished performing LF-mapping for all the BWT characters in the random order.");
    INFO_MSG("Total fast forward: " + std::to_string(ff_count_tot));
    TIMING_MSG(total_elapsed, "LF-mapping for all the BWT characters in the random order");
}

void MoveStructure::verify_lfs() {
    // This function only works if the LF function works correctly and
    // that depends on populating the occs bitvector correctly.
    // This function is only used for debugging purposes.
    uint64_t not_matched = 0;
    for (uint64_t i = 0; i < all_p.size(); i++) {
        std::uint64_t end_ = (i < all_p.size() - 1) ? all_p[i + 1] : length;
        for (uint64_t j = all_p[i]; j < end_; j++) {
            uint64_t offset_ = j - all_p[i];
            uint64_t idx_ = i;
            uint64_t lf = 0;
            if (i != end_bwt_idx) {
                uint64_t alphabet_index = rlbwt[idx_].get_c();
                lf = LF(j, alphabet_index);
            } else {
                WARNING_MSG("end_run = " + std::to_string(i) + " len: " + std::to_string(rlbwt[i].get_n ()));
            }
            LF_move(offset_, idx_);
            uint64_t lf_move = all_p[idx_] + offset_;
            if (lf != lf_move) {
                not_matched += 1;
                DEBUG_MSG("j\t" + std::to_string(j) + "\n");
                DEBUG_MSG("idx\t" + std::to_string(i) + "\n");
                DEBUG_MSG("offset\t" + std::to_string(j - all_p[i]) + "\n");
                DEBUG_MSG("rlbwt[idx].get_id\t" + std::to_string(get_id(i)) + "\n");
                DEBUG_MSG("get_offset(i)\t" + std::to_string(get_offset(i)) + "\n");
                for (uint64_t k = 0; k <= i; k++) {
                    DEBUG_MSG(rlbwt[k].get_n() + " ");
                }
                DEBUG_MSG("\n\n");

                DEBUG_MSG("lf\t" + std::to_string(lf) + "\n");
                DEBUG_MSG("lf_move\t" + std::to_string(lf_move) + "\n");
                DEBUG_MSG("idx_\t" + std::to_string(idx_) + "\n");
                DEBUG_MSG("offset_\t" + std::to_string(offset_) + "\n");
                DEBUG_MSG("all_p[idx_]\t" + std::to_string(all_p[idx_]) + "\n");
                DEBUG_MSG("\n\n\n");
            }
        }
    }
    if (not_matched == 0) {
        SUCCESS_MSG("All the LF_move operations are correct.");
    } else {
        original_r = 0;
        WARNING_MSG("There are " + std::to_string(not_matched) + " LF_move operations that failed to match the true lf results.");
    }
}

void MoveStructure::verify_lf_loop() {
    // This function only works if the LF function works correctly and
    // It verifies that n LF_move operations loops back to the same BWT offeset.
    uint64_t idx = end_bwt_idx;
    uint64_t offset = 0;

    std::vector<std::vector<uint64_t>> bwt_offsets(r);
    for (uint64_t i = 0; i < r; i++) {
        bwt_offsets[i].resize(rlbwt[i].get_n(), 0);
    }

    INFO_MSG("Verifying that all the LF_move operations loop back to the same BWT offset.");

    for (uint64_t i = 0; i < length; i++) {
        LF_move(offset, idx);
        bwt_offsets[idx][offset] = 1;
    }

    uint64_t visited_offsets = 0;
    for (uint64_t i = 0; i < r; i++) {
        for (uint64_t j = 0; j < rlbwt[i].get_n(); j++) {
            if (bwt_offsets[i][j] == 1) {
                visited_offsets += 1;
            }
        }
    }

    if (idx == end_bwt_idx and offset == 0 and visited_offsets == length) {
        SUCCESS_MSG("All the LF_move operations are correct.");
    } else {
        INFO_MSG("\tlast idx: " + std::to_string(idx) + " last offset: " + std::to_string(offset));
        INFO_MSG("\tNumber of visited offsets: " + std::to_string(visited_offsets));
        INFO_MSG("\tlength of the BWT: " + std::to_string(length));
        throw std::runtime_error(ERROR_MSG("[verify lf_loop] LF_move operations failed to loop back to the same BWT offset."));
    }
}

uint64_t MoveStructure::reposition_up(uint64_t idx, char c, uint64_t& scan_count) {
    if (idx == 0)
        return r;
    char row_c = alphabet[rlbwt[idx].get_c()];

    while (idx > 0 and row_c != c) {
        scan_count += 1;
        idx -= 1;
        row_c = alphabet[rlbwt[idx].get_c()];
    }

    /* if (logs) {
        if (repositions.find(scan_count) != repositions.end())
            repositions[scan_count] += 1;
        else
            repositions[scan_count] = 1;
    } */

    if (movi_options->is_debug())
        DEBUG_MSG("\tidx after the while in the reposition up: " + std::to_string(idx));
    return (row_c == c) ? idx : r;
}

uint64_t MoveStructure::reposition_down(uint64_t idx, char c, uint64_t& scan_count) {
    if (idx == r - 1)
        return r;
    char row_c = alphabet[rlbwt[idx].get_c()];

    while (idx < r - 1 && row_c != c) {
        scan_count += 1;
        idx += 1;
        row_c = alphabet[rlbwt[idx].get_c()];
    }

    /* if (logs) {
        if (repositions.find(scan_count) != repositions.end())
            repositions[scan_count] += 1;
        else
            repositions[scan_count] = 1;
    } */

    if (movi_options->is_debug())
        DEBUG_MSG("\tidx after the while in the reposition down: " + std::to_string(idx) + " " + c + " " + row_c);
    return (row_c == c) ? idx : r;
}


#if USE_THRESHOLDS
void MoveStructure::handle_reposition_up(uint64_t& idx, uint64_t saved_idx, char r_char, uint16_t next_up, uint64_t& scan_count) {
#if USE_NEXT_POINTERS
    if (constant) {
        scan_count += 1;
        if (next_up == std::numeric_limits<uint16_t>::max())
            idx = r;
        else
            idx = saved_idx - next_up;
    } else {
        throw std::runtime_error(ERROR_MSG("[reposition thresholds] MODE is set to " + std::to_string(MODE) +
                                           ", but the constant variable is false."));
   }
#endif
#if THRESHOLDS_WITHOUT_NEXTS
    idx = reposition_up(saved_idx, r_char, scan_count);
#endif
}

void MoveStructure::handle_reposition_down(uint64_t& idx, uint64_t saved_idx, char r_char, uint16_t next_down, uint64_t& scan_count) {
#if USE_NEXT_POINTERS
    if (constant) {
        scan_count += 1;
        if (next_down == std::numeric_limits<uint16_t>::max())
            idx = r;
        else
            idx = saved_idx + next_down;
    } else {
         throw std::runtime_error(ERROR_MSG("[reposition thresholds] MODE is set to " + std::to_string(MODE) +
                                            ", but the constant variable is false."));
    }
#endif
#if THRESHOLDS_WITHOUT_NEXTS
    idx = reposition_down(saved_idx, r_char, scan_count);
#endif
}

bool MoveStructure::reposition_thresholds(uint64_t& idx, uint64_t offset, char r_char, uint64_t& scan_count) {
    // If offset is greather than or equal to the threshold, reposition down
    // otherwise, reposition up
    uint64_t saved_idx = idx;
    uint64_t alphabet_index = alphamap[static_cast<uint64_t>(r_char)];
    if (use_separator()) {
        if (alphabet_index == 0) {
            throw std::runtime_error(ERROR_MSG("[reposition thresholds] the alphabet index equal to 0 should not happen with separators."));
        }
        alphabet_index -= 1;
    }
    scan_count = 0;

    char rlbwt_char = alphabet[rlbwt[idx].get_c()];

    uint64_t threshold_value = 0;

    // Used in constant mode where we store explicit pointers for repositioning
    uint16_t next_down = 0;
    uint16_t next_up = 0;

    if (idx == end_bwt_idx) {
        threshold_value = end_bwt_idx_thresholds[alphabet_index];
#if USE_NEXT_POINTERS
        next_down = end_bwt_idx_next_down[alphabet_index];
        next_up = end_bwt_idx_next_up[alphabet_index];
#endif
    } else if (use_separator() and rlbwt_char == SEPARATOR) {
        threshold_value = separators_thresholds[separators_thresholds_map[idx]].values[alphabet_index];
        // TODO: Next pointers with separators are not implemented yet
        // Constant mode doesn't store next pointers for separators
    } else {
        if (use_separator()) {
            if (movi_options->is_debug()) {
                DEBUG_MSG("Use separators: rlbwt_char = " + std::string(1, rlbwt_char) +
                          " alphabet_index before applying alphamap_3: " + std::to_string(alphabet_index));
            }
            alphabet_index = alphamap_3[alphamap[rlbwt_char] - 1][alphabet_index];
        } else {
            if (movi_options->is_debug()) {
                DEBUG_MSG("No separators: rlbwt_char = " + std::string(1, rlbwt_char) + " alphabet_index: " + std::to_string(alphabet_index));
            }
            alphabet_index = alphamap_3[alphamap[rlbwt_char]][alphabet_index];
        }

        if (alphabet_index == 3) {
            throw std::runtime_error(ERROR_MSG("[reposition thresholds] alphamap_3 is incorrect, alphabet_index = " + std::to_string(alphabet_index)));
        }
        threshold_value = get_thresholds(idx, alphabet_index);
#if USE_NEXT_POINTERS
        next_down = rlbwt[idx].get_next_down(alphabet_index);
        next_up = rlbwt[idx].get_next_up(alphabet_index);
#endif
    }

    if (movi_options->is_debug()) {
        DEBUG_MSG("[reposition_thresholds] alphabet_index: " + std::to_string(alphabet_index) +
                                         " r_char:" + std::to_string(r_char) + " rlbwt_char:" + std::to_string(rlbwt_char) +
                                         " idx:" + std::to_string(idx) + " offset: " + std::to_string(offset) +
                                         " threshold_value: " + std::to_string(threshold_value));
    }

    if (offset >= threshold_value) {

        if (movi_options->is_debug())
            DEBUG_MSG("[reposition_thresholds] Repositioning down with thresholds..:");

        handle_reposition_down(idx, saved_idx, r_char, next_down, scan_count);

        if (r_char != alphabet[rlbwt[idx].get_c()]) {
            throw std::runtime_error(ERROR_MSG("[reposition thresholds] r_char != alphabet[rlbwt[idx].get_c()], r_char: " + std::to_string(r_char) +
                                               ", alphabet[rlbwt[idx].get_c()]: " + std::to_string(alphabet[rlbwt[idx].get_c()])));
        }
        return false;
    } else {

        if (movi_options->is_debug())
            DEBUG_MSG("[reposition_thresholds] Repositioning up with thresholds..");

        handle_reposition_up(idx, saved_idx, r_char, next_up, scan_count);

        if (r_char != alphabet[rlbwt[idx].get_c()]) {
            throw std::runtime_error(ERROR_MSG("[reposition thresholds] r_char != alphabet[rlbwt[idx].get_c()], r_char: " + std::to_string(r_char) +
                                               ", alphabet[rlbwt[idx].get_c()]: " + std::to_string(alphabet[rlbwt[idx].get_c()])));

        }
        return true;
    }
}
#endif

bool MoveStructure::reposition_randomly(uint64_t& idx, uint64_t& offset, char r_char, uint64_t& scan_count) {
    uint64_t saved_idx = idx;
    thread_local ThreadRandom random_generator;
    // uint16_t reposition_direction = random_generator.get_random() % 2;
    uint16_t reposition_direction =  offset * 2 < get_n(idx) ? 1 : 0;
    bool up = false;
    scan_count = 0;
    if (movi_options->is_debug())
        DEBUG_MSG("idx before repositioning: " + std::to_string(idx));


    // Selecting the right direction for trivial cases at the beginning and end of the move table
    if (idx == r - 1) {
        reposition_direction = 1;
    }
    if (idx == 0) {
        reposition_direction = 0;
    }

    if ( (reposition_direction == 1 && idx > 0) or idx == r - 1) {

        // repositioning up
        up = true;

        if (movi_options->is_debug())
            DEBUG_MSG("Repositioning up randomly:");

        idx = reposition_up(saved_idx, r_char, scan_count);
        if (movi_options->is_debug())
            DEBUG_MSG("idx after repositioning up: " + std::to_string(idx));

        if (idx >= r) {
            if (movi_options->is_debug())
                DEBUG_MSG("Up didn't work, try repositioning down:");

            // repositioning down
            up = false;
            idx = reposition_down(saved_idx, r_char, scan_count);
            if (movi_options->is_debug())
                DEBUG_MSG("idx after repositioning down: " + std::to_string(idx));
            if (idx == r) {
                // TODO
                throw std::runtime_error(ERROR_MSG("[reposition randomly] Neither up or down repositioning works.\n" +
                                    "The character does not exist in the index.\n"));
                throw std::runtime_error(ERROR_MSG("[reposition randomly] Neither up or down repositioning works.\n"));
            }
        }
    } else {

        // repositioning down
        up = false;

        if (movi_options->is_debug())
            DEBUG_MSG("Repositioning down randomly:");

        idx = reposition_down(saved_idx, r_char, scan_count);
        if (movi_options->is_debug())
            DEBUG_MSG("idx after repositioning down: " + std::to_string(idx));

        if (idx >= r) {
            if (movi_options->is_debug())
                DEBUG_MSG("Down didn't work, try repositioning up:");

            // repositioning up
            up = true;
            idx = reposition_up(saved_idx, r_char, scan_count);
            if (movi_options->is_debug())
                DEBUG_MSG("idx after repositioning up: " + std::to_string(idx));
            if (idx == r) {
                // TODO
                throw std::runtime_error(ERROR_MSG("[reposition randomly] Neither up or down repositioning works.\n" +
                                                   "The character does not exist in the index.\n"));
            }
        }
    }

    // sanity check
    char c = alphabet[rlbwt[idx].get_c()];
    if (c != r_char or idx == r) {
        throw std::runtime_error(ERROR_MSG("[reposition randomly] This should never happen.""c: " + c + "\n" +
                                           "\tidx: " + std::to_string(idx) + "\n"));
    }

    return up;
}

