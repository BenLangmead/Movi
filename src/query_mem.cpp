#include <sys/stat.h> 

#include "move_structure.hpp"
#include "utils.hpp"

// IMPORTANT: For the bidirectional search, the index has to be built with the reverse complement of the reference
uint64_t MoveStructure::query_mems(MoveQuery& mq) {
    int32_t min_mem_length = movi_options->get_min_mem_length();
    if (min_mem_length <= 1) return query_all_mems(mq);

    int32_t pos_on_r = 0; // MEM finding goes left to right
    std::string& query_seq = mq.query();

    size_t ftab_k = movi_options->get_ftab_k();
    // At the time, the initialization works if ftab with ftab_k exists only
    // And the multi-ftab strategy must be turned off
    movi_options->set_multi_ftab(false);
    
    uint64_t mems_found = 0;
    do {
        mems_found += query_mem_bml(mq, pos_on_r, min_mem_length, query_seq, ftab_k);
    } while (pos_on_r < query_seq.length());

    return mems_found;
}

// We search for a MEM in P[pos_on_r..m-1] starting at pos_on_r
// Returns true if a MEM is found, false otherwise, and updates pos_on_r
bool MoveStructure::query_mem_bml(MoveQuery& mq, int32_t& pos_on_r, int32_t& min_mem_length, std::string& query_seq, size_t& ftab_k) {
    // Reached end of query sequence
    if (pos_on_r + min_mem_length > query_seq.length()) {
        pos_on_r = query_seq.length();
        return false;
    }

    uint64_t match_len = 0;
    int32_t init_pos = pos_on_r + min_mem_length - 1;
    MoveBiInterval bi_interval;
    bi_interval = initialize_bidirectional_search(mq, init_pos, match_len);
    // If min_mem_length is at least ftab_k, we skip using BML
    bool ftab_skip = match_len <= 1 && ftab_k <= min_mem_length;
    --init_pos; // Move to next character left of initial match range

    if (ftab_skip) {
        // If BML skip is known from ftab, find next left end with only backward extension (extend_left slow without ftab)
        MoveInterval fw_interval = bi_interval.fw_interval;
        // Signed counter: the bound goes negative when the ftab covers the whole seed,
        // and a size_t one would wrap and loop forever.
        for (int32_t j = 0; j <= init_pos - pos_on_r; ++j) {
            if (!backward_search_step(query_seq[init_pos - j], fw_interval)) {
                pos_on_r = init_pos - j + 1;
                return false;
            }
            ++match_len;
        }
        // Should never reach here
        throw std::runtime_error("Extended past failed ftab");
    }
    else {
        // Backward extension to find if left end permits a sufficiently long MEM, should equal pos_on_r if true.
        // Signed counter, same reason as above.
        for (int32_t j = 0; j <= init_pos - pos_on_r; ++j) {
            if (!extend_left(query_seq[init_pos - j], bi_interval)) {
                pos_on_r = init_pos - j + 1;
                return false;
            }
            ++match_len;
        }
    }

    // Forward extension to find right end of MEM (exclusive)
    MoveInterval rc_interval = bi_interval.rc_interval; // We don't need to keep the fw_interval
    MoveInterval rc_interval_before_extension = rc_interval;
    size_t i;
    for (i = pos_on_r + min_mem_length; i < query_seq.length(); ++i) {
        rc_interval_before_extension = rc_interval;
        if (!forward_search_step(query_seq[i], rc_interval)) {
            rc_interval = rc_interval_before_extension;
            break;
        }
        ++match_len;
    }

    mq.add_mem(pos_on_r, i, rc_interval.count(rlbwt));

    // Backward extension to find next left end of candidate MEM (inclusive)
    // TODO: if this extension is at least L, we can skip the first BML step in next recursion
    size_t end_pos_on_r = i;
    size_t j;
    if (end_pos_on_r < query_seq.length()) {
        init_pos = end_pos_on_r;
        match_len = 0;
        MoveInterval fw_interval = initialize_backward_search(mq, init_pos, match_len);
        ++match_len; // TODO: this is still off by one for initialize_backward_search
        --init_pos; // Move to next character left of ftab k-mer or character left of initial range if ftab fails
        // Signed comparison, same reason as above.
        for (i = 0; static_cast<int32_t>(i) <= init_pos - (pos_on_r + 1); ++i) {
            if (!backward_search_step(query_seq[init_pos - i], fw_interval)) {
                break;
            }
            ++match_len;
        }
    }

    pos_on_r = (init_pos - i) + 1;
    return true;
}

uint64_t MoveStructure::query_all_mems(MoveQuery& mq) {
    std::string& query_seq = mq.query();
    const size_t m = query_seq.length();

    // At the time, the initialization works if ftab with ftab_k exists only
    // And the multi-ftab strategy must be turned off
    movi_options->set_multi_ftab(false);

    // bi_interval holds the occurrences of Q[s, s + match_len). match_len == 0 means no
    // match is held, and the next MEM is seeded from Q[s] alone.
    size_t s = 0;
    uint64_t match_len = 0;
    MoveBiInterval bi_interval;

    while (s < m) {
        if (match_len == 0) {
            // A character outside the alphabet is in no MEM.
            if (!check_alphabet(query_seq[s])) {
                ++s;
                continue;
            }
            int32_t init_pos = s;
            bi_interval = initialize_bidirectional_search(mq, init_pos, match_len);
            // s is 0 or follows an illegal character, so no ftab k-mer ending at s is
            // usable and the seed is Q[s] alone; init_pos is taken back in case the
            // seed ever covers characters left of s.
            s = init_pos;
        }

        // Extend right until mismatch. A failed extension leaves the intervals
        // modified, so the last good pair is restored.
        while (s + match_len < m) {
            MoveBiInterval bi_interval_before_extension = bi_interval;
            if (!extend_right(query_seq[s + match_len], bi_interval)) {
                bi_interval = bi_interval_before_extension;
                break;
            }
            ++match_len;
        }
        size_t e = s + match_len;
        mq.add_mem(s, e, bi_interval.fw_interval.count(rlbwt));
        if (e >= m) break;

        // The next MEM is the longest match ending at e. It starts after s, because
        // Q[s, e] does not occur (the MEM above stopped short of e).
        if (!check_alphabet(query_seq[e])) {
            s = e + 1;
            match_len = 0;
            continue;
        }
        int32_t left = e;
        uint64_t seed_len = 0;
        bi_interval = initialize_bidirectional_search(mq, left, seed_len);
        while (left - 1 > static_cast<int64_t>(s)) {
            MoveBiInterval bi_interval_before_extension = bi_interval;
            if (!extend_left(query_seq[left - 1], bi_interval)) {
                bi_interval = bi_interval_before_extension;
                break;
            }
            --left;
        }
        s = left;
        match_len = e - left + 1;
    }

    return mq.get_mems().size();
}
