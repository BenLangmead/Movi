#include "move_structure.hpp"

uint64_t MoveStructure::query_zml(MoveQuery& mq) {
    auto& query_seq = mq.query();
    int32_t pos_on_r = query_seq.length() - 1;
    uint64_t match_len = 0;
    uint64_t ff_count_tot = 0;

    while (!check_alphabet(query_seq[pos_on_r]) and pos_on_r >= 0) {
        mq.add_ml(0, movi_options->is_stdout());
        pos_on_r -= 1;
    }

    if (pos_on_r < 0) {
        // Special case where no character in the read exists in the index.
        return 0;
    }

    // Multi-class classification
    if (movi_options->is_multi_classify()) {
        for (uint16_t i = 0; i < num_species; i++) {
            classify_cnts[i] = 0;
        }
    }

    MoveInterval interval = initialize_backward_search(mq, pos_on_r, match_len);
    while (pos_on_r > 0) {
        MoveInterval prev_interval = interval;
        ff_count_tot += backward_search_step(query_seq, pos_on_r, interval);
        if (!interval.is_empty()) {
            mq.add_ml(match_len, movi_options->is_stdout());
            pos_on_r -= 1;
            match_len += 1;
        } else {
            // Classification based on maximal matching
            if (movi_options->is_multi_classify() && match_len >= movi_options->get_min_match_len()) {
                if (prev_interval.run_start == prev_interval.run_end) {
                    for (uint64_t i = prev_interval.offset_start; i <= prev_interval.offset_end; i++) {
                        uint64_t full_ind = run_offsets[prev_interval.run_start] + i;
                        uint16_t cur_doc = doc_pats[full_ind];
                        classify_cnts[cur_doc] += match_len;
                    }   
                } else {
                    for (uint64_t i = prev_interval.offset_start; i < get_n(prev_interval.run_start); i++) {
                        uint64_t full_ind = run_offsets[prev_interval.run_start] + i;
                        uint16_t cur_doc = doc_pats[full_ind];
                        classify_cnts[cur_doc] += match_len;
                    }    
                    for (uint64_t r_ind = prev_interval.run_start + 1; r_ind < prev_interval.run_end; r_ind++) {
                        for (uint64_t i = 0; i < get_n(r_ind); i++) {
                            uint64_t full_ind = run_offsets[r_ind] + i;
                            uint16_t cur_doc = doc_pats[full_ind];
                            classify_cnts[cur_doc] += match_len;
                        }  
                    } 
                    for (uint64_t i = 0; i <= prev_interval.offset_end; i++) {
                        uint64_t full_ind = run_offsets[prev_interval.run_end] + i;
                        uint16_t cur_doc = doc_pats[full_ind];
                        classify_cnts[cur_doc] += match_len;
                    }           
                }
            }

            mq.add_ml(match_len, movi_options->is_stdout());
            pos_on_r -= 1;
            match_len = 0;
            while (!check_alphabet(query_seq[pos_on_r]) and pos_on_r > 0) {
                mq.add_ml(match_len, movi_options->is_stdout());
                pos_on_r -= 1;
            }
            // Special case where the character at position 0 of the read does not exist in the index.
            if (check_alphabet(query_seq[pos_on_r]))
                interval = initialize_backward_search(mq, pos_on_r, match_len);
        }
    }
    if (interval.is_empty()) {
        match_len = 0;
    }
    mq.add_ml(match_len, movi_options->is_stdout());

    // Document occuring the most is the genotype we think the query is from.
    if (movi_options->is_multi_classify()) {
        uint16_t best_doc = 0;
        for (uint16_t i = 1; i < num_species; i++) {
            if (classify_cnts[i] > classify_cnts[best_doc]) {
                best_doc = i;
            }
        }
        
        // Document occuring the most is the genotype we think the query is from.
        output_files->out_file << to_taxon_id[best_doc] << " ";
        for (uint16_t i = 0; i < num_species; i++) {
        //    output_files->out_file << classify_cnts[i] << " ";
        }
        output_files->out_file << "\n";
    }

    return ff_count_tot;
}
