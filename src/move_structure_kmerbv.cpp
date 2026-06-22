#include "move_structure.hpp"

#include <stack>
#include <unordered_map>

void MoveStructure::rebuild_all_p_if_needed() {
    if (!all_p.empty()) return;
    all_p.resize(r + 1, 0);
    for (uint64_t i = 0; i < r; ++i) {
        all_p[i + 1] = all_p[i] + get_n(i);
    }
}

void MoveStructure::build_kmerbv(const std::vector<uint32_t>& ks) {
    if (ks.empty()) return;

    // ------------------------------------------------------------------ //
    // B_k marks, in Movi's NATIVE BWT-row order, the first row of each
    // distinct real (separator-free) k-mer's interval.  rank(lb) over B_k
    // is then the dense lexicographic id of that k-mer among all distinct
    // real k-mers (an SSHash-style MPHF), and it is collision-free because
    // every distinct k-mer interval has a distinct first row.
    //
    // We obtain the native group-start rows directly from Movi's own RLBWT
    // (no second suffix-array build): a DFS that grows a k-mer one base at a
    // time to the left via backward_search_step, branching only over A/C/G/T
    // (never the '%' separator, which Movi treats as a per-record
    // terminator).  Every distinct d-mer present in the index is visited
    // exactly once as a depth-d node; we set B_d[lb] = 1 there.  Because the
    // bitvector is built and queried in the same native row order, there is
    // no ordering to disagree (the libsais SA-order approach failed exactly
    // here: libsais reads suffixes through '%' while Movi terminates at it).
    // ------------------------------------------------------------------ //
    rebuild_all_p_if_needed();

    uint32_t max_k = 0;
    std::unordered_map<uint32_t, sdsl::bit_vector> bvs;
    for (uint32_t k : ks) {
        max_k = std::max(max_k, k);
        bvs.emplace(k, sdsl::bit_vector(length, 0));
    }
    // Fast depth -> bitvector lookup (avoids a hash probe per visited node).
    std::vector<sdsl::bit_vector*> bv_at(max_k + 1, nullptr);
    for (uint32_t k : ks) bv_at[k] = &bvs.at(k);

    // The four ACGT characters (the '%' separator and sentinel are excluded:
    // real k-mers never span them).
    const char bases[4] = {'A', 'C', 'G', 'T'};

    // Atomically set bit `lb` of `bv` (different threads may touch the same
    // 64-bit word; distinct k-mers always have distinct lb, so no bit is set
    // by two threads, but the read-modify-write on the word must be atomic).
    auto mark = [&](uint32_t depth, const MoveInterval& iv) {
        sdsl::bit_vector* bv = bv_at[depth];
        if (bv == nullptr) return;
        uint64_t lb = all_p[iv.run_start] + iv.offset_start;
        __atomic_fetch_or(&bv->data()[lb >> 6], 1ULL << (lb & 63ULL), __ATOMIC_RELAXED);
    };

    struct Node { MoveInterval interval; uint32_t depth; };

    // ------------------------------------------------------------------ //
    // Phase 1 (serial): DFS down to `seed_depth`, marking every requested k
    // in [1, seed_depth-1] and collecting the depth-`seed_depth` intervals as
    // independent roots for the parallel phase.  Seeding is cheap (there are
    // at most 4^seed_depth distinct prefixes).
    // ------------------------------------------------------------------ //
    const uint32_t seed_depth = std::min<uint32_t>(max_k, 6);
    std::vector<MoveInterval> seeds;
    {
        std::stack<Node> stk;
        for (char c : bases) {
            if (!check_alphabet(c)) continue;
            uint64_t ci = alphamap[static_cast<uint64_t>(c)] + 1;
            MoveInterval iv(first_runs[ci], first_offsets[ci], last_runs[ci], last_offsets[ci]);
            if (iv.is_empty()) continue;
            stk.push({iv, 1});
        }
        while (!stk.empty()) {
            Node node = stk.top();
            stk.pop();
            if (node.depth == seed_depth) {
                seeds.push_back(node.interval);  // marked in the parallel phase
                continue;
            }
            mark(node.depth, node.interval);
            for (char c : bases) {
                MoveInterval child = node.interval;
                if (backward_search_step(c, child))  // prepend c (left-extend)
                    stk.push({child, node.depth + 1});
            }
        }
    }

    INFO_MSG("Enumerating distinct k-mers from the RLBWT (max_k=" +
             std::to_string(max_k) + ", " + std::to_string(seeds.size()) +
             " parallel seeds at depth " + std::to_string(seed_depth) + ")...");

    // ------------------------------------------------------------------ //
    // Phase 2 (parallel): each seed subtree is independent.  A thread DFSes
    // its seed from `seed_depth` to `max_k`, marking every requested k.
    // ------------------------------------------------------------------ //
    #pragma omp parallel for schedule(dynamic)
    for (size_t s = 0; s < seeds.size(); ++s) {
        std::stack<Node> stk;
        stk.push({seeds[s], seed_depth});
        while (!stk.empty()) {
            Node node = stk.top();
            stk.pop();
            mark(node.depth, node.interval);
            if (node.depth == max_k) continue;
            for (char c : bases) {
                MoveInterval child = node.interval;
                if (backward_search_step(c, child))
                    stk.push({child, node.depth + 1});
            }
        }
    }

    // ------------------------------------------------------------------ //
    // Build rank support and serialize each requested k.
    // ------------------------------------------------------------------ //
    std::string index_dir = movi_options->get_index_dir();
    for (uint32_t k : ks) {
        sdsl::bit_vector& bv = bvs.at(k);
        sdsl::rank_support_v<> rank_bv(&bv);
        uint64_t ones = rank_bv(bv.size());

        std::string bv_path   = index_dir + "/kmerbv." + std::to_string(k) + ".bv";
        std::string rank_path = index_dir + "/kmerbv." + std::to_string(k) + ".rank";
        sdsl::store_to_file(bv,      bv_path);
        sdsl::store_to_file(rank_bv, rank_path);

        INFO_MSG("k=" + std::to_string(k) + ": " + std::to_string(ones) +
                 " distinct k-mers; bitvector written to " + bv_path);
    }
}

void MoveStructure::load_kmerbv(uint32_t k) {
    std::string index_dir = movi_options->get_index_dir();
    std::string bv_path   = index_dir + "/kmerbv." + std::to_string(k) + ".bv";
    std::string rank_path = index_dir + "/kmerbv." + std::to_string(k) + ".rank";

    if (!std::filesystem::exists(bv_path)) {
        throw std::runtime_error(
            ERROR_MSG("K-mer bitvector for k=" + std::to_string(k) +
                      " not found at " + bv_path +
                      ". Run 'movi build-kmerbv --index <dir> --kmer-lengths " +
                      std::to_string(k) + "' first."));
    }

    sdsl::load_from_file(kmerbv,      bv_path);
    sdsl::load_from_file(kmerbv_rank, rank_path);
    kmerbv_rank.set_vector(&kmerbv);

    // Rebuild all_p from the RLBWT so we can compute absolute BWT row numbers
    // from (run_start, offset_start) intervals during queries.
    rebuild_all_p_if_needed();

    INFO_MSG("Loaded k-mer bitvector for k=" + std::to_string(k));
}
