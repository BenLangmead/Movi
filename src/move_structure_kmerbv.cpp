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
    // Dense-id bitvector B_k: one mark at each real k-mer's first row lb.
    std::unordered_map<uint32_t, sdsl::bit_vector> bvs;
    // Count bitvector C_k: marks at lb AND rb+1 (both boundaries of each real
    // k-mer's interval). A k-mer's interval [lb,rb] is a single group, so a
    // predecessor query from any row p in [lb,rb] returns lb and a successor
    // query returns rb+1, giving count = succ(p)-pred(p). This is what the
    // positive-skip count walk uses: at each backward-search step the current
    // (shrinking) range sits inside the leftmost k-mer's interval, and pred/succ
    // snap out to that k-mer's exact group boundaries. C_k is consumed only to
    // derive the compact run-local structure below; it is never serialized.
    // Size length+1 so the final group's rb+1 == length is representable.
    std::unordered_map<uint32_t, sdsl::bit_vector> cbvs;
    for (uint32_t k : ks) {
        max_k = std::max(max_k, k);
        bvs.emplace(k, sdsl::bit_vector(length, 0));
        cbvs.emplace(k, sdsl::bit_vector(length + 1, 0));
    }
    // Fast depth -> bitvector lookup (avoids a hash probe per visited node).
    std::vector<sdsl::bit_vector*> bv_at(max_k + 1, nullptr);
    std::vector<sdsl::bit_vector*> cbv_at(max_k + 1, nullptr);
    for (uint32_t k : ks) { bv_at[k] = &bvs.at(k); cbv_at[k] = &cbvs.at(k); }

    // The four ACGT characters (the '%' separator and sentinel are excluded:
    // real k-mers never span them).
    const char bases[4] = {'A', 'C', 'G', 'T'};

    auto set_bit = [](sdsl::bit_vector* bv, uint64_t pos) {
        // Atomic OR: different threads may touch the same 64-bit word (e.g. one
        // k-mer's rb+1 coincides with an adjacent k-mer's lb), so the
        // read-modify-write must be atomic; setting the same bit twice is fine.
        __atomic_fetch_or(&bv->data()[pos >> 6], 1ULL << (pos & 63ULL), __ATOMIC_RELAXED);
    };
    auto mark = [&](uint32_t depth, const MoveInterval& iv) {
        if (bv_at[depth] == nullptr) return;
        uint64_t lb = all_p[iv.run_start] + iv.offset_start;
        uint64_t rb = all_p[iv.run_end] + iv.offset_end;
        set_bit(bv_at[depth], lb);          // dense-id bitvector: lb only
        set_bit(cbv_at[depth], lb);         // count bitvector: lb ...
        set_bit(cbv_at[depth], rb + 1);     //              ... and rb+1
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
    // Serialize each requested k: the dense-id bitvector B_k (+rank) and the
    // run-local count structure derived from C_k.
    // ------------------------------------------------------------------ //
    std::string index_dir = movi_options->get_index_dir();
    for (uint32_t k : ks) {
        std::string ks_ = std::to_string(k);
        // Dense-id bitvector (lb only) + rank.
        sdsl::bit_vector& bv = bvs.at(k);
        sdsl::rank_support_v<> rank_bv(&bv);
        uint64_t ones = rank_bv(bv.size());
        sdsl::store_to_file(bv,      index_dir + "/kmerbv." + ks_ + ".bv");
        sdsl::store_to_file(rank_bv, index_dir + "/kmerbv." + ks_ + ".rank");

        // Run-local (all_p-free) count structure derived from C_k. Decompose the
        // group-starts relative to Movi's run heads (which come for free from the
        // move structure): for each row p,
        //   - an interior group-start (marked, not a run head) -> record its
        //     within-run offset (per-run, delta-coded) and flag its run in hi;
        //   - a run head that is NOT a group-start -> flag it in ex.
        // Counting then needs only run lengths + these per-run offsets: no
        // absolute positions and no n-space count bitvector. (Earlier revisions
        // also built plain/sd_vector/rrr_vector and a run-heads+Gi+E hybrid to map
        // the space/speed frontier; this run-local form was both smallest and
        // fastest, so it is the only representation kept.)
        sdsl::bit_vector& cbv = cbvs.at(k);
        sdsl::bit_vector ex(r, 0);   // ex[run] = 1 iff that run head is not a group-start
        sdsl::bit_vector hi(r, 0);   // hi[run] = 1 iff that run contains an interior group-start
        std::vector<uint64_t> gioff_vec;   // within-run offset of each interior, in (run,offset) order
        std::vector<bool>     gimark_vec;  // 1 at the FIRST interior of each run
        {
            uint64_t prev_run = std::numeric_limits<uint64_t>::max();
            uint64_t run = 0;
            for (uint64_t p = 0; p < length; ++p) {
                while (run + 1 < all_p.size() and all_p[run + 1] <= p) run++;
                bool is_head = (all_p[run] == p);
                if (cbv[p] && !is_head) {
                    hi[run] = 1;
                    gioff_vec.push_back(p - all_p[run]);
                    gimark_vec.push_back(run != prev_run);
                    prev_run = run;
                }
                if (is_head && !cbv[p]) ex[run] = 1;
            }
        }
        sdsl::store_to_file(hi, index_dir + "/kmerbv." + ks_ + ".cnt.hi");
        sdsl::store_to_file(ex, index_dir + "/kmerbv." + ks_ + ".cnt.ex");

        // Re-key each interior to the delta from the previous interior in its run
        // (absolute for the run's first interior), then variable-length code the
        // result with a dac_vector<2> (O(1) random read). Within-run offsets are
        // small and strictly increasing, so this is ~2.5x smaller than fixed-width
        // while staying exact.
        sdsl::int_vector<> gioff_dt(gioff_vec.size(), 0, 64);
        { uint64_t prev = 0; for (size_t i = 0; i < gioff_vec.size(); ++i) {
              uint64_t v = gioff_vec[i]; gioff_dt[i] = gimark_vec[i] ? v : (v - prev); prev = v; } }
        sdsl::util::bit_compress(gioff_dt);
        sdsl::dac_vector<2> gioff_dac(gioff_dt);
        sdsl::store_to_file(gioff_dac, index_dir + "/kmerbv." + ks_ + ".cnt.gioffdac");
        sdsl::bit_vector gimark(gimark_vec.size(), 0);
        for (size_t i = 0; i < gimark_vec.size(); ++i) if (gimark_vec[i]) gimark[i] = 1;
        sdsl::store_to_file(gimark, index_dir + "/kmerbv." + ks_ + ".cnt.gimark");

        INFO_MSG("k=" + ks_ + ": " + std::to_string(ones) +
                 " distinct k-mers; bitvectors written to " + index_dir +
                 "/kmerbv." + ks_ + ".{bv,cnt.*}");
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

    // Run-local (all_p-free) count structure: per-run interior offsets + the
    // run-head exception/has-interior flags. kmer_count_from_bv resolves a
    // k-mer's count from these plus run lengths alone.
    std::string pfx = index_dir + "/kmerbv." + std::to_string(k) + ".cnt";
    sdsl::load_from_file(kmerbv_ex,     pfx + ".ex");
    sdsl::load_from_file(kmerbv_hi,     pfx + ".hi");
    sdsl::load_from_file(kmerbv_gioff,  pfx + ".gioffdac");   // delta-transformed, dac_vector<2>
    sdsl::load_from_file(kmerbv_gimark, pfx + ".gimark");
    kmerbv_hi_rank = sdsl::rank_support_v<>(&kmerbv_hi);
    kmerbv_hi_ones = kmerbv_hi_rank(kmerbv_hi.size());
    kmerbv_gimark_sel = sdsl::select_support_mcl<1>(&kmerbv_gimark);

    // The MPHF id (--kmer --kmer-bv) forms lb = kmerbv_all_p[run] + offset to do
    // rank(lb) on B_k, so it needs a per-run absolute-position table; the
    // run-local count needs no absolute positions. Build the table only for id
    // queries (not --kmer-count) -- that is the all_p elimination on the count
    // path. The table is bit-packed (width ceil(log2 n)), not the 8-byte/run all_p.
    if (!movi_options->is_kmer_count()) {
        uint64_t w = (length <= 1) ? 1 : (sdsl::bits::hi(length) + 1);
        kmerbv_all_p = sdsl::int_vector<>(r + 1, 0, w);
        uint64_t acc = 0;
        for (uint64_t i = 0; i < r; ++i) { kmerbv_all_p[i] = acc; acc += get_n(i); }
        kmerbv_all_p[r] = length;
    }

    INFO_MSG("Loaded k-mer bitvector for k=" + std::to_string(k));
}
