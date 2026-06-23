#include "move_structure.hpp"

#include <stack>
#include <unordered_map>
#include <fstream>
#include <cstdlib>

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
    // bases[] is in A<C<G<T order, so its index doubles as the 2-bit code; the
    // complement of code x is 3-x (A<->T, C<->G).
    const char bases[4] = {'A', 'C', 'G', 'T'};

    // Canonical-id mode (SSHash-compatible): mark the dense-id bitvector B_k only
    // at canonical k-mers (s == min(s, rc(s))), so rank(lb) is a dense id in
    // [0, #canonical k-mers) matching SSHash's num_kmers. The count bitvector C_k
    // is ALWAYS marked for every real k-mer: counting is orientation-agnostic (a
    // query searches the literal k-mer and reads its interval size = the canonical
    // count occ(x)+occ(rc(x))), so the non-canonical orientation's interval must
    // still be resolvable.
    const bool canonical = movi_options->is_kmerbv_canonical();

    // 2-bit k-mer packed as it is grown left-to-right down the DFS: each
    // left-extension prepends the new (text-leftmost) base into a higher 2-bit
    // group. group j (bits [2j,2j+2)) holds the (j+1)-th base added, i.e. text
    // position depth-1-j, so the integer value orders lexicographically by the
    // real text and `packed <= rc(packed)` tests canonicality. k<=63 fits 126 bits.
    auto kmer_rc = [](unsigned __int128 x, uint32_t d) -> unsigned __int128 {
        unsigned __int128 r = 0;
        for (uint32_t i = 0; i < d; ++i) { r = (r << 2) | (3u - (unsigned)(x & 3u)); x >>= 2; }
        return r;
    };

    auto set_bit = [](sdsl::bit_vector* bv, uint64_t pos) {
        // Atomic OR: different threads may touch the same 64-bit word (e.g. one
        // k-mer's rb+1 coincides with an adjacent k-mer's lb), so the
        // read-modify-write must be atomic; setting the same bit twice is fine.
        __atomic_fetch_or(&bv->data()[pos >> 6], 1ULL << (pos & 63ULL), __ATOMIC_RELAXED);
    };
    auto mark = [&](uint32_t depth, const MoveInterval& iv, unsigned __int128 packed) {
        if (bv_at[depth] == nullptr) return;
        uint64_t lb = all_p[iv.run_start] + iv.offset_start;
        uint64_t rb = all_p[iv.run_end] + iv.offset_end;
        if (!canonical || packed <= kmer_rc(packed, depth))
            set_bit(bv_at[depth], lb);      // dense-id bitvector: lb only (canonical k-mers)
        set_bit(cbv_at[depth], lb);         // count bitvector: lb ...
        set_bit(cbv_at[depth], rb + 1);     //              ... and rb+1 (every real k-mer)
    };

    struct Node { MoveInterval interval; uint32_t depth; unsigned __int128 packed; };

    // ------------------------------------------------------------------ //
    // Phase 1 (serial): DFS down to `seed_depth`, marking every requested k
    // in [1, seed_depth-1] and collecting the depth-`seed_depth` intervals as
    // independent roots for the parallel phase.  Seeding is cheap (there are
    // at most 4^seed_depth distinct prefixes).
    // ------------------------------------------------------------------ //
    const uint32_t seed_depth = std::min<uint32_t>(max_k, 6);
    struct Seed { MoveInterval interval; unsigned __int128 packed; };
    std::vector<Seed> seeds;
    {
        std::stack<Node> stk;
        for (int i = 0; i < 4; ++i) {
            char c = bases[i];
            if (!check_alphabet(c)) continue;
            uint64_t ci = alphamap[static_cast<uint64_t>(c)] + 1;
            MoveInterval iv(first_runs[ci], first_offsets[ci], last_runs[ci], last_offsets[ci]);
            if (iv.is_empty()) continue;
            stk.push({iv, 1, (unsigned __int128)i});  // group 0 = first (rightmost) base
        }
        while (!stk.empty()) {
            Node node = stk.top();
            stk.pop();
            if (node.depth == seed_depth) {
                seeds.push_back({node.interval, node.packed});  // marked in the parallel phase
                continue;
            }
            mark(node.depth, node.interval, node.packed);
            for (int i = 0; i < 4; ++i) {
                MoveInterval child = node.interval;
                if (backward_search_step(bases[i], child))  // prepend base (left-extend)
                    stk.push({child, node.depth + 1,
                              node.packed | ((unsigned __int128)i << (2 * node.depth))});
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
        stk.push({seeds[s].interval, seed_depth, seeds[s].packed});
        while (!stk.empty()) {
            Node node = stk.top();
            stk.pop();
            mark(node.depth, node.interval, node.packed);
            if (node.depth == max_k) continue;
            for (int i = 0; i < 4; ++i) {
                MoveInterval child = node.interval;
                if (backward_search_step(bases[i], child))
                    stk.push({child, node.depth + 1,
                              node.packed | ((unsigned __int128)i << (2 * node.depth))});
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
        // MPHF-id structure for B_k. The sparse Elias-Fano sd_vector is the default
        // (far smaller than dense bv+rank since B_k is sparse, and competitive with
        // SSHash). The dense bit_vector+rank is ~6% faster at query but ~4-5x larger;
        // it is only written with --id-bv-dense (for benchmarking / max speed).
        sdsl::bit_vector& bv = bvs.at(k);
        sdsl::rank_support_v<> rank_bv(&bv);
        uint64_t ones = rank_bv(bv.size());          // = num_kmers
        // Sparse (Elias-Fano) B_k over ABSOLUTE BWT rows: universe = n (length),
        // each mark set at its absolute row p. This is the smallest on-disk rep; the
        // per-run base position needed to address it at query (lb = all_p[run]+offset)
        // is reconstructed from sampled checkpoints + run lengths at load time (see
        // reconstruct_allp / load_kmerbv), so the full all_p table is never stored.
        // Built via sd_vector_builder since p is strictly increasing.
        sdsl::sd_vector_builder builder(length, ones);
        for (uint64_t p = 0; p < length; ++p)
            if (bv[p]) builder.set(p);
        sdsl::sd_vector<> sdb(builder);
        sdsl::store_to_file(sdb, index_dir + "/kmerbv." + ks_ + ".sd");
        uint64_t sz_sd    = sdsl::size_in_bytes(sdb);
        uint64_t sz_dense = sdsl::size_in_bytes(bv) + sdsl::size_in_bytes(rank_bv);
        if (movi_options->is_kmerbv_id_dense()) {
            // Dense rep keeps absolute-row addressing (needs all_p at query).
            sdsl::store_to_file(bv,      index_dir + "/kmerbv." + ks_ + ".bv");
            sdsl::store_to_file(rank_bv, index_dir + "/kmerbv." + ks_ + ".rank");
        }

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

        // Sidecar metadata so the query / access side knows the id space:
        //   line 1: canonical flag (0|1);  line 2: num_kmers (= ones).
        {
            std::ofstream meta(index_dir + "/kmerbv." + ks_ + ".meta");
            meta << (canonical ? 1 : 0) << "\n" << ones << "\n";
        }

        INFO_MSG("k=" + ks_ + ": num_kmers=" + std::to_string(ones) + " distinct " +
                 (canonical ? "CANONICAL " : "") + "k-mers; B_k id sd_vector=" +
                 std::to_string(sz_sd / 1048576) + "MB (" +
                 std::to_string(sz_dense ? (100 * sz_sd / sz_dense) : 0) +
                 "% of dense bv+rank=" + std::to_string(sz_dense / 1048576) + "MB" +
                 (movi_options->is_kmerbv_id_dense() ? ", dense also written" : "") +
                 "); k-mer files at " + index_dir + "/kmerbv." + ks_ + ".*");
    }
}

// Full backward search of k-mer string s over the (doubled) BWT; returns s's
// BWT interval, or an empty interval if s does not occur. Mirrors the per-base
// left-extension the build DFS uses. Used to locate rc(x)'s interval for the
// canonical-id lookup (id = rank(min(lb(x), lb(rc(x))))).
MoveInterval MoveStructure::search_kmer_interval(std::string s) {
    MoveInterval iv;
    int n = static_cast<int>(s.size());
    if (n == 0 || !check_alphabet(s[n - 1])) { iv.make_empty(); return iv; }
    uint64_t ci = alphamap[static_cast<uint64_t>(s[n - 1])] + 1;
    iv = MoveInterval(first_runs[ci], first_offsets[ci], last_runs[ci], last_offsets[ci]);
    for (int i = n - 2; i >= 0; --i) {
        if (!backward_search_step(s[i], iv)) { iv.make_empty(); return iv; }
    }
    return iv;
}

void MoveStructure::load_kmerbv(uint32_t k) {
    std::string index_dir = movi_options->get_index_dir();
    std::string ks_       = std::to_string(k);
    std::string bv_path   = index_dir + "/kmerbv." + ks_ + ".bv";
    std::string rank_path = index_dir + "/kmerbv." + ks_ + ".rank";
    std::string sd_path   = index_dir + "/kmerbv." + ks_ + ".sd";

    // MPHF-id B_k representation. Default: sparse Elias-Fano sd_vector (small,
    // SSHash-competitive). MOVI_ID_BV=dense forces the dense bit_vector+rank rep
    // (~6% faster query, ~4-5x larger; only present if built with --id-bv-dense).
    const char* idbv_env = std::getenv("MOVI_ID_BV");
    bool want_dense = (idbv_env != nullptr && std::string(idbv_env) == "dense");
    bool sd_ok = std::filesystem::exists(sd_path);
    bool dense_ok = std::filesystem::exists(bv_path);
    if (!sd_ok && !dense_ok) {
        throw std::runtime_error(
            ERROR_MSG("K-mer id bitvector for k=" + ks_ + " not found at " + sd_path +
                      ". Run 'movi build-kmerbv --index <dir> --kmer-lengths " + ks_ + "' first."));
    }
    kmerbv_use_sd = sd_ok && !(want_dense && dense_ok);
    if (kmerbv_use_sd) {
        sdsl::load_from_file(kmerbv_sd, sd_path);
        kmerbv_sd_rank = sdsl::sd_vector<>::rank_1_type(&kmerbv_sd);
        kmerbv_sd_sel  = sdsl::sd_vector<>::select_1_type(&kmerbv_sd);
    } else {
        sdsl::load_from_file(kmerbv,      bv_path);
        sdsl::load_from_file(kmerbv_rank, rank_path);
        kmerbv_rank.set_vector(&kmerbv);
    }

    // Read the id-space metadata (canonical flag + num_kmers). In canonical mode
    // the query canonicalizes each k-mer before the id lookup. (Older indexes wrote
    // a 3rd line, a now-unused sd stride; it is simply ignored.)
    kmerbv_is_canonical = false; kmerbv_num_kmers = 0;
    std::ifstream meta(index_dir + "/kmerbv." + ks_ + ".meta");
    if (meta) { int c = 0; meta >> c >> kmerbv_num_kmers; kmerbv_is_canonical = (c == 1); }

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

    // Sampled per-run absolute-position checkpoints for the MPHF id query
    // (--kmer --kmer-bv). Both the sd and dense B_k reps address ABSOLUTE BWT rows
    // (lb = all_p[run]+offset), so instead of storing all r+1 positions we keep one
    // checkpoint every kmerbv_allp_S runs and reconstruct exact positions from the
    // move structure's run lengths at query (reconstruct_allp). The stride defaults
    // to 16 and is overridable via MOVI_ALLP_SAMPLE for space/speed benchmarking.
    // --kmer-count is run-local and needs none.
    if (!movi_options->is_kmer_count()) {
        kmerbv_allp_S = 16;
        if (const char* s = std::getenv("MOVI_ALLP_SAMPLE")) {
            uint64_t v = std::strtoull(s, nullptr, 10);
            if (v >= 1) kmerbv_allp_S = v;
        }
        uint64_t w = (length <= 1) ? 1 : (sdsl::bits::hi(length) + 1);
        uint64_t nck = r / kmerbv_allp_S + 1;
        kmerbv_allp_ckpt = sdsl::int_vector<>(nck, 0, w);
        uint64_t acc = 0;
        for (uint64_t i = 0; i < r; ++i) {
            if (i % kmerbv_allp_S == 0) kmerbv_allp_ckpt[i / kmerbv_allp_S] = acc;
            acc += get_n(i);
        }
    }

    INFO_MSG("Loaded k-mer bitvector for k=" + std::to_string(k));
}
