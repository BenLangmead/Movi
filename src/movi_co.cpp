// ============================================================================
// movi-co: coroutine-based latency hiding (query implementation "style (c)")
//
// Movi runs each hot query in up to three styles; this file is style (c):
//   (a) sequential, no latency hiding:
//         move_structure_query.cpp (query_pml/query_zml),
//         query_kmer.cpp (query_all_kmers/query_kmers_from),
//         query_mem.cpp (query_mems/query_mem_bml).
//   (b) manual "strand" state machine: prefetch, then hand-switch among
//         in-flight reads at each prefetch point:
//         read_processor.cpp (process_latency_hiding, kmer_search_latency_hiding),
//         struct Strand in read_processor.hpp.
//   (c) coroutines (this file): prefetch, then co_yield at each prefetch point;
//         the scheduler resumes another coroutine while the cache line fills.
//
// The pattern in every query_*_coroutine below: compute the LF destination row
// id (step_prep), prefetch it, co_yield, then dereference it on resume
// (step_finish). Only long-range non-sequential move-structure accesses are
// hidden; sequential run scans rely on the hardware prefetcher.
// ============================================================================

/**
 * Author: Ben Langmead
 * Date: Oct 4, 2025
 *
 * Simple program to compute Pseudo Matching Lengths (PMLs) for reads.  Meant
 * to be a starting point for work on using co-routines.
 *
 * Distributed under the GPL3 license.
 * See accompanying LICENSE or https://opensource.org/license/gpl-3-0
 */

#include <cstdint>
#include <stdio.h>
#include <cstdio>
#include <chrono>
#include <cstddef>
#include <unistd.h>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <algorithm>
#include <variant>
#include <stdexcept>
#include <mutex>
#include <cstring>
#include <map>

#include <sdsl/int_vector.hpp>
#include <zlib.h>

#include "utils.hpp"
#include "move_structure.hpp"
#include "move_query.hpp"
#include "movi_options.hpp"
#include <coroutine>

using std::suspend_always;
using std::coroutine_handle;
using std::string;
using std::numeric_limits;
using std::monostate;
using std::exception;
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

// Global debug flag
static bool debug_enabled = false;
// When set, write binary BPF output (like mainline Movi) instead of text.
static bool g_bpf_output = false;
// When true (--reorder), emit reads in input order via the reorder buffer;
// default is completion order (lower overhead, matches mainline Movi's order).
static bool g_ordered_output = false;
// MEM query mode (movi-co also supports --mem) and its parameters.
static bool g_mem_query = false;
static int  g_ftab_k = 0;
static int  g_min_mem_length = 0;
// k-mer query mode (movi-co also supports --kmer / --kmer-count) and its length.
static bool g_kmer_query = false;
static bool g_kmer_count = false;
static bool g_kmer_bv = false;   // use the count bitvector (pred/succ) for counts
static int  g_k = 0;

// Setup buffered stdout and stderr
static FILE* stdout_buf = nullptr;
static FILE* stderr_buf = nullptr;

// Initialize buffered I/O
static void init_buffered_io() {
    stdout_buf = stdout;
    stderr_buf = stderr;
    // Use full buffering for stdout (better performance)
    setvbuf(stdout_buf, nullptr, _IOFBF, 8192);
    // Use line buffering for stderr (immediate error visibility)
    setvbuf(stderr_buf, nullptr, _IOLBF, 0);
}

// Custom function to write unsigned integer to buffer, digit by digit
// Returns number of bytes written
static size_t write_uint_to_buffer(char* buffer, size_t buffer_size, uint64_t value) {
    if (buffer_size == 0) return 0;
    
    // Handle zero case
    if (value == 0) {
        if (buffer_size > 0) {
            buffer[0] = '0';
            return 1;
        }
        return 0;
    }
    
    // Calculate number of digits
    uint64_t temp = value;
    int num_digits = 0;
    while (temp > 0) {
        num_digits++;
        temp /= 10;
    }
    
    // Write digits from least to most significant, then reverse
    if (static_cast<size_t>(num_digits) > buffer_size) {
        num_digits = static_cast<int>(buffer_size);
    }
    
    char* digits = buffer;
    temp = value;
    for (int i = 0; i < num_digits; i++) {
        digits[i] = '0' + (temp % 10);
        temp /= 10;
    }
    
    // Reverse the digits
    for (int i = 0; i < num_digits / 2; i++) {
        char tmp = digits[i];
        digits[i] = digits[num_digits - 1 - i];
        digits[num_digits - 1 - i] = tmp;
    }
    
    return static_cast<size_t>(num_digits);
}

// Custom function to write string to buffer
static size_t write_string_to_buffer(char* buffer, size_t buffer_size, const char* str, size_t str_len) {
    size_t copy_len = (str_len < buffer_size) ? str_len : buffer_size;
    if (copy_len > 0) {
        std::memcpy(buffer, str, copy_len);
    }
    return copy_len;
}

// DEBUG macro that only prints when NDEBUG is not defined and debug flag is set
// Uses fprintf-style format strings
#ifndef NDEBUG
#define DEBUG_MSG_CO(...) do { if (debug_enabled) { fprintf(stderr_buf, __VA_ARGS__); } } while(0)
#else
#define DEBUG_MSG_CO(...) ((void)0)
#endif

// Mutex for thread-safe output from concurrent coroutines
static std::mutex output_mutex;

// Flush interval for periodic flushing
static constexpr int FLUSH_INTERVAL = 100;  // Flush every 100 lines

// In-order, chunked output emitter. Coroutines finish reads out of order; this
// buffers completed lines and writes them in input (read) order, so downstream
// tools see the same order as the input. Single-threaded (the scheduler resumes
// one coroutine at a time), so no locking is needed.
static constexpr size_t OUTPUT_BUFFER_SIZE = 65536;
struct OrderedEmitter {
    uint64_t next_emit = 0;
    bool ordered = false;         // when true, buffer + write in input order; else completion order
    std::map<uint64_t, std::string> pending;  // completed lines awaiting their turn
    char buf[OUTPUT_BUFFER_SIZE];
    size_t pos = 0;

    void flush_chunk() {
        if (pos > 0) { fwrite(buf, 1, pos, stdout_buf); pos = 0; }
    }
    void write_line(const std::string& line) {
        if (line.size() > OUTPUT_BUFFER_SIZE) {
            flush_chunk();
            fwrite(line.data(), 1, line.size(), stdout_buf);
        } else {
            if (pos + line.size() > OUTPUT_BUFFER_SIZE) flush_chunk();
            std::memcpy(buf + pos, line.data(), line.size());
            pos += line.size();
        }
    }
    // Emit the line for read 'seq'. Writes immediately if it is the next read in
    // input order, then drains any consecutive buffered lines; otherwise buffers
    // a copy until the earlier reads complete.
    void emit(uint64_t seq, const std::string& line) {
        if (!ordered) { write_line(line); return; }  // completion order (reorder disabled)
        if (seq == next_emit) {
            write_line(line);
            ++next_emit;
            auto it = pending.find(next_emit);
            while (it != pending.end()) {
                write_line(it->second);
                pending.erase(it);
                ++next_emit;
                it = pending.find(next_emit);
            }
        } else {
            pending.emplace(seq, line);  // copy: only for out-of-order completions
        }
    }
    void finish() {
        for (auto& kv : pending) write_line(kv.second);
        pending.clear();
        flush_chunk();
        fflush(stdout_buf);
    }
};
static OrderedEmitter g_emitter;

// Shared reader that parses FASTQ/A records in batches and hands them to
// coroutines through per-coroutine slots. The hot path neither parses one
// record at a time nor copies each read: the scheduler moves a record from the
// current batch into the awaiting coroutine's slot, and the coroutine uses it
// in place (zero per-read copy). Single-threaded -- no background I/O thread --
// so the query still runs on exactly one CPU.
class SharedFastqReader {
public:
    struct ReadData {
        bool valid = false;
        string name;
        string sequence;
        uint64_t seq = 0;          // input-order index, assigned when served
    };

private:
    gzFile fp;
    kseq_t* seq;
    std::vector<ReadData> batch;   // current batch of parsed records
    size_t batch_pos = 0;          // next unserved record within batch
    std::vector<ReadData> slots;   // one slot per coroutine (zero-copy handoff)
    uint64_t served_ = 0;          // monotonic read index, for in-order output
    static constexpr size_t BATCH_N = 256;

    // Parse up to BATCH_N records in a tight loop. This amortizes the kseq
    // parse / gzip-decompression call overhead and keeps it off the per-read
    // critical path (one refill stall every BATCH_N reads instead of a parse
    // interleaved with every read's compute).
    void refill() {
        batch.clear();
        batch_pos = 0;
        for (size_t k = 0; k < BATCH_N; k++) {
            int l = kseq_read(seq);
            if (l < 0) break;  // EOF (-1) or error (-2/-3)
            ReadData rd;
            rd.valid = true;
            rd.name.assign(seq->name.s, seq->name.l);
            rd.sequence.assign(seq->seq.s, seq->seq.l);
            batch.push_back(std::move(rd));
        }
    }

public:
    SharedFastqReader(const string& fastq_file, int concurrency) : fp(nullptr), seq(nullptr) {
        fp = gzopen(fastq_file.c_str(), "r");
        if (!fp) {
            throw std::runtime_error("Cannot open FASTQ file: " + fastq_file);
        }
        seq = kseq_init(fp);
        if (!seq) {
            gzclose(fp);
            throw std::runtime_error("Cannot initialize kseq");
        }
        slots.resize(concurrency);
    }

    ~SharedFastqReader() {
        if (seq) kseq_destroy(seq);
        if (fp) gzclose(fp);
    }

    // Move the next record into coroutine id's slot. Returns false at EOF
    // (the slot is marked invalid).
    bool serve_next(int id) {
        if (batch_pos >= batch.size()) refill();
        if (batch_pos >= batch.size()) { slots[id].valid = false; return false; }
        slots[id] = std::move(batch[batch_pos++]);
        slots[id].seq = served_++;
        return true;
    }

    ReadData& slot(int id) { return slots[id]; }
};

// Awaitable: suspends the coroutine until the scheduler serves a read into this
// coroutine's slot, then hands back a reference to that slot (no copy).
struct read_awaitable {
    SharedFastqReader* reader;
    coroutine_handle<>* stored_handle;  // Where to store the coroutine handle
    int coroutine_id;

    bool await_ready() const noexcept { return false; }

    template<typename Promise>
    void await_suspend(coroutine_handle<Promise> h) {
        *stored_handle = h;
    }

    SharedFastqReader::ReadData& await_resume() {
        return reader->slot(coroutine_id);
    }
};

// Helper to create awaitable
read_awaitable get_next_read(SharedFastqReader& reader, coroutine_handle<>& handle_storage, int coroutine_id) {
    return read_awaitable{&reader, &handle_storage, coroutine_id};
}

// Coroutine return type for query_pml_coroutine
struct MoveStructure::coroutine_task {
    struct promise_type {
        coroutine_task get_return_object() {
            return coroutine_task{coroutine_handle<promise_type>::from_promise(*this)};
        }
        suspend_always initial_suspend() { return {}; }
        suspend_always final_suspend() noexcept { return {}; }
        void unhandled_exception() {
            DEBUG_MSG_CO("DEBUG: Coroutine exception caught in unhandled_exception!\n");
            try {
                throw;  // Re-throw to see what it is
            } catch (const std::exception& e) {
                DEBUG_MSG_CO("DEBUG: Exception type: std::exception, what(): %s\n", e.what());
            } catch (...) {
                DEBUG_MSG_CO("DEBUG: Exception type: unknown\n");
            }
            // For now, still swallow it to prevent crash, but we'll know it happened
        }
        void return_void() {}
        suspend_always yield_value(monostate) noexcept { return {}; }
    };
    
    coroutine_handle<promise_type> coro;
    
    coroutine_task() : coro(nullptr) {}
    coroutine_task(coroutine_handle<promise_type> h) : coro(h) {}
    
    // Copy constructor - delete it since coroutine handles shouldn't be copied
    coroutine_task(const coroutine_task& other) = delete;
    
    // Move constructor
    coroutine_task(coroutine_task&& other) noexcept 
        : coro(other.coro) {
        other.coro = nullptr;
    }
    
    // Copy assignment - delete it
    coroutine_task& operator=(const coroutine_task& other) = delete;
    
    // Move assignment
    coroutine_task& operator=(coroutine_task&& other) noexcept {
        if (this != &other) {
            if (coro) {
                if (!coro.done()) {
                    coro.destroy();
                }
            }
            coro = other.coro;
            other.coro = nullptr;
        }
        return *this;
    }
    
    ~coroutine_task() { 
        if (coro) {
            if (!coro.done()) {
                coro.destroy();
            }
        }
    }
    
    // Helper to check if coroutine is valid and done
    bool is_done() const {
        return !coro || coro.done();
    }
};

// Coroutine implementation that requests reads from shared reader
MoveStructure::coroutine_task MoveStructure::query_pml_coroutine(
    SharedFastqReader& reader, 
    coroutine_handle<>& my_handle_storage,
    int coroutine_id) {
    
    int read_count = 0;
    int64_t total_bases = 0;
    const bool use_separator = this->use_separator();
    const int sep_adjust = use_separator ? 1 : 0;
    
    // Per-coroutine work buffers reused across reads (cleared, not reallocated,
    // each read) so there is no per-read heap churn.
    std::vector<uint16_t> matching_lens;
    std::string out_line;
    char numbuf[24];

    while (true) {
        DEBUG_MSG_CO("DEBUG: Coroutine %d about to await next read\n", coroutine_id);
        // Request next read - this suspends until scheduler reads it
        auto& read_data = co_await get_next_read(reader, my_handle_storage, coroutine_id);
        
        DEBUG_MSG_CO("DEBUG: Coroutine %d resumed with read, valid=%d\n", coroutine_id, read_data.valid ? 1 : 0);
        if (!read_data.valid) break;  // EOF
        
        read_count++;
        int read_bases = read_data.sequence.length();
        total_bases += read_bases;
        
        // Process the read in place: reverse this coroutine's own copy of the
        // sequence and use it directly as the query (no extra copies, no
        // per-read MoveQuery allocation).
        std::string& R = read_data.sequence;
        std::reverse(R.begin(), R.end());
        matching_lens.clear();
        int32_t roff = R.length() - 1;    // offset in read
        uint64_t idx = r - 1;             // row index
        uint64_t match_len = 0;           // match length (consecurive case 1s) so far
        uint64_t ff_count_tot = 0, scan_count = 0;
        uint64_t offset = get_n(idx) - 1;
        assert(idx < rlbwt.size());
        assert(offset < get_n(idx));

        // NOTE: the PML inner loop below is intentionally an INLINE reimplementation
        // of the repositioning and the LF/fast-forward step, rather than calls to the
        // shared reposition_thresholds / get_id / LF_move primitives (the way the MEM
        // and k-mer coroutines reuse them). A reuse version was written and measured
        // byte-identical but materially slower, and that was confirmed on a dedicated
        // c5.metal box (7 reps, tight): at W=8 the hand-inlined loop is ~109 ns/base
        // vs ~130 for the primitive-reusing version (~19%). PML's per-step work is so
        // small that once latency hiding is engaged, the un-inlined/un-fused primitive
        // calls dominate. LTO (which inlined LF_move) and even __attribute__((flatten))
        // on this function (which force-inlines the whole reposition chain) recovered
        // only ~5 ns of the ~21 ns gap (~126 ns/base) -- the rest is structural, not
        // just un-inlined calls. So the inline form is kept deliberately; it was
        // corrected to match the engine exactly (the four bugs noted in the foundation
        // commit). If you edit it, re-validate byte-for-byte against production
        // query_pml (--pml --reverse).
        while (roff > -1) {
            DEBUG_MSG_CO("DEBUG: Coroutine %d processing, roff=%d\n", coroutine_id, roff);
            char row_c_id = static_cast<char>((rlbwt[idx].n & (~mask_c)) >> SHIFT_C);
            char row_c = alphabet[row_c_id];
            if (!check_alphabet(R[roff])) {  // char doens't exist in reference
                match_len = 0;
            } else if (row_c == R[roff]) {   // Case 1: Match
                match_len++;
            } else {                         // Case 2: Reposition
                const char r_ch = R[roff];
                uint64_t r_ch_id = alphamap[static_cast<uint64_t>(r_ch)];
                uint64_t saved_idx = idx;
                bool up = false;
                if (idx == end_bwt_idx) {
                    up = offset < end_bwt_idx_thresholds[r_ch_id];
                    if (up) { // reposition up
                        if (saved_idx == 0) {
                            idx = r;
                        } else {
                            char row_c_up_id = static_cast<char>((rlbwt[saved_idx].n & (~mask_c)) >> SHIFT_C);
                            uint64_t repo_idx = saved_idx;
                            while (repo_idx > 0 && row_c_up_id != r_ch_id) {
                                repo_idx--;
                                row_c_up_id = static_cast<char>((rlbwt[repo_idx].n & (~mask_c)) >> SHIFT_C);
                            }
                            idx = (row_c_up_id == r_ch_id) ? repo_idx : r;
                        }
                        assert(idx < saved_idx);
                    } else { // reposition down
                        if (saved_idx == r - 1) {
                            idx = r;
                        } else {
                            char row_c_dn_id = static_cast<char>((rlbwt[saved_idx].n & (~mask_c)) >> SHIFT_C);
                            uint64_t repo_idx = saved_idx;
                            while (repo_idx < r - 1 && row_c_dn_id != r_ch_id) {
                                repo_idx++;
                                row_c_dn_id = static_cast<char>((rlbwt[repo_idx].n & (~mask_c)) >> SHIFT_C);
                            }
                            idx = (row_c_dn_id == r_ch_id) ? repo_idx : r;
                        }
                        assert(idx > saved_idx);
                    }
                } else {
                    // Handle both separator and non-separator cases
                    assert(!use_separator || r_ch_id != 0);
                    uint64_t r_ch_id_adj = r_ch_id - sep_adjust;
                    uint64_t threshold_value;
                    if (use_separator && row_c == SEPARATOR) {
                        // Handle separator case - use separators_thresholds
                        threshold_value = separators_thresholds[separators_thresholds_map[idx]].values[r_ch_id_adj];
                    } else {
                        // Regular case - use alphamap_3 (with or without separator adjustment)
                        r_ch_id_adj = alphamap_3[row_c_id - sep_adjust][r_ch_id_adj];
                        assert(r_ch_id_adj != 3);
                        threshold_value = get_thresholds(idx, r_ch_id_adj);
                    }
                    
                    if (offset < threshold_value) {
                        up = true;
                        char row_c_up_id = static_cast<char>((rlbwt[saved_idx].n & (~mask_c)) >> SHIFT_C);
                        uint64_t repo_idx = saved_idx;
                        while (repo_idx > 0 && row_c_up_id != r_ch_id) {
                            repo_idx--;
                            row_c_up_id = static_cast<char>((rlbwt[repo_idx].n & (~mask_c)) >> SHIFT_C);
                        }
                        idx = (row_c_up_id == r_ch_id) ? repo_idx : r;
                        assert(idx < saved_idx);
                    } else {
                        char row_c_dn_id = static_cast<char>((rlbwt[saved_idx].n & (~mask_c)) >> SHIFT_C);
                        uint64_t repo_idx = saved_idx;
                        while (repo_idx < r - 1 && row_c_dn_id != r_ch_id) {
                            repo_idx++;
                            row_c_dn_id = static_cast<char>((rlbwt[repo_idx].n & (~mask_c)) >> SHIFT_C);
                        }
                        idx = (row_c_dn_id == r_ch_id) ? repo_idx : r;
                        assert(idx > saved_idx);
                    }
                }
                match_len = 0;
                assert(idx < r);
                assert(idx < rlbwt.size());
                assert(alphabet[rlbwt[idx].get_c()] == R[roff] && "Repositioning failed - character mismatch");
                offset = up ? (get_n(idx) - 1) : 0;
            }
            // At this point, if repositioning was needed it has been done
            // 'idx' and 'offset' are the new move-structure row index and offset
            assert(offset < get_n(idx));
            assert(idx < rlbwt.size());
            
            // Record the PML (clamped to uint16_t), reusing the vector capacity.
            matching_lens.push_back(match_len > numeric_limits<uint16_t>::max()
                                    ? numeric_limits<uint16_t>::max()
                                    : static_cast<uint16_t>(match_len));
            roff--; // move left by 1 on query
            
            // Begin LF step.  Start by computing index of next row id where this
            // row maps
            auto& row = rlbwt[idx];
            uint64_t new_idx;
            if (row.offset >= (1U << 12)) {
                uint64_t res = static_cast<uint64_t>((row.offset & ~0x0FFF) >> 12) << 32;
                uint64_t rowid = static_cast<uint64_t>(row.id);
                rowid = rowid | res;
                new_idx = rowid;
            } else {
                new_idx = static_cast<uint64_t>(row.id);
            }
            my_prefetch_r((void*)(&(rlbwt[0]) + new_idx));
            offset = get_offset(idx) + offset;
            DEBUG_MSG_CO("DEBUG: Coroutine %d about to co_yield (prefetch)\n", coroutine_id);
            co_yield monostate{}; // wait for prefetch
            DEBUG_MSG_CO("DEBUG: Coroutine %d resumed after co_yield\n", coroutine_id);
            uint64_t n = get_n(new_idx);
            if (new_idx < r - 1 && offset >= n) {
                uint64_t niter = 0;
                while (new_idx < r - 1 && offset >= n) {
                    offset -= n;
                    niter++;
                    new_idx++;
                    n = get_n(new_idx);
                }
                assert(niter < numeric_limits<uint16_t>::max());
                ff_count_tot += static_cast<uint16_t>(niter);
            }
            idx = new_idx;
        }
        
        // Output results for this read into the reused line buffer (cleared
        // each read so capacity is retained). Building the whole line first
        // means arbitrarily long reads can never overflow/truncate the shared
        // buffer.
        out_line.clear();
        if (g_bpf_output) {
            // Binary BPF record, byte-compatible with mainline Movi's
            // output_base_stats: uint16 name length, name bytes, uint64 count,
            // then the raw uint16 matching-length array.
            uint16_t name_len = static_cast<uint16_t>(read_data.name.size());
            out_line.append(reinterpret_cast<const char*>(&name_len), sizeof(name_len));
            out_line.append(read_data.name);
            uint64_t cnt = matching_lens.size();
            out_line.append(reinterpret_cast<const char*>(&cnt), sizeof(cnt));
            if (cnt) out_line.append(reinterpret_cast<const char*>(matching_lens.data()), cnt * sizeof(uint16_t));
        } else {
            out_line.append(read_data.name);
            out_line.push_back(' ');
            for (size_t i = 0; i < matching_lens.size(); i++) {
                size_t w = write_uint_to_buffer(numbuf, sizeof(numbuf), matching_lens[i]);
                out_line.append(numbuf, w);
                if (i + 1 < matching_lens.size()) out_line.push_back(' ');
            }
            out_line.push_back('\n');
        }

        // Hand the completed record to the in-order emitter: it writes now if
        // this is the next read in input order, otherwise buffers a copy until
        // the earlier reads finish. (out_line is reused across reads.)
        g_emitter.emit(read_data.seq, out_line);
        DEBUG_MSG_CO("DEBUG: Coroutine %d finished processing read %d, looping back\n", coroutine_id, read_count);
    }
    
    co_return;
}


// Coroutine MEM query with latency hiding. Mirrors MoveStructure::query_mems /
// query_mem_bml (min_mem_length > 1 path). Yields only at the LF destination
// jump inside each backward/forward search step, and before ftab lookups.
// extend_bidirectional is inlined so the single co_yield can be placed at the
// fw LF_move; the subsequent sequential work (skip scan, rc walk) gets no yield
// because the hardware prefetcher covers sequential rlbwt row access.
// The skip scan also accumulates fw_count (= count of c in fw_before), replacing
// the separate fw_interval.count(rlbwt) scan for the rc-end walk.
MoveStructure::coroutine_task MoveStructure::query_mem_coroutine(
    SharedFastqReader& reader,
    coroutine_handle<>& my_handle_storage,
    int coroutine_id) {

    const int32_t min_mem_length = movi_options->get_min_mem_length();
    const size_t  ftab_k = movi_options->get_ftab_k();
    char numbuf[24];
    std::string out_line;

    // Split backward_search_step at the LF destination jump — the only cold-miss
    // access in search. step_prep does update_interval (sequential, hw-prefetchable)
    // then computes and prefetches the LF destination rows. step_finish does the
    // LF_move using the precomputed ids so the prefetched rows are exactly used.
    uint64_t pf_id_s = 0, pf_id_e = 0;
    auto step_prep = [&](MoveInterval& iv, char c) -> bool {
        if (!check_alphabet(c)) { iv.make_empty(); return false; }
        update_interval(iv, c);
        if (iv.is_empty()) return false;
        pf_id_s = get_id(iv.run_start);
        pf_id_e = get_id(iv.run_end);
        my_prefetch_r((void*)(&(rlbwt[0]) + pf_id_s));
        my_prefetch_r((void*)(&(rlbwt[0]) + pf_id_e));
        return true;
    };
    auto step_finish = [&](MoveInterval& iv) {
        LF_move(iv.offset_start, iv.run_start, pf_id_s);
        LF_move(iv.offset_end, iv.run_end, pf_id_e);
    };

    // Prefetch both fw and rc ftab entries before initialize_bidirectional_search.
    auto ftab_pf_bidi = [&](std::string& qs, int32_t pos) -> bool {
        if (ftab_k <= 1 || pos < static_cast<int32_t>(ftab_k) - 1) return false;
        int32_t start = pos - static_cast<int32_t>(ftab_k) + 1;
        uint64_t fw = kmer_to_number(ftab_k, qs, start, alphamap, false);
        uint64_t rc = kmer_to_number(ftab_k, qs, start, alphamap, true);
        bool any = false;
        if (fw != std::numeric_limits<uint64_t>::max()) { my_prefetch_r((void*)(&(ftab[0]) + fw)); any = true; }
        if (rc != std::numeric_limits<uint64_t>::max()) { my_prefetch_r((void*)(&(ftab[0]) + rc)); any = true; }
        return any;
    };
    // Prefetch fw ftab entry before initialize_backward_search (next MEM start).
    auto ftab_pf = [&](std::string& qs, int32_t pos) -> bool {
        if (ftab_k <= 1 || pos < static_cast<int32_t>(ftab_k) - 1) return false;
        uint64_t code = kmer_to_number(ftab_k, qs, pos - static_cast<int32_t>(ftab_k) + 1, alphamap);
        if (code == std::numeric_limits<uint64_t>::max()) return false;
        my_prefetch_r((void*)(&(ftab[0]) + code));
        return true;
    };

    while (true) {
        auto& read_data = co_await get_next_read(reader, my_handle_storage, coroutine_id);
        if (!read_data.valid) break;

        MoveQuery mq(std::move(read_data.sequence));
        mq.set_query_id(read_data.name);
        std::string& query_seq = mq.query();
        const int32_t qlen = static_cast<int32_t>(query_seq.length());

        int32_t pos_on_r = 0;
        // Cast to size_t replicates production's `pos_on_r < query_seq.length()`
        // (int32_t < size_t): when pos_on_r goes negative after a MEM that reaches
        // end-of-read, the cast wraps to SIZE_MAX > qlen and the loop exits.
        while (static_cast<size_t>(pos_on_r) < static_cast<size_t>(qlen)) {
            // ---- query_mem_bml(pos_on_r) inlined with prefetch + co_yield ----
            if (pos_on_r + min_mem_length > qlen) { pos_on_r = qlen; break; }

            uint64_t match_len = 0;
            int32_t init_pos = pos_on_r + min_mem_length - 1;

            // Prefetch both ftab entries before the cold bidirectional initialization.
            if (ftab_pf_bidi(query_seq, init_pos)) co_yield monostate{};
            MoveBiInterval bi_interval = initialize_bidirectional_search(mq, init_pos, match_len);
            bool ftab_skip = match_len <= 1 && ftab_k <= static_cast<size_t>(min_mem_length);
            --init_pos;

            // Prefetch fw_interval entry-point rows so the first update_interval
            // in the left-extension loop doesn't stall on a cold ftab-seeded interval.
            if (!ftab_skip && match_len > 1 && !bi_interval.fw_interval.is_empty()) {
                my_prefetch_r((void*)(&(rlbwt[0]) + bi_interval.fw_interval.run_start));
                my_prefetch_r((void*)(&(rlbwt[0]) + bi_interval.fw_interval.run_end));
                co_yield monostate{};
            }

            bool failed = false;
            if (ftab_skip) {
                // ftab can't seed bidirectional; pure backward search to find pos_on_r.
                MoveInterval fw_interval = bi_interval.fw_interval;
                for (size_t j = 0; j <= static_cast<size_t>(init_pos - pos_on_r); ++j) {
                    if (!step_prep(fw_interval, query_seq[init_pos - j])) {
                        pos_on_r = init_pos - j + 1; failed = true; break;
                    }
                    co_yield monostate{};
                    step_finish(fw_interval);
                    ++match_len;
                }
                if (!failed) throw std::runtime_error("Extended past failed ftab");
            } else {
                // Left extension: extend_bidirectional inlined so co_yield is placed
                // only at the LF jump on fw_interval. Sequential scan and rc walk
                // follow step_finish with no yield.
                for (size_t j = 0; j <= static_cast<size_t>(init_pos - pos_on_r); ++j) {
                    char c_ = query_seq[init_pos - j];
                    char c_comp = complement(c_);

                    MoveInterval fw_before = bi_interval.fw_interval;

                    if (!step_prep(bi_interval.fw_interval, c_)) {
                        pos_on_r = init_pos - j + 1; failed = true; break;
                    }
                    co_yield monostate{};   // hide LF destination miss on fw
                    step_finish(bi_interval.fw_interval);

                    // Combined skip-count + fw_count scan over fw_before (one pass).
                    // skip  → how far to walk rc_interval.start forward.
                    // fw_count → count of c_ in fw_before = count of new fw_interval
                    //            (bidirectional BWT invariant), used for rc end walk
                    //            instead of a separate fw_interval.count(rlbwt) call.
                    uint64_t skip = 0, fw_count = 0;
                    uint64_t cur = fw_before.run_start;
                    uint64_t cur_off = fw_before.offset_start;
                    while (cur <= fw_before.run_end) {
                        if (cur != end_bwt_idx) {
                            char ch = get_char(cur);
                            uint64_t n = (cur != fw_before.run_end) ?
                                get_n(cur) - cur_off :
                                fw_before.offset_end - cur_off + 1;
                            if (complement(ch) < c_comp) skip += n;
                            if (ch == c_) fw_count += n;
                        } else {
                            skip += 1;
                        }
                        ++cur;
                        cur_off = 0;
                    }

                    // Walk rc_interval.start forward by skip.
                    MoveInterval& rc = bi_interval.rc_interval;
                    while (skip != 0) {
                        uint64_t rows_after = get_n(rc.run_start) - 1 - rc.offset_start;
                        if (rows_after >= skip) {
                            rc.offset_start += skip;
                            skip = 0;
                        } else {
                            rc.run_start += 1;
                            rc.offset_start = 0;
                            skip -= rows_after + 1;
                        }
                    }

                    // Walk rc_interval.end to rc_interval.start + (fw_count - 1).
                    rc.run_end = rc.run_start;
                    rc.offset_end = rc.offset_start;
                    uint64_t fw_skip = fw_count - 1;
                    while (fw_skip != 0) {
                        uint64_t rows_after = get_n(rc.run_end) - 1 - rc.offset_end;
                        if (rows_after >= fw_skip) {
                            rc.offset_end += fw_skip;
                            fw_skip = 0;
                        } else {
                            rc.run_end += 1;
                            rc.offset_end = 0;
                            fw_skip -= rows_after + 1;
                        }
                    }

                    bi_interval.match_len += 1;
                    ++match_len;
                }
            }
            if (failed) continue;

            // Forward extension: find right end of MEM (exclusive).
            // forward_search_step(c, iv) == backward_search_step(complement(c), iv).
            MoveInterval rc_interval = bi_interval.rc_interval;
            MoveInterval rc_before = rc_interval;
            size_t i;
            for (i = static_cast<size_t>(pos_on_r) + min_mem_length;
                 i < static_cast<size_t>(qlen); ++i) {
                rc_before = rc_interval;
                if (!step_prep(rc_interval, complement(query_seq[i]))) {
                    rc_interval = rc_before; break;
                }
                co_yield monostate{};       // hide LF destination miss on rc
                step_finish(rc_interval);
                ++match_len;
            }

            mq.add_mem(pos_on_r, i, rc_interval.count(rlbwt));

            // Backward extension to find next candidate MEM's left end.
            size_t end_pos_on_r = i;
            if (end_pos_on_r < static_cast<size_t>(qlen)) {
                init_pos = static_cast<int32_t>(end_pos_on_r);
                match_len = 0;
                if (ftab_pf(query_seq, init_pos)) co_yield monostate{};
                MoveInterval fw_interval = initialize_backward_search(mq, init_pos, match_len);
                ++match_len;
                --init_pos;
                for (i = 0; i <= static_cast<size_t>(init_pos - (pos_on_r + 1)); ++i) {
                    if (!step_prep(fw_interval, query_seq[init_pos - i])) break;
                    co_yield monostate{};   // hide LF destination miss on fw
                    step_finish(fw_interval);
                    ++match_len;
                }
            }
            // When end_pos_on_r >= qlen the else branch is skipped (matching
            // the reference), leaving i == qlen from the forward loop so that
            // pos_on_r = init_pos - qlen + 1 goes negative; the size_t cast in
            // the outer while condition then terminates the loop correctly.
            pos_on_r = (init_pos - static_cast<int32_t>(i)) + 1;
        }

        // Emit this read's MEMs (id \t start \t end \t count), one per line.
        out_line.clear();
        for (auto& m : mq.get_mems()) {
            out_line.append(mq.get_query_id());
            out_line.push_back('\t');
            out_line.append(numbuf, write_uint_to_buffer(numbuf, sizeof(numbuf), m.start));
            out_line.push_back('\t');
            out_line.append(numbuf, write_uint_to_buffer(numbuf, sizeof(numbuf), m.end));
            out_line.push_back('\t');
            out_line.append(numbuf, write_uint_to_buffer(numbuf, sizeof(numbuf), m.count));
            out_line.push_back('\n');
        }
        g_emitter.emit(read_data.seq, out_line);
    }
    co_return;
}


// Coroutine k-mer query with latency hiding. Mirrors MoveStructure::query_all_kmers
// (together with query_kmers_from and look_ahead_backward_search) exactly -- both
// skip heuristics included -- so output is byte-identical to production
// `movi query --kmer [--kmer-count]`. The only change versus production is a
// prefetch + co_yield issued before each backward_search_step (the cache-missing
// LF_move). Those two inner backward-search loops are inlined here because suspend
// points cannot live in called functions. One coroutine serves both modes via
// count_mode. k-mer search uses plain (not bidirectional) backward search, so no
// reverse-complement index is needed and the read is NOT reversed.
MoveStructure::coroutine_task MoveStructure::query_kmer_coroutine(
    SharedFastqReader& reader,
    coroutine_handle<>& my_handle_storage,
    int coroutine_id) {

    const size_t  ftab_k = movi_options->get_ftab_k();
    const int32_t k = static_cast<int32_t>(movi_options->get_k());
    const bool    count_mode = movi_options->is_kmer_count();
    std::string out_line;
    char numbuf[24];

    // One backward_search_step, split so the prefetch + co_yield sit at the only
    // long-range jump in MODE-6 interval search: the LF destination row. step_prep
    // does the sequential part (update_interval's run scan needs no prefetch),
    // then -- the instant we learn the LF destinations get_id(run_start/run_end) --
    // prefetches those move-structure rows and returns. The caller co_yields so the
    // prefetch lands while other coroutines run, then step_finish does the LF_move
    // that actually dereferences the destination (get_n on the prefetched row).
    // The precomputed ids are handed to LF_move so the prefetched row is exactly
    // the one accessed. Returns false (no jump) when the step empties the interval.
    uint64_t pf_id_s = 0, pf_id_e = 0;
    auto step_prep = [&](MoveInterval& iv, char c) -> bool {
        if (!check_alphabet(c)) { iv.make_empty(); return false; }
        update_interval(iv, c);
        if (iv.is_empty()) return false;
        pf_id_s = get_id(iv.run_start);
        pf_id_e = get_id(iv.run_end);
        my_prefetch_r((void*)(&(rlbwt[0]) + pf_id_s));
        my_prefetch_r((void*)(&(rlbwt[0]) + pf_id_e));
        return true;
    };
    auto step_finish = [&](MoveInterval& iv) {
        LF_move(iv.offset_start, iv.run_start, pf_id_s);
        LF_move(iv.offset_end, iv.run_end, pf_id_e);
    };

    // The other long-range jump: initialize_backward_search's ftab lookup. The
    // ftab index (kmer_code) is computed cheaply from the query characters
    // (sequential, cached), but ftab[kmer_code] is a random access into a large
    // array (tens to hundreds of MB). Prefetch that entry; the caller co_yields so
    // it lands before initialize_backward_search dereferences it. Returns false
    // when no ftab access will happen (ftab_k<=1, too few chars, or illegal kmer)
    // so the caller can skip the (then-pointless) yield -- no overhead in no-ftab.
    auto ftab_pf = [&](std::string& qs, int32_t pos) -> bool {
        if (ftab_k <= 1 || pos < static_cast<int32_t>(ftab_k) - 1) return false;
        uint64_t code = kmer_to_number(ftab_k, qs, pos - static_cast<int32_t>(ftab_k) + 1, alphamap);
        if (code == std::numeric_limits<uint64_t>::max()) return false;
        my_prefetch_r((void*)(&(ftab[0]) + code));
        return true;
    };

    while (true) {
        auto& read_data = co_await get_next_read(reader, my_handle_storage, coroutine_id);
        if (!read_data.valid) break;

        // Own the sequence (cheap move) so the shared primitives, which take a
        // MoveQuery&, can be reused directly.
        MoveQuery mq(std::move(read_data.sequence));
        mq.set_query_id(read_data.name);
        std::string& query_seq = mq.query();
        const int32_t qlen = static_cast<int32_t>(query_seq.length());
        // Match production's denominator exactly (movi.cpp output_kmers call),
        // including the unsigned underflow when the read is shorter than k.
        const uint64_t all_kmer_count = mq.length() - movi_options->get_k() + 1;

        int32_t pos_on_r = qlen - 1;

        // k == 1 is a no-op for per-read kmer output in production query_all_kmers
        // (it only bumps a global stat, never calls add_kmer).
        if (k != 1 && qlen >= k) {
            while (pos_on_r >= 0 && !check_alphabet(query_seq[pos_on_r])) pos_on_r -= 1;

            int32_t step = k / 3;
            if (k - step < static_cast<int32_t>(ftab_k)) step = k - static_cast<int32_t>(ftab_k) - 1;

            while (pos_on_r >= k - 1) {
                bool did_search = true;

                if (pos_on_r >= k - 1 + step) {
                    // ---- look_ahead_backward_search(mq, pos_on_r, step) inlined ----
                    uint64_t la_match_len = 0;
                    int32_t la_pos = pos_on_r - step;
                    if (ftab_pf(query_seq, la_pos)) co_yield monostate{};
                    MoveInterval la_iv = initialize_backward_search(mq, la_pos, la_match_len);
                    const int32_t la_max = k - step - static_cast<int32_t>(la_match_len);
                    // inlined backward_search(query_seq, la_pos, la_iv, la_max)
                    const int32_t la_saved = la_pos;
                    while (la_pos > 0 && !la_iv.is_empty()) {
                        if (step_prep(la_iv, query_seq[la_pos - 1])) {
                            co_yield monostate{};
                            step_finish(la_iv);
                            la_pos -= 1;
                        }
                        if (la_saved - la_pos > la_max) break;
                    }
                    if (pos_on_r - la_pos < k - 1) {     // look-ahead says "skip"
                        pos_on_r = pos_on_r - step - 1;
                        did_search = false;
                    }
                }

                if (did_search) {
                    // ---- query_kmers_from(mq, pos_on_r, single=count_mode, &kc) inlined ----
                    const int32_t kmer_end = pos_on_r;   // captured before pos_on_r moves
                    int32_t pos_saved = pos_on_r;
                    uint64_t match_len = 0;
                    MoveInterval init_iv;
                    do {
                        if (ftab_pf(query_seq, pos_on_r)) co_yield monostate{};
                        init_iv = initialize_backward_search(mq, pos_on_r, match_len);
                        if (match_len == 0 && ftab_k > 1) {
                            pos_on_r -= 1;
                            pos_saved = pos_on_r;
                        }
                    } while (match_len == 0 && pos_on_r >= k - 1 && ftab_k > 1);

                    if (count_mode && movi_options->is_kmer_bv()) {
                        // ---- keep-going count via the bitvector (W1+W2; inlined
                        // query_kmers_count_bv). Walk left at presence speed; at every
                        // position whose covered match is >= k, resolve the leftmost
                        // k-mer's count from the count bitvector. The shrinking interval
                        // is always a subset of that k-mer's group, so pred(lb)/succ(rb)
                        // give its exact count. The LF step is the same prefetch+co_yield
                        // point as presence; the pred/succ probe rides the same pipeline.
                        const int32_t anchor = pos_saved;
                        if (init_iv.is_empty()) {
                            pos_on_r = pos_saved - 1;
                        } else {
                            MoveInterval iv = init_iv;
                            uint64_t kmers_found = 0;
                            while (true) {
                                if (anchor - pos_on_r + 1 >= k) {
                                    uint64_t cnt = kmer_count_from_bv(iv);
                                    mq.add_kmer(pos_on_r, /*present=*/1,
                                                std::numeric_limits<uint64_t>::max(), cnt);
                                    kmers_found += 1;
                                }
                                if (pos_on_r <= 0) break;
                                if (!step_prep(iv, query_seq[pos_on_r - 1])) break;
                                co_yield monostate{};
                                step_finish(iv);
                                pos_on_r -= 1;
                            }
                            pos_on_r = (kmers_found > 0) ? (pos_on_r + k - 2) : (pos_saved - 1);
                        }
                    } else {
                    const int32_t bs_max = count_mode
                        ? (k - static_cast<int32_t>(match_len) - 2)
                        : std::numeric_limits<int32_t>::max();
                    // inlined backward_search(query_seq, pos_on_r, init_iv, bs_max)
                    MoveInterval bs_iv = init_iv;
                    MoveInterval bs_prev = bs_iv;
                    const int32_t bs_saved = pos_on_r;
                    while (pos_on_r > 0 && !bs_iv.is_empty()) {
                        bs_prev = bs_iv;
                        if (step_prep(bs_iv, query_seq[pos_on_r - 1])) {
                            co_yield monostate{};
                            step_finish(bs_iv);
                            pos_on_r -= 1;
                        }
                        if (bs_saved - pos_on_r > bs_max) break;
                    }
                    MoveInterval bs_result = bs_iv.is_empty() ? bs_prev : bs_iv;

                    uint64_t found = 0, kc = 0;
                    if (bs_result.is_empty()) {
                        pos_on_r = pos_saved - 1;
                    } else if (pos_saved - pos_on_r >= k - 1) {
                        found = static_cast<uint64_t>(pos_saved - pos_on_r - k + 2);
                        if (count_mode) kc = bs_result.count(rlbwt);
                        pos_on_r = pos_on_r + k - 2;
                    } else {
                        pos_on_r = pos_saved - 1;
                    }

                    if (count_mode) {
                        // present=1 per k-mer keeps the found/all header correct;
                        // kc is the displayed multiplicity.
                        mq.add_kmer(kmer_end - k + 1, /*present=*/1,
                                    std::numeric_limits<uint64_t>::max(), kc);
                    } else {
                        mq.add_kmer(pos_on_r + 2 - k, found);
                    }
                    }
                }

                while (pos_on_r >= 0 && !check_alphabet(query_seq[pos_on_r])) pos_on_r -= 1;
            }
        }

        // Build the output line exactly like output_kmers(): id \t found/all \t poslist
        out_line.clear();
        out_line.append(mq.get_query_id());
        out_line.push_back('\t');
        out_line.append(numbuf, write_uint_to_buffer(numbuf, sizeof(numbuf), mq.found_kmer_count));
        out_line.push_back('/');
        out_line.append(numbuf, write_uint_to_buffer(numbuf, sizeof(numbuf), all_kmer_count));
        out_line.push_back('\t');
        out_line.append(mq.get_matching_lengths_string());
        out_line.push_back('\n');
        g_emitter.emit(read_data.seq, out_line);
    }
    co_return;
}


void process_fastq(const string& fastq_file, const string& index_dir, int concurrency) {
    MoviOptions movi_options;
    movi_options.set_index_dir(index_dir);
    if (g_mem_query) {
        movi_options.set_mem();
        movi_options.set_ftab_k(static_cast<uint32_t>(g_ftab_k));
        movi_options.set_min_mem_length(static_cast<uint32_t>(g_min_mem_length));
        movi_options.set_multi_ftab(false);
    } else if (g_kmer_query) {
        movi_options.set_kmer();
        movi_options.set_kmer_count(g_kmer_count);
        movi_options.set_kmer_bv(g_kmer_bv);
        movi_options.set_k(static_cast<uint32_t>(g_k));
        movi_options.set_ftab_k(static_cast<uint32_t>(g_ftab_k));
        movi_options.set_multi_ftab(false);
    } else {
        movi_options.set_pml(); // Enable PML mode
    }

    MoveStructure mv(&movi_options);
    mv.deserialize();
    if ((g_mem_query && g_ftab_k > 0) || (g_kmer_query && g_ftab_k > 1)) mv.read_ftab();
    // Run-local k-mer count structure (and the MPHF id bitvector): only needed for
    // the --kmer-bv path; loaded once after the move structure.
    if (g_kmer_query && g_kmer_bv) mv.load_kmerbv(static_cast<uint32_t>(g_k));
    fprintf(stderr_buf, "Successfully loaded Movi index from: %s\n", index_dir.c_str());

    // PML/text/BPF header only applies to the PML output; MEM and k-mer emit
    // their own tab-delimited lines with no header.
    if (!g_mem_query && !g_kmer_query) {
        if (g_bpf_output) {
            BPFHeader header;
            header.init(16);  // matching lengths are stored as 16-bit values
            fwrite(&header, sizeof(header), 1, stdout_buf);
        } else {
            fprintf(stdout_buf, "# Read_ID PML_Values\n");
        }
    }
    g_emitter.ordered = g_ordered_output;
    
    auto start_time = steady_clock::now();
    int64_t total_bases = 0;
    
    // Single shared reader for all coroutines
    SharedFastqReader reader(fastq_file, concurrency);
    
    // Storage for coroutine handles waiting for reads
    std::vector<coroutine_handle<>> waiting_handles(concurrency);
    std::vector<MoveStructure::coroutine_task> coroutines;
    coroutines.reserve(concurrency);
    
    // Initialize all coroutines - use move semantics
    for (int i = 0; i < concurrency; ++i) {
        if (g_mem_query)
            coroutines.push_back(std::move(mv.query_mem_coroutine(reader, waiting_handles[i], i)));
        else if (g_kmer_query)
            coroutines.push_back(std::move(mv.query_kmer_coroutine(reader, waiting_handles[i], i)));
        else
            coroutines.push_back(std::move(mv.query_pml_coroutine(reader, waiting_handles[i], i)));
        // They start suspended, so resume to begin execution.
        if (coroutines[i].coro) {
            coroutines[i].coro.resume();
        }
    }
    
    // Master scheduler loop: coordinate reads and coroutine execution
    int active_coroutines = concurrency;
    int iteration = 0;
    int total_reads_served = 0;
    
    while (active_coroutines > 0) {
        iteration++;
        active_coroutines = 0;
        
        for (int i = 0; i < concurrency; ++i) {
            if (coroutines[i].is_done()) {
                DEBUG_MSG_CO("DEBUG: Coroutine %d is done\n", i);
                continue;
            }
            
            // Check if this coroutine is waiting for a read
            if (waiting_handles[i]) {
                total_reads_served++;
                // Move the next batched record into this coroutine's slot.
                if (reader.serve_next(i)) {
                    total_bases += reader.slot(i).sequence.length();
                }
                // Clear the slot BEFORE resuming. If the coroutine finishes this
                // read without ever hitting co_yield (e.g. a read shorter than k,
                // which skips the whole search), it will loop straight back to
                // co_await and re-register its handle here during the resume();
                // clearing after the resume would clobber that fresh await and the
                // coroutine would be resumed next iteration with no new read served
                // -- re-processing the moved-from (empty) slot as a phantom record.
                auto h = waiting_handles[i];
                waiting_handles[i] = nullptr;
                h.resume();
                if (!coroutines[i].is_done()) {
                    active_coroutines++;
                }
            } else {
                DEBUG_MSG_CO("DEBUG: Scheduler iteration %d, coroutine %d not waiting (handle is %s)\n",
                     iteration, i, (waiting_handles[i] ? "set" : "null"));
                // Not waiting for read - just resume (handles co_yield for voluntary suspension)
                if (coroutines[i].coro) {
                    coroutines[i].coro.resume();
                    if (!coroutines[i].is_done()) {
                        active_coroutines++;
                    }
                }
            }
        }
    }
    
    fprintf(stderr_buf, "Scheduler completed: %d iterations, %d reads served to coroutines\n", 
            iteration, total_reads_served);
    
    // Emit any buffered output in input order and flush (part of the timed region).
    g_emitter.finish();

    // Calculate final statistics
    auto end_time = steady_clock::now();
    double total_elapsed = duration_cast<duration<double>>(end_time - start_time).count();
    
    // Calculate nanoseconds per base
    double nanoseconds_per_base = 0.0;
    if (total_bases > 0) {
        double total_elapsed_ns = total_elapsed * 1e9; // Convert seconds to nanoseconds
        nanoseconds_per_base = total_elapsed_ns / total_bases;
    }
    
    fprintf(stderr_buf, "Completed processing (elapsed: %.2f sec, %.1f ns/base, bases=%ld)\n",
            total_elapsed, nanoseconds_per_base, static_cast<long>(total_bases));
}

void print_usage(const char* program_name) {
    fprintf(stderr_buf, "Usage: %s [--debug] [--bpf] [--reorder] [query opts] <fastq_file> <index_dir> [concurrency]\n", program_name);
    fprintf(stderr_buf, "  --debug:     Enable debug output\n");
    fprintf(stderr_buf, "  --bpf:       Write binary BPF output (like mainline Movi) instead of text\n");
    fprintf(stderr_buf, "  --reorder:   Emit reads in input order (enable the reorder buffer; default is completion order)\n");
    fprintf(stderr_buf, "  --mem --ftab-k N --min-mem-length N:  MEM query mode\n");
    fprintf(stderr_buf, "  --kmer | --kmer-count, -k/--k-length N [--ftab-k N]:  k-mer query mode\n");
    fprintf(stderr_buf, "               (default query mode is PML)\n");
    fprintf(stderr_buf, "  fastq_file: Input FASTQ file with reads\n");
    fprintf(stderr_buf, "  index_dir:  Directory containing the Movi index\n");
    fprintf(stderr_buf, "  concurrency: Number of concurrent coroutines (default: 1)\n");
    fprintf(stderr_buf, "\n");
    fprintf(stderr_buf, "This program computes Pseudo Matching Lengths (PMLs) for reads\n");
    fprintf(stderr_buf, "using a Movi regular-thresholds index (MODE=6).\n");
}

int main(int argc, char* argv[]) {
    // Initialize buffered I/O
    init_buffered_io();
    
    int arg_idx = 1;

    // Parse leading flags (--debug, --bpf) in any order.
    // Accept both long (--foo) and short (-k) option flags; a positional (fastq
    // path or concurrency) never starts with '-'.
    while (arg_idx < argc && argv[arg_idx][0] == '-' && string(argv[arg_idx]) != "-") {
        string flag = argv[arg_idx];
        if (flag == "--debug") { debug_enabled = true; }
        else if (flag == "--bpf") { g_bpf_output = true; }
        else if (flag == "--reorder") { g_ordered_output = true; }
        else if (flag == "--mem") { g_mem_query = true; }
        else if (flag == "--kmer") { g_kmer_query = true; }
        else if (flag == "--kmer-count") { g_kmer_query = true; g_kmer_count = true; }
        else if (flag == "--kmer-bv") { g_kmer_bv = true; }
        else if (flag == "-k" || flag == "--k-length") { if (arg_idx + 1 >= argc) { print_usage(argv[0]); return 1; } g_k = std::stoi(argv[++arg_idx]); }
        else if (flag == "--ftab-k") { if (arg_idx + 1 >= argc) { print_usage(argv[0]); return 1; } g_ftab_k = std::stoi(argv[++arg_idx]); }
        else if (flag == "--min-mem-length") { if (arg_idx + 1 >= argc) { print_usage(argv[0]); return 1; } g_min_mem_length = std::stoi(argv[++arg_idx]); }
        else { fprintf(stderr_buf, "Unknown flag: %s\n", flag.c_str()); print_usage(argv[0]); return 1; }
        arg_idx++;
    }

    // Validate query-mode flag combinations.
    if (g_mem_query && g_kmer_query) {
        fprintf(stderr_buf, "Error: --mem and --kmer/--kmer-count are mutually exclusive\n");
        return 1;
    }
    if (g_kmer_query && g_k < 1) {
        fprintf(stderr_buf, "Error: --kmer/--kmer-count requires -k/--k-length N (N >= 1)\n");
        return 1;
    }

    // Check remaining arguments
    if (argc - arg_idx < 2 || argc - arg_idx > 3) {
        fprintf(stderr_buf, "argc = %d\n", argc);
        print_usage(argv[0]);
        return 1;
    }
    
    string fastq_file{argv[arg_idx]}, index_dir{argv[arg_idx + 1]};
    int concurrency = 1;
    
    if (argc - arg_idx == 3) {
        try {
            concurrency = std::stoi(argv[arg_idx + 2]);
            if (concurrency < 1) {
                fprintf(stderr_buf, "Error: Concurrency must be at least 1\n");
                return 1;
            }
        } catch (const exception& e) {
            fprintf(stderr_buf, "Error: Invalid concurrency value: %s\n", argv[arg_idx + 2]);
            return 1;
        }
    }
    
    try {
        process_fastq(fastq_file, index_dir, concurrency);
    } catch (const exception& e) {
        fprintf(stderr_buf, "Error: %s\n", e.what());
        return 1;
    }
    return 0;
}


