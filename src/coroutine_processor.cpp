// coroutine_processor.cpp: coroutine-based latency hiding -- query style "(c)" -- the engine.
//
// Movi runs each hot query in up to three latency-hiding styles (docs/latency_hiding.md
// is the full map); this file is style (c):
//   (a) sequential, no latency hiding: query_pml.cpp, query_zml.cpp, query_kmer.cpp,
//         query_kmer_bv.cpp, query_mem.cpp (shared primitives in move_structure_search.cpp).
//   (b) manual "strand" state machine: read_processor.cpp (process_latency_hiding),
//         struct Strand in read_processor.hpp.
//   (c) coroutines (this file): prefetch, then co_yield at each prefetch point; the
//         scheduler resumes another coroutine while the cache line fills. Reached from
//         the main movi binary via `movi query --coroutine` (dispatch in movi.cpp).
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
 * The C++20-coroutine latency-hiding query engine for Movi (PML, MEM and k-mer
 * queries), linked into the main movi binary and driven by `movi query --coroutine`.
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
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include <omp.h>
#include <mutex>
#include <cstring>
#include <map>
#include <sstream>
#include <iostream>

#include <sdsl/int_vector.hpp>
#include <zlib.h>

#include "utils.hpp"
#include "move_structure.hpp"
#include "move_query.hpp"
#include "movi_options.hpp"
#include "coroutine_processor.hpp"
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

// Coroutine latency hiding is implemented for the threshold index modes -- regular-
// thresholds (6), sampled-thresholds (7), and blocked-thresholds (8) -- in non-color
// builds. Those are the modes that carry the thresholds PML/MEM/k-mer queries need,
// and the engine uses the mode-portable move-structure primitives throughout. Under
// any other MODE/COLOR_MODE this translation unit compiles to nothing, so it is safe
// to keep in the shared MOVI_SOURCES; the `movi query --coroutine` dispatch in
// movi.cpp is guarded by the same condition.
#if (MODE == 6 || MODE == 7 || MODE == 8) && COLOR_MODE == 0

// The only cross-cutting global the query bodies need: the debug flag read by the
// DEBUG_MSG_CO macro. All per-query parameters (query mode, k, ftab_k, min-mem,
// count/bv flags) are read directly from the MoveStructure's MoviOptions in the
// bodies, and the scheduler reads them from the CoroutineQueryOptions argument, so
// no per-query file-scope globals are needed.
static bool debug_enabled = false;

// Setup buffered stdout and stderr
static FILE* stdout_buf = nullptr;
static FILE* stderr_buf = nullptr;

// Initialize buffered I/O
static void init_buffered_io() {
    if (stdout_buf != nullptr) return;  // idempotent: safe to call more than once
    stdout_buf = stdout;
    stderr_buf = stderr;
    // Use full buffering for stdout (better performance)
    setvbuf(stdout_buf, nullptr, _IOFBF, 8192);
    // Use line buffering for stderr (immediate error visibility)
    setvbuf(stderr_buf, nullptr, _IOLBF, 0);
}

// DEBUG macro that only prints when NDEBUG is not defined and debug flag is set
// Uses fprintf-style format strings
#ifndef NDEBUG
#define DEBUG_MSG_CO(...) do { if (debug_enabled) { fprintf(stderr_buf, __VA_ARGS__); } } while(0)
#else
#define DEBUG_MSG_CO(...) ((void)0)
#endif

// In-order, chunked output emitter. Coroutines finish reads out of order; this
// buffers completed lines and writes them in input (read) order, so downstream
// tools see the same order as the input. Single-threaded (the scheduler resumes
// one coroutine at a time), so no locking is needed.
struct OrderedEmitter {
    uint64_t next_emit = 0;
    bool ordered = false;         // when true, buffer + write in input order; else completion order
    bool threaded = false;        // when true, guard emit() with the mutex (many worker threads)
    std::mutex mu;                // serializes emit()/finish() across worker threads
    std::map<uint64_t, std::string> pending;  // completed lines awaiting their turn
    // Output destination, set per run by run_coroutine_query to the appropriate
    // output_files stream (or std::cout for --stdout) so coroutine output lands in
    // the same place as the sequential path.
    std::ostream* out = &std::cout;

    void write_line(const std::string& line) {
        out->write(line.data(), line.size());
    }
    // Emit the line for read 'seq'. Writes immediately if it is the next read in
    // input order, then drains any consecutive buffered lines; otherwise buffers
    // a copy until the earlier reads complete. When more than one worker thread is
    // running the body is taken under the mutex; a single-thread run skips the lock
    // entirely, so the sequential fast path pays no locking overhead.
    void emit(uint64_t seq, const std::string& line) {
        if (threaded) {
            std::lock_guard<std::mutex> lk(mu);
            emit_body(seq, line);
        } else {
            emit_body(seq, line);
        }
    }
    void emit_body(uint64_t seq, const std::string& line) {
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
        out->flush();
        next_emit = 0;  // reset for any subsequent run reusing this static emitter
    }
};
static OrderedEmitter g_emitter;

// One parsed FASTQ/A record. `seq` is the global input-order index, assigned by
// FastqSource under its lock so it is gap-free and monotonic across every worker
// thread -- which is what lets the OrderedEmitter reassemble input order.
struct CoReadData {
    bool valid = false;
    string name;
    string sequence;
    uint64_t seq = 0;
};

// The single shared input for a run, in one of two modes:
//
//  - MMAP (raw, uncompressed, regular files): the file is memory-mapped, and the only
//    work done under the lock is claiming a record-aligned byte range and counting its
//    records -- a newline scan, no allocation. The expensive part, building each read's
//    name/sequence strings, runs OFF the lock in the worker thread straight from the
//    mapped bytes. Keeping the serialized section tiny is what stops the single input
//    from capping multi-thread throughput.
//  - STREAM (gzip, a pipe, or stdin): falls back to one kseq parser run entirely under
//    the lock, batch_n records at a time -- the original behavior.
//
// One instance per run, shared by every worker thread; locks engage only when threaded.
class FastqSource {
    std::mutex mu;
    bool threaded;
    uint64_t served_ = 0;          // monotonic read index, for in-order output
    char fmt = 0;                  // '@' FASTQ or '>' FASTA

    // MMAP mode
    bool use_mmap = false;
    const char* data = nullptr;    // mapped file bytes
    size_t data_size = 0;
    size_t offset = 0;             // next unclaimed byte (guarded by mu)
    static constexpr size_t BLOCK = 1u << 18;  // ~256 KB target claim (record-aligned)

    // STREAM mode
    gzFile fp = nullptr;
    kseq_t* seq = nullptr;
    size_t batch_n;                // records per refill (see MOVI_CO_BATCH)


    // STREAM: parse up to batch_n records with kseq, under the lock.
    size_t fill(std::vector<CoReadData>& out) {
        out.clear();
        for (size_t k = 0; k < batch_n; k++) {
            int l = kseq_read(seq);
            if (l < 0) break;  // EOF (-1) or error (-2/-3)
            CoReadData rd;
            rd.valid = true;
            rd.name.assign(seq->name.s, seq->name.l);
            rd.sequence.assign(seq->seq.s, seq->seq.l);
            rd.seq = served_++;
            out.push_back(std::move(rd));
        }
        return out.size();
    }

    // MMAP, under the lock: from `offset`, advance to a record boundary at least BLOCK
    // bytes ahead (or EOF), returning the byte range [p,e) and the number of complete
    // records in it, bumping served_ and offset. Only a newline scan -- no allocation.
    // The record rule matches parse_range exactly (4 lines = one FASTQ record; each
    // line-initial '>' = one FASTA record), so the count equals what parse_range builds.
    bool claim(const char*& p, const char*& e, uint64_t& base, size_t& count) {
        if (offset >= data_size) return false;
        size_t start = offset, pos = offset, target = offset + BLOCK, n = 0;
        if (fmt == '@') {
            size_t lines = 0;
            while (pos < data_size) {
                const char* nl = static_cast<const char*>(memchr(data + pos, '\n', data_size - pos));
                pos = nl ? static_cast<size_t>(nl - data) + 1 : data_size;
                if (++lines % 4 == 0) { ++n; if (pos >= target) break; }
                if (!nl) break;
            }
        } else {  // FASTA: each record starts at a line-initial '>'
            while (pos < data_size) {
                ++n;                       // pos is at this record's '>'
                size_t q = pos + 1;
                while (q < data_size) {
                    const char* nl = static_cast<const char*>(memchr(data + q, '\n', data_size - q));
                    if (!nl) { q = data_size; break; }
                    size_t nq = static_cast<size_t>(nl - data) + 1;
                    if (nq < data_size && data[nq] == '>') { q = nq; break; }
                    q = nq;
                }
                pos = q;
                if (pos >= target) break;
            }
        }
        p = data + start; e = data + pos;
        base = served_; served_ += n; offset = pos; count = n;
        return n > 0;
    }

    // MMAP, off the lock: parse [p,e) into `out`, matching kseq's fields (name = bytes
    // after '@'/'>' up to the first whitespace; sequence = the sequence line(s), any
    // trailing '\r' stripped). Produces exactly the record count claim() reported.
    void parse_range(const char* p, const char* e, uint64_t base, std::vector<CoReadData>& out) {
        out.clear();
        uint64_t idx = base;
        if (fmt == '@') {
            while (p < e && *p == '@') {
                const char* n0 = static_cast<const char*>(memchr(p, '\n', e - p));      if (!n0) break;
                const char* s1 = n0 + 1;
                const char* n1 = static_cast<const char*>(memchr(s1, '\n', e - s1));    if (!n1) break;
                const char* s2 = n1 + 1;
                const char* n2 = static_cast<const char*>(memchr(s2, '\n', e - s2));    if (!n2) break;
                const char* s3 = n2 + 1;
                const char* n3 = static_cast<const char*>(memchr(s3, '\n', e - s3));    // may be null at EOF
                const char* nameB = p + 1, *nameE = nameB;
                while (nameE < n0 && *nameE != ' ' && *nameE != '\t' && *nameE != '\r') ++nameE;
                const char* se = n1; if (se > s1 && se[-1] == '\r') --se;
                CoReadData rd;
                rd.valid = true;
                rd.name.assign(nameB, nameE - nameB);
                rd.sequence.assign(s1, se - s1);
                rd.seq = idx++;
                out.push_back(std::move(rd));
                if (!n3) break;
                p = n3 + 1;
            }
        } else {  // FASTA (possibly multi-line sequence)
            while (p < e && *p == '>') {
                const char* n0 = static_cast<const char*>(memchr(p, '\n', e - p));      if (!n0) break;
                const char* nameB = p + 1, *nameE = nameB;
                while (nameE < n0 && *nameE != ' ' && *nameE != '\t' && *nameE != '\r') ++nameE;
                CoReadData rd;
                rd.valid = true;
                rd.name.assign(nameB, nameE - nameB);
                const char* s = n0 + 1;
                std::string sb;
                while (s < e && *s != '>') {
                    const char* nl = static_cast<const char*>(memchr(s, '\n', e - s));
                    const char* le = nl ? nl : e; if (le > s && le[-1] == '\r') --le;
                    sb.append(s, le - s);
                    if (!nl) { s = e; break; }
                    s = nl + 1;
                }
                rd.sequence = std::move(sb);
                rd.seq = idx++;
                out.push_back(std::move(rd));
                p = s;
            }
        }
    }

public:
    FastqSource(const string& fastq_file, bool threaded_, size_t batch_n_ = 32)
        : threaded(threaded_), batch_n(batch_n_ ? batch_n_ : 32) {
        // Prefer mmap for a raw, regular, uncompressed FASTQ/FASTA file; anything else
        // (gzip, stdin, a pipe) falls back to the kseq stream path.
        int fd = (fastq_file != "-") ? ::open(fastq_file.c_str(), O_RDONLY) : -1;
        if (fd >= 0) {
            struct stat st;
            if (fstat(fd, &st) == 0 && S_ISREG(st.st_mode) && st.st_size > 0) {
                void* m = mmap(nullptr, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
                if (m != MAP_FAILED) {
                    const unsigned char* b = static_cast<const unsigned char*>(m);
                    bool gz = (st.st_size >= 2 && b[0] == 0x1f && b[1] == 0x8b);
                    size_t i = 0;
                    while (i < static_cast<size_t>(st.st_size) &&
                           (b[i] == '\n' || b[i] == '\r' || b[i] == ' ' || b[i] == '\t')) ++i;
                    char f = (i < static_cast<size_t>(st.st_size)) ? static_cast<char>(b[i]) : 0;
                    if (!gz && (f == '@' || f == '>')) {
                        use_mmap = true;
                        data = static_cast<const char*>(m);
                        data_size = st.st_size;
                        offset = i;
                        fmt = f;
                        madvise(m, st.st_size, MADV_SEQUENTIAL);
                        ::close(fd);
                        return;
                    }
                    munmap(m, st.st_size);
                }
            }
            ::close(fd);
        }
        fp = gzopen(fastq_file.c_str(), "r");
        if (!fp) throw std::runtime_error("Cannot open FASTQ file: " + fastq_file);
        seq = kseq_init(fp);
        if (!seq) { gzclose(fp); throw std::runtime_error("Cannot initialize kseq"); }
    }

    ~FastqSource() {
        if (seq) kseq_destroy(seq);
        if (fp) gzclose(fp);
        if (use_mmap && data) munmap(const_cast<char*>(data), data_size);
    }

    // Refill `out` with the next batch. In MMAP mode only the claim (a record-boundary
    // scan) is serialized; the parse runs after the lock is released. Global read indices
    // stay contiguous because ranges are claimed in lock order.
    size_t next_batch(std::vector<CoReadData>& out) {
        if (use_mmap) {
            const char* p; const char* e; uint64_t base; size_t count;
            bool got;
            if (!threaded) {
                got = claim(p, e, base, count);
                if (!got) { out.clear(); return 0; }
                parse_range(p, e, base, out);
                return out.size();
            }
            {
                std::lock_guard<std::mutex> lk(mu);
                got = claim(p, e, base, count);
            }
            if (!got) { out.clear(); return 0; }
            parse_range(p, e, base, out);   // off the lock
            return out.size();
        }
        // STREAM: kseq under the lock.
        if (!threaded) {
            size_t n = fill(out);
            return n;
        }
        std::lock_guard<std::mutex> lk(mu);
        size_t n = fill(out);
        return n;
    }

};

// Per-thread view over the shared FastqSource. Holds this thread's current batch and
// its per-coroutine slots; it refills from the source (locked once per batch) and then
// hands reads to coroutines with zero per-read copy. Each worker thread owns one of
// these, so the per-read path is entirely thread-local. Coroutine bodies take this by
// reference, unchanged.
class SharedFastqReader {
public:
    using ReadData = CoReadData;

private:
    FastqSource* src;
    std::vector<CoReadData> batch;   // this thread's current batch of parsed records
    size_t batch_pos = 0;            // next unserved record within batch
    std::vector<CoReadData> slots;   // one slot per coroutine (zero-copy handoff)

public:
    SharedFastqReader(FastqSource& src_, int concurrency) : src(&src_) {
        slots.resize(concurrency);
    }

    // Move the next record into coroutine id's slot. Returns false at EOF
    // (the slot is marked invalid).
    bool serve_next(int id) {
        if (batch_pos >= batch.size()) { src->next_batch(batch); batch_pos = 0; }
        if (batch_pos >= batch.size()) { slots[id].valid = false; return false; }
        slots[id] = std::move(batch[batch_pos++]);
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
    
    int read_count = 0;  // only used by DEBUG_MSG_CO tracing below

    // Per-coroutine work buffer reused across reads (cleared, not reallocated,
    // each read) so there is no per-read heap churn.
    std::vector<uint16_t> matching_lens;

    while (true) {
        DEBUG_MSG_CO("DEBUG: Coroutine %d about to await next read\n", coroutine_id);
        // Request next read - this suspends until scheduler reads it
        auto& read_data = co_await get_next_read(reader, my_handle_storage, coroutine_id);
        
        DEBUG_MSG_CO("DEBUG: Coroutine %d resumed with read, valid=%d\n", coroutine_id, read_data.valid ? 1 : 0);
        if (!read_data.valid) break;  // EOF
        
        read_count++;

        // Process the read in place, without reversing it: query_pml walks the read
        // right-to-left (pos from length-1 down to 0), producing forward-read PMLs.
        // Reversing the read would instead produce the reversed-read PMLs (the
        // `--pml --reverse` result), which is a different query.
        std::string& R = read_data.sequence;
        matching_lens.clear();
        int32_t roff = R.length() - 1;    // offset in read
        uint64_t idx = r - 1;             // row index
        uint64_t match_len = 0;           // match length (consecurive case 1s) so far
        uint64_t offset = get_n(idx) - 1;
        assert(idx < rlbwt.size());
        assert(offset < get_n(idx));

        // The PML loop mirrors the sequential query_pml, using the same mode-portable
        // primitives (get_c / reposition_thresholds / get_id / get_n / get_offset) so it
        // works in every threshold index mode (regular-, sampled-, and blocked-thresholds),
        // not just mode 6. The single latency-hiding point is the LF destination row: its
        // id (get_id) is prefetched, then co_yield lets other coroutines run while that
        // cache line fills. Repositioning is a sequential run scan with no long-range
        // jump, so it stays a plain call with no yield. If you edit this loop, re-validate
        // byte-for-byte against `movi query --pml` (the
        // tests/regression_coroutine_separators.sh check does exactly this).
        uint64_t scan_count = 0;   // reposition scan length (needed by the reposition API)
        while (roff > -1) {
            DEBUG_MSG_CO("DEBUG: Coroutine %d processing, roff=%d\n", coroutine_id, roff);
            char row_c = alphabet[rlbwt[idx].get_c()];
            if (!check_alphabet(R[roff])) {  // char doens't exist in reference
                match_len = 0;
            } else if (row_c == R[roff]) {   // Case 1: Match
                match_len++;
            } else {                         // Case 2: reposition (up or down)
                // Repositioning is a sequential run scan with no long-range jump, so it
                // calls the shared, mode-portable primitive instead of an inlined mode-6
                // copy -- there is no prefetch/co_yield here. This is exactly what the
                // sequential query_pml does, so the PMLs stay byte-identical.
                bool up = movi_options->is_random_repositioning()
                            ? reposition_randomly(idx, offset, R[roff], scan_count)
                            : reposition_thresholds(idx, offset, R[roff], scan_count);
                match_len = 0;
                assert(idx < r && idx < rlbwt.size());
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
            
            // Begin LF step. The destination row id is the one long-range jump; get_id
            // is the mode-portable accessor (a direct field in mode 6, a block or tally
            // lookup in the blocked/sampled-threshold modes). Prefetch it, then yield.
            uint64_t new_idx = get_id(idx);
            my_prefetch_r((void*)(&(rlbwt[0]) + new_idx));
            offset = get_offset(idx) + offset;
            DEBUG_MSG_CO("DEBUG: Coroutine %d about to co_yield (prefetch)\n", coroutine_id);
            co_yield monostate{}; // wait for prefetch
            DEBUG_MSG_CO("DEBUG: Coroutine %d resumed after co_yield\n", coroutine_id);
            uint64_t n = get_n(new_idx);
            // Fast-forward past exhausted runs (the sequential query_pml does the same).
            while (new_idx < r - 1 && offset >= n) {
                offset -= n;
                new_idx++;
                n = get_n(new_idx);
            }
            idx = new_idx;
        }
        
        // Emit via the shared mainline output_base_stats so coroutine PML is
        // byte-identical to the sequential query_pml path (values, order, and format):
        // both walk the read right-to-left, so matching_lens is already in the order
        // query_pml records via add_ml. output_base_stats then formats it identically
        // (binary st_length+name+lengths to a file, or the reversed-string text form
        // to stdout).
        MoveQuery mq_out{std::string()};
        mq_out.set_query_id(read_data.name);
        const bool to_stdout = movi_options->is_stdout();
        for (uint16_t v : matching_lens)
            mq_out.add_ml(static_cast<uint64_t>(v), to_stdout);
        std::ostringstream oss;
        output_base_stats(DataType::match_length, to_stdout, oss, mq_out);
        g_emitter.emit(read_data.seq, oss.str());
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
                // Signed counter, as in query_mem_bml: the bound is negative when the
                // ftab already covers the whole seed (min-mem-length == ftab-k, the
                // zero-extend case), and a size_t one wraps to SIZE_MAX instead of
                // skipping the loop.
                for (int32_t j = 0; j <= init_pos - pos_on_r; ++j) {
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
                // Signed counter, same reason as the ftab_skip loop above.
                for (int32_t j = 0; j <= init_pos - pos_on_r; ++j) {
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

        // Format via the shared mainline output_mems so coroutine output is
        // identical to the sequential/strand path by construction; the reorder
        // buffer still emits in input order.
        std::ostringstream oss;
        output_mems(false, oss, mq);
        g_emitter.emit(read_data.seq, oss.str());
    }
    co_return;
}


// Coroutine k-mer query with latency hiding. Mirrors MoveStructure::query_all_kmers
// (together with query_kmers_from and look_ahead_backward_search), including both
// skip heuristics, so its output is byte-identical to `movi query --kmer
// [--kmer-count]`. It adds a prefetch + co_yield before each backward_search_step
// (the cache-missing LF_move) to hide memory latency. Those two inner backward-search
// loops are inlined here because suspend points cannot live in called functions. One
// coroutine serves both modes via count_mode. k-mer search uses plain (not
// bidirectional) backward search, so no reverse-complement index is needed and the
// read is not reversed.
MoveStructure::coroutine_task MoveStructure::query_kmer_coroutine(
    SharedFastqReader& reader,
    coroutine_handle<>& my_handle_storage,
    int coroutine_id) {

    const size_t  ftab_k = movi_options->get_ftab_k();
    const int32_t k = static_cast<int32_t>(movi_options->get_k());
    const bool    count_mode = movi_options->is_kmer_count();

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
            if (step < 0) step = 0;   // ftab_k >= k: no look-ahead room (mirrors query_all_kmers)

            while (pos_on_r >= 0 and static_cast<size_t>(pos_on_r) + 1 >= k) {
                bool did_search = true;

                if (step > 0 && static_cast<size_t>(pos_on_r) + 1 >= k + step) {
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
                    } while (match_len == 0 && pos_on_r >= 0
                             && static_cast<size_t>(pos_on_r) + 1 >= k && ftab_k > 1);

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
                        // Presence is `found` (0 or 1), matching the sequential count
                        // path (query_all_kmers): an absent k-mer passes count 0, which
                        // add_kmer treats as a no-op (no poslist entry, no found tally),
                        // while a present k-mer records kc as its occurrence multiplicity.
                        // The position is the leftmost k-mer of the found run, pos_on_r +
                        // 2 - k (pos_on_r has been advanced to that run's left end + k - 2),
                        // the same expression the non-count branch and the sequential
                        // query_kmers_from caller use.
                        mq.add_kmer(pos_on_r + 2 - k, found,
                                    std::numeric_limits<uint64_t>::max(), kc);
                    } else {
                        mq.add_kmer(pos_on_r + 2 - k, found);
                    }
                    }
                }

                while (pos_on_r >= 0 && !check_alphabet(query_seq[pos_on_r])) pos_on_r -= 1;
            }
        }

        // Format via the shared mainline output_kmers so coroutine output is
        // identical to the sequential/strand path by construction; the reorder
        // buffer still emits in input order.
        std::ostringstream oss;
        output_kmers(false, oss, all_kmer_count, mq, *movi_options);
        g_emitter.emit(read_data.seq, oss.str());
    }
    co_return;
}


// Entry point for the engine: movi.cpp's query() dispatch calls this on an
// already-deserialized index. The engine does not own index construction -- the
// caller loads the index and sets up ftab/kmerbv before handing it over here.
void run_coroutine_query(MoveStructure& mv, const string& fastq_file, int concurrency,
                         int nthreads, const CoroutineQueryOptions& opts, std::ostream& out) {
    init_buffered_io();
    if (nthreads < 1) nthreads = 1;
    if (concurrency < 1) concurrency = 1;
    const bool threaded = nthreads > 1;

    g_emitter.out = &out;
    // No per-run PML header is written here: PML output goes through the shared
    // output_base_stats (per-read). When writing to a file, the one-time BPFHeader is
    // written by open_output_files; stdout PML has no header. This matches the
    // sequential path, so coroutine PML output is byte-identical to it.
    g_emitter.ordered = opts.ordered_output;
    g_emitter.threaded = threaded;
    g_emitter.next_emit = 0;
    debug_enabled = opts.debug;

    auto start_time = steady_clock::now();

    // Per-thread work-claim / parse granularity, chosen by thread count. The batch is
    // both the parse unit and the unit of work a thread claims at once, so the two ends
    // pull opposite ways: a large batch amortizes the parse/refill (which a single thread
    // cannot hide behind other threads), while a small batch spreads work evenly across
    // threads (with N reads and B-read batches, at most N/B threads ever get work). A
    // single-thread run has no balancing to do, so it keeps the large batch and its lower
    // overhead; a multi-thread run takes the small batch, which a sweep found gives the
    // same peak throughput as any batch up to 64 while 256 costs ~17% at the optimum and
    // starves threads outright on small inputs. MOVI_CO_BATCH overrides both.
    size_t co_batch = (nthreads <= 1) ? 256 : 32;
    if (const char* e = std::getenv("MOVI_CO_BATCH")) {
        long v = std::atol(e);
        if (v > 0) co_batch = static_cast<size_t>(v);
    }

    // The one shared input; every worker thread pulls batches from it. Output ordering
    // is preserved because it assigns a single global input-order index per read.
    FastqSource source(fastq_file, threaded, co_batch);

    std::atomic<int64_t> total_bases{0};
    std::atomic<long> total_iterations{0};
    std::atomic<long> total_reads_served{0};

    // Each thread runs an independent scheduler over its own coroutines and its own
    // view of the input; the coroutine bodies read the shared, read-only index and emit
    // through the shared OrderedEmitter. `concurrency` is the per-thread in-flight
    // coroutine count (latency hiding within a thread), orthogonal to nthreads (cores).
    // The scheduler is query-agnostic -- it only resumes handles -- so PML, MEM and the
    // k-mer queries are all parallelized here by the same loop.
    #pragma omp parallel num_threads(nthreads)
    {
        SharedFastqReader reader(source, concurrency);
        std::vector<coroutine_handle<>> waiting_handles(concurrency);
        std::vector<MoveStructure::coroutine_task> coroutines;
        coroutines.reserve(concurrency);

        // Initialize all coroutines - use move semantics
        for (int i = 0; i < concurrency; ++i) {
            if (opts.mem_query)
                coroutines.push_back(std::move(mv.query_mem_coroutine(reader, waiting_handles[i], i)));
            else if (opts.kmer_query)
                coroutines.push_back(std::move(mv.query_kmer_coroutine(reader, waiting_handles[i], i)));
            else
                coroutines.push_back(std::move(mv.query_pml_coroutine(reader, waiting_handles[i], i)));
            // They start suspended, so resume to begin execution.
            if (coroutines[i].coro) {
                coroutines[i].coro.resume();
            }
        }

        // Per-thread scheduler loop: coordinate reads and coroutine execution
        int active_coroutines = concurrency;
        long iteration = 0;
        long reads_served = 0;
        int64_t bases = 0;

        while (active_coroutines > 0) {
            iteration++;
            active_coroutines = 0;

            for (int i = 0; i < concurrency; ++i) {
                if (coroutines[i].is_done()) {
                    continue;
                }

                // Check if this coroutine is waiting for a read
                if (waiting_handles[i]) {
                    reads_served++;
                    // Move the next batched record into this coroutine's slot.
                    if (reader.serve_next(i)) {
                        bases += reader.slot(i).sequence.length();
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

        total_bases += bases;
        total_iterations += iteration;
        total_reads_served += reads_served;
    }  // implicit barrier: all threads have finished before we drain the emitter

    // Emit any buffered output in input order and flush (part of the timed region).
    g_emitter.finish();

    // Calculate final statistics
    auto end_time = steady_clock::now();
    double total_elapsed = duration_cast<duration<double>>(end_time - start_time).count();

    // Calculate nanoseconds per base
    double nanoseconds_per_base = 0.0;
    if (total_bases.load() > 0) {
        double total_elapsed_ns = total_elapsed * 1e9; // Convert seconds to nanoseconds
        nanoseconds_per_base = total_elapsed_ns / static_cast<double>(total_bases.load());
    }

    fprintf(stderr_buf, "Scheduler completed: %ld iterations (summed over %d thread(s)), %ld reads served to coroutines\n",
            total_iterations.load(), nthreads, total_reads_served.load());
    fprintf(stderr_buf, "Completed processing (elapsed: %.2f sec, %.1f ns/base, bases=%ld)\n",
            total_elapsed, nanoseconds_per_base, static_cast<long>(total_bases.load()));

}


#endif  // MODE == 6 && COLOR_MODE == 0
