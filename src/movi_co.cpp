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
    };

private:
    gzFile fp;
    kseq_t* seq;
    std::vector<ReadData> batch;   // current batch of parsed records
    size_t batch_pos = 0;          // next unserved record within batch
    std::vector<ReadData> slots;   // one slot per coroutine (zero-copy handoff)
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
struct MoveStructure::query_pml_coroutine_return_type {
    struct promise_type {
        query_pml_coroutine_return_type get_return_object() {
            return query_pml_coroutine_return_type{coroutine_handle<promise_type>::from_promise(*this)};
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
    
    query_pml_coroutine_return_type() : coro(nullptr) {}
    query_pml_coroutine_return_type(coroutine_handle<promise_type> h) : coro(h) {}
    
    // Copy constructor - delete it since coroutine handles shouldn't be copied
    query_pml_coroutine_return_type(const query_pml_coroutine_return_type& other) = delete;
    
    // Move constructor
    query_pml_coroutine_return_type(query_pml_coroutine_return_type&& other) noexcept 
        : coro(other.coro) {
        other.coro = nullptr;
    }
    
    // Copy assignment - delete it
    query_pml_coroutine_return_type& operator=(const query_pml_coroutine_return_type& other) = delete;
    
    // Move assignment
    query_pml_coroutine_return_type& operator=(query_pml_coroutine_return_type&& other) noexcept {
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
    
    ~query_pml_coroutine_return_type() { 
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
MoveStructure::query_pml_coroutine_return_type MoveStructure::query_pml_coroutine(
    SharedFastqReader& reader, 
    coroutine_handle<>& my_handle_storage,
    int coroutine_id) {
    
    int read_count = 0;
    int64_t total_bases = 0;
    const bool use_separator = this->use_separator();
    const int sep_adjust = use_separator ? 1 : 0;
    
    // Per-coroutine output buffer (64KB should be enough for many lines)
    static constexpr size_t OUTPUT_BUFFER_SIZE = 65536;
    char output_buffer[OUTPUT_BUFFER_SIZE];
    size_t buffer_pos = 0;
    int buffer_line_count = 0;

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
        out_line.append(read_data.name);
        out_line.push_back(' ');
        for (size_t i = 0; i < matching_lens.size(); i++) {
            size_t w = write_uint_to_buffer(numbuf, sizeof(numbuf), matching_lens[i]);
            out_line.append(numbuf, w);
            if (i + 1 < matching_lens.size()) out_line.push_back(' ');
        }
        out_line.push_back('\n');

        // Emit: flush the shared buffer if the line wouldn't fit; write
        // oversized lines directly so nothing is ever truncated.
        if (out_line.size() > OUTPUT_BUFFER_SIZE) {
            std::lock_guard<std::mutex> lock(output_mutex);
            if (buffer_pos > 0) { fwrite(output_buffer, 1, buffer_pos, stdout_buf); buffer_pos = 0; }
            fwrite(out_line.data(), 1, out_line.size(), stdout_buf);
        } else {
            if (buffer_pos + out_line.size() > OUTPUT_BUFFER_SIZE) {
                std::lock_guard<std::mutex> lock(output_mutex);
                fwrite(output_buffer, 1, buffer_pos, stdout_buf);
                buffer_pos = 0;
            }
            std::memcpy(output_buffer + buffer_pos, out_line.data(), out_line.size());
            buffer_pos += out_line.size();
        }
        DEBUG_MSG_CO("DEBUG: Coroutine %d finished processing read %d, looping back\n", coroutine_id, read_count);
    }
    
    // Final flush of any remaining buffered output before returning
    if (buffer_pos > 0) {
        std::lock_guard<std::mutex> lock(output_mutex);
        fwrite(output_buffer, 1, buffer_pos, stdout_buf);
        fflush(stdout_buf);
    }
    
    co_return;
}


void process_fastq(const string& fastq_file, const string& index_dir, int concurrency) {
    MoviOptions movi_options;
    movi_options.set_index_dir(index_dir);
    movi_options.set_pml(); // Enable PML mode
    
    MoveStructure mv(&movi_options);
    mv.deserialize();
    fprintf(stderr_buf, "Successfully loaded Movi index from: %s\n", index_dir.c_str());
    
    fprintf(stdout_buf, "# Read_ID PML_Values\n");
    
    auto start_time = steady_clock::now();
    int64_t total_bases = 0;
    
    // Single shared reader for all coroutines
    SharedFastqReader reader(fastq_file, concurrency);
    
    // Storage for coroutine handles waiting for reads
    std::vector<coroutine_handle<>> waiting_handles(concurrency);
    std::vector<MoveStructure::query_pml_coroutine_return_type> coroutines;
    coroutines.reserve(concurrency);
    
    // Initialize all coroutines - use move semantics
    for (int i = 0; i < concurrency; ++i) {
        coroutines.push_back(std::move(mv.query_pml_coroutine(reader, waiting_handles[i], i)));
        // Resume coroutines initially to get them started
        // They start suspended, so we need to resume them to begin execution
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
                waiting_handles[i].resume();
                waiting_handles[i] = nullptr;
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
    fprintf(stderr_buf, "Usage: %s [--debug] <fastq_file> <index_dir> [concurrency]\n", program_name);
    fprintf(stderr_buf, "  --debug:     Enable debug output\n");
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
    
    // Parse --debug flag if present
    if (argc > 1 && string(argv[arg_idx]) == "--debug") {
        debug_enabled = true;
        arg_idx++;
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


