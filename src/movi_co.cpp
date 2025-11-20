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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <variant>
#include <stdexcept>

#include <sdsl/int_vector.hpp>
#include <zlib.h>

#include "utils.hpp"
#include "move_structure.hpp"
#include "move_query.hpp"
#include "movi_options.hpp"
#include <coroutine>

using std::suspend_always;
using std::coroutine_handle;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::numeric_limits;
using std::monostate;
using std::exception;
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

// Shared reader that manages single gzFile and kseq_t
class SharedFastqReader {
public:
    struct ReadData {
        bool valid;
        string name;
        string sequence;
    };
    
private:
    gzFile fp;
    kseq_t* seq;
    ReadData pending_read;
    int read_call_count = 0;  // Add this
    
public:
    SharedFastqReader(const string& fastq_file) : fp(nullptr), seq(nullptr) {
        fp = gzopen(fastq_file.c_str(), "r");
        if (!fp) {
            throw std::runtime_error("Cannot open FASTQ file: " + fastq_file);
        }
        seq = kseq_init(fp);
        if (!seq) {
            gzclose(fp);
            throw std::runtime_error("Cannot initialize kseq");
        }
    }
    
    ~SharedFastqReader() {
        if (seq) kseq_destroy(seq);
        if (fp) gzclose(fp);
    }
    
    // Read next sequence into pending_read
    bool read_next() {
        read_call_count++;
        cerr << "DEBUG: read_next() called (call #" << read_call_count << ")" << endl;
        int l = kseq_read(seq);
        if (l < 0) {
            cerr << "DEBUG: kseq_read returned " << l << " on call #" << read_call_count << endl;
            // kseq_read returns:
            // -1: EOF
            // -2: truncated quality string  
            // -3: error reading stream
            pending_read.valid = false;
            return false;
        }
        cerr << "DEBUG: kseq_read returned length " << l << " on call #" << read_call_count << endl;
        pending_read.valid = true;
        pending_read.name = string(seq->name.s);
        pending_read.sequence = string(seq->seq.s);
        return true;
    }
    
    const ReadData& get_pending_read() const { return pending_read; }
};

// Awaitable for requesting a read from the shared reader
struct read_awaitable {
    SharedFastqReader* reader;
    coroutine_handle<>* stored_handle;  // Where to store the coroutine handle
    
    bool await_ready() const noexcept { return false; }
    
    template<typename Promise>
    void await_suspend(coroutine_handle<Promise> h) {
        cerr << "DEBUG: await_suspend storing handle for coroutine" << endl;
        *stored_handle = h;
        cerr << "DEBUG: handle stored, stored_handle pointer is " 
             << (stored_handle ? "non-null" : "null") << ", handle value is "
             << (h ? "non-null" : "null") << endl;
    }
    
    SharedFastqReader::ReadData await_resume() {
        return reader->get_pending_read();
    }
};

// Helper to create awaitable
read_awaitable get_next_read(SharedFastqReader& reader, coroutine_handle<>& handle_storage) {
    read_awaitable awaiter;
    awaiter.reader = &reader;
    awaiter.stored_handle = &handle_storage;
    return awaiter;
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
            cerr << "DEBUG: Coroutine exception caught in unhandled_exception!" << endl;
            try {
                throw;  // Re-throw to see what it is
            } catch (const std::exception& e) {
                cerr << "DEBUG: Exception type: std::exception, what(): " << e.what() << endl;
            } catch (...) {
                cerr << "DEBUG: Exception type: unknown" << endl;
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
    
    while (true) {
        cerr << "DEBUG: Coroutine " << coroutine_id << " about to await next read" << endl;
        // Request next read - this suspends until scheduler reads it
        auto read_data = co_await get_next_read(reader, my_handle_storage);
        
        cerr << "DEBUG: Coroutine " << coroutine_id << " resumed with read, valid=" 
             << read_data.valid << endl;
        if (!read_data.valid) break;  // EOF
        
        read_count++;
        int read_bases = read_data.sequence.length();
        total_bases += read_bases;
        
        // Process the read
        string query_seq = read_data.sequence;
        std::reverse(query_seq.begin(), query_seq.end());
        MoveQuery mq(query_seq);

        auto& R = mq.query();
        int32_t roff = R.length() - 1;    // offset in read
        uint64_t idx = r - 1;             // row index
        uint64_t offset = get_n(idx) - 1; // offset into row
        uint64_t match_len = 0;           // match length (consecurive case 1s) so far
        uint64_t ff_count_tot = 0, scan_count = 0;
        
        while (roff > -1) {
            cerr << "DEBUG: Coroutine " << coroutine_id << " processing, roff=" << roff << endl;
            char row_c = alphabet[rlbwt[idx].get_c()];
            if (!check_alphabet(R[roff])) {  // char doens't exist in reference
                match_len = 0;
            } else if (row_c == R[roff]) {   // Case 1: Match
                match_len++;
            } else {                         // Case 2: Reposition
                uint64_t saved_idx = idx;
                uint64_t alphabet_index = alphamap[static_cast<uint64_t>(R[roff])];
                bool up = false;
                if (idx == end_bwt_idx) {
                    up = offset < end_bwt_idx_thresholds[alphabet_index];
                    if (up) {
                        idx = reposition_up(saved_idx, R[roff], scan_count);
                        assert(idx < saved_idx);
                    } else {
                        idx = reposition_down(saved_idx, R[roff], scan_count);
                        assert(idx > saved_idx);
                    }
                } else {
                    alphabet_index = alphamap_3[alphamap[row_c]][alphabet_index];                
                    up = offset < get_thresholds(idx, alphabet_index);
                    if (up) {
                        idx = reposition_up(saved_idx, R[roff], scan_count);
                        assert(idx < saved_idx);
                    } else {
                        idx = reposition_down(saved_idx, R[roff], scan_count);
                        assert(idx > saved_idx);
                    }
                }
                match_len = 0;
                assert(alphabet[rlbwt[idx].get_c()] == R[roff] && "Repositioning failed - character mismatch");
                offset = up ? get_n(idx) - 1 : 0;
            }
            // At this point, if repositioning was needed it has been done
            // 'idx' and 'offset' are the new move-structure row index and offset
            
            // Add matching length to query
            mq.add_ml(match_len, movi_options->is_stdout());
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
            cerr << "DEBUG: Coroutine " << coroutine_id << " about to co_yield (prefetch)" << endl;
            co_yield monostate{}; // wait for prefetch
            cerr << "DEBUG: Coroutine " << coroutine_id << " resumed after co_yield" << endl;
            if (new_idx < r - 1 && offset >= get_n(new_idx)) {
                uint64_t niter = 0;
                while (new_idx < r - 1 && offset >= get_n(new_idx)) {
                    offset -= get_n(new_idx);
                    niter++;
                    new_idx++;
                }
                assert(niter < numeric_limits<uint16_t>::max());
                ff_count_tot += static_cast<uint16_t>(niter);
            }
            idx = new_idx;
        }
        
        // Output results for this read
        cout << read_data.name << " ";
        auto& matching_lengths = mq.get_matching_lengths();
        for (size_t i = 0; i < matching_lengths.size(); i++) {
            cout << matching_lengths[i];
            if (i < matching_lengths.size() - 1) {
                cout << " ";
            }
        }
        cout << endl;
        cerr << "DEBUG: Coroutine " << coroutine_id << " finished processing read " 
             << read_count << ", looping back" << endl;
    }
    
    co_return;
}


void process_fastq(const string& fastq_file, const string& index_dir, int concurrency) {
    MoviOptions movi_options;
    movi_options.set_index_dir(index_dir);
    movi_options.set_pml(); // Enable PML mode
    
    MoveStructure mv(&movi_options);
    mv.deserialize();
    cerr << "Successfully loaded Movi index from: " << index_dir << endl;
    
    cout << "# Read_ID PML_Values" << endl;
    
    auto start_time = steady_clock::now();
    int64_t total_bases = 0;
    
    // Single shared reader for all coroutines
    SharedFastqReader reader(fastq_file);
    
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
                cerr << "DEBUG: Coroutine " << i << " is done" << endl;
                continue;
            }
            
            // Check if this coroutine is waiting for a read
            if (waiting_handles[i]) {
                cerr << "DEBUG: Scheduler detected waiting_handles[" << i << "] is set" << endl;
                total_reads_served++;
                // Read next sequence and resume the coroutine
                bool read_success = reader.read_next();
                if (read_success) {
                    // Track total bases processed
                    const auto& read_data = reader.get_pending_read();
                    total_bases += read_data.sequence.length();
                    
                    waiting_handles[i].resume();
                    waiting_handles[i] = nullptr;  // Clear the handle
                    if (!coroutines[i].is_done()) {
                        active_coroutines++;
                    }
                } else {
                    // EOF - resume to let coroutine know there are no more reads
                    cerr << "DEBUG: EOF detected after serving " << total_reads_served 
                         << " reads to coroutines" << endl;
                    waiting_handles[i].resume();
                    waiting_handles[i] = nullptr;
                    if (!coroutines[i].is_done()) {
                        active_coroutines++;
                    }
                }
            } else {
                cerr << "DEBUG: Scheduler iteration " << iteration << ", coroutine " << i 
                     << " not waiting (handle is " << (waiting_handles[i] ? "set" : "null") << ")" << endl;
                // Not waiting for read - just resume (handles co_yield for voluntary suspension)
                if (coroutines[i].coro) {
                    coroutines[i].coro.resume();
                    if (!coroutines[i].is_done()) {
                        active_coroutines++;
                    }
                }
            }
        }
        
        // Safety check
        if (iteration > 100000 && total_reads_served < 2) {
            cerr << "WARNING: After " << iteration << " iterations, only " 
                 << total_reads_served << " reads served. Something is wrong." << endl;
            cerr << "waiting_handles: ";
            for (int j = 0; j < concurrency; ++j) {
                cerr << (waiting_handles[j] ? "1" : "0") << " ";
            }
            cerr << endl;
            break;
        }
    }
    
    cerr << "Scheduler completed: " << iteration << " iterations, " 
         << total_reads_served << " reads served to coroutines" << endl;
    
    // Calculate final statistics
    auto end_time = steady_clock::now();
    double total_elapsed = duration_cast<duration<double>>(end_time - start_time).count();
    
    // Calculate nanoseconds per base
    double nanoseconds_per_base = 0.0;
    if (total_bases > 0) {
        double total_elapsed_ns = total_elapsed * 1e9; // Convert seconds to nanoseconds
        nanoseconds_per_base = total_elapsed_ns / total_bases;
    }
    
    cerr << "Completed processing (elapsed: " << std::fixed << std::setprecision(2)
         << total_elapsed << " sec, " << std::setprecision(1) << nanoseconds_per_base
         << " ns/base, bases=" << total_bases << ")" << endl;
}

void print_usage(const char* program_name) {
    cerr << "Usage: " << program_name << " <fastq_file> <index_dir> [concurrency]" << endl;
    cerr << "  fastq_file: Input FASTQ file with reads" << endl;
    cerr << "  index_dir:  Directory containing the Movi index" << endl;
    cerr << "  concurrency: Number of concurrent coroutines (default: 1)" << endl;
    cerr << endl;
    cerr << "This program computes Pseudo Matching Lengths (PMLs) for reads" << endl;
    cerr << "using a Movi regular-thresholds index (MODE=6)." << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 4) {
        cerr << "argc = " << argc << endl;
        print_usage(argv[0]);
        return 1;
    }
    
    string fastq_file{argv[1]}, index_dir{argv[2]};
    int concurrency = 1;
    
    if (argc == 4) {
        try {
            concurrency = std::stoi(argv[3]);
            if (concurrency < 1) {
                cerr << "Error: Concurrency must be at least 1" << endl;
                return 1;
            }
        } catch (const exception& e) {
            cerr << "Error: Invalid concurrency value: " << argv[3] << endl;
            return 1;
        }
    }
    
    try {
        process_fastq(fastq_file, index_dir, concurrency);
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
