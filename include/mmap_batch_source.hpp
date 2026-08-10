// mmap_batch_source.hpp: the shared, memory-mapped input for the strand and sequential
// query paths.
//
// Both of those paths drive a BatchLoader, which used to read straight from one shared
// std::ifstream under an `omp critical` section. That put the whole parse -- getline per
// line, plus building the batch string -- inside the serialized region, so at high thread
// counts the single reader, not the index, set the ceiling.
//
// This class splits that work in two. Under the lock it only *claims* a record-aligned
// byte range out of the mapped file: a newline scan with no allocation. The claiming
// thread then hands that range to BatchLoader::loadBatchFromRange off the lock, so the
// expensive part runs in the worker. The coroutine path already works this way (see
// FastqSource in coroutine_processor.cpp); this brings the same split to every other
// query type, notably ZML, whole-read exact count, and MPHF-id.
//
// Only raw, regular, uncompressed FASTA/FASTQ files are mapped. Gzip, stdin and pipes
// report is_mmap() == false, and the caller keeps its original under-lock ifstream path.

#ifndef MMAP_BATCH_SOURCE_H
#define MMAP_BATCH_SOURCE_H

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

class MmapBatchSource {
public:
    explicit MmapBatchSource(const std::string& read_file) {
        if (read_file == "-") return;                 // stdin: caller keeps the stream path
        int fd = ::open(read_file.c_str(), O_RDONLY);
        if (fd < 0) return;
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
                    use_mmap_ = true;
                    data_ = static_cast<const char*>(m);
                    data_size_ = static_cast<size_t>(st.st_size);
                    offset_ = i;
                    fmt_ = f;
                    madvise(m, st.st_size, MADV_SEQUENTIAL);
                    ::close(fd);
                    return;
                }
                munmap(m, st.st_size);
            }
        }
        ::close(fd);
    }

    ~MmapBatchSource() {
        if (use_mmap_ && data_) munmap(const_cast<char*>(data_), data_size_);
    }

    MmapBatchSource(const MmapBatchSource&) = delete;
    MmapBatchSource& operator=(const MmapBatchSource&) = delete;

    bool is_mmap() const { return use_mmap_; }
    char format() const { return fmt_; }

    // Claim the next record-aligned byte range [p, e) and advance the shared offset.
    // Returns false at end of input. The caller serializes this; the parse that follows
    // does not need to be.
    //
    // The stopping rule reproduces BatchLoader::loadBatch's own batch-size heuristic --
    // keep taking records until both the base target and the read floor are met, where a
    // FASTQ record contributes half its bytes as an estimate of its sequence length and a
    // FASTA record contributes all of them -- so claim granularity, and therefore the
    // per-thread work-claim size that MOVI_STRAND_BATCH tunes, matches the stream path.
    bool claim(const char*& p_out, const char*& e_out, size_t num_bases, size_t min_reads) {
        if (offset_ >= data_size_) return false;
        const size_t start = offset_;
        size_t pos = offset_, bases = 0, reads = 0, record_size = 0, lines = 0;
        while (pos < data_size_ && (bases < num_bases || reads < min_reads)) {
            const char* nl = static_cast<const char*>(memchr(data_ + pos, '\n', data_size_ - pos));
            const size_t line_end = nl ? static_cast<size_t>(nl - data_) : data_size_;
            record_size += line_end - pos;              // line length, excluding the newline
            pos = nl ? line_end + 1 : data_size_;
            ++lines;
            if (fmt_ == '@') {
                if (lines % 4 == 0) { bases += record_size / 2; record_size = 0; ++reads; }
            } else {  // FASTA: a record ends where the next line-initial '>' begins
                if (pos >= data_size_ || data_[pos] == '>') { bases += record_size; record_size = 0; ++reads; }
            }
            if (!nl) break;
        }
        p_out = data_ + start;
        e_out = data_ + pos;
        offset_ = pos;
        return pos > start;
    }

private:
    bool use_mmap_ = false;
    const char* data_ = nullptr;
    size_t data_size_ = 0;
    size_t offset_ = 0;      // next unclaimed byte; guarded by the caller's lock
    char fmt_ = 0;           // '@' FASTQ or '>' FASTA
};

#endif /* End of MMAP_BATCH_SOURCE_H */
