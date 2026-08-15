// The following code is borrowed from SPUMONI by Omar Y. Ahmed with minor modifications
// https://github.com/oma219/spumoni/blob/main/src/batch_loader.cpp

/*
  * File: batch_loader.cpp
  * Description: Loads in a block of reads that can be then 
  *              processed by SPUMONI
  *             
  * Start Date: April 21, 2022
  *
  * Note: The BatchLoader class is inspired by the implementation
  *       in the Kraken2 repository at this link 
  *       (https://github.com/DerrickWood/kraken2/blob/master/src/seqreader.cpp). It uses
  *       different logic to determine the batch size to load balance between
  *       multiple reader threads.
  */

#include <batch_loader.hpp>

BatchLoader::BatchLoader() {
    /* Main constructor for BatchLoader */
    input_format = NOT_CLEAR;
    batch_buffer.reserve(8192);
}

bool BatchLoader::loadBatch(std::ifstream& input, size_t num_bases, size_t min_number_of_reads) {
    /* Attempts to load a batch of reads that have at least num_bases in them */

    // figure out the type of input file (if needed)
    if (input_format == NOT_CLEAR) {
        if (!input.good()) return false;
        switch (input.peek()) {
            case '>': input_format = FA; break;
            case '@': input_format = FQ; break;
            case EOF: return false; break;
            default: FATAL_ERROR("unrecognized input query file type - expects FASTA or FASTQ.");
        }
    }

    batch_stream.clear();
    batch_stream.str("");

    // reads in a batch of reads until we hit the required number of bases
    size_t num_bases_covered = 0;
    size_t num_reads_covered = 0;
    size_t lines_covered = 0;
    size_t record_size = 0;
    bool valid_batch = false;

    while (input.good() && (num_bases_covered < num_bases || num_reads_covered < min_number_of_reads)) {
        
        if (!std::getline(input, batch_buffer)) {
            // Ensure the last batch is read
            if (input_format == FQ and lines_covered % 4 == 0) {
                return true;
            } else if (input_format == FA and lines_covered % 2 == 0) {
                return true;
            } else {
                return false;
            }
        }

        lines_covered++;

        record_size += batch_buffer.size();
        valid_batch = true;

        // Use heuristic to determine number of bases read so far. Half of the
        // length of FASTQ records is roughly length of read, and full length
        // of FASTA record is roughly the length of read. Assuming the header
        // line is adding negligible length to record length.

        if (input_format == FQ) {
            if (lines_covered % 4 == 0) {
                num_bases_covered += record_size/2;
                record_size = 0;
                num_reads_covered += 1;
            }
        } else {
            if (input.peek() == '>') {
                num_bases_covered += record_size;
                record_size = 0;
                num_reads_covered += 1;
            }
        }
        batch_stream << batch_buffer << '\n';
    }
    return valid_batch;
}

bool BatchLoader::loadBatchFromRange(const char* begin, const char* end, char fmt) {
    /* Loads a batch from a record-aligned byte range claimed out of a mapped file.
       The range already ends on a record boundary, so unlike loadBatch there is no
       scanning to do here: the bytes are handed to the same batch_stream that
       grabNextRead consumes. */

    if (input_format == NOT_CLEAR) {
        switch (fmt) {
            case '>': input_format = FA; break;
            case '@': input_format = FQ; break;
            default: FATAL_ERROR("unrecognized input query file type - expects FASTA or FASTQ.");
        }
    }

    if (begin >= end) return false;

    batch_stream.clear();
    batch_stream.str(std::string(begin, end));
    return true;
}

bool BatchLoader::grabNextRead(Read& curr_read) {
    /* Grabs a read from the batch loaded into the BatchLoader object */

    ASSERT((input_format == FQ || input_format == FA), 
           "input query file type is not recognized when processing.");

    // grab the id for the read
    curr_read.read_format = this->input_format;
    if (!std::getline(batch_stream, batch_buffer)) return false;
    if (!batch_buffer.size()) return false; // an empty line 

    if (this->input_format == FQ) {
        if (batch_buffer[0] != '@') {
            FATAL_ERROR("Incorrect FASTQ entry, it should start with '@' but found %c", 
                         batch_buffer[0]);}
    } else { // FASTA
        if (batch_buffer[0] != '>') {
            FATAL_ERROR("Incorrect FASTA entry, it should start with '>' but found %c",
                        batch_buffer[0]);}
    }
    curr_read.header_line.assign(batch_buffer);

    // grab the id from the header line
    if (batch_buffer.size() <= 2) 
        FATAL_ERROR("header line is missing an id. invalid query cannot be processed.");

    // The id is the header up to the first whitespace, without the leading sigil.
    // id_end is an index, so the count passed to assign is one less than it: taking
    // id_end characters from position 1 would carry the delimiter into the id, and
    // the id is written into the query output.
    size_t id_end = batch_buffer.find_first_of(" \t\r", 1);
    if (id_end == std::string::npos) id_end = batch_buffer.size();
    curr_read.id.assign(batch_buffer, 1, id_end - 1);

    // grab the sequence
    if (this->input_format == FQ) {
        // get seq, and clean of trailing white space
        if (!std::getline(batch_stream, batch_buffer)) return false;
        while (std::isspace(batch_buffer.back())) batch_buffer.pop_back();
        curr_read.seq.assign(batch_buffer);
        
        // skip '+' line, and get qualities
        if (!std::getline(batch_stream, batch_buffer)) return false;
        if (!std::getline(batch_stream, batch_buffer)) return false;
        while (std::isspace(batch_buffer.back())) batch_buffer.pop_back();
        curr_read.qual.assign(batch_buffer);

    } else { // FASTA
        curr_read.qual.assign("");
        curr_read.seq.assign("");
        while (batch_stream.good() && batch_stream.peek() != '>') {
            if (!std::getline(batch_stream, batch_buffer)) return curr_read.seq.size();
            while (std::isspace(batch_buffer.back())) batch_buffer.pop_back();
            curr_read.seq.append(batch_buffer);
        }
    }    
    return true;
}