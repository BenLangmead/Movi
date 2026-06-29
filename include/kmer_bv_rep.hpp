#ifndef KMER_BV_REP_HPP
#define KMER_BV_REP_HPP

// Compile-time indirection for the succinct bitvector that backs the k-mer-id
// boundary structure B_k (built by build_kmerbv, queried by the --kmer-bv id
// path). B_k is query-time critical, so its representation is isolated behind
// these type aliases: swapping the representation for the whole engine is a
// one-line change here, with no edits to the build or query code.
//
// The default is an Elias-Fano sd_vector. B_k is sparse -- only canonical
// k-mers are marked -- so the sd_vector is several times smaller than a dense
// bit_vector + rank, and carries its own rank/select. To evaluate a different
// representation, redefine these three aliases (e.g. to sdsl::bit_vector with
// rank_support_v / select_support_mcl) and rebuild; the rest of the kmerbv
// subsystem is written against the aliases alone.

#include "sdsl_wrapper.hpp"

namespace movi {

using KmerIdBv     = sdsl::sd_vector<>;
using KmerIdBvRank = sdsl::sd_vector<>::rank_1_type;
using KmerIdBvSel  = sdsl::sd_vector<>::select_1_type;

} // namespace movi

#endif
