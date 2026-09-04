#ifndef MOVE_QUERY_HPP
#define MOVE_QUERY_HPP

#include <algorithm>
#include <vector>
#include <utility>

class MoveQuery {
    public:
        MoveQuery() {}
        MoveQuery(std::string q) : query_string(q) {}

        std::string& query() { return query_string; }
        uint64_t length() { return query_string.length(); }

        // One record per present query k-mer, in the order the search found them.
        // The search walks a read right to left, so positions arrive descending; a
        // consumer that needs input order scatters the records by pos.
        struct kmer_hit_t {
            int32_t pos;     // start position of the k-mer in the read
            uint64_t count;  // the value shown in the count field
            uint64_t id;     // MPHF id, or absent_id when no bitvector id was resolved
        };
        static constexpr uint64_t absent_id = std::numeric_limits<uint64_t>::max();

        // kmer_count drives the read-level presence tally (found_kmer_count =
        // number of present k-mers).  display_count, when provided, is what is
        // printed in the per-k-mer count field instead of kmer_count — used by
        // the --kmer-bv path to show the true occurrence multiplicity while
        // still counting presence one k-mer at a time.
        // orientation: -1 = unknown/not-applicable (per-orientation ids); 0 = the
        // query k-mer is the canonical (forward) one; 1 = the query k-mer is the
        // reverse complement of the canonical k-mer that owns the id. When set, it
        // is emitted as a 4th field (':f' / ':r'), matching SSHash's
        // lookup_result.kmer_orientation.
        void add_kmer(int32_t pos_on_r, uint64_t kmer_count,
                      uint64_t kmer_id = std::numeric_limits<uint64_t>::max(),
                      uint64_t display_count = std::numeric_limits<uint64_t>::max(),
                      int orientation = -1) {
            if (kmer_count > 0) {
                found_kmer_count += kmer_count;
                uint64_t shown = (display_count != std::numeric_limits<uint64_t>::max())
                                 ? display_count : kmer_count;
                // Structured record for the non-native output formats: the k-mer's
                // start position, its shown count, and its id where one was resolved.
                kmer_hits.push_back(kmer_hit_t{pos_on_r, shown, kmer_id});
                if (kmer_id != std::numeric_limits<uint64_t>::max()) {
                    matching_lengths_string += std::to_string(pos_on_r) + ":" +
                                               std::to_string(shown) + ":" +
                                               std::to_string(kmer_id);
                    if (orientation >= 0)
                        matching_lengths_string += (orientation == 0 ? ":f" : ":r");
                    matching_lengths_string += " ";
                } else {
                    matching_lengths_string += std::to_string(pos_on_r) + ":" +
                                               std::to_string(shown) + " ";
                }
            }
        }

        void add_color(uint64_t color_id) {
            // Check if the color id is ever larger than 2^32?

            matching_colors.push_back(color_id);
        }

        void add_ml(uint64_t len_, bool to_stdout) {
            uint16_t len = static_cast<uint16_t>(len_);
            if (len_ > std::numeric_limits<uint16_t>::max()) {
                len = std::numeric_limits<uint16_t>::max();
            }

            matching_lens.push_back(len);
            if (to_stdout) {
                std::string len_str = std::to_string(len);
                std::reverse(len_str.begin(), len_str.end());
                matching_lengths_string += " " + len_str;
            }
        }
        void add_sa_entries(uint64_t sa_entry) {
            sa_entries.push_back(sa_entry);
        }
        struct mem_t {
            uint32_t start;
            uint32_t end;
            // Occurrences of the MEM in the indexed text, as returned by
            // MoveInterval::count, which is 64-bit. A conserved MEM occurs once per
            // sequence in the collection, so on a large bacterial collection this
            // exceeds any 16- or 32-bit bound; the field matches the source type so
            // the value cannot wrap on the way in.
            uint64_t count;
        };
        void add_mem(uint32_t start, uint32_t end, uint64_t count) {
            mems.push_back(mem_t(start, end, count));
        }

        void add_cost(std::chrono::nanoseconds cost) { costs.push_back(cost); }
        void add_scan(uint64_t scan) { scans.push_back(scan); }
        void add_fastforward(uint64_t fastforward) { fastforwards.push_back(fastforward); }
        std::vector<uint16_t>& get_matching_lengths() { return matching_lens; }
        std::vector<uint64_t>& get_matching_colors() { return matching_colors; }
        std::vector<uint64_t>& get_sa_entries() { return sa_entries; }
        std::string& get_matching_lengths_string() { return matching_lengths_string; }
        std::vector<kmer_hit_t>& get_kmer_hits() { return kmer_hits; }
        std::vector<mem_t>& get_mems() { return mems; }
        std::vector<uint16_t>& get_scans() { return scans; }
        std::vector<uint16_t>& get_fastforwards() { return fastforwards; }
        std::vector<std::chrono::nanoseconds>& get_costs() { return costs; }

        void set_query_id(std::string id) { query_id = id; }
        std::string& get_query_id() { return query_id; }

        friend std::ostream& operator<<(std::ostream& output, const MoveQuery& dt) {
            // output << "The matching statistics are:\n";
            for (int64_t i = dt.matching_lens.size() - 1; i >= 0; i--) {
                auto pml = dt.matching_lens[i];
                output << pml << " ";
            }
            return output;
        }
        uint64_t found_kmer_count = 0;
    private:
        std::string query_id = "";
        std::string query_string = "";
        std::string matching_lengths_string = "";
        std::vector<kmer_hit_t> kmer_hits;  // one entry per present k-mer
        std::vector<uint64_t> sa_entries;
        std::vector<uint16_t> matching_lens;
        std::vector<uint64_t> matching_colors;
        std::vector<mem_t> mems;

        std::vector<uint16_t> scans;
        std::vector<uint16_t> fastforwards;
        std::vector<std::chrono::nanoseconds> costs;
};

#endif
