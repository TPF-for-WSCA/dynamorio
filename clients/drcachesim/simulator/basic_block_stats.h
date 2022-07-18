#ifndef _BB_ANALYZER_H_
#define _BB_ANALYZER_H_

#define MAX_X86_INSTR_SIZE 15

#include <cstdint>
#include <string>
#include <unordered_map>
#include <set>
#include <vector>
#include <deque>

#include "caching_device_stats.h"
#include "CTikz.hpp"
#include "CException.hpp"

struct BasicBlock {
    addr_t starting_addr;
    addr_t end_addr;
    size_t instr_size;
    size_t byte_size;
    mutable size_t hits = 1;
    mutable bool miss = true;
};

/**
 * @brief hitbytes describes a pair [counter, accessed bytes],
 * denoting how often the number of bytes were accessed
 *
 */
typedef std::pair<size_t, uint8_t> hitbytes;

void
print_type(trace_type_t type);

bool
operator==(const BasicBlock &lhs, const BasicBlock &rhs);

struct BBLessThanBasic {
    inline bool
    operator()(const BasicBlock &lhs, const BasicBlock &rhs) const
    {
        return lhs.end_addr <
            rhs.end_addr; // assuming that we always exit a bb at a branch
    }
};

namespace std {
template <> struct hash<BasicBlock> {
    size_t
    operator()(const BasicBlock &b) const
    {
        using std::hash;
        using std::size_t;
        // Hack to not cancel out if start and end address are the same (using xor)
        return (hash<addr_t>()(b.starting_addr) ^ hash<addr_t>()(b.end_addr << 1) >> 1);
    }
};
}

void
insert_bbcount(std::vector<int_least64_t> &vec, int_least64_t count);

class basic_block_stats_t : public caching_device_stats_t {
public:
    explicit basic_block_stats_t(int block_size, const std::string &miss_file = "",
                                 const std::string &output_dir = "",
                                 bool warmup_enabled = false, bool is_coherent = false);
    void
    access(const memref_t &memref, bool hit,
           caching_device_block_t *cache_block) override;

    void
    print_stats(std::string prefix) override;

    void
    reset() override;

protected:
    // NOTE: Cacheline presence means for every time this address range has been loaded
    // into the cache we start a new counter
    std::unordered_map<addr_t, std::vector<std::vector<uint64_t>>>
        bytes_accessed_per_presence_per_cacheline;
    std::vector<size_t> count_per_basic_block_instr_size_;
    std::vector<size_t> count_per_basic_block_byte_size_;
    std::vector<size_t> basic_block_size_history;
    std::unordered_map<BasicBlock, size_t> basic_blocks_hit_count;
    std::unordered_map<addr_t, std::vector<BasicBlock>> number_of_bytes_accessed;

private:
    void
    record_block(const BasicBlock &bb);
    bool
    handle_instr(const memref_t &memeref, bool hit);
    void
    handle_branch(const memref_t &memref);
    void
    handle_interrupt(const memref_t &memref, bool hit);

    void
    reset_current_cacheline_block();

    void
    reset_current_cacheline_block(const memref_t &memref, bool hit);

    void
    record_cacheblock(const caching_device_block_t &t);

    void
    print_last_n_cacheblocks(int n);

    void
    record_memref(memref_t m, bool hit);

    void
    print_last_n_memrefs(int n);

    std::vector<uint64_t>
    aggregate_byte_accesses_in_cacheline_presence(
        const std::vector<uint64_t> &access_masks);

    std::vector<uint64_t>
    aggregate_byte_accesses_in_cacheline(
        const std::vector<std::vector<uint64_t>> &cacheline_presences);

    void
    track_cacheline_access(const memref_t &memref);

    uint64_t
    bytes_accessed_by_block(const addr_t &cacheline_base, BasicBlock &block);

    std::pair<uint8_t, std::vector<hitbytes>>
    bytes_accessed(const addr_t &cacheline_base,
                   std::vector<BasicBlock> &blocks_contained);

    void
    print_bytes_accessed();

    size_t
    create_histogram_of_cachelineaccesses(std::vector<uint64_t> &histogram,
                                          std::vector<hitbytes> &accesses);

    size_t max_cacheline_bb = 0;
    size_t handled_instructions = 0;
    size_t num_interrupts = 0;
    const std::string output_dir;
    size_t max_instr_size = 0;
    BasicBlock current_block;
    BasicBlock current_block_cacheline_constrained;
    static const size_t cache_line_address_mask = 0xFFFFFFFFFFFFFFC0;
    std::deque<BasicBlock> basic_block_buffer;
    std::deque<caching_device_block_t> cache_block_buffer;
    // TODO: Implement counter function for this (mapping of current access cacheline size to set of distinct addresses)
    std::vector<std::set<addr_t>> distinct_cachelines_by_size;
};

#endif