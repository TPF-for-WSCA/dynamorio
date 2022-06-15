#ifndef _BB_ANALYZER_H_
#define _BB_ANALYZER_H_

#define MAX_X86_INSTR_SIZE 15

#include <cstdint>
#include <string>
#include <unordered_map>
#include <set>
#include <vector>

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
    std::vector<size_t> count_per_basic_block_instr_size_;
    std::vector<size_t> count_per_basic_block_byte_size_;
    std::vector<size_t> basic_block_size_history;
    std::unordered_map<BasicBlock, size_t> basic_block_cache_aligned_hit_count;
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
    track_cacheline_access(const memref_t &memref);

    uint64_t
    bytes_accessed_by_block(const addr_t &cacheline_base, BasicBlock &block);

    std::pair<uint8_t, std::vector<uint8_t>>
    bytes_accessed(const addr_t &cacheline_base,
                   std::vector<BasicBlock> &blocks_contained);

    void
    print_bytes_accessed();

    void
    create_histogram_of_cachelineaccesses(std::vector<uint64_t> &histogram,
                                          std::vector<uint8_t> &accesses);

    size_t max_cacheline_bb = 0;
    const std::string output_dir;
    size_t max_instr_size = 0;
    BasicBlock current_block;
    size_t max_memory_consumption = 0;
    static const size_t cache_block_address_mask = 0xFFFFFFFFFFFFFFC0;
};

#endif