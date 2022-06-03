#ifndef _BB_ANALYZER_H_
#define _BB_ANALYZER_H_

#define MAX_X86_INSTR_SIZE 13

#include <cstdint>
#include <string>
#include <unordered_map>
#include <set>
#include <vector>

#include "caching_device_stats.h"

struct BasicBlock {
    addr_t starting_addr;
    addr_t end_addr;
    size_t instr_size;
    size_t byte_size;
};

void
print_type(trace_type_t type);

bool
operator==(const BasicBlock &lhs, const BasicBlock &rhs);

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
    std::unordered_map<BasicBlock, size_t> basic_blocks_hit_count;

private:
    void
    record_block(const BasicBlock &bb);

    size_t max_instr_size = 0;
    BasicBlock current_block;
};

#endif