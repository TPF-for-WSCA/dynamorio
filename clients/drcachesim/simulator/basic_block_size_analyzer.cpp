#include <iostream>
#include <iomanip>
#include "basic_block_stats.h"

void
insert_bbcount(std::vector<int_least64_t> &vec, int_least64_t count)
{
    if (vec.size() <= count) {
        vec.resize(2 * count, 0);
    }
    vec[count]++;
}

basic_block_stats_t::basic_block_stats_t(int block_size, const std::string &miss_file,
                                         bool warmup_enabled, bool is_coherent)
    : caching_device_stats_t(miss_file, block_size, warmup_enabled, is_coherent)
    , count_per_basic_block_instr_size_(block_size)
    , count_per_basic_block_byte_size_(block_size)
{
    stats_map_.emplace(metric_name_t::BASIC_BLOCK_COUNTS, count_per_basic_block_size);
}

void
basic_block_stats_t::access(const memref_t &memref, bool hit,
                            caching_device_block_t *cache_block)
{
    if (type_is_instr_branch(memref.data.type)) {
        // TODO: Fetch address from branching information file.
        bb_end_addr = memref.data.addr;
        if (bb_start_addr != 0) {
            insert_bbcount(count_per_basic_block_byte_size_, bb_end_addr - bb_start_addr);
        }
    } else if (type_is_instr(memref.data.type) && bb_start_addr == 0 &&
               bb_end_addr != 0) {
        // We're the first address getting fetched after the first jump
        bb_start_addr = memref.data.addr;
    }
}