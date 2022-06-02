#include <iostream>
#include <iomanip>
#include <cassert>
#include "basic_block_stats.h"

#define assertm(exp, msg) assert(((void)msg, exp))

void
insert_bbcount(std::vector<size_t> &vec, size_t count)
{
    try {
        if (vec.size() <= count) {
            vec.resize(2 * count, 0);
        }
        vec[count] += 1;
    } catch (std::length_error &e) {
        printf("FAILED TO ADD DATA TO VECTOR.\n");
        printf("Current vector size: %lu\n", vec.size());
        printf("count              : %lu\n", count);
        std::cout << e.what() << std::endl;
    }
}

basic_block_stats_t::basic_block_stats_t(int block_size, const std::string &miss_file,
                                         bool warmup_enabled, bool is_coherent)
    : caching_device_stats_t(miss_file, block_size, warmup_enabled, is_coherent)
    , count_per_basic_block_instr_size_(block_size)
    , count_per_basic_block_byte_size_(block_size)
{
    // TODO: stats map is expecting long ints. We might not need this
    //    stats_map_.emplace(metric_name_t::BASIC_BLOCK_COUNTS,
    //                       count_per_basic_block_byte_size_);
}

void
basic_block_stats_t::access(const memref_t &memref, bool hit,
                            caching_device_block_t *cache_block)
{
    if (type_is_prefetch(memref.data.type)) {
        return; // we only care about basic blocks
    }

    bool is_instr_ = type_is_instr(memref.data.type);

    // count number of instructions in basic block
    if (is_instr_) {
        bb_size_instr++;
    }

    // calculate number of bytes in instruction based on addresses.
    if (type_is_instr_branch(memref.data.type)) {
        // TODO: Fetch address from branching information file.
        bb_end_addr = memref.data.addr;
        if (bb_end_addr < bb_start_addr) {
            printf("PANIC\n");
            printf("bbsize: %lu\n", bb_end_addr - bb_start_addr);
            printf("bb_size_instr: %lu\n", bb_size_instr);
            printf("IGNORING\n");

            bb_start_addr = 0;
            bb_size_instr = 0;
            return;
        }

        if (bb_start_addr != 0) {
            size_t bbsize = bb_end_addr - bb_start_addr;
            // assuming x86, instructions in basic blocks have monotoneously increasing
            // addresses
            insert_bbcount(count_per_basic_block_byte_size_, bbsize);
            insert_bbcount(count_per_basic_block_instr_size_, bb_size_instr);
        }
        bb_start_addr = 0;
        bb_size_instr = 0;
    } else if (is_instr_ && bb_start_addr == 0 && bb_end_addr != 0) {
        // We're the first address getting fetched after the first jump
        bb_start_addr = memref.data.addr;
    }
}

void
trim_vector(std::vector<size_t> &vec)
{
    int i;
    for (i = vec.size(); i >= 0 && vec[i] == 0; --i)
        ;
    vec.resize(i + 1);
}

void
basic_block_stats_t::print_stats(std::string prefix)
{
    trim_vector(count_per_basic_block_byte_size_);
    trim_vector(count_per_basic_block_instr_size_);

    std::cout << prefix << "bb byte size:" << std::endl;
    for (auto it = count_per_basic_block_byte_size_.begin();
         it != count_per_basic_block_byte_size_.end(); it++) {
        std::cout << prefix << "    " << it - count_per_basic_block_byte_size_.begin()
                  << ": " << *it << std::endl;
    }

    std::cout << prefix << "bb instr size:" << std::endl;
    for (auto it = count_per_basic_block_instr_size_.begin();
         it != count_per_basic_block_instr_size_.end(); it++) {
        std::cout << prefix << "    " << it - count_per_basic_block_instr_size_.begin()
                  << ": " << *it << std::endl;
    }
}

void
basic_block_stats_t::reset()
{
    num_hits_at_reset_ = num_hits_;
    num_misses_at_reset_ = num_misses_;
    num_child_hits_at_reset_ = num_child_hits_;
    num_hits_ = 0;
    num_misses_ = 0;
    num_compulsory_misses_ = 0;
    num_child_hits_ = 0;
    num_inclusive_invalidates_ = 0;
    num_coherence_invalidates_ = 0;
    count_per_basic_block_byte_size_.clear();
    count_per_basic_block_instr_size_.clear();
    bb_end_addr = __INT64_MAX__;
    bb_start_addr = 0;
    bb_size_instr = 0;
}

// void
// basic_block_stats_t::flush(const memref_t &memref)
//{
//
// }