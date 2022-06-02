#ifndef _BB_ANALYZER_H_
#define _BB_ANALYZER_H_

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

// #include "cache_simulator.h"
#include "caching_device_stats.h"
// #include "../common/memref.h"

void
insert_bbcount(std::vector<int_least64_t> &vec, int_least64_t count);

class basic_block_stats_t : public caching_device_stats_t {
public:
    explicit basic_block_stats_t(int block_size, const std::string &miss_file = "",
                                 bool warmup_enabled = false, bool is_coherent = false);
    void
    access(const memref_t &memref, bool hit,
           caching_device_block_t *cache_block) override;

    //    virtual void
    //    flush(const memref_t &memref);

    void
    print_stats(std::string prefix) override;

    void
    reset() override;

protected:
    // void
    // print_counts(std::string prefix) override;

    std::vector<size_t> count_per_basic_block_instr_size_;
    std::vector<size_t> count_per_basic_block_byte_size_;

private:
    size_t bb_start_addr = 0;
    size_t bb_end_addr = __UINT64_MAX__;
    size_t bb_size_instr = 0;
};

#endif