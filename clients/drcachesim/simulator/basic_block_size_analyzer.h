#ifndef _BB_ANALYZER_H_
#define _BB_ANALYZER_H_

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "cache_simulator.h"
#include "caching_device_stats.h"
#include "../common/memref.h"

class basic_block_stats_t : public caching_device_stats_t {
public:
    explicit basic_block_stats_t(int block_size, const std::string &branch_file = "",
                                 bool warmup_enabled = false, bool is_coherent = false);
    void
    access(const memref_t &memref, bool hit,
           caching_device_block_t *cache_block) override;

    virtual void
    flush(const memref_t &memref);

    void
    reset() override;

protected:
    void
    print_counts(std::string prefix) override;

    std::vector<int_least64_t> count_per_basic_block_size_;
};

#endif