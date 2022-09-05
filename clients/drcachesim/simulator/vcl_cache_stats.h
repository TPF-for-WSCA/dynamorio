#ifndef _VCL_ANALYZER_H_
#define _VCL_ANALYZER_H_

#define MAX_X86_INSTR_SIZE 15

#include <cstdint>
#include <string>
#include <unordered_map>
#include <set>
#include <vector>
#include <deque>

#include "basic_block_stats.h"
#include "vcl_caching_device_block.h"
#include "CTikz.hpp"
#include "CException.hpp"

class vcl_cache_stats_t : public basic_block_stats_t {
public:
    explicit vcl_cache_stats_t(int block_size, const std::string &miss_file = "",
                               const std::string &output_dir = "",
                               bool warmup_enbled = false, bool is_coherent = false);
    void
    access(const memref_t &memref, bool hit,
           caching_device_block_t *cache_block) override;

    void
    print_stats(std::string prefixt) override;

    void
    reset() override;

protected:
private:
};

#endif // _VCL_ANALYZER_H_