#include "vcl_cache_stats.h"

vcl_cache_stats_t::vcl_cache_stats_t(int block_size, const std::string &miss_file,
                                     const std::string &output_dir, bool warmup_enbled,
                                     bool is_coherent)
    : basic_block_stats_t(block_size, miss_file, output_dir, warmup_enbled, is_coherent)
{
}

void
vcl_cache_stats_t::access(const memref_t &memref, bool hit,
                          caching_device_block_t *cache_block)
{
}

void
vcl_cache_stats_t::print_stats(std::string prefixt)
{
}

void
vcl_cache_stats_t::reset()
{
}