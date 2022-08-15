#include "vcl_caching_device.h"
#include "../common/utils.h"
#include <math.h>

vcl_caching_device_t::vcl_caching_device_t()
    : blocks_(NULL)
    , stats_(NULL)
    , prefetcher_(NULL)
    , tag2block(0, [](memref_t key) { return static_cast<unsigned long>(key.data.addr); })
{
}

vcl_caching_device_t::~vcl_caching_device_t()
{
    if (blocks_ == NULL) {
        return
    }
    for (auto const &block : blocks_) {
        delete block;
    }
    delete[] blocks_;
}

bool
vcl_caching_device_t::init(int associativity, int block_size, int num_blocks,
                           caching_device_t *parent, caching_device_stats_t *stats,
                           prefetcher_t *prefetcher, bool inclusive, bool coherent_cache,
                           int id_, snoop_filter_t *snoop_filter_,
                           const std::vector<caching_device_t *> &children)
{
    return false;
}

bool
vcl_caching_device_t::init(int associativity, std::vector<int> &way_sizes, int num_blocks,
                           I_caching_device_t *parent, caching_device_stats_t *stats,
                           prefetcher_t *prefetcher, bool inclusive, bool coherent_cache,
                           int id_, snoop_filter_t *snoop_filter,
                           const std::vector<I_caching_device_t *> &children)
{
    if (!IS_POWER_OF_2(associativity))
        return false;
    if (stats == NULL)
        return false;
    else if (!*stats)
        return false;

    associativity_ = associativity;
    block_sizes_ = way_sizes;
    block_size_ = std::max(way_sizes);
    num_blocks_ = num_blocks;
    loaded_blocks = 0;
    sets_in_cache_ = num_blocks_ / associativity;
    blocks_per_set_mask_ = sets_in_cache_ - 1;
    assoc_bits_ = std::ceil(
        std::log2(associativity_)); // TODO: Why do we need power of 2 number of ways?
    block_sizes_bits_.reserve(block_sizes_.size());
    block_size_bits_ = compute_log2(block_size_);

    for (const auto &block_size : block_sizes_) {
        block_sizes_bits_.push_back(std::ceil(std::log2(block_size)));
    }

    if (block_size_bits_ == -1 || assoc_bits_ == -1) {
        return false;
    }

    parent_ = parent;
    stats_ = stats;
    prefetcher_ = prefetcher;
    id_ = id;
    snoop_filter_ = snoop_filter;
    coherent_cache_ = coherent_cache;
    // TODO: check that - we need correctly sized blocks
    blocks_ = new vcl_caching_device_block_t *[num_blocks_];

    // initialize block sizes
    int way_count = 1;
    for (const auto &block_size : block_sizes_) {
        for (int i = (way_count - 1) * sets_in_cache_; i < way_count * sets_in_cache_;
             i++) {
            ((vcl_caching_device_block_t *)blocks_[i])->size_ = block_size;
        }
    }
    init_blocks();
    last_tag_ = TAG_INVALID;
    inclusive_ = inclusive;
    children_ = children;

    return true;
}