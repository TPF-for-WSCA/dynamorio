#include "vcl_caching_device.h"
#include "../common/utils.h"
#include <math.h>

vcl_caching_device_t::vcl_caching_device_t()
    : I_caching_device_t()
{
}

vcl_caching_device_t::~vcl_caching_device_t()
{
    if (blocks_ == NULL) {
        return;
    }
    for (int i = 0; i < num_blocks_; i++) {
        delete blocks_[i];
    }
    delete[] blocks_;
}

void
vcl_caching_device_t::init_blocks()
{
    for (int i = 0; i < num_blocks_; i++) {
        blocks_[i] = new vcl_caching_device_block_t;
    }
}

bool
vcl_caching_device_t::init(int associativity, int block_size, int num_blocks,
                           I_caching_device_t *parent, caching_device_stats_t *stats,
                           prefetcher_t *prefetcher, bool inclusive, bool coherent_cache,
                           int id, snoop_filter_t *snoop_filter_,
                           const std::vector<I_caching_device_t *> &children)
{
    return false;
}

bool
vcl_caching_device_t::init(int associativity, std::vector<int> &way_sizes, int num_blocks,
                           I_caching_device_t *parent, caching_device_stats_t *stats,
                           prefetcher_t *prefetcher, bool inclusive, bool coherent_cache,
                           int id, snoop_filter_t *snoop_filter,
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
    block_size_ = *(std::max_element(way_sizes.begin(), way_sizes.end()));
    num_blocks_ = num_blocks;
    loaded_blocks_ = 0;
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
    blocks_ = (caching_device_block_t **)new vcl_caching_device_block_t *[num_blocks_];

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

void
vcl_caching_device_t::request(_memref_t const &memref)
{
}

void
vcl_caching_device_t::invalidate(addr_t tag, invalidation_type_t invalidation_type_)
{
}

bool
vcl_caching_device_t::contains_tag(addr_t tag)
{
    return false;
}

void
vcl_caching_device_t::propagate_eviction(addr_t tag, const I_caching_device_t *requester)
{
}

void
vcl_caching_device_t::propagate_write(addr_t tag, const I_caching_device_t *requester)
{
}

void
vcl_caching_device_t::access_update(int block_idx, int way)
{
}

int
vcl_caching_device_t::replace_which_way(int block_idx)
{
    return -1;
}

int
vcl_caching_device_t::get_next_way_to_replace(const int block_idx) const
{
    return -1;
}

// TODO: Type double that function for the vcl blocks r unify block interface
void
vcl_caching_device_t::record_access_stats(const _memref_t &memref, bool hit,
                                          caching_device_block_t *cache_block)
{
}

std::pair<caching_device_block_t *, int>
vcl_caching_device_t::find_caching_device_block(addr_t tag)
{
    return std::make_pair(new caching_device_block_t, -1);
}
