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

bool
vcl_caching_device_t::vcl_enabled()
{
    return true;
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
    num_sets_ = num_blocks_ / associativity;
    set_idx_mask_ = sets_in_cache_ - 1;
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
vcl_caching_device_t::request(_memref_t const &memref_in)
{
    memref_t memref;
    addr_t final_addr = memref_in.data.addr + memref_in.data.size - 1;
    addr_t final_tag = compute_tag(final_addr);
    addr_t tag = compute_tag(memref_in.data.addr);

    memref = memref_in;
    for (; tag <= final_tag; ++tag) {
        int way = associativity_;
        int block_idx = compute_block_idx(tag);

        bool missed = false;
        if (tag + 1 <= final_tag) {
            memref.data.size = ((tag + 1) << block_size_bits_) -
                memref.data.addr; // is this good enough?
        }

        auto block_way = find_caching_device_block(tag);
        if (block_way.first != nullptr) {
            // found
        } else {
            // miss
            way = replace_which_way(block_idx, )
        }
    }
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
vcl_caching_device_t::replace_which_way(int block_idx, int size)
{
    int first_viable_way;
    for (int i = 0; i < block_sizes_.size(); ++i) {
        if (block_sizes_[i] > size) {
            first_viable_way = i;
            break;
        }
    }

    int min_way = 0;
    int min_counter = 0;
    // We assume that the ways are sorted by size, such that we can assume that the
    // smalles one fitting comes first.
    std::vector<std::pair<int, int>>
        candidate_sizes_and_idx; // We use it as <size, idx>. All INVALID lines are
                                 // candidates + LRU candidate
    for (int way = first_viable_way; way < associativity_; ++way) {
        if (get_caching_device_block(block_idx, way).tag_ == TAG_INVALID) {
            candidate_sizes_and_idx.push_back(std::pair(block_sizes_[way], way));
        }
        if (way == first_viable_way ||
            get_caching_device_block(block_idx, way).counter_ < min_counter) {
            min_counter = get_caching_device_block(block_idx, way).counter_;
            min_way = way;
        }
    }
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
