#ifndef _I_CACHING_DEVICE_H_
#define _I_CACHING_DEVICE_H_ 1

#include <functional>
#include <unordered_map>
#include <vector>

#include "caching_device_block.h"
#include "vcl_caching_device_block.h"
#include "caching_device_stats.h"
#include "basic_block_stats.h"
#include "memref.h"
#include "prefetcher.h"

class snoop_filter_t;

class I_caching_device_t {
public:
    virtual ~I_caching_device_t()
    {
    }

    virtual bool
    init(int associativity, int block_size, int num_blocks, I_caching_device_t *parent,
         caching_device_stats_t *stats, prefetcher_t *prefetcher = nullptr,
         bool inclusive = false, bool coherent_cache = false, int id_ = -1,
         snoop_filter_t *snoop_filter_ = nullptr,
         const std::vector<I_caching_device_t *> &children = {}) = 0;

    virtual bool
    init(int associativity, std::vector<int> &way_sizes, int num_blocks,
         I_caching_device_t *parent, caching_device_stats_t *stats,
         prefetcher_t *prefetcher = nullptr, bool inclusive = false,
         bool coherent_cache = false, int id_ = -1,
         snoop_filter_t *snoop_filter_ = nullptr,
         const std::vector<I_caching_device_t *> &children = {}) = 0;

    virtual void
    request(const memref_t &memref) = 0;
    virtual void
    invalidate(addr_t tag, invalidation_type_t invalidation_type_) = 0;
    virtual bool
    contains_tag(addr_t tag) = 0;
    virtual void
    propagate_eviction(addr_t tag, const I_caching_device_t *requester) = 0;
    virtual void
    propagate_write(addr_t tag, const I_caching_device_t *requester) = 0;

    virtual caching_device_stats_t *
    get_stats() const = 0;
    virtual void
    set_stats(caching_device_stats_t *stats) = 0;
    virtual prefetcher_t *
    get_prefetcher() const = 0;
    virtual I_caching_device_t *
    get_parent() const = 0;
    virtual inline double
    get_loaded_fraction() const = 0;
    // Must be called prior to any call to request().
    virtual inline void
    set_hashtable_use(bool use_hashtable) = 0;
    virtual int
    get_block_index(const addr_t addr) const = 0;

protected:
    virtual void
    access_update(int block_idx, int way) {};
    virtual int
    replace_which_way(int block_idx)
    {
        return 0;
    };
    virtual int
    get_next_way_to_replace(const int block_idx) const
    {
        return 0;
    };
    virtual void
    record_access_stats(const memref_t &memref, bool hit,
                        caching_device_block_t *cache_block) {};

    virtual inline addr_t
    compute_tag(addr_t addr) const;
    virtual inline int
    compute_block_idx(addr_t tag) const;
    virtual inline caching_device_block_t &
    get_caching_device_block(int block_idx, int way) const;

    virtual inline void
    invalidate_caching_device_block(caching_device_block_t *block) = 0;

    virtual inline void
    update_tag(caching_device_block_t *block, int way, addr_t new_tag) = 0;

    // Returns the block (and its way) whose tag equals `tag`.
    // Returns <nullptr,0> if there is no such block.
    virtual std::pair<caching_device_block_t *, int>
    find_caching_device_block(addr_t tag);

    // a pure virtual function for subclasses to initialize their own block array
    virtual void
    init_blocks() = 0;
};

#endif