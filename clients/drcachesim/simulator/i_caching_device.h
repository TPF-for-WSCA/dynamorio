#ifndef _I_CACHING_DEVICE_H_
#define _I_CACHING_DEVICE_H_ 1

#include <functional>
#include <unordered_map>
#include <vector>

#include "caching_device_block.h"
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
    init(int associativity, int block_size, int num_blocks, caching_device_t *parent,
         caching_device_stats_t *stats, prefetcher_t *prefetcher = nullptr,
         bool inclusive = false, bool coherent_cache = false, int id_ = -1,
         snoop_filter_t *snoop_filter_ = nullptr,
         const std::vector<caching_device_t *> &children = {}) = 0;

    virtual void
    request(const memref_t &memref) = 0;
    virtual void
    invalidate(addr_t tag, invalidation_type_t invalidation_type_) = 0;
    virtual bool
    contains_tag(addr_t tag) = 0;
    virtual void
    propagate_eviction(addr_t tag, const caching_device_t *requester) = 0;
    virtual void
    propagate_write(addr_t tag, const caching_device_t *requester) = 0;

    virtual caching_device_stats_t *
    get_stats() const = 0;
    virtual void
    set_stats(caching_device_stats_t *stats) = 0;
    virtual prefetcher_t *
    get_prefetcher() const = 0;
    virtual caching_device_t *
    get_parent() const = 0;
    virtual inline double
    get_loaded_fraction() const = 0;
    // Must be called prior to any call to request().
    virtual inline void
    set_hashtable_use(bool use_hashtable) = 0;
    virtual int
    get_block_index(const addr_t addr) const = 0;
};

#endif