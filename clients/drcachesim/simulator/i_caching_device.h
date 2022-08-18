#ifndef _I_CACHING_DEVICE_H_
#define _I_CACHING_DEVICE_H_ 1

#include <functional>
#include <unordered_map>
#include <vector>

#include "caching_device_block.h"
#include "basic_block_stats.h"
#include "memref.h"
#include "prefetcher.h"

class snoop_filter_t;

class I_caching_device_t {
public:
    I_caching_device_t()
        : blocks_(NULL)
        , stats_(NULL)
        , prefetcher_(NULL)
        , cache_(NULL)
        // The tag being hashed is already right-shifted to the cache line and
        // an identity hash is plenty good enough and nice and fast.
        // We set the size and load factor only if being used, in set_hashtable_use().
        , tag2block(0, [](addr_t key) { return static_cast<unsigned long>(key); })
    {
    }

    virtual ~I_caching_device_t()
    {
    }

    virtual bool vcl_enabled() = 0;

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

    virtual addr_t
    compute_tag(addr_t addr) const = 0;
    virtual int
    compute_block_idx(addr_t tag) const = 0;
    virtual caching_device_block_t &
    get_caching_device_block(int block_idx, int way) const = 0;

    virtual inline void
    invalidate_caching_device_block(caching_device_block_t *block) = 0;

    virtual inline void
    update_tag(caching_device_block_t *block, int way, addr_t new_tag) = 0;

    // Returns the block (and its way) whose tag equals `tag`.
    // Returns <nullptr,0> if there is no such block.
    virtual std::pair<caching_device_block_t *, int>
    find_caching_device_block(addr_t tag) = 0;

    // a pure virtual function for subclasses to initialize their own block array
    virtual void
    init_blocks() = 0;

    int num_blocks_;
    bool coherent_cache_;
    // This is an index into snoop filter's array of caches.
    int id_;

    // Current valid blocks in the cache
    int loaded_blocks_;

    // Pointers to the caching device's parent and children devices.
    I_caching_device_t *parent_;
    std::vector<I_caching_device_t *> children_;

    snoop_filter_t *snoop_filter_;

    // If true, this device is inclusive of its children.
    bool inclusive_;

    // This should be an array of caching_device_block_t pointers, otherwise
    // an extended block class which has its own member variables cannot be indexed
    // correctly by base class pointers.
    caching_device_block_t **blocks_;
    int blocks_per_set_;
    int sets_in_cache_;
    int associativity_;
    int block_size_;
    // len of block_sizes_ needs to be equal associativity
    std::vector<int> block_sizes_;
    caching_device_stats_t *stats_;
    // Optimization fields for fast bit operations
    int blocks_per_set_mask_;
    int assoc_bits_;
    std::vector<int> block_sizes_bits_;
    int block_size_bits_;

    prefetcher_t *prefetcher_;

    // Optimization: remember last tag
    addr_t last_tag_;
    int last_way_;
    int last_block_idx_;

    void *cache_;

    // Optimization: keep a hashtable for quick lookup of {block,way}
    // given a tag, if using a large cache hierarchy where serial
    // walks over the associativity end up as bottlenecks.
    // We can't easily remove the blocks_ array and replace with just
    // the hashtable as replace_which_way(), etc. want quick access to
    // every way for a given line index.
    std::unordered_map<addr_t, std::pair<caching_device_block_t *, int>,
                       std::function<unsigned long(addr_t)>>
        tag2block;
    bool use_tag2block_table_ = false;
};

#endif