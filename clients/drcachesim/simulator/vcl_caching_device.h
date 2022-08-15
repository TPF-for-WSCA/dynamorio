#ifndef _VCL_CACHING_DEVICE_H_
#define _VCL_CACHING_DEVICE_H_

#include "i_caching_device.h"

class vcl_caching_device_t : public I_caching_device_t {

public:
    vcl_caching_device_t();
    using I_caching_device_t::init;
    virtual bool
    init(int associativity, std::vector<int> &way_sizes, int num_blocks,
         I_caching_device_t *parent, caching_device_stats_t *stats,
         prefetcher_t *prefetcher = nullptr, bool inclusive = false,
         bool coherent_cache = false, int id_ = -1,
         snoop_filter_t *snoop_filter_ = nullptr,
         const std::vector<I_caching_device_t *> &children = {}) override;
    virtual ~vcl_caching_device_t() override;
    virtual void
    request(const memref_t &memref) override;
    virtual void
    invalidate(addr_t tag, invalidation_type_t invalidation_type_) override;
    bool
    contains_tag(addr_t tag) override;
    void
    propagate_eviction(addr_t tag, const I_caching_device_t *requester) override;
    void
    propagate_write(addr_t tag, const I_caching_device_t *requester) override;

    caching_device_stats_t *
    get_stats() const override
    {
        return stats_;
    }

    void
    set_stats(caching_device_stats_t *stats) override
    {
        this->stats_ = stats;
    }

    prefetcher_t *
    get_prefetcher() const override
    {
        return prefetcher_;
    }
    caching_device_t *
    get_parent() const override
    {
        return parent_;
    }
    inline double
    get_loaded_fraction() const override
    {
        return double(loaded_blocks_) / num_blocks_;
    }
    // Must be called prior to any call to request().
    virtual inline void
    set_hashtable_use(bool use_hashtable) override
    {
        if (!use_tag2block_table_ && use_hashtable) {
            // Resizing from an initial small table causes noticeable overhead, so we
            // start with a relatively large table.
            tag2block.reserve(1 << 16);
            // Even with the large initial size, for large caches we want to keep the
            // load factor small.
            tag2block.max_load_factor(0.5);
        }
        use_tag2block_table_ = use_hashtable;
    }
    int
    get_block_index(const addr_t addr) const override
    {
        addr_t tag = compute_tag(addr);
        int block_idx = compute_block_idx(tag);
        return block_idx;
    }

protected:
    virtual void
    access_update(int block_idx, int way) override;
    virtual int
    replace_which_way(int block_idx) override;
    virtual int
    get_next_way_to_replace(const int block_idx) const override;
    virtual void
    record_access_stats(const memref_t &memref, bool hit,
                        caching_device_block_t *cache_block) override;

    inline addr_t
    compute_tag(addr_t addr) const override
    {
        // look up next branch
        return 0;
    }
    inline int
    compute_block_idx(addr_t tag) const override
    {
        // TODO:;
        return 0;
    }
    inline caching_device_block_t &
    get_caching_device_block(int block_idx, int way) const override
    {
        // TODO;
        return *this->blocks_[0];
    }

    inline void
    invalidate_caching_device_block(caching_device_block_t *block) override
    {
        // TODO
    }

    inline void
    update_tag(caching_device_block_t *block, int way, addr_t new_tag) override
    {
        // TODO
    }

    // Returns the block (and its way) whose tag equals `tag`.
    // Returns <nullptr,0> if there is no such block.
    virtual std::pair<caching_device_block_t *, int>
    find_caching_device_block(addr_t tag) override;

    // a pure virtual function for subclasses to initialize their own block array
    virtual void
    init_blocks() override;

    int associativity_;
    // len of block_sizes_ needs to be equal associativity
    std::vector<int> block_sizes_;
    int num_blocks_;
    bool coherent_cache_;
    // This is an index into snoop filter's array of caches.
    int id_;

    // Current valid blocks in the cache
    int loaded_blocks_;

    // Pointers to the caching device's parent and children devices.
    caching_device_t *parent_;
    std::vector<caching_device_t *> children_;

    snoop_filter_t *snoop_filter_;

    // If true, this device is inclusive of its children.
    bool inclusive_;

    // This should be an array of caching_device_block_t pointers, otherwise
    // an extended block class which has its own member variables cannot be indexed
    // correctly by base class pointers.
    caching_device_block_t **blocks_;
    int sets_in_cache_;
    // Optimization fields for fast bit operations
    int blocks_per_set_mask_;
    int assoc_bits_;
    std::vector<int> block_sizes_bits_;
    int block_size_bits_;

    caching_device_stats_t *stats_;
    prefetcher_t *prefetcher_;

    // Optimization: remember last tag
    addr_t last_tag_;
    int last_way_;
    int last_block_idx_;
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
#endif /* _VCL_CACHING_DEVICE_H_ */