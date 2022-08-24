#ifndef _VCL_CACHING_DEVICE_H_
#define _VCL_CACHING_DEVICE_H_

// #define _VCL_SETS_EXTENSION
#define _VCL_WAY_EXTENSION

#include "i_caching_device.h"
#include <set>
#include <fstream>
#include "vcl_caching_device_block.h"
struct AddrBlockCmp {
    bool
    operator()(const std::pair<addr_t, addr_t> &lhs,
               const std::pair<addr_t, addr_t> &rhs) const
    {
        bool r = lhs.first < rhs.first;
        if (lhs.first == rhs.first)
            r = lhs.second < rhs.second;
        return r;
    }
};
class vcl_caching_device_t : public I_caching_device_t {

public:
    vcl_caching_device_t(std::string perfect_fetch_file);
    using I_caching_device_t::init;

    virtual bool
    vcl_enabled() override;

    virtual bool
    init(int associativity, std::vector<int> &way_sizes, int num_blocks,
         I_caching_device_t *parent, caching_device_stats_t *stats,
         prefetcher_t *prefetcher = nullptr, bool inclusive = false,
         bool coherent_cache = false, int id_ = -1,
         snoop_filter_t *snoop_filter_ = nullptr,
         const std::vector<I_caching_device_t *> &children = {}) override;
    virtual bool
    init(int associativity, int block_size, int num_blocks, I_caching_device_t *parent,
         caching_device_stats_t *stats, prefetcher_t *prefetcher = nullptr,
         bool inclusive = false, bool coherent_cache = false, int id_ = -1,
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
    I_caching_device_t *
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
    virtual int
    replace_which_way(int block_idx) override;
    int
    replace_which_way(int block_idx, int size);
    virtual int
    get_next_way_to_replace(const int block_idx) const override;
    virtual void
    record_access_stats(const memref_t &memref, bool hit,
                        caching_device_block_t *cache_block) override;

    // TODO: How to fix tag calculation if we dont yet know size of the way
    inline addr_t
    compute_tag(addr_t addr) const override
    {
        return addr >> block_size_bits_;
    }
    inline int
    compute_block_idx(addr_t tag) const override
    {
#ifdef _VCL_WAY_EXTENSION
        return (tag & set_idx_mask_) * block_sizes_.size();
#else // _VCL_SET_EXTENSION
        return -1; // todo
#endif
    }
    inline caching_device_block_t &
    get_caching_device_block(int block_idx, int way) const override
    {
        return *(blocks_[block_idx + way]);
    }

    inline void
    invalidate_caching_device_block(caching_device_block_t *block) override
    {
        if (use_tag2block_table_) {
            tag2block.erase(block->tag_);
        }

        block->tag_ = TAG_INVALID;
        block->counter_ = 0;
    }

    inline int
    get_smallest_possible_way(int size)
    {
        int first_viable_way;
        for (size_t i = 0; i < block_sizes_.size(); ++i) {
            if (block_sizes_[i] >= size) {
                first_viable_way = i;
                break;
            }
        }
        return first_viable_way;
    }

    inline void
    update_tag(caching_device_block_t *block, int way, addr_t new_tag) override
    {
        if (use_tag2block_table_) {
            if (block->tag_ != TAG_INVALID) {
                tag2block.erase(block->tag_);
            }
            tag2block[new_tag] = std::make_pair(block, way);
        }
        block->tag_ = new_tag;
        block->validity_ = true;
    }

    // Returns the block (and its way) whose tag equals `tag`.
    // Returns <nullptr,0> if there is no such block.
    virtual std::pair<caching_device_block_t *, int>
    find_caching_device_block(addr_t tag) override;

    virtual std::vector<std::pair<caching_device_block_t *, int>>
    find_all_caching_device_blocks(addr_t addr, bool only_one = false);

    // a pure virtual function for subclasses to initialize their own block array
    virtual void
    init_blocks() override;

private:
    std::pair<int, int>
    start_and_end_oracle(addr_t address);
    bool
    read_n_oracle_lines(size_t n);
    std::map<addr_t, std::set<std::pair<addr_t, addr_t>, AddrBlockCmp>>
        base_address_to_blocks_mapping;
    std::ifstream perfect_loading_decision_fh_;
};
#endif /* _VCL_CACHING_DEVICE_H_ */