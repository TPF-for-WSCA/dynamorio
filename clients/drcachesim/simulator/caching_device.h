/* **********************************************************
 * Copyright (c) 2015-2021 Google, Inc.  All rights reserved.
 * **********************************************************/

/*
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * * Neither the name of Google, Inc. nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without
 *   specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL VMWARE, INC. OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 */

/* caching_device: represents a hardware caching device.
 */

#ifndef _CACHING_DEVICE_H_
#define _CACHING_DEVICE_H_ 1

#include "i_caching_device.h"

// Statistics collection is abstracted out into the caching_device_stats_t class.

// Different replacement policies are expected to be implemented by
// subclassing caching_device_t.

// We assume we're only invoked from a single thread of control and do
// not need to synchronize data access.

class caching_device_t : public I_caching_device_t {
public:
    caching_device_t();
    bool
    init(int associativity, int block_size, int num_blocks, I_caching_device_t *parent,
         caching_device_stats_t *stats, prefetcher_t *prefetcher = nullptr,
         bool inclusive = false, bool coherent_cache = false, int id_ = -1,
         snoop_filter_t *snoop_filter_ = nullptr,
         const std::vector<I_caching_device_t *> &children = {}) override;

    bool
    init(int associativity, std::vector<int> &way_sizes, int num_blocks,
         I_caching_device_t *parent, caching_device_stats_t *stats,
         prefetcher_t *prefetcher = nullptr, bool inclusive = false,
         bool coherent_cache = false, int id_ = -1,
         snoop_filter_t *snoop_filter_ = nullptr,
         const std::vector<I_caching_device_t *> &children = {}) override;

    virtual bool
    vcl_enabled() override;

    virtual ~caching_device_t();
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
        stats_ = stats;
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

    virtual int
    replace_which_way(int block_idx) override;
    virtual int
    get_next_way_to_replace(const int block_idx) const override;
    virtual void
    record_access_stats(const memref_t &memref, bool hit,
                        caching_device_block_t *cache_block) override;

    virtual inline addr_t
    compute_tag(addr_t addr) const override
    {
        return addr >> block_size_bits_;
    }
    virtual inline int
    compute_block_idx(addr_t tag) const override
    {
        return (tag & set_idx_mask_) << assoc_bits_;
    }
    virtual inline caching_device_block_t &
    get_caching_device_block(int block_idx, int way) const override
    {
        return *(blocks_[block_idx + way]);
    }

    virtual inline void
    invalidate_caching_device_block(caching_device_block_t *block) override
    {
        if (use_tag2block_table_)
            tag2block.erase(block->tag_);
        block->tag_ = TAG_INVALID;
        // Xref cache_block_t constructor about why we set counter to 0.
        block->counter_ = 0;
    }

    virtual inline void
    update_tag(caching_device_block_t *block, int way, addr_t new_tag) override
    {
        if (use_tag2block_table_) {
            if (block->tag_ != TAG_INVALID)
                tag2block.erase(block->tag_);
            tag2block[new_tag] = std::make_pair(block, way);
        }
        block->tag_ = new_tag;
    }

    // Returns the block (and its way) whose tag equals `tag`.
    // Returns <nullptr,0> if there is no such block.
    virtual std::pair<caching_device_block_t *, int>
    find_caching_device_block(addr_t tag) override;

    // a pure virtual function for subclasses to initialize their own block array
    virtual void
    init_blocks() override;
};

#endif /* _CACHING_DEVICE_H_ */
