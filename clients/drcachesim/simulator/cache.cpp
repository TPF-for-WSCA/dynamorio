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

#include "cache.h"
#include "../common/utils.h"
#include <assert.h>

cache_t::cache_t(I_caching_device_t *self)
    : self_(self)
{
}

cache_t::~cache_t()
{
    delete self_;
    return;
}

bool
cache_t::init(int associativity, int line_size, int total_size,
              I_caching_device_t *parent, caching_device_stats_t *stats,
              prefetcher_t *prefetcher, bool inclusive, bool coherent_cache, int id,
              snoop_filter_t *snoop_filter,
              const std::vector<I_caching_device_t *> &children)
{
    // convert total_size to num_blocks to fit for self_->init
    int num_lines = total_size / line_size;

    return self_->init(associativity, line_size, num_lines, parent, stats, prefetcher,
                       inclusive, coherent_cache, id, snoop_filter, children);
}

void
cache_t::init_blocks()
{
    for (int i = 0; i < self_->num_blocks_; i++) {
        self_->blocks_[i] = new cache_line_t;
    }
}

void
cache_t::request(const memref_t &memref)
{
    self_->request(memref);
}

void
cache_t::flush(const memref_t &memref)
{
    addr_t tag = self_->compute_tag(memref.flush.addr);
    addr_t final_tag =
        self_->compute_tag(memref.flush.addr + memref.flush.size - 1 /*no overflow*/);
    self_->last_tag_ = TAG_INVALID;
    for (; tag <= final_tag; ++tag) {
        auto block_way = self_->find_caching_device_block(tag);
        if (block_way.first == nullptr)
            continue;
        self_->invalidate_caching_device_block(block_way.first);
    }
    // We flush parent_'s code cache here.
    // XXX: should L1 data cache be flushed when L1 instr cache is flushed?
    if (self_->parent_ != NULL)
        ((cache_t *)(self_->parent_->cache_))->flush(memref);
    if (self_->stats_ != NULL)
        ((cache_stats_t *)self_->stats_)->flush(memref);
}

void
cache_t::access_update(int block_idx, int way)
{
    self_->access_update(block_idx, way);
};

int
cache_t::replace_which_way(int block_idx)
{
    return self_->replace_which_way(block_idx);
};
int
cache_t::get_next_way_to_replace(const int block_idx) const
{
    return self_->get_next_way_to_replace(block_idx);
};
