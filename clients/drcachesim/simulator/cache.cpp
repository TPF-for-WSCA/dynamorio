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

bool
cache_t::init(int associativity, int line_size, int total_size, caching_device_t *parent,
              caching_device_stats_t *stats, prefetcher_t *prefetcher, bool inclusive,
              bool coherent_cache, int id, snoop_filter_t *snoop_filter,
              const std::vector<caching_device_t *> &children)
{
    // convert total_size to num_blocks to fit for caching_device_t::init
    int num_lines = total_size / line_size;

    return caching_device_t::init(associativity, line_size, num_lines,
                                  (caching_device_t *)parent, stats, prefetcher,
                                  inclusive, coherent_cache, id, snoop_filter,
                                  static_cast<std::vector<caching_device_t *>>(children));
}

void
cache_t::init_blocks()
{
    for (int i = 0; i < caching_device_t::num_blocks_; i++) {
        caching_device_t::blocks_[i] = new cache_line_t;
    }
}

void
cache_t::request(const memref_t &memref)
{
    caching_device_t::request(memref);
}

void
cache_t::flush(const memref_t &memref)
{
    addr_t tag = caching_device_t::compute_tag(memref.flush.addr);
    addr_t final_tag = caching_device_t::compute_tag(
        memref.flush.addr + memref.flush.size - 1 /*no overflow*/);
    caching_device_t::last_tag_ = TAG_INVALID;
    for (; tag <= final_tag; ++tag) {
        auto block_way = caching_device_t::find_caching_device_block(tag);
        if (block_way.first == nullptr)
            continue;
        caching_device_t::invalidate_caching_device_block(block_way.first);
    }
    // We flush parent_'s code cache here.
    // XXX: should L1 data cache be flushed when L1 instr cache is flushed?
    if (caching_device_t::parent_ != NULL)
        ((cache_t *)caching_device_t::parent_)->flush(memref);
    if (caching_device_t::stats_ != NULL)
        ((cache_stats_t *)caching_device_t::stats_)->flush(memref);
}
