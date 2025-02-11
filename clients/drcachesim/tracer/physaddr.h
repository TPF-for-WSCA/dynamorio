/* **********************************************************
 * Copyright (c) 2015-2022 Google, Inc.  All rights reserved.
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

/* physaddr: translates virtual addresses to physical addresses.
 */

#ifndef _PHYSADDR_H_
#define _PHYSADDR_H_ 1

#include <fstream>
#include <unordered_map>
#include "dr_api.h"
#include "hashtable.h"
#include "../common/trace_entry.h"

// This class is not thread-safe: the caller should create a separate instance
// per thread.
class physaddr_t {
public:
    physaddr_t();
    ~physaddr_t();
    bool
    init();
    addr_t
    virtual2physical(addr_t virt);

private:
#ifdef LINUX
    addr_t last_vpage_;
    addr_t last_ppage_;
    // TODO i#4014: An app with thousands of threads might hit open file limits,
    // and even a hundred threads will use up DR's private FD limit and push
    // other files into potential app conflicts.
    // Sharing the descriptor would require locks, however.  Evaluating
    // how best to do that (maybe the caching will reduce the contention enough)
    // is future work.
    int fd_;
    // We would use std::unordered_map, but that is not compatible with
    // statically linking drmemtrace into an app.
    hashtable_t v2p_;
    // With hashtable_t nullptr is how non-existence is shown, so we store
    // an actual 0 address (can happen for physical) as this sentinel.
    static constexpr addr_t ZERO_ADDR_PAYLOAD = 1;
    unsigned int count_;
#endif
};

#endif /* _PHYSADDR_H_ */
