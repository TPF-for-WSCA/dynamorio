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

#include <iostream>
#include <sstream>
#ifdef LINUX
#    include <sys/types.h>
#    include <unistd.h>
#    include <sys/stat.h>
#    include <fcntl.h>
#endif
#include "physaddr.h"
#include "../common/options.h"
#include "../common/utils.h"

#if defined(X86_32) || defined(ARM_32)
#    define IF_X64_ELSE(x, y) y
#else
#    define IF_X64_ELSE(x, y) x
#endif

// XXX: can we share w/ core DR?
#define TESTALL(mask, var) (((mask) & (var)) == (mask))
#define TESTANY(mask, var) (((mask) & (var)) != 0)

#ifdef LINUX
#    define PAGEMAP_VALID 0x8000000000000000
#    define PAGEMAP_SWAP 0x4000000000000000
#    define PAGEMAP_PFN 0x007fffffffffffff
#    define PAGE_BITS 12 // XXX i#1703: handle large pages
#    define PAGE_START(addr) ((addr) & (~((1 << PAGE_BITS) - 1)))
#    define PAGE_OFFS(addr) ((addr) & ((1 << PAGE_BITS) - 1))
static const addr_t PAGE_INVALID = (addr_t)-1;
#endif

physaddr_t::physaddr_t()
#ifdef LINUX
    : last_vpage_(PAGE_INVALID)
    , last_ppage_(PAGE_INVALID)
    , fd_(-1)
    , count_(0)
#endif
{
#ifdef LINUX
    v2p_.table = nullptr;
#endif
}

physaddr_t::~physaddr_t()
{
#ifdef LINUX
    if (v2p_.table != nullptr)
        hashtable_delete(&v2p_);
#endif
}

bool
physaddr_t::init()
{
#ifdef LINUX
    // Some threads may not do much, so start out small.
    constexpr int V2P_INITIAL_BITS = 9;
    hashtable_init_ex(&v2p_, V2P_INITIAL_BITS, HASH_INTPTR, /*strdup=*/false,
                      /*synch=*/false, nullptr, nullptr, nullptr);

    // We avoid std::ostringstream to avoid malloc use for static linking.
    constexpr int MAX_PAGEMAP_FNAME = 64;
    char fname[MAX_PAGEMAP_FNAME];
    dr_snprintf(fname, BUFFER_SIZE_ELEMENTS(fname), "/proc/%d/pagemap", getpid());
    NULL_TERMINATE_BUFFER(fname);
    // We can't read pagemap with any buffered i/o, like ifstream, as we'll
    // get EINVAL on any non-8-aligned size, and ifstream at least likes to
    // read buffers of non-aligned sizes.
    fd_ = open(fname, O_RDONLY);
    // Accessing /proc/pid/pagemap requires privileges on some distributions,
    // such as Fedora with recent kernels.  We have no choice but to fail there.
    return (fd_ != -1);
#else
    // i#1727: we assume this is not possible on Windows.  If it is we
    // may want to split into physaddr_linux.cpp vs others.
    return false;
#endif
}

addr_t
physaddr_t::virtual2physical(addr_t virt)
{
#ifdef LINUX
    addr_t vpage = PAGE_START(virt);
    bool use_cache = true;
    if (op_virt2phys_freq.get_value() > 0 && ++count_ >= op_virt2phys_freq.get_value()) {
        // Flush the cache and re-sync with the kernel
        use_cache = false;
        last_vpage_ = PAGE_INVALID;
        hashtable_clear(&v2p_);
        count_ = 0;
    }
    if (use_cache) {
        // Use cached values on the assumption that the kernel hasn't re-mapped
        // this virtual page.
        if (vpage == last_vpage_)
            return last_ppage_ + PAGE_OFFS(virt);
        // XXX i#1703: add (debug-build-only) internal stats here and
        // on cache_t::request() fastpath.
        void *lookup = hashtable_lookup(&v2p_, reinterpret_cast<void *>(vpage));
        if (lookup != nullptr) {
            addr_t ppage = reinterpret_cast<addr_t>(lookup);
            // Restore a 0 payload.
            if (ppage == ZERO_ADDR_PAYLOAD)
                ppage = 0;
            last_vpage_ = vpage;
            last_ppage_ = ppage;
            return last_ppage_ + PAGE_OFFS(virt);
        }
    }
    // Not cached, or forced to re-sync, so we have to read from the file.
    if (fd_ == -1)
        return 0;
    // The pagemap file contains one 64-bit int per page, which we assume
    // here is 4096 bytes.
    // (XXX i#1703: handle large pages)
    // Thus we want offset:
    //   (addr / 4096 * 8) == ((addr >> 12) << 3) == addr >> 9
    if (lseek64(fd_, vpage >> 9, SEEK_SET) < 0)
        return 0;
    unsigned long long entry;
    if (read(fd_, (char *)&entry, sizeof(entry)) != sizeof(entry))
        return 0;
    if (!TESTALL(PAGEMAP_VALID, entry) || TESTANY(PAGEMAP_SWAP, entry))
        return 0;
    last_ppage_ = (addr_t)((entry & PAGEMAP_PFN) << PAGE_BITS);
    if (op_verbose.get_value() >= 2) {
        std::cerr << "virtual " << virt << " => physical "
                  << (last_ppage_ + PAGE_OFFS(virt)) << std::endl;
    }
    // Store 0 as a sentinel since 0 means no entry.
    hashtable_add(
        &v2p_, reinterpret_cast<void *>(vpage),
        reinterpret_cast<void *>(last_ppage_ == 0 ? ZERO_ADDR_PAYLOAD : last_ppage_));
    last_vpage_ = vpage;
    return last_ppage_ + PAGE_OFFS(virt);
#else
    return 0;
#endif
}
