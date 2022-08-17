/**
 * @file basic_block_stats.cpp
 * @author Roman Kaspar Brunner (roman.k.brunner@ntnu.no)
 * @brief Basic block analyzer for instruction specific cache.
 * @version 0.1
 * @date 2022-06-13
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <execution>
#include <sys/stat.h>

#include "basic_block_stats.h"
#include "../common/options.h"
// DEBUGGING IMPORTS: HOW MUCH MEMORY DO WE USE
#include "stdlib.h"
#include "string.h"

// DEBUGGING ABORT
// #define ANALYSED_INSTRUCTIONS_PER_ITERATION 50000000

/**
 * @brief Parses status line
 *
 * @param line The textual output of status file to be parsed
 * @return int The amount of resources used
 */
int
parseLine(char *line)
{
    int i = strlen(line);
    const char *p = line;
    while (*p < '0' || *p > '9')
        p++;
    line[i - 3] = '\0';
    i = atoi(p);
    return i;
}

/**
 * @brief Get the Physical RAM Usage Value
 *
 * @return int The amount of physical RAM used in KB
 */
int
getPhysicalRAMUsageValue()
{
    FILE *file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

void
trim_vector(std::vector<size_t> &vec)
{
    int i;
    for (i = vec.size(); i >= 0 && vec[i] == 0; --i)
        ;
    vec.resize(i + 1);
}

uint64_t
set_accessed(uint64_t mask, uint8_t lower, uint8_t upper)
{
    if (upper > 64) {
        upper = 64;
    }

    uint64_t bitmask = 1;
    for (uint8_t i = 0; i < lower; i++) {
        bitmask = bitmask << 1;
    }

    for (uint8_t i = lower; i < upper; i++) {
        mask |= bitmask;
        bitmask = bitmask << 1;
    }
    return mask;
}

uint64_t
get_total_mask_for_presence(const std::vector<uint64_t> &masks)
{
    uint64_t current_mask = 0;
    for (auto const &mask : masks) {
        current_mask |= mask;
    }
    return current_mask;
}

uint8_t
get_total_access_from_masks(const std::vector<uint64_t> &masks)
{
    return (uint8_t)__builtin_popcountll(get_total_mask_for_presence(masks));
}

/**
 * @brief Parses a bit mask to detect holes and blocks in a cacheline
 *
 * @param the mask of the cacheline under test
 * @return std::vector<int> A vector containing B H B H B sizes (Block and hole sizes),
 * where it alwazs needs to return an odd number of entries, as it alwazs has to start
 * and end with a block and in between tw blcks there is always exclusively one single
 * hole.
 */
std::vector<int>
count_holes_in_masks(const uint64_t &mask)
{
    std::vector<int> holes;
    uint8_t prev_bit = 0;
    uint8_t count_hole_size = 0;
    uint8_t count_block_size = 0;
    bool trailing = true;
    for (int byte = 0; byte < 64; byte++) {
        uint8_t current_bit = ((mask >> byte) & 0x1);
        if (current_bit == 1 && trailing) {
            trailing = false;
            prev_bit = current_bit;
            count_block_size = 1;
            continue;
        } else if (current_bit == 0 && trailing) {
            continue;
        }

        // posedge/negedge kind of analysis
        if (current_bit == 0 && prev_bit == 1) {
            holes.push_back(count_block_size);
            count_block_size = 0;
        } else if (current_bit == 1 && prev_bit == 0) {
            holes.push_back(count_hole_size);
            count_hole_size = 0;
        }

        if (current_bit == 0) {
            count_hole_size += 1;
        } else {
            count_block_size += 1;
        }
        prev_bit = current_bit;
    }

    if (prev_bit == 1 && count_block_size > 0) {
        holes.push_back(count_block_size);
    }

    if (holes.size() % 2 != 1) {
        std::cout << "WHAT IS WRONG WITH YOU???" << std::endl;
    }
    return holes;
}

/**
 * @brief assert with a message
 *
 */
#define assertm(exp, msg) assert(((void)msg, exp))

bool
operator==(const BasicBlock &lhs, const BasicBlock &rhs)
{
    return (rhs.starting_addr <= lhs.starting_addr && lhs.end_addr <= rhs.end_addr) ||
        (lhs.starting_addr <= rhs.starting_addr && rhs.end_addr <= lhs.end_addr);
}

/**
 * @brief Adjust the count in the vector, where the index denotes how large the
 * basic block is
 *
 * @param vec Vector contains the counts for the different sizes of basic blocks
 * @param size The size of the basic block which we want to count
 */
void
insert_bbcount(std::vector<size_t> &vec, size_t size)
{
    try {
        if (vec.size() <= size) {
            vec.resize(2 * size, 0);
        }
        vec[size] += 1;
    } catch (std::length_error &e) {
        printf("FAILED TO ADD DATA TO VECTOR.\n");
        printf("Current vector size: %lu\n", vec.size());
        printf("size              : %lu\n", size);
        std::cout << e.what() << std::endl;
    }
}

basic_block_stats_t::basic_block_stats_t(int block_size, const std::string &miss_file,
                                         const std::string &output_dir,
                                         bool warmup_enabled, bool is_coherent)
    : caching_device_stats_t(miss_file, block_size, warmup_enabled, is_coherent)
    , count_per_basic_block_instr_size_(block_size)
    , count_per_basic_block_byte_size_(block_size)
    , basic_block_size_history(0)
    , basic_blocks_hit_count(0)
    , output_dir(output_dir)
    , current_block({ 0, 0, 0, 0 })
    , basic_block_buffer()
    , distinct_cachelines_by_size(65)
{
    basic_block_size_history.reserve(10000);
    // basic_blocks_hit_count.reserve(10000);

    std::cout << "creating output dir: " << output_dir << std::endl;
    std::string cmd = "/bin/mkdir -p " + output_dir;
    const int dir_err = system(cmd.c_str());
    if (0 != dir_err) {
        std::cout << "Sth is wrong with the output dir: " << output_dir << std::endl;
        exit(1);
    }
}

void
basic_block_stats_t::record_cacheblock(const caching_device_block_t &t)
{
    if (cache_block_buffer.size() == 10) {
        cache_block_buffer.pop_front();
    }

    cache_block_buffer.push_back(t);
}

void
basic_block_stats_t::print_last_n_cacheblocks(int n)
{
    int start = cache_block_buffer.size() - 1;
    n = start - n;
    if (n < 0)
        n = 0;

    for (int i = start; i >= n; i++) {
        auto curr = cache_block_buffer[i];
        std::cout << "cacheblock: " << i << "\n";
        std::cout << "\tTag: " << curr.tag_ << "\n";
        std::cout << "\tCounter: " << curr.counter_ << "\n";
    }
    std::cout << std::flush;
}

void
basic_block_stats_t::record_memref(memref_t mem, bool hit)
{
    if (basic_block_buffer.size() == 10) {
        basic_block_buffer.pop_front();
    }

    basic_block_buffer.push_back(
        BasicBlock { .starting_addr = mem.instr.addr,
                     .end_addr = (mem.instr.addr + mem.instr.size - 1),
                     .instr_size = 1,
                     .byte_size = mem.instr.size,
                     .hits = 1,
                     .miss = !hit });
}

void
basic_block_stats_t::print_last_n_memrefs(int n)
{
    int start = basic_block_buffer.size() - 1;
    n = start - n;
    if (n < 0) {
        n = 0;
    }

    for (int i = start; i >= n; --i) {
        auto curr = basic_block_buffer[i];
        std::cout << "memref -" << (i + 1) << "\n";
        std::cout << "start addr: " << curr.starting_addr << "\n";
        std::cout << "size: " << curr.byte_size << "\n";
        // print_type(curr.instr.type);  // TODO: If required adjust bb struct
        std::cout << "miss: " << curr.miss << "\n";
    }
    std::cout << std::flush;
}

bool
basic_block_stats_t::handle_instr(const memref_t &memref, bool hit)
{
    bool is_interrupt_ = false;
    if (current_block.starting_addr == 0) {
        current_block.miss = !hit;
        current_block.starting_addr = memref.data.addr;
        current_block.end_addr = memref.data.addr;
    }

    if (current_block_cacheline_constrained.starting_addr == 0) {
        reset_current_cacheline_block(memref, hit);
    }

    // TODO: Given we want to track actual hit rates we need to adjust for cache lines
    // here -> Two lists:
    //      One per basic block, without recording hits/misses
    //      One per block per cacheline with recorded hits/misses, misses only recorded
    //      when we access the first byte of the block

    // NOTE: When an instruction overlaps into the next cacheline we will see an
    // additional access to the overlapping bytes
    auto prev_cacheline_base_address =
        current_block_cacheline_constrained.starting_addr & cache_line_address_mask;
    auto curr_cacheline_base_address = memref.instr.addr & cache_line_address_mask;

    // go by current unconstrained basic block as cacheline constraints lead to more
    // frequent updates
    bool is_adjacent_instr = (memref.data.addr ==
                              (current_block_cacheline_constrained.starting_addr +
                               current_block_cacheline_constrained.byte_size));

    bool aliasing_eviction = (is_adjacent_instr && !hit && prev_hit);

    // We do not orient ourselves at branches but only at cacheline base_addresses
    // and corresponding jumps
    if (prev_cacheline_base_address != curr_cacheline_base_address ||
        !is_adjacent_instr || aliasing_eviction) {
        // TODO: Track new cacheline constrained block
        uint8_t lower = (current_block_cacheline_constrained.starting_addr -
                         prev_cacheline_base_address);
        uint8_t upper =
            (current_block_cacheline_constrained.starting_addr +
             current_block_cacheline_constrained.byte_size - prev_cacheline_base_address);
        uint64_t mask = set_accessed(0, lower, upper);
        if (!current_block_cacheline_constrained.miss &&
            !bytes_accessed_per_presence_per_cacheline[prev_cacheline_base_address]
                 .empty()) {
            // still present - add current access
            bytes_accessed_per_presence_per_cacheline[prev_cacheline_base_address]
                .back()           // get last vector = current present line counter
                .push_back(mask); // add current access
        } else {
            // newly loaded into cache - add instance
            // or we are after warmup - we start counting from 0
            bytes_accessed_per_presence_per_cacheline[prev_cacheline_base_address]
                .push_back({ mask });
        }

        reset_current_cacheline_block(memref, hit);
        // We don't need to fixup this if it is only a cacheline overrun, as this is fixed
        // in belows general update logic
        if (!is_adjacent_instr || aliasing_eviction) {
            current_block_cacheline_constrained.byte_size = memref.instr.size;
            current_block_cacheline_constrained.instr_size = 1;
        }
    }

    if (is_adjacent_instr && !aliasing_eviction) {
        // If we jumped addresses (without seeing a branch), then we have encountered
        // an interrupt/scheduling event
        current_block.instr_size++;
        current_block.byte_size += memref.data.size;
        current_block.end_addr = memref.data.addr;

        current_block_cacheline_constrained.instr_size++;
        current_block_cacheline_constrained.byte_size += memref.data.size;
        current_block_cacheline_constrained.end_addr = memref.data.addr;
        if (current_block_cacheline_constrained.byte_size > 64) {
            std::cout << "WTF is going on" << std::endl;
        }

    } else if (aliasing_eviction && is_adjacent_instr) {
        // if not adjacent: we were and interrupt/
        // TODO: count alias evictions
        // do we need a seperate return value? or anything else?
    } else {
        // We moved more than the size of the previous x86 instruction - we must have been
        // interrupted
        num_interrupts++;
        is_interrupt_ = true;
    }
    if ((MAX_X86_INSTR_SIZE <
         (current_block.byte_size / (float)current_block.instr_size))) {
        printf("WARN: SUCH LARGE x86 INSTRUCTIONS DO NOT EXIST: %10.2f\n",
               current_block.byte_size / (float)current_block.instr_size);
    }

    if (MAX_X86_INSTR_SIZE < memref.data.size) {
        printf("WARN: SUCH LARGE x86 INSTRUCTIONS DO NOT EXIST: %lu\n", memref.data.size);
    }

    if (max_instr_size < memref.data.size) {
        max_instr_size = memref.data.size;
    }

    prev_hit = hit;
    return is_interrupt_;
}

void
basic_block_stats_t::track_cacheline_access(const memref_t &memref)
{
    addr_t cacheline_start_address =
        cache_line_address_mask & current_block.starting_addr;
    addr_t cacheline_end_address =
        cache_line_address_mask & (current_block.starting_addr + current_block.byte_size);
    cacheline_end_address = (cacheline_end_address == cacheline_start_address)
        ? cacheline_end_address + 64
        : cacheline_end_address;

    if (cacheline_end_address - cacheline_start_address > max_cacheline_bb) {
        max_cacheline_bb = cacheline_end_address - cacheline_start_address;
        std::cout << "Current largest overlap of a basic block: " << max_cacheline_bb / 64
                  << " lines" << std::endl;
        std::cout << "\tBB START: " << current_block.starting_addr << "\n\tBB END: "
                  << current_block.starting_addr + current_block.byte_size
                  << "\n\tBB INSTR COUNT: " << current_block.instr_size << std::endl;
    }

    for (; cacheline_start_address < cacheline_end_address;
         cacheline_start_address += 64) {
        auto basic_block_vec = &number_of_bytes_accessed[cacheline_start_address];
        // if it spans more than one line it might be that the vec is still empty
        // if it spans more than one line it might be that the miss handling per basic
        // block is overly pessimistic
        if (current_block.miss || basic_block_vec->size() == 0) {
            basic_block_vec->push_back(current_block);
        } else {
            auto is_same_block = [this](const BasicBlock &bb) {
                return (bb.starting_addr <= current_block.starting_addr) &&
                    (bb.end_addr >= current_block.end_addr);
            };
            auto bb_entry = std::find_if(basic_block_vec->rbegin(),
                                         basic_block_vec->rend(), is_same_block);

            if (bb_entry == basic_block_vec->rend()) {
                basic_block_vec->push_back(current_block);
            } else {
                bb_entry->hits += 1;
            }
        }
    }
}

void
basic_block_stats_t::handle_branch(const memref_t &memref)
{
    // assuming x86, instructions in basic blocks have monotoneously increasing
    // addresses
    // Deactivate for now: only tracing cachelines
    // if (current_block.starting_addr <= current_block.end_addr &&
    //     current_block.starting_addr != 0) {
    //     record_block(current_block);
    // } else {
    //     printf("WARN: ADDRESS MISMATCH, BLOCK END IS SMALLER THAN BLOCK START\n");
    // }
    // insert_bbcount(count_per_basic_block_byte_size_, current_block.byte_size);
    // insert_bbcount(count_per_basic_block_instr_size_, current_block.instr_size);

    //  basic_blocks_hit_count[current_block] += 1;
    //  basic_block_size_history.push_back(current_block.byte_size);

    // track_cacheline_access(memref);

    if (basic_block_size_history.size() == basic_block_size_history.capacity()) {
        basic_block_size_history.reserve(basic_block_size_history.capacity() * 2);
        size_t current_ram_consumption = getPhysicalRAMUsageValue() >> 10;
        // how much memory do I consume?
        std::cout << "CURRENT RAM: " << current_ram_consumption << " MB" << std::endl;
    }
    // current_block = { .starting_addr = 0, .end_addr = 0 };
    current_block.byte_size = 0;
    current_block.instr_size = 0;
    current_block.starting_addr = 0;
    current_block.end_addr = 0;
}

void
basic_block_stats_t::handle_interrupt(const memref_t &memref, bool hit)
{
    current_block.starting_addr = memref.data.addr;
    current_block.end_addr = memref.data.addr;
    current_block.byte_size = memref.data.size;
    current_block.instr_size = 1;
    current_block.miss = !hit;
    current_block_cacheline_constrained.starting_addr = memref.data.addr;
    current_block_cacheline_constrained.end_addr = memref.data.addr;
    current_block_cacheline_constrained.byte_size = memref.data.size;
    current_block_cacheline_constrained.instr_size = 1;
    current_block_cacheline_constrained.miss = !hit;
}

void
basic_block_stats_t::access(const memref_t &memref, bool hit,
                            caching_device_block_t *cache_block)
{
    // TODO: Handle misses even if we are adjacent instructions
    // record_cacheblock(*cache_block);
    caching_device_stats_t::access(memref, hit, cache_block);

    // size_t current_ram_consumption = getPhysicalRAMUsageValue() >> 10;
    // if (max_memory_consumption < current_ram_consumption) {
    //     std::cout << "CURRENT RAM: " << current_ram_consumption << " MB" << std::endl;
    // }

    //    double progress =
    //        ((double)handled_instructions / ANALYSED_INSTRUCTIONS_PER_ITERATION) *
    //        100.0;
    if (handled_instructions % 1000000 == 0) {
        std::cout << "PROGRESS " << handled_instructions << " instructions handled\n";

        size_t current_ram_consumption = getPhysicalRAMUsageValue() >> 10;
        // how much memory do I consume?
        std::cout << "CURRENT RAM: " << current_ram_consumption << " MB\n";
        std::cout << std::flush;
    }

    if (type_is_prefetch(memref.data.type)) {
        std::cout << "YES SAW A PREFETCH\n";
        return; // we only care about executed instrs
    }

    if (memref.data.type == TRACE_TYPE_THREAD ||
        memref.data.type == TRACE_TYPE_THREAD_EXIT) {
        printf("WARN: we do not expect thread entries in this part of the reader. Still "
               "we saw one. Ignoring\n");
        return;
    }

    bool is_instr_ = type_is_instr(memref.data.type);
    bool is_branch_ = type_is_instr_branch(memref.data.type);

    if (!is_branch_ && !is_instr_) {
        print_type(memref.data.type);
        return; // we do not care about anything that is not an instruction and not a
                // branch
    }

    handled_instructions++;

    bool is_interrupt_ = false;

    // count number of instructions in basic block
    if (is_instr_) {
        is_interrupt_ = handle_instr(memref, hit);
        is_branch_ = is_branch_ || is_interrupt_;
    }
    // record_memref(memref, hit);

    // calculate number of bytes in instruction based on addresses.
    if (is_branch_) {
        //        if (current_block.byte_size == 1)
        //            print_last_n_memrefs(1);
        handle_branch(memref);
    }

    if (is_interrupt_) {
        // if (current_block.byte_size == 1)
        //     print_last_n_memrefs(9);
        handle_interrupt(memref, hit);
    }

    if (!hit) {
        auto base_addr = memref.instr.addr & cache_line_address_mask;
        auto it = eviction_to_fetch_index_map.find(base_addr);
        int prev_idx = -1;
        if (it != eviction_to_fetch_index_map.end()) {
            prev_idx = it->second;
        }

        if (prev_idx != -1) {
            auto presences = bytes_accessed_per_presence_per_cacheline[base_addr];
            auto last_presence = presences.back();
            auto total_accesses = get_total_mask_for_presence(last_presence);
            perfect_fetch_history[prev_idx] = std::pair(base_addr, total_accesses);
        }
        // ELSE: This is a compulsory miss

        perfect_fetch_history.push_back(std::pair(base_addr, 0));
        eviction_to_fetch_index_map[base_addr] = perfect_fetch_history.size() - 1;
    }
    // if (handled_instructions % ANALYSED_INSTRUCTIONS_PER_ITERATION == 0) {
    //     std::ostringstream oss;
    //     oss << "Stats after " << handled_instructions << " instructions:\n";
    //     print_stats(oss.str());
    //     oss.str("");
    //     oss.clear();
    //     // short hack to keep instr count
    //     auto prev = handled_instructions;
    //     reset();
    //     handled_instructions = prev;
    // }
}

/**
 * @brief Aggregate access counts per contiguous blocks with the same access patterns.
 * The current implementation could lead to falsely counting too long blocks.
 *
 * @param access_masks Vector of masks indicating which bytes were accessed contiguously
 * @return std::vector<uint64_t> Vector of size 64, counting the number of accesses to
 * blocks of size (idx + 1).
 */
std::vector<uint64_t>
basic_block_stats_t::aggregate_byte_accesses_in_cacheline_presence(
    const std::vector<uint64_t> &access_masks)
{
    std::vector<int> byte_bucket_indices(64);
    std::vector<uint64_t> byte_bucket(64, 0);
    std::iota(std::begin(byte_bucket_indices), std::end(byte_bucket_indices), 0);
    std::for_each(std::execution::par, std::begin(byte_bucket_indices),
                  std::end(byte_bucket_indices),
                  [&access_masks, &byte_bucket](auto &&byte) {
                      for (const auto &mask : access_masks) {
                          if (((mask >> byte) & 0x1) == 1) {
                              byte_bucket[byte]++;
                          }
                      }
                  });

    std::vector<uint64_t> result(65, 0);
    uint64_t access_count = 0;
    uint64_t byte_count = 0;
    for (auto &loc_access_count : byte_bucket) {
        if (access_count == loc_access_count) {
            byte_count += 1;
        } else {
            result[byte_count] += access_count;
            access_count = loc_access_count;
            byte_count = 1;
        }
    }
    if (access_count != 0) {
        result[byte_count] += access_count;
    }

    result.erase(result.begin());
    return result;
}

/**
 * @brief Aggregate access counts over multiple presences of a single memory region, split
 * by cache line.
 *
 * @param cacheline_presences Vector containing all presences of current cacheline.
 * @return std::vector<uint64_t> Vector of size 64 indicating the number of accesses to
 * blocks of size (idx + 1) within the current cacheline.
 */
std::vector<uint64_t>
basic_block_stats_t::aggregate_byte_accesses_in_cacheline(
    const std::vector<std::vector<uint64_t>> &cacheline_presences)
{
    std::vector<uint64_t> result(64, 0);
    for (auto &masks : cacheline_presences) {
        auto aggregated_masks = aggregate_byte_accesses_in_cacheline_presence(masks);
        std::transform(result.begin(), result.end(), aggregated_masks.begin(),
                       result.begin(), std::plus<uint64_t>());
    }
    return result;
}

void
print_type(trace_type_t type)
{
    switch (type) {
    case TRACE_TYPE_READ: printf("READ\n"); break;
    case TRACE_TYPE_WRITE: printf("WRITE\n"); break;
    case TRACE_TYPE_PREFETCH: printf("PREFETCH\n"); break;
    // X86 specific prefetch
    case TRACE_TYPE_PREFETCHT0: printf("PREFETCHT0\n"); break;
    case TRACE_TYPE_PREFETCHT1: printf("PREFETCHT1\n"); break;
    case TRACE_TYPE_PREFETCHT2: printf("PREFETCHT2\n"); break;
    case TRACE_TYPE_PREFETCHNTA: printf("PREFETCHNTA\n"); break;
    case TRACE_TYPE_PREFETCH_READ: printf("PREFETCH_READ\n"); break;
    case TRACE_TYPE_PREFETCH_WRITE: printf("PREFETCH_WRITE\n"); break;
    case TRACE_TYPE_PREFETCH_INSTR: printf("PREFETCH_INSTR\n"); break;
    case TRACE_TYPE_INSTR: printf("INSTR\n"); break;
    case TRACE_TYPE_INSTR_DIRECT_JUMP: printf("INSTR_DIRECT_JUMP\n"); break;
    case TRACE_TYPE_INSTR_INDIRECT_JUMP: printf("INSTR_INDIRECT_JUMP\n"); break;
    case TRACE_TYPE_INSTR_CONDITIONAL_JUMP: printf("INSTR_CONDITIONAL_JUMP\n"); break;
    case TRACE_TYPE_INSTR_DIRECT_CALL: printf("INSTR_DIRECT_CALL\n"); break;
    case TRACE_TYPE_INSTR_INDIRECT_CALL: printf("INSTR_INDIRECT_CALL\n"); break;
    case TRACE_TYPE_INSTR_RETURN: printf("INSTR_RETURN\n"); break;
    case TRACE_TYPE_INSTR_BUNDLE: printf("INSTR_BUNDLE\n"); break;
    case TRACE_TYPE_INSTR_FLUSH: printf("INSTR_FLUSH\n"); break;
    case TRACE_TYPE_INSTR_FLUSH_END: printf("INSTR_FLUSH_END\n"); break;
    case TRACE_TYPE_DATA_FLUSH: printf("DATA_FLUSH\n"); break;
    case TRACE_TYPE_DATA_FLUSH_END: printf("DATA_FLUSH_END\n"); break;
    case TRACE_TYPE_THREAD: printf("THREAD\n"); break;
    case TRACE_TYPE_MARKER: printf("MARKER\n"); break;
    case TRACE_TYPE_THREAD_EXIT: printf("THREAD_EXIT\n"); break;
    case TRACE_TYPE_PID: printf("PID\n"); break;
    case TRACE_TYPE_HEADER: printf("HEADER\n"); break;
    case TRACE_TYPE_FOOTER: printf("FOOTER\n"); break;
    case TRACE_TYPE_HARDWARE_PREFETCH: printf("HARDWARE_PREFETCH\n"); break;
    case TRACE_TYPE_INSTR_MAYBE_FETCH: printf("INSTR_MAYBE_FETCH\n"); break;
    case TRACE_TYPE_INSTR_SYSENTER: printf("INSTR_SYSENTER\n"); break;
    case TRACE_TYPE_PREFETCH_READ_L1_NT: printf("PREFETCH_READ_L1_NT\n"); break;
    case TRACE_TYPE_PREFETCH_READ_L2_NT: printf("PREFETCH_READ_L2_NT\n"); break;
    case TRACE_TYPE_PREFETCH_READ_L3_NT: printf("PREFETCH_READ_L3_NT\n"); break;
    case TRACE_TYPE_PREFETCH_INSTR_L1_NT: printf("PREFETCH_INSTR_L1_NT\n"); break;
    case TRACE_TYPE_PREFETCH_INSTR_L2: printf("PREFETCH_INSTR_L2\n"); break;
    case TRACE_TYPE_PREFETCH_INSTR_L2_NT: printf("PREFETCH_INSTR_L2_NT\n"); break;
    case TRACE_TYPE_PREFETCH_INSTR_L3: printf("PREFETCH_INSTR_L3\n"); break;
    case TRACE_TYPE_PREFETCH_INSTR_L3_NT: printf("PREFETCH_INSTR_L3_NT\n"); break;
    case TRACE_TYPE_PREFETCH_WRITE_L1_NT: printf("PREFETCH_WRITE_L1_NT\n"); break;
    case TRACE_TYPE_PREFETCH_WRITE_L2: printf("PREFETCH_WRITE_L2\n"); break;
    case TRACE_TYPE_PREFETCH_WRITE_L2_NT: printf("PREFETCH_WRITE_L2_NT\n"); break;
    case TRACE_TYPE_PREFETCH_WRITE_L3: printf("PREFETCH_WRITE_L3\n"); break;
    case TRACE_TYPE_PREFETCH_WRITE_L3_NT: printf("PREFETCH_WRITE_L3_NT\n"); break;
    default: printf("NOT HANDLED CASE OF %d\n", type); break;
    }
}

/*
void
basic_block_stats_t::record_block(const BasicBlock &bb)
{
    auto target_bb_node = basic_blocks_hit_count.find(bb);
    if (target_bb_node == basic_blocks_hit_count.end()) {
        basic_blocks_hit_count.insert({ bb, 1 });
        return;
    }

    // key and value both return references
    auto target_bb = target_bb_node->first;
    BasicBlock merged_bb = target_bb;
    if (bb.starting_addr < target_bb.starting_addr) {
        merged_bb.starting_addr = bb.starting_addr;
    }

    if (bb.end_addr > target_bb.end_addr) {
        printf("WARN: BLOCK ENDS EXTENDS - SHOULD NEVER HAPPEN AS WE ALWAYS FINISH BB AT "
               "BRANCHING INSTR\n");
        target_bb.end_addr = bb.end_addr;
    }

    target_bb_node->second += 1;
}
*/

uint64_t
basic_block_stats_t::bytes_accessed_by_block(const addr_t &cacheline_base,
                                             BasicBlock &block)
{
    addr_t start_block = block.starting_addr;
    addr_t end_block = start_block + block.byte_size;
    if (block.starting_addr < cacheline_base) {
        start_block = cacheline_base;
    }
    if ((cacheline_base + 64) < block.end_addr) {
        end_block = cacheline_base + 64;
    }
    start_block -= cacheline_base;
    end_block -= cacheline_base;

    auto result = set_accessed(0, start_block, end_block);
    return result;
}

std::pair<uint8_t, std::vector<hitbytes>>
basic_block_stats_t::bytes_accessed(const addr_t &cacheline_base,
                                    std::vector<BasicBlock> &blocks_contained)
{
    uint64_t mask = 0;
    std::vector<hitbytes> accesses;
    accesses.reserve(blocks_contained.size());
    for (auto it = blocks_contained.begin(); it != blocks_contained.end(); it++) {
        auto tmp_mask = bytes_accessed_by_block(cacheline_base, *it);
        auto tmp_count = (uint8_t)__builtin_popcountll(tmp_mask);
        accesses.push_back(hitbytes(it->hits, tmp_count));
        mask |= tmp_mask;
    }

    return std::pair<uint8_t, std::vector<hitbytes>> {
        (uint8_t)__builtin_popcountll(mask), accesses
    };
}

std::vector<size_t>
get_moving_average_block_size(std::vector<size_t> &history, int window_size)
{
    std::vector<size_t> result;
    result.reserve(history.size());
    for (size_t i = 0; i < history.size(); i += 1) {
        auto it = history.begin() + i;
        size_t sum = std::accumulate(it, it + window_size, 0);
        result.push_back(sum / window_size);
    }
    return result;
}

std::vector<size_t>
get_local_average_block_size(std::vector<size_t> &history, int window_size)
{
    std::vector<size_t> result;
    result.reserve(history.size() / window_size + 1);
    for (size_t i = 0; i < history.size(); i += window_size) {
        if (history.size() < (size_t)window_size + i - 2)
            break;
        auto it = history.begin() + i;
        size_t sum = std::accumulate(it, it + window_size, 0);
        result.push_back(sum / window_size);
    }
    return result;
}

void
tikz_print_series()
{
}

// overload fnct for stacked relative vs series
void
tikz_print_histo()
{
}

void
basic_block_stats_t::print_bytes_accessed()
{
    std::vector<uint64_t> histogram(65, 0);
    std::vector<double> overall_accessed_histgram(65, 0);
    uint64_t total_allocations = 0;

    for (auto it = number_of_bytes_accessed.begin(); it != number_of_bytes_accessed.end();
         it++) {
        std::vector<BasicBlock> vec = (*it).second;
        total_allocations += vec.size();
        const addr_t base_addr = (*it).first;
        auto accessed = bytes_accessed(base_addr, vec);
        total_allocations +=
            create_histogram_of_cachelineaccesses(histogram, accessed.second);
        int total_accessed = (int)accessed.first;
        overall_accessed_histgram[total_accessed]++;
        // TODO: Integrate this into counter
        // std::cout << base_addr << ": " << total_accessed << std::endl;
    }

    std::vector<double> relative_histogram;
    relative_histogram.reserve(65);

    for (auto const &bucket : histogram) {
        relative_histogram.push_back((long double)bucket /
                                     (long double)total_allocations);
    }

    std::vector<std::vector<size_t>>
        count_accesses_of_accessed_bytes_per_total_accessed_bytes_of_cacheline(
            65, std::vector<size_t>(65, 0));

    std::map<size_t, size_t> distribution_of_presences;

    //    std::vector<double> count_accessed_bytes_per_cacheline(65, 0.0);
    std::vector<uint64_t> access_sizes_to_cache(64, 0);
    std::vector<uint64_t> num_lines_with_accesses_of_size(64, 0);
    std::vector<int> holes(64);
    std::vector<std::pair<int, double>>
        hole_sizes_and_relative; // (hole size, relative part of the local block (B H B))
    // TODO: Replace by parallel foreach loops
    for (auto const &[cacheline_baseaddress, accesses_per_presence] :
         bytes_accessed_per_presence_per_cacheline) {
        auto tmp_result = aggregate_byte_accesses_in_cacheline(accesses_per_presence);
        std::transform(access_sizes_to_cache.begin(), access_sizes_to_cache.end(),
                       tmp_result.begin(), access_sizes_to_cache.begin(),
                       std::plus<uint64_t>());
        distribution_of_presences[accesses_per_presence.size()] += 1;
        for (auto const &accesses : accesses_per_presence) {
            uint8_t total_accessed_bytes = get_total_access_from_masks(accesses);
            auto curr_presence_holes =
                count_holes_in_masks(get_total_mask_for_presence(accesses));
            int num_holes = curr_presence_holes.size() / 2;
            holes[num_holes] += 1;
            for (size_t i = 1; i < curr_presence_holes.size(); i += 2) {
                int curr_hole_size = curr_presence_holes[i];
                hole_sizes_and_relative.push_back(
                    std::pair(curr_hole_size,
                              (double)curr_hole_size /
                                  (curr_presence_holes[i - 1] +
                                   curr_presence_holes[i + 1] + curr_hole_size)));
            }
            distinct_cachelines_by_size[total_accessed_bytes].insert(
                cacheline_baseaddress);
            //            count_accessed_bytes_per_cacheline[total_accessed_bytes]++;
            for (auto const &mask : accesses) {
                // TODO: Assert that no mask is non-contiguous
                // count and insert
                auto accessed_bytes = __builtin_popcountll(mask);
                if (accessed_bytes > total_accessed_bytes) {
                    std::cout << "OH OH THATS A HUGE NO NO" << std::endl;
                    exit(-1);
                }
                count_accesses_of_accessed_bytes_per_total_accessed_bytes_of_cacheline
                    [total_accessed_bytes][accessed_bytes]++;
            }
        }

        for (size_t i = 0; i < tmp_result.size(); i++) {
            if (tmp_result[i] == 0)
                continue;
            num_lines_with_accesses_of_size[i] += 1;
        }
    }

    std::ofstream perfect_loading_decisions_csv;
    std::cout << "Write file to " << output_dir << std::endl;
    perfect_loading_decisions_csv.open((output_dir + "/perfect_loading.csv"),
                                       std::ios::out);
    for (auto const &[base_addr, idx] : eviction_to_fetch_index_map) {
        auto presences = bytes_accessed_per_presence_per_cacheline[base_addr];
        auto last_presence = presences.back();
        auto total_accesses = get_total_mask_for_presence(last_presence);
        perfect_fetch_history[idx] = std::pair(base_addr, total_accesses);
    }
    if (!perfect_loading_decisions_csv) {
        std::cerr << "COULD NOT CREATE/WRITE FILE " << output_dir
                  << "/perfect_loading.csv\n";
        std::cout << "==== PERFECT LOADING DECISIONS ====\n";
        for (const auto &addr_mask_pair : perfect_fetch_history) {
            std::cout << addr_mask_pair.first << "; " << addr_mask_pair.second << "\n";
        }
    } else {
        perfect_loading_decisions_csv << "sep=;\n";
        perfect_loading_decisions_csv << "64b aligned addr;mask\n";
        for (const auto &addr_mask_pair : perfect_fetch_history) {
            perfect_loading_decisions_csv << addr_mask_pair.first << "; "
                                          << addr_mask_pair.second << "\n";
        }
        perfect_loading_decisions_csv.flush();
        perfect_loading_decisions_csv.close();
    }

    std::ofstream distinct_cachelines_csv;
    std::cout << "Write file to " << output_dir << std::endl;
    distinct_cachelines_csv.open((output_dir + "/distinct_cachelines.csv"),
                                 std::ios::out);
    if (!distinct_cachelines_csv) {
        std::cerr << "COULD NOT CREATE/WRITE FILE " << output_dir
                  << "/distinct_cachelines.csv\n";

        std::cout << "==== DISTINCT COUNT OF CACHELINES OF ACCESSED REGION SIZE ====\n";
        for (size_t i = 0; i < distinct_cachelines_by_size.size(); i++) {
            std::cout << i << ": " << distinct_cachelines_by_size[i].size() << "\n";
        }
    } else {
        distinct_cachelines_csv << "sep=;\n";
        distinct_cachelines_csv << "ACCESSED BYTES;COUNT OF DISTINCT CACHELINES\n";
        for (size_t i = 0; i < distinct_cachelines_by_size.size(); i++) {
            distinct_cachelines_csv << i << ";" << distinct_cachelines_by_size[i].size()
                                    << "\n";
        }
        distinct_cachelines_csv.flush();
        distinct_cachelines_csv.close();
    }

    std::ofstream holes_in_cachelines_csv;
    std::ofstream num_cl_with_num_holes_csv;
    holes_in_cachelines_csv.open((output_dir + "/holes_in_cachelines.csv"),
                                 std::ios::out);
    num_cl_with_num_holes_csv.open(
        (output_dir + "/number_of_cachelines_with_number_of_holes.csv"), std::ios::out);

    if (!holes_in_cachelines_csv || !num_cl_with_num_holes_csv) {
        std::stringstream holes_in_cachelines;
        std::stringstream num_cl_with_num_holes;
        std::cerr << "COULD NOT CREATE/WRITE FILE " << output_dir
                  << "/holes_in_cachelines.csv v "
                     "number_of_cachelines_with_number_of_holes.csv\n";
        for (size_t i = 0; i < holes.size(); i++) {
            num_cl_with_num_holes << i << ";" << holes[i] << "\n";
        }
        for (const auto &hole_information : hole_sizes_and_relative) {
            holes_in_cachelines << hole_information.first << "; "
                                << hole_information.second << "\n";
        }
        std::cout << "==== SIZE OF HOLES IN CACHELINES ====\n";
        std::cout << holes_in_cachelines.rdbuf();
        std::cout << "==== #CACHELINES WITH #HOLES ====\n";
        std::cout << num_cl_with_num_holes.rdbuf();
        std::cout << std::flush;
    } else {
        holes_in_cachelines_csv
            << "SIZE OF HOLES IN CACHELINES;RELATIVE CONTRIBUTION TO LOCAL BLOCK\n";
        for (size_t i = 0; i < holes.size(); i++) {
            num_cl_with_num_holes_csv << i << ";" << holes[i] << "\n";
        }
        for (const auto &hole_information : hole_sizes_and_relative) {
            holes_in_cachelines_csv << hole_information.first << "; "
                                    << hole_information.second << "\n";
        }
        holes_in_cachelines_csv.flush();
        holes_in_cachelines_csv.close();
        num_cl_with_num_holes_csv.flush();
        num_cl_with_num_holes_csv.close();
    }

    // for (size_t i = 0; i < count_accessed_bytes_per_cacheline.size(); i++) {
    //     std::cout << i << ": " << count_accessed_bytes_per_cacheline[i] <<
    //     std::endl;
    // }

    std::ofstream cacheline_load_counts;
    cacheline_load_counts.open((output_dir + "/cacheline_load_counts.csv"),
                               std::ios::out);
    if (!cacheline_load_counts) {
        std::cerr << "COULD NOT CREATE/WRITE FILE " << output_dir
                  << "/cacheline_load_counts.csv\n";

        std::cout << "==== NUMBER OF CACHELINES SEEN X TIMES ====\n";
        for (auto const &[loads, count] : distribution_of_presences) {
            std::cout << loads << ": " << count << "\n";
        }
    } else {
        cacheline_load_counts << "sep=;\n";
        cacheline_load_counts << "#LOADS;#CACHELINES\n";
        for (auto const &[loads, count] : distribution_of_presences) {
            cacheline_load_counts << loads << ";" << count << "\n";
        }
        cacheline_load_counts.flush();
        cacheline_load_counts.close();
    }

    std::ofstream access_count_to_cl_of_size;
    access_count_to_cl_of_size.open((output_dir + "/access_count_to_cl_of_size.csv"),
                                    std::ios::out);
    if (!access_count_to_cl_of_size) {
        std::cerr << "COULD NOT CREATE/WRITE FILE " << output_dir
                  << "/access_count_to_cl_of_size.csv\n";

        std::cout << "==== NUMBER OF ACCESSES TO LINE OF SIZE ====\n";
        for (size_t i = 0; i < access_sizes_to_cache.size(); i++) {
            std::cout << (i + 1) << ": " << access_sizes_to_cache[i] << "\n";
        }
    } else {
        access_count_to_cl_of_size << "sep=;\n";
        access_count_to_cl_of_size << "ACCESS SIZE;#ACCESSES\n";
        for (size_t i = 0; i < access_sizes_to_cache.size(); i++) {
            access_count_to_cl_of_size << (i + 1) << ";" << access_sizes_to_cache[i]
                                       << "\n";
        }
        access_count_to_cl_of_size.flush();
        access_count_to_cl_of_size.close();
    }

    std::ofstream number_of_cachelines_size;
    number_of_cachelines_size.open((output_dir + "/number_of_cachelines_size.csv"),
                                   std::ios::out);
    if (!number_of_cachelines_size) {
        std::cerr << "COULD NOT CREATE/WRITE FILE " << output_dir
                  << "/number_of_cachelines_size.csv\n";

        std::cout << "==== NUMBER OF CACHELINES OF SIZE ====\n";
        for (size_t i = 0; i < num_lines_with_accesses_of_size.size(); i++) {
            std::cout << (i + 1) << ": " << num_lines_with_accesses_of_size[i] << "\n";
        }

    } else {
        number_of_cachelines_size << "sep=;\n";
        number_of_cachelines_size << "ACCESS SIZE;#CACHELINES\n";
        for (size_t i = 0; i < num_lines_with_accesses_of_size.size(); i++) {
            number_of_cachelines_size << (i + 1) << ";"
                                      << num_lines_with_accesses_of_size[i] << "\n";
        }
        number_of_cachelines_size.flush();
        number_of_cachelines_size.close();
    }

    std::ofstream count_of_access_sizes_per_total_accessed_bytes_per_presence;
    count_of_access_sizes_per_total_accessed_bytes_per_presence.open(
        (output_dir + "/count_of_access_sizes_per_total_accessed_bytes_per_presence.csv"),
        std::ios::out);
    if (count_of_access_sizes_per_total_accessed_bytes_per_presence) {
        count_of_access_sizes_per_total_accessed_bytes_per_presence << "sep=;\n";
        for (size_t i = 0;
             i < count_accesses_of_accessed_bytes_per_total_accessed_bytes_of_cacheline
                     .size();
             ++i) {
            auto accesses =
                count_accesses_of_accessed_bytes_per_total_accessed_bytes_of_cacheline[i];
            count_of_access_sizes_per_total_accessed_bytes_per_presence << i << ";";
            for (size_t j = 1; j < accesses.size(); ++j) {
                if (i == 0) {
                    count_of_access_sizes_per_total_accessed_bytes_per_presence << j
                                                                                << ";";
                    continue;
                }
                count_of_access_sizes_per_total_accessed_bytes_per_presence << accesses[j]
                                                                            << ";";
            }
            count_of_access_sizes_per_total_accessed_bytes_per_presence << "\n";
        }
        count_of_access_sizes_per_total_accessed_bytes_per_presence.flush();
        count_of_access_sizes_per_total_accessed_bytes_per_presence.close();
    }

    // std::vector<std::vector<double>>
    //     accessed_bytes_per_access_per_total_accesses_per_cacheline(
    //         65, std::vector<double>(65, 0.0));
    //
    // for (int i = 0; i < 66; i++) {
    //     uint8_t total_bytes_accessed_per_cacheline_presence = 0;
    //     auto sub_accesses =
    //     count_of_access_sizes_per_bytes_accessed_per_cacheline[i]; for (int j = 0;
    //     j < 66; j++) {
    //         auto total_accesses_of_size_j = histogram[j];
    //         relative_accessed_bytes_per_total_accesses_per_cacheline[i][j] =
    //             ((long double)sub_accesses[j] / (long
    //             double)total_accesses_of_size_j);
    //     }
    // }
    //
    // std::cout << relative_accessed_bytes_per_total_accesses_per_cacheline[5][32]
    //           << std::endl;

    /*
        for (int i = 100; i < 10001;) {
            // auto moving_average_block_size =
            //     get_moving_average_block_size(basic_block_size_history, i);
            auto local_average_block_size =
                get_local_average_block_size(basic_block_size_history, i);
            std::cout << "#################### WINDOW SIZE " << i << "
       ####################"
                      << std::endl;
            std::cout << "#################### MOVING AVERAGE ####################"
                      << std::endl;
            // for (int j = 0; j < moving_average_block_size.size(); ++j) {
            //     std::cout << j << ": " << moving_average_block_size[j] << std::endl;
            // }
            std::cout << "#################### LOCAL AVERAGE ####################"
                      << std::endl;
            for (int j = 0; j < local_average_block_size.size(); ++j) {
                std::cout << j << ": " << local_average_block_size[j] << std::endl;
            }

            if (i < 1000)
                i += 100;
            else
                i += 1000;
        }*/

    // TODO: Extract to drawing functions
    // CTikz tikz_canvas;
    // CTikz tikz_global_bb;
    // CTikz tikz_presence_cacheline;
    // CTikz stacked_accesses;
    // try {
    //     std::vector<double> idx(65);
    //     std::iota(std::begin(idx), std::end(idx), 0);
    //     tikz_canvas.addData_vd(idx, relative_histogram,
    //                            "Bytes Accessed in Basic Block/Cacheline", "blue");
    //     tikz_global_bb.addData_vd(idx, overall_accessed_histgram, "#Bytes/Line
    //     (overall)",
    //                               "red");
    //     // tikz_presence_cacheline.addData_vd(idx, count_accessed_bytes_per_cacheline,
    //     //                                   "\#Bytes Accessed/Cacheline", "green");
    //     // tikz_canvas.createTikzPdfHist_vd(file, 64, 1, 64);
    //     // file += "norm.tikz";
    //     tikz_canvas.setXlabel_vd("\#Useful bytes in Cacheline");
    //     tikz_canvas.setYlabel_vd("Relative frequency");
    //     tikz_canvas.setLegendStyle_vd("draw=none");
    //     tikz_canvas.addAdditionalSettings_vd("legend pos={outer north east}");
    //     // tikz_presence_cacheline.setXlabel_vd("\#Accessed bytes in Cacheline");
    //     // tikz_presence_cacheline.setYlabel_vd("Count");
    //     // tikz_presence_cacheline.setLegendStyle_vd("draw=none");
    //     // tikz_presence_cacheline.addAdditionalSettings_vd("legend pos={outer north
    //     // east}");
    //
    // } catch (const CException &e) {
    //     std::cerr << e.what() << '\n';
    //     exit(-1);
    // }
    // bool created = false;
    // int counter = 0;
    // while (!created) {
    //     try {
    //         std::ostringstream oss;
    //         oss << output_dir << "numberofbytespercacheline_instr_"
    //             << handled_instructions << "." << counter << ".tikz";
    //         std::string file = oss.str();
    //         oss.str("");
    //         oss.clear();
    //         oss << output_dir << "numberofbytespercacheline_global_instr_"
    //             << handled_instructions << "." << counter << ".tikz";
    //         std::string global_file = oss.str();
    //         oss.str("");
    //         oss.clear();
    //         oss << output_dir << "numberofaccessedinstr_per_cacheline_"
    //             << handled_instructions << "." << counter << ".tikz";
    //         std::string cacheline_local_file = oss.str();
    //         tikz_canvas.createTikzPdf_vd(file);
    //         tikz_global_bb.createTikzPdf_vd(global_file);
    //         // tikz_presence_cacheline.createTikzPdf_vd(cacheline_local_file);
    //         created = true;
    //     } catch (const CException &e) {
    //         std::cerr << e.what() << '\n';
    //         counter++;
    //         if (counter == 100) {
    //             created = true;
    //             std::cout << "WARNING: COULDN'T CREATE GRAPH FILES" << std::endl;
    //         }
    //     }
    // }

    // std::cout << "Number of total cachelines: " << number_of_bytes_accessed.size()
    //           << std::endl;
    // std::cout << "Number of useful bytes : Number of cachelines" << std::endl;
    // for (size_t i = 0; i < histogram.size(); i++) {
    //     std::cout << i << ": " << histogram[i] << "\n";
    // }
}

/**
 * @brief Count how many occurances of the same number of byte accesses happened
 *
 * @param histogram The histogram vector containing the aggregated count
 * @param accesses A vector containing the history of cachelines accessed bit counts,
 * measured over one presence in L1-I
 *
 * @result The total of all accesses encoded in accesses.
 */
size_t
basic_block_stats_t::create_histogram_of_cachelineaccesses(
    std::vector<uint64_t> &histogram, std::vector<hitbytes> &accesses)
{
    size_t total_accesses = 0;
    for (auto it = accesses.begin(); it != accesses.end(); it++) {
        histogram[it->second] += it->first;
        total_accesses += it->first;
    }

    return total_accesses;
}

void
basic_block_stats_t::print_stats(std::string prefix)
{
    // trim_vector(count_per_basic_block_byte_size_);
    // trim_vector(count_per_basic_block_instr_size_);
    //
    // std::cout << "Max instr size: " << max_instr_size << std::endl;
    // std::cout << "Num basic blocks:" << basic_block_size_history.size() << std::endl;
    //
    // std::cout << prefix << "bb byte size:" << std::endl;
    // for (auto it = count_per_basic_block_byte_size_.begin()++;
    //     it != count_per_basic_block_byte_size_.end(); it++) {
    //    std::cout << "    " << it - count_per_basic_block_byte_size_.begin() << ": "
    //              << *it << std::endl;
    //}
    //
    // std::cout << prefix << "bb instr size:" << std::endl;
    // for (auto it = count_per_basic_block_instr_size_.begin();
    //     it != count_per_basic_block_instr_size_.end(); it++) {
    //    std::cout << "    " << it - count_per_basic_block_instr_size_.begin() << ": "
    //              << *it << std::endl;
    //}
    //
    print_bytes_accessed();
    /*
        for (auto it = number_of_bytes_accessed.begin(); it !=
       number_of_bytes_accessed.end(); it++) { const addr_t base_addr = (*it).first;
            std::vector<BasicBlock> vec = (*it).second;
            std::cout << base_addr << ": " << bytes_accessed(base_addr, vec) <<
       std::endl;
        }
    */

    std::cout << "Total L1-I Instructions: " << handled_instructions << std::endl;
    std::cout << "\tInterrupts:\t" << num_interrupts << std::endl;
    caching_device_stats_t::print_stats(prefix);
}

void
basic_block_stats_t::reset_current_cacheline_block()
{
    current_block_cacheline_constrained = {
        .starting_addr = 0, .end_addr = 0, .instr_size = 0, .byte_size = 0
    };
}

void
basic_block_stats_t::reset_current_cacheline_block(const memref_t &memref, bool hit)
{
    current_block_cacheline_constrained.miss = !hit;
    current_block_cacheline_constrained.starting_addr = memref.data.addr;
    current_block_cacheline_constrained.end_addr = memref.data.addr;
    current_block_cacheline_constrained.instr_size = 0;
    current_block_cacheline_constrained.byte_size = 0;
}

void
basic_block_stats_t::reset()
{
    caching_device_stats_t::reset();
    // TODO: Fixup missing variables
    count_per_basic_block_byte_size_.clear();
    count_per_basic_block_instr_size_.clear();
    eviction_to_fetch_index_map.clear();
    perfect_fetch_history.clear();
    current_block = {
        .starting_addr = 0, .end_addr = 0, .instr_size = 0, .byte_size = 0
    };
    reset_current_cacheline_block();
    basic_block_size_history.clear();
    basic_blocks_hit_count.clear();
    number_of_bytes_accessed.clear();
    handled_instructions = 0;
    num_interrupts = 0;
    bytes_accessed_per_presence_per_cacheline.clear();
    prev_hit = false;
}
