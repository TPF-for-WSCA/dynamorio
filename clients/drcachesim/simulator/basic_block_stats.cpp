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
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <numeric>

#include "basic_block_stats.h"
// DEBUGGING IMPORTS: HOW MUCH MEMORY DO WE USE
#include "stdlib.h"
#include "string.h"

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
{
    basic_block_size_history.reserve(10000);
    basic_blocks_hit_count.reserve(10000);
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

    if ((memref.data.addr - current_block.end_addr) <= MAX_X86_INSTR_SIZE) {
        // x86 instrs can take up a max of 15 bytes, hence we assume a single basic
        // block if the address does not differ by more than this
        current_block.instr_size++;
        current_block.byte_size += memref.data.size;
        current_block.end_addr = memref.data.addr;
    } else {
        // We moved more than the size of any allowed x86 instruction - we must have been
        // interrupted
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
    return is_interrupt_;
}

void
basic_block_stats_t::track_cacheline_access(const memref_t &memref)
{
    addr_t cacheline_start_address =
        cache_block_address_mask & current_block.starting_addr;
    addr_t cacheline_end_address = cache_block_address_mask &
        (current_block.starting_addr + current_block.byte_size);
    cacheline_end_address = (cacheline_end_address == cacheline_start_address)
        ? cacheline_end_address + 64
        : cacheline_end_address;

    if (cacheline_end_address - cacheline_start_address > max_cacheline_bb) {
        max_cacheline_bb = cacheline_end_address - cacheline_start_address;
        std::cout << "Current largest overlap of a basic block: " << max_cacheline_bb / 64
                  << " lines" << std::endl;
    }

    for (; cacheline_start_address < cacheline_end_address;
         cacheline_start_address += 64) {
        auto basic_block_vec = number_of_bytes_accessed[cacheline_start_address];
        // if it spans more than one line it might be that the vec is still empty
        if (current_block.miss || basic_block_vec.size() == 0) {
            basic_block_vec.push_back(current_block);
        } else {
            auto is_same_block = [this](const BasicBlock &bb) {
                return (bb.starting_addr <= current_block.starting_addr) &&
                    (bb.end_addr >= current_block.end_addr);
            };
            auto bb_entry = std::find_if(basic_block_vec.rbegin(), basic_block_vec.rend(),
                                         is_same_block);

            if (bb_entry == basic_block_vec.rend()) {
                current_block.access_sizes[current_block.byte_size]++;
                basic_block_vec.push_back(current_block);
            } else {
                bb_entry->access_sizes[current_block.byte_size]++;
                bb_entry->hits += 1;
            }
        }
        number_of_bytes_accessed[cacheline_start_address] = basic_block_vec;
    }
}

void
basic_block_stats_t::handle_branch(const memref_t &memref)
{
    // assuming x86, instructions in basic blocks have monotoneously increasing
    // addresses
    if (current_block.starting_addr <= current_block.end_addr &&
        current_block.starting_addr != 0) {
        record_block(current_block);
    } else {
        printf("WARN: ADDRESS MISMATCH, BLOCK END IS SMALLER THAN BLOCK START\n");
    }
    insert_bbcount(count_per_basic_block_byte_size_, current_block.byte_size);
    insert_bbcount(count_per_basic_block_instr_size_, current_block.instr_size);

    basic_blocks_hit_count[current_block] += 1;
    basic_block_size_history.push_back(current_block.byte_size);

    track_cacheline_access(memref);

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
    current_block.access_sizes = std::vector<size_t>(65, 0);
}

void
basic_block_stats_t::handle_interrupt(const memref_t &memref, bool hit)
{
    current_block.starting_addr = memref.data.addr;
    current_block.end_addr = memref.data.addr;
    current_block.byte_size = memref.data.size;
    current_block.instr_size = 1;
    current_block.miss = !hit;
}

void
basic_block_stats_t::access(const memref_t &memref, bool hit,
                            caching_device_block_t *cache_block)
{
    caching_device_stats_t::access(memref, hit, cache_block);

    // size_t current_ram_consumption = getPhysicalRAMUsageValue() >> 10;
    // if (max_memory_consumption < current_ram_consumption) {
    //     std::cout << "CURRENT RAM: " << current_ram_consumption << " MB" << std::endl;
    // }

    if (type_is_prefetch(memref.data.type)) {
        std::cout << "YES SAW A PREFETCH" << std::endl;
        return; // we only care about basic blocks
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

    bool is_interrupt_ = false;

    // count number of instructions in basic block
    if (is_instr_) {
        is_interrupt_ = handle_instr(memref, hit);
        is_branch_ = is_branch_ || is_interrupt_;
    }

    // calculate number of bytes in instruction based on addresses.
    if (is_branch_) {
        handle_branch(memref);
    }

    if (is_interrupt_) {
        handle_interrupt(memref, hit);
    }
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

void
basic_block_stats_t::record_block(const BasicBlock &bb)
{
    auto target_bb_node = basic_blocks_hit_count.extract(bb);
    if (target_bb_node.empty()) {
        basic_blocks_hit_count.insert({ bb, 1 });
        return;
    }

    // key and value both return references
    auto target_bb = target_bb_node.key();
    BasicBlock merged_bb = target_bb;
    if (bb.starting_addr < target_bb.starting_addr) {
        merged_bb.starting_addr = bb.starting_addr;
    }

    if (bb.end_addr > target_bb.end_addr) {
        printf("WARN: BLOCK ENDS EXTENDS - SHOULD NEVER HAPPEN AS WE ALWAYS FINISH BB AT "
               "BRANCHING INSTR\n");
        target_bb.end_addr = bb.end_addr;
    }

    target_bb_node.mapped() += 1;

    basic_blocks_hit_count.insert(move(target_bb_node));
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
    uint64_t bitmask = 1;
    for (int i = 0; i < lower; i++) {
        bitmask = bitmask << 1;
    }

    for (int i = lower; i < upper; i++) {
        mask |= bitmask;
        bitmask = bitmask << 1;
    }
    return mask;
}

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

std::pair<uint8_t, std::vector<uint8_t>>
basic_block_stats_t::bytes_accessed(
    const addr_t &cacheline_base, std::vector<BasicBlock> &blocks_contained,
    std::vector<std::vector<uint64_t>>
        count_of_access_sizes_per_bytes_accessed_per_cacheline)
{
    uint64_t mask = 0;
    std::vector<uint8_t> accesses;
    accesses.reserve(blocks_contained.size());
    for (auto it = blocks_contained.begin(); it != blocks_contained.end(); it++) {
        auto tmp_mask = bytes_accessed_by_block(cacheline_base, *it);
        auto tmp_count = (uint8_t)__builtin_popcountll(tmp_mask);
        for (auto const &accessed_bytes : it->access_sizes) {
            count_of_access_sizes_per_bytes_accessed_per_cacheline[accessed_bytes]
                                                                  [tmp_count] += 1;
        }
        accesses.push_back(tmp_count);
        mask |= tmp_mask;
    }

    return std::pair<uint8_t, std::vector<uint8_t>> { (uint8_t)__builtin_popcountll(mask),
                                                      accesses };
}

void
basic_block_stats_t::print_bytes_accessed()
{
    std::vector<uint64_t> histogram(65, 0);
    std::vector<double> overall_accessed_histgram(65, 0);
    std::vector<std::vector<uint64_t>>
        count_of_access_sizes_per_bytes_accessed_per_cacheline(
            65, std::vector<uint64_t>(65, 0));
    uint64_t total_allocations = 0;

    for (auto it = number_of_bytes_accessed.begin(); it != number_of_bytes_accessed.end();
         it++) {
        std::vector<BasicBlock> vec = (*it).second;
        total_allocations += vec.size();
        const addr_t base_addr = (*it).first;
        auto accessed = bytes_accessed(
            base_addr, vec, count_of_access_sizes_per_bytes_accessed_per_cacheline);
        create_histogram_of_cachelineaccesses(histogram, accessed.second);
        int total_accessed = (int)accessed.first;
        overall_accessed_histgram[total_accessed]++;
        std::cout << base_addr << ": " << total_accessed << std::endl;
    }

    std::vector<double> relative_histogram;
    relative_histogram.reserve(65);

    for (auto const &bucket : histogram) {
        relative_histogram.push_back((long double)bucket /
                                     (long double)total_allocations);
    }

    std::vector<std::vector<double>>
        relative_accessed_bytes_per_total_accesses_per_cacheline(
            65, std::vector<double>(65, 0.0));

    for (int i = 0; i < 66; i++) {
        auto sub_accesses = relative_accessed_bytes_per_total_accesses_per_cacheline[i];
        for (int j = 0; j < 66; j++) {
            auto total_accesses_of_size_j = histogram[j];
            relative_accessed_bytes_per_total_accesses_per_cacheline[i][j] =
                ((long double)sub_accesses[j] / (long double)total_accesses_of_size_j);
        }
    }

    std::cout << relative_accessed_bytes_per_total_accesses_per_cacheline[5][32]
              << std::endl;

    try {

        CTikz tikz_canvas;
        CTikz tikz_global_bb;
        CTikz stacked_accesses;
        std::vector<double> idx(65);
        std::iota(std::begin(idx), std::end(idx), 0);
        tikz_canvas.addData_vd(idx, relative_histogram, "Bytes Accessed/Cacheline",
                               "blue");
        tikz_global_bb.addData_vd(idx, overall_accessed_histgram, "#Bytes/Line (overall)",
                                  "red");
        std::string file = output_dir + "numberofbytespercacheline.tikz";
        std::string global_file = output_dir + "numberofbytespercacheline_global.tikz";
        // tikz_canvas.createTikzPdfHist_vd(file, 64, 1, 64);
        // file += "norm.tikz";
        tikz_canvas.setXlabel_vd("\\#Useful bytes in Cacheline");
        tikz_canvas.setYlabel_vd("Relative frequency");
        tikz_canvas.setLegendStyle_vd("draw=none");
        tikz_canvas.addAdditionalSettings_vd("legend pos={outer north east}");
        tikz_canvas.createTikzPdf_vd(file);
        tikz_global_bb.createTikzPdf_vd(global_file);
    } catch (const CException &e) {
        std::cerr << e.what() << '\n';
        exit(-1);
    }

    std::cout << "Number of total cachelines: " << number_of_bytes_accessed.size()
              << std::endl;
    std::cout << "Number of useful bytes : Number of cachelines" << std::endl;
    for (size_t i = 0; i < histogram.size(); i++) {
        std::cout << i << ": " << histogram[i] << std::endl;
    }
}

/**
 * @brief Count how many occurances of the same number of byte accesses happened
 *
 * @param histogram The histogram vector containing the aggregated count
 * @param accesses A vector containing the history of cachelines accessed bit counts,
 * measured over one presence in L1-I
 */
void
basic_block_stats_t::create_histogram_of_cachelineaccesses(
    std::vector<uint64_t> &histogram, std::vector<uint8_t> &accesses)
{
    for (auto it = accesses.begin(); it != accesses.end(); it++) {
        histogram[*it] += 1;
    }
}

void
basic_block_stats_t::print_stats(std::string prefix)
{
    trim_vector(count_per_basic_block_byte_size_);
    trim_vector(count_per_basic_block_instr_size_);

    std::cout << "Max instr size: " << max_instr_size << std::endl;
    std::cout << "Num instructions:" << basic_block_size_history.size() << std::endl;

    std::cout << prefix << "bb byte size:" << std::endl;
    for (auto it = count_per_basic_block_byte_size_.begin()++;
         it != count_per_basic_block_byte_size_.end(); it++) {
        std::cout << prefix << "    " << it - count_per_basic_block_byte_size_.begin()
                  << ": " << *it << std::endl;
    }

    std::cout << prefix << "bb instr size:" << std::endl;
    for (auto it = count_per_basic_block_instr_size_.begin();
         it != count_per_basic_block_instr_size_.end(); it++) {
        std::cout << prefix << "    " << it - count_per_basic_block_instr_size_.begin()
                  << ": " << *it << std::endl;
    }

    print_bytes_accessed();
    /*
        for (auto it = number_of_bytes_accessed.begin(); it !=
       number_of_bytes_accessed.end(); it++) { const addr_t base_addr = (*it).first;
            std::vector<BasicBlock> vec = (*it).second;
            std::cout << base_addr << ": " << bytes_accessed(base_addr, vec) <<
       std::endl;
        }
    */
    caching_device_stats_t::print_stats(prefix);
}

void
basic_block_stats_t::reset()
{
    caching_device_stats_t::reset();
    // TODO: Fixup missing variables
    count_per_basic_block_byte_size_.clear();
    count_per_basic_block_instr_size_.clear();
    current_block = { .starting_addr = 0,
                      .end_addr = 0,
                      .instr_size = 0,
                      .byte_size = 0,
                      .access_sizes = std::vector<size_t>(65, 0) };
    basic_block_size_history.clear();
    basic_blocks_hit_count.clear();
    number_of_bytes_accessed.clear();
}

// void
// basic_block_stats_t::flush(const memref_t &memref)
//{
//
//}