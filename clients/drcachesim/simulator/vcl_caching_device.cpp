#include "vcl_caching_device.h"
#include "cacheline_access_bitmask_helpers.h"
#include "../common/utils.h"
#include "snoop_filter.h"
#include <bits/basic_string.h>
#include <math.h>
#include <iterator>
#include "cache.h"
#include <filesystem>

vcl_caching_device_t::vcl_caching_device_t(std::string perfect_fetch_file)
    : I_caching_device_t()
{
    perfect_prefetching_path = perfect_fetch_file;
}

// TODO: Implement more sophisticated splitting behaviour here
std::set<std::pair<addr_t, addr_t>, AddrBlockCmp>
get_blocks_base_and_size(uint64_t mask, addr_t basic_block_addr)
{
    std::set<std::pair<addr_t, addr_t>, AddrBlockCmp> result;
    int prev_bit = 0;
    int curr_bit = 0;
    uint8_t start = 0;
    for (int i = 0; i < 64; ++i) {
        curr_bit = mask >> i & 1;
        if (curr_bit == 1 && prev_bit == 0) {
            start = i;
        } else if (curr_bit == 0 && prev_bit == 1) {
            result.insert(std::pair(basic_block_addr + start, basic_block_addr + i - 1));
        }
        prev_bit = curr_bit;
    }
    if (curr_bit == 1) {
        result.insert(std::pair(basic_block_addr + start, basic_block_addr + 63));
    }
    return result;
}

bool
vcl_caching_device_t::read_n_oracle_lines(size_t n)
{
    std::string line;
    bool endoffile = true;
    while (n > 0) {
        endoffile = std::getline(perfect_loading_decision_fh_, line).fail();
        if (endoffile) {
            return endoffile;
        }
        size_t idx = line.find(";");
        if (idx == std::string::npos) {
            std::cerr << "Unexpected input for perfect_loading_decision_file"
                      << std::endl;
            return true;
        }
        uint64_t base_address = std::stoull(line.substr(0, idx));
        line.erase(0, idx + 1);
        uint64_t mask = std::stoull(line);
        base_address_to_blocks_mapping[base_address].merge(
            get_blocks_base_and_size(mask, base_address));
        n--;
    }
    return endoffile;
}

std::pair<addr_t, addr_t>
get_candidate(const addr_t &address,
              const std::set<std::pair<addr_t, addr_t>, AddrBlockCmp> &blocks)
{
    // currently going for min enclosing block
    std::pair<addr_t, addr_t> candidate(0, 0);
    size_t min_block = 65;
    for (const auto &block : blocks) {
        if (block.first <= address && address <= block.second &&
            block.second - block.first < min_block) {
            candidate = block;
            min_block = block.second - block.first + 1;
        }
    }
    return candidate;
}

std::pair<int, int>
vcl_caching_device_t::start_and_end_oracle(addr_t address)
{
    addr_t base_addr = address & _CACHELINE_BASEADDRESS_MASK;
    auto add_to_block_mapping =
        base_address_to_blocks_mapping.end(); // we want to eventually read the whole file
    int64_t exp_backoff = 2;
    bool endoffile = false;
    auto candidate = std::pair<addr_t, addr_t>(0, 0);
    while (add_to_block_mapping == base_address_to_blocks_mapping.end() ||
           !(candidate.first <= address && address <= candidate.second)) {
        if (endoffile) {
            return std::pair<int, int>(
                0, 63); // we are at the end - not found - we should return 64
        }
        endoffile = read_n_oracle_lines(exp_backoff);
        add_to_block_mapping = base_address_to_blocks_mapping.find(base_addr);
        exp_backoff *= 1.5; // read larger shares if we miss something
        if (add_to_block_mapping != base_address_to_blocks_mapping.end())
            candidate = get_candidate(address, add_to_block_mapping->second);
    }

    return std::pair<int, int>(candidate.first - base_addr, candidate.second - base_addr);
}

vcl_caching_device_t::~vcl_caching_device_t()
{
    if (blocks_ == NULL) {
        return;
    }
    for (int i = 0; i < num_blocks_; i++) {
        delete blocks_[i];
    }
    delete[] blocks_;
}

bool
vcl_caching_device_t::vcl_enabled()
{
    return true;
}

void
vcl_caching_device_t::init_blocks()
{
    for (int i = 0; i < num_blocks_; i++) {
        blocks_[i] = new vcl_caching_device_block_t;
    }
}

bool
vcl_caching_device_t::init(int associativity, int block_size, int num_blocks,
                           I_caching_device_t *parent, caching_device_stats_t *stats,
                           prefetcher_t *prefetcher, bool inclusive, bool coherent_cache,
                           int id, snoop_filter_t *snoop_filter_,
                           const std::vector<I_caching_device_t *> &children)
{
    return false;
}

bool
vcl_caching_device_t::init(int associativity, std::vector<int> &way_sizes, int num_sets,
                           I_caching_device_t *parent, caching_device_stats_t *stats,
                           prefetcher_t *prefetcher, bool inclusive, bool coherent_cache,
                           int id, snoop_filter_t *snoop_filter,
                           const std::vector<I_caching_device_t *> &children)
{
    if (!IS_POWER_OF_2(num_sets))
        return false;
    if (stats == NULL)
        return false;
    else if (!*stats)
        return false;

    int total_size =
        num_sets * 8 * 64; // default so far - we tested with 8 way 64 byte lines

    std::filesystem::path root_dir(perfect_prefetching_path);
    std::filesystem::path size_dir(std::to_string(total_size));
    std::filesystem::path file("perfect_loading.dat");
    std::filesystem::path fullpath = root_dir / size_dir / file;
    perfect_loading_decision_fh_.open(fullpath.string(), std::ios::in | std::ios::binary);

    associativity_ = associativity;
    block_sizes_ = way_sizes;
    block_size_ = *(std::max_element(way_sizes.begin(), way_sizes.end()));
    num_sets_ = num_sets;
    num_blocks_ = num_sets * associativity;
    loaded_blocks_ = 0;
    set_idx_mask_ = num_sets_ - 1;
    assoc_bits_ = std::ceil(
        std::log2(associativity_)); // TODO: Why do we need power of 2 number of ways?
    block_sizes_bits_.reserve(block_sizes_.size());
    block_size_bits_ = compute_log2(block_size_);
    // For any reasonable cache this keeps us in a range between 1 and 10 for the buffer
    // size (L1-I, for data this can go up quite a bit higher, often ~ 20 buffer entries)
    int buffer_size = num_sets;
    fifo_buffer.init(std::max(1, buffer_size));

    for (const auto &block_size : block_sizes_) {
        block_sizes_bits_.push_back(std::ceil(std::log2(block_size)));
    }

    if (block_size_bits_ == -1 || assoc_bits_ == -1) {
        return false;
    }

    parent_ = parent;
    stats_ = stats;
    prefetcher_ = prefetcher;
    id_ = id;
    snoop_filter_ = snoop_filter;
    coherent_cache_ = coherent_cache;
    blocks_ = (caching_device_block_t **)new vcl_caching_device_block_t *[num_blocks_];

    init_blocks();
    // initialize block sizes
    for (int curr_set = 0; curr_set < num_sets_; curr_set++) {
        for (size_t curr_way = 0; curr_way < block_sizes_.size(); curr_way++) {
            ((vcl_caching_device_block_t *)
                 blocks_[curr_set * block_sizes_.size() + curr_way])
                ->size_ = block_sizes_[curr_way];
        }
    }
    last_tag_ = TAG_INVALID;
    inclusive_ = inclusive;
    children_ = children;

    return true;
}

void
vcl_caching_device_t::insert_cacheblock(vcl_caching_device_block_t *cache_block)
{
    addr_t victim_tag = cache_block->tag_;
    if (victim_tag == TAG_INVALID) {
        loaded_blocks_++;
        return;
    }
    if (!children_.empty() && inclusive_) {
        for (auto &child : children_) {
            child->invalidate(victim_tag, INVALIDATION_INCLUSIVE);
        }
    }
    if (!coherent_cache_) {
        return;
    }
    bool child_holds_tag = false;
    if (!children_.empty()) {
        for (auto &child : children_) {
            if (child->contains_tag(victim_tag)) {
                child_holds_tag = true;
                break;
            }
        }
    }
    if (child_holds_tag) {
        return;
    }
    if (snoop_filter_ != NULL) {
        snoop_filter_->snoop_eviction(victim_tag, id_);
    } else if (parent_ != NULL) {
        parent_->propagate_eviction(victim_tag, this);
    }
}

bool
vcl_caching_device_t::check_vcl_buffer(const memref_t &memref)
{
    auto base_address = memref.data.addr & _CACHELINE_BASEADDRESS_MASK;
    auto result = std::find_if(fifo_buffer.begin(), fifo_buffer.end(),
                               [base_address](const std::pair<addr_t, uint64_t> &val) {
                                   return val.first == base_address;
                               });
    if (result == fifo_buffer.end()) {
        return true;
    } else {
        // nullptr indicating hit in buffer
        uint8_t start = memref.data.addr - base_address;
        uint8_t end =
            start + memref.data.size; // points to the one after, but end is exclusive
        result->second = set_accessed(result->second, start, end); // [start, end[
        return false;
    }
}

void
vcl_caching_device_t::handle_miss(const memref_t &memref, const addr_t &tag)
{
    // TODO: put out into separate func - check which params are needed

    int block_idx = compute_block_idx(tag);
    auto base_address = memref.data.addr & _CACHELINE_BASEADDRESS_MASK;
    uint8_t start = memref.data.addr - base_address;
    uint8_t end =
        start + memref.data.size; // points to the one after, but end is exclusive
    std::pair<addr_t, uint64_t> *to_insert = nullptr;
    if (!fifo_buffer.empty()) {
        to_insert = &fifo_buffer.front();
    }

    // Fetch before inserting into buffer
    if (parent_ != NULL) {
        parent_->request(memref);
    }
    if (snoop_filter_ != NULL) {
        snoop_filter_->snoop(tag, id_, (memref.data.type == TRACE_TYPE_WRITE));
    }

    // TODO: Keep data if we eagerly evict
    auto bitmask = set_accessed(0, start, end);
    fifo_buffer.push(std::pair<addr_t, uint64_t>(base_address, bitmask));

    // we have a candidate to insert into the cache
    if (to_insert != nullptr && (*to_insert) != fifo_buffer.front()) {
        vcl_caching_device_block_t *cache_block = nullptr;
        auto blocks = get_start_end_of_bitmask(to_insert->second);
        // TODO: handle holes ...
        /* +1 for size, data.addr to adjust for overlapping */
        int size = blocks.back().second - blocks.front().first + 1;
        auto way = replace_which_way(block_idx, size);
        cache_block =
            (vcl_caching_device_block_t *)get_caching_device_block(block_idx, way);
        // we have something to insert...
        insert_cacheblock(cache_block);
        update_tag(cache_block, way, tag);
        uint8_t max_offset = 64 - cache_block->size_;
        // we will not start later in the cacheline than that
        cache_block->offset_ = std::min(start, max_offset);
        ((cache_t *)cache_)->access_update(block_idx, way);
    }
}

bool
vcl_caching_device_t::handle_candidate_blocks(
    std::vector<std::pair<caching_device_block_t *, int>> &block_ways,
    const memref_t &memref, vcl_caching_device_block_t **containing_block,
    const addr_t &tag)
{
    int block_idx = compute_block_idx(tag);
    auto base_address = memref.data.addr & _CACHELINE_BASEADDRESS_MASK;
    uint8_t start = memref.data.addr - base_address;
    uint8_t end =
        start + memref.data.size; // points to the one after, but end is exclusive

    for (auto &block_way : block_ways) {
        // We can ever hit in one - we ensure this by not inserting overlapping
        // regions
        // Update start and end if we hit in overlapping regions
        auto cache_block = (vcl_caching_device_block_t *)
                               block_way.first; // we only store vcl blocks in vcl cache
        auto way = block_way.second;
        uint8_t next_block = cache_block->size_ + cache_block->offset_;
        if (cache_block->validity_ && cache_block->offset_ <= start &&
            end <= next_block) {
            // [DONE] hit
            ((cache_t *)cache_)->access_update(block_idx, way);
            *containing_block = cache_block;
            return false;
        }
        // Old impossible implementation kept for reference (for now)
        // else if (cache_block->validity_ && cache_block->offset_ <= start &&
        //            start < next_block) {
        //     // we overlap partially - handle non-stored part [DONE]
        //     // NOTE: It is possible that we have a part of one instruction in
        //     one
        //     // cacheline and the other part in the next one - this is perfectly
        //     // fine for x86 and has been that way even for default cache sizes
        //     memref.data.addr =
        //         baseaddr + cache_block->offset_ + cache_block->size_;
        //     memref.data.size = memref.data.size -
        //         (cache_block->offset_ + cache_block->size_ - start);
        //     start = memref.data.addr - baseaddr;
        //     if (memref.data.size > MAX_X86_INSTR_SIZE) {
        //         std::cout << "WHAT THE FCK HAPPENED" << std::endl;
        //     }
        //     missed = true;
        //     partial_match = true;
        // } else if (cache_block->validity_ && cache_block->offset_ <= end &&
        //            end < next_block) {
        //
        //    // Handle overlapping in the end instead of the beginning. [DONE]
        //    memref.data.size -= (end - cache_block->offset_ + 1);
        //    partial_match = true;
        //
        //} else if (cache_block->validity_ && start <= cache_block->offset_ &&
        //           (next_block - 1) <= end) {
        //    // Handle case where cache block contains proper subset of memref
        //    TODO std::cout
        //        << "CACHEBLOCK RIPPING HOLE INTO MEMREF -- DOES THIS EVEN
        //        HAPPEN"
        //        << std::endl;
        //}
    }

    // We are not contained - check vcl buffer
    return check_vcl_buffer(memref);
}

void
vcl_caching_device_t::request(_memref_t const &memref_in)
{
    memref_t memref;
    addr_t final_addr = memref_in.data.addr + memref_in.data.size - 1;
    addr_t final_tag = compute_tag(final_addr);
    addr_t tag = compute_tag(memref_in.data.addr);
    // Note: the tag still contains the index bits here
    // if (tag < final_tag) {
    //     std::cout << "tag: " << tag << "; final tag: " << final_tag << std::endl;
    // }

    memref = memref_in;
    for (; tag <= final_tag; ++tag) {
        bool missed = false;
        vcl_caching_device_block_t *cache_block = nullptr;

        if (tag + 1 <= final_tag) {
            memref.data.size = ((tag + 1) << block_size_bits_) -
                memref.data.addr; // is this good enough?
        }

        if (memref.data.size > MAX_X86_INSTR_SIZE) {
            std::cout << "WTF" << std::endl;
        }

        auto block_ways = find_all_caching_device_blocks(memref.data.addr, false, true);
        if (block_ways.size() > 0) {
            missed = handle_candidate_blocks(block_ways, memref, &cache_block, tag);
        } else {
            // check buffer
            missed = check_vcl_buffer(memref);
        }

        if (missed) {
            handle_miss(memref, tag);
        }

        if (missed && !type_is_prefetch(memref.data.type) && prefetcher_ != nullptr) {
            prefetcher_->prefetch(this, memref);
        }

        if (tag + 1 <= final_tag) {
            addr_t next_addr = (tag + 1) << block_size_bits_;
            memref.data.addr = next_addr;
            memref.data.size = final_addr - next_addr + 1;
        }
        if (memref.data.size > MAX_X86_INSTR_SIZE) {
            std::cout << "WTF" << std::endl;
        }
        record_access_stats(memref, !missed, cache_block);
    }
}

void
vcl_caching_device_t::invalidate(addr_t tag, invalidation_type_t invalidation_type_)
{
    auto all_affected_blocks_way =
        find_all_caching_device_blocks((tag << block_size_bits_), false, false);
    if (all_affected_blocks_way.size() > 0) {
        for (auto &block_way : all_affected_blocks_way) {
            invalidate_caching_device_block(block_way.first);
            stats_->invalidate(
                invalidation_type_); // should we count this per invalidation request or
                                     // per invalidated block?
        }
        if (invalidation_type_ == INVALIDATION_INCLUSIVE && inclusive_ &&
            !children_.empty()) {
            for (auto &child : children_) {
                child->invalidate(tag, invalidation_type_);
            }
        }
    }
}

bool
vcl_caching_device_t::contains_tag(addr_t tag)
{
    throw std::runtime_error("contains_tag not yet implemented");
    return false;
}

void
vcl_caching_device_t::propagate_eviction(addr_t tag, const I_caching_device_t *requester)
{
    throw std::runtime_error("propagate_eviction not yet implemented");
}

void
vcl_caching_device_t::propagate_write(addr_t tag, const I_caching_device_t *requester)
{
    throw std::runtime_error("propagate_write not yet implemented");
}

int
vcl_caching_device_t::replace_which_way(int block_idx)
{
    return replace_which_way(
        block_idx, 64); // We behave like a smaller cache but the same in terms of storage
}

int
vcl_caching_device_t::replace_which_way(int block_idx, int size)
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
    int first_viable_way = 0;
    first_viable_way = get_smallest_possible_way(size);

    // For now we simply do lru over the large enough entries
    int max_way = 0;
    int max_counter = 0;
    // double best_weighted_relation = 0; // This could be an optimization that takes
    // wasted bytes into account
    // We assume that the ways are sorted by size, such that we
    // can assume that the smalles one fitting comes first. std::vector<std::pair<int,
    // int>>
    //     candidate_sizes_and_idx; // We use it as <size, idx>. All INVALID lines are
    //                              // candidates + LRU candidate
    for (int way = first_viable_way; way < associativity_; ++way) {
        if (get_caching_device_block(block_idx, way)->tag_ == TAG_INVALID) {
            max_way = way; // we found an empty way - use this
            break;
        }
        if (way == first_viable_way ||
            get_caching_device_block(block_idx, way)->counter_ > max_counter) {
            max_counter = get_caching_device_block(block_idx, way)->counter_;
            max_way = way;
        }
    }
    return max_way;
#pragma GCC diagnostic pop
}

int
vcl_caching_device_t::get_next_way_to_replace(const int block_idx) const
{
    throw std::runtime_error("get_next_way_to_replace not yet implemented");
    return -1;
}

// TODO: Type double that function for the vcl blocks r unify block interface
void
vcl_caching_device_t::record_access_stats(const _memref_t &memref, bool hit,
                                          caching_device_block_t *cache_block)
{
    stats_->access(memref, hit, cache_block);

    if (hit) {
        for (I_caching_device_t *up = parent_; up != nullptr; up = up->parent_) {
            up->stats_->child_access(memref, hit, cache_block);
        }
    } else if (parent_ != nullptr) {
        parent_->stats_->child_access(memref, hit, cache_block);
    }
}
std::vector<std::pair<caching_device_block_t *, int>>
vcl_caching_device_t::find_all_caching_device_blocks(addr_t addr, bool only_one,
                                                     bool check_inlier)
{
    addr_t tag = compute_tag(addr);
    addr_t blockidx = compute_block_idx(tag);
    std::vector<std::pair<caching_device_block_t *, int>> candidate_blocks;
    for (int way = 0; way < associativity_; ++way) {
        vcl_caching_device_block_t *block =
            (vcl_caching_device_block_t *)get_caching_device_block(blockidx, way);
        // TODO: Check if < or <=
        addr_t startaddr = (addr & _CACHELINE_BASEADDRESS_MASK) + block->offset_;
        addr_t endaddr = startaddr + block->size_;
        if (block->tag_ == tag &&
            ((startaddr < addr && addr <= endaddr) || !check_inlier)) {
            candidate_blocks.push_back(std::make_pair(block, way));
            if (only_one && check_inlier)
                return candidate_blocks;
        }
    }
    return candidate_blocks; // its a miss if empty
}

std::pair<caching_device_block_t *, int>
vcl_caching_device_t::find_caching_device_block(addr_t addr)
{
    auto result = find_all_caching_device_blocks(addr, true);
    if (result.size() == 0) {
        return std::make_pair(nullptr, 0);
    }
    return result[0]; // can at most have len 1
}
