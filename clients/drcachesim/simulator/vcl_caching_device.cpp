#include "vcl_caching_device.h"
#include "../common/utils.h"
#include "snoop_filter.h"
#include <bits/basic_string.h>
#include <math.h>
#include <iterator>

vcl_caching_device_t::vcl_caching_device_t(std::string perfect_fetch_file)
    : I_caching_device_t()
{
    perfect_loading_decision_fh_.open(perfect_fetch_file,
                                      std::ios::in | std::ios::binary);
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
        endoffile = !std::getline(perfect_loading_decision_fh_, line).fail();
        if (endoffile) {
            return endoffile;
        }
        size_t idx = line.find(";");
        if (idx == std::string::npos) {
            std::cerr << "Unexpected input for perfect_loading_decision_file"
                      << std::endl;
            return false;
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

    associativity_ = associativity;
    block_sizes_ = way_sizes;
    block_size_ = *(std::max_element(way_sizes.begin(), way_sizes.end()));
    sets_in_cache_ = num_sets;
    num_blocks_ = num_sets * associativity;
    loaded_blocks_ = 0;
    set_idx_mask_ = sets_in_cache_ - 1;
    assoc_bits_ = std::ceil(
        std::log2(associativity_)); // TODO: Why do we need power of 2 number of ways?
    block_sizes_bits_.reserve(block_sizes_.size());
    block_size_bits_ = compute_log2(block_size_);

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
    int way_count = 1;
    for (const auto &block_size : block_sizes_) {
        for (int i = (way_count - 1) * sets_in_cache_; i < way_count * sets_in_cache_;
             i++) {
            ((vcl_caching_device_block_t *)blocks_[i])->size_ = block_size;
        }
        way_count += 1;
    }
    last_tag_ = TAG_INVALID;
    inclusive_ = inclusive;
    children_ = children;

    return true;
}

void
vcl_caching_device_t::request(_memref_t const &memref_in)
{
    memref_t memref;
    addr_t final_addr = memref_in.data.addr + memref_in.data.size - 1;
    addr_t final_tag = compute_tag(final_addr);
    addr_t tag = compute_tag(memref_in.data.addr);

    memref = memref_in;
    for (; tag <= final_tag; ++tag) {
        int way = associativity_;
        int block_idx = compute_block_idx(tag);

        bool missed = false;
        if (tag + 1 <= final_tag) {
            memref.data.size = ((tag + 1) << block_size_bits_) -
                memref.data.addr; // is this good enough?
        }

        auto block_ways = find_all_caching_device_blocks(memref.data.addr);
        if (block_ways.size() > 0) {
            // found - check if in range
            for (auto &block_way : block_ways) {
                // We can ever hit in one - we ensure this by not inserting overlapping
                // regions
                vcl_caching_device_block_t *cache_block =
                    (vcl_caching_device_block_t *)
                        block_way.first; // we only store vcl blocks in vcl cache
                way = block_way.second;
                addr_t baseaddr = memref.data.addr & _CACHELINE_BASEADDRESS_MASK;
                uint8_t start = memref.data.addr - baseaddr;
                uint8_t end = start + memref.data.size;
                if (cache_block->validity_ && cache_block->offset_ <= start &&
                    end <= (cache_block->size_ + cache_block->offset_)) {
                    // hit
                } else if (cache_block->validity_ && cache_block->offset_ <= start &&
                           start <= (cache_block->size_ + cache_block->offset_)) {
                    // we overlap partially - handle non-stored part
                    memref.data.addr =
                        baseaddr + cache_block->offset_ + cache_block->size_;
                    memref.data.size = memref.data.size -
                        (cache_block->offset_ + cache_block->size_ - start);
                    missed = true;
                } else {
                    // we do not overlap - handle a miss
                    missed = true;
                }
            }
        } else {
            // miss
            missed = true;
        }

        if (missed) {
            // TODO: put out into separate func - check which params are needed
            auto predicted_block = start_and_end_oracle(memref.data.addr);
            way = replace_which_way(block_idx,
                                    /* +1 for size */
                                    predicted_block.second - predicted_block.first + 1);
            vcl_caching_device_block_t *cache_block =
                (vcl_caching_device_block_t *)&get_caching_device_block(block_idx, way);
            record_access_stats(memref, false, cache_block);
            if (parent_ != NULL) {
                parent_->request(memref);
            }
            if (snoop_filter_ != NULL) {
                snoop_filter_->snoop(tag, id_, (memref.data.type == TRACE_TYPE_WRITE));
            }

            addr_t victim_tag = cache_block->tag_;
            if (victim_tag == TAG_INVALID) {
                loaded_blocks_++;
            } else {
                if (!children_.empty() && inclusive_) {
                    for (auto &child : children_) {
                        child->invalidate(victim_tag, INVALIDATION_INCLUSIVE);
                    }
                }
                if (coherent_cache_) {
                    bool child_holds_tag = false;
                    if (!children_.empty()) {
                        for (auto &child : children_) {
                            if (child->contains_tag(victim_tag)) {
                                child_holds_tag = true;
                                break;
                            }
                        }
                    }
                    if (!child_holds_tag) {
                        if (snoop_filter_ != NULL) {
                            snoop_filter_->snoop_eviction(victim_tag, id_);
                        } else if (parent_ != NULL) {
                            parent_->propagate_eviction(victim_tag, this);
                        }
                    }
                }
            }
            update_tag(cache_block, way, tag);
        }

        access_update(block_idx, way);

        if (missed && !type_is_prefetch(memref.data.type) && prefetcher_ != nullptr) {
            prefetcher_->prefetch(this, memref);
        }

        if (tag + 1 <= final_tag) {
            addr_t next_addr = (tag + 1) << block_size_bits_;
            memref.data.addr = next_addr;
            memref.data.size = final_addr - next_addr + 1;
        }
    }
}

void
vcl_caching_device_t::invalidate(addr_t tag, invalidation_type_t invalidation_type_)
{
}

bool
vcl_caching_device_t::contains_tag(addr_t tag)
{
    return false;
}

void
vcl_caching_device_t::propagate_eviction(addr_t tag, const I_caching_device_t *requester)
{
}

void
vcl_caching_device_t::propagate_write(addr_t tag, const I_caching_device_t *requester)
{
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
    int first_viable_way = get_smallest_possible_way(size);

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
        if (get_caching_device_block(block_idx, way).tag_ == TAG_INVALID) {
            max_way = way; // we found an empty way - use this
            break;
        }
        if (way == first_viable_way ||
            get_caching_device_block(block_idx, way).counter_ > max_counter) {
            max_counter = get_caching_device_block(block_idx, way).counter_;
            max_way = way;
        }
    }
    return max_way;
}

int
vcl_caching_device_t::get_next_way_to_replace(const int block_idx) const
{
    return -1;
}

// TODO: Type double that function for the vcl blocks r unify block interface
void
vcl_caching_device_t::record_access_stats(const _memref_t &memref, bool hit,
                                          caching_device_block_t *cache_block)
{
}
std::vector<std::pair<caching_device_block_t *, int>>
vcl_caching_device_t::find_all_caching_device_blocks(addr_t addr, bool only_one)
{
    addr_t tag = compute_tag(addr);
    addr_t blockidx = compute_block_idx(tag);
    std::vector<std::pair<caching_device_block_t *, int>> candidate_blocks;
    for (int way = 0; way < associativity_; ++way) {
        vcl_caching_device_block_t *block =
            (vcl_caching_device_block_t *)&get_caching_device_block(blockidx, way);
        // TODO: Check if < or <=
        addr_t startaddr = (addr & _CACHELINE_BASEADDRESS_MASK) + block->offset_;
        addr_t endaddr = startaddr + block->size_;
        if (block->tag_ == tag && startaddr < addr && addr <= endaddr) {
            candidate_blocks.push_back(std::make_pair(block, way));
            if (only_one)
                return candidate_blocks; // we only get here once before returning
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
