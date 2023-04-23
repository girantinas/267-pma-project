#pragma once

#include <upcxx/upcxx.hpp>
// #include "butil.hpp"
#include <algorithm>
#include <climits>
#include <cmath>
#include <chrono>
#include <iostream>
#include <vector>
#include <deque>
#include <utility>
#include <functional>
#include <unordered_set>
#include "set_pma.hpp"

using namespace std;

class DistPCSR {
    public:
        typedef dist_range_t = upcxx::dist_object<vector<range_t>>;
        typedef range_t = tuple<chrono::milliseconds, uint64_t, uint64_t>
        DistPCSR(uint32_t size, uint64_t max_val);
        uint32_t size();
        uint32_t num_vertices; // Restrict this to power of 2
        /* Inserts an edge from one vertex to another. Fails if either vertex does not exist. */
        void insert_edge_local(uint32_t from, uint32_t to);
        /* Queries if an edge exists. */
        bool query_edge(uint32_t from, uint32_t to);
        /* Edge list of a vertex. If the vertex does not exist, the list is empty. */
        vector<uint32_t> edges(uint32_t from);
        /* Returns the entire adjacency list representation of the graph. */
        vector<pair<uint32_t, vector<uint32_t>>> adjacency_lists();
    private:
        static constexpr double INIT_LEAF_MAX = 0.2;
        static uint64_t make_edge_tuple(uint32_t from, uint32_t to);
        static pair<uint32_t, uint32_t> get_edge_tuple(uint64_t edge);
        uint32_t target_rank(uint32_t from);
        SetPMA spma;
        const uint32_t TO_ONES = 0xFFFFFFFF;
        void acquire_redis_lock();
        void release_redis_lock();

        pair<uint64_t, uint64_t> my_range;
        dist_range_t cached_ranges;

        bool redistributing;
        upcxx::global_ptr<int64_t> redis_lock; 
        upcxx::atomic_domain<int64_t> ad;

        struct Insert {
            uint32_t edge_from;
            uint32_t edge_to;
        }

        deque<Insert> rq_queue;
        vector<upcxx::team> teams; // teams[0] = you + neighbor, teams[1] = you + 3 adjacent ranks, etc. Only 2p teams in total, should be manageable

};

uint32_t MSSB(uint32_t x) {
    uint32_t i = 0;
    while (x != 0) {
        x = x >> 1;
        ++i;
    }
    return i - 1;
}

DistPCSR::DistPCSR(uint32_t size, uint64_t max_val) : spma(SetPMA(size, INIT_LEAF_MAX)) {
    num_vertices = max_val;
    auto lower_from = (max_val / upcxx::rank_n()) * upcxx::rank_me();
    auto upper_from = (max_val / upcxx::rank_n()) * (upcxx::rank_me() + 1);
    auto lower_tup = make_edge_tuple(lower_from, 0);
    auto upper_tup = make_edge_tuple(upper_from, TO_ONES);
    my_range = make_pair(lower_tup, upper_tup);
    vector<pair<uint32_t, uint32_t>> ranges;
    for (int i = 0; i < upcxx::rank_n(); i += 1) {
        auto lower = (max_val / upcxx::rank_n()) * i;
        auto upper = (max_val / upcxx::rank_n()) * (i + 1);
        auto lower_tup = make_edge_tuple(lower_from, 0);
        auto upper_tup = make_edge_tuple(upper_from, TO_ONES);
        ranges.push_back(make_pair(chrono::milliseconds::zero(), lower_tup, upper_tup));
    }
    cached_from_ranges = dist_range_t(ranges);

    redistributing = false;
    if (upcxx::rank_me() == 0) {
        redis_lock = new_<int64_t>(0);
    }
    redis_lock = upcxx::broadcast(redis_lock, 0).wait();
    ad = upcxx::atomic_domain<int64_t>({upcxx::atomic_op::compare_exchange, upcxx::atomic_op::store});
    
    uint32_t n_procs_in_team = 1;
    uint32_t proc = upcxx::rank_me();
    uint32_t total_procs = upcxx::rank_n();
    upcxx::team& world_team = upcxx::world();
    for (int level = 1; level < MSSB(upcxx::rank_n()); level += 1) {
        n_procs_in_team *= 2;
        uint32_t color = total_procs / n_procs;
        uint32_t key = total_procs % n_procs;
        upcxx::team new_team = world_team.split(color, key);
        teams.push_back(new_team);
    }
}

void DistPCSR::acquire_redis_lock() {
    while (ad.compare_exchange(redis_lock, 0, 1).wait());
}

void DistPCSR::release_redis_lock() {
    ad.store(redis_lock, 0).wait();
}

uint64_t DistPCSR::make_edge_tuple(uint32_t from, uint32_t to) { return (((uint64_t) from) << 32) | to; }
pair<uint32_t, uint32_t> DistPCSR::get_edge_tuple(uint64_t edge) { return make_pair((uint32_t)(edge >> 32), (uint32_t)((edge << 32) >> 32)); } /* doesn't work */

// binary search in our cached values of ranges; could be outdated
uint32_t DistPCSR::target_rank(uint32_t from, uint32_t to) {
    uint64_t edge_tuple = make_edge_tuple(from, to);

    // first check local rank, will always be up to date.
    if (edge_tuple >= my_range.first && edge_tuple < my_range.second) {
        return upcxx::rank_me();
    }

    auto x = range_t(0, edge_tuple, 0); // dummy object to ensure types match
    auto target = std::lower_bound(cached_ranges->begin(), cached_ranges->end(), edge_tuple, [](const range_t& a, const range_t& b) { return std::get<1>(a) < std::get<1>(b); });
    auto index = target - cached_ranges->begin();
    return index;
}

void flush_queue(upcxx::dist_object<DistPCSR>& pcsr) {
    if (pcsr->redistributing) {
        return;
    }
    while (pcsr->rq_queue.size() != 0) {
        auto insert = pcsr->rq_queue.front();
        pcsr->rq_queue.pop_front();
        pcsr->insert_edge_local(insert.from, insert.to);        
    }
}

void DistPCSR::insert_edge_local(uint32_t from, uint32_t to) {
    if (spma._num_elements / spma.size() >= spma._leaf_max) {
        redistribute();
    }
    if (target_rank(from, to) != upcxx::rank_me()) {
        insert_edge(pcsr, from, to);
    }
    uint64_t edge = make_edge_tuple(from, to);
    spma.insert(edge);
}

upcxx::future<> insert_edge(upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, uint32_t to) {
    uint32_t rank = target_rank(from, to);
    
    if (rank == upcxx::rank_me()) {
        pcsr->rq_queue.push_back(Insert(from, to));
        return upcxx::make_future();
    }
    
    return upcxx::rpc(rank,
        [](upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, uint32_t to) {
            return insert_edge(pcsr, from, to);
        },
        pcsr, from, to
    );
}

bool Dist_PCSR::query_edge(uint32_t from, uint32_t to) {
    uint64_t edge = make_edge_tuple(from, to);
    return spma.query(edge);
}

upcxx::future<bool> query_edge(upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, uint32_t to) {
    uint32_t rank = target_rank(from);
    if (rank == upcxx::rank_me()) {
        return upcxx::make_future(pcsr->query_edge(from, to));
    }
    return upcxx::rpc(rank, 
        [](upcxx::dist_object<DistPCSR>& local_pcsr, uint32_t from, uint32_t to) {
            return query_edge(local_pcsr, from, to);
        },
        pcsr, from, to
    );
}

uint32_t redis_count(
    upcxx::dist_object<DistPCSR>& pcsr, 
    uint32_t start_proc, 
    uint32_t num_procs_contact,
    vector<uint32_t>& element_counts
) {
    uint32_t count = 0;
    
    upcxx::future<> count_future = upcxx::make_future();
    for (int i = 0; i < num_procs_contact; i += 1) {
        upcxx::future<> fut = upcxx::rpc(
            start_proc + i,
            [](upcxx::dist_object<DistPCSR>& pcsr) {
                pcsr->redistributing = true;
                return make_pair(pcsr->spma._num_elements, upcxx::rank_me());
            },
            pcsr
        ).then(
            [&count, &element_counts] (uint32_t num_elements, uint32_t rank) {
                count += num_elements;
                element_counts[rank] = num_elements;
            }
        )
        count_future = upcxx::when_all(count_future, fut);
    }

    count_future.wait();
    return count;
}
// creates a list of (remote start index, remote end index, processor id to retrieve from), one for each retriever.
// for all processors participating in redistribute. each processor is responsible for acquiring data
// it needs from other processors
vector<pair<uint32_t, vector<tuple<uint32_t, uint32_t, uint32_t>>>> get_redistribute_commands(
    vector<uint32_t>& element_counts, 
    uint32_t start_proc, 
    uint32_t num_procs, 
    uint32_t total_elements
) {
    vector<pair<uint32_t, vector<tuple<uint32_t, uint32_t, uint32_t>>>> command_lists;
    uint32_t elements_per_proc_floor = total_elements / num_procs;
    uint32_t current_local_position = 0;
    uint32_t remote_proc = start_proc;
    for (int i = 0; i < num_procs; i += 1) {
        uint32_t retrieving_proc = start_proc + i;
        uint32_t elements = element_counts + (retrieving_proc < total_elements % num_procs);
        vector<tuple<uint32_t, uint32_t, uint32_t>> command_list;
        while (elements > 0) {
            uint32_t remaining_elems = element_counts[remote_proc] - current_local_position;
            if (remaining_elems > elements) {
                command_lists.push_back(make_tuple(current_local_position, current_local_position + elements, remote_proc));
                current_local_position += elements;
                elements = 0;
            } 
            else {
                command_lists.push_back(make_tuple(current_local_position, element_counts[remote_proc], remote_proc));
                elements -= remaining_elem
                remote_proc += 1;
                current_local_position = 0;
            }
        }
        command_lists.push_back(make_pair(retrieving_proc, std::move(command_list)));
    }
    return command_lists;
}

upcxx::future<> retrieve_command(upcxx::dist_object<DistPCSR>& pcsr, tuple<uint32_t, uint32_t, uint32_t> cmd, vector<uint64_t> temp) {
    uint32_t remote_start = std::get<0>(cmd);
    uint32_t remote_end = std::get<1>(cmd);
    uint32_t target_proc = std::get<2>(cmd);
    return upcxx::rpc(target_proc,
        [remote_start, remote_end](upcxx::dist_object<DistPCSR>& pcsr){
            return pcsr->spma.get_min_range(remote_start, remote_end);
        }, pcsr).then(
            [](vector<uint64_t>& retrieved_values) {
                temp.insert(temp.end(), retrieved_values.begin(), retrieved_values.end());
            }
        );
}

vector<uint64_t> temp;

void redistribute(upcxx::dist_object<DistPCSR>& pcsr) {
    pcsr->redistributing = true;
    pcsr->acquire_redis_lock();

    // if someone else redistributed with us while we were waiting for lock,
    // making a redistribute unnecessary, then don't redistribute
    if (!pcsr->redistributing) {
        pcsr->release_redis_lock();
        return;
    }

    uint32_t total_elements = pcsr->spma._num_elements;
    vector<uint32_t> element_counts(upcxx::rank_n(), 0);
    uint32_t n_procs = 1;
    uint32_t start_proc = upcxx::rank_me();
    element_counts[start_proc] = total_elements;
    uint32_t level = 0;
    while (total_elements >= (uint32_t) (pcsr->spma._leaf_max * n_procs * size()) && (n_procs < upcxx::rank_n())) {
        n_procs *= 2;
        level++;
        uint32_t new_proc = (start_proc / n_procs) * n_procs;
        if (new_proc < start_proc) {
            total_elements += redis_count(pcsr, new_proc, n_procs / 2, element_counts);
        } else {
            total_elements += redis_count(pcsr, new_proc + n_procs / 2, n_procs / 2, element_counts);
        }
        start_proc = new_proc;
    }
    level -= 1;
    
    upcxx::team redis_team = teams[level];

    if (n_procs == upcxx::rank_n() && total_elements >= (uint32_t) (pcsr->spma._leaf_max * n_procs)) {
        cerr << "PMA too full" << endl;
        exit(0);
    }

    auto commands = get_redistribute_commands(element_counts, start_proc, n_procs, total_elements);
    upcxx::future<> all_team_finished = upcxx::make_future();
    
    for (int i = 0; i < commands.size(); ++i) {
        uint32_t p;
        vector<tuple<uint32_t, uint32_t, uint32_t>> p_commands;
        tie(p, p_commands) = commands[i];
        auto f = upcxx::rpc(p, [level](upcxx::dist_object<DistPCSR>& pcsr) {
            temp.clear();
            upcxx::future<> retrieve_future = upcxx::make_future();
            for (auto command : p_commands) {
                retrieve_future = upcxx::when_all(retrieve_future, retrieve_command(pcsr, command, temp));
            }
            return retrieve_future.then(
                [&pcsr]() {
                    std::sort(temp.begin(), temp.end());
                    pcsr->my_range = make_pair(temp.front(), temp.back());
                    pcsr->cached_ranges[upcxx::rank_me()] = range_t(duration_cast< milliseconds >(system_clock::now().time_since_epoch()), temp.front(), temp.back());
                    return pcsr->cached_ranges[upcxx::rank_me()];
                }
            );
        }, pcsr);
        all_team_finished = upcxx::when_all(all_team_finished, f.then(
            [i]() {
                pcsr->cached_ranges[i + start_proc] = range;
            }
        ));
    }

    all_team_finished.wait();

    upcxx::future<> all_team_swapped;
    for (int i = 0; i < commands.size(); ++i) {
        uint32_t p;
        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> p_commands;
        std::tie(p, p_commands) = commands[i];
        auto f = upcxx::rpc(p, [](upcxx::dist_object<DistPCSR>& pcsr, vector<range_t>& updated_ranges){
            for (int i = 0; i < updated_ranges.size(); i += 1) {
                if (std::get<0>(pcsr->cached_ranges[i]) < std::get<0>(updated_ranges[i])) {
                    pcsr->cached_ranges[i] = updated_ranges[i];
                }
            }
            pcsr->spma.swap_data(temp);
            pcsr->redistributing = false;
        }, pcsr, pcsr->cached_ranges);
        all_team_swapped = upcxx::when_all(all_team_swapped, f);
    }
    
    all_team_swapped.wait();
    for (int i = 0; i < upcxx::rank_n(); i += 1) {
        if (i >= start_proc && i < start_proc + n_procs) {
            continue;
        }
        upcxx::rpc(i, [](upcxx::dist_object<DistPCSR>& pcsr, vector<range_t>& updated_ranges){
            for (int i = 0; i < updated_ranges.size(); i += 1) {
                if (std::get<0>(pcsr->cached_ranges[i]) < std::get<0>(updated_ranges[i])) {
                    pcsr->cached_ranges[i] = updated_ranges[i];
                }
            }
        }, pcsr, pcsr->cached_ranges);
    }

    pcsr->release_redis_lock();
    pcsr->redistributing = false;
    pcsr->spma._leaf_max = min(spma._leaf_max + 0.1, 0.8);
}

vector<pair<uint32_t, vector<uint32_t>>> DistPCSR::adjacency_lists() { // fix
    std::vector<upcxx::future<vector<pair<uint32_t, vector<uint32_t>>>>> futures;
    for (int rank = 0; rank < upcxx::rank_n(); rank++) {
        futures.push_back(upcxx::rpc(
            rank,
            [](upcxx::dist_object<SetPMA>& dist_spma) {
                vector<pair<uint32_t, vector<uint32_t>>> adjacencies;
                uint32_t curr_source = UINT32_MAX;
                uint32_t* curr_source_ptr = &curr_source;
                vector<uint32_t> adjacency;
                auto f = SetPMA::range_func([&adjacencies, &adjacency, curr_source_ptr](uint64_t v) {
                    uint32_t source, dest;
                    source = v >> 32;
                    dest = v;
                    
                    if (source != *curr_source_ptr) {
                        if (*curr_source_ptr != UINT32_MAX) {
                            adjacencies.push_back(
                                std::make_pair(*curr_source_ptr, std::move(adjacency))
                            );
                            adjacency.clear();
                        }
                        *curr_source_ptr = source;
                        adjacency.push_back(dest);
                    } else {
                        adjacency.push_back(dest);
                    }
                });
                dist_spma->range(0, UINT64_MAX - 1, f);
                    
                // handle tail
                if (curr_source != UINT32_MAX) {
                    adjacencies.push_back(
                        std::make_pair(curr_source, std::move(adjacency))
                    );
                }
                return adjacencies;
            }, dist_spma
        ));
    }
    vector<pair<uint32_t, vector<uint32_t>>> combined_adj_list;
    for (const auto& future : futures) {
        const auto& adj_list = future.wait();
        combined_adj_list.insert(combined_adj_list.end(), adj_list.begin(), adj_list.end());
    }
    return combined_adj_list;
}

vector<uint32_t> DistPCSR::edges(uint32_t from) {
    uint32_t rank = target_rank(from);

    return upcxx::rpc(rank,
        [](upcxx::dist_object<SetPMA> &dist_spma, uint32_t from) {
            uint64_t from_padded = ((uint64_t)from) << 32;
            vector<uint32_t> edge_list;
            auto edge_populater = SetPMA::range_func([&edge_list](uint64_t v) {
                edge_list.push_back((uint32_t)v);
            });

            dist_spma->range(from_padded, from_padded | 0xFFFFFFFF, edge_populater);
            return edge_list;
        },
        dist_spma, from
    ).wait();
}