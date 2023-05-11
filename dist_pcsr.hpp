#pragma once

#include <upcxx/upcxx.hpp>
// #include "butil.hpp"
#include <climits>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <deque>
#include "set_pma.hpp"

using namespace std;

// struct Insert {
//     uint32_t from;
//     uint32_t to;
//     Insert(uint32_t from, uint32_t to) : from(from), to(to) {}
// };

class DistPCSR {
    public:
        DistPCSR(uint32_t size, uint64_t max_val);
        uint32_t size();
        uint32_t num_vertices; // Restrict this to power of 2
        /* Inserts an edge from one vertex to another. Fails if either vertex does not exist. */
        void insert_edge(uint32_t from, uint32_t to);
        /* Queries if an edge exists. */
        bool query_edge(uint32_t from, uint32_t to);
        /* Edge list of a vertex. If the vertex does not exist, the list is empty. */
        vector<uint32_t> edges(uint32_t from);
        /* Returns the entire adjacency list representation of the graph. */
        vector<pair<uint32_t, vector<uint32_t>>> adjacency_lists();
        // deque<Insert> rq_queue;
        uint32_t num_elements();
        static uint64_t make_edge_tuple(uint32_t from, uint32_t to);
        static pair<uint32_t, uint32_t> get_edge_tuple(uint64_t edge);
        uint32_t target_rank(uint32_t from);
        upcxx::dist_object<SetPMA> dist_spma;
        const uint32_t TO_ONES = 0xFFFFFFFF;
};


DistPCSR::DistPCSR(uint32_t size, uint64_t max_val) : dist_spma(SetPMA(size)) {
    num_vertices = max_val;
}

uint64_t DistPCSR::make_edge_tuple(uint32_t from, uint32_t to) { return (((uint64_t) from) << 32) | to; }
pair<uint32_t, uint32_t> DistPCSR::get_edge_tuple(uint64_t edge) { return std::make_pair((uint32_t)(edge >> 32), (uint32_t)((edge << 32) >> 32)); } /* doesn't work */

uint32_t DistPCSR::target_rank(uint32_t from) {
    return from / (num_vertices / upcxx::rank_n());
}

uint32_t DistPCSR::num_elements() {
    return dist_spma->_num_elements;
}

void DistPCSR::insert_edge(uint32_t from, uint32_t to) {
    uint32_t rank = target_rank(from);
    uint64_t edge = make_edge_tuple(from, to);
    if (upcxx::rank_me() == rank) {
        dist_spma->insert(edge);
        return;
    }
    upcxx::rpc(rank,
        [](upcxx::dist_object<SetPMA> &dist_spma, uint64_t edge) {
            dist_spma->insert(edge);
        },
        dist_spma, edge
    ).wait();
}

bool DistPCSR::query_edge(uint32_t from, uint32_t to) {
    uint32_t rank = target_rank(from);
    uint64_t edge = make_edge_tuple(from, to);
    return upcxx::rpc(rank, 
        [](upcxx::dist_object<SetPMA>& dist_spma, uint64_t edge) {
            return dist_spma->query(edge);
        },
        dist_spma, edge
    ).wait();
}

vector<uint32_t> DistPCSR::edges(uint32_t from) {
    uint32_t rank = target_rank(from);
    if (rank == upcxx::rank_me()) {
        uint64_t from_padded = ((uint64_t)from) << 32;
        vector<uint32_t> edge_list;
        auto edge_populater = SetPMA::range_func([&edge_list](uint64_t v) {
            edge_list.push_back((uint32_t)v);
        });

        dist_spma->range(from_padded, from_padded | 0xFFFFFFFF, edge_populater);
        return edge_list;
    }
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

vector<pair<uint32_t, vector<uint32_t>>> DistPCSR::adjacency_lists() {
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
