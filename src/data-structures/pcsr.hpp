#pragma once

// #include <upcxx/upcxx.hpp>
// #include "butil.hpp"
#include <climits>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include "set_pma.hpp"

using namespace std;

class PCSR {
    public:
        PCSR(uint32_t size);
        uint32_t size();

        /* Inserts an edge from one vertex to another. Fails if either vertex does not exist. */
        void insert_edge(uint32_t from, uint32_t to);
        /* Queries if an edge exists. */
        bool query_edge(uint32_t from, uint32_t to);
        /* Edge list of a vertex. If the vertex does not exist, the list is empty. */
        vector<uint32_t> edges(uint32_t from);
        /* Returns the entire adjacency list representation of the graph. */
        vector<pair<uint32_t, vector<uint32_t>>> adjacency_lists();
    private:
        static uint64_t make_edge_tuple(uint32_t from, uint32_t to);
        static pair<uint32_t, uint32_t> get_edge_tuple(uint64_t edge);
        SetPMA spma;
        const uint32_t TO_ONES = 0xFFFFFFFF;
};


PCSR::PCSR(uint32_t size) : spma(size) {}

uint64_t PCSR::make_edge_tuple(uint32_t from, uint32_t to) { return (((uint64_t) from) << 32) | to; }
pair<uint32_t, uint32_t> PCSR::get_edge_tuple(uint64_t edge) { return std::make_pair((uint32_t)(edge >> 32), (uint32_t)((edge << 32) >> 32)); } /* doesn't work */

void PCSR::insert_edge(uint32_t from, uint32_t to) {
    uint64_t edge = make_edge_tuple(from, to);
    spma.insert(edge);
}

bool PCSR::query_edge(uint32_t from, uint32_t to) {
    uint64_t edge = make_edge_tuple(from, to);
    return spma.query(edge);
}

vector<uint32_t> PCSR::edges(uint32_t from) {
    uint64_t from_padded = ((uint64_t) from) << 32;
    vector<uint32_t> edge_list;
    auto edge_populater = SetPMA::range_func([&edge_list](uint64_t v) { 
        edge_list.push_back((uint32_t) v); 
    });
    spma.range(from_padded, from_padded | TO_ONES, edge_populater);
    return edge_list;
}

vector<pair<uint32_t, vector<uint32_t>>> PCSR::adjacency_lists() {
    vector<pair<uint32_t, vector<uint32_t>>> adjacencies;
    uint32_t curr_source = UINT32_MAX; // don't use this as a vertex
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
        }
        else {
            adjacency.push_back(dest);
        }
    });
    spma.range(0, UINT64_MAX - 1, f);
    
    // handle tail
    if (curr_source != UINT32_MAX) {
        adjacencies.push_back(
            std::make_pair(curr_source, std::move(adjacency))
        );
    }

    return adjacencies;
}