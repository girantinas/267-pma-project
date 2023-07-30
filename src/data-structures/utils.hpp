#pragma once
#include <cstdint>
#include <utility>

typedef uint32_t vertex_t;
typedef uint64_t edge_t; // upper 32 bits are source, lower 32 bits are destination

uint64_t MSSB(uint64_t x) {
    uint64_t i = 0;
    while (x != 0) {
        x >>= 1;
        ++i;
    }
    return i - 1;
}

uint64_t next_power_of_2(uint64_t x) {
    return 1 << (MSSB(x) + 1);
}

edge_t make_edge_tuple(vertex_t from, vertex_t to) { 
    return (((edge_t) from) << 32) | to; 
}

std::pair<vertex_t, vertex_t> get_edge_tuple(edge_t edge) { 
    return std::make_pair((vertex_t)(edge >> 32), (vertex_t)((edge << 32) >> 32)); 
}
