#include <random>
#include <set>

class Graph {
    public:
        void insert(uint32_t from, uint32_t to);
        bool query_edge(uint32_t from, uint32_t to);
        pair<uint32_t, uint32_t> random_edge(mt19937& rng);
        pair<uint32_t, uint32_t> random_non_edge(mt19937& rng);
    private:
        std::set<uint64_t> data;
};

static uint64_t make_edge_tuple(uint32_t from, uint32_t to) { return (((uint64_t) from) << 32) | to; }
static pair<uint32_t, uint32_t> get_edge_tuple(uint64_t edge) { return std::make_pair((uint32_t)(edge >> 32), (uint32_t)((edge << 32) >> 32)); }

void Graph::insert(uint32_t from, uint32_t to) {
    data.insert(make_edge_tuple(from, to));
}

bool Graph::query_edge(uint32_t from, uint32_t to) {
    uint64_t edge = make_edge_tuple(from, to);
    return data.find(edge) != data.end();
}

pair<uint32_t, uint32_t> Graph::random_edge(mt19937& rng) {
    assert(data.size() > 0);
    auto it = data.end();
    while (it == data.end()) {
        uint64_t larger = rng();
        it = data.lower_bound(larger);
    }

    uint64_t value_in_set = *it;
    return get_edge_tuple(value_in_set);
}

pair<uint32_t, uint32_t> Graph::random_non_edge(mt19937& rng) {
    uint64_t value_not_in_set;
    auto it = data.begin();
    while (it != data.end()) {
        value_not_in_set = rng();
        it = data.find(value_not_in_set);
    }

    return get_edge_tuple(value_not_in_set);
}