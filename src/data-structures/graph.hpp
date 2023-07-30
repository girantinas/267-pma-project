#include "utils.hpp"
#include <random>
#include <set>

class Graph {
    public:
        void insert(vertex_t from, vertex_t to);
        bool query_edge(vertex_t from, vertex_t to);
        pair<vertex_t, vertex_t> random_edge(mt19937& rng);
        pair<vertex_t, vertex_t> random_non_edge(mt19937& rng);
    private:
        std::set<edge_t> data;
};

void Graph::insert(vertex_t from, vertex_t to) {
    data.insert(make_edge_tuple(from, to));
}

bool Graph::query_edge(vertex_t from, vertex_t to) {
    edge_t edge = make_edge_tuple(from, to);
    return data.find(edge) != data.end();
}

pair<vertex_t, vertex_t> Graph::random_edge(mt19937& rng) {
    assert(data.size() > 0);
    auto it = data.end();
    while (it == data.end()) {
        edge_t larger = rng();
        it = data.lower_bound(larger);
    }

    edge_t value_in_set = *it;
    return get_edge_tuple(value_in_set);
}

pair<vertex_t, vertex_t> Graph::random_non_edge(mt19937& rng) {
    edge_t value_not_in_set;
    auto it = data.begin();
    while (it != data.end()) {
        value_not_in_set = rng();
        it = data.find(value_not_in_set);
    }

    return get_edge_tuple(value_not_in_set);
}