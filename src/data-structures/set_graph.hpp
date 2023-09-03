#include "utils.hpp"
#include <random>
#include <set>
#include <unordered_map>
#include <utility>
#include <deque>

class SetGraph {
    public:
        void insert(vertex_t from, vertex_t to);
        bool query_edge(vertex_t from, vertex_t to);

        unordered_map<vertex_t, int> bfs(vertex_t source);
        vector<double> pagerank();

        pair<vertex_t, vertex_t> random_edge(mt19937& rng);
        pair<vertex_t, vertex_t> random_non_edge(mt19937& rng);
    private:
        std::set<edge_t> data;
        std::unordered_map<vertex_t, std::set<vertex_t>> adj_set;
};

void SetGraph::insert(vertex_t from, vertex_t to) {
    data.insert(make_edge_tuple(from, to));
    if (adj_set.find(from) == adj_set.end()) {
        adj_set[from] = std::set<vertex_t>();
    }
    adj_set[from].insert(to);   
}

bool SetGraph::query_edge(vertex_t from, vertex_t to) {
    edge_t edge = make_edge_tuple(from, to);
    return data.find(edge) != data.end();
}

pair<vertex_t, vertex_t> SetGraph::random_edge(mt19937& rng) {
    assert(data.size() > 0);
    auto it = data.end();
    while (it == data.end()) {
        edge_t larger = rng();
        it = data.lower_bound(larger);
    }

    edge_t value_in_set = *it;
    return get_edge_tuple(value_in_set);
}

pair<vertex_t, vertex_t> SetGraph::random_non_edge(mt19937& rng) {
    edge_t value_not_in_set;
    auto it = data.begin();
    while (it != data.end()) {
        value_not_in_set = rng();
        it = data.find(value_not_in_set);
    }

    return get_edge_tuple(value_not_in_set);
}

std::unordered_map<vertex_t, int> SetGraph::bfs(vertex_t source) {
    std::deque<std::pair<vertex_t, int>> bfs_queue;
    std::set<vertex_t> visited;
    std::unordered_map<vertex_t, int> distances;
    
    bfs_queue.push_back(std::make_pair(source, 0));
    visited.insert(source);
    distances[source] = 0;
    while (!bfs_queue.empty()) {
        vertex_t v;
        int d;

        std::tie(v, d) = bfs_queue.front();
        bfs_queue.pop_front();
        distances[v] = d;

        auto lo = data.lower_bound(make_edge_tuple(v, 0));
        auto hi = data.upper_bound(make_edge_tuple(v, UINT32_MAX));
        while (lo != hi) {
            vertex_t neighbor = get_edge_tuple(*lo).second;
            if (visited.find(neighbor) == visited.end()) {
                bfs_queue.push_back(std::make_pair(neighbor, d + 1));
                visited.insert(neighbor);
            }
            lo++;
        }
    }

    return distances;
}

// requires num_vertices to be correct, i.e. # vertices is correctly known before inserts
vector<double> SetGraph::pagerank() {
    
}