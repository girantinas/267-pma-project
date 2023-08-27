#include "utils.hpp"
#include <random>

class Graph {
    public:
        virtual void insert(vertex_t from, vertex_t to) = 0;
        virtual bool query_edge(vertex_t from, vertex_t to) = 0;

        virtual vector<int> bfs(vertex_t source) = 0;
        virtual vector<double> pagerank() = 0;
    private:
    // all the helpers
};
