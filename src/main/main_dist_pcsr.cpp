#include "../data-structures/dist_pcsr2.hpp"
#include "../data-structures/graph.hpp"

#include <upcxx/upcxx.hpp>
#include <chrono>
#include <iostream>
#include <random>

mt19937 rng(42);

template <typename T>
void assertEqual(T actual, T expected) {
    if (actual != expected) {
        cout << "Expected " << expected << " got " << actual << endl;
        exit(-1);
    }
}

int main() {
    upcxx::init();
    upcxx::dist_object<DistPCSR> pcsr = DistPCSR::make_dist_pcsr(1 << 10);

    /*
    if (upcxx::rank_me() == 0) {
        pcsr->insert_edge(0, 1);
        pcsr->insert_edge(6, 3);
    }
    else {
        pcsr->insert_edge(1, 7);
        pcsr->insert_edge(8, 9);
    }

    upcxx::barrier();
    pcsr->flush_queue();

    if (upcxx::rank_me() == 0) {
        assertEqual(pcsr->query_edge(0, 1).wait(), true);
        assertEqual(pcsr->query_edge(1, 7).wait(), true);
        assertEqual(pcsr->query_edge(6, 3).wait(), true);
        assertEqual(pcsr->query_edge(8, 9).wait(), true);
    }
    */


    Graph reference;
    const int NUM_ITERATIONS = 1e3;
    const int INSERT_PHASE_SIZE = 1e2;
    const int QUERY_PHASE_SIZE = 1e2;

    std::uniform_int_distribution<int> random_server(0, upcxx::rank_n());
    std::uniform_int_distribution<int> random_vertex(0, 1 << 10);
    std::uniform_real_distribution<> existsDis(0, 1);

    for (int i = 0; i < NUM_ITERATIONS; i++) {
        upcxx::barrier();
        if (i % 2 == 0) {
            for (int j = 0; j < INSERT_PHASE_SIZE; j++) {
                // Generate random server to be called on
                int server = random_server(rng);
                // Generate random edge
                uint32_t from = random_vertex(rng);
                uint32_t to = random_vertex(rng);
                if (upcxx::rank_me() != server) { continue; }
                reference.insert(from, to);
                pcsr->insert_edge(from, to);
            }
            for (int stupid = 0; stupid < 5; stupid++) {
                pcsr->flush_queue();
                upcxx::barrier();
            }
            // sync_dist_pcsr(pcsr);
        }
        else {
            for (int j = 0; j < QUERY_PHASE_SIZE; j++) {
                // Generate random server to perform query
                int server = random_server(rng);
                // Generate random edge
                uint32_t from;
                uint32_t to;
                if (existsDis(rng) < 0.5) {
                    tie(from, to) = reference.random_edge(rng);
                }
                else {
                    tie(from, to) = reference.random_non_edge(rng);
                }
    
                if (upcxx::rank_me() != server) { continue; }
                // Check correctness
                assert(reference.query_edge(from, to) == pcsr.query_edge(from, to).wait());
            }
        }

        if (i % 10000 == 0) {
            cout << i << endl;
        }
    }

    upcxx::finalize();
}