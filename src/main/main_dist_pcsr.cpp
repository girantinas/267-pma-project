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
    upcxx::dist_object<DistPCSR> pcsr = DistPCSR::make_dist_pcsr(1 << 20);

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
    const int INSERT_PHASE_SIZE = 1e4;
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
            
            for (int i = 0; i < 5; i += 1) {
                pcsr->flush_queue();
                upcxx::barrier();
            }

            // A, B, C
            // inactivity = no outgoing RPCs, not redistributing, nothing in insert queue
            // 1. Flush queue
            // 2. Upon inactivity, write 1 to your inactivity flag in processor 0 and spin (making user-level progress)
            // (For processor 0), upon inactivity, check for any inactivity flag being 0 
            // 3. Spin (making user-level progress) until all inactivity flags are 1 (i.e. you think no processor is still active)
            ////////////////////
            //      RMA all processors for their completeness flag (length of flush_queue is 0) again (and stop them from making user-level progress by spinning)
            //      If any flag is 0
            //         Set all flags to 0
            //         RMA to all processors that they should exit spin, redo steps 1 and 2
            //      If all flags are 1
            //         RMA to all processors that they are done with insert phase
            
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