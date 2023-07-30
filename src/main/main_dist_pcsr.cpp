#include "../data-structures/dist_pcsr.hpp"
#include "../data-structures/graph.hpp"

#include <upcxx/upcxx.hpp>
#include <chrono>
#include <iostream>
#include <random>

mt19937 rng(42);

// fully synchronizes a distributed PCSR
void sync_dist_pcsr(upcxx::dist_object<DistPCSR>& pcsr) {
    // finish up all your inserts
    int loop_count = 0;
    auto watchdog = std::chrono::high_resolution_clock::now();
    bool printed = false;

    while (pcsr->outstanding_rpcs != 0 || !pcsr->rq_queue.empty() || pcsr->redistributing) {
        /*if (upcxx::rank_me() == 0) {
            if (pcsr->redistributing) {
                cout << "redistributing..." << endl;
            }
            else if (pcsr->outstanding_rpcs != 0) {
                cout << "waiting on rpcs, " << pcsr->outstanding_rpcs << " remaining" << endl;
            }
            else {
                cout << "request queue empty" << endl;
            }
        }*/
        flush_queue(pcsr);
        /*if (upcxx::rank_me() == 0) {
            cout << "flushed queue" << endl;
        }*/
        upcxx::progress();
        /*if (upcxx::rank_me() == 0) {
            cout << "finished progress" << endl;
        }*/
    }
    
    // at this point, all of OUR insert futures are satisfied, OUR queue is empty, and we are not redistributing
    auto finished_inserts = upcxx::barrier_async();
    while (pcsr->outstanding_rpcs != 0 || !finished_inserts.ready() || !pcsr->rq_queue.empty() || pcsr->redistributing) {
        upcxx::progress();
        flush_queue(pcsr);
        if (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - watchdog).count() > 15.0 && !printed) {
            cout << "rank " << upcxx::rank_me() << "seems stuck" << endl;
            cout << "outstanding: " << pcsr->outstanding_rpcs << endl;
            cout << "request_queue: " << pcsr->rq_queue.size() << endl;
            cout << "redistributing: " << pcsr->redistributing << endl; 
            printed = true;
        }
    }
    finished_inserts.wait();
    // if (upcxx::rank_me() == 0) cout << "reached first barrier in sync_dist_pcsr" << endl;
    finished_inserts = upcxx::barrier_async();
    while (pcsr->outstanding_rpcs != 0 || !finished_inserts.ready() || !pcsr->rq_queue.empty() || pcsr->redistributing) {
        upcxx::progress();
        flush_queue(pcsr);
        if (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - watchdog).count() > 15.0 && !printed) {
            cout << "rank " << upcxx::rank_me() << "seems stuck" << endl;
            cout << "outstanding: " << pcsr->outstanding_rpcs << endl;
            cout << "request_queue: " << pcsr->rq_queue.size() << endl;
            cout << "redistributing: " << pcsr->redistributing << endl; 
            printed = true;
        }
    }
    finished_inserts.wait();
    // if (upcxx::rank_me() == 0) cout << "reached second barrier in sync_dist_pcsr" << endl;
}

int main() {
    DistPCSR pcsr = DistPCSR::make_dist_pcsr(1 << 10);
    Graph reference;
    const int NUM_ITERATIONS = 1e4;
    const int INSERT_PHASE_SIZE = 1e4;
    const int QUERY_PHASE_SIZE = 1e2;

    std::uniform_int_distribution<int> random_server(0, upcxx::rank_n());
    std::uniform_int_distribution<int> random_vertex(0, 1 << 20);
    std::uniform_real_distribution<> existsDis(0, 1);

    for (int i = 0; i < NUM_ITERATIONS; i++) {
        if (i % 2 == 0) {
            for (int j = 0; j < INSERT_PHASE_SIZE; j++) {
                // Generate random server to be called on
                int server = random_server(rng);
                if (upcxx::rank_me() != server) { continue; }
                // Generate random edge
                uint32_t from = random_vertex(rng);
                uint32_t to = random_vertex(rng);
                reference.insert(from, to);
                pcsr.insert(from, to);
            }
            // sync_dist_pcsr(pcsr);
        }
        else {
            for (int j = 0; j < QUERY_PHASE_SIZE; j++) {
                // Generate random server to perform query
                int server = random_server(rng);
                if (upcxx::rank_me() != server) { continue; }
                
                uint32_t from;
                uint32_t to;
                if (existsDis(rng) < 0.5) {
                    tie(from, to) = reference.random_edge(rng);
                }
                else {
                    tie(from, to) = reference.random_non_edge(rng);
                }
    
                // Check correctness
                assert(reference.query_edge(from, to) == pcsr.query_edge(from, to).wait());
            }
        }

        if (i % 10000 == 0) {
            cout << i << endl;
        }
    }
}