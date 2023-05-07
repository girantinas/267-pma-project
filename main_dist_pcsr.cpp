#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <thread>
#include <upcxx/upcxx.hpp>
#include <vector>
#include <fstream> 
#include <sstream>
#include <iostream>

#include "dist_pcsr.hpp"
// #include "butil.hpp"

using namespace std;

// fully synchronizes a distributed PCSR
void sync_dist_pcsr(upcxx::dist_object<DistPCSR>& pcsr) {
    // finish up all your inserts
    int loop_count = 0;
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
        loop_count += 1;
    }
    
    // at this point, all of OUR insert futures are satisfied, OUR queue is empty, and we are not redistributing
    auto finished_inserts = upcxx::barrier_async();
    while (pcsr->outstanding_rpcs != 0 || !finished_inserts.ready() || !pcsr->rq_queue.empty() || pcsr->redistributing) {
        upcxx::progress();
        flush_queue(pcsr);
    }
    finished_inserts.wait();
    // if (upcxx::rank_me() == 0) cout << "reached first barrier in sync_dist_pcsr" << endl;
    finished_inserts = upcxx::barrier_async();
    while (pcsr->outstanding_rpcs != 0 || !finished_inserts.ready() || !pcsr->rq_queue.empty() || pcsr->redistributing) {
        upcxx::progress();
        flush_queue(pcsr);
    }
    finished_inserts.wait();
    // if (upcxx::rank_me() == 0) cout << "reached second barrier in sync_dist_pcsr" << endl;
}

int main(int argc, char** argv) {
    if (std::strcmp(std::getenv("GASNET_OFI_RECEIVE_BUFF_SIZE"), "single")) {
        cout << "Critical error: GASNET_OFI_RECEIVE_BUFF_SIZE workaround not detected (value=" << std::getenv("GASNET_OFI_RECEIVE_BUFF_SIZE") << ")" << endl;
        exit(-1);
    }
    cout << "rank " << upcxx::rank_me() << " alive" << endl;
    upcxx::init();
    cout << "rank " << upcxx::rank_me() << " initialized" << endl;
    auto start = std::chrono::high_resolution_clock::now();

    ifstream infile;
    int num_files = 16;
    int ranks_per_file = upcxx::rank_n() / num_files; // Make sure both of these are a power of 2.
    infile.open("../rmat-tests/rmat-inserts-" + std::to_string(upcxx::rank_me() / ranks_per_file) + ".txt");
    if (!infile.is_open()) {
        cerr << "Couldn't open LiveJournal files" << endl;
        return -1;
    }

    string line;

    if (upcxx::rank_me() == 0) cout << "(rank 0) finished setting up file I/O at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
    int num_edges = 1 << 29;
    int num_ranks = upcxx::rank_n();

    upcxx::dist_object<DistPCSR> pcsr(DistPCSR((num_edges / num_ranks) * 16, (1 << 24) + 1000));
    if (upcxx::rank_me() == 0) cout << "pma depth" << pcsr->spma.depth() << endl;
    if (upcxx::rank_me() == 0) cout << "pma size" << pcsr->spma.size() << endl; 
    if (upcxx::rank_me() == 0) cout << "finished distpcsr construction at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
    upcxx::barrier();
    auto commands_start = std::chrono::high_resolution_clock::now();
    if (upcxx::rank_me() == 0) cout << "starting commands at " << std::chrono::duration<double>(commands_start - start).count() << endl;
    int line_number = 0;
    upcxx::future<> my_futures;
    bool inserting = false;

    int total_inserts_local = 0;
    while (getline(infile, line)) {
        istringstream iss(line);
        string command;
        iss >> command;
        if (line_number % ranks_per_file != upcxx::rank_me() % ranks_per_file) {
        } else if (command == "BFS") {
            
        } else {
            uint64_t edge = stoull(command);
            uint32_t u, v;
            // cout << "inserting edge " << u << " " << v << endl;
            tie(u, v) = DistPCSR::get_edge_tuple(edge);
            insert_edge_batched(pcsr, u, v);
            total_inserts_local += 1;
        }
        
        if (line_number % 100000 == 0) {
            if (upcxx::rank_me() == 0 && line_number % 1000000 == 0) {
                cout << "processed 10^6 lines, at " << line_number << endl;
            }
            sync_dist_pcsr(pcsr);
            if (upcxx::rank_me() == 0) {
                // cout << "processed all batches" << endl;
            }
        } 
        line_number++;
    }

    infile.close();
    clear_batch_buffers(pcsr);
    sync_dist_pcsr(pcsr);
    upcxx::barrier();
    uint64_t total_inserts = upcxx::reduce_one(total_inserts_local, upcxx::op_fast_add, 0).wait();
    uint64_t total_elements = upcxx::reduce_one(pcsr->num_elements(), upcxx::op_fast_add, 0).wait();
    if (upcxx::rank_me() == 0) {
        cout << "pma total inserts: " << total_inserts << endl;
        cout << "pma total elements: " << total_elements << endl;
    }
    // cout << "rank " << upcxx::rank_me() << " total elements: " << pcsr->num_elements() << endl;
    
    
    
    uint32_t rank = upcxx::rank_me();
    // cout << "rank " << rank << " finished at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
    
    upcxx::finalize();
    if (rank == 0) {
        cout << "all finished at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
        cout << "time spent: " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - commands_start).count() << endl;
    }
    return 0;
}

/*
vector<int> bfs(DistPCSR &pcsr, uint32_t source) {
    uint32_t start_vertex = pcsr.num_vertices / upcxx::rank_n() * upcxx::rank_me();
    uint32_t end_vertex = start_vertex + pcsr.num_vertices / upcxx::rank_n();
    uint32_t local_num_vertices = end_vertex - start_vertex;

    vector<int> local_distances(local_num_vertices, -1);

    if (source >= start_vertex && source < end_vertex) {
        local_distances[source - start_vertex] = 0;
    }

    deque<uint32_t> frontiers;
    if (upcxx::rank_me() == pcsr.target_rank(source)) {
        frontiers.push_back(source);
    }
    int level = 0;

    upcxx::dist_object<vector<uint32_t>> next_set_recv{vector<uint32_t>()};

    upcxx::barrier();

    while (true) {

        vector<uint32_t> next_set;
        bool is_frontier_empty = frontiers.empty();
        bool all_frontiers_empty = upcxx::reduce_all(is_frontier_empty, upcxx::op_fast_bit_and).wait();
        if (all_frontiers_empty) {
            break;
        }

        for (uint32_t vertex : frontiers) {
            vector<uint32_t> neighbors = pcsr.edges(vertex);
            next_set.insert(next_set.end(), neighbors.begin(), neighbors.end());
        }

        vector<vector<uint32_t>> vertex_destinations;
        vertex_destinations.resize(upcxx::rank_n());

        for (uint32_t next_vertex : next_set) {
            int target_rank = pcsr.target_rank(next_vertex);
            vertex_destinations[target_rank].push_back(next_vertex);
        }

        (*next_set_recv).clear();

        upcxx::barrier();

        upcxx::future<> all_futures = upcxx::make_future();

        for (uint32_t target_rank = 0; target_rank < upcxx::rank_n(); target_rank++) {
            if (target_rank == upcxx::rank_me()) {
                (*next_set_recv).insert((*next_set_recv).end(), vertex_destinations[target_rank].begin(), vertex_destinations[target_rank].end());
                continue;
            }
            all_futures = upcxx::when_all(all_futures, upcxx::rpc(
                target_rank, 
                [](upcxx::dist_object<vector<uint32_t>> &next_set_recv, const vector<uint32_t>& vertices) {
                    (*next_set_recv).insert((*next_set_recv).end(), vertices.begin(), vertices.end());
                },
                next_set_recv, vertex_destinations[target_rank] // Should I capture this instead?
            ));

        }

        all_futures.wait();

        upcxx::barrier();
        frontiers.clear();

        for (uint32_t vertex : *next_set_recv) {
            int val = local_distances[vertex - start_vertex];
            if (local_distances[vertex - start_vertex] == -1) {
                local_distances[vertex - start_vertex] = level + 1;
                frontiers.push_back(vertex);
            }
        }

        level++;
    }
    return local_distances;
}
*/
/*
vector<int> bfs_serial(DistPCSR &pcsr, uint32_t source) {
    if (upcxx::rank_me() == 0) cout << "BFS Start" << endl;
    vector<int> distances(pcsr.num_vertices, -1);
    deque<uint32_t> queue;

    distances[source] = 0;
    queue.push_back(source);

    while (!queue.empty()) {
        uint32_t current = queue.front();
        queue.pop_front();

        vector<uint32_t> neighbors = pcsr.edges(current).wait();
        for (uint32_t neighbor : neighbors) {
            if (distances[neighbor] == -1) {
                distances[neighbor] = distances[current] + 1;
                queue.push_back(neighbor);
            }
        }
    }
    if (upcxx::rank_me() == 0) cout << "BFS Finish" << endl;
    return distances;
}
*/