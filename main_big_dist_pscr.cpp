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
    while (pcsr->outstanding_rpcs != 0 || !pcsr->rq_queue.empty() || pcsr->redistributing) {
        upcxx::progress();
        flush_queue(pcsr);
    }
    
    // at this point, all of OUR insert futures are satisfied, OUR queue is empty, and we are not redistributing
    auto finished_inserts = upcxx::barrier_async();
    while (pcsr->outstanding_rpcs != 0 || !finished_inserts.ready() || !pcsr->rq_queue.empty() || pcsr->redistributing) {
        upcxx::progress();
        flush_queue(pcsr);
    }
    finished_inserts.wait();

    finished_inserts = upcxx::barrier_async();
    while (pcsr->outstanding_rpcs != 0 || !finished_inserts.ready() || !pcsr->rq_queue.empty() || pcsr->redistributing) {
        upcxx::progress();
        flush_queue(pcsr);
    }
    finished_inserts.wait();
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
    ofstream outfile;
    
    if (argc < 2) {
        infile.open("???");
        if (!infile.is_open()) {
            cerr << "Couldn't open pcsr_inserts.txt" << endl;
            return -1;
        }
    } else {
        infile.open(argv[1]);
        if (!infile.is_open()) {
            cerr << "Couldn't open " << argv[1] << endl;
            return -1;
        }
    }

    if (argc >= 3) {
        outfile.open(std::string(argv[2]) + "_" + std::to_string(upcxx::rank_me()) + ".dat");
        if (!outfile.is_open()) {
            cerr << "Couldn't open " << argv[2] << endl;
            return -1; 
        }
    }

    bool write_output = outfile.is_open();

    string line;

    if (upcxx::rank_me() == 0) cout << "(rank 0) finished setting up file I/O at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
    upcxx::dist_object<DistPCSR> pcsr(DistPCSR(1 << 8, 16000));
    if (upcxx::rank_me() == 0) cout << "finished distpcsr construction at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
    upcxx::barrier();
    if (upcxx::rank_me() == 0) cout << "starting commands at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
    int line_number = 0;
    upcxx::future<> my_futures;
    while (getline(infile, line)) {
        istringstream iss(line);
        string command;
        iss >> command;
        if (command == "START_INSERTS") {
            upcxx::barrier();
            if (upcxx::rank_me() == 0) cout << "starting insert phase at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
            my_futures = upcxx::make_future();
        } else if (command == "START_QUERIES") {
            // finish up all your inserts
            while (pcsr->outstanding_rpcs != 0 || !pcsr->rq_queue.empty() || pcsr->redistributing) {
                upcxx::progress();
                flush_queue(pcsr);
            }
            
            // at this point, all of OUR insert futures are satisfied, OUR queue is empty, and we are not redistributing
            auto finished_inserts = upcxx::barrier_async();
            while (pcsr->outstanding_rpcs != 0 || !finished_inserts.ready() || !pcsr->rq_queue.empty() || pcsr->redistributing) {
                upcxx::progress();
                flush_queue(pcsr);
            }
            finished_inserts.wait();

            finished_inserts = upcxx::barrier_async();
            while (pcsr->outstanding_rpcs != 0 || !finished_inserts.ready() || !pcsr->rq_queue.empty() || pcsr->redistributing) {
                upcxx::progress();
                flush_queue(pcsr);
            }
            finished_inserts.wait();

            if (upcxx::rank_me() == 0) cout << "starting query phase at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
        } else if (line_number % upcxx::rank_n() != upcxx::rank_me()) {
        } else if (command == "PUT_EDGE") {
            uint32_t u, v;
            iss >> u >> v;
            insert_edge(pcsr, u, v);
        } else if (command == "QUERY_EDGE") {
            uint32_t u, v;
            iss >> u >> v;
            bool exists = query_edge(pcsr, u, v).wait(); // TODO: batch these   
            if (exists) {
                outfile << u << " " << v << " True" << endl;
            } else {
                outfile << u << " " << v << " False" << endl;
            }
            // outfile << endl;
        } else if (command == "GET_OUT_EDGES") {
            uint32_t vertex;
            iss >> vertex;
            vector<uint32_t> out_edges;
            edges(pcsr, vertex, out_edges).wait();
            std::sort(out_edges.begin(), out_edges.end());
            if (write_output) {
                if (out_edges.empty()) {
                    outfile << "V: " << vertex << ", E: ";
                } else {
                    outfile << "V: " << vertex << ", E:";
                    for (uint32_t dest_vertex : out_edges) {
                        outfile << " " << dest_vertex;
                    }
                }
                outfile << endl;
            }
        } else if (command == "QUERY_ALL_GRAPH") {
            /*
            vector<pair<uint32_t, vector<uint32_t>>> result = adjacency_lists(pcsr);
            if (write_output) {
                for (auto& entry : result) {
                    uint32_t& vertex = entry.first;
                    outfile << "V: " << vertex << ", E:";
                    for (int dest_vertex : entry.second) {
                        outfile << " " << dest_vertex;
                    }
                    outfile << endl;
                }
                outfile << endl;
            }
            */
        } else {
            cerr << "Received unsupported command: " << command << endl;
        }
        line_number++;
    }

    infile.close();
    outfile.close();

    pcsr->print_dist_pcsr();
    
    uint32_t rank = upcxx::rank_me();
    cout << "rank " << rank << " finished at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
    
    upcxx::finalize();
    if (rank == 0) cout << "all finished at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
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