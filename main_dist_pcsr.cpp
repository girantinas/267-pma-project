#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>
#include <fstream> 
#include <sstream>
#include <iostream>
#include <deque>

#include "dist_pcsr.hpp"
// #include "butil.hpp"

using namespace std;

vector<int> bfs(DistPCSR &pcsr, uint32_t source);
vector<double> pagerank(DistPCSR &pcsr);

int main(int argc, char** argv) {
    upcxx::init();
    ifstream infile;
    ofstream outfile;
    
    if (argc < 2) {
        infile.open("dist_pcsr_inserts.txt");
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

    DistPCSR pcsr(1 << 4, 1000);
    int line_number = 0;

    chrono::high_resolution_clock::time_point start_time, end_time;
    std::chrono::duration<double> insert_duration(0), query_duration(0);

    while (getline(infile, line)) {
        istringstream iss(line);
        string command;
        iss >> command;
        if (command == "START_INSERTS") {
            if (line_number != 0) {
                end_time = std::chrono::high_resolution_clock::now();
                query_duration += end_time - start_time;
            }
            start_time = std::chrono::high_resolution_clock::now();
            upcxx::barrier();
            line_number++;
            continue;
        } else if (command == "START_QUERIES") {
            if (line_number != 0) {
                end_time = std::chrono::high_resolution_clock::now();
                insert_duration += end_time - start_time;
            }
            start_time = std::chrono::high_resolution_clock::now();
            upcxx::barrier();
            line_number++;
            continue;
        }
        if (line_number % upcxx::rank_n() != upcxx::rank_me()) {
            line_number++;
            continue;
        }
        if (command == "PUT_EDGE") {
            uint32_t u, v;
            iss >> u >> v;
            pcsr.insert_edge(u, v);
        } else if (command == "QUERY_EDGE") {
            uint32_t u, v;
            iss >> u >> v;
            bool exists = pcsr.query_edge(u, v);
            if (exists) {
                outfile << u << " " << v << " True" << endl;
            } else {
                outfile << u << " " << v << " False" << endl;
            }
            // outfile << endl;
        } else if (command == "GET_OUT_EDGES") {
            uint32_t vertex;
            iss >> vertex;
            vector<uint32_t> out_edges = pcsr.edges(vertex);
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
            // outfile << endl;
        } else if (command == "QUERY_ALL_GRAPH") {
            vector<pair<uint32_t, vector<uint32_t>>> result = pcsr.adjacency_lists();
            if (write_output) {
                for (auto& entry : result) {
                    uint32_t& vertex = entry.first;
                    outfile << "V: " << vertex << ", E:";
                    for (int dest_vertex : entry.second) {
                        outfile << " " << dest_vertex;
                    }
                    outfile << endl;
                }
                // outfile << endl;
            }
        } else {
            cerr << "Received unsupported command: " << command << endl;
        }
        line_number++;
    }
    end_time = std::chrono::high_resolution_clock::now();
    query_duration += end_time - start_time;
    upcxx::barrier();
    if (upcxx::rank_me() == 0) {
        cout << "Inserts total time: " << std::chrono::duration<double>(insert_duration).count() << endl;
        cout << "Queries total time: " << std::chrono::duration<double>(query_duration).count() << endl;
    }
    infile.close();
    outfile.close();
    upcxx::barrier();
    start_time = std::chrono::high_resolution_clock::now();
    vector<int> distances = bfs(pcsr, 0);
    upcxx::barrier();
    end_time = std::chrono::high_resolution_clock::now();
    if (upcxx::rank_me() == 0) {
        cout << "BFS total time: " << std::chrono::duration<double>(end_time - start_time).count() << endl;
    }
    if (write_output) {
        outfile.open(std::string(argv[2]) + "_" + std::to_string(upcxx::rank_me()) + ".bfs");
        uint32_t start_vertex = pcsr.num_vertices / upcxx::rank_n() * upcxx::rank_me();
        for (uint32_t i = 0; i < distances.size(); i++) {
            uint32_t vertex = start_vertex + i;
            outfile << vertex << " " << distances[i] << endl;
        }
        outfile.close();
    }

    upcxx::barrier();
    start_time = std::chrono::high_resolution_clock::now();
    vector<double> pr_values = pagerank(pcsr);
    upcxx::barrier();
    end_time = std::chrono::high_resolution_clock::now();
    if (upcxx::rank_me() == 0) {
        cout << "PageRank total time: " << std::chrono::duration<double>(end_time - start_time).count() << endl;
    }
    if (write_output) {
        outfile.open(std::string(argv[2]) + "_" + std::to_string(upcxx::rank_me()) + ".pr");
        uint32_t start_vertex = pcsr.num_vertices / upcxx::rank_n() * upcxx::rank_me();
        for (uint32_t i = 0; i < pr_values.size(); i++) {
            uint32_t vertex = start_vertex + i;
            outfile << vertex << " " << pr_values[i] << endl;
        }
        outfile.close();
    }
    upcxx::finalize();
    return 0;
}

vector<int> bfs(DistPCSR &pcsr, uint32_t source) {
    uint32_t start_vertex = pcsr.num_vertices / upcxx::rank_n() * upcxx::rank_me();
    uint32_t end_vertex = start_vertex + pcsr.num_vertices / upcxx::rank_n();
    uint32_t local_num_vertices = end_vertex - start_vertex;

    vector<int> local_distances(local_num_vertices, -1);

    // cout << "Checkpoint 1" << endl;

    if (source >= start_vertex && source < end_vertex) {
        local_distances[source - start_vertex] = 0;
    }

    deque<uint32_t> frontiers;
    if (upcxx::rank_me() == pcsr.target_rank(source)) {
        frontiers.push_back(source);
    }
    int level = 0;

    upcxx::dist_object<vector<uint32_t>> next_set_recv{vector<uint32_t>()};

    // cout << "Checkpoint 2" << endl;
    upcxx::barrier();

    while (true) {
        // cout << "Checkpoint 2.1" << endl;

        vector<uint32_t> next_set;
        // cout << "Checkpoint 2.2" << endl;
        bool is_frontier_empty = frontiers.empty();
        // cout << "Checkpoint 2.3" << endl;
        bool all_frontiers_empty = upcxx::reduce_all(is_frontier_empty, upcxx::op_fast_bit_and).wait();
        // cout << "Checkpoint 2.4" << endl;
        if (all_frontiers_empty) {
            // cout << "Checkpoint 2.5" << endl;
            break;
        }

        // cout << "Checkpoint 3" << endl;

        for (uint32_t vertex : frontiers) {
            

            vector<uint32_t> neighbors = pcsr.edges(vertex);
            next_set.insert(next_set.end(), neighbors.begin(), neighbors.end());

            // for (uint32_t neighbor : neighbors) {
            //     if (neighbor >= start_vertex && neighbor < end_vertex) {
            //         if ((*local_distances)[neighbor - start_vertex] == -1) {
            //             (*local_distances)[neighbor - start_vertex] = current_distance + 1;
            //             global_queue.push_back(neighbor);
            //         }
            //     } else {
            //         int target_rank = pcsr.target_rank(neighbor);
            //         int dest_rank_local_index = neighbor - pcsr.target_rank(neighbor) * local_num_vertices
            //         int neighbor_distance = local_distances->fetch().wait()
            //     }
            // }
        }

        // cout << "Checkpoint 4" << endl;

        vector<vector<uint32_t>> vertex_destinations;
        vertex_destinations.resize(upcxx::rank_n());

        for (uint32_t next_vertex : next_set) {
            int target_rank = pcsr.target_rank(next_vertex);
            vertex_destinations[target_rank].push_back(next_vertex);
        }

        (*next_set_recv).clear();

        upcxx::barrier();

        // cout << "Checkpoint 5" << endl;

        for (uint32_t target_rank = 0; target_rank < upcxx::rank_n(); target_rank++) {
            if (target_rank == upcxx::rank_me()) {
                (*next_set_recv).insert((*next_set_recv).end(), vertex_destinations[target_rank].begin(), vertex_destinations[target_rank].end());
                continue;
            }

            upcxx::rpc(
                target_rank, 
                [](upcxx::dist_object<vector<uint32_t>> &next_set_recv, const vector<uint32_t>& vertices) {
                    (*next_set_recv).insert((*next_set_recv).end(), vertices.begin(), vertices.end());
                },
                next_set_recv, vertex_destinations[target_rank] // Should I capture this instead?
            ).wait();

        }

        upcxx::barrier();
        frontiers.clear();

        // cout << "Checkpoint 6" << endl;
        for (uint32_t vertex : *next_set_recv) {
            // cout << "Another checkpoint 1" << endl;
            int val = local_distances[vertex - start_vertex];
            // if (level == 0)
            // cout << "Another checkpoint 2. Val is " << val << endl;
            // cout << "vertex - start_vertex is" << vertex - start_vertex << endl;
            if (local_distances[vertex - start_vertex] == -1) {
                local_distances[vertex - start_vertex] = level + 1;
                // cout << "Another checkpoint 3" << endl;
                frontiers.push_back(vertex);
                // cout << "Another checkpoint 4" << endl;
            }
        }

        level++;
    }
    return local_distances;
    // cout << "Checkpoint inf" << endl;
}

vector<double> pagerank(DistPCSR &pcsr) {
    double damping_factor = 0.85;
    int max_iterations = 100;
    uint32_t start_vertex = pcsr.num_vertices / upcxx::rank_n() * upcxx::rank_me();
    uint32_t end_vertex = start_vertex + pcsr.num_vertices / upcxx::rank_n();
    uint32_t local_num_vertices = end_vertex - start_vertex;

    vector<double> local_pagerank(local_num_vertices, 1.0 / pcsr.num_vertices);
    vector<double> local_prev_pagerank(local_num_vertices, 0.0);

    upcxx::dist_object<vector<pair<uint32_t, double>>> contributions_recv{vector<pair<uint32_t, double>>()};

    int iteration = 0;

    upcxx::barrier();

    while (iteration < max_iterations) {
        local_prev_pagerank = local_pagerank;
        local_pagerank.assign(local_num_vertices, (1.0 - damping_factor) / pcsr.num_vertices);

        vector<vector<pair<uint32_t, double>>> contributions_destinations;
        contributions_destinations.resize(upcxx::rank_n());

        for (uint32_t vertex = start_vertex; vertex < end_vertex; vertex++) {
            vector<uint32_t> out_edges = pcsr.edges(vertex);
            double contribution = local_prev_pagerank[vertex - start_vertex] / out_edges.size();
            for (uint32_t neighbor : out_edges) {
                int target_rank = pcsr.target_rank(neighbor);
                contributions_destinations[target_rank].push_back(make_pair(neighbor, contribution));
            }
        }
        
        (*contributions_recv).clear();

        upcxx::barrier();

        for (int target_rank = 0; target_rank < upcxx::rank_n(); target_rank++) {
            if (target_rank == upcxx::rank_me()) {
                (*contributions_recv).insert((*contributions_recv).end(), contributions_destinations[target_rank].begin(), contributions_destinations[target_rank].end());
            } else {
                upcxx::rpc(
                    target_rank, [](auto& contributions_recv, const vector<pair<uint32_t, double>>& contributions) {
                        (*contributions_recv).insert((*contributions_recv).end(), contributions.begin(), contributions.end());
                    },
                    contributions_recv, contributions_destinations[target_rank] // Should I capture this instead?
                ).wait(); // Join these together then wait on all at once if we feel like 
            }
        }

        upcxx::barrier();

        for (pair<uint32_t, double> contribution_pair : (*contributions_recv)) {
            uint32_t vertex = contribution_pair.first;
            double contribution = contribution_pair.second;
            local_pagerank[vertex - start_vertex] += contribution * damping_factor;
        }
        iteration++;
        // cout << "Iteration: " << iteration << endl;
    }
    return local_pagerank;
}
