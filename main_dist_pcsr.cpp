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
vector<int> bfs_serial(DistPCSR &pcsr, uint32_t source);
vector<int> bfs(DistPCSR &pcsr, uint32_t source);
vector<double> pagerank(DistPCSR &pcsr);

int main_2(int argc, char** argv) {
    upcxx::init();
    if (upcxx::rank_me() == 0) {
        cout << "alive" << endl;
        unordered_set<uint64_t> edges;
        ifstream infile;
        string line;
        uint32_t num_dupes = 0;
        for (int i = 0; i < 16; i++) {
            infile.open(std::string(argv[1]) + "-" + std::to_string(i) + ".txt");
            
            while (getline(infile, line)) {
                istringstream iss(line);
                string command;
                iss >> command;
                // cout << command << endl;
                if (command == "BFS") {
                } else {
                    uint64_t edge = stoull(command);
                    if (edges.find(edge) != edges.end()) {
                        num_dupes++;
                    } else {
                        edges.insert(edge);
                    }
                }
                
            }
            infile.close();
            cout << "finished a file" << endl;
        }
        cout << "Num dupes in insert files: " << num_dupes << endl;
    }
    upcxx::finalize();
    return 0;
}

int main(int argc, char** argv) {
    upcxx::init();
    int num_files = 16;
    int ranks_per_file = upcxx::rank_n() / num_files; // Make sure both of these are a power of 2.
    // cout << "Checkpoint 1" << endl;
    ifstream infile;
    ofstream outfile;
    if (argc < 2) {
        infile.open("dist_pcsr_inserts.txt");
        if (!infile.is_open()) {
            cerr << "Couldn't open pcsr_inserts.txt" << endl;
            return -1;
        }
    // cout << "Checkpoint 2" << endl;
    } else {
        // infile.open(argv[1]);
        infile.open(std::string(argv[1]) + "-" + std::to_string(upcxx::rank_me() / ranks_per_file) + ".txt");
        if (!infile.is_open()) {
            cerr << "Couldn't open " << std::string(argv[1]) + "-" + std::to_string(upcxx::rank_me() / ranks_per_file) + ".txt" << endl;
            return -1;
        }
    }
    // cout << "Checkpoint 3" << endl;
    if (argc >= 3) {
        outfile.open(std::string(argv[2]) + "_" + std::to_string(upcxx::rank_me()) + ".dat");
        if (!outfile.is_open()) {
            cerr << "Couldn't open " << argv[2] << endl;
            return -1; 
        }
    }

    bool write_output = outfile.is_open();
    // cout << "Checkpoint 4" << endl;
    string line;
    int num_edges = 1 << 29;
    DistPCSR pcsr(1 << 27, (1 << 24) + 1000);
    // DistPCSR pcsr(num_edges / upcxx::rank_n() * 4, (1 << 24) + 1000);
    int line_number = 0;
    int num_inserts = 0;
    int num_bfs = 0;
    // cout << "Checkpoint 5" << endl;
    chrono::high_resolution_clock::time_point start_time, end_time;
    std::chrono::duration<double> insert_duration(0), query_duration(0), bfs_duration(0), pagerank_duration(0);
    start_time = std::chrono::high_resolution_clock::now();
    while (getline(infile, line)) {
        // cout << "Line number " << line_number << endl;
        // if (line_number % 10000 == 0 && upcxx::rank_me() == 0) {

        //     double insert_time = std::chrono::duration<double>(insert_duration).count();
        //     double bfs_time = std::chrono::duration<double>(bfs_duration).count();
        //     cout << "Processed " << line_number << " lines.\nTotal insert duration so far: " << insert_time << "\nTotal BFS duration so far: " << bfs_time << "\n Inserts completed: " << num_inserts << endl;

        // }
        istringstream iss(line);
        string command;
        iss >> command;
        if (line_number % ranks_per_file != upcxx::rank_me() % ranks_per_file) {
        // if (line_number % upcxx::rank_n() != upcxx::rank_me()) {
            if (command == "BFS") {
                num_bfs++;
                end_time = std::chrono::high_resolution_clock::now();
                upcxx::barrier();
                insert_duration += (end_time - start_time);
                start_time = std::chrono::high_resolution_clock::now();
                uint32_t vertex;
                iss >> vertex;
                vector<int> distances = bfs(pcsr, vertex);
                end_time = std::chrono::high_resolution_clock::now();
                bfs_duration += (end_time - start_time) / pcsr.num_elements();
                upcxx::barrier();
                start_time = std::chrono::high_resolution_clock::now();
            }
            line_number++;
            continue;
        }
        if (command == "BFS") {
            num_bfs++;
            end_time = std::chrono::high_resolution_clock::now();
            upcxx::barrier();
            insert_duration += (end_time - start_time);
            start_time = std::chrono::high_resolution_clock::now();
            uint32_t vertex;
            iss >> vertex;
            vector<int> distances = bfs(pcsr, vertex);
            end_time = std::chrono::high_resolution_clock::now();
            bfs_duration += (end_time - start_time) / pcsr.num_elements();
            upcxx::barrier();
            start_time = std::chrono::high_resolution_clock::now();
        } else {
            num_inserts++;
            uint64_t edge = stoull(command);
            uint32_t u, v;
            tie(u, v) = DistPCSR::get_edge_tuple(edge);
            // start_time = std::chrono::high_resolution_clock::now();
            // cout << "Inserting edge " << u << " " << v << endl;
            pcsr.insert_edge(u, v);
            // end_time = std::chrono::high_resolution_clock::now();
            // insert_duration += (end_time - start_time) / pcsr.num_elements;
        }
        line_number++;
    }
    end_time = std::chrono::high_resolution_clock::now();
    insert_duration += end_time - start_time;
    double insert_time = std::chrono::duration<double>(insert_duration).count();
    double bfs_time = std::chrono::duration<double>(bfs_duration).count();
    uint32_t num_elements = pcsr.num_elements();
    upcxx::barrier();
    
    double inserts_average_time = upcxx::reduce_one(insert_time, upcxx::op_fast_add, 0).wait() / upcxx::rank_n();
    uint32_t total_inserts = upcxx::reduce_one(num_inserts, upcxx::op_fast_add, 0).wait();
    double inserts_throughput = total_inserts / inserts_average_time;
    int total_bfs_calls = upcxx::reduce_one(num_bfs, upcxx::op_fast_add, 0).wait();
    double bfs_average_time = upcxx::reduce_one(bfs_time, upcxx::op_fast_add, 0).wait() / total_bfs_calls * pcsr.num_elements();
    if (upcxx::rank_me() == 0) {
        cout << "Insert time: " << inserts_average_time << endl;
        cout << "Insert throughput: " << inserts_throughput << endl;
        cout << "bfs average time scaled to full graph size: " << bfs_average_time << endl;
    }


    // uint32_t total_num_inserts = upcxx::reduce_one(num_inserts, upcxx::op_fast_add, 0).wait();
    // uint32_t total_num_elements = upcxx::reduce_one(num_elements, upcxx::op_fast_add, 0).wait();
    // cout << "pma total inserts: " << total_num_inserts << endl; 
    // cout << "pma total elements: " << total_num_elements << endl; 
    // cout << "Inserts total time: " << upcxx::reduce_all(insert_time, upcxx::op_fast_add).wait() / upcxx::rank_n() << endl;
    // cout << "Queries total time: " << upcxx::reduce_all(query_time, upcxx::op_fast_add).wait() / upcxx::rank_n() << endl;
    // cout << "BFS total time: " << upcxx::reduce_all(bfs_time, upcxx::op_fast_add).wait() / upcxx::rank_n() << endl;
    // cout << "PageRank total time: " << upcxx::reduce_all(pagerank_time, upcxx::op_fast_add).wait() / upcxx::rank_n() << endl;

    infile.close();
    outfile.close();
    // cout << "Rank " << upcxx::rank_me() << " finished" << endl;
    upcxx::barrier();
    start_time = std::chrono::high_resolution_clock::now();
    vector<double> pr_values = pagerank(pcsr);
    end_time = std::chrono::high_resolution_clock::now();
    upcxx::barrier();
    if (upcxx::rank_me() == 0) {
        cout << "PageRank total time: " << std::chrono::duration<double>(end_time - start_time).count() << endl;
    }
    upcxx::finalize();
    return 0;
}

vector<int> bfs_serial(DistPCSR &pcsr, uint32_t source) {
    // if (upcxx::rank_me() == 0) cout << "BFS Start" << endl;
    vector<int> distances(pcsr.num_vertices, -1);
    deque<uint32_t> queue;

    distances[source] = 0;
    queue.push_back(source);

    while (!queue.empty()) {
        uint32_t current = queue.front();
        queue.pop_front();

        vector<uint32_t> neighbors = pcsr.edges(current);
        for (uint32_t neighbor : neighbors) {
            if (distances[neighbor] == -1) {
                distances[neighbor] = distances[current] + 1;
                queue.push_back(neighbor);
            }
        }
    }
    // if (upcxx::rank_me() == 0) cout << "BFS Finish" << endl;
    return distances;
}

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
            // if (vertex == source) {
                // cout << "Degree of source vertex is: " << neighbors.size() << endl;
            // }
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

        upcxx::future<> all_futures = upcxx::make_future();

        for (int target_rank = 0; target_rank < upcxx::rank_n(); target_rank++) {
            if (target_rank == upcxx::rank_me()) {
                (*contributions_recv).insert((*contributions_recv).end(), contributions_destinations[target_rank].begin(), contributions_destinations[target_rank].end());
            } else {
                all_futures = upcxx::when_all(all_futures, 
                upcxx::rpc(
                    target_rank, [](auto& contributions_recv, const vector<pair<uint32_t, double>>& contributions) {
                        (*contributions_recv).insert((*contributions_recv).end(), contributions.begin(), contributions.end());
                    },
                    contributions_recv, contributions_destinations[target_rank] // Should I capture this instead?
                )); // Join these together then wait on all at once if we feel like 
            }
        }

        all_futures.wait();

        upcxx::barrier();

        for (pair<uint32_t, double> contribution_pair : (*contributions_recv)) {
            uint32_t vertex = contribution_pair.first;
            double contribution = contribution_pair.second;
            local_pagerank[vertex - start_vertex] += contribution * damping_factor;
        }
        iteration++;
    }
    return local_pagerank;
}
