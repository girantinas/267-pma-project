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
#include <cmath>

#include "dist_pcsr.hpp"
// #include "butil.hpp"

using namespace std;

vector<int> bfs(upcxx::dist_object<DistPCSR> &pcsr, uint32_t source);
vector<double> pagerank(upcxx::dist_object<DistPCSR> &pcsr);

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

vector<int> bfs_serial(upcxx::dist_object<DistPCSR> &pcsr, uint32_t source);

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
    infile.open("../rmat-tests/rmat-inserts-shuffled-" + std::to_string(upcxx::rank_me() / ranks_per_file) + ".txt");
    if (!infile.is_open()) {
        cerr << "Couldn't open LiveJournal files" << endl;
        return -1;
    }

    string line;

    if (upcxx::rank_me() == 0) cout << "(rank 0) finished setting up file I/O at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
    int num_edges = 1 << 29;
    int num_ranks = upcxx::rank_n();

    upcxx::dist_object<DistPCSR> pcsr(DistPCSR((num_edges / num_ranks) * 4, (1 << 23) + 1000));
    if (upcxx::rank_me() == 0) cout << "pma depth" << pcsr->spma.depth() << endl;
    if (upcxx::rank_me() == 0) cout << "pma size" << pcsr->spma.size() << endl; 
    if (upcxx::rank_me() == 0) cout << "finished distpcsr construction at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
    upcxx::barrier();
    auto commands_start = std::chrono::high_resolution_clock::now();
    if (upcxx::rank_me() == 0) cout << "starting commands at " << std::chrono::duration<double>(commands_start - start).count() << endl;
    int line_number = 0;
    upcxx::future<> my_futures;
    bool inserting = false;
    int num_inserts = 0;
    chrono::high_resolution_clock::time_point start_time, end_time;
    std::chrono::duration<double> insert_duration(0), query_duration(0), bfs_duration(0), pagerank_duration(0);
    start_time = std::chrono::high_resolution_clock::now();
    int total_inserts_local = 0;
    while (getline(infile, line)) {
        istringstream iss(line);
        string command;
        iss >> command;
        if (line_number % ranks_per_file != upcxx::rank_me() % ranks_per_file) {
        } else if (command == "BFS") {
        } else {
            num_inserts++;
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
    clear_batch_buffers(pcsr);
    sync_dist_pcsr(pcsr);
    end_time = std::chrono::high_resolution_clock::now();
    insert_duration += end_time - start_time;
    upcxx::barrier();
    uint64_t total_inserts = upcxx::reduce_one(total_inserts_local, upcxx::op_fast_add, 0).wait();
    uint64_t total_elements = upcxx::reduce_one(pcsr->num_elements(), upcxx::op_fast_add, 0).wait();
    if (upcxx::rank_me() == 0) {
        cout << "pma total inserts: " << total_inserts << endl;
        cout << "pma total elements: " << total_elements << endl;
        cout << "all inserts finished at " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << endl;
        cout << "time spent: " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - commands_start).count() << endl;
    }


    uint32_t vertex = 0;
    upcxx::barrier();
    start_time = std::chrono::high_resolution_clock::now();
    vector<int> distances = bfs(pcsr, vertex);
    end_time = std::chrono::high_resolution_clock::now();
    upcxx::barrier();
    
    bfs_duration += (end_time - start_time);
    double insert_time = std::chrono::duration<double>(insert_duration).count();
    double bfs_time = std::chrono::duration<double>(bfs_duration).count();
    uint32_t num_elements = pcsr->num_elements();
    upcxx::barrier();

    double inserts_average_time = upcxx::reduce_one(insert_time, upcxx::op_fast_add, 0).wait() / upcxx::rank_n();
    double inserts_throughput = total_inserts / inserts_average_time;

    if (upcxx::rank_me() == 0) {
        cout << "Insert time: " << inserts_average_time << endl;
        cout << "Insert throughput: " << inserts_throughput << endl;
        cout << "bfs average time: " << bfs_time << endl;
    }

    infile.close();
    
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

vector<uint32_t> make_vertex_ranges(upcxx::dist_object<DistPCSR> &pcsr) {
    vector<uint32_t> vertex_ranges;
    uint32_t next_range_start = 0;
    for (tuple<uint32_t, uint64_t, uint64_t> range : pcsr->cached_ranges) {
        pair<uint32_t, uint32_t> range_end = DistPCSR::get_edge_tuple(std::get<2>(range));
        vertex_ranges.push_back(next_range_start);
        next_range_start = range_end.first + 1;
    }
    vertex_ranges.push_back(next_range_start);
    return vertex_ranges;
}

uint32_t vertex_target_rank(vector<uint32_t>& vertex_ranges, uint32_t v) {
    auto target = std::lower_bound(
        vertex_ranges.begin(),
        vertex_ranges.end(),
        v
    );
    uint32_t index;
    if (target != vertex_ranges.end() && *target == v) {
        index = target - vertex_ranges.begin();
    } else {
        if (target == vertex_ranges.begin()) {
            // if this is smaller than all known intervals, send to rank 0
            index = 0;
        }
        else {
            index = target - vertex_ranges.begin() - 1;
        }
    }
    return index;
}

vector<int> bfs(upcxx::dist_object<DistPCSR> &pcsr, uint32_t source) {
    vector<uint32_t> vertex_ranges = make_vertex_ranges(pcsr);
    uint32_t start_vertex = vertex_ranges[upcxx::rank_me()];
    uint32_t end_vertex = vertex_ranges[upcxx::rank_me() + 1];

    uint32_t local_num_vertices = end_vertex - start_vertex;

    vector<int> local_distances(local_num_vertices, -1);

    if (source >= start_vertex && source < end_vertex) {
        local_distances[source - start_vertex] = 0;
    }

    deque<uint32_t> frontier;
    /*
    if (upcxx::rank_me() == 0) {
        cout << "target rank for source vertex: " << vertex_target_rank(vertex_ranges, source) << endl;
        cout << "vertex ranges: " << endl;
        for (int i = 0; i < upcxx::rank_n(); i += 1) {
            cout << "proc " << i << " [" << vertex_ranges[i] << ", " << vertex_ranges[i+1] << ")" << endl;
        }
    }
    */
    if (upcxx::rank_me() == vertex_target_rank(vertex_ranges, source)) {
        frontier.push_back(source);
    }
    int level = 0;

    upcxx::dist_object<vector<uint32_t>> next_set_recv{vector<uint32_t>()};
    vector<vector<uint32_t>> vertex_destinations;
    vertex_destinations.resize(upcxx::rank_n());
    
    unordered_set<uint32_t> next_set;

    upcxx::barrier();

    while (true) {
        if (upcxx::rank_me() == 0) cout << "level " << level << endl;
        next_set.clear();
        bool is_frontier_empty = frontier.empty();
        bool all_frontiers_empty = upcxx::reduce_all(is_frontier_empty, upcxx::op_fast_bit_and).wait();
        if (all_frontiers_empty) {
            break;
        }

        upcxx::future<> all_futures = upcxx::make_future();
        for (uint32_t vertex : frontier) {
            upcxx::future<> fut = add_neighbors_to_set(pcsr, vertex, next_set);
            all_futures = upcxx::when_all(all_futures, fut);
        }   
        all_futures.wait();

        for (int i = 0; i < upcxx::rank_n(); i += 1) {
            vertex_destinations[i].clear();
        }

        for (uint32_t next_vertex : next_set) {
            int target_rank = vertex_target_rank(vertex_ranges, next_vertex);
            vertex_destinations[target_rank].push_back(next_vertex);
        }

        (*next_set_recv).clear();

        upcxx::barrier();

        all_futures = upcxx::make_future();

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
        frontier.clear();

        for (uint32_t vertex : *next_set_recv) {
            int val = local_distances[vertex - start_vertex];
            if (local_distances[vertex - start_vertex] == -1) {
                local_distances[vertex - start_vertex] = level + 1;
                frontier.push_back(vertex);
            }
        }

        level++;
    }
    return local_distances;
}

// requires num_vertices to be correct, i.e. # vertices is correctly known before inserts
vector<double> pagerank(upcxx::dist_object<DistPCSR> &pcsr) {
    double epsilon = 0.0000001;
    double damping_factor = 0.85;
    int max_iterations = 100;
    vector<uint32_t> vertex_ranges = make_vertex_ranges(pcsr);
    uint32_t start_vertex = vertex_ranges[upcxx::rank_me()];
    uint32_t end_vertex = vertex_ranges[upcxx::rank_me() + 1];

    uint32_t local_num_vertices = end_vertex - start_vertex;

    vector<double> local_pagerank(local_num_vertices, 1.0 / pcsr->num_vertices);
    vector<double> local_prev_pagerank(local_num_vertices, 0.0);

    upcxx::dist_object<vector<pair<uint32_t, double>>> contributions_recv{vector<pair<uint32_t, double>>()};

    int iteration = 0;

    upcxx::barrier();

    vector<vector<pair<uint32_t, double>>> contributions_destinations;
    contributions_destinations.resize(upcxx::rank_n());

    while (iteration < max_iterations) {
        if (upcxx::rank_me() == 0) { cout << "pagerank iteration " << iteration << endl; }
        double l1_norm = 0;
        for (int i = 0; i < local_pagerank.size(); i++) {
            l1_norm += std::abs(local_pagerank[i] - local_prev_pagerank[i]);
        }
        double all_l1_norm = upcxx::reduce_all(l1_norm, upcxx::op_fast_add).wait();
        if (all_l1_norm < epsilon) {
            break;
        }

        vector<double> tmp = std::move(local_prev_pagerank);
        local_prev_pagerank = std::move(local_pagerank);
        local_pagerank = std::move(tmp);
        
        local_pagerank.assign(local_num_vertices, (1.0 - damping_factor) / pcsr->num_vertices);

        for (int i = 0; i < upcxx::rank_n(); i += 1) {
            contributions_destinations[i].clear();
        }
        

        for (uint32_t vertex = start_vertex; vertex < end_vertex; vertex++) {
            vector<uint32_t> out_edges;
            edges(pcsr, vertex, out_edges).wait();
            double contribution = local_prev_pagerank[vertex - start_vertex] / out_edges.size();
            for (uint32_t neighbor : out_edges) {
                int target_rank = vertex_target_rank(vertex_ranges, neighbor);
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
