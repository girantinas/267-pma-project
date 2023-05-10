#pragma once

#include <upcxx/upcxx.hpp>
// #include "butil.hpp"
#include <atomic>
#include <algorithm>
#include <climits>
#include <cmath>
#include <chrono>
#include <iostream>
#include <vector>
#include <deque>
#include <utility>
#include <functional>
#include <sstream>
#include <unordered_set>
#include "set_pma.hpp"

using namespace std;
using namespace std::chrono;
#define timestamp_endl "    " << (duration_cast<milliseconds>(system_clock::now().time_since_epoch()) - reference_time).count() << endl;
milliseconds reference_time;
ofstream redistribute_log;

struct Insert {
    uint32_t from;
    uint32_t to;
    Insert(uint32_t from, uint32_t to) : from(from), to(to) {}
};
typedef tuple<uint32_t, uint64_t, uint64_t> range_t;
const int BATCH_SIZE = 4096;
class DistPCSR {
    public:
        DistPCSR(uint32_t size, uint64_t max_val);
        uint32_t size();
        uint32_t num_elements(); // number of nonempty elements in this rank's local PCSR
        uint32_t num_vertices; // Restrict this to power of 2
        /* Inserts an edge from one vertex to another. Fails if either vertex does not exist. */
        void insert_edge_local(upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, uint32_t to);
        /* Queries if an edge exists. */
        bool query_edge(uint32_t from, uint32_t to);
        /* Edge list of a vertex. If the vertex does not exist, the list is empty. */
        vector<uint32_t> edges_local(uint32_t from);
        /* Returns the entire adjacency list representation of the graph. */
        vector<pair<uint32_t, vector<uint32_t>>> adjacency_lists();
        /* Prints all attributes of this DistPCSR for debugging */
        void print_dist_pcsr();

        deque<Insert> rq_queue;
        static constexpr double INIT_LEAF_MAX = 0.9;
        static uint64_t make_edge_tuple(uint32_t from, uint32_t to);
        static pair<uint32_t, uint32_t> get_edge_tuple(uint64_t edge);
        uint32_t target_rank(uint32_t from, uint32_t to);
        SetPMA spma;
        const uint32_t TO_ONES = 0xFFFFFFFF;
        void acquire_redis_lock();
        void release_redis_lock();

        pair<uint64_t, uint64_t> my_range;
        vector<range_t> cached_ranges;

        bool redistributing;
        upcxx::global_ptr<int64_t> redis_lock; 
        upcxx::atomic_domain<int64_t> ad;

        int outstanding_rpcs;

        vector<vector<Insert>> batches;
};

std::ostream& operator<<(std::ostream& _os, const range_t& _p) {
    uint32_t timestamp;
    uint64_t low;
    uint64_t high;
    tie(timestamp, low, high) = _p;
    uint32_t low_from, low_to;
    uint32_t high_from, high_to;
    tie(low_from, low_to) = DistPCSR::get_edge_tuple(low);
    tie(high_from, high_to) = DistPCSR::get_edge_tuple(high);
    
    _os << "[(from=" << low_from << ", to=" << low_to << "), (from=" << high_from << ", to=" << high_to << ")) (version=" << timestamp << ")";
    return _os;
}

void insert_edge(upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, uint32_t to);
void insert_edge_batched(upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, uint32_t to);
upcxx::future<bool> query_edge(upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, uint32_t to);
uint32_t redis_count(
    upcxx::dist_object<DistPCSR>& pcsr, 
    uint32_t start_proc, 
    uint32_t num_procs_contact,
    vector<uint32_t>& element_counts
);
vector<pair<uint32_t, vector<uint32_t>>> get_redistribute_commands(
    vector<uint32_t>& element_counts, 
    uint32_t start_proc, 
    uint32_t num_procs, 
    uint32_t total_elements
);
upcxx::future<> retrieve_command(upcxx::dist_object<DistPCSR>& pcsr, tuple<uint32_t, uint32_t, uint32_t> cmd, vector<uint64_t> temp);
void redistribute(upcxx::dist_object<DistPCSR>& pcsr);


uint32_t DistPCSR::size() {
    return spma.size();
}

uint32_t DistPCSR::num_elements() {
    return spma._num_elements;
}

DistPCSR::DistPCSR(uint32_t size, uint64_t max_val) : spma(SetPMA(size, INIT_LEAF_MAX, false)) {
    num_vertices = max_val;
    auto lower_from = (max_val / upcxx::rank_n()) * upcxx::rank_me();
    auto upper_from = (max_val / upcxx::rank_n()) * (upcxx::rank_me() + 1);
    auto lower_tup = make_edge_tuple(lower_from, 0);
    auto upper_tup = make_edge_tuple(upper_from, TO_ONES);
    my_range = make_pair(lower_tup, upper_tup);
    for (uint32_t i = 0; i < upcxx::rank_n(); i += 1) {
        auto lower_from = (max_val / upcxx::rank_n()) * i;
        auto upper_from = (max_val / upcxx::rank_n()) * (i + 1);
        auto lower_tup = make_edge_tuple(lower_from, 0);
        auto upper_tup = make_edge_tuple(upper_from, TO_ONES);
        cached_ranges.push_back(make_tuple(0, lower_tup, upper_tup));
    }

    outstanding_rpcs = 0;
    redistributing = false;
    if (upcxx::rank_me() == 0) {
        redis_lock = upcxx::new_<int64_t>(0);
    }
    redis_lock = upcxx::broadcast(redis_lock, 0).wait();
    ad = upcxx::atomic_domain<int64_t>({upcxx::atomic_op::compare_exchange, upcxx::atomic_op::store});
    
    if (upcxx::rank_me() == 0) {
        reference_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    }
    reference_time = upcxx::broadcast(reference_time, 0).wait();

    batches.resize(upcxx::rank_n());
    for (auto& insert_vec : batches) {
        insert_vec.reserve(BATCH_SIZE);
    }

    // redistribute_log.open("redistribute_" + std::to_string(upcxx::rank_me()) + ".txt");
}

void DistPCSR::acquire_redis_lock() {
    while (ad.compare_exchange(redis_lock, 0, 1, std::memory_order_seq_cst).wait()) {
        upcxx::progress();
    };
}

void DistPCSR::release_redis_lock() {
    ad.store(redis_lock, 0, std::memory_order_seq_cst).wait();
}

uint64_t DistPCSR::make_edge_tuple(uint32_t from, uint32_t to) { return (((uint64_t) from) << 32) | to; }
pair<uint32_t, uint32_t> DistPCSR::get_edge_tuple(uint64_t edge) { return make_pair(static_cast<uint32_t>(edge >> 32), static_cast<uint32_t>(edge)); }

// binary search in our cached values of ranges; could be outdated
uint32_t DistPCSR::target_rank(uint32_t from, uint32_t to) {
    uint64_t edge_tuple = make_edge_tuple(from, to);

    // first check local rank, will always be up to date.
    if (edge_tuple >= my_range.first && edge_tuple < my_range.second) {
        return upcxx::rank_me();
    }

    auto target = std::lower_bound(
        cached_ranges.begin(),
        cached_ranges.end(),
        edge_tuple, 
        [](const range_t& a, const uint64_t b) { return std::get<1>(a) < b; }
    );

    uint32_t index;
    if (target != cached_ranges.end() && std::get<1>(*target) == edge_tuple) {
        index = target - cached_ranges.begin(); // need this because if its actually equal to low end of interval of some range, don't want smaller proc. basically, discrepancy between < and <=
    } else {
        if (target == cached_ranges.begin()) {
            // if this is smaller than all known intervals, send to rank 0
            index = 0;
        }
        else {
            index = target - cached_ranges.begin() - 1;
        }
    }
    return index;
}

void flush_queue(upcxx::dist_object<DistPCSR>& pcsr) {
    if (upcxx::rank_me() == 0) {
        // cout << "queue size " << pcsr->rq_queue.size() << endl;
    }
    while (pcsr->rq_queue.size() != 0 && !pcsr->redistributing) {
        auto insert = pcsr->rq_queue.front();
        pcsr->rq_queue.pop_front();
        pcsr->insert_edge_local(pcsr, insert.from, insert.to);
        if (upcxx::rank_me() == 0) {
            // cout << "queue size " << pcsr->rq_queue.size() << endl;
        }     
    }
}

void DistPCSR::insert_edge_local(upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, uint32_t to) {
    // forwarding
    if (target_rank(from, to) != upcxx::rank_me()) {
        if (upcxx::rank_me() == 0) {
            // cout << "forwarding edge " << from << " " << to << " from " << upcxx::rank_me() << " to " << target_rank(from, to) << endl;
        }
        insert_edge(pcsr, from, to);
    }
    else {
        uint64_t edge = make_edge_tuple(from, to);
        spma.insert(edge);
        if (spma._num_elements > (uint32_t) ((spma._leaf_max - spma.depth() * 0.01) * spma.size())) {
            cout << "rank " << upcxx::rank_me() << " redistributing!" << endl;
            redistribute(pcsr);
        }
    }
}

void insert_edge_batch(upcxx::dist_object<DistPCSR>& pcsr, uint32_t target_rank, const vector<Insert>& edge_batch) {
    pcsr->outstanding_rpcs += 1;
    upcxx::rpc(
        target_rank,
        [](upcxx::dist_object<DistPCSR>& pcsr, const vector<Insert>& edge_batch) {
            pcsr->rq_queue.insert(pcsr->rq_queue.end(), edge_batch.begin(), edge_batch.end());
            // cout << "received edge batch" << endl;
        },
        pcsr, edge_batch
    ).then(
        [&pcsr]() {
            pcsr->outstanding_rpcs -= 1;
        }
    );
}

void insert_edge_batched(upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, uint32_t to) {
    uint32_t rank = pcsr->target_rank(from, to);
    pcsr->batches[rank].push_back(Insert(from, to));
    if (pcsr->batches[rank].size() >= BATCH_SIZE) {
        insert_edge_batch(pcsr, rank, pcsr->batches[rank]);
        pcsr->batches[rank].clear();
    }
}

void clear_batch_buffers(upcxx::dist_object<DistPCSR>& pcsr) {
    for (int rank = 0; rank < pcsr->batches.size(); rank += 1) {
        if (pcsr->batches[rank].size() == 0) {
            continue;
        }
        insert_edge_batch(pcsr, rank, pcsr->batches[rank]);
        pcsr->batches[rank].clear();
    }
}

void insert_edge(upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, uint32_t to) {
    uint32_t rank = pcsr->target_rank(from, to);
    if (rank == upcxx::rank_me()) {
        pcsr->rq_queue.push_back(Insert(from, to));
        return;
    }
    
    pcsr->outstanding_rpcs += 1;
    upcxx::rpc(rank,
        [](upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, uint32_t to) {
            insert_edge(pcsr, from, to);
        },
        pcsr, from, to
    ).then(
        [&pcsr]() {
            pcsr->outstanding_rpcs -= 1;
        }
    );
}

bool DistPCSR::query_edge(uint32_t from, uint32_t to) {
    uint64_t edge = make_edge_tuple(from, to);
    return spma.query(edge);
}

upcxx::future<bool> query_edge(upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, uint32_t to) {
    uint32_t rank = pcsr->target_rank(from, to);
    if (pcsr->redistributing) {
        cout << "why are we query edge during redistribute" << timestamp_endl;
        cout << from << " " << to << timestamp_endl;
        cout << "sending rpc to rank " << rank << " from " << upcxx::rank_me() << timestamp_endl;
        exit(-1);
    }
    if (rank == upcxx::rank_me()) {
        bool result = pcsr->query_edge(from, to);
        if (result == false) {
            /*cout << "critical: rank " << upcxx::rank_me() << " failed to find " << from << " " << to << endl;
            cout << "my_range: " << make_tuple(0, pcsr->my_range.first, pcsr->my_range.second) << endl;
            int i = 0;
            for (auto r : pcsr->cached_ranges) {
                cout << "rank " << i << " cached range: " << r << endl;
                i += 1;
            }*/
        }
        return upcxx::make_future(result);
    }
    
    // cout << upcxx::rank_me() << " sending edge " << from << "," << to << " to " << rank << timestamp_endl;
    // pcsr->print_dist_pcsr();
    
    return upcxx::rpc(rank, 
        [](upcxx::dist_object<DistPCSR>& local_pcsr, uint32_t from, uint32_t to) {
            return query_edge(local_pcsr, from, to);
        },
        pcsr, from, to
    );
}

uint32_t redis_count(
    upcxx::dist_object<DistPCSR>& pcsr, 
    uint32_t start_proc, 
    uint32_t num_procs_contact,
    vector<uint32_t>& element_counts
) {
    uint32_t count = 0;
    
    upcxx::future<> count_future = upcxx::make_future();
    for (uint32_t i = 0; i < num_procs_contact; i += 1) {
        upcxx::future<> fut = upcxx::rpc(
            start_proc + i,
            [](upcxx::dist_object<DistPCSR>& pcsr) {
                pcsr->redistributing = true;
                // pcsr->spma.print_pma(redistribute_log);
                return make_pair(pcsr->spma._num_elements, upcxx::rank_me());
            },
            pcsr
        ).then(
            [&count, &element_counts] (pair<uint32_t, upcxx::intrank_t> element_count) {
                uint32_t num_elements;
                upcxx::intrank_t rank;
                tie(num_elements, rank) = element_count;
                count += num_elements;
                element_counts[rank] = num_elements;
            }
        );
        count_future = upcxx::when_all(count_future, fut);
    }

    count_future.wait();
    return count;
}
// creates a list of (remote start index, remote end index, processor id to retrieve from), one for each retriever.
// for all processors participating in redistribute. each processor is responsible for acquiring data
// it needs from other processors
vector<pair<uint32_t, vector<uint32_t>>> get_redistribute_commands(
    vector<uint32_t>& element_counts, 
    uint32_t start_proc, 
    uint32_t num_procs, 
    uint32_t total_elements
) {
    
    vector<pair<uint32_t, vector<uint32_t>>> command_lists;
    uint32_t elements_per_proc_floor = total_elements / num_procs;
    uint32_t current_local_position = 0;
    uint32_t remote_proc = start_proc;
    for (uint32_t i = 0; i < num_procs; i += 1) {
        uint32_t retrieving_proc = start_proc + i;
        uint32_t elements = elements_per_proc_floor + (i < total_elements % num_procs);
        
        vector<uint32_t> command_list;
        while (elements > 0) {
            uint32_t remaining_elems = element_counts[remote_proc] - current_local_position;
            if (remaining_elems == 0) {
                
                remote_proc += 1;
                current_local_position = 0;
            } else if (remaining_elems >= elements) {
                command_list.push_back(current_local_position);
                command_list.push_back(current_local_position + elements);
                command_list.push_back(remote_proc);
                current_local_position += elements;
                
                elements = 0;
            } else {
                command_list.push_back(current_local_position);
                command_list.push_back(element_counts[remote_proc]);
                command_list.push_back(remote_proc);
                elements -= remaining_elems;
                remote_proc += 1;
                
                current_local_position = 0;
            }
        }
        command_lists.push_back(make_pair(retrieving_proc, std::move(command_list)));
    }

    for (auto& command_list : command_lists) {
        uint32_t retrieving_proc = command_list.first;
        vector<uint32_t>& flattened = command_list.second;
        
        for (int i = 0; i < flattened.size(); i += 3) {
            
        }
    }
    return command_lists;
}

vector<uint64_t> temp;

upcxx::future<> retrieve_command(upcxx::dist_object<DistPCSR>& pcsr, tuple<uint32_t, uint32_t, uint32_t> cmd) {
    uint32_t remote_start = std::get<0>(cmd);
    uint32_t remote_end = std::get<1>(cmd);
    uint32_t target_proc = std::get<2>(cmd);
    
    return upcxx::rpc(target_proc,
        [remote_start, remote_end](upcxx::dist_object<DistPCSR>& pcsr){
            
            return pcsr->spma.get_min_range(remote_start, remote_end);
        }, pcsr).then(
            [target_proc](const vector<uint64_t>& retrieved_values) {
                for (uint64_t retrieved : retrieved_values) {
                    if (retrieved == SetPMA::INT_NULL) {
                        cerr << "proc " << upcxx::rank_me() << " received an INT_NULL in values retrieved from target_proc=" << target_proc << endl;
                        exit(-1);
                    }
                }
                temp.insert(temp.end(), retrieved_values.begin(), retrieved_values.end());
            }
        );
}

void redistribute(upcxx::dist_object<DistPCSR>& pcsr) {
    pcsr->redistributing = true;
    pcsr->acquire_redis_lock();

    // if someone else redistributed with us while we were waiting for lock,
    // making a redistribute unnecessary, then don't redistribute
    if (!pcsr->redistributing) {
        
        pcsr->release_redis_lock();
        return;
    }

    
    uint32_t total_elements = pcsr->spma._num_elements;
    vector<uint32_t> element_counts(upcxx::rank_n(), 0);
    uint32_t n_procs = 1;
    uint32_t start_proc = upcxx::rank_me();
    element_counts[start_proc] = total_elements;
    // need a buffer of at least 1 in each processor (arguably probably more)
    int level = 0;
    while (total_elements > (uint32_t) ((pcsr->spma._leaf_max - pcsr->spma.depth() * 0.01 - level * 0.02) * n_procs * pcsr->size()) && (n_procs < upcxx::rank_n())) {
        n_procs *= 2;
        
        uint32_t new_proc = (start_proc / n_procs) * n_procs;
        if (new_proc < start_proc) {
            total_elements += redis_count(pcsr, new_proc, n_procs / 2, element_counts);
        } else {
            total_elements += redis_count(pcsr, new_proc + n_procs / 2, n_procs / 2, element_counts);
        }
        start_proc = new_proc;
        level += 1;
    }
    
    if (n_procs == upcxx::rank_n() && total_elements >= (uint32_t) ((pcsr->spma._leaf_max - level * 0.01) * n_procs * pcsr->size())) {
        cerr << "total elements " << total_elements << timestamp_endl;
        cerr << "total capacity " << (uint32_t) (pcsr->spma._leaf_max * n_procs * pcsr->size()) << timestamp_endl;
        cerr << "PMA too full" << timestamp_endl;
        exit(0);
    }

    auto commands = get_redistribute_commands(element_counts, start_proc, n_procs, total_elements);
    upcxx::future<> all_team_finished = upcxx::make_future();
    

    for (uint32_t i = 0; i < commands.size(); ++i) {
        uint32_t p;
        vector<uint32_t> flattened_commands;
        tie(p, flattened_commands) = commands[i];
        auto lambda = [](upcxx::dist_object<DistPCSR>& pcsr, const vector<uint32_t>& flattened_commands) -> upcxx::future<range_t> {
            temp.clear();
            upcxx::future<> retrieve_future = upcxx::make_future();
            

            for (uint32_t c = 0; c < flattened_commands.size(); c += 3) {
                
                tuple<uint32_t, uint32_t, uint32_t> command = make_tuple(flattened_commands[c], flattened_commands[c + 1], flattened_commands[c + 2]);
                retrieve_future = upcxx::when_all(retrieve_future, retrieve_command(pcsr, command));
            }

            upcxx::future<range_t> return_future = retrieve_future.then(
                [&pcsr]() -> range_t {
                    std::sort(temp.begin(), temp.end());
                    pcsr->my_range = make_pair(temp.front(), temp.back());
                    uint32_t current_version = std::get<0>(pcsr->cached_ranges[upcxx::rank_me()]);
                    pcsr->cached_ranges[upcxx::rank_me()] = (range_t) make_tuple(current_version + 1, temp.front(), temp.back());
                    

                    return pcsr->cached_ranges[upcxx::rank_me()];
                }
            );

            return return_future;
        };

        upcxx::future<range_t> f = upcxx::rpc(p, lambda, pcsr, flattened_commands);

        all_team_finished = upcxx::when_all(all_team_finished, f.then(
            [p, &pcsr](range_t range) {
                pcsr->cached_ranges[p] = range;
            }
        ));
    }

    all_team_finished.wait();


    upcxx::future<> all_team_swapped = upcxx::make_future();
    for (uint32_t i = 0; i < commands.size(); ++i) {
        uint32_t p;
        std::vector<uint32_t> p_commands;
        std::tie(p, p_commands) = commands[i];
        auto f = upcxx::rpc(p, [root_p=upcxx::rank_me()](upcxx::dist_object<DistPCSR>& pcsr, const vector<range_t>& updated_ranges){
            for (uint32_t i = 0; i < updated_ranges.size(); i += 1) {
                // check if updated_ranges is newer than cached_ranges
                if (std::get<0>(pcsr->cached_ranges[i]) < std::get<0>(updated_ranges[i])) {
                    pcsr->cached_ranges[i] = updated_ranges[i];
                }
                else {

                }
            }
            pcsr->spma.swap_data(temp);
            // pcsr->spma.print_pma(redistribute_log);
            if (upcxx::rank_me() != root_p) { // team leader will set their own redistributing flag to false later
                pcsr->redistributing = false;
            }
        }, pcsr, pcsr->cached_ranges);
        all_team_swapped = upcxx::when_all(all_team_swapped, f);
    }
    
    all_team_swapped.wait();

    
    for (uint32_t i = 0; i < upcxx::rank_n(); i += 1) {
        if (i >= start_proc && i < start_proc + n_procs) {
            continue;
        }
        pcsr->outstanding_rpcs += 1;
        upcxx::rpc(i, [](upcxx::dist_object<DistPCSR>& pcsr, const vector<range_t>& updated_ranges){
            for (uint32_t i = 0; i < updated_ranges.size(); i += 1) {
                if (std::get<0>(pcsr->cached_ranges[i]) < std::get<0>(updated_ranges[i])) {
                    pcsr->cached_ranges[i] = updated_ranges[i];
                }
            }
        }, pcsr, pcsr->cached_ranges).then(
            [&pcsr]() {
                pcsr->outstanding_rpcs -= 1;
            }
        );
    }

    pcsr->redistributing = false;
    pcsr->release_redis_lock();

    // pcsr->spma._leaf_max = min(pcsr->spma._leaf_max + 0.1, 0.8);
}

void DistPCSR::print_dist_pcsr() {
    ofstream outfile;
    outfile.open("final_" + std::to_string(upcxx::rank_me()) + ".txt");
    outfile << "Proc " << upcxx::rank_me() << " ends with " << spma._num_elements << " elements, has range " << \
        get<0>(my_range) << "," << get<1>(my_range) << ")" << timestamp_endl;
    int i = 0;
    for (auto range : cached_ranges) {
        outfile << "proc " << upcxx::rank_me() << " thinks proc " << i << " has range of " << range << timestamp_endl; 
        i += 1;
    }
    
    std::ostringstream ss;
    ss << "Proc " << upcxx::rank_me() << " ending contents: [";
    auto printer = SetPMA::range_func(
        [&ss](uint64_t item) { 
            uint32_t from, to;
            tie(from, to) = DistPCSR::get_edge_tuple(item);
            ss << "(" << from << "," << to << "), "; 
        }
    ); 
    spma.range(0, UINT64_MAX - 1, printer);
    ss << "]";
    outfile << ss.str() << timestamp_endl;
    outfile.close();
    // redistribute_log.close();
}

vector<uint32_t> DistPCSR::edges_local(uint32_t from) {
    vector<uint32_t> edge_list;
    auto edge_adder = SetPMA::range_func([&edge_list](uint64_t e){
        edge_list.push_back(get_edge_tuple(e).second);
    });
    spma.range(make_edge_tuple(from, 0), make_edge_tuple(from, UINT32_MAX), edge_adder);
    return edge_list;
}

upcxx::future<> edges(upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, vector<uint32_t>& dest) {
    uint32_t begin_rank = pcsr->target_rank(from, 0);
    uint32_t end_rank = pcsr->target_rank(from, UINT32_MAX);

    upcxx::future<> finish_future = upcxx::make_future();
    for (int rank = begin_rank; rank <= end_rank; ++rank) {
        if (rank == upcxx::rank_me()) {
            for (uint32_t v : pcsr->edges_local(from)) {
                dest.push_back(v);
            };
        } else {
            auto f = upcxx::rpc(rank, 
                [from](upcxx::dist_object<DistPCSR>& local_pcsr) {
                    return local_pcsr->edges_local(from);
                }, pcsr
            ).then([&dest](vector<uint32_t> v) {
                dest.insert(dest.end(), v.begin(), v.end());
            });
            finish_future = upcxx::when_all(f, finish_future);
        }
    }
    return finish_future;
}

upcxx::future<> add_neighbors_to_set(upcxx::dist_object<DistPCSR>& pcsr, uint32_t from, unordered_set<uint32_t>& dest) {
    uint32_t begin_rank = pcsr->target_rank(from, 0);
    uint32_t end_rank = pcsr->target_rank(from, UINT32_MAX);

    upcxx::future<> finish_future = upcxx::make_future();
    for (int rank = begin_rank; rank <= end_rank; ++rank) {
        if (rank == upcxx::rank_me()) {
            for (uint32_t v : pcsr->edges_local(from)) {
                dest.insert(v);
            };
        } else {
            auto f = upcxx::rpc(rank, 
                [from](upcxx::dist_object<DistPCSR>& local_pcsr) {
                    return local_pcsr->edges_local(from);
                }, pcsr
            ).then([&dest](vector<uint32_t> v) {
                dest.insert(v.begin(), v.end());
            });
            finish_future = upcxx::when_all(f, finish_future);
        }
    }
    return finish_future;
}