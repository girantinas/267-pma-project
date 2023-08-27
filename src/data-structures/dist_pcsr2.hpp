#include "utils.hpp"
#include "set_pma.hpp"
#include <deque>
#include <upcxx/upcxx.hpp>
#include <iostream>

typedef tuple<uint32_t, edge_t> range_t; // timestamp, start
using namespace std;


struct Command {
    uint64_t start_index;
    uint64_t end_index;
    uint32_t target_rank;
};
        
class DistPCSR {
    public:
        static upcxx::dist_object<DistPCSR> make_dist_pcsr(int initial_capacity);
        int capacity();
        int size();
        int depth();

        void insert_edge(edge_t edge);
        void insert_edge(vertex_t from, vertex_t to);

        upcxx::future<bool> query_edge(edge_t edge);
        upcxx::future<bool> query_edge(vertex_t from, vertex_t to);

        void flush_queue();
    private:
        DistPCSR(int initial_capacity);
        uint32_t target_rank(edge_t edge);
        uint32_t target_rank(vertex_t from, vertex_t to);

        // Gets the upper density bound for a subtree of servers
        double get_upper_density_bound(int level);
        double _root_density = 0.67;

        // Processing insert queue
        std::deque<edge_t> insert_queue;
        void insert_edge_local(edge_t edge);
        
        upcxx::dist_object<DistPCSR> *dist_pcsr_obj;
        SetPMA pma;
        bool redistributing = false;
        int outstanding_inserts;
        
        pair<uint64_t, uint64_t> my_range;
        std::vector<range_t> edge_ranges;

        /* Redistribution functions */
        upcxx::global_ptr<int64_t> redis_lock; 
        upcxx::atomic_domain<int64_t> ad;
        void acquire_redis_lock();
        void release_redis_lock();
        void redistribute();
        uint32_t redis_count(uint32_t start_proc, uint32_t num_procs_contact, vector<uint64_t>& element_counts);
        std::tuple<uint64_t, uint32_t, uint32_t> gather_redis_counts(vector<uint64_t>& element_counts);
        void send_all_commands(uint32_t lowest_proc, vector<vector<Command>>& command_lists);
        vector<vector<Command>> get_redistribute_commands(vector<uint64_t>& element_counts, uint32_t start_proc, uint32_t num_procs, uint64_t total_elements);
        upcxx::future<> retrieve_command(Command command);
        void swap_all_data(int lowest_proc, int num_procs, vector<range_t>& updated_ranges);
};



/* Constructors */                                              // root density
DistPCSR::DistPCSR(int initial_capacity) : pma(initial_capacity, 0.75, false) {
    if (upcxx::rank_me() == 0) {
        cout << "# Ranks: " << upcxx::rank_n() << endl;
    }

    for (int rank = 0; rank < upcxx::rank_n(); rank++) {
        // edge_ranges.push_back(std::make_pair(0, rank * (10 * UINT32_MAX / upcxx::rank_n()))); // FOR TESTING
        edge_ranges.push_back(std::make_pair(0, rank * (UINT64_MAX / upcxx::rank_n()))); 
    }

    if (upcxx::rank_me() == 0) {
        redis_lock = upcxx::new_<int64_t>(0);
    }
    redis_lock = upcxx::broadcast(redis_lock, 0).wait();
    ad = upcxx::atomic_domain<int64_t>({upcxx::atomic_op::compare_exchange, upcxx::atomic_op::store});
}

upcxx::dist_object<DistPCSR> DistPCSR::make_dist_pcsr(int initial_capacity) {
    upcxx::dist_object<DistPCSR> pcsr(initial_capacity);
    pcsr->dist_pcsr_obj = &pcsr;
    return pcsr;
}

/* Accessors */
int DistPCSR::capacity() {
    return pma.capacity();
}

int DistPCSR::size() {
    return pma.size();
}

int DistPCSR::depth() {
    return MSSB(upcxx::rank_n());
}

double DistPCSR::get_upper_density_bound(int level) {
    double max_upper_limit = pma.get_upper_density_bound(pma.depth());
    return max_upper_limit - (double) level / depth() * (max_upper_limit - _root_density); // mx + b
}

/* Gets the edge for a */
uint32_t DistPCSR::target_rank(vertex_t from, vertex_t to) {
    edge_t edge = make_edge_tuple(from, to);
    return target_rank(edge);
}

uint32_t DistPCSR::target_rank(edge_t edge) {
    int current_argmin = 0;
    while(current_argmin < edge_ranges.size() && !(get<1>(edge_ranges[current_argmin]) > edge)) {
        current_argmin++;
    }
    return max((current_argmin - 1), 0); // left of 0 pushed to 0
}

upcxx::future<bool> DistPCSR::query_edge(vertex_t from, vertex_t to) {
    edge_t edge = make_edge_tuple(from, to);
    return query_edge(edge);
}

upcxx::future<bool> DistPCSR::query_edge(edge_t edge) {
    uint32_t rank = target_rank(edge);
    if (rank == upcxx::rank_me()) {
        bool result = pma.query(edge);
        cout << "query is: " << result << endl;
        return upcxx::make_future(result);
    }

    return upcxx::rpc(rank, 
        [](upcxx::dist_object<DistPCSR>& local_pcsr, edge_t e) {
            cout << "rpc running" << endl;
            return local_pcsr->query_edge(e);
        },
        *dist_pcsr_obj, edge
    );
}

void DistPCSR::insert_edge(vertex_t from, vertex_t to) {
    edge_t edge = make_edge_tuple(from, to);
    return insert_edge(edge);
}

void DistPCSR::insert_edge(edge_t edge) {
    uint32_t rank = target_rank(edge);
    if (rank == upcxx::rank_me()) {
        insert_queue.push_back(edge);
        return;
    }

    outstanding_inserts += 1;
    upcxx::rpc(rank,
        [](upcxx::dist_object<DistPCSR>& local_pcsr, edge_t edge) {
            local_pcsr->insert_edge(edge);
        },
        *dist_pcsr_obj, edge
    ).then(
        [this]() {
            outstanding_inserts -= 1;
        }
    );
}

void DistPCSR::insert_edge_local(edge_t edge) {
    if (target_rank(edge) != upcxx::rank_me()) {
        insert_edge(edge);
    } else {
        // cout << "inserting edge " << edge << endl;
        pma.insert(edge);
        if (pma.size() >= pma.get_upper_density_bound(pma.depth()) * pma.capacity()) {
            // cout << "oh no we're too big: size: " << pma.size() << ", density bound: " << pma.get_upper_density_bound(pma.depth()) << endl;
            // cout << "Capacity: " << pma.capacity() << endl;
            // cout << "Depth: " << pma.depth() << endl;
            redistribute();
        }
    }
}

void DistPCSR::flush_queue() {
    // cout << "flush" << endl;
    while (insert_queue.size() > 0 && !redistributing) {
        edge_t edge = insert_queue.front();
        insert_queue.pop_front();
        insert_edge_local(edge);
    }
}

void DistPCSR::acquire_redis_lock() {
    bool got_lock = false;
    while (!got_lock) {
        upcxx::progress();
        upcxx::future<int64_t> lock_req = ad.compare_exchange(redis_lock, 0, 1, std::memory_order_seq_cst);
        while(!lock_req.ready()) {
            // may want to call upcxx::progress()
        }
        got_lock = !lock_req.result();
    }
}

void DistPCSR::release_redis_lock() {
    ad.store(redis_lock, 0, std::memory_order_seq_cst).wait();
}

void DistPCSR::redistribute() {
    redistributing = true;
    acquire_redis_lock();

    cout << "Redistributing" << endl;

    // if someone redistributed with us, redistribute is no longer required (they should have reduced number of elements held in us)
    if (!redistributing) {
        release_redis_lock();
        return;
    }

    // 1. Grab density count of each processor, keep going up the tree until the total number of elements falls below some density bound
    // 2. WE have the total number of elements within each subtree, and we want to redistribute them
    //      - Generate commands/actions for each processor to take to do this
    //      - proc 1: n - 50 items, proc 2: n + 50 items, generate command for processor 1 to steal 50 items off the beginning of proc 2's list
    //      - All-to-all data movement within themselves
    //      - Redistribute leader waits for this to finish, and while all this is happening,
    // 3. Individual processors update their ranges and send it to the redistribute leader, who then redistributes the ranges AMONG US
    //      - Each redistributing processor has a completely updated range by the end of the redistribute, if a query/insert is sent to another processor,
    //      it will have the guarantee of eventually hitting the right destination (?)
    //      - The union of all the ranges of the redistributing processors is the same, therefore ranges are only changed within them. So if something belongs to a range
    //      that is not up-to-date, it will be forwarded once to one of the formerly redistributing processors, which will then correctly route it to the right destination
    //      within one more jump.
    //      
    uint64_t total_elements;
    uint32_t redis_n_procs, lowest_proc;
    vector<uint64_t> element_counts(upcxx::rank_n(), 0);
    tie(total_elements, redis_n_procs, lowest_proc) = gather_redis_counts(element_counts);
    
    cout << "redistributing with [" << lowest_proc << ", " << lowest_proc + redis_n_procs << ")" << endl; 
    
    auto command_lists = get_redistribute_commands(element_counts, lowest_proc, redis_n_procs, total_elements);

    send_all_commands(lowest_proc, command_lists);    
    swap_all_data(lowest_proc, redis_n_procs, edge_ranges);

    upcxx::future<> update_ranges_fut = upcxx::make_future();
    for (uint32_t i = 0; i < upcxx::rank_n(); i += 1) {
        if (i >= lowest_proc && i < lowest_proc + redis_n_procs) {
            continue;
        }
        // pcsr->outstanding_rpcs += 1;
        auto f = upcxx::rpc(i, [](upcxx::dist_object<DistPCSR>& pcsr, const vector<range_t>& updated_ranges){
            for (uint32_t i = 0; i < updated_ranges.size(); i += 1) {
                if (std::get<0>(pcsr->edge_ranges[i]) < std::get<0>(updated_ranges[i])) {
                    pcsr->edge_ranges[i] = updated_ranges[i];
                }
            }
        }, *dist_pcsr_obj, edge_ranges); /*.then( // TODO: Make this not a wait
            [&pcsr]() {
                pcsr->outstanding_rpcs -= 1;
            }
        );*/
        update_ranges_fut = upcxx::when_all(update_ranges_fut, f);
    }

    update_ranges_fut.wait();

    redistributing = false;
    release_redis_lock();
}

uint32_t DistPCSR::redis_count(
    uint32_t start_proc, 
    uint32_t num_procs_contact,
    vector<uint64_t>& element_counts
) {
    uint32_t count = 0;
    
    upcxx::future<> count_future = upcxx::make_future();
    for (uint32_t i = 0; i < num_procs_contact; i += 1) {
        cout << "sp: " << start_proc << endl;
        if (start_proc + i == upcxx::rank_me()) {
            cout << "We're rpc'ing ourselves" << endl;
        }
        upcxx::future<> fut = upcxx::rpc(
            start_proc + i,
            [](upcxx::dist_object<DistPCSR>& pcsr) {
                /*if (pcsr->redistributing) {
                    cerr << "Already redistributing when redistribute requested: " << upcxx::rank_me() << endl;
                    exit(-1);
                }*/
                cout << "asserting redistributing: " << upcxx::rank_me() << endl;
                pcsr->redistributing = true;
                // pcsr->pma.print_pma(redistribute_log);
                return make_pair(pcsr->pma.size(), upcxx::rank_me());
            },
            *dist_pcsr_obj
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

std::tuple<uint64_t, uint32_t, uint32_t> DistPCSR::gather_redis_counts(vector<uint64_t>& element_counts) {
    uint64_t total_elements = pma.size();
    uint32_t redis_n_procs = 1;
    uint32_t lowest_proc = upcxx::rank_me();
    element_counts[upcxx::rank_me()] = total_elements;
    int level = 0;
    
    while (total_elements >= get_upper_density_bound(level) * redis_n_procs * pma.capacity() && redis_n_procs < upcxx::rank_n()) {
        redis_n_procs *= 2;
        uint32_t new_proc = (lowest_proc / redis_n_procs) * redis_n_procs;
        cout << "np:" << new_proc << "lp: " << lowest_proc << endl;
        if (new_proc < lowest_proc) {
            total_elements += redis_count(new_proc, redis_n_procs / 2, element_counts);
        } else { // new_proc == lowest_proc
            cout << "?" << new_proc + redis_n_procs / 2 << endl;
            cout << "rp:" << redis_n_procs << endl;
            total_elements += redis_count(new_proc + redis_n_procs / 2, redis_n_procs / 2, element_counts);
        }
        lowest_proc = new_proc;
        level += 1;
    }

    if (redis_n_procs == upcxx::rank_n() && total_elements > get_upper_density_bound(level) * redis_n_procs * pma.capacity()) {
        cerr << "total elements " << total_elements << endl;
        cerr << "total capacity " << redis_n_procs * pma.capacity() << endl;
        cerr << "PMA too full" << endl;
        exit(0);
    }

    return std::make_tuple(total_elements, redis_n_procs, lowest_proc);
}

vector<vector<Command>> DistPCSR::get_redistribute_commands(
    vector<uint64_t>& element_counts, 
    uint32_t start_proc, 
    uint32_t num_procs, 
    uint64_t total_elements
) {
    
    vector<vector<Command>> command_lists;
    uint32_t elements_per_proc_floor = total_elements / num_procs;
    uint32_t current_local_position = 0;
    uint32_t remote_proc = start_proc;
    for (uint32_t i = 0; i < num_procs; i += 1) {
        uint32_t elements = elements_per_proc_floor + (i < total_elements % num_procs);
        
        vector<Command> command_list;
        while (elements > 0) {
            uint32_t remaining_elems = element_counts[remote_proc] - current_local_position;
            if (remaining_elems == 0) {
                remote_proc += 1;
                current_local_position = 0;
            } else if (remaining_elems >= elements) { // == case handled on next outer loop iteration
                Command command = {current_local_position, current_local_position + elements, remote_proc};
                command_list.push_back(command);
                current_local_position += elements;
                elements = 0;
            } else {
                Command command = {current_local_position, element_counts[remote_proc], remote_proc};
                command_list.push_back(command);
                elements -= remaining_elems;
                remote_proc += 1;
                current_local_position = 0;
            }
        }
        command_lists.push_back(std::move(command_list));
    }

    return command_lists;
}

vector<uint64_t> temp;
void DistPCSR::send_all_commands(uint32_t lowest_proc, vector<vector<Command>>& command_lists) {
    upcxx::future<> all_finished = upcxx::make_future();

    for (int i = 0; i < command_lists.size(); ++i) {
        auto send_command_list = [](upcxx::dist_object<DistPCSR>& pcsr, const vector<Command>& command_list) -> upcxx::future<range_t> {
            temp.clear();
            upcxx::future<> retrieve_future = upcxx::make_future();
            for (Command command : command_list) {
                auto fut = pcsr->retrieve_command(command);
                retrieve_future = upcxx::when_all(retrieve_future, fut);
            }
            return retrieve_future.then(
                [&pcsr]() -> range_t {
                    std::sort(temp.begin(), temp.end());
                    pcsr->my_range = make_pair(temp.front(), temp.back());
                    uint32_t current_version = std::get<0>(pcsr->edge_ranges[upcxx::rank_me()]);
                    pcsr->edge_ranges[upcxx::rank_me()] = (range_t) std::make_tuple(current_version + 1, temp.front());
                    return pcsr->edge_ranges[upcxx::rank_me()];
                }
            );
        };

        upcxx::future<> f = upcxx::rpc(lowest_proc + i, send_command_list, *dist_pcsr_obj, command_lists[i]).then(
            [lowest_proc, i, this] (range_t range) { edge_ranges[lowest_proc + i] = range; }
        );

        all_finished = upcxx::when_all(all_finished, f);
    }

    all_finished.wait();
}

upcxx::future<> DistPCSR::retrieve_command(Command command) {
    uint32_t target_proc = command.target_rank;
    uint64_t remote_start = command.start_index;
    uint64_t remote_end = command.end_index;
    return 
    upcxx::rpc(
        target_proc,
        [remote_start, remote_end](upcxx::dist_object<DistPCSR>& pcsr){
            return pcsr->pma.get_min_range(remote_start, remote_end);
        }, 
        *dist_pcsr_obj)
    .then(
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

void DistPCSR::swap_all_data(int lowest_proc, int num_procs, vector<range_t>& updated_ranges) {
    upcxx::future<> all_team_swapped = upcxx::make_future();
    int team_leader = upcxx::rank_me();
    for (uint32_t i = lowest_proc; i < lowest_proc + num_procs; ++i) {
        auto swap_local = upcxx::rpc(i, [team_leader](upcxx::dist_object<DistPCSR>& pcsr, const vector<range_t>& updated_ranges){
            for (uint32_t i = 0; i < updated_ranges.size(); i += 1) {
                // check if updated_ranges is newer than edge_ranges
                if (std::get<0>(pcsr->edge_ranges[i]) < std::get<0>(updated_ranges[i])) {
                    pcsr->edge_ranges[i] = updated_ranges[i];
                }
            }
            pcsr->pma.swap_data(temp);
            // pcsr->pma.print_pma(redistribute_log);
            cout << "team leader: " << team_leader << endl;
            if (upcxx::rank_me() != team_leader) { // team leader will set their own redistributing flag to false later
                cout << "deasserting redistributing: " << upcxx::rank_me() << endl;
                pcsr->redistributing = false;
            }
        }, *dist_pcsr_obj, updated_ranges);
        all_team_swapped = upcxx::when_all(all_team_swapped, swap_local);
    }
    
    all_team_swapped.wait();
}