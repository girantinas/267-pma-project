#include "utils.hpp"
#include "set_pma.hpp"
#include <deque>
#include <upcxx/upcxx.hpp>
#include <iostream>

typedef tuple<uint32_t, edge_t> range_t; // timestamp, start
using namespace std;
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
        int outstanding_rpcs;
        std::vector<range_t> edge_ranges;
        
};

/* Constructors */
DistPCSR::DistPCSR(int initial_capacity) : pma(initial_capacity) {
    if (upcxx::rank_me() == 0) {
        cout << "# Ranks: " << upcxx::rank_n() << endl;
    }

    for (int rank = 0; rank < upcxx::rank_n(); rank++) {
        edge_ranges.push_back(std::make_pair(0, rank * (10 * UINT32_MAX / upcxx::rank_n())));
        // edge_ranges.push_back(std::make_pair(0, rank * (UINT64_MAX / upcxx::rank_n()))); 
    }
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

    upcxx::rpc(rank,
        [](upcxx::dist_object<DistPCSR>& local_pcsr, edge_t edge) {
            local_pcsr->insert_edge(edge);
        },
        *dist_pcsr_obj, edge
    ).wait();
}

void DistPCSR::insert_edge_local(edge_t edge) {
    if (target_rank(edge) != upcxx::rank_me()) {
        insert_edge(edge);
    } else {
        cout << "inserting edge " << edge << endl;
        pma.insert(edge);
        if (pma.size() >= (uint32_t) pma.get_upper_density_bound(pma.depth()) * pma.capacity()) {
            // TODO: Redistribute
        }
    }
}

void DistPCSR::flush_queue() {
    cout << "flush" << endl;
    while (insert_queue.size() > 0 && !redistributing) {
        edge_t edge = insert_queue.front();
        insert_queue.pop_front();
        insert_edge_local(edge);
    }
}