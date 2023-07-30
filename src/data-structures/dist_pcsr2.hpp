#include "set_pma.hpp"
#include <upcxx/upcxx.hpp>

class DistPCSR {
    public:
        static DistPCSR make_dist_pcsr(int initial_capacity);
        int capacity();
        int size();

        void insert_edge(uint32_t from, uint32_t to);
        upcxx::future<bool> query_edge(uint32_t from, uint32_t to);

    private:
        DistPCSR(int initial_capacity);
        upcxx::dist_object<DistPCSR> *dist_pcsr_obj;
        SetPMA pma;
        bool redistributing;
        int outstanding_rpcs;
};

/* Constructors */
DistPCSR::DistPCSR(int initial_capacity) : pma(initial_capacity) {}

DistPCSR DistPCSR::make_dist_pcsr(int initial_capacity) {
    upcxx::dist_object<DistPCSR> pcsr(initial_capacity);
    pcsr->dist_pcsr_obj = &pcsr;
    return *pcsr;
}

/* Accessors */
int DistPCSR::capacity() {
    return pma.capacity();
}

int DistPCSR::size() {
    return pma.size();
}

/* Gets the edge for a */
uint32_t DistPCSR::target_rank(uint32_t from, uint32_t to) {
    return ;
}

upcxx::future<bool> query_edge(uint32_t from, uint32_t to) {
    
}