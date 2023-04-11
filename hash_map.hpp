#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>
#include "butil.hpp"

#define PROC_OF(x) ((x) % num_procs)
#define BATCH_LEN 32 // max is 32

struct kmer_batch {
    kmer_pair kmers[BATCH_LEN];
    int len = 0;
};

// C++ globals are rpc-accessible
kmer_pair* data; // own
int* used; // own
size_t my_size;
size_t my_rank;
size_t num_procs;

struct HashMap {
    upcxx::global_ptr<kmer_pair>* global_datas;
    upcxx::global_ptr<int>* global_useds;

    kmer_batch* outgoing;

    size_t size() const noexcept;

    HashMap(size_t size);

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    void insert(const kmer_pair& kmer);
    upcxx::future<> pending_inserts = upcxx::make_future();
    void finish_inserts();
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);

    // Helper functions

    // Write and read to a logical data slot in the table.
    void write_slot(uint64_t slot, const kmer_pair& kmer);
    kmer_pair read_slot(uint64_t slot);
    upcxx::future<kmer_pair> async_read_remote_slot(uint64_t slot, int proc);

    // Request a slot or check if it's already used.
    bool request_slot(uint64_t slot);
    bool slot_used(uint64_t slot);

    // non-blocking
    upcxx::future<bool> async_find(const pkmer_t& key_kmer, kmer_pair& val_kmer);
    upcxx::future<bool> async_remote_find(const pkmer_t& key_kmer, kmer_pair& val_kmer, uint64_t probe=0);
};

HashMap::HashMap(size_t size) {
    my_size = size;
    my_rank = upcxx::rank_me();
    num_procs = upcxx::rank_n();
    auto global_data = upcxx::new_array<kmer_pair>(size);
    auto global_used = upcxx::new_array<int>(size);
    data = global_data.local();
    used = global_used.local();
    global_datas = new upcxx::global_ptr<kmer_pair>[num_procs];
    global_useds = new upcxx::global_ptr<int>[num_procs];
    outgoing = new kmer_batch[num_procs];
    for (int i = 0; i < num_procs; ++i) {
        global_datas[i] = upcxx::broadcast(global_data, i).wait();
        global_useds[i] = upcxx::broadcast(global_used, i).wait();
    }
}

bool local_insert(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    for (uint64_t probe = 0; probe < my_size; probe++) {
        uint64_t slot = (hash + probe) % my_size;
        if (used[slot] == 0) {
            used[slot] = 1;
            data[slot] = kmer;
            return true;
        }
    }
    return false;
}

bool batch_insert(const kmer_batch& batch) {
    for (int i = 0; i < batch.len; i++) {
        if (!local_insert(batch.kmers[i])) {
            return false;
        }
    }
    return true;
}

void HashMap::insert(const kmer_pair& kmer) {
    int proc = PROC_OF(kmer.hash());
    kmer_batch& batch = outgoing[proc];
    batch.kmers[batch.len++] = kmer;
    if (batch.len == BATCH_LEN) {
        auto fut = upcxx::rpc(proc, batch_insert, batch).then([](bool success) {
            if (!success) {
                throw std::runtime_error("Error: HashMap is full!");
            }
        });
        pending_inserts = upcxx::when_all(pending_inserts, fut);
        batch.len = 0;
    }
}

void HashMap::finish_inserts() {
    for (int proc = 0; proc < num_procs; proc++) {
        if (outgoing[proc].len) {
            auto fut = upcxx::rpc(proc, batch_insert, outgoing[proc]).then([](bool success) {
                if (!success) {
                    throw std::runtime_error("Error: HashMap is full!");
                }
            });
            pending_inserts = upcxx::when_all(pending_inserts, fut);
        }
    }
    pending_inserts.wait();
}

bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    uint64_t hash = key_kmer.hash();
    if (PROC_OF(hash) == my_rank) {
        for (int probe = 0; probe < my_size; probe++) {
            uint64_t slot = (hash + probe) % my_size;
            val_kmer = read_slot(slot);
            if (val_kmer.kmer == key_kmer) {
                return true;
            }
        }
        return false;
    }
    for (int probe = 0; probe < my_size; probe++) {
        uint64_t slot = (hash + probe) % my_size;
        val_kmer = async_read_remote_slot(slot, hash % num_procs).wait();
        if (val_kmer.kmer == key_kmer) {
            return true;
        }
    }
    return false;
}

upcxx::future<bool> HashMap::async_remote_find(const pkmer_t& key_kmer, kmer_pair& val_kmer, uint64_t probe) {
    if (val_kmer.kmer == key_kmer) {
        return upcxx::make_future(true);
    }
    if (probe >= size()) {
        throw std::runtime_error("===========----------==REMOTE FIND FAIL");
        return upcxx::make_future(false);
    }
    uint64_t hash = key_kmer.hash();
    uint64_t slot = (hash + probe++) % size();
    return async_read_remote_slot(slot, hash % num_procs).then(
        [this, &key_kmer, &val_kmer, probe] (kmer_pair candidate) {
            val_kmer = candidate;
            return async_remote_find(key_kmer, val_kmer, probe);
        }
    );
}

upcxx::future<bool> HashMap::async_find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    uint64_t hash = key_kmer.hash();
    if (PROC_OF(hash) == my_rank) {
        for (int probe = 0; probe < my_size; probe++) {
            uint64_t slot = (hash + probe) % my_size;
            val_kmer = read_slot(slot);
            if (val_kmer.kmer == key_kmer) {
                return upcxx::make_future(true);
            }
        }
        throw std::runtime_error("===========----------==LOCAL FIND FAIL");
        return upcxx::make_future(false);
    }
    return async_remote_find(key_kmer, val_kmer);
}

bool HashMap::slot_used(uint64_t slot) { return used[slot] != 0; }

void HashMap::write_slot(uint64_t slot, const kmer_pair& kmer) { data[slot] = kmer; }

kmer_pair HashMap::read_slot(uint64_t slot) { return data[slot]; }

upcxx::future<kmer_pair> HashMap::async_read_remote_slot(uint64_t slot, int proc) { return rget(global_datas[proc] + slot); }

bool HashMap::request_slot(uint64_t slot) {
    if (used[slot] != 0) {
        return false;
    } else {
        used[slot] = 1;
        return true;
    }
}

size_t HashMap::size() const noexcept { return my_size; }
