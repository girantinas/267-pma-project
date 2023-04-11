#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>

#include "hash_map.hpp"
#include "kmer_t.hpp"
#include "read_kmers.hpp"

#include "butil.hpp"

int main(int argc, char** argv) {
    upcxx::init();

    if (argc < 2) {
        BUtil::print("usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        exit(1);
    }

    std::string kmer_fname = std::string(argv[1]);
    std::string run_type = "";

    if (argc >= 3) {
        run_type = std::string(argv[2]);
    }

    std::string test_prefix = "test";
    if (run_type == "test" && argc >= 4) {
        test_prefix = std::string(argv[3]);
    }

    int ks = kmer_size(kmer_fname);

    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: " + kmer_fname + " contains " + std::to_string(ks) +
                                 "-mers, while this binary is compiled for " +
                                 std::to_string(KMER_LEN) +
                                 "-mers.  Modify packing.hpp and recompile.");
    }

    size_t n_kmers = line_count(kmer_fname);

    // Load factor of 0.5
    size_t hash_table_size = n_kmers * (1.0 / 0.5) / upcxx::rank_n(); 
    HashMap hashmap(hash_table_size);

    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %d for %d kmers.\n", hash_table_size,
                     n_kmers);
    }

    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());

    if (run_type == "verbose") {
        BUtil::print("Finished reading kmers.\n");
    }

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<kmer_pair> start_nodes;

    int progress_time = 2 * upcxx::rank_n();
    int i = 0;
    for (auto& kmer : kmers) {
        hashmap.insert(kmer);
        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
        }
        if (i % progress_time == 0) { upcxx::progress(); }
        ++i;
    }
    hashmap.finish_inserts();
    auto end_insert = std::chrono::high_resolution_clock::now();
    upcxx::barrier();

    double insert_time = std::chrono::duration<double>(end_insert - start).count();
    if (run_type != "test") {
        BUtil::print("Finished inserting in %lf\n", insert_time);
    }
    upcxx::barrier();

    auto start_read = std::chrono::high_resolution_clock::now();

    #define ASSEMBLY_BATCHING 1
    #define TESTING_async_find 0
    #if ASSEMBLY_BATCHING == 0
        std::list<std::list<kmer_pair>> contigs;
        for (const auto& start_kmer : start_nodes) {
            std::list<kmer_pair> contig;
            contig.push_back(start_kmer);
            while (contig.back().forwardExt() != 'F') {
                kmer_pair kmer;
                #if TESTING_async_find == 1
                bool success = hashmap.async_find(contig.back().next_kmer(), kmer).wait();
                #else
                bool success = hashmap.find(contig.back().next_kmer(), kmer);
                #endif
                if (!success) {
                    throw std::runtime_error("Error: k-mer not found in hashmap.");
                }
                contig.push_back(kmer);
            }
            contigs.push_back(contig);
        }
    #elif ASSEMBLY_BATCHING == 1
        // not compatible with async_find
        std::list<std::list<kmer_pair>> contigs;
        upcxx::future<> fut_all_starts = upcxx::make_future();
        i = 0;
        for (const auto& start_kmer : start_nodes) {
            upcxx::future<> fut = upcxx::make_future();
            fut.then([&start_kmer, &contigs, &hashmap]() {
                std::list<kmer_pair> contig;
                contig.push_back(start_kmer);
                while (contig.back().forwardExt() != 'F') {
                    kmer_pair kmer;
                    bool success = hashmap.find(contig.back().next_kmer(), kmer);
                    if (!success) {
                        throw std::runtime_error("Error: k-mer not found in hashmap.");
                    }
                    contig.push_back(kmer);
                }
                contigs.push_back(contig);
            });
            fut_all_starts = upcxx::when_all(fut_all_starts, fut);
            if (i % progress_time == 0) { upcxx::progress(); }
            ++i;
        }
        // conjoin all futures and wait on them
        fut_all_starts.wait();
    #else
        #define WORKERS 10
        std::vector<std::list<kmer_pair>> contigs(start_nodes.size());
        for (auto& start_node : start_nodes) {
            std::list<kmer_pair> contig;
            contig.push_back(start_node);
            contigs.push_back(contig);
        }
        bool done = false;
        std::list<kmer_pair>* growing_contigs[WORKERS];
        upcxx::future<bool> futures[WORKERS];
        kmer_pair incoming[WORKERS];
        for (i = 0; i < WORKERS && i < contigs.size(); i++) {
            growing_contigs[i] = &contigs[i];
            futures[i] = hashmap.async_find(growing_contigs[i]->back().next_kmer(), incoming[i]);
        }
        int next_contig = i;
        if (next_contig == 0) { done = true; }
        for (; i < WORKERS; i++) {
            growing_contigs[i] = NULL;
        }
        while (!done) {
            done = true;
            for (i = 0; i < WORKERS; i++) {
                if (growing_contigs[i] == NULL) {
                    continue;
                }
                if (!futures[i].wait()) {
                    throw std::runtime_error("Error: k-mer not found in hashmap.");
                }
                growing_contigs[i]->push_back(incoming[i]);
                if (growing_contigs[i]->back().forwardExt() != 'F') {
                    futures[i] = hashmap.async_find(growing_contigs[i]->back().next_kmer(), incoming[i]);
                    done = false;
                } else if (next_contig < contigs.size()) {
                    growing_contigs[i] = &contigs[next_contig++];
                    futures[i] = hashmap.async_find(growing_contigs[i]->back().next_kmer(), incoming[i]);
                    done = false;
                } else {
                    growing_contigs[i] = NULL;
                }
                upcxx::progress();
            }
            upcxx::progress();
        }
    #endif

    auto end_read = std::chrono::high_resolution_clock::now();
    upcxx::barrier();
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> read = end_read - start_read;
    std::chrono::duration<double> insert = end_insert - start;
    std::chrono::duration<double> total = end - start;

    int numKmers = std::accumulate(
        contigs.begin(), contigs.end(), 0,
        [](int sum, const std::list<kmer_pair>& contig) { return sum + contig.size(); });

    if (run_type != "test") {
        BUtil::print("Assembled in %lf total\n", total.count());
        BUtil::print("Read in %lf\n", read.count());
    }

    if (run_type == "verbose") {
        printf("Rank %d reconstructed %d contigs with %d nodes from %d start nodes."
               " (%lf read, %lf insert, %lf total)\n",
               upcxx::rank_me(), contigs.size(), numKmers, start_nodes.size(), read.count(),
               insert.count(), total.count());
    }

    if (run_type == "test") {
        std::ofstream fout(test_prefix + "_" + std::to_string(upcxx::rank_me()) + ".dat");
        for (const auto& contig : contigs) {
            fout << extract_contig(contig) << std::endl;
        }
        fout.close();
    }

    upcxx::finalize();
    return 0;
}
