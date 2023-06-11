#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
// #include <upcxx/upcxx.hpp>
#include <vector>
#include <fstream> 
#include <sstream>
#include <random>

#include "set_pma.hpp"
// #include "butil.hpp"

using namespace std;

mt19937 rng(4);

bool insert(SetPMA& pma, set<uint64_t>& reference) {
    uint64_t value = rng();
    bool resized = pma.insert(value);
    reference.insert(value);
    return resized;
}

uint64_t find_value_in_set(set<uint64_t>& reference) {
    assert(reference.size() > 0);
    auto it = reference.end();
    while (it == reference.end()) {
        uint64_t larger = rng();
        it = reference.lower_bound(larger);
    }

    uint64_t value_in_set = *it;
    return value_in_set;
}

uint64_t find_value_not_in_set(set<uint64_t>& reference) {
    uint64_t value_not_in_set;
    auto it = reference.begin();
    while (it != reference.end()) {
        value_not_in_set = rng();
        it = reference.find(value_not_in_set);
    }

    return value_not_in_set;
}

void test_query_hit(SetPMA& pma, set<uint64_t>& reference) {
    uint64_t value_in_set = find_value_in_set(reference);
    assert(reference.find(value_in_set) != reference.end());
    assert(pma.query(value_in_set));
}

void test_query_miss(SetPMA& pma, set<uint64_t>& reference) {
    uint64_t value_not_in_set = find_value_not_in_set(reference);
    assert(reference.find(value_not_in_set) == reference.end());
    assert(!pma.query(value_not_in_set));
}

void test_range_sum(SetPMA& pma, set<uint64_t>& reference) {
    uint64_t range_start = 0;
    uint64_t range_end = -2;
    if (range_start > range_end) {
        std::swap(range_start, range_end);
    }

    uint64_t pmaResult = pma.range_sum(range_start, range_end);
    auto set_it_start = reference.lower_bound(range_start);
    auto set_it_end = reference.upper_bound(range_end);
    uint64_t setSum = 0;
    for (auto it = set_it_start; it != set_it_end; it++) {
        setSum += *it;
    }

    if (pmaResult != setSum) {
        cerr << "Mismatch in range query: PMA returned " << pmaResult
                << " but ordered set returned " << setSum << endl;
        assert(false);
    }
}

void test_contents(SetPMA& pma, set<uint64_t>& reference) {
    vector<uint64_t> all_values = pma.get_min_range(0, pma.size());
    assert(all_values.size() == reference.size()); 
    assert(std::equal(all_values.begin(), all_values.end(), reference.begin(), reference.end()));
}

int main(int argc, char** argv) {
    const int NUM_ITERATIONS = 5e6;
    
    std::uniform_real_distribution<> actionDis(0, 1);
    std::uniform_real_distribution<> existsDis(0, 1);
    
    SetPMA pma(1 << 8);
    set<uint64_t> reference;
    
    for (int i = 0 ; i < NUM_ITERATIONS; i++) {
        double action = actionDis(rng);
        if (action < 0.6) { // Insert with probability 0.6
            bool resized = insert(pma, reference);
            if (resized) {
                cout << "Testing contents after resize" << endl;
                test_range_sum(pma, reference);
                test_contents(pma, reference);
            }
            
        } else if (action < 0.8) { // Point query (hit) with probability 0.1
            if (existsDis(rng) < 0.5 && !reference.empty()) {
                test_query_hit(pma, reference);
            } else { // Point query (miss) with probability 0.1
                test_query_miss(pma, reference);
            }
        } else { // Range query with probability 0.2
            test_range_sum(pma, reference);
        }

        if (i % 10000 == 0) {
            cout << i << endl;
        }
    }
    
    test_contents(pma, reference);
}