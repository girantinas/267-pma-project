#pragma once

// #include <upcxx/upcxx.hpp>
// #include "butil.hpp"
#include <climits>
#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include <cstring>
#include <cassert>
#include <unordered_set>

using namespace std;

class SetPMA {
    public:
        typedef std::function<void(uint64_t)> range_func;
        static constexpr uint64_t INT_NULL = UINT64_MAX;
        static constexpr uint32_t INVALID_IDX = UINT32_MAX;
        SetPMA(uint32_t size);
        SetPMA(uint32_t size, double leaf_max);
        SetPMA(uint32_t size, double leaf_max, bool resize_allowed);
        
        uint32_t capacity();
        uint32_t size();
        SetPMA() = default;

        bool insert(uint64_t i);
        bool query(uint64_t i);
        /* Sum keys in [left, right) */
        uint64_t range_sum(uint64_t left, uint64_t right);
        vector<uint64_t> get_min_range(uint32_t left, uint32_t right);
        void swap_data(std::vector<uint64_t>& temp);
        void range(uint64_t left, uint64_t right, range_func& op);
        double _leaf_max = 0.75;

        std::vector<uint64_t> data;
        uint32_t _capacity;
        uint32_t search(uint64_t i);
        uint32_t logN();
        uint32_t loglogN();
        uint32_t leaf_index(uint32_t index);
        uint32_t next_leaf(uint32_t index);
        uint32_t leaf_number(uint32_t index);
        uint32_t leaf_position(uint32_t leaf_num);
        uint32_t num_leaves();
        uint32_t depth();
        uint32_t count_nonempty(uint32_t index, uint32_t len);
        void gather(uint32_t index, uint32_t len, uint32_t density_count);
        void distribute(uint32_t index, uint32_t len, uint32_t density_count);
        void redistribute(uint32_t index, uint32_t len, uint32_t density_count);
        void resize();
        void slide_right(uint32_t index);
        void print_pma(ostream& stream = cout);
        uint32_t _size;
        bool _resize_allowed;
        uint32_t _max_index;
        vector<uint64_t> _temp;
};

SetPMA::SetPMA(uint32_t size) {
    data.resize(size, INT_NULL);
    _capacity = size;
    _size = 0;
    _resize_allowed = true;
    _max_index = INVALID_IDX;
}

SetPMA::SetPMA(uint32_t size, double leaf_max) {
    data.resize(size, INT_NULL);
    _capacity = size;
    _size = 0;
    _leaf_max = leaf_max;
    _resize_allowed = true;
    _max_index = INVALID_IDX;
}

SetPMA::SetPMA(uint32_t size, double leaf_max, bool resize_allowed) {
    data.resize(size, INT_NULL);
    _capacity = size;
    _size = 0;
    _leaf_max = leaf_max;
    _resize_allowed = resize_allowed;
    _max_index = INVALID_IDX;
}

// Gets the range [left-th minimum, right-th minimum)
// left-th item to the right-th item, not including nulls.
vector<uint64_t> SetPMA::get_min_range(uint32_t left, uint32_t right) {
    uint32_t min_count = 0;
    assert (left < data.size());
    vector<uint64_t> range;
    range.reserve(right - left);

    // Scan through the PMA, save all non-nulls.
    for (int i = 0; i < data.size();) {
        if (data[i] != INT_NULL) {
            if (min_count >= left && min_count < right) {
                range.push_back(data[i]);
            }
            min_count += 1;
            i += 1;
        } else {
            // Reached a null, can skip the rest of the values of the leaf.
            i = next_leaf(i);
        }
        if (min_count >= right) {
            break;
        }
    }
    assert ((right - left) == range.size());
    return range;
}

void SetPMA::swap_data(std::vector<uint64_t>& temp) {
    uint32_t density_count = temp.size();
    temp.resize(_capacity, INT_NULL);
    std::swap(data, temp);
    _size = density_count;
    // cout << "swapping data..." << endl;
    // cout << density_count << endl;
    redistribute(0, _capacity, density_count);
    // cout << _max_index << endl;
    // cout << data[_max_index] << endl;
}

uint32_t MSSB(uint32_t x) {
    uint32_t i = 0;
    while (x != 0) {
        x = x >> 1;
        ++i;
    }
    return i - 1;
}

uint32_t next_power_of_2(uint32_t x) {
    return 1 << (MSSB(x) + 1);
}

uint32_t SetPMA::capacity() { return _capacity; } // N
uint32_t SetPMA::size() { return _size; }
uint32_t SetPMA::logN() { return next_power_of_2((uint32_t) log2(_capacity)); }
uint32_t SetPMA::loglogN() { return (uint32_t) log2(logN()); }
uint32_t SetPMA::leaf_index(uint32_t index) { return (index & ~(logN() - 1)); }
uint32_t SetPMA::next_leaf(uint32_t index) { return leaf_index(index + logN()); }
uint32_t SetPMA::leaf_number(uint32_t index) { return leaf_index(index) >> loglogN(); }
uint32_t SetPMA::leaf_position(uint32_t leaf_num) { return leaf_num << loglogN(); }
uint32_t SetPMA::num_leaves() { return _capacity / logN(); }
uint32_t SetPMA::depth() { return MSSB(num_leaves()); }
uint32_t SetPMA::count_nonempty(uint32_t index, uint32_t len)  { 
    uint32_t full = 0;
    for (uint32_t i = index; i < index + len;) {
        if (data[i] != INT_NULL) {
            ++full;
            i += 1;
        }
        else {
            i = next_leaf(i);
        }
    }
    return full;
}

bool SetPMA::query(uint64_t key) {
    uint32_t idx = search(key);
    if (idx == INVALID_IDX) {
        return false;
    }
    return data[idx] == key;
}

// finds index of key in the PMA, if it is present
// if it is not present, returns index of largest key in PMA smaller than key
uint32_t SetPMA::search(uint64_t key) {
    uint32_t low = 0;
    uint32_t high = leaf_index(_capacity - 1);
    // find minimum value in PMA
    uint64_t min_key = data[low];
    if (key == min_key) {
        return low;
    }
    else if (key < min_key) {
        // smaller than minimum element
        return INVALID_IDX;
    }
    
    uint32_t argmax = _max_index;
    uint64_t max_key = data[argmax];
    high = leaf_index(argmax);
    
    if (key >= max_key) {
        return argmax;
    }

    while (low < high) {
        uint32_t mid = (low + high) / 2;
        uint32_t mid_leaf = leaf_index(mid);
        if (data[mid_leaf] == key) {
            return mid_leaf;
        }
        else if (data[mid_leaf] > key) {
            high = mid_leaf - logN();
        }
        else {
            if (mid_leaf == low) {
                break; // avoid infinite loop, can only occur if high = low + 1 leaf
            }
            low = mid_leaf;
        }
    }

    // if low != high, then key could be either in low or high. Check first element of high to find out which
    if (data[high] == key) {
        return high;
    }
    if (data[high] < key) {
        low = high;
    }
    // search leaf
    uint32_t leaf = leaf_index(low);
    for (uint32_t i = 0; i < logN(); i += 1) {
        if (data[leaf + i] == INT_NULL) {
            return leaf + i - 1;
        }
        if (data[leaf + i] == key) {
            return leaf + i;
        }
        if (i != (logN() - 1) && data[leaf + i] < key && data[leaf + i + 1] > key) {
            return leaf + i;
        }
    }

    cerr << "Leaf of PMA too dense (is the PMA properly left packed?)" << endl;
    assert(false); // in theory this should never happen, since that would mean leaf is completely full
    return -1;
}

void SetPMA::slide_right(uint32_t index) {
    uint64_t right;
    uint64_t left = data[index];
    for (uint32_t i = index; i < leaf_position(leaf_number(index) + 1); i++) {
        right = data[i + 1];
        data[i + 1] = left;
        left = right;
        if (i == _max_index) {
            _max_index = i + 1;
        }
        if (left == INT_NULL) {
            break;
        }
    }
    // for (uint32_t i = leaf_position(leaf_number(index) + 1) - 1; i > index; --i) {
    //     if (i - 1 == _max_index) {
    //         _max_index = i;
    //     }
    //     data[i] = data[i - 1];
    // }
}

bool SetPMA::insert(uint64_t key) {
    // print_pma();
    // auto start = std::chrono::high_resolution_clock::now();
    bool resized = false;
    uint32_t index = search(key);
    // auto search_end = std::chrono::high_resolution_clock::now();
    if (index != INVALID_IDX && data[index] == key) {
        return resized;
    }
    _size += 1;
    if (index == _max_index) {
        _max_index = index + 1;
    }
    index = index + 1;

    // always deposit on the left
    if (data[index] != INT_NULL) {
        slide_right(index);
    }

    data[index] = key;

    // get density of the leaf you are in

    uint32_t len = logN();
    uint32_t node_index = leaf_index(index);
    uint32_t density_count = count_nonempty(node_index, len);
    // while density too high, go up the implicit tree
    // go up to the biggest node above the density bound
    int level = 0;
    while (density_count > (uint32_t) ((_leaf_max - 0.01 * level) * len) && (len < _capacity)) {
        len *= 2;
        uint32_t new_node_index = (node_index / len) * len;

        if (new_node_index < node_index) {
            density_count += count_nonempty(new_node_index, len / 2);
        } else {
            density_count += count_nonempty(new_node_index + len / 2, len / 2);
        }
        node_index = new_node_index;
        level += 1;
    }

    if (len == _capacity && density_count > (uint32_t) (_leaf_max * len)) {
        // need to double PMA, will disallow this
        if (_resize_allowed) {
            resized = true;
            resize();
        }
    }
    else if (len > logN()) {
        redistribute(node_index, len, density_count);
    }
    
    return resized;
}

void SetPMA::gather(uint32_t index, uint32_t len, uint32_t density_count) {
    _temp.reserve(density_count); // t - s, t = index
    _temp.clear();
    for (uint32_t i = index; i < len + index;) {
        if (data[i] != INT_NULL) {
            _temp.push_back(data[i]);
            data[i] = INT_NULL;
            i += 1;
        }
        else {
            i = next_leaf(i); // doesn't work if capacity just doubled. 
        }
    }

    // Fill with nulls, because all of our data is in temp
    if (index == 0 && len == capacity()) {
        cout << "Temp size: " << _temp.size() << endl;
        std::fill(data.begin(), data.end(), INT_NULL);
    }
}

void SetPMA::distribute(uint32_t index, uint32_t len, uint32_t density_count) {
    uint32_t nl = len / logN();
    uint32_t elems_per_leaf = density_count / nl;
    uint32_t x = 0;
    for (uint32_t leaf = 0; leaf < nl; ++leaf) {
        uint32_t num_elems_to_copy = elems_per_leaf + (leaf < density_count % nl);
        memcpy(&data[index + leaf * logN()], &_temp[x], num_elems_to_copy * sizeof(uint64_t));
        x += num_elems_to_copy;
        // if we are on final iteration and _max_index is within the redistribute, update _max_index
        // _max_index can be INVALID_IDX if there have never been any inserts (e.g. we are redistributed to as first operation of this PMA)
        if (leaf == nl - 1 && (_max_index >= index && _max_index < index + len || _max_index == INVALID_IDX)) {
            if (num_elems_to_copy == 0) {
                cout << "???" << endl;
                exit(-1);
            }
            _max_index = index + leaf * logN() + (num_elems_to_copy - 1);
        }
    }
}

void SetPMA::redistribute(uint32_t index, uint32_t len, uint32_t density_count) {
    gather(index, len, density_count);
    distribute(index, len, density_count);
}


// void SetPMA::debug_resize() {
//     unordered_set<uint64_t> before;
//     unordered_set<uint64_t> after;
//     cout << "doing resize" << endl;
//     for (uint64_t element : data) {
//         if (element != INT_NULL) {
//             before.insert(element);
//         }
//     }
    
//     _capacity *= 2;
//     data.resize(_capacity, INT_NULL);
//     redistribute(0, _capacity, _size);
//     for (uint64_t element : data) {
//         if (element != INT_NULL) {
//             after.insert(element);
//         }
//     }
    
//     // ensure no duplicates
//     cout << "Before size: " << before.size() << endl;
//     cout << "After size: " << after.size() << endl;
//     assert(before.size() == after.size());
//     assert(before == after);
// }


void SetPMA::resize() {
    gather(0, _capacity, _size);
    _capacity *= 2;
    data.resize(_capacity, INT_NULL);
    distribute(0, _capacity, _size);
}

void SetPMA::print_pma(ostream& stream) {
    stream << "[";
    for (auto element : data) {
        if (element == INT_NULL) {
            stream << "null" << " ";
        }
        else {
            stream << element << " ";
        }
    }
    stream << "]" << endl;
}

uint64_t SetPMA::range_sum(uint64_t left, uint64_t right) {
    uint64_t sum = 0;
    range_func summer([&sum](uint64_t v){ sum += v; });
    range(left, right, summer);
    return sum;
}

// left inclusive, right exclusive. calls op with value of every element within range
// potentially return whether this found any values in range?
void SetPMA::range(uint64_t left, uint64_t right, range_func& op) {
    // print_pma();
    uint32_t left_index = search(left);
    uint32_t right_index = search(right);
    
    if (right_index == INVALID_IDX) {
        return;
    }

    if (left_index == INVALID_IDX) {
        left_index = 0;
    }


    if (data[left_index] < left || data[left_index] == INT_NULL) {
        if (left_index == right_index) {
            return;
        } else {
            left_index += 1;
        }
    }

    right_index += 1;
    for (uint32_t i = left_index; i < right_index; ) {
        if (data[i] != INT_NULL) {
            op(data[i]);
            i += 1;
        }
        else {
            i = next_leaf(i);
        }
    }
}

