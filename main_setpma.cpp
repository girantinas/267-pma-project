#include <chrono>
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

#include "set_pma.hpp"
// #include "butil.hpp"

using namespace std;

int main(int argc, char** argv) {
    ifstream infile;
    ofstream outfile;
    if (argc < 2) {
        infile.open("benchmark_inserts.txt");
        if (!infile.is_open()) {
            cerr << "Couldn't open benchmark_inserts.txt" << endl;
            return -1;
        }
    } else {
        infile.open(argv[1]);
        if (!infile.is_open()) {
            cerr << "Couldn't open " << argv[1] << endl;
            return -1;
        }
    }

    if (argc >= 3) {
        outfile.open(argv[2]); // Use "benchmark_output.txt" or something like that 
        if (!outfile.is_open()) {
            cerr << "Couldn't open " << argv[2] << endl;
            return -1;
        }
    }

    bool write_output = outfile.is_open();
    
    string line;

    SetPMA pma(1 << 22);
    int num_range_queries = 1;
    int line_num = 0;
    while (getline(infile, line)) {
        istringstream iss(line);
        string command;
        iss >> command;
        if (command == "PUT") {
            uint64_t insert_num;
            iss >> insert_num;
            // cout << "put(" << insert_num << ")" << endl;
            pma.insert(insert_num);
        } else if (command == "RANGE_QUERY") {
            uint64_t range_start, range_end;
            iss >> range_start >> range_end;
            // cout << "range query(" << range_start << "," << range_end << ")" << ", number " << num_range_queries << endl;
            uint64_t result = pma.range_sum(range_start, range_end);
            if (write_output) {
                outfile << result << endl;
            }
            num_range_queries++;
        } else if (command == "QUERY_ALL") {
            cout << "query all" << endl;
            uint64_t result = pma.range_sum(0, INT_MAX - 1);
            if (write_output) {
                outfile << result << endl;
            }
        } else {
            cerr << "Received unsupported command" << endl;
        }
        
        if (line_num % 100000 == 0 || line_num > 2700000 && line_num % 10000 == 0) {
            cout << "line " << line_num << endl;
        }
        line_num += 1;
    }
    cout << "done!" << endl;
    infile.close();
    outfile.close();

    return 0;
}
