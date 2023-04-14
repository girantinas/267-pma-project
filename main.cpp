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

    SetPMA pma(1 << 4);
    int num_range_queries = 1;

    while (getline(infile, line)) {
        istringstream iss(line);
        string command;
        iss >> command;
        if (command == "PUT") {
            int insert_num;
            iss >> insert_num;
            cout << "put(" << insert_num << ")" << endl;
            pma.insert(insert_num);
        } else if (command == "RANGE_QUERY") {
            int range_start, range_end;
            iss >> range_start >> range_end;
            cout << "range query(" << range_start << "," << range_end << ")" << ", number " << num_range_queries << endl;
            int result = pma.range(range_start, range_end);
            if (write_output) {
                outfile << result << endl;
            }
            num_range_queries++;
        } else if (command == "QUERY_ALL") {
            cout << "query all" << endl;
            int result = pma.range(0, INT_MAX - 1);
            if (write_output) {
                outfile << result << endl;
            }
        } else {
            cerr << "Received unsupported command" << endl;
        }
    }
    infile.close();
    outfile.close();

    return 0;
}
