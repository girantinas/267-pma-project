#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>
#include <fstream> 
#include <sstream>
#include <iostream>

#include "dist_pcsr.hpp"
// #include "butil.hpp"

using namespace std;

int main(int argc, char** argv) {
    upcxx::init();
    ifstream infile;
    
    if (argc < 2) {
        infile.open("dist_pcsr_inserts.txt");
        if (!infile.is_open()) {
            cerr << "Couldn't open pcsr_inserts.txt" << endl;
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
        outfile.open(argv[2] + "_" + std::to_string(upcxx::rank_me()) + ".dat");
        if (!outfile.is_open()) {
            cerr << "Couldn't open " << argv[2] << endl;
            return -1; 
        }
    }

    bool write_output = outfile.is_open();

    string line;

    upcxx::dist_object<DistPCSR> pcsr(1 << 4, 16384);

    upcxx::barrier();
    line_number = 0;
    while (getline(infile, line)) {
        istringstream iss(line);
        string command;
        iss >> command;
        if (command == "START_INSERTS" || command == "START_QUERIES") {
            upcxx::barrier();
            line_number++;
            continue;
        }
        if (line_number % upcxx::rank_n() != upcxx::rank_me()) {
            line_number++;
            continue;
        }
        if (command == "PUT_EDGE") {
            uint32_t u, v;
            iss >> u >> v;
            insert_edge(pcsr, u, v);
        } else if (command == "QUERY_EDGE") {
            uint32_t u, v;
            iss >> u >> v;
            bool exists = query_edge(pcsr, u, v); // TODO: batch these
            if (exists) {
                outfile << "True" << endl;
            } else {
                outfile << "False" << endl;
            }
            outfile << endl;
        } else if (command == "GET_OUT_EDGES") {
            uint32_t vertex;
            iss >> vertex;
            vector<uint32_t> out_edges = edges(pcsr, vertex);
            if (write_output) {
                if (out_edges.empty()) {
                    outfile << "V: " << vertex << ", E: ";
                } else {
                    outfile << "V: " << vertex << ", E:";
                    for (uint32_t dest_vertex : out_edges) {
                        outfile << " " << dest_vertex;
                    }
                }
                outfile << endl;
            }
            outfile << endl;
        } else if (command == "QUERY_ALL_GRAPH") {
            vector<pair<uint32_t, vector<uint32_t>>> result = adjacency_lists(pcsr);
            if (write_output) {
                for (auto& entry : result) {
                    uint32_t& vertex = entry.first;
                    outfile << "V: " << vertex << ", E:";
                    for (int dest_vertex : entry.second) {
                        outfile << " " << dest_vertex;
                    }
                    outfile << endl;
                }
                outfile << endl;
            }
        } else {
            cerr << "Received unsupported command: " << command << endl;
        }
    }

    infile.close();
    outfile.close();
    upcxx::finalize();
    return 0;
}