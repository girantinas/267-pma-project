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
#include <deque>

using namespace std;

int main(int argc, char** argv) {
    upcxx::init();
    if (upcxx::rank_me() == 0) {
        cout << "alive" << endl;
        unordered_set<uint64_t> edges;
        ifstream infile;
        string line;
        uint32_t num_dupes = 0;
        for (int i = 0; i < 16; i++) {
            infile.open("../rmat-tests/rmat-inserts-" + std::to_string(i) + ".txt");
            
            while (getline(infile, line)) {
                istringstream iss(line);
                string command;
                iss >> command;
                // cout << command << endl;
                if (command == "BFS") {
                } else {
                    uint64_t edge = stoull(command);
                    if (edges.find(edge) != edges.end()) {
                        num_dupes++;
                    } else {
                        edges.insert(edge);
                    }
                }
                
            }
            infile.close();
            cout << "finished a file" << endl;
        }
        cout << "Num dupes in insert files: " << num_dupes << endl;
    }
    upcxx::finalize();
}