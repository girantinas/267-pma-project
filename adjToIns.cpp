#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <omp.h>
#include <sstream>
#include <string>

using namespace std;

inline uint64_t make_edge_tuple(uint32_t from, uint32_t to) { return (((uint64_t) from) << 32) | to; }
inline bool isSpace(char c) {
    switch (c)  {
        case '\r': 
        case '\t': 
        case '\n': 
        case 0:
        case ' ' : return true;
        default : return false;
    }
}

int main(int argc, char* argv[]) {
    // Open files
    std::string infileName = "../rmat_ligra.adj";
    std::string outfileFormat = "./rmat-tests/rmat-inserts-";
    if (argc >= 3) {
        infileName = argv[1];
        outfileFormat = argv[2];
    }

    ifstream infile (infileName, ios::in | ios::binary | ios::ate);
    if (!infile.is_open()) {
        cout << "Unable to open file: " << infileName << endl;
        return -1;
    }

    // Read the file
    long end = infile.tellg();
    infile.seekg(0, ios::beg);
    long n = end - infile.tellg();
    char* bytes = (char*) malloc(sizeof(char) * (n + 1));
    bytes[n] = '\0';
    infile.read(bytes, n);
    infile.close();

    // Separate into words
    #pragma omp parallel for
    for (long i = 0; i < n; ++i) {
        if (isSpace(bytes[i])) { bytes[i] = 0; }
    }

    vector<long> wordStarts = {0};
    
    for (long i = 1; i < n; ++i) {
        if (bytes[i] && !bytes[i - 1]) {
            wordStarts.push_back(i);
        }
    }

    char** strings = (char**) malloc(sizeof(char*) * wordStarts.size());
    #pragma omp parallel for
    for (long j = 0; j < wordStarts.size(); ++j) {
        strings[j] = bytes + wordStarts[j];
    }

    //
    long len = wordStarts.size() - 1;
    uint32_t* In = (uint32_t*) malloc(sizeof(uint32_t) * len);
    #pragma omp parallel for
    for (long i = 0; i < len; ++i) {
        In[i] = static_cast<uint32_t>(stoul(strings[i + 1]));
    }
    n = In[0];
    long m = In[1];

    if (len != n + m + 2) {
        cout << "Bad input file" << endl;
        return -1;
    }

    uint32_t* offsets = In + 2;
    uint32_t* edges = In + 2 + n;

    int errorFlag = 0;

    #pragma omp parallel 
    {
        int numThreads = omp_get_num_threads();
        int thisThread = omp_get_thread_num();
        ofstream outfile;

        outfile.open(outfileFormat + std::to_string(thisThread) + ".txt");
        if (!outfile.is_open()) {
            cerr << "Couldn't open " << argv[2] << endl;
            errorFlag = -1; 
        }

        if (errorFlag == 0) {
            uint32_t currInVertex = 0;
            for (long j = thisThread; j < m; j += numThreads) {
                while (offsets[currInVertex + 1] <= j) {
                    currInVertex++;
                }
                outfile << make_edge_tuple(currInVertex, edges[j]) << '\n';
            }
        }
    }
    return errorFlag;
}