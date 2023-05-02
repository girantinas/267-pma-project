#include <cstdint>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <omp.h>

using namespace std;

inline uint64_t make_edge_tuple(uint32_t from, uint32_t to) { return (((uint64_t) from) << 32) | to; }
inline pair<uint32_t, uint32_t> get_edge_tuple(uint64_t edge) { return make_pair(static_cast<uint32_t>(edge >> 32), static_cast<uint32_t>(edge)); }

inline uint32_t hashInt(uint32_t a) {
   a = (a+0x7ed55d16) + (a<<12);
   a = (a^0xc761c23c) ^ (a>>19);
   a = (a+0x165667b1) + (a<<5);
   a = (a+0xd3a2646c) ^ (a<<9);
   a = (a+0xfd7046c5) + (a<<3);
   a = (a^0xb55a4f09) ^ (a>>16);
   return a;
}

inline uint64_t hashInt(uint64_t a) {
   a = (a+0x7ed55d166bef7a1d) + (a<<12);
   a = (a^0xc761c23c510fa2dd) ^ (a>>9);
   a = (a+0x165667b183a9c0e1) + (a<<59);
   a = (a+0xd3a2646cab3487e3) ^ (a<<49);
   a = (a+0xfd7046c5ef9ab54c) + (a<<3);
   a = (a^0xb55a4f090dd4a67b) ^ (a>>32);
   return a;
}

double hashDouble(uint32_t i) {
  return ((double) (hashInt((uint32_t)i))/((double) UINT32_MAX));
}

struct rMat {
  double a, ab, abc;
  uint32_t n; 
  uint32_t h;
  ofstream outfile;
  rMat(uint32_t _n, uint32_t _seed, 
       double _a, double _b, double _c) {
    n = _n; a = _a; ab = _a + _b; abc = _a+_b+_c;
    h = hashInt(_seed);
  }

  pair<uint32_t, uint32_t> rMatRec(uint32_t nn, uint32_t randStart, uint32_t randStride) {
    if (nn==1) return make_pair(0, 0);
    else {
      pair<uint32_t, uint32_t> x = rMatRec(nn/2, randStart + randStride, randStride);
      double r = hashDouble(randStart); 
      if (r < a) { return x; }
      if (r < ab) {
        return make_pair(x.first,x.second+nn/2);
      }
      if (r < abc) { 
        return make_pair(x.first+nn/2, x.second); 
      }
      return make_pair(x.first+nn/2, x.second+nn/2);
    }
  }

  uint64_t operator() (uint64_t i) {
    uint32_t randStart = hashInt((2*i)*h);
    uint32_t randStride = hashInt((2*i+1)*h);
    auto x = rMatRec(n, randStart, randStride);
    return make_edge_tuple(x.first, x.second);
  }
};

void edgeRmat(uint32_t logn, uint64_t m, uint32_t seed, 
		   float a, float b, float c) {
  uint32_t nn = (1 << logn);
  rMat g(nn, seed, a, b, c);
  ofstream outfile;
  outfile.open("tests/rmat-inserts.txt");
  if (!outfile.is_open()) {
      cerr << "Couldn't open outfile" << endl;
      exit(-1);
  }  
  
  #pragma omp parallel
  {
    uint32_t bfs_vertex;
    int thread_num = omp_get_thread_num();
    int total_threads = omp_get_num_threads();

    for (uint64_t i = 0; i < m / total_threads; i++) {
      uint64_t edge = g(i + (m / total_threads) * thread_num);
      outfile << edge << '\n';

      if (i % 10000 == 6666) { // lol
        bfs_vertex = get_edge_tuple(edge).first;
      } else if (i % 10000 == 10000 - 1) {
        outfile << "BFS " << bfs_vertex << endl;
      }
    }
  }
}

int main(int argc, char* argv[]) {
  uint32_t logn = 24;
  double a = 0.5;
  double b = 0.1;
  double c = 0.1;
  uint64_t m = 32 * (uint64_t) (1 << logn); // 500 million edges
  uint32_t seed = 5;
  edgeRmat(logn, m, seed, a, b, c);

  return 0;
}
