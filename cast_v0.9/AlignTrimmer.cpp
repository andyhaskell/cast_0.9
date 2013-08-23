#include <utility>
#include <cstring>

using namespace std;

class AlignTrimmer {

public:
  AlignTrimmer(void) {
    memset(to_int_lookup, 0, sizeof(to_int_lookup));
    to_int_lookup['A'] = 0;
    to_int_lookup['C'] = 1;
    to_int_lookup['G'] = 2;
    to_int_lookup['T'] = 3;
    memset(k_mer_locs, -1, sizeof(k_mer_locs));
  }
  
  pair <int, int> trim(char *qdna, int qsize, char *tdna, int tsize) {

    pair <int, int> trim_lengths = make_pair(-100, -100);
    
    // set k_mer_locs table
    int k_mer = make_k_mer(qdna);
    for (int i = 0; /*i < NUM_QUERY_K_MERS*/; i++) {
      k_mers[i] = k_mer;
      k_mer_locs[k_mer] = i;
      if (i == NUM_QUERY_K_MERS-1) break;
      int base = to_int_lookup[qdna[TRIM_K+i]];
      k_mer = ((k_mer<<2)|base)&K_MASK;
    }

    // go through tdna
    k_mer = make_k_mer(tdna);
    int stop_i = min(tsize-TRIM_K, MAX_TARGET_K_MERS);
    for (int i = 0; /*i < min(MAX_TARGET_K_MERS, tsize-TRIM_K)*/; i++) {
      //cerr << "k_mer: " << k_mer << endl;
      if (k_mer_locs[k_mer] != -1) { // found
	trim_lengths.first = i - k_mer_locs[k_mer];
	//cerr << "found " << i << " " << k_mer_locs[k_mer] << " " << k_mer << endl;
	break;
      }
      if (i == stop_i) break;
      int base = to_int_lookup[tdna[TRIM_K+i]];
      k_mer = ((k_mer<<2)|base)&K_MASK;
    }

    // unset k_mer_locs_table
    for (int i = 0; i < NUM_QUERY_K_MERS; i++)
      k_mer_locs[k_mers[i]] = -1;

    return trim_lengths;
  }

private:
  static const int TRIM_K = 5;
  static const int K_MASK = (1<<(2*TRIM_K))-1;
  static const int NUM_QUERY_K_MERS = 50; // hash this many k-mers at each end of the query
  static const int MAX_TARGET_K_MERS = 75; // scan distance from each end of the target
  int to_int_lookup[256];
  int k_mer_locs[1<<(2*TRIM_K)];
  int k_mers[NUM_QUERY_K_MERS];

  int make_k_mer(char *seq) {
    int k_mer = 0;
    for (int i = 0; i < TRIM_K; i++) {
      int base = to_int_lookup[seq[i]];
      k_mer = (k_mer<<2)|base;
    }
    return k_mer;
  }
};
