#include <cassert>
#include <string>
#include <vector>
#include <cmath>
#include <cctype>
#include <queue>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <numeric>

using namespace std;

const int NUM_MATCH_THRESHES = 6;
const int MATCH_THRESHES[NUM_MATCH_THRESHES] = {20, 40, 60, 80, 100, 120};

/*
 0 matches int unsigned ,           # Number of bases that match that aren't repeats
 1 misMatches int unsigned ,        # Number of bases that don't match
 2 repMatches int unsigned ,        # Number of bases that match but are part of repeats
 3 nCount int unsigned ,            # Number of 'N' bases
 4 qNumInsert int unsigned ,        # Number of inserts in query
 5 qBaseInsert int unsigned ,       # Number of bases inserted in query
 6 tNumInsert int unsigned ,        # Number of inserts in target
 7 tBaseInsert int unsigned ,       # Number of bases inserted in target
 8 strand char(2) ,                 # + or - for query strand, opt. followed by + or . for tar
 9 qName varchar(255) ,             # Query sequence name
10 qSize int unsigned ,             # Query sequence size
11 qStart int unsigned ,            # Alignment start position in query
12 qEnd int unsigned ,              # Alignment end position in query
13 tName varchar(255) ,             # Target sequence name
14 tSize int unsigned ,             # Target sequence size
15 tStart int unsigned ,            # Alignment start position in target
16 tEnd int unsigned ,              # Alignment end position in target
17 blockCount int unsigned ,        # Number of blocks in alignment. A block contains no gaps.
18 blockSizes longblob ,            # Size of each block in a comma separated list
19 qStarts longblob ,               # Start of each block in query in a comma separated list
20 tStarts longblob ,               # Start of each block in target in a comma separated list
*/

struct Hit {
  string fullLine;
  string qName, tName;
  int matches, misMatches, qId, qStart, qEnd, tStart, tEnd;
};

map < int, vector <Hit> > read_psl(const char *filename) {
  map < int, vector <Hit> > hits;

  ifstream fin(filename);
  string line;
  while (getline(fin, line)) {
    vector <int> field_starts(1, 0);
    for (int i = 0; i < (int) line.length(); i++)
      if (line[i] == '\t')
	field_starts.push_back(i+1);
    Hit hit;
    hit.fullLine = line;
    sscanf(line.c_str()+field_starts[0], "%d", &hit.matches);
    sscanf(line.c_str()+field_starts[1], "%d", &hit.misMatches);
    hit.qName = line.substr(field_starts[9], field_starts[10]-1-field_starts[9]);
    int qIdPos = field_starts[9];
    while (!isdigit(line[qIdPos])) qIdPos++;
    sscanf(line.c_str()+qIdPos, "%d", &hit.qId);
    hit.tName = line.substr(field_starts[13], field_starts[14]-1-field_starts[13]);
    sscanf(line.c_str()+field_starts[11], "%d", &hit.qStart);
    sscanf(line.c_str()+field_starts[12], "%d", &hit.qEnd);
    sscanf(line.c_str()+field_starts[15], "%d", &hit.tStart);
    sscanf(line.c_str()+field_starts[16], "%d", &hit.tEnd);
    if (hit.tEnd - hit.tStart <= 1000) // ignore hits to disparate chunks
      hits[hit.qId].push_back(hit);
  }
  return hits;
}

bool good_overlap(int start1, int end1, int start2, int end2) {
  // overlap at least 50% of each range
  return 2 * (min(end1, end2) - max(start1, start2)) >= max(end1-start1, end2-start2);
    
}

bool is_present(const Hit &hit, const vector <Hit> &hits_vec) {
  for (int i = 0; i < (int) hits_vec.size(); i++) {
    if (/*hit.matches == hits_vec[i].matches &&
	  hit.misMatches == hits_vec[i].misMatches &&*/
	good_overlap(hit.tStart, hit.tEnd, hits_vec[i].tStart, hits_vec[i].tEnd) &&
	hit.tName.substr(0, 20) == hits_vec[i].tName.substr(0, 20))
      return true;
  }
  return false;
}

void write_array(const int *arr, int len) {
  for (int i = 0; i < len; i++)
    fprintf(stderr, "%7d", arr[i]);
  fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {

  map < int, vector <Hit> > blat_hits = read_psl(argv[1]);
  map < int, vector <Hit> > cablat_hits = read_psl(argv[2]);

  set <int> qIds;
  for (__typeof(blat_hits.begin()) it = blat_hits.begin(); it != blat_hits.end(); it++)
    qIds.insert(it->first);
  for (__typeof(cablat_hits.begin()) it = cablat_hits.begin(); it != cablat_hits.end(); it++)
    qIds.insert(it->first);

  vector <int> num_found(NUM_MATCH_THRESHES), num_missed(NUM_MATCH_THRESHES),
    num_extra(NUM_MATCH_THRESHES);
  
  for (set <int>::iterator it = qIds.begin(); it != qIds.end(); it++) {
    int qId = *it;
    printf("---- query %d ----\n", qId);
    for (int i = 0; i < (int) blat_hits[qId].size(); i++) {
      const Hit &hit = blat_hits[qId][i];

      if (is_present(hit, cablat_hits[qId])) {
	printf("found: match %d, misMatch %d, query %d-%d, target %d-%d %s\n",
	       hit.matches, hit.misMatches, hit.qStart, hit.qEnd, hit.tStart, hit.tEnd,
	       hit.tName.substr(0, 20).c_str());
	for (int i = 0; i < NUM_MATCH_THRESHES; i++)
	  if (hit.matches >= MATCH_THRESHES[i])
	    num_found[i]++;
      }
      else {
	printf("*** missed ***: match %d, misMatch %d, query %d-%d, target %d-%d %s\n",
	       hit.matches, hit.misMatches, hit.qStart, hit.qEnd, hit.tStart, hit.tEnd,
	       hit.tName.substr(0, 20).c_str());
	cout << hit.fullLine << endl;
	for (int i = 0; i < NUM_MATCH_THRESHES; i++)
	  if (hit.matches >= MATCH_THRESHES[i])
	    num_missed[i]++;
      }
    }
    for (int i = 0; i < (int) cablat_hits[qId].size(); i++) {
      const Hit &hit = cablat_hits[qId][i];

      // check for double-counts
      bool already_seen = false;
      for (int j = 0; j < i; j++)
	if (cablat_hits[qId][j].tStart <= hit.tEnd &&
	    cablat_hits[qId][j].tEnd >= hit.tStart) {
	  cout << "ignoring double count" << endl;
	  already_seen = true;
	  break;
	}
      if (already_seen) continue;

      if (is_present(hit, blat_hits[qId])) {
	;
      }
      else {
	printf("extra: match %d, misMatch %d, query %d-%d, target %d-%d %s\n",
	       hit.matches, hit.misMatches, hit.qStart, hit.qEnd, hit.tStart, hit.tEnd,
	       hit.tName.substr(0, 20).c_str());
	cout << hit.fullLine << endl;
	for (int i = 0; i < NUM_MATCH_THRESHES; i++)
	  if (hit.matches >= MATCH_THRESHES[i])
	    num_extra[i]++;
      }
    }
  }
  
  fprintf(stderr, "\n");
  fprintf(stderr, "thresh: "); write_array(&MATCH_THRESHES[0], NUM_MATCH_THRESHES);
  fprintf(stderr, "found:  "); write_array(&num_found[0], NUM_MATCH_THRESHES);
  fprintf(stderr, "missed: "); write_array(&num_missed[0], NUM_MATCH_THRESHES);
  fprintf(stderr, "extra:  "); write_array(&num_extra[0], NUM_MATCH_THRESHES);

  return 0;
}
