// example usage:
// ./gen_compressed_db uniques.fasta links_compressed.dat edit_scripts.dat smushed < test.fasta | tee gencomp_stdout.txt

/*
 * input args:
 * 1 - compressed db filename: uniques (e.g., uniques_all_seqs.fasta)
 * 2 - compressed db filename: links (e.g., links_compressed.dat)
 * 3 - compressed db filename: edit scripts (e.g., edit_scripts.dat)
 * 4 - uniques file format: 'smushed' OR 'separate' sequences
 * 5 - (optional) window identity threshold for compression (default 85)
 * reads sequences in fasta format from stdin
 * writes log to stderr
 */

#include <cstring>
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

#include "HalfByteWriter.cpp"
#include "LinkInfo.cpp"
#include "DNAUtils.cpp"
#include "Timer.cpp"
#include "IOUtils.cpp"

using namespace std;

Timer timer;

//#define PRINT_MATCHES // to get examples for testing compression
//#define CHECK_AXT
//#define DEBUG_OUTPUT
//#define PRINT_OLD_LINKS_TXT

#define GLOBAL_STATS
#ifdef GLOBAL_STATS
int tot_processed, tot_len;
#endif

// input arguments
const int UNIQ_ALL_ARG_NUM = 1;
const int LINKS_ARG_NUM = 2;
const int EDITS_ARG_NUM = 3;
const int SEPARATE_SEQS_ARG_NUM = 4;
const int WINDOW_THRESH_ARG_NUM = 5;



const int BLAST_K = 10;
const int MAX_K_MER_FREQ = 500;
const int MAX_CONSEC_MISMATCH = 3;
const int OVERLAP = 100;
const int MIN_MATCH_LEN = 3*OVERLAP;
const int MAX_CHUNK_SIZE = 10000;

int max_k_mer_freq = 0, tot_bases_read = 0; // max_k_mer_freq adapts

struct loc {
  int pos;
  int id;
  loc(int _pos=0, int _id=0) : pos(_pos), id(_id) {};
};
vector <loc> k_mer_locs[1<<(2*BLAST_K)];

int k_mer_freqs[1<<(2*BLAST_K)];

typedef unsigned short ushort;

struct match {
  int start0, stop0, start1, stop1;
  int id; // for second sequence in repository
  char dir; // +/-
  int start_script;
  match(int _start0=-1, int _stop0=-1, int _start1=-1, int _stop1=-1,
	int _id=-1, char _dir='-', int _start_script=-1) :
    start0(_start0), stop0(_stop0), start1(_start1), stop1(_stop1), id(_id),
    dir(_dir), start_script(_start_script) {};
};


void fill_buf(char *buf, int &buflen, istream &fin,
	      bool &remaining_record_data, bool &getline_success,
	      string &fasta_header) {
  while (buflen < 2*MAX_CHUNK_SIZE) {
    if (!remaining_record_data) return;
    string line;
    getline_success = getline(fin, line);
    if (!getline_success || line[0] == '>') {
      remaining_record_data = false;
      fasta_header = line;
      return;
    }
    strcpy(buf+buflen, line.c_str());
    buflen += line.length();
  }
}

string to_octal_str(int i) {
  char buf[16];
  sprintf(buf, "%o", i);
  return string(buf);
}

string make_edit_script(const vector <string> &axt, char dir) {
  string edit_script = (dir == '+' ? "0" : "1");
  const string &str = axt[0]; // strings are aligned and include gaps
  const string &ref = axt[2];
  bool insert_open = false, subdel_open = false;
  int last_edit = 0;
  for (int i = 0; i < (int) ref.length(); i++) {
    if (str[i] == ref[i]) {
      insert_open = false;
      subdel_open = false;
    }
    else { // mismatch
      if (ref[i] == '-') { // insertion in str relative to ref (i.e., gap in ref)
	subdel_open = false;
	if (!insert_open) { // indicate start of insertion
	  insert_open = true;
	  edit_script += "i" + to_octal_str(i - last_edit);
	  last_edit = i;
	}
	edit_script += str[i];
      }
      else { // substitution or deletion in str (represented in script by '-') relative to ref
	insert_open = false;
	if (!subdel_open) { // indicate start of subdel
	  subdel_open = true;
	  edit_script += "s" + to_octal_str(i - last_edit);
	  last_edit = i;
	}
	edit_script += str[i];
      }
    }
  }
  return edit_script;
}

struct edit_info {
  bool is_subdel; // 0 if insert, 1 if subdel
  int last_dist;
  string str;
};

// reads one edit worth of info (or fails if current decoded char is digit)
// moves pos accordingly
bool next_edit(const string &edit_script, int &pos, edit_info &edit) {
  if (isdigit(edit_script[pos]))
    return false;
  edit.is_subdel = edit_script[pos++] == 's';
  edit.last_dist = 0; edit.str = "";
  while (isdigit(edit_script[pos])) {
    edit.last_dist *= 8; // octal encoding
    edit.last_dist += edit_script[pos++] - '0';
  }
  while (isupper(edit_script[pos]) || edit_script[pos] == '-')
    edit.str += edit_script[pos++];
  return true;
}

string read_edit_script(const string &edit_script, const string &orig) {
  string str;
  int script_pos = 0;
  edit_info edit;
  int orig_pos = 0, last_edit_str_len = 0; // length of last edit str
  while (next_edit(edit_script, script_pos, edit)) {
    // chunk after previous edit
    str += orig.substr(orig_pos, edit.last_dist - last_edit_str_len);

    // update position in original string
    orig_pos += edit.last_dist - last_edit_str_len;

    // append replacement string in edit script; get rid of dashes
    for (int i = 0; i < (int) edit.str.length(); i++)
      if (edit.str[i] != '-')
	str += edit.str[i];

    // skip subdel along original string
    if (edit.is_subdel) orig_pos += edit.str.length();

    last_edit_str_len = edit.str.length();

    /*
    check += (edit.is_subdel ? 's' : 'i');
    check += to_octal_str(edit.last_dist);
    check += edit.str;
    */
  }
  str += orig.substr(orig_pos); // chunk after last edit
  return str;
}

match copy_chunk(vector <string> &seqs, vector <int> &seqs_cum_lengths,
		 char *buf, int buf_copied, int offset) {
  int id = seqs.size();
  seqs.push_back(string(buf, buf_copied));
  seqs_cum_lengths.push_back(seqs_cum_lengths.back() + buf_copied);
#ifdef GLOBAL_STATS
  tot_len += buf_copied;
#endif
  // update k_mer table
  int k_mask = (1<<(2*BLAST_K)) - 1;
  int k_mer = 0;
  int last_bad = 0;
  for (int i = 0; i < buf_copied; i++) {
    int seq_i = to_int(buf[i]);
    k_mer <<= 2;
    if (seq_i == -1) {
      last_bad = 0;
    }
    else {
      last_bad++;
      k_mer |= seq_i;
      tot_bases_read++;
    }
    k_mer &= k_mask;
    if (last_bad >= BLAST_K) {
      if (k_mer_freqs[k_mer] < MAX_K_MER_FREQ) {
	k_mer_locs[k_mer][k_mer_freqs[k_mer]] = loc(i, id);
	k_mer_freqs[k_mer]++;
      }
    }
  }
  max_k_mer_freq = (tot_bases_read>>(2*BLAST_K))+2;

  return match(offset, offset+buf_copied-1, 0, buf_copied-1, id, '+', -1);
}

int mismatch_score(char c1, char c2, char dir_prod) {
  /* may want to rethink this... note that 100 will break the backtrack below
  if (c1 == '.' || c2 == '.') return 100;
  else if (c1 == 'N' || c2 == 'N') return 1;
  else */
  if (c1 == 'N' || c2 == 'N') return 1;
  else return dir_prod == 1 ? c1 != c2 : to_int(c1) + to_int(c2) != 3;
}



const int MIN_PROGRESS = 50; // trigger gapped extension
const int CONSEC_MATCH_CLUMP_SIZE = 4;
const int MAX_DP_LEN = 25;
const int WINDOW_SIZE = 100;
int WINDOW_IDENT_THRESH = 85; // default window identity threshold
const int BTWN_MATCH_MIN_DIST_CHECK = 10;
const double BTWN_MATCH_IDENT_THRESH = 0.5;
char ismatch[2*MAX_CHUNK_SIZE];
int dp_score[MAX_DP_LEN+1][MAX_DP_LEN+1];
int dp_from[MAX_DP_LEN+1][MAX_DP_LEN+1];

int max_dp_len(int i, int dir, int len) {
  return dir == 1 ? min(MAX_DP_LEN, len-i) : min(MAX_DP_LEN, i+1);
}

bool check_and_update(int &last_consec_match_pos, int ismatch_pos,
		      int &num_matches_in_window) {
  int &p = last_consec_match_pos;
  for (; p < ismatch_pos; p++) {
    num_matches_in_window += ismatch[p] - ismatch[p-WINDOW_SIZE];
    if (num_matches_in_window < WINDOW_IDENT_THRESH)
      return false; // check... last_consec_match_pos may be wrong? (includes bad region?)
  }
  return true;
}

void reverse(string &s) {
  int len = s.length();
  for (int i = 0; 2*i < len-1; i++)
    swap(s[i], s[len-1-i]);
}

vector <string> attempt_gapped_ext(int &i1, const int dir1, const char *s1, int len1,
				   int &i2, const int dir2, const string &s2) {
  
  vector <string> axt(3);
  
#ifdef DEBUG_OUTPUT
  cerr << "starting attempt_gapped_ext: i1=" << i1 << (dir1==1?'+':'-')
       << " i2=" << i2 << (dir2==1?'+':'-') << endl;
#endif
  const int dir_prod = dir1*dir2;
  assert(mismatch_score(s1[i1], s2[i2], dir_prod) == 0);
  i1 += dir1; i2 += dir2;
  int len2 = s2.length();
  
  int ismatch_pos = WINDOW_SIZE;
  int num_matches_in_window = WINDOW_SIZE;
  int last_consec_match_pos = ismatch_pos;
  bool found_bad_window = false;

  while (true) {

    int old_ismatch_pos = ismatch_pos;
    int old_i1 = i1, old_i2 = i2;
    
    // ungapped extension until unable to find 4-mer match
    // with 50% identity in between
    int consec_matches = CONSEC_MATCH_CLUMP_SIZE;
    int matches_since_last_consec = 0;
    while (i1 >= 0 && i1 < len1 && i2 >= 0 && i2 < len2) {
      char cur_ismatch = ismatch[ismatch_pos++] =
	mismatch_score(s1[i1], s2[i2], dir_prod) == 0;
      i1 += dir1; i2 += dir2;
      if (cur_ismatch) { // current bases match
	consec_matches++;
	if (consec_matches >= CONSEC_MATCH_CLUMP_SIZE) { // if 4-mer match, check window condition
	  if (!check_and_update(last_consec_match_pos, ismatch_pos,
				num_matches_in_window)) {
	    int pos_past_bad_window = ismatch_pos - last_consec_match_pos;
	    i1 -= dir1 * pos_past_bad_window;
	    i2 -= dir2 * pos_past_bad_window;
	    ismatch_pos = last_consec_match_pos;
	    found_bad_window = true; // window test failed; full stop (will exit outer loop)
#ifdef DEBUG_OUTPUT
	    cerr << "found bad window: i1 = " << i1 << " i2 = " << i2 << endl;
#endif
	    break;
	  }
	  matches_since_last_consec = 0;
	}
	else
	  matches_since_last_consec++;
      }
      else { // mismatch
	consec_matches = 0;
	int pos_since_last_consec = ismatch_pos - last_consec_match_pos;
	if (pos_since_last_consec >= BTWN_MATCH_MIN_DIST_CHECK) { // if 10 bases since last 4-mer
	  if (matches_since_last_consec <
	      pos_since_last_consec * BTWN_MATCH_IDENT_THRESH) { // check for 50% identity
	    i1 -= dir1 * pos_since_last_consec;
	    i2 -= dir2 * pos_since_last_consec;
	    ismatch_pos = last_consec_match_pos;
	    break; // trigger DP
	  }
	}
      }
    }
#ifdef DEBUG_OUTPUT
    assert(ismatch_pos - old_ismatch_pos == (i1 - old_i1)*dir1);
    assert(ismatch_pos - old_ismatch_pos == (i2 - old_i2)*dir2);
#endif
    // create alignment strings corresponding to non-gap region
    string subs1, smatch, subs2;
    for (int x = 0; x < ismatch_pos - old_ismatch_pos; x++) {
      subs1 += dir_prod == 1 ? s1[old_i1+x*dir1] :
	base_pair_comp(s1[old_i1+x*dir1]);
      subs2 += s2[old_i2+x*dir2];
      assert((subs1[x] == subs2[x] && subs1[x] != 'N') == ismatch[old_ismatch_pos+x]);
      smatch += ismatch[old_ismatch_pos+x] ? '|' : ' ';
    }
    axt[0] += subs1; // will reverse them all at the end if dir2 == -1
    axt[1] += smatch;
    axt[2] += subs2;

#ifdef DEBUG_OUTPUT
    for (int x = 0; x < ismatch_pos - old_ismatch_pos; x += 80) {
      cerr << subs1.substr(x, 80) << endl;
      cerr << smatch.substr(x, 80) << endl;
      cerr << subs2.substr(x, 80) << endl;
    }
#endif

    if (found_bad_window) break; // stop align

    // DP on small square to find gapped extension to next 4-mer match
    int dp_len1 = max_dp_len(i1, dir1, len1);
    int dp_len2 = max_dp_len(i2, dir2, len2);
#ifdef DEBUG_OUTPUT
    cerr << "dp_len1: " << dp_len1 << " ";
    for (int x = 0; x < dp_len1; x++)
      cerr << (dir1 == 1 ? s1[i1+x*dir1] : base_pair_comp(s1[i1+x*dir1]));
    cerr << endl;
    cerr << "dp_len2: " << dp_len2 << " ";
    for (int x = 0; x < dp_len2; x++)
      cerr << (dir2 == 1 ? s2[i2+x*dir2] : base_pair_comp(s2[i2+x*dir2]));
    cerr << endl;
#endif
    for (int j2 = 0; j2 <= dp_len2; j2++) {
      dp_score[0][j2] = -3*j2;
      dp_from[0][j2] = 2;
    }
    for (int j1 = 1; j1 <= dp_len1; j1++) {
      dp_score[j1][0] = -3*j1;
      dp_from[j1][0] = 1;
      for (int j2 = 1; j2 <= dp_len2; j2++) {
	int score0 = dp_score[j1-1][j2-1] +
	  (mismatch_score(s1[i1+dir1*(j1-1)], s2[i2+dir2*(j2-1)], dir_prod) ?
	   -3 : 1);
	int score1 = dp_score[j1-1][j2] - 3; // advance 1; gap in 2
	int score2 = dp_score[j1][j2-1] - 3; // advance 2; gap in 1
	if (score0 >= score1 && score0 >= score2) {
	  dp_score[j1][j2] = score0;
	  dp_from[j1][j2] = 0;
	}
	else if (score1 >= score2) {
	  dp_score[j1][j2] = score1;
	  dp_from[j1][j2] = 1;
	}
	else {
	  dp_score[j1][j2] = score2;
	  dp_from[j1][j2] = 2;
	}
      }
    }

    // find edge cell on bottom-right border with best score
    int max_j1 = 0, max_j2 = 0, max_dp_score = -1000;
    for (int j1 = 0; j1 <= dp_len1; j1++)
      if (dp_score[j1][dp_len2] >= max_dp_score) {
	max_dp_score = dp_score[j1][dp_len2];
	max_j1 = j1; max_j2 = dp_len2;
      }
    for (int j2 = 0; j2 <= dp_len2; j2++)
      if (dp_score[dp_len1][j2] >= max_dp_score) {
	max_dp_score = dp_score[dp_len1][j2];
	max_j1 = dp_len1; max_j2 = j2;
      }
    if (max_dp_score == -1000) break;

    // backtrack to last 4-mer match on path
    int cur_j1 = max_j1, cur_j2 = max_j2;
    consec_matches = 0;
    while (!(cur_j1 == 0 && cur_j2 == 0)) {
      if (consec_matches == CONSEC_MATCH_CLUMP_SIZE) { // found chunk; stop
	cur_j1 += CONSEC_MATCH_CLUMP_SIZE;
	cur_j2 += CONSEC_MATCH_CLUMP_SIZE;
	break;
      }
      int prev_j1, prev_j2;
      switch (dp_from[cur_j1][cur_j2]) { // backtrack to previous cell
      case 0: prev_j1 = cur_j1-1; prev_j2 = cur_j2-1; break;
      case 1: prev_j1 = cur_j1-1; prev_j2 = cur_j2; break;
      default: prev_j1 = cur_j1; prev_j2 = cur_j2-1;
      }
      if (dp_from[cur_j1][cur_j2] == 0) { // match or substitution
	if (dp_score[cur_j1][cur_j2] > dp_score[prev_j1][prev_j2]) // match
	  consec_matches++;
	else
	  consec_matches = 0;
	cur_j1--; cur_j2--;
      }
      else { // gap
	consec_matches = 0;
	if (dp_from[cur_j1][cur_j2] == 1)
	  cur_j1--;
	else
	  cur_j2--;
      }
    }
    if (consec_matches < CONSEC_MATCH_CLUMP_SIZE) // failed to find 4-mer; quit
      break;

#ifdef DEBUG_OUTPUT
    cerr << "max_j1: " << max_j1 << " max_j2: " << max_j2 << endl;
    cerr << "cur_j1: " << cur_j1 << " cur_j2: " << cur_j2 << endl;
#endif

    string subs1_dp, smatch_dp, subs2_dp;

    max_j1 = cur_j1; max_j2 = cur_j2;
    // append path and update matches
    int num_steps = 0;
    while (!(cur_j1 == 0 && cur_j2 == 0)) {
      int prev_j1, prev_j2;
      switch (dp_from[cur_j1][cur_j2]) {
	char c1, c2;
      case 0: prev_j1 = cur_j1-1; prev_j2 = cur_j2-1; // match or substitution
	c1 = s1[i1+dir1*prev_j1]; // comp if antisense
	if (dir_prod == -1) c1 = base_pair_comp(c1);
	c2 = s2[i2+dir2*prev_j2];
	subs1_dp += c1;
	smatch_dp += c1 == c2 ? '|' : ' ';
	subs2_dp += c2;
	break;
      case 1: prev_j1 = cur_j1-1; prev_j2 = cur_j2; // advance 1; gap in 2
	c1 = s1[i1+dir1*prev_j1];
	if (dir_prod == -1) c1 = base_pair_comp(c1); // comp if antisense
	subs1_dp += c1;
	smatch_dp += ' ';
	subs2_dp += '-';
	break;
      default: prev_j1 = cur_j1; prev_j2 = cur_j2-1; // advance 2; gap in 1
	c2 = s2[i2+dir2*prev_j2];
	subs1_dp += '-';
	smatch_dp += ' ';
	subs2_dp += c2;
      }
      ismatch[ismatch_pos++] = // note: need to flip order
	dp_score[cur_j1][cur_j2] > dp_score[prev_j1][prev_j2];
      num_steps++;
      cur_j1 = prev_j1; cur_j2 = prev_j2;
    }

    for (int i = 0; i < num_steps/2; i++) // flip order
      swap(ismatch[ismatch_pos-1-i], ismatch[last_consec_match_pos+i]);
    if (!check_and_update(last_consec_match_pos, ismatch_pos,
			  num_matches_in_window)) // failed window threshold; quit
      break;
    else {
      i1 += dir1*max_j1; i2 += dir2*max_j2; // move positions to end of match
      // update axt alignment with dp region
      reverse(subs1_dp);
      reverse(smatch_dp);
      reverse(subs2_dp);
      axt[0] += subs1_dp; // will reverse them all at the end if dir2 == -1
      axt[1] += smatch_dp;
      axt[2] += subs2_dp;

#ifdef DEBUG_OUTPUT
      cerr << subs1_dp << endl;
      cerr << smatch_dp << endl;
      cerr << subs2_dp << endl;
#endif
    }
  }

  i1 -= dir1; i2 -= dir2; // interval is inclusive

  if (dir2 == -1) { // choose direction so reference (s2) is always going forward
    reverse(axt[0]);
    reverse(axt[1]);
    reverse(axt[2]);
  }

  return axt;
}

// returns extension distance (but does not move pointers)
int attempt_ext(int i1, const int dir1, const char *s1, int len1,
		int i2, const int dir2, const string &s2) {
  const int dir_prod = dir1*dir2;
  i1 += dir1; i2 += dir2;
  int len2 = s2.length();
  int progress = 0, consec_mismatch = 0;
  while (consec_mismatch < MAX_CONSEC_MISMATCH &&
	 i1 >= 0 && i1 < len1 && i2 >= 0 && i2 < len2) {
    if (mismatch_score(s1[i1], s2[i2], dir_prod))
      consec_mismatch++;
    else
      consec_mismatch = 0;
    i1 += dir1; i2 += dir2;
    progress++;
  }
  return progress;
}

match find_match(const vector <string> &seqs, const char *buf,
		 const int buflen, HalfByteWriter &script_writer) {
  int k_mer = 0, k_mer_revcomp = 0, k_mask = (1<<(2*BLAST_K))-1;
  int last_bad = 0;
  for (int i = 0; i < min(buflen, MAX_CHUNK_SIZE); i++) {
    int buf_i = to_int(buf[i]);
    k_mer <<= 2; k_mer_revcomp >>= 2;
    if (buf_i == -1) {
      last_bad = 0;
    }
    else {
      last_bad++;
      k_mer |= buf_i;
      k_mer_revcomp |= (3-buf_i)<<(2*BLAST_K-2);
    }
    k_mer &= k_mask;
    
    if (i % BLAST_K != 0) continue;

    if (last_bad >= BLAST_K) {
      if (k_mer_freqs[k_mer] < max_k_mer_freq) { // same dir
	for (int j = 0; j < k_mer_freqs[k_mer]; j++) {
	  int stop0 = i, start0 = stop0 - BLAST_K + 1;
	  int stop1 = k_mer_locs[k_mer][j].pos, start1 = stop1 - BLAST_K + 1,
	    id = k_mer_locs[k_mer][j].id;

#ifdef DEBUG_OUTPUT
	  cerr << "seed k-mer: ";
	  for (int x = start1; x <= stop1; x++)
	    cerr << seqs[id][x];
	  cerr << endl;
#endif
	  vector <string> axt_seed(3);
	  axt_seed[0] = seqs[id].substr(start1, BLAST_K);
	  axt_seed[1] = string(BLAST_K, '|');
	  axt_seed[2] = axt_seed[0];

	  if (attempt_ext(start0, -1, buf, buflen, start1, -1, seqs[id]) +
	      attempt_ext(stop0, 1, buf, buflen, stop1, 1, seqs[id])
	      > MIN_PROGRESS) {
	    vector <string> axt_prev = attempt_gapped_ext(start0, -1, buf, buflen, start1, -1, seqs[id]);
	    vector <string> axt_next = attempt_gapped_ext(stop0, 1, buf, buflen, stop1, 1, seqs[id]);
#ifdef DEBUG_OUTPUT
	    cerr << "aligned: " << stop0 - start0 << endl;
#endif
	    if (stop0 - start0 >= MIN_MATCH_LEN) {

	      vector <string> axt(3);
	      for (int a = 0; a < 3; a++)
		axt[a] = axt_prev[a] + axt_seed[a] + axt_next[a];

	      string edit_script = make_edit_script(axt, '+');
	      script_writer.write(edit_script);

#ifdef PRINT_MATCHES
	      /*if ((rand() & 0xfff) == 0)*/ {
		for (int i = start0; i <= stop0; i++)
		  cerr << buf[i];
		cerr << endl;
		for (int i = start1; i <= stop1; i++)
		  cerr << seqs[id][i];
		cerr << endl;
		cerr << '+' << " " << stop0-start0 << endl;

		cerr << axt[0] << endl << axt[1] << endl << axt[2] << endl;

		cerr << "edit script: " << edit_script << endl;
		cerr << "length: " << edit_script.length() << "/" << stop0-start0+1 << endl;
	      }
#endif

#ifdef CHECK_AXT
	      // 1. axt1 = axt0==axt2
	      for (int x = 0; x < (int) axt[0].length(); x++)
		assert((axt[0][x] == axt[2][x] ? '|' : ' ') == axt[1][x]);
	      // 2. smush(axt0) = s1sub, smush(axt2) = s2sub (revcomp if necessary)
	      int check_pos;
	      check_pos = start0;
	      for (int x = 0; x < (int) axt[0].length(); x++)
		if (axt[0][x] != '-') {
		  assert(axt[0][x] == buf[check_pos]);
		  check_pos++;
		}
	      check_pos = start1;
	      for (int x = 0; x < (int) axt[2].length(); x++)
		if (axt[2][x] != '-') {
		  assert(axt[2][x] == seqs[id][check_pos]);
		  check_pos++;
		}
#endif
	      return match(start0, stop0, start1, stop1, id, '+',
			   script_writer.get_tot_script_chars() - edit_script.length());
	    }
	  }
	}
      }
      if (k_mer_freqs[k_mer_revcomp] < max_k_mer_freq) { // opp dir
	for (int j = 0; j < k_mer_freqs[k_mer_revcomp]; j++) {
	  int stop0 = i, start0 = stop0 - BLAST_K + 1;
	  int stop1 = k_mer_locs[k_mer_revcomp][j].pos,
	    start1 = stop1-BLAST_K+1, id = k_mer_locs[k_mer_revcomp][j].id;

#ifdef DEBUG_OUTPUT
	  cerr << "seed k-mer: ";
	  for (int x = start1; x <= stop1; x++)
	    cerr << seqs[id][x];
	  cerr << endl;
#endif
	  vector <string> axt_seed(3);
	  axt_seed[0] = seqs[id].substr(start1, BLAST_K);
	  axt_seed[1] = string(BLAST_K, '|');
	  axt_seed[2] = axt_seed[0];

	  if (attempt_ext(start0, -1, buf, buflen, stop1, 1, seqs[id]) +
	      attempt_ext(stop0, 1, buf, buflen, start1, -1, seqs[id])
	      > MIN_PROGRESS) {
	    vector <string> axt_prev = attempt_gapped_ext(stop0, 1, buf, buflen, start1, -1, seqs[id]);
	    vector <string> axt_next = attempt_gapped_ext(start0, -1, buf, buflen, stop1, 1, seqs[id]);
#ifdef DEBUG_OUTPUT
	    cerr << "aligned: " << stop0 - start0 << " (reverse)" << endl;
#endif
	    if (stop0 - start0 >= MIN_MATCH_LEN) {

	      // combine alignment strings
	      vector <string> axt(3);
	      for (int a = 0; a < 3; a++)
		axt[a] = axt_prev[a] + axt_seed[a] + axt_next[a];

	      string edit_script = make_edit_script(axt, '-');
	      script_writer.write(edit_script);

#ifdef PRINT_MATCHES
	      /*if ((rand() & 0xfff) == 0)*/ {
		for (int i = start0; i <= stop0; i++)
		  cerr << buf[i];
		cerr << endl;
		for (int i = start1; i <= stop1; i++)
		  cerr << seqs[id][i];
		cerr << endl;
		cerr << '-' << " " << stop0-start0 << endl;

		cerr << axt[0] << endl << axt[1] << endl << axt[2] << endl;

		cerr << "edit script: " << edit_script << endl;
		cerr << "length: " << edit_script.length() << "/" << stop0-start0+1 << endl;
	      }
#endif

#ifdef CHECK_AXT
	      // 1. axt1 = axt0==axt2
	      for (int x = 0; x < (int) axt[0].length(); x++)
		assert((axt[0][x] == axt[2][x] ? '|' : ' ') == axt[1][x]);
	      // 2. smush(axt0) = s1sub *** revcomp ***, smush(axt2) = s2sub
	      int check_pos;
	      check_pos = stop0;
	      for (int x = 0; x < (int) axt[0].length(); x++)
		if (axt[0][x] != '-') {
		  assert(axt[0][x] == base_pair_comp(buf[check_pos]));
		  check_pos--;
		}
	      check_pos = start1;
	      for (int x = 0; x < (int) axt[2].length(); x++)
		if (axt[2][x] != '-') {
		  assert(axt[2][x] == seqs[id][check_pos]);
		  check_pos++;
		}
#endif
	      return match(start0, stop0, start1, stop1, id, '-',
			   script_writer.get_tot_script_chars() - edit_script.length());
	    }
	  }
	}
      }
    }
  }
  return match();
}

void read_FASTA(const char *links_compressed_file, vector <string> &seqs,
		vector < vector <match> > &links, HalfByteWriter &script_writer) {
  vector <int> seqs_cum_lengths(1, 0);
  string fasta_header;
  bool getline_success = getline(cin, fasta_header);
  if (fasta_header[0] != '>') {
    cerr << "Error: First line of input file is missing FASTA header." << endl;
    return;
  }
  
  FILE *links_compressed_ptr = fopen(links_compressed_file, "wb"); // mustOpen would be better

  char buffer[3*MAX_CHUNK_SIZE];
  while (getline_success) { // last line should be fasta header; empty = EOF

    // links_compressed.txt: output length of fasta header, then header
    int header_len = fasta_header.length();
    fwrite(&header_len, sizeof(header_len), 1, links_compressed_ptr);
    fwrite(fasta_header.c_str(), sizeof(char), header_len, links_compressed_ptr);

#ifdef PRINT_OLD_LINKS_TXT
    cout << fasta_header << endl;
#endif
    vector <match> cur_record_matches;
    int buflen = 0, offset = 0; // number of bases processed so far
    int last_offset_printed = 0; timer.update_time();
    bool remaining_record_data = true;
    fill_buf(buffer, buflen, cin, remaining_record_data, getline_success,
	     fasta_header);
    while (buflen) {
      if (offset > last_offset_printed + 1000000) {
#ifndef PRINT_MATCHES
	cerr << offset << " " << timer.update_time() << endl;
#endif
	last_offset_printed = offset;
      }
      // look for matches until MAX_CHUNK_SIZE or end of buflen
      // note: the returned match will have no offset
      match cur_match = find_match(seqs, buffer, buflen, script_writer);

      // if match found:
      // - copy out pre-match chunk (+ OVERLAP) to unique seqs; save link
      // - save match link (append to cur_record_matches vector)
      // - move remnant of buffer (+ OVERLAP at end of match) to front
      // - refill buffer

      // else (no matches found):
      // - copy out MAX_CHUNK_SIZE (or entire buffer if at end); save link
      // - move remnant of buffer (+ OVERLAP at end) to front
      // - refill buffer
      
      int buf_used;
      if (cur_match.start0 != -1) { // found a match
	if (cur_match.start0 != 0) { // need to copy over pre-match chunk
	  int buf_copied = cur_match.start0 + OVERLAP;
	  cur_record_matches.push_back(copy_chunk(seqs, seqs_cum_lengths,
						  buffer, buf_copied, offset));
	}
	buf_used = cur_match.stop0 + 1;
	cur_match.start0 += offset; cur_match.stop0 += offset;
	cur_record_matches.push_back(cur_match);
      }
      else {
	buf_used = min(buflen, MAX_CHUNK_SIZE);
	cur_record_matches.push_back(copy_chunk(seqs, seqs_cum_lengths,
						buffer, buf_used, offset));
      }

      if (buf_used == buflen) {
	buflen = 0; // should only happen if no remaining record data
	offset += buf_used;
	assert(remaining_record_data == false);
      }
      else {
	buf_used -= OVERLAP; // back up by OVERLAP
	for (int i = buf_used; i < buflen; i++)
	  buffer[i-buf_used] = buffer[i];
	buflen -= buf_used;
	offset += buf_used;
      }
      fill_buf(buffer, buflen, cin, remaining_record_data, getline_success,
	       fasta_header);
    }
    links.push_back(cur_record_matches);

    // links_compressed.txt: output number of links, then link_info
    int num_links = cur_record_matches.size();
    fwrite(&num_links, sizeof(num_links), 1, links_compressed_ptr);

    for (int i = 0; i < (int) cur_record_matches.size(); i++) {
      const match &m = cur_record_matches[i];
      int cum_offset = seqs_cum_lengths[m.id];
      // outputs start1, stop1 in uniques_all_seqs concatenated coords
#ifdef PRINT_OLD_LINKS_TXT
      printf("%d-%d, %d-%d %c\n", m.start0, m.stop0, m.start1 + cum_offset,
	     m.stop1 + cum_offset, m.dir);
#endif
      // links_compressed.txt: output link_info
      LinkInfo link(m.start0, m.start1+cum_offset, m.start_script,
		    m.stop0-m.start0, m.stop1-m.start1);
      fwrite(&link, sizeof(link), 1, links_compressed_ptr);
    }

#ifdef GLOBAL_STATS
    tot_processed += offset;
    cerr << "*** total processed: " << tot_processed << endl;
    cerr << "    seqs len: " << tot_len << endl;
#endif
  }
  fclose(links_compressed_ptr);

}

void print_usage() {
  cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  cerr << "Arguments:" << endl;
  cerr << UNIQ_ALL_ARG_NUM << " - compressed db filename: uniques (e.g., uniques_all_seqs.fasta)" << endl;
  cerr << LINKS_ARG_NUM << " - compressed db filename: links (e.g., links_compressed.dat)" << endl;
  cerr << EDITS_ARG_NUM << " - compressed db filename: edit scripts (e.g., edit_scripts.dat)" << endl;
  cerr << SEPARATE_SEQS_ARG_NUM << " - uniques file format: 'smushed' OR 'separate' sequences" << endl;
  cerr << WINDOW_THRESH_ARG_NUM << " - window identity threshold for compression (default 85)" << endl;
  cerr << endl;
  cerr << "Reads sequences in fasta format from stdin." << endl;
  cerr << "Writes log to stdout." << endl;
  cerr << endl;
  cerr << "Example usage: ./gen_compressed_db uniques.fasta links.dat" << endl;
  cerr << "                   edit_scripts.dat smushed < test.fasta" << endl;
  cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  cerr << endl;
}

int main(int argc, char *argv[]) {

  if (argc <= SEPARATE_SEQS_ARG_NUM) {
    cerr << "not enough input args" << endl;
    print_usage();
    exit(1);
  }
  if (argc > WINDOW_THRESH_ARG_NUM) {
    sscanf(argv[WINDOW_THRESH_ARG_NUM], "%d", &WINDOW_IDENT_THRESH);
    if (!(70 <= WINDOW_IDENT_THRESH && WINDOW_IDENT_THRESH <= 100)) {
      cerr << "invalid window identity threshold (must be 70-100)" << endl;
      print_usage();
      exit(1);
    }
  }

  bool smush_uniques;
  if (strcmp(argv[SEPARATE_SEQS_ARG_NUM], "smushed") == 0)
    smush_uniques = true;
  else if (strcmp(argv[SEPARATE_SEQS_ARG_NUM], "separate") == 0)
    smush_uniques = false;
  else {
    cerr << "argument " << SEPARATE_SEQS_ARG_NUM << " has incorrect format" << endl;
    print_usage();
    exit(1);
  }

  memset(ismatch, 1, WINDOW_SIZE);
  
  for (int k_mer = 0; k_mer < (1<<(2*BLAST_K)); k_mer++) {
    k_mer_locs[k_mer].resize(MAX_K_MER_FREQ);
  }

  vector < vector <match> > links;

  timer.update_time(); double starttime = timer.curtime;

  HalfByteWriter script_writer(argv[EDITS_ARG_NUM]);

  vector <string> seqs;
  read_FASTA(argv[LINKS_ARG_NUM], seqs, links, script_writer);

  script_writer.flush_and_close();
  cerr << "total script chars: " << script_writer.get_tot_script_chars() << endl;
  
  timer.update_time();
  cerr << "time taken: " << timer.curtime - starttime << endl;

  fstream fout;
  mustOpenStream(fout, argv[UNIQ_ALL_ARG_NUM], ios::out);
  
  for (int i = 0; i < (int) seqs.size(); i++) {
    if (i == 0 || !smush_uniques)
      fout << ">" << i << endl; // fasta header
    fout << seqs[i] << endl;
  }
  fout.close();
  return 0;
}
