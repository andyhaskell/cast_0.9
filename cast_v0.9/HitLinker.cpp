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

//#include "gfClientLib.h" // for mustOpen

#include "CharTranslation.cpp"
#include "LinkInfo.cpp"
#include "DNAUtils.cpp"
#include "IOUtils.cpp"
#include "Timer.cpp"

//#define UNCOMPRESSED_CHECK

using namespace std;

const int CHUNK_SIZE = 30000;
const int MAX_EDIT_LEN = 10000;

struct rmatch {
  int start0, stop0, start1, stop1; // 0: str aligned to 1: ref in uniques db
  int id0; // which fasta entry this match comes from (>seq0, >seq1, >seq2, ...)
  char dir; // +/- alignment
  rmatch(int _start0=-1, int _stop0=-1, int _start1=-1, int _stop1=-1,
	 int _id0=-1, char _dir='-') :
    start0(_start0), stop0(_stop0), start1(_start1), stop1(_stop1), id0(_id0),
    dir(_dir) {};
};

struct LinkInfoID {
  LinkInfo link;
  int id0;
};

struct edit_info {
  bool is_subdel; // 0 if insert, 1 if subdel
  int last_dist;
  //string str;
  char str[MAX_EDIT_LEN];
  int strlen;
};

struct HitExpansion {
  string header_plus_dna;
  int offset, orig_length;
  char strand;
};

class HitLinker {
public:

#ifdef UNCOMPRESSED_CHECK
  void old_init(void) {
    ifstream fasta_in("Bacteria_BLAT_tests/Brucella/all_seqs3.fasta");
    //ifstream fasta_in("cat_flies.fasta");
    string line;
    ostringstream oss;
    bool first = true;
    while (getline(fasta_in, line)) {
      if (line[0] == '>') { // FASTA header
	//m_fasta_headers.push_back(line);
	if (!first) {
	  m_orig_seqs.push_back(oss.str());
	  oss.str("");
	}
	first = false;
      }
      else {
	oss << line;
      }
    }
    m_orig_seqs.push_back(oss.str());
    fasta_in.close();

    // read in the links
    ifstream link_in("links3.txt");
    //ifstream link_in("links.txt");
    rmatch r; r.id0 = -1;
    while (getline(link_in, line)) {
      if (line[0] == '>') { // FASTA header
	r.id0++;
      }
      else {
	sscanf(line.c_str(), "%d-%d, %d-%d %c\n", &r.start0, &r.stop0,
	       &r.start1, &r.stop1, &r.dir);
	int chunk_start = r.start1 / CHUNK_SIZE,
	  chunk_stop = r.stop1 / CHUNK_SIZE;
	for (int chunk = chunk_start; chunk <= chunk_stop; chunk++) {
	  if ((int) m_links_by_chunk.size() <= chunk)
	    m_links_by_chunk.resize(chunk+1);
	  m_links_by_chunk[chunk].push_back(r);
	}
      }
    }
    link_in.close();
  }
#endif

  HitLinker(const char *uniques_file, const char *links_file, const char *scripts_file,
	    int _hit_pad_length, char _header_dna_delim) {

    hit_pad_length = _hit_pad_length;
    header_dna_delim = _header_dna_delim;

#ifdef UNCOMPRESSED_CHECK
    old_init();
#endif

    Timer timer; timer.update_time();

    /***** read uniques fasta file *****/

    fstream uniques_in;
    mustOpenStream(uniques_in, uniques_file, ios::in);
    string line;
    while (getline(uniques_in, line)) {
      if (line[0] == '>')
	uniq_all_offsets.push_back(uniq_all.length());
      else
	uniq_all += line;
    }
    uniques_in.close();
#ifdef REPORT_TIMES
    cerr << "time for reading uniques: " << timer.update_time() << endl;
#endif

    /***** read compressed edit scripts file *****/

    fstream scripts_in;
    mustOpenStream(scripts_in, scripts_file, ios::in|ios::binary|ios::ate);
    int size = scripts_in.tellg();
    edit_scripts_packed = new char[size];
    scripts_in.seekg(0, ios::beg);
    scripts_in.read(edit_scripts_packed, size);
    scripts_in.close();
#ifdef REPORT_TIMES
    cerr << "time for reading edit scripts: " << timer.update_time() << endl;
#endif

    /***** read compressed links file *****/

    FILE *links_compressed_ptr = mustOpenFile(links_file, "rb");
    int header_len;
    LinkInfoID link_with_ID; link_with_ID.id0 = 0;
    while (fread(&header_len, sizeof(header_len), 1, links_compressed_ptr)) {
      //cout << "header_len: " << header_len << endl;
      char header_buf[header_len+1];
      fread(header_buf, sizeof(header_buf[0]), header_len, links_compressed_ptr);
      header_buf[header_len] = '\0';
      string header_str(header_buf);
      m_fasta_headers.push_back(header_str);
      //cout << header_str << endl;
      
      int num_links;
      fread(&num_links, sizeof(num_links), 1, links_compressed_ptr);
      for (int i = 0; i < num_links; i++) {
	fread(&link_with_ID.link, sizeof(link_with_ID.link), 1, links_compressed_ptr);
	int link_ind = m_links.size();
	m_links.push_back(link_with_ID);

	int chunk_start = link_with_ID.link.start1 / CHUNK_SIZE,
	  chunk_stop = (link_with_ID.link.start1 + link_with_ID.link.dist1) / CHUNK_SIZE;
	for (int chunk = chunk_start; chunk <= chunk_stop; chunk++) {
	  if ((int) m_link_inds_by_chunk.size() <= chunk)
	    m_link_inds_by_chunk.resize(chunk+1);
	  m_link_inds_by_chunk[chunk].push_back(link_ind);
	}
	/*
	printf("%d-%d, %d-%d, start %d\n", link.start0, link.start0+link.dist0,
	       link.start1, link.start1+link.dist1, link.start_script);
	string dest(link.dist0+1, '?');
	decode_edit_script(dest, link.start0, link);
	*/	
      }
      // last stop0 is length-1
      m_orig_lengths.push_back(link_with_ID.link.start0 + link_with_ID.link.dist0 + 1);
      //cerr << "original sequence length: " << m_orig_lengths.back() << endl;
      link_with_ID.id0++;
    }
    fclose(links_compressed_ptr);
#ifdef REPORT_TIMES
    cerr << "time for reading links: " << timer.update_time() << endl;
#endif
  }

  ~HitLinker() {
    delete[] edit_scripts_packed;
  }

#ifdef UNCOMPRESSED_CHECK
  vector < pair <string, int> > expand_hits(int id1, int start1, int stop1) {
    vector < pair <string, int> > ret;
    int chunk_start = start1 / CHUNK_SIZE, chunk_stop = stop1 / CHUNK_SIZE;
    for (int chunk = chunk_start; chunk <= chunk_stop; chunk++)
      for (int i = 0; i < (int) m_links_by_chunk[chunk].size(); i++) {
	rmatch &r = m_links_by_chunk[chunk][i];
	if (r.start1 <= stop1 && r.stop1 >= start1) { // intervals intersect

	  //printf("%d-%d, %d-%d %c\n", r.start0, r.stop0, r.start1, r.stop1, r.dir);

	  // can reduce this to the relevant interval -- unnecessary to look at the entire hit
	  int start0 = max(0, (r.dir == '+' ?
			       min(start1 + r.start0-r.start1, start1 + r.stop0-r.stop1) :
			       min(r.start0 + r.stop1 - stop1, r.stop0 - (stop1 - r.start1)))
			   - hit_pad_length);
	  int stop0 = (r.dir == '+' ?
		       max(stop1 + r.start0-r.start1, stop1 + r.stop0-r.stop1) :
		       max(r.stop0 - (start1 - r.start1), r.start0 + r.stop1 - start1))
	    + hit_pad_length;
	  // note this isn't capped because substr() doesn't mind overshooting
	  ret.push_back(make_pair(m_fasta_headers[r.id0] + '\n' +
				  m_orig_seqs[r.id0].substr(start0, stop0-start0+1),
				  start0));
	}
      }
    return ret;
  };
#endif

  /*
   * input:
   * - id1: which entry of uniques_all_seqs we're expanding (currently there's only one)
   * - start1, stop1: range to expand
   *
   * output: vector of entries corresponding to expanded links
   * - first: string containing sequence header + '\0' + dna
   * - second: offset (sequence coordinate corresponding to start of dna)
   *           total length of original sequence
   */
  vector <HitExpansion> expand_hits_packed(int id1, int start1, int stop1) {
    //cerr << "start1: " << start1 << " stop1: " << stop1 << endl;

    // convert to uniq_all coordinates
    start1 += uniq_all_offsets[id1];
    stop1 += uniq_all_offsets[id1];

    vector <HitExpansion> ret;
    int chunk_start = start1 / CHUNK_SIZE, chunk_stop = stop1 / CHUNK_SIZE;
    for (int chunk = chunk_start; chunk <= chunk_stop; chunk++)
      for (int i = 0; i < (int) m_link_inds_by_chunk[chunk].size(); i++) {
	int base_ind = m_link_inds_by_chunk[chunk][i];
	const LinkInfo &r = m_links[base_ind].link;
	int id0 = m_links[base_ind].id0;

	if (r.start1 <= stop1 && r.start1+r.dist1 >= start1) { // intervals intersect

	  int r_dir = get_link_dir(r);
	  /*
	  fprintf(stderr, "%d-%d, %d-%d %c\n", r.start0, r.start0+r.dist0,
		  r.start1, r.start1+r.dist1, r_dir == 1 ? '+' : '-');
	  */

	  // compute extent of the expanded region we want from the original sequence
	  int start0 = max(0, (r_dir == 1 ?
			       min(start1 + r.start0-r.start1,
				   start1 + r.start0+r.dist0-(r.start1+r.dist1)) :
			       min(r.start0 + r.start1+r.dist1 - stop1,
				   r.start0+r.dist0 - (stop1 - r.start1)))
			   - hit_pad_length);
	  int stop0 = min((r_dir == 1 ?
			   max(stop1 + r.start0-r.start1,
			       stop1 + r.start0+r.dist0-(r.start1+r.dist1)) :
			   max(r.start0+r.dist0 - (start1 - r.start1),
			       r.start0 + r.start1+r.dist1 - start1))
			  + hit_pad_length, m_orig_lengths[id0]-1);

	  string orig_str(stop0-start0+1, '?');
	  for (int ind = base_ind; ind >= 0 && m_links[ind].id0==id0; ind--) {
	    const LinkInfo &cur_link = m_links[ind].link;
	    if (!(cur_link.start0 <= stop0 && cur_link.start0+cur_link.dist0 >= start0)) break;
	    decode_edit_script(orig_str, start0, cur_link);
	  }
	  for (int ind = base_ind+1; ind < (int) m_links.size() && m_links[ind].id0==id0; ind++) {
	    const LinkInfo &cur_link = m_links[ind].link;
	    if (!(cur_link.start0 <= stop0 && cur_link.start0+cur_link.dist0 >= start0)) break;
	    decode_edit_script(orig_str, start0, cur_link);
	  }
	  HitExpansion expansion;
	  expansion.header_plus_dna = m_fasta_headers[id0] + header_dna_delim + orig_str;
	  expansion.offset = start0;
	  expansion.orig_length = m_orig_lengths[id0];
	  expansion.strand = (r_dir == 1 ? '+' : '-');
	  ret.push_back(expansion);
	}
      }
    return ret;
  };

private:
  int hit_pad_length;
  char header_dna_delim;

#ifdef UNCOMPRESSED_CHECK
  vector <string> m_orig_seqs; // OLD
  vector < vector <rmatch> > m_links_by_chunk; // OLD
#endif

  // info about original sequences
  vector <string> m_fasta_headers;
  vector <int> m_orig_lengths;

  // info about uniques and links
  string uniq_all; // uniques, all smushed together
  vector <int> uniq_all_offsets; // starts of unique seqs >0, >1, >2, ... in uniq_all
  char *edit_scripts_packed; // edit scripts, all smushed together
  
  vector <LinkInfoID> m_links; // in consecutive order!!
  // rmatch contains link_info along with id0 indexed into list of fasta headers
  vector < vector <int> > m_link_inds_by_chunk;

  edit_info edit; // contains edit string buffer (only allocate once)


  char get_script_char(int pos) {
    if (pos&1) // odd
      return half_byte_to_char[edit_scripts_packed[pos>>1]&0xf];
    else
      return half_byte_to_char[((unsigned char) edit_scripts_packed[pos>>1])>>4];
  }

  int get_link_dir(const LinkInfo &link) {
    if (link.start_script == -1) return 1;
    return get_script_char(link.start_script) == '0' ? 1 : -1;    
  }

  // reads one edit worth of info (or fails if current decoded char is digit)
  // moves pos accordingly and returns info in edit param
  bool next_edit(int &pos) {
    if (isdigit(get_script_char(pos)))
      return false;
    edit.is_subdel = get_script_char(pos++) == 's';
    edit.last_dist = 0;
    edit.strlen = 0;
    while (isdigit(get_script_char(pos))) {
      edit.last_dist *= 8; // octal encoding
      edit.last_dist += get_script_char(pos++) - '0';
    }
    while (isupper(get_script_char(pos)) || get_script_char(pos) == '-')
      edit.str[edit.strlen++] = get_script_char(pos++);
    return true;
  }

  char comp_base_if_needed(char base, int dir) {
    return dir == 1 ? base : base_pair_comp(base);
  }

  bool intervals_isect(int start1, int end1, int start2, int end2) {
    return start1 <= end2 && start2 <= end1;
  }

  void decode_edit_script(string &dest, int dest0_coord, const LinkInfo &link) {

    int dest_len = dest.length();
    int i0, dir0;

    if (link.start_script == -1) { // no edits; verbatim copy
      i0 = link.start0 - dest0_coord;
      for (int i1 = link.start1; i1 <= link.start1+link.dist1; i0++, i1++)
	if (0 <= i0 && i0 < dest_len)
	  dest[i0] = uniq_all[i1];
      return;
    }
     
    int script_pos = link.start_script;
    int orig_pos = link.start1; // position in original string... "original"=str1?
    int last_edit_str_len = 0; // length of last edit str
    if (get_script_char(script_pos++) == '0') {
      dir0 = 1;
      i0 = link.start0 - dest0_coord;

      while (next_edit(script_pos)) {
	// chunk after previous edit
	//str += uniq_all.substr(orig_pos, edit.last_dist - last_edit_str_len);
	int xmin = -i0;
	int xmax = dest_len-i0;
	for (int x = max(0, xmin); x < min(edit.last_dist - last_edit_str_len, xmax); x++)
	  dest[i0+x] = uniq_all[x+orig_pos];
	i0 += edit.last_dist - last_edit_str_len;

	// update position in original string
	orig_pos += edit.last_dist - last_edit_str_len;

	// append replacement string in edit script; get rid of dashes
	//str += edit.str[i];
	for (int i = 0; i < edit.strlen; i++)
	  if (edit.str[i] != '-') {
	    if (0 <= i0 && i0 < dest_len)
	      dest[i0] = edit.str[i];
	    i0++;
	  }

	// skip subdel along original string
	if (edit.is_subdel) orig_pos += edit.strlen;

	last_edit_str_len = edit.strlen;

	if (i0 >= dest_len) return; // passed end of region of interest: done!
      }

    }
    else {
      dir0 = -1;
      i0 = link.start0+link.dist0 - dest0_coord;

      while (next_edit(script_pos)) {
	// chunk after previous edit
	//str += uniq_all.substr(orig_pos, edit.last_dist - last_edit_str_len);
	int xmin = i0-dest_len+1;
	int xmax = i0+1;
	for (int x = max(0, xmin); x < min(edit.last_dist - last_edit_str_len, xmax); x++)
	  dest[i0-x] = base_pair_comp(uniq_all[x+orig_pos]);
	i0 -= edit.last_dist - last_edit_str_len;

	// update position in original string
	orig_pos += edit.last_dist - last_edit_str_len;

	// append replacement string in edit script; get rid of dashes
	//str += edit.str[i];
	for (int i = 0; i < edit.strlen; i++)
	  if (edit.str[i] != '-') {
	    if (0 <= i0 && i0 < dest_len)
	      dest[i0] = base_pair_comp(edit.str[i]);
	    i0--;
	  }

	// skip subdel along original string
	if (edit.is_subdel) orig_pos += edit.strlen;

	last_edit_str_len = edit.strlen;

	if (i0 < 0) return; // passed end of region of interest: done!
      }

    }

    // chunk after last edit
    //str += uniq_all.substr(orig_pos, link.start1+link.dist1+1-orig_pos);
    if (dir0 == 1 && i0 < dest_len || dir0 == -1 && i0 >= 0) {
      for (int i1 = orig_pos; i1 <= link.start1+link.dist1; i1++) {
	if (0 <= i0 && i0 < dest_len)
	  dest[i0] = comp_base_if_needed(uniq_all[i1], dir0);
	i0 += dir0;
      }
    }
  }
};
