/*
 * args:
 * 1 - file containing coarse BLAT results (psl)
 * 2 - compressed db file: uniques (fasta)
 * 3 - compressed db file: links (dat)
 * 4 - compressed db file: edit scripts (dat)
 * 5 - query file (fasta)
 * 6 - final output file (psl)
 * ? - additional args such as minScore?
 */

#include "common.h"
#include "memalloc.h"
#include "linefile.h"
#include "bits.h"
#include "hash.h"
#include "dnautil.h"
#include "dnaseq.h"
#include "fa.h"
#include "nib.h"
#include "twoBit.h"
#include "psl.h"
#include "sig.h"
#include "options.h"
#include "obscure.h"
#include "genoFind.h"
#include "trans3.h"
#include "gfClientLib.h"

#include "HitLinker.cpp"
#include "AlignTrimmer.cpp"

//#define CABLAT_DEBUG_OUTPUT

using namespace std;

// input arguments
const int COARSE_OUT_ARG_NUM = 1;
const int UNIQ_ALL_ARG_NUM = 2;
const int LINKS_ARG_NUM = 3;
const int EDITS_ARG_NUM = 4;
const int QUERY_ARG_NUM = 5;
const int OUT_ARG_NUM = 6;


/********* default global vars copied from blat.c *********/

struct gfOutput *gvo;		/* Overall output controller */

/* Variables shared with other modules.  Set in this module, read only
 * elsewhere. */
char *databaseName;		/* File name of database. */
int databaseSeqCount = 0;	/* Number of sequences in database. */
unsigned long databaseLetters = 0;	/* Number of bases in database. */

/* Variables that can be set from command line. */
int tileSize = 11;
int stepSize = 0;	/* Default (same as tileSize) */
int minMatch = 2;
int minScore = 30;
int maxGap = 2;
int repMatch = 1024*4;
int dotEvery = 0;
boolean oneOff = FALSE;
boolean noHead = FALSE;
boolean trimA = FALSE;
boolean trimHardA = FALSE;
boolean trimT = FALSE;
boolean fastMap = FALSE;
char *makeOoc = NULL;
char *ooc = NULL;
enum gfType qType = gftDna;
enum gfType tType = gftDna;
char *mask = NULL;
char *repeats = NULL;
char *qMask = NULL;
double minRepDivergence = 15;
double minIdentity = 90;
char outputFormat[] = "psl";

/********* from jkOwnLib/gfBlatLib.c *********/

int scoreAli(struct ffAli *ali, boolean isProt, 
	enum ffStringency stringency, 
	struct dnaSeq *tSeq, struct trans3 *t3List)
/* Score alignment. */
{
int (*scoreFunc)(char *a, char *b, int size);
struct ffAli *ff, *nextFf;
int score = 0;
if (isProt) 
    scoreFunc = aaScoreMatch;
else
    scoreFunc = dnaScoreMatch;
for (ff = ali; ff != NULL; ff = nextFf)
    {
    nextFf = ff->right;
    score += scoreFunc(ff->nStart, ff->hStart, ff->nEnd-ff->nStart);
    if (nextFf != NULL)
        {
	int nhStart = trans3GenoPos(nextFf->hStart, tSeq, t3List, FALSE);
	int ohEnd = trans3GenoPos(ff->hEnd, tSeq, t3List, TRUE);
	int hGap = nhStart - ohEnd;
	int nGap = nextFf->nStart - ff->nEnd;
	score -= ffCalcGapPenalty(hGap, nGap, stringency);
	}
    }
return score;
}





bool check_usage(int argc, char *argv[]) {
  return argc == 7; // validate existence of files? additional args?
}

void print_usage(void) {
  cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  cerr << "Arguments: " << endl;
  cerr << COARSE_OUT_ARG_NUM << " - file containing coarse BLAT results (psl)" << endl;
  cerr << UNIQ_ALL_ARG_NUM << " - compressed db file: uniques (fasta)" << endl;
  cerr << LINKS_ARG_NUM << " - compressed db file: links (dat)" << endl;
  cerr << EDITS_ARG_NUM << " - compressed db file: edit scripts (dat)" << endl;
  cerr << QUERY_ARG_NUM << " - query file (fasta)" << endl;
  cerr << OUT_ARG_NUM << " - final output file (psl)" << endl;
  cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  //cerr << "? - additional args such as minScore?" << endl;
}

struct CoarseHit {
  unsigned int qInd, qStart, qEnd, tInd, tStart, tEnd;
  char strand;
};


void process_seq(string &header, string &seq, int cur_ind, vector <string> &query_names_with_null,
		 vector <string> &query_strings, map <string, int> &query_names_to_inds) {
  string name = header.substr(1);
  query_names_with_null.push_back(name + '\0'); // add null-termination
  query_strings.push_back(seq);
  query_names_to_inds[name] = cur_ind;
}

void read_query_file(const char *query_file, vector <string> &query_names_with_null,
		     vector <string> &query_strings, map <string, int> &query_names_to_inds) {
  fstream query_in;
  mustOpenStream(query_in, query_file, ios::in);
  string line, header, seq;
  int cur_ind = -1;
  while (getline(query_in, line)) {
    if (line[0] == '>') {
      if (cur_ind != -1) {
	process_seq(header, seq, cur_ind, query_names_with_null, query_strings,
		    query_names_to_inds);
	seq = "";
      }
      header = line;
      cur_ind++;
    }
    else
      seq += line;
  }
  process_seq(header, seq, cur_ind, query_names_with_null, query_strings, query_names_to_inds);
  query_in.close();
}
  
vector <dnaSeq> make_qSeqs(vector <string> &query_names_with_null,
			   vector <string> &query_strings) {
  vector <dnaSeq> qSeqs;
  for (int i = 0; i < (int) query_strings.size(); i++) {
    dnaSeq qSeq;
    qSeq.name = &query_names_with_null[i][0];
    qSeq.dna = &query_strings[i][0];
    qSeq.size = query_strings[i].length();
    qSeqs.push_back(qSeq);
  }
  return qSeqs;
}


vector <CoarseHit> read_coarse_output(const char *coarse_output_file,
				      map <string, int> &query_names_to_inds) {

  // relevant fields of psl file format (0-based):
  //  9: qName
  // 11: qStart
  // 12: qEnd
  // 15: tStart
  // 16: tEnd

  vector <CoarseHit> coarse_hits;

  fstream coarse_in;
  mustOpenStream(coarse_in, coarse_output_file, ios::in);
  string line;
  while (getline(coarse_in, line)) {
    vector <int> field_starts(1, 0);
    for (int i = 0; i < (int) line.length(); i++)
      if (line[i] == '\t')
	field_starts.push_back(i+1);
    CoarseHit hit;
    hit.qInd = query_names_to_inds[line.substr(field_starts[9],
					       field_starts[10]-1-field_starts[9])];
    hit.strand = line[field_starts[8]];
    sscanf(line.c_str()+field_starts[11], "%d", &hit.qStart);
    sscanf(line.c_str()+field_starts[12], "%d", &hit.qEnd);
    sscanf(line.c_str()+field_starts[13], "%d", &hit.tInd); // assumes target header is just '>###'
    sscanf(line.c_str()+field_starts[15], "%d", &hit.tStart);
    sscanf(line.c_str()+field_starts[16], "%d", &hit.tEnd);
    coarse_hits.push_back(hit);
  }
  return coarse_hits;
}

int main(int argc, char *argv[]) {

  if (!check_usage(argc, argv)) {
    print_usage();
    exit(1);
  }

  /********** set up output controller **********/

  databaseName = "testdb"; // todo... do this and databaseSeqCount, databaseLetters matter?
  boolean tIsProt = FALSE;
  boolean qIsProt = FALSE;

  FILE *f = mustOpen(argv[OUT_ARG_NUM], "w");

  gvo = gfOutputAny(outputFormat, minIdentity*10, qIsProt, tIsProt, noHead, 
		    databaseName, databaseSeqCount, databaseLetters, minIdentity, f);

  /********** read query file and build map from query names to indices **********/

  vector <string> query_names_with_null;
  vector <string> query_strings;
  map <string, int> query_names_to_inds;
  read_query_file(argv[QUERY_ARG_NUM], query_names_with_null, query_strings, query_names_to_inds);
  // kind of unsafe... the pointers reference data in query_names and query_strings
  vector <dnaSeq> qSeqs = make_qSeqs(query_names_with_null, query_strings);

  /********** read coarse output **********/

  vector <CoarseHit> coarse_hits = read_coarse_output(argv[COARSE_OUT_ARG_NUM],
						      query_names_to_inds);

  /********** set up HitLinker, expand hits, and perform fine search with blat **********/

  HitLinker hit_linker(argv[UNIQ_ALL_ARG_NUM], argv[LINKS_ARG_NUM], argv[EDITS_ARG_NUM], 10, '\0');

  AlignTrimmer alignTrimmer;

  int num_ffFind_calls = 0;

  // iterate through coarse hits
  for (int i = 0; i < (int) coarse_hits.size(); i++) {
    CoarseHit &hit = coarse_hits[i];
    dnaSeq &qSeq = qSeqs[hit.qInd];

#ifdef CABLAT_DEBUG_OUTPUT
    cerr << "hit " << i << " query" << hit.qInd << endl;
#endif
    if (hit.tEnd - hit.tStart > 1000) {
#ifdef CABLAT_DEBUG_OUTPUT
      cerr << "ignoring hit " << i << " query" << hit.qInd << ": too long" << endl;
#endif
      continue;
    }

    // expand hits to target regions of original sequences
    vector <HitExpansion> hit_expansions =
      hit_linker.expand_hits_packed(hit.tInd, hit.tStart, hit.tEnd);

    // iterate through hit expansions (decompressed target regions)
    for (int j = 0; j < (int) hit_expansions.size(); j++) {
      HitExpansion &target = hit_expansions[j];
      int null_pos = target.header_plus_dna.find('\0'); // termination of target name
      int offset = target.offset;
      int orig_length = target.orig_length;

      boolean qIsRc = FALSE;
      if (hit.strand != target.strand) {
#ifdef CABLAT_DEBUG_OUTPUT
	cerr << "reverse complementing query" << endl;
#endif
	qIsRc = TRUE;
	reverseComplement(qSeq.dna, qSeq.size);
      }

      dnaSeq tSeq;
      tSeq.name = &target.header_plus_dna[1];
      tSeq.dna = &target.header_plus_dna[null_pos+1];
      tSeq.size = target.header_plus_dna.length() - (null_pos+1);
      boolean tIsRc = FALSE;

      // might be better style to use:
      // - struct dnaSeq *newDnaSeq(DNA *dna, int size, char *name);
      // - void freeDnaSeq(struct dnaSeq **pSeq);
      // ... but some unnecessary copying this way

      pair <int, int> trim_lengths;// = alignTrimmer.trim(&qSeq.dna[0], qSeq.size, &tSeq.dna[0], tSeq.size);

      if (trim_lengths.first != -100)  {

      //cerr << "trim: " << trim_lengths.first << endl;
      trim_lengths.first = max(trim_lengths.first, 0);
      trim_lengths.second = max(trim_lengths.second, 0);

      char *tSeqStart = &tSeq.dna[0] + trim_lengths.first;
      char *tSeqEnd = &tSeq.dna[0] + tSeq.size - trim_lengths.second;
      ffAli *ff = ffFind(&qSeq.dna[0], &qSeq.dna[0] + qSeq.size,
			 tSeqStart, tSeqEnd, ffCdna);
      num_ffFind_calls++;

      int score = scoreAli(ff, FALSE, ffCdna, &tSeq, NULL);

#ifdef CABLAT_DEBUG_OUTPUT
      cerr << j << " offset: " << offset << " length: " << tSeq.size << " score: " << score;
      fprintf(stderr, "\t%s\n", tSeq.name);
#endif      
      // todo: does gvo->out test this?
      // gfOut.c:savePslx seems to do output including thresholding... if gvo->out calls this,
      //   then no need to do scoreAli... but minScore doesn't seem to have much effect      
      if (score >= minScore) {
	gvo->out(tSeq.name, orig_length, offset/*+trim_lengths.first*/, ff, &tSeq, NULL, &qSeq,
		 qIsRc, tIsRc, ffCdna, minScore, gvo);
	/*
	if (trim_lengths.first == -100)  {
	  cerr << "hit " << i << " query" << hit.qInd << endl;
	  cerr << "query:  ";
	  fwrite(&qSeq.dna[0], 1, qSeq.size, stderr);
	  cerr << endl;
	  cerr << "target: ";
	  fwrite(&tSeq.dna[0], 1, tSeq.size, stderr);
	  cerr << endl;
	}
	*/
      }

      ffFreeAli(&ff);
      }

      if (qIsRc) // reset qSeq string
	reverseComplement(qSeq.dna, qSeq.size);
    }
  }

  cerr << "num_ffFind_calls: " << num_ffFind_calls << endl;

  carefulClose(&f);

  return 0;
}
