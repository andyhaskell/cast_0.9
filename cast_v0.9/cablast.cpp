/*
input:
------

-program   blastn (others not tested)
-fulldb    original uncompressed blast database created by makeblastdb; no file ext
           used as a reference to run blast on and compare to cablast run
-db        uniques database created by makeblastdb on uniques from gen_compressed_db; no file ext
-uniques   uniques FASTA file from gen_compressed_db: unique sequence chunks concatenated together
-links     links DAT file from gen_compressed_db
-edits     edits DAT file from gen_compressed_db
-in        query FASTA file
-evalue    evalue for final output/scoring
-coarse_evalue   evalue for coarse search

output:
-------

summary of cablast vs. blast comparison to stdout

*/


//#define VERBOSE
#define COMPARE_RESULTS
//#define PRINT_MISSES

#include <cstring>

#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <corelib/ncbienv.hpp>
#include <corelib/ncbiargs.hpp>

#include <objmgr/object_manager.hpp>

#include <objects/seqalign/Seq_align_set.hpp>

#include <algo/blast/api/sseqloc.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/bl2seq.hpp> // added
#include <algo/blast/api/uniform_search.hpp>
#include <algo/blast/api/blast_types.hpp>
#include <algo/blast/api/blast_aux.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>
#include <algo/blast/api/blast_prot_options.hpp>

#include <algo/blast/blastinput/blast_input.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>

#include "Timer.cpp"
#include "HitLinker.cpp"

USING_NCBI_SCOPE;
USING_SCOPE(blast);


/*
 * Boilerplate command-line argument processing code automatically generated
 * by new_project.sh script in NCBI C++ Toolkit
 *
 *     Authors:  Denis Vakatov, Vladimir Ivanov)
 */


/////////////////////////////////////////////////////////////////////////////
//  CCablastApplication::


class CCablastApplication : public CNcbiApplication
{
private:
    Timer timer;
    virtual void Init(void);
    virtual int  Run(void);
    virtual void Exit(void);
    CSearchResultSet full_blast(TSeqLocVector query_loc, CRef<CBlastOptionsHandle> opts_handle, Uint8 *full_db_size);
    void ProcessCommandLineArgs(CRef<CBlastOptionsHandle> opts_handle);

};


/////////////////////////////////////////////////////////////////////////////
//  Init test for all different types of arguments


void CCablastApplication::Init(void)
{
    // Create command-line argument descriptions class
    auto_ptr<CArgDescriptions> arg_desc(new CArgDescriptions);

    // Specify USAGE context
    arg_desc->SetUsageContext(GetArguments().GetProgramBasename(), "CaBLAST vs. BLAST comparison program");

    arg_desc->AddDefaultKey
        ("program", "ProgramName",
         "One of blastn, megablast, disc_megablast, blastp, blastx, tblastn, tblastx, rpsblast",
         CArgDescriptions::eString, "blastn");
    arg_desc->SetConstraint
        ("program", &(*new CArgAllow_Strings,
                "blastn", "megablast", "disc_megablast", "blastp", "blastx", "tblastn", "tblastx", "rpsblast"));

    arg_desc->AddKey
        ("fulldb", "FullDatabase",
         "Name of full (uncompressed) database created by makeblastdb, without file extension",
         CArgDescriptions::eString);

    arg_desc->AddKey
        ("db", "Database",
         "Name of database created by running makeblastdb on uniques (from gen_compressed_db)",
         CArgDescriptions::eString);

    arg_desc->AddKey
        ("uniques", "Uniques",
         "FASTA file containing unique sequence chunks (concatenated together)",
         CArgDescriptions::eString);

    arg_desc->AddKey
        ("links", "Links",
         "DAT file containing link information (from gen_compressed_db)",
         CArgDescriptions::eString);

    arg_desc->AddKey
        ("edits", "EditScripts",
         "DAT file containing edit scripts (from gen_compressed_db)",
         CArgDescriptions::eString);

    arg_desc->AddKey("in", "QueryFile",
		     "FASTA file containing queries", CArgDescriptions::eInputFile);
    /*
    arg_desc->AddKey("out", "Outputfile",
		     "Output file", CArgDescriptions::eOutputFile);
    */
    arg_desc->AddDefaultKey("evalue", "evalue",
                        "E-value threshold for final hits", CArgDescriptions::eDouble, "1e-30");

    arg_desc->AddDefaultKey("coarse_evalue", "coarse_evalue",
                        "E-value threshold for coarse search against uniques database", CArgDescriptions::eDouble, "1e-20");

    // Setup arg.descriptions for this application
    SetupArgDescriptions(arg_desc.release());
}

/// Modify BLAST options from defaults based upon command-line args.
///
/// @param opts_handle already created CBlastOptionsHandle to modify [in]
void CCablastApplication::ProcessCommandLineArgs(CRef<CBlastOptionsHandle> opts_handle)
{
	CArgs args = GetArgs();

        // Expect value is a supported option for all flavors of BLAST.
        if(args["evalue"].AsDouble()) {
	  cout << "Full BLAST evalue: " << args["evalue"].AsDouble() << endl;
          opts_handle->SetEvalueThreshold(args["evalue"].AsDouble());
	}
        
        // The first branch is used if the program is blastn or a flavor of megablast
        // as reward and penalty is a valid option.
        //
        // The second branch is used for all other programs except rpsblast as matrix
        // is a valid option for blastp and other programs that perform protein-protein
        // comparisons.
        //
	/*
        if (CBlastNucleotideOptionsHandle* nucl_handle =
              dynamic_cast<CBlastNucleotideOptionsHandle*>(&*opts_handle)) {

              if (args["reward"].AsInteger())
                nucl_handle->SetMatchReward(args["reward"].AsInteger());
            
              if (args["penalty"].AsInteger())
                nucl_handle->SetMismatchPenalty(args["penalty"].AsInteger());
        }
        else if (CBlastProteinOptionsHandle* prot_handle =
               dynamic_cast<CBlastProteinOptionsHandle*>(&*opts_handle)) {
              if (args["matrix"]) 
                prot_handle->SetMatrixName(args["matrix"].AsString().c_str());
        }
	*/

        return;
}


CSearchResultSet CCablastApplication::full_blast(TSeqLocVector query_loc, CRef<CBlastOptionsHandle> opts, Uint8 *full_db_size) {

    const CArgs& args = GetArgs();

    EProgram program = ProgramNameToEnum(args["program"].AsString());

    bool db_is_aa = (program == eBlastp || program == eBlastx ||
                     program == eRPSBlast || program == eRPSTblastn);

    const CSearchDatabase target_db(args["fulldb"].AsString(),
        db_is_aa ? CSearchDatabase::eBlastDbIsProtein : CSearchDatabase::eBlastDbIsNucleotide);

    *full_db_size = target_db.GetSeqDb()->GetTotalLength();
    cout << "Full db size: " << *full_db_size << endl;

    CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(query_loc));
    
    timer.update_time();
    CLocalBlast blaster(query_factory, opts, target_db);
    cout << "Time spent initializing full BLAST: " << timer.update_time() << endl;

    timer.update_time();
    CSearchResultSet results = *blaster.Run();
    cout << "Full BLAST search time: " << timer.update_time() << endl;

    /* print hits
    for (unsigned int i = 0; i < results.GetNumResults(); i++) {
         CConstRef<CSeq_align_set> sas = results[i].GetSeqAlign();
         //cout << MSerial_AsnText << *sas;

	 const list < CRef <CSeq_align> > &seqAlignList = sas->Get();
	 ITERATE(list < CRef <CSeq_align> >, seqAlign_it, seqAlignList) {
	   string label;
	   //(*seqAlign_it)->GetSeq_id(1).GetLabel(&label);
	   cout << "full blast: " << (*seqAlign_it)->GetSeq_id(1).GetSeqIdString() << " hit at loc " << (*seqAlign_it)->GetSeqStart(1)
		<< "-" << (*seqAlign_it)->GetSeqStop(1) << " " << (*seqAlign_it)->GetSeq_id(1).AsFastaString() << endl;
	 }
    }
    */
    return results;

}



/////////////////////////////////////////////////////////////////////////////
//  Run test (printout arguments obtained from command-line)


int CCablastApplication::Run(void)
{
    // Get arguments
    const CArgs& args = GetArgs();

    EProgram program = ProgramNameToEnum(args["program"].AsString());

    bool db_is_aa = (program == eBlastp || program == eBlastx ||
                     program == eRPSBlast || program == eRPSTblastn);

    CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));

    ProcessCommandLineArgs(opts);
    
    opts->Validate();  // Can throw CBlastException::eInvalidOptions for invalid option.

    // This will dump the options to stderr.
    // opts->GetOptions().DebugDumpText(cerr, "opts", 1);

    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
         throw std::runtime_error("Could not initialize object manager");
    }

    const bool is_protein =
        !!Blast_QueryIsProtein(opts->GetOptions().GetProgramType());
    SDataLoaderConfig dlconfig(is_protein);
    CBlastInputSourceConfig iconfig(dlconfig);
    CBlastFastaInputSource fasta_input(args["in"].AsInputFile(), iconfig);
    CScope scope(*objmgr);

    CBlastInput blast_input(&fasta_input);

    TSeqLocVector query_loc = blast_input.GetAllSeqLocs(scope);

    Uint8 full_db_size;
    CSearchResultSet full_results = full_blast(query_loc, opts, &full_db_size);

    // output number of full BLAST hits
    int num_full_results = 0;
    for (unsigned int i = 0; i < full_results.GetNumResults(); i++) {
      const list < CRef <CSeq_align> > &full_seqAlignList =
	full_results[i].GetSeqAlign()->Get();
      num_full_results += full_seqAlignList.size();
    }
    cout << "num_full_results: " << num_full_results << endl;
    cout << endl;

    cout << "---------- Coarse BLAST ----------" << endl;

    if(args["coarse_evalue"].AsDouble()) {
      cout << "Coarse evalue: " << args["coarse_evalue"].AsDouble() << endl;
      opts->SetEvalueThreshold(args["coarse_evalue"].AsDouble());
    }

    CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(query_loc));

    const CSearchDatabase target_db(args["db"].AsString(),
        db_is_aa ? CSearchDatabase::eBlastDbIsProtein : CSearchDatabase::eBlastDbIsNucleotide);

    cout << "Uniques db size: " << target_db.GetSeqDb()->GetTotalLength() << endl;

    opts->SetDbLength(full_db_size);
    cout << "Setting DbLength parameter to full db size: " << full_db_size << endl;

    timer.update_time();
    CLocalBlast blaster(query_factory, opts, target_db);
    cout << "Time spent initializing coarse BLAST: " << timer.update_time() << endl;

    timer.update_time();
    CSearchResultSet results = *blaster.Run();
    cout << "Coarse BLAST search time: " << timer.update_time() << endl;

    // Get warning messages.
    for (unsigned int i = 0; i < results.GetNumResults(); i++) 
    {
        TQueryMessages messages = results[i].GetErrors(eBlastSevWarning);
        if (messages.size() > 0)
        {
            CConstRef<CSeq_id> seq_id = results[i].GetSeqId();
            if (seq_id.NotEmpty())
                cerr << "ID: " << seq_id->AsFastaString() << endl;
            else
                cerr << "ID: " << "Unknown" << endl;

            ITERATE(vector<CRef<CSearchMessage> >, it, messages)
                cerr << (*it)->GetMessage() << endl;
        }
    }
    
    //CNcbiOstream& out = args["out"].AsOutputFile();

    cout << endl;
    cout << "---------- Fine BLAST ----------" << endl;

    timer.update_time();
    HitLinker HL(args["uniques"].AsString().c_str(), args["links"].AsString().c_str(),
		 args["edits"].AsString().c_str(), 50, '\n');
    cout << "Preparing for link-tracing... HitLinker setup time: " << timer.update_time() << endl;


    if (args["evalue"].AsDouble()) {
      cout << "Resetting evalue to " << args["evalue"].AsDouble()
	   << " for fine BLAST" << endl;
      opts->SetEvalueThreshold(args["evalue"].AsDouble());
    }

    int full_ctr = 0, missed_ctr = 0;
    bool punted = false;
    // should be the same /size/ as full_results
    for (unsigned int i = 0; i < results.GetNumResults(); i++) {

      /*
      cout << "starting refinement for query " << i
	   << " (" << full_results[i].GetSeqAlign()->Get().size()
	   << " full BLAST hits)" << endl;
      */
      //timer.update_time();

      // get coarse BLAST results for query i
      CConstRef<CSeq_align_set> sas = results[i].GetSeqAlign();
      const list < CRef <CSeq_align> > &seqAlignList = sas->Get();

      // punt if tons of hits (for some reason BLAST is slow on many sequences)
      if (seqAlignList.size() >= 1000) {
	cout << "WARNING: skipping query " << i << ": "
	     << seqAlignList.size() << " coarse BLAST hits ("
	     << full_results[i].GetSeqAlign()->Get().size()
	     << " full BLAST hits)" << endl;
	punted = true;
	//continue;
      }

      // if no coarse BLAST results, record all full BLAST hits as missed
      if (seqAlignList.empty()) {
	const list < CRef <CSeq_align> > &full_seqAlignList =
	  full_results[i].GetSeqAlign()->Get();
#ifdef COMPARE_RESULTS
	if (full_seqAlignList.size()) {
	  
	  ITERATE(list < CRef <CSeq_align> >, seqAlign_it, full_seqAlignList) {
	    full_ctr++; missed_ctr++;
#ifdef PRINT_MISSES
	    string label;
	    string fasta_prefix = (*seqAlign_it)->GetSeq_id(1).AsFastaString();
	    int start = (*seqAlign_it)->GetSeqStart(1);
	    int end = (*seqAlign_it)->GetSeqStop(1);
	    cout << "full blast: hit from " << fasta_prefix
		 << " at pos " << start << "-" << end << " ";
	    cout << "Missed. ";
	    double evalue;
	    (*seqAlign_it)->GetNamedScore(CSeq_align::eScore_EValue, evalue);
	    cout << evalue << endl;
#endif
	  }
#ifdef PRINT_MISSES
	  cout << "missed all " << full_seqAlignList.size()
	       << " from query " << i << endl;
#endif
	}
#endif      
	continue;
      }

      // follow links to create expanded list of candidate hits
      vector <HitExpansion> all_expansions;
      // iterate through coarse BLAST hits
      ITERATE(list < CRef <CSeq_align> >, seqAlign_it, seqAlignList) {
	/*
	  cout << "start0-stop0 (query): id "
	  << (*seqAlign_it)->GetSeq_id(0).GetSeqIdString()
	  << " " << (*seqAlign_it)->GetSeqStart(0)
	  << "-" << (*seqAlign_it)->GetSeqStop(0) << endl;
	*/
#ifdef VERBOSE
	cout << "start1-stop1 (target): id "
	     << (*seqAlign_it)->GetSeq_id(1).GetSeqIdString()
	     << " " << (*seqAlign_it)->GetSeqStart(1)
	     << "-" << (*seqAlign_it)->GetSeqStop(1) << endl;
#endif
	
	
	int start1 = (*seqAlign_it)->GetSeqStart(1);
	int stop1 = (*seqAlign_it)->GetSeqStop(1);
	int id1 = 0;
	/*
	sscanf((*seqAlign_it)->GetSeq_id(1).GetSeqIdString().c_str(),
	       "%d", &id1);
	*/
	vector <HitExpansion> expansions = HL.expand_hits_packed(id1, start1, stop1);
	/*
	for (int i = 0; i < vp.size(); i++)
	  cout << vp[i].first << endl << vp[i].second << endl;
	*/
	all_expansions.insert(all_expansions.end(), expansions.begin(), expansions.end());
      }
#ifdef VERY_VERBOSE
      out << MSerial_AsnText << *sas;
#endif
      /*
      if (all_expansions.size())
	cout << "number of expansions for query " << i << ": "
	     << all_expansions.size() << endl;
      */

      string user_input;
      for (int k = 0; k < all_expansions.size(); k++)
	user_input += all_expansions[k].header_plus_dna + '\n';
      CBlastFastaInputSource fasta_input_from_str(user_input, iconfig);
      CBlastInput blast_input_from_str(&fasta_input_from_str);
      
      
      CRef<CObjectManager> objmgr2 = CObjectManager::GetInstance();
      if (!objmgr2) {
	throw std::runtime_error("Could not initialize object manager2");
      }
      CScope scope2(*objmgr2);
      
      TSeqLocVector target_from_str = blast_input_from_str.GetAllSeqLocs(scope2);
      
      //cerr << "beginning re-blast... ";
      
      // blast just the current (i-th) query against the hit expansions
      CBl2Seq blaster_from_str(query_loc[i], target_from_str, *opts);
      TSeqAlignVector results_from_str(blaster_from_str.Run());
      
      //cerr << "finished re-blast" << endl;
      
      set < pair <string, pair <int, int> > > cablast_hits;
      
      // one result (seq align set) for each target sequence
      for (unsigned int j = 0; j < results_from_str.size(); j++) {
#ifdef VERBOSE
	cerr << "results_from_str[" << j << "]" << endl;
#endif
	CConstRef<CSeq_align_set> sas = results_from_str[j];
#ifdef VERY_VERBOSE
	out << MSerial_AsnText << *sas;
#endif
	const list < CRef <CSeq_align> > &seqAlignList = sas->Get();
	ITERATE(list < CRef <CSeq_align> >, seqAlign_it, seqAlignList) {
	  int expansions_ind;
	  sscanf((*seqAlign_it)->GetSeq_id(1).GetSeqIdString().c_str(),
		 "%d", &expansions_ind);
	  const string &target_fasta =
	    all_expansions[expansions_ind-1].header_plus_dna;
	  int second_pipe = target_fasta.find('|');
	  second_pipe = target_fasta.find('|', second_pipe+1);
	  string fasta_prefix =
	    target_fasta.substr(1, min(min(target_fasta.find('\n')-1,
					   target_fasta.find(' ')-1),
				       (size_t) second_pipe-1));
	  int start = (*seqAlign_it)->GetSeqStart(1)
	    + all_expansions[expansions_ind-1].offset;
	  int end = (*seqAlign_it)->GetSeqStop(1)
	    + all_expansions[expansions_ind-1].offset;
	  cablast_hits.insert(make_pair(fasta_prefix,
					make_pair(start, end)));
#ifdef VERBOSE
	  cerr << "hit from " << fasta_prefix
	       << " at pos " << start << "-" << end << endl;
#endif
	}
      }
#ifdef COMPARE_RESULTS
      // compare to full_results
      const list < CRef <CSeq_align> > &full_seqAlignList =
	full_results[i].GetSeqAlign()->Get();
      bool missed_at_least_one = false;
      
      ITERATE(list < CRef <CSeq_align> >, seqAlign_it, full_seqAlignList) {
	full_ctr++;
	string label;
	//(*seqAlign_it)->GetSeq_id(1).GetLabel(&label);
	/*
	string fasta_prefix = (*seqAlign_it)->GetSeq_id(1).AsFastaString();
	int start = (*seqAlign_it)->GetSeqStart(1);
	int end = (*seqAlign_it)->GetSeqStop(1);
	cout << "full blast: hit from " << fasta_prefix
	     << " at pos " << start << "-" << end << " ";
	if (cablast_hits.count(make_pair(fasta_prefix, make_pair(start, end))))
	  cout << "FOUND! ";
	else
	  cout << "Missed. ";
	double evalue;
	(*seqAlign_it)->GetNamedScore(CSeq_align::eScore_EValue, evalue);
	cout << evalue << endl;
	*/
	string fasta_prefix = (*seqAlign_it)->GetSeq_id(1).AsFastaString();
	int start = (*seqAlign_it)->GetSeqStart(1);
	int end = (*seqAlign_it)->GetSeqStop(1);

	if (fasta_prefix.substr(0, 4) == "lcl|") fasta_prefix = fasta_prefix.substr(4);

	if (!cablast_hits.count(make_pair(fasta_prefix,
					  make_pair(start, end)))) {
	  bool found_imperfect = false;

	  // look for hits that aren't quite right
	  for (set < pair <string, pair <int, int> > >::iterator it =
		 cablast_hits.begin(); it != cablast_hits.end(); it++) {
	    if (it->first == fasta_prefix) {
	      int hit_start = it->second.first;
	      int hit_end = it->second.second;
	      if (hit_start <= end && hit_end >= start) {
		found_imperfect = true;
#ifdef PRINT_MISSES
		cout << "found imperfect: " << hit_start << " " << hit_end << endl;
#endif	  
		break;
	      }
	    }
	  }
	  
	  if (!found_imperfect) {
#ifdef PRINT_MISSES
	  cout << "full blast: hit from " << fasta_prefix
	       << " at pos " << start << "-" << end << " ";
	  cout << "Missed. ";
	  double evalue;
	  (*seqAlign_it)->GetNamedScore(CSeq_align::eScore_EValue, evalue);
	  cout << evalue << endl;
#endif	  
	    missed_ctr++;
	    missed_at_least_one = true;
	  }
	}
      }
#ifdef PRINT_MISSES
      if (missed_at_least_one) cout << "above from query " << i << endl;
#endif
#endif
    }
    cout << "Fine BLAST time: " << timer.update_time() << endl;
    
    cout << "Overall accuracy: " << 1-missed_ctr/(double) full_ctr << endl;
    cout << "Missed " << missed_ctr << "/" << full_ctr << endl;

    if (punted) {
      cout << "WARNING: ignored one or more queries that created an anomalously large number of coarse BLAST hits (which would skew results)" << endl;
    }

    /*
    for (unsigned int i = 0; i < results_from_str.size(); i++) {
         CConstRef<CSeq_align_set> sas = results_from_str[i];
         out << MSerial_AsnText << *sas;
    }
    */
    return 0;
}


/////////////////////////////////////////////////////////////////////////////
//  Cleanup


void CCablastApplication::Exit(void)
{
    SetDiagStream(0);
}


/////////////////////////////////////////////////////////////////////////////
//  MAIN


#ifndef SKIP_DOXYGEN_PROCESSING
int main(int argc, const char* argv[])
{
    // Execute main application function
    return CCablastApplication().AppMain(argc, argv, 0, eDS_Default, 0);
}
#endif /* SKIP_DOXYGEN_PROCESSING */
