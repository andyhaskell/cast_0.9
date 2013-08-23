***************************************************************
***** Compression-accelerated BLAST and BLAT, version 0.9 *****
*****    (pre-release for inclusion with publication)     *****
*****                   June 21, 2012                     *****
*****                                                     *****
*****      Po-Ru Loh, Michael Baym and Bonnie Berger      *****
***************************************************************

This software contains implementations of compression-accelerated
BLAST and BLAT described in:

  Loh P-R, Baym M, Berger B.  Compressive genomics.  Nature
    Biotechnology, Volume 30 Number 7, July 2012.

These algorithms serve as proof-of-concept that computationally-aware
compression not only reduces storage space but also accelerates
analysis (in this case, sequence search).  Note that the
implementations contained herein are prototypes, not yet sufficiently
tuned, tested, or interfaced for inclusion in computational pipelines.
We anticipate that to achieve optimal performance, the code will need
to be tailored to match the particular engineering trade-offs that
arise in real-world applications.  The code in its current form is
intended primarily as a resource for developers interested in adapting
it for specific applications (or working with us to build it into
"industrial-strength" software for general use by practitioners).  We
provide it for academic and non-commercial use only.


****************************************
***** Before proceeding further... *****
****************************************

Please check our website:

  http://cast.csail.mit.edu

to make sure you have the latest version!  We expect to be updating
this package with bug fixes, added functionality and improved
documentation as we receive feedback.


********************
***** Contents *****
********************

1. Installing the software
   a. Compiling the preprocessing (compression) phase
   b. Installing CaBLAST
   c. Installing CaBLAT
2. Running the software
   a. Creating a compressed database
   b. Running CaBLAST search
   c. Running CaBLAT search
3. Known issues, limitations, etc.
4. License



~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
|                   1. Installing the software                       |
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


***************************************************************
***** 1a. Compiling the preprocessing (compression) phase *****
***************************************************************

Just run:

  g++ -O2 gen_compressed_db.cpp -o gen_compressed_db

to create the executable 'gen_compressed_db'.


**********************************
***** 1b. Installing CaBLAST *****
**********************************

- Download and build the NCBI C++ Toolkit.  We tested CaBLAST using
  version 7.0.0 of the C++ Toolkit (June 2011 release) and recommend
  compiling CaBLAST using version 7; CaBLAST's accuracy seems to
  decrease when built with version 9.0.0 (May 2012 release) and we
  have not yet had a chance to investigate this behavior.  Version
  7.0.0 of the toolkit is available at:

    ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools++/ARCHIVE/7_0_0/ncbi_cxx--7_0_0.tar.gz

  On our machine, we configured the toolkit with:

    ./configure --with-dll --with-64 --without-debug


- Create a new toolkit project 'cablast'.  This can be done by running
  the new_project script (scripts/common/new_project.sh in the toolkit
  source tree):

    bash <path_to>/new_project.sh cablast app [builddir]

  This will create a subdirectory of the current directory called
  'cablast' with template project files.

  For documentation on project creation with the C++ Toolkit, see
  e.g. Chapter 6 of the NCBI C++ Toolkit Book, "Project Creation and
  Management":

    http://www.ncbi.nlm.nih.gov/books/NBK7171/
    -> new_project: Starting a New Project outside the C++ Toolkit Tree


- Check that the 'Makefile.builddir' file contains the path to the
  build directory of your installation of the NCBI C++ Toolkit.

- Modify the set of additional libraries to link to by adding the
  following lines to 'Makefile.cablast.app':

    LIB_ = $(BLAST_INPUT_LIBS) $(BLAST_LIBS) $(OBJMGR_LIBS)
    LIB = $(LIB_:%=%$(STATIC))

  (These are copied from 'src/app/blast/Makefile.blastn.app', the
  corresponding file for the blastn project.)


- Copy the *.cpp files from this package to the cablast project
  directory.

- Run 'make' within the cablast project directory to create the
  'cablast' executable.


*********************************
***** 1c. Installing CaBLAT *****
*********************************

- Download the BLAT source and follow its installation instructions.
  We tested CaBLAT using version 34 of the BLAT source:

    http://users.soe.ucsc.edu/~kent/src/blatSrc34.zip

  To successfully compile BLAT, you may need to define and export the
  MACHTYPE environment variable:

    MACHTYPE=x86_64   <-- replace with your machine configuration
    export MACHTYPE

  You may also need to disable the -Werror compile flag by removing it
  from line 17 of the 'inc/common.mk' file.


- Create a subdirectory 'cablat' of the 'blatSrc' directory and copy
  in the contents of this package.

- Within the cablat directory, run:

    ./add_extern_C.sh

  to allow C++ linking to the BLAT library (jkOwnLib.a).

- Run 'make' within cablat directory to create the 'cablat'
  executable.



~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
|                    2. Running the software                         |
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


**********************************************
***** 2a. Creating a compressed database *****
**********************************************

The tool 'gen_compressed_db' creates a compressed database from a
FASTA file; for usage information, run it without any input arguments.


**************************************
***** 2b. Running CaBLAST search *****
**************************************

The CaBLAST code currently produces a 'cablast' executable that runs a
performance comparison between BLAST and CaBLAST.  Run 'cablast -help'
for usage information.  Note that 'cablast' compares results but does
not print the actual hits; outputting hits will require modifying the
'cablast.cpp' source accordingly and recompiling.

   
*************************************
***** 2c. Running CaBLAT search *****
*************************************

We provide a script 'run_cablat.sh' that runs the coarse search step
(directly calling the blat executable) followed by fine search (using
our cablat executable on the results of the coarse step).  For usage
information, run 'run_cablat.sh' without any input arguments.  Note
that you will need to modify the paths to the blat and cablat
executables in the script in order for it to work.



~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
|                 3. Known issues, limitations, etc.                 |
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As mentioned in the preface to this readme, the code provided here is
proof-of-concept software and as such has various limitations,
including the following:

- The code that constructs the compressed database (gen_compressed_db)
  currently stores sequence positions and link pointers with 32-bit
  signed integers; as such, it will fail for data sizes greater than
  about 2 billion bases.

- The CaBLAST and CaBLAT fine search steps perform searches against
  potential hit regions individually and from scratch (i.e., without
  using information from the coarse search and link-tracing -- which
  should already have done most of the work).  As such, they are quite
  inefficient: individual searches incur overhead (especially in
  BLAST) and re-alignment takes work.  In extreme cases where a query
  matches to a very large number of hits, CaBLAST/CaBLAT can even end
  up taking longer than BLAST/BLAT.  In particular, our current
  implementations are not suitable for searching for weak hits.

- In the parameter ranges we tested, CaBLAST and CaBLAT produce
  results similar to BLAST and BLAT but not exactly the same, as
  discussed in detail in our paper.  Accuracy will depend on many
  parameters, including in particular the window identity threshold
  during preprocessing and the coarse and fine E-value/minIdentity
  thresholds during search.  (Increasing the preprocessing window
  identity threshold will tend to increase search accuracy at the
  expense of search speed and compression ratio; relaxing the coarse
  E-value/minIdentity search threshold will tend to increase accuracy
  at the expense of speed.)

- The CaBLAST code currently runs BLAST and CaBLAST and compares the
  results without explicitly writing the search results.

- The search results written by our CaBLAT code currently contain some
  duplicates (caused by overlapping link regions); we eliminate these
  during post-processing in our analyses but the duplicates are
  present in output files.

- The CaBLAT code does not currently take coarse and fine minIdentity
  thresholds as input parameters; they are set to default values of 80
  and 90, respectively.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
|                             4. License                             |
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This software is made available for academic and non-commercial use
only.