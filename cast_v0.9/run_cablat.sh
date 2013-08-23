#!/bin/bash

if [ $# -ne 5 ]
then
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Script to run CaBLAT, which consists of two steps:

- Coarse search, done by running the blat executable directly
- Fine search, done by running the cablat executable
  (passing it the coarse results file from the previous step).

Arguments:

1 - compressed db file: uniques (fasta)
2 - compressed db file: links (dat)
3 - compressed db file: edit scripts (dat)
4 - query file (fasta)
5 - output file for writing results (psl)
? - (maybe in a later version) additional args such as minScore?

Example usage: ./run_cablat.sh uniques.fasta links.dat edit_scripts.dat query.fasta output.psl

Make sure this script contains the correct paths to the blat and
cablat executables!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    exit
fi

if [ -z "$MACHTYPE" ]
then
    echo "error: MACHTYPE is undefined or empty"
    exit 1
fi

coarse_out_file="coarse_output.psl"
blat_exec_path="$HOME/bin/$MACHTYPE/blat" # replace with path to blat executable
cablat_exec_path="$HOME/cablast/BLAT/blatSrc/cablat/cablat" # replace with path to cablat executable

echo "coarse_out_file: $coarse_out_file"
echo "blat_exec_path: $blat_exec_path"
echo "cablat_exec_path: $cablat_exec_path"

########## run blat on uniques: produce psl output ##########

time $blat_exec_path $1 $4 $coarse_out_file -noHead -minIdentity=80


########## run cablat: read psl output from uniques; do fine search; make final output  ##########

time $cablat_exec_path $coarse_out_file $1 $2 $3 $4 $5
