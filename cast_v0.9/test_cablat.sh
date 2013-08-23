#!/bin/bash

if [ $# -ne 8 ]
then
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Script to compare BLAT and CaBLAT results.  You will first need to
compile compare_results:

  g++ compare_results.cpp -o compare_results

Arguments:
1 - compressed db: uniques (fasta)
2 - compressed db: links (dat)
3 - compressed db: edit scripts (dat)
4 - query (fasta)
5 - output (psl)
6 - full db (fasta)
7 - full output (psl)
8 - compare output file
? - additional args such as minScore?

Example usage:
./test_cablat.sh uniques.fasta links.dat edit_scripts.dat query.fasta output.psl all_seqs.fasta output_check.psl
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    exit
fi

blat_exec_path="$HOME/bin/$MACHTYPE/blat" # replace with path to blat executable
cablat_dir="$HOME/cablast/BLAT/blatSrc/cablat" # replace with path to cablat directory
cablat_exec_path="$cablat_dir/cablat" 

coarse_out_file="coarse_output.psl"

# run original BLAT

( time -p ( $blat_exec_path $6 $4 $7 -noHead ) ) 2>&1

# run coarse BLAT

( time -p ( $blat_exec_path $1 $4 $coarse_out_file -noHead -minIdentity=80 ) ) 2>&1

# run fine BLAT (cablat executable)

( time -p ( $cablat_exec_path $coarse_out_file $1 $2 $3 $4 $5 ) ) 2>&1

# compare results

$cablat_dir/compare_results $7 $5 > $8
