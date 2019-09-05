#!/usr/bin/bash

#
# This script computes the alignment using TMscore
# bash compute_alignment.sh <multi-fasta file> <template_list> <Pfam family name> 
#


fasta=$1
template=$2
pfam=$3


cd out/aln
cp ../../*.pdb .
t_coffee -seq $fasta -template_file $template -method sap_pair TMalign_pair
t_coffee -other_pg seq_reformat -in *.aln -output phylip_aln -out ${pfam}_tmalign.ph
t_coffee -other_pg seq_reformat -in *.aln -output fasta_aln -out ${pfam}_tmalign.fa
rm *.pdb

 
