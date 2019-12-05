#!/usr/bin/bash

#
# This script computes the alignment using TMscore
# bash compute_alignment.sh <multi-fasta file> <template_list> <Pfam family name> 
#


fasta=$1
template=$2
pfam=$3
thread=$4


cd out/aln
cp ../../*.pdb .

# Align MSA (in .aln format)
t_coffee -seq $fasta -template_file $template -method sap_pair TMalign_pair -thread $thread
a=${fasta##*/}
mv ${a%.fa}.aln ${pfam}.aln

# Change sequence names (to shorter ones) 
t_coffee -other_pg seq_reformat -in ${pfam}.aln -output code_name > table
t_coffee -other_pg seq_reformat -in ${pfam}.aln -code table > ${pfam}_coded.aln

# Reformat to .fa & .ph
t_coffee -other_pg seq_reformat -in ${pfam}.aln -output fasta_aln -out ${pfam}.fa 
t_coffee -other_pg seq_reformat -in ${pfam}_coded.aln -output fasta_aln -out ${pfam}_coded.fa 
t_coffee -other_pg seq_reformat -in ${pfam}_coded.aln -output phylip_aln -out ${pfam}_coded.ph

rm *.pdb

 
