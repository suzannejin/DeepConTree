#!/usr/bin/bash

#
# This script computes the structure based tree
# It calls a first script in order to extract the distance matrix,
# then it calls the second script to determine the 3d-tree
# and a third script to compute the boostrap values
#
# bash compute_3dtree.sh <Pfam family name Eg. PF00127> <template_list> <MSA method> <T-Coffee MODE> <T-Coffee MODE $EXP (2=square; 3=cubic)>


pfam=$1
template=$2
msa=$3
mode=$4
mexp=$5

dir=$(dirname "$0")
script1=${dir}/RUN_phylo3d_unweighted_d4.pl
script2=${dir}/script_extr_fastme_per_family_D4_unweighted.sh
script3=${dir}/bootrees.R



cd out/tree3d
echo $pfam > list
cp ../aln/${pfam}_tmalign.fa ${pfam}/.
cp $template ${pfam}/${pfam}_ref.template_list2
cp ../../*.pdb ${pfam}/.


# Extract distance matrix
outbash1=${pfam}_${msa}_phylo3D_unw_d4-${mode}-${mexp}.sh
perl $script1 list $msa $mode $mexp  # It creates $outbash script
bash $outbash1


# Compute 3D tree
outbash2=${pfam}_extract_mat_run_FASTME_boottrees_D4_unw.sh
bash $script2
bash $outbash2

# Compute bootstrap values
Rscript $script3 list unw_d4


tree=${pfam}/${pfam}_unw_d4_${msa}${mode}-${mexp}.mat_1.txt.nwk
cp $tree .
rm list
cd ${pfam}
rm *_tmalign.fa *_ref.template_list2 *.pdb 
 
  
