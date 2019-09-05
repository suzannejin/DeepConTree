#!/usr/bin/bash
################################################################################
#This script generates a script called $family'_extract_mat_run_FASTME_boottrees_D4_unw.sh'
# that divides the 101 matrices that are in the matrices file in a folder called
#'boot3d-d4-unw/' then it generates a tree for each matrix.
# It then copies the matrix and the tree in the family folder for easy access!
################################################################################


main_dir=$1
pfam=$2

outbash=${main_dir}/out/tree3d/${pfam}_extract_mat_run_FASTME_boottrees_D4_unw.sh

echo cd ${main_dir}/out/tree3d > $outbash ;
echo mkdir boot3d-d4-unw/ >> $outbash ;
echo cp ${pfam}.matrices boot3d-d4-unw/. >> $outbash ;
#echo 'cp ../fastme boot3d-d4-unw/.' >> $outbash ;
echo cd boot3d-d4-unw/ >> $outbash ;
echo "awk -v RS= '{print > (\""$pfam"_unw_d4_tmalign4-2.mat_"\" NR "\".txt\")}' ${pfam}.matrices"  >> $outbash ;
echo 'for j in `ls|grep '.txt'`; do fastme -i $j -o ${j}.nwk -g 1.0 -s -n -z 5; done' >> $outbash ;
echo cd ../../../  >> $outbash ;
echo mv out/tree3d/boot3d-d4-unw/${pfam}_unw_d4_tmalign4-2.mat_1.txt out/tree3d/. >> $outbash ;
echo mv out/tree3d/boot3d-d4-unw/${pfam}_unw_d4_tmalign4-2.mat_1.txt_fastme_stat.txt out/tree3d/. >> $outbash ;
echo mv out/tree3d/boot3d-d4-unw/${pfam}_unw_d4_tmalign4-2.mat_1.txt.nwk out/tree3d/.  >> $outbash ;
echo 'for k in `ls '$main_dir'/out/tree3d/boot3d-d4-unw/*txt.nwk`; do echo $k ; cat $k; done > '$main_dir'/out/tree3d/boot3d-d4-unw/'$pfam'_bootrees' >> $outbash

