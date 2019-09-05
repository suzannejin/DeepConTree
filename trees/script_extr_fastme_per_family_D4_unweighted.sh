#!/usr/bin/bash
################################################################################
#This script generates a script for each family called $family'_extract_mat_run_FASTME_boottrees_D4_unw.sh'
# that divides the 101 matrices that are in the matrices file in a folder called
#'boot3d-d4-unw/' then it generates a tree for each matrix.
# It then copies the matrix and the tree in the family folder for easy access!
################################################################################
for i in `cat list`; do echo cd $i/ > $i'_extract_mat_run_FASTME_boottrees_D4_unw.sh';
 echo 'mkdir boot3d-d4-unw/' >> $i'_extract_mat_run_FASTME_boottrees_D4_unw.sh';
 echo cp $i'_unw_d4_tmalign4-2.matrices' boot3d-d4-unw/. >> $i'_extract_mat_run_FASTME_boottrees_D4_unw.sh' ;
 #echo cp ../fastme boot3d-d4-unw/. >> $i'_extract_mat_run_FASTME_boottrees_D4_unw.sh' ;
 echo cd boot3d-d4-unw/ >> $i'_extract_mat_run_FASTME_boottrees_D4_unw.sh' ;
 echo "awk -v RS= '{print > (\""$i"_unw_d4_tmalign4-2.mat_"\" NR "\".txt\")}'" $i'_unw_d4_tmalign4-2.matrices'  >> $i'_extract_mat_run_FASTME_boottrees_D4_unw.sh';
 echo 'for j in `ls|grep '.txt'`;do fastme -i $j -o $j".nwk" -g 1.0 -s -n -z 5; done' >> $i'_extract_mat_run_FASTME_boottrees_D4_unw.sh';
 echo cd ../../  >> $i'_extract_mat_run_FASTME_boottrees_D4_unw.sh';
 echo mv $i/boot3d-d4-unw/$i'_unw_d4_tmalign4-2.mat_1.txt' $i/. >> $i'_extract_mat_run_FASTME_boottrees_D4_unw.sh' ;
 echo mv $i/boot3d-d4-unw/$i'_unw_d4_tmalign4-2.mat_1.txt_fastme_stat.txt' $i/. >> $i'_extract_mat_run_FASTME_boottrees_D4_unw.sh';
 echo mv $i/boot3d-d4-unw/$i'_unw_d4_tmalign4-2.mat_1.txt.nwk' $i/.  >> $i'_extract_mat_run_FASTME_boottrees_D4_unw.sh';
 echo 'for k in `ls '$i'/boot3d-d4-unw/*txt.nwk`; do echo $k ; cat $k; done > '$i'/boot3d-d4-unw/'$i'_bootrees' >> $i'_extract_mat_run_FASTME_boottrees_D4_unw.sh';
done
