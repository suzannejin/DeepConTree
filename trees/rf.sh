
#
# This script compares each tree to each other and calculates a normalized RF score for each pair of trees
# and writes an output file with a pairwise matrix of RF score
#

pfam=$1

msa=$2
mode=$3
mexp=$4

dir=$(dirname "$0")
script=${dir}/rf.R
env=${dir}/../set_envs.sh


cd trees

## Create a text file with two columns: 1) path to tree.nwk  and  2) tree name
#touch trees.txt
#echo "experimental/out/tree3d/${pfam}_unw_d4_${msa}${mode}-${mexp}.mat_1.txt.nwk experimental" >> trees.txt
#trees=*/*/out/tree3d/${pfam}_unw_d4_${msa}${mode}-${mexp}.mat_1.txt.nwk
#for tree in $trees; do
#    first=${tree%%/*}
#    second=$(echo $tree | awk '{split($0,a,"/");print(a[2])}')
#    name=${first}_${second}
#    echo "$tree $name" >> trees.txt
#done



# Merge tree.nwk files
[[ -s trees.nwk ]] && rm trees.nwk
touch trees.nwk
files=$(cut -f1 -d" " tree_list)
for file in $files; do
    cat $file >> trees.nwk
done

# Calculate RF score
source $env tree && Rscript $script tree_list trees.nwk rf.txt rf.xlsx

