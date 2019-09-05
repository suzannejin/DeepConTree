
#
# This script calculates the Cscore and filters the contacts
# Note: better give the absolute path
# 


fa=$1   # Alignment in fasta format, with PDB ID as sequence name
con=$2  # Contactbench file. It should be saved in the same folder as the contactbench files of the remaining proteins
out=$3  # Output prefix
p=$4    # Pscore threshold
c=$5    # Cscore threshold
w=$6    # Consider "yes/no" the Voronoid weight




dir=$(dirname "$0")
script=${dir}/msa2cscore.py


python $script $fa $con ${out}.cscore $p 0 $w
awk -v p=$p -v c=$c 'NR!=1&&$5>=c{printf "%s %s 0 8 %s\n",$2,$3,$4}' ${out}.cscore | sort -k5gr > ${out}.con


