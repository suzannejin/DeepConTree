
#
# This scripts uses Spider3 to predict secondary structure (SSE)
#



spiderf=$1  # Path to Spider3
fasta=$2  # Fasta file
id=$3  # Protein ID
prefix=$4  # Directory + prot ID


dir=$(dirname "$0")
env=${dir}/../set_envs.sh


cd $spiderf

cp $fasta $id
cp ${prefix}.pssm ${id}.pssm
cp ${prefix}.hhm ${id}.hhm


echo ${id} ${id}.pssm ${id}.hhm > tlist
source $env spider && bash ./scripts/impute_script.sh tlist


mv ${id}.spd33 ${prefix}.spd33
rm ${id}* tlist
