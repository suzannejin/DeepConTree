

pdb1=$1
pdb2=$2

command="TMscore $pdb1 $pdb2"
rmsd=$($command | grep -E "^RMSD" | sed 's/   */ /g' | cut -f6 -d" ")
tmscore=$($command | grep -E "^TM-score" | sed 's/  */ /g' | cut -f3 -d" ")
echo "$rmsd $tmscore"
