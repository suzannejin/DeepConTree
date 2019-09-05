
prediction=$1
experimental=$2
out=$3

awk '$1=="HEADER"||$1=="COMPND"||$1=="REMARK"||$1=="SEQRES"{print $0}' $experimental > $out
chain=$(awk '$1=="COMPND"{print $4}' $experimental | sed 's/;//g')
awk -v chain=$chain '$1=="ATOM"{printf"%s%+7s  %-4s%s %s%+4s%+12s%+8s%+8s  %s\n",$1,$2,$3,$4,chain,$6,$7,$8,$9,substr($0,57,9)}' $prediction >> $out
