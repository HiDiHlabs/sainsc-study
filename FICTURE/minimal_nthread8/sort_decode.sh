#!/bin/bash
input=$1
output=$2
coor=$3
model_id=$4
bsize=$5
scale=$6
topk=$7
bgzip=$8
tabix=$9

K=$( echo $model_id | sed 's/nF\([0-9]\{1,\}\)\..*/\1/' )
while IFS=$'\t' read -r r_key r_val; do
    export "${r_key}"="${r_val}"
done < ${coor}
echo -e "${xmin}, ${xmax}; ${ymin}, ${ymax}"

offsetx=${xmin}
offsety=${ymin}
rangex=$( echo "(${xmax} - ${xmin} + 0.5)/1+1" | bc )
rangey=$( echo "(${ymax} - ${ymin} + 0.5)/1+1" | bc )
bsize=2000
scale=100
header="##K=${K};TOPK=${topk}\n##BLOCK_SIZE=${bsize};BLOCK_AXIS=X;INDEX_AXIS=Y\n##OFFSET_X=${offsetx};OFFSET_Y=${offsety};SIZE_X=${rangex};SIZE_Y=${rangey};SCALE=${scale}\n#BLOCK\tX\tY\tK1\tK2\tK3\tP1\tP2\tP3"

(echo -e "${header}" && gzip -cd "${input}" | tail -n +2 | perl -slane '$F[0]=int(($F[1]-$offx)/$bsize) * $bsize; $F[1]=int(($F[1]-$offx)*$scale); $F[1]=($F[1]>=0)?$F[1]:0; $F[2]=int(($F[2]-$offy)*$scale); $F[2]=($F[2]>=0)?$F[2]:0; print join("\t", @F);' -- -bsize="${bsize}" -scale="${scale}" -offx="${offsetx}" -offy="${offsety}" | sort -S 1G -k1,1g -k3,3g ) | ${bgzip} -c > ${output}

${tabix} -f -s1 -b3 -e3 ${output}
rm ${input}
