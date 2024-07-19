#!/bin/bash
input=$1
output=$2
mu_scale=$3
gzip -cd ${input} | awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++){if($i=="X")x=i;if($i=="Y")y=i}print $x,$y;next}{print $x,$y}' | perl -slane 'print join("\t",$F[0]/${mu_scale},$F[1]/${mu_scale})' -- -mu_scale="${mu_scale}" | awk -F'\t' ' BEGIN { min1 = "undef"; max1 = "undef"; min2 = "undef"; max2 = "undef"; } { if (NR == 2 || $1 < min1) min1 = $1; if (NR == 2 || $1 > max1) max1 = $1; if (NR == 2 || $2 < min2) min2 = $2; if (NR == 2 || $2 > max2) max2 = $2; } END { print "xmin\t", min1; print "xmax\t", max1; print "ymin\t", min2; print "ymax\t", max2; }' > ${output}
