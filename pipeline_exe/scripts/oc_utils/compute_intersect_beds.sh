#!/bin/bash

file_a=$1
file_b=$2
basename=$3

tmp_overlap="$file_a".overlap.tmp

bedtools intersect -wao -a $file_a -b $file_b  | awk '{a[$1"\t"$2"\t"$3"\t"$4]+=$12}END{for (i in a){print i,a[i]}}' > $tmp_overlap

overlap_out="$file_a".overlap
awk '{print($0"\t"($5/($3-$2)))}' $tmp_overlap | awk '{$1=$1}1' OFS="\t" | cut -f4,6 | awk -v var="$basename" '{split($1,a,"%"); print a[1]"\t"a[2]"_"var"\t"$2}'  > $overlap_out 

#rm -f $tmp_overlap
