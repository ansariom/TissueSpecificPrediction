#!/bin/bash


if [ $# -lt 4 ]; then
	echo "./compute_oc_overlap.sh [ranges.bed] [oc_peaks.bed] [base] [outfile] "
	exit 
fi

base=$3
range_file=`readlink -f $1`
oc_file=`readlink -f $2`
outfile=$4

feature_bed=$range_file.$base.bed
awk '{print $2"\t"$3"\t"$4"\t"$1"\t.\t+"}' $range_file > $feature_bed

#echo $feature_bed
tmp_overlap="$feature_bed".overlap.tmp
bedtools intersect -wao -a $feature_bed  -b $oc_file  | awk '{a[$1"\t"$2"\t"$3"\t"$4]+=$12}END{for (i in a){print i,a[i]}}' > $tmp_overlap

awk '{print($0"\t"($5/($3-$2)))}' $tmp_overlap | awk '{$1=$1}1' OFS="\t" | cut -f4,6 | awk -v var="$base" '{print $1"\tOC_P_OVERALL_"var"\t"$2}'  > $outfile 

rm -f $tmp_overlap
rm -f $feature_bed

