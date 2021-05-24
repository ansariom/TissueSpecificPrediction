#!/bin/bash


if [ $# -lt 6 ]; then
	echo "./compute_oc_overlap.sh [feature_map] [leaf_oc_peaks.bed] [outbasename] [outfile] [outbasedir]  [nlines_perfile] [ncpu] "
	exit 
fi

outbasename=$3
outdir=`readlink -f $5`/oc_$outbasename
if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

feature_map=`readlink -f $1`
leaf_oc=`readlink -f $2`
outbase=`readlink -f $outdir/$outbasename`

outfile=$4
ncpu=$7
#leaf_overlap_out=features_oc_long_leaf.txt

## Split input features and run in parallel
nfeatures=$6
fprefix="$outbase".

echo "split -a 5 -d -l $nfeatures $feature_map $fprefix"
split -a 5 -d -l $nfeatures $feature_map $fprefix

echo $fprefix
# Remove header from first file
first_file="$fprefix"00000
tail -n +2 "$first_file" > "$first_file.tmp" && mv "$first_file.tmp" "$first_file"

echo "Start computing OC overlaps ... "
count=0
total=0
for feature_file in `ls $fprefix*`; do
	if [ $count -gt $ncpu ]; then
		echo "$total files compeleted!"
		wait
		let "total=total+count"
		count=0
	else
		let "count=count+1"
	fi

	feature_bed="$feature_file".bed
	awk '{print $8"\t"$9"\t"$10"\t"$1"%"$2"\t.\t+"}' $feature_file > $feature_bed

	#echo $feature_bed
	scripts/oc_utils/compute_intersect_beds.sh $feature_bed $leaf_oc $outbasename &	
done

wait
cat $outbase.*.overlap > $outfile
echo "OC overlaps computed!"

echo "Cleaning up!"
rm -f $outbase.*.overlap
rm -f $outbase.*.bed
rm -f $fprefix*


