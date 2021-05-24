#!/bin/bash

if [ "$#" -lt 4 ]; then
    echo "insufficient arguments: [roe_input_dir] [tile_input_dir] [src_dir] [outdir] "
    exit
fi

#src_dir=../src/
src_dir=$3

#outdir=paper_utils/aug2019/
outdir=$4

#nums=("2" "2.5" "3"  "3.5")
nums=("3")
#nums=( "2.5" "3"  "3.5")
#exp_levels=( "med_low" "med_high")
exp_levels=( "med_high")
model_type=main
#model_type="tfbs_only"
#model_type="oc_only"

#indirs=(ibdc_roe-only ibdc_tile-only)
indirs=($1 $2)
name=(ROE_Model Tiled_Model)


outfile1=${indirs[0]}/$model_type.performance.txt
outfile2=${indirs[1]}/$model_type.performance.txt
rm -f $outfile1
rm -f $outfile2

pattern1=$model_type
pattern2=$model_type
if [ $model_type == "main" ]; then
	pattern1="roe_only"
	pattern2="tile_only"
fi

for e in ${exp_levels[@]}
do
	for i in ${nums[@]}
	do
		head ${indirs[0]}/$e/$i/l1_logreg/[1-9]*/$pattern1*held* | grep -v auroc | grep -v '==' | awk -v l=$e -v f=$i '{print $0"\t"l"\t"f}' >> $outfile1
		head ${indirs[1]}/$e/$i/l1_logreg/[1-9]*/$pattern2*held* | grep -v auroc | grep -v '==' | awk -v l=$e -v f=$i '{print $0"\t"l"\t"f}' >> $outfile2
		echo $i
	done
done
$src_dir/performance/plot_performances_models.R $outfile1 $outfile2 "TEP-ROE" "TEP-Tiled" "TEP" $outdir
