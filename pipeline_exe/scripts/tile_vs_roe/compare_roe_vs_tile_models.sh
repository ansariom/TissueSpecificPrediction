#!/bin/bash

roe_main_indir=$1
tile_main_indir=$2
outdir=$3


echo $outdir
roe_feature_map=$roe_main_indir/features_map.txt
foldchange=3
model_type=l1_logreg
roe_coef_basedir=$roe_main_indir/med_high/$foldchange/$model_type
tile_coef_basedir=$tile_main_indir/med_high/$foldchange/$model_type

if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

roe_coords_out=$outdir/roe_win_coords.txt
tiles_dir=$outdir/tiles
roes_dir=$outdir/roes

if [ ! -d $tiles_dir ]; then
	mkdir -p $tiles_dir
fi

if [ ! -d $roes_dir ]; then
	mkdir -p $roes_dir
fi

# pick one tss name from features_map file to narrow down the search 
#ex_tss="AT1G01010.1_Chr1_650_+_0"
cat $roe_feature_map | head -n 30000 | grep "AT1G01010.1_Chr1_3650_+_0" | cut -f2-4 > $roe_coords_out

w=20
from=-200
to=40
for i in `seq 1 12`; do
	s=$(((i-1) * w + from))
	e=$((s+20))
	echo -e TAContent$i'\t'$s'\t'$e >> $roe_coords_out
done

echo -e GCContent'\t-100\t0' >> $roe_coords_out
echo -e CAContent'\t-100\t0' >> $roe_coords_out
echo -e OC_P_OVERALL_LEAF'\t-100\t0' >> $roe_coords_out
echo -e OC_P_OVERALL_ROOT'\t-100\t0' >> $roe_coords_out

#grep $ex_tss $roe_feature_map | cut -f2-4 > $roe_coords_out
cp -r --backup=t  $tile_coef_basedir/*/tile*coef* $tiles_dir
cp -r --backup=t $roe_coef_basedir/*/roe*coef* $roes_dir

