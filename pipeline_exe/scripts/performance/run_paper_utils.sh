#!/bin/bash

# output is in paper_util directory
outdir=/nfs0/BPP/Megraw_Lab/mitra/Projects/IDBC/github/IBDC/exe/supplemetal_scripts/aug2019/
roe_main_idir=../ibdc_roe-only_aug2019/
tile_main_indir=../ibdc_tile-only_aug2019/
src_dir=../../src/

rm -f -r logs
mkdir logs

if [ ! -d $outdir ]; then
	mkdir $outdir
fi

# model stability analsysis + roe_tile performance comparision
echo "gather all 30 model performances for ROE and Tile and plot them .."
#../../src//performance/gather_model_performances.sh $roe_main_idir $tile_main_indir $src_dir $outdir
echo "done! output in: $outdir"

# compare ROE vs. Tile features
echo "Comapre ROE vs. Tiled feature overlaps .."
#../../src/tile_vs_roe/compare_roe_vs_tile_models.sh $roe_main_idir $tile_main_indir $outdir/tile_vs_roe/
echo "done! output in $outdir/tile_vs_roe/"

echo ""
echo ../../src/utils_paper/get_oc_regions_for_promoters.sh 3000 3000 $roe_main_idir/peaks_3000_region.fa $outdir/oc $src_dir $roe_main_idir/oc_peaks_root.bed $roe_main_idir/oc_peaks_leaf.bed | SGE_Array -m 100G -q *@megraw* -r logs/logs_oc_tss_stats_plots 

