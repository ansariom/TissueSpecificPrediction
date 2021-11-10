#!/bin/bash

# output is in paper_util directory
outdir=.
roe_main_idir=../pipeline_output_roe/
tile_main_indir=../pipeline_output_tile/
src_dir=../scripts/

rm -f -r logs
mkdir logs

if [ ! -d $outdir ]; then
	mkdir $outdir
fi

leaf_oc_promoter=$outdir/oc/leaf_promoter_open_regions_3000-3000.txt
root_oc_promoter=$outdir/oc/root_promoter_open_regions_3000-3000.txt
diff_expr_genes=../external_tools/RNASeq_run1199/rsem_output/ath_root_leaf_rsem_deseq_diff_expr_results_filtered.txt
tss_peaks=$roe_main_dir/aligned.peaks.annotated.capped.filtered

function common {
        echo $src_dir/utils_paper/get_oc_regions_for_promoters.sh 3000 3000 $roe_main_idir/peaks_3000_region.fa $outdir/oc $src_dir $roe_main_idir/oc_peaks_root.bed $roe_main_idir/oc_peaks_leaf.bed | SGE_Array -m 100G -q *@megraw* -r logs/logs_oc_promoter_dist
}

function fig2A {
	echo $src_dir/utils_paper/oc_plots_for_region.R $outdir/oc/ $outdir/oc/  | SGE_Array -m 100G -q *@megraw* -r logs/logs_oc_plots_for_region --hold_names logs/logs_oc_promoter_dist
}

function fig2B {
	echo echo $src_dir/basic_data_stats.R $leaf_oc_promoter $root_oc_promoter $tss_peaks $diff_expr_genes | SGE_Array -m 100G -q *@megraw* -r logs/log_basic_stats_tss_oc --hold_names logs/logs_oc_promoter_dist
}
function suppFig9{
	echo ../scripts/performance/gather_model_performances.sh $roe_main_idir $tile_main_indir $src_dir $outdir/model_performance | SGE_Array -m 100G -q *@megraw* -r logs/logs_gather_model_performances_plots
}
function suppFig10_11 {
	$src_dir/tile_vs_roe/compare_roe_vs_tile_models.sh $roe_main_idir $tile_main_indir tile_vs_roe
	$src_dir/utils_paper/model_coef_variability.R tile_vs_roe $roe_main_idir/roe_only_out/med_high/3/l1_logreg/810376/featureInfo_hardCodedSoftCoded.rdat $tile_main_indir/tile_only_out/med_high/3/l1_logreg/67086/featureInfo_hardCodedSoftCoded.rdat
}
#common
#fig2A
#fig2B
#suppFig9


# compare ROE vs. Tile features
echo "Comapre ROE vs. Tiled feature overlaps .."
#../scripts/tile_vs_roe/compare_roe_vs_tile_models.sh $roe_main_idir $tile_main_indir $outdir/tile_vs_roe/
# run "" in ../scripts/tile_vs_roe/heatmap_top_features_tile_roe_multiRuns.R"
echo "done! output in $outdir/tile_vs_roe/"

echo ""
#echo $src_dir/utils_paper/get_oc_regions_for_promoters.sh 3000 3000 $roe_main_idir/peaks_3000_region.fa $outdir/oc $src_dir $roe_main_idir/oc_peaks_root.bed $roe_main_idir/oc_peaks_leaf.bed | SGE_Array -m 100G -q *@megraw* -r logs/logs_oc_tss_stats_plots 

