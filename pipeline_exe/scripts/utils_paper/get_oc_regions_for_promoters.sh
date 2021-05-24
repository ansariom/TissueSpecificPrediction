#!/bin/bash

up=$1
down=$2
in_fasta=$3
outdir=$4
src_dir=$5

if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

promoter_bed=$outdir/root_promoter.bed
outfile=$outdir/root_promoter_open_regions_$up-$down".txt"

oc_peaks_root=$6
oc_peaks_leaf=$7

$src_dir/seq_utils/get_tss_region.sh $up $down $up $down $in_fasta > $promoter_bed 
$src_dir/oc_utils/get_oc_promoter_intersects.sh $promoter_bed $oc_peaks_root $outfile

promoter_bed=$outdir/leaf_promoter.bed
outfile=$outdir/leaf_promoter_open_regions_$up-$down".txt"

$src_dir/seq_utils/get_tss_region.sh $up $down $up $down $in_fasta > $promoter_bed
$src_dir/oc_utils/get_oc_promoter_intersects.sh $promoter_bed $oc_peaks_leaf $outfile
