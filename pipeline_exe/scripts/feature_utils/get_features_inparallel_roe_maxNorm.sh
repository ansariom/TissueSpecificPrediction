#!/bin/bash

ncpu=50

input_fasta=$1
outdir=$2
final_outfile=$3
final_mapfile=$4
roe_fwd=$5
roe_rev=$6
psudocount=$7
pwm=$8
maxScoresFile=$9

basename=seq
seqs_outdir=$outdir/seqs_features

if [ ! -d $seqs_outdir ]; then
	mkdir -p $seqs_outdir
fi

scripts/seq_utils/split_fasta_file.sh $input_fasta 200 $seqs_outdir $basename

outdir_tmp=$outdir/features
rm -f -r $outdir_tmp

if [ ! -d $outdir_tmp ]; then
	mkdir -p $outdir_tmp
fi

outdir_maps_tmp=$outdir/maps
if [ ! -d $outdir_maps_tmp ]; then
        mkdir -p $outdir_maps_tmp
fi

count=1
for seqfile in `ls $seqs_outdir/$basename.*`; do
	echo $seqfile $outfile $maps_out
	outfile=$outdir_tmp/`basename $seqfile`.Rdat
        maps_out=$outdir_maps_tmp/`basename $seqfile`.maps.txt
	#java -Xms10G -Xmx20G -jar software/tfbs_scan.jar GenFeatures --pseudoCounts $psudocount --mapFile $maps_out -n 5 $roe_fwd $roe_rev $seqfile $pwm $outfile &
	java -Xms10G -Xmx20G -jar external_tools/tfbs_scan.jar GenFeaturesNormByMax --pseudoCounts $psudocount --mapFile $maps_out -n 5 $roe_fwd $roe_rev $seqfile $pwm $outfile --maxScores $maxScoresFile &

        if [ $count -lt $ncpu ]; then
		let "count=count+1"
	else
		echo "wait for batch to complete! (submitted = $count)"
                wait
                count=1
	fi
done
wait

echo "All finished!!!"

rm -f $final_outfile

cat $outdir_tmp/*.Rdat > $outdir_tmp/all.Rdat.tmp
head -n 1 $outdir_tmp/all.Rdat.tmp > $final_outfile
grep -v "_FWD_" $outdir_tmp/all.Rdat.tmp >> $final_outfile

cat $outdir_maps_tmp/*.maps.txt > $outdir_maps_tmp/all_maps.tmp
head -n 1 $outdir_maps_tmp/all_maps.tmp > $final_mapfile
grep -v "seq_id" $outdir_maps_tmp/all_maps.tmp >> $final_mapfile

# Cleaning up
