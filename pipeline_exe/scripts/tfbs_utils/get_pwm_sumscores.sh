#!/bin/bash

pwm_name=$1
scores_dir=$2
tmpdir=$3
tmpout_dir=$4
outdir=$5

for strand in 'FWD' 'REV'; do
        for score_file in `ls $scores_dir/*.$strand.cumscores`; do
#		echo $score_file
                sbase=`basename $score_file`
                tmp_file=$tmpdir/$pwm_name.$sbase.tmp
                grep $pwm_name $score_file > $tmp_file
                outbase_scoreext=`mktemp --tmpdir=$tmpdir`
                scripts/tfbs_utils/flip_row_rm_header.pl $tmp_file > $outbase_scoreext

                outbase_locext=`mktemp --tmpdir=$tmpdir`
                filebase=`echo $sbase | cut -d '.' -f1-2`
                locs_file=$scores_dir/$filebase.$strand.locs
                grep $pwm_name $locs_file > $tmp_file
                scripts/tfbs_utils//flip_row_rm_header.pl $tmp_file > $outbase_locext

                paste $outbase_locext $outbase_scoreext > $tmpout_dir/$strand.$pwm_name.$filebase
                rm -f $tmp_file $outbase_scoreext $outbase_locext
        done
        # start gathering PWM info
        cat $tmpout_dir/$strand.$pwm_name.* > $tmpout_dir/$pwm_name.$strand.all.txt
        rm -f $tmpout_dir/$strand.$pwm_name.*

        # call R to sum up all values in the relative location
        outfile=$outdir/$strand/$pwm_name.dist
        scripts/tfbs_utils/sum_dists.R $tmpout_dir/$pwm_name.$strand.all.txt $outfile
        rm -f $tmpout_dir/$pwm_name.$strand.all.txt
        echo $pwm_name " - " $strand
done # end of strand
