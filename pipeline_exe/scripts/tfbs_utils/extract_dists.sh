#!/bin/bash

in_pwms=$1
scores_dir=$2 # cumulitive scores that are gathered from parallel exe
outdir=$3

fwd_outdir=$outdir/FWD
rev_outdir=$outdir/REV

if [ ! -d $fwd_outdir ]; then
	mkdir -p $fwd_outdir
fi

if [ ! -d $rev_outdir ]; then
        mkdir -p $rev_outdir
fi

# get the pwm labels into array
all_pwms_file=pwm_labels.tmp.txt
grep '>' $in_pwms  | cut -c 3- > $all_pwms_file
mapfile -t pwm_names < $all_pwms_file

tmpdir=`mktemp -d --tmpdir=.`
tmpout_dir=`mktemp -d --tmpdir=.`

count=1
ncpu=40

for pwm_name in ${pwm_names[@]}; do
	echo $pwm_name
	if [ $count -lt $ncpu ]; then
#		echo software/get_pwm_sumscores.sh $pwm_name $scores_dir $tmpdir $tmpout_dir $outdir 
		scripts/tfbs_utils/get_pwm_sumscores.sh $pwm_name $scores_dir $tmpdir $tmpout_dir $outdir &
		let "count=count+1"
	else
		wait
		scripts/tfbs_utils/get_pwm_sumscores.sh $pwm_name $scores_dir $tmpdir $tmpout_dir $outdir &
		count=1
	fi
done
wait

rm -f -r $tmpdir $tmpout_dir
