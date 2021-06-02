#!/bin/bash

#this script collects all PWM matrices from HUghes data and make a single PWM file for pipeline use

ath_indir=/dfs/Megraw_Lab/data/public_db/HughesPWM/Arabidopsis_thaliana_2016_08_30_5-11_pm/pwms_all_motifs/
human_indir=/dfs/Megraw_Lab/data/public_db/HughesPWM/Homo_sapiens_2021_05_21_3_48_pm/pwms_all_motifs/

indir=$human_indir
outfile=hughes_pwms.mat

rm -f $outfile

for f in `ls $indir`; do
	pwm=$indir/$f
	l=`wc -l $pwm | cut -d " " -f1`
	echo $l
	if [ $l -gt 4 ]; then
		mname=$(echo $f | cut -f 1-2 -d '.')
		awk -v mname=$mname '{if (NR == 1){print "> "mname} else {print "=\t"$2"\t"$3"\t"$4"\t"$5}}' $pwm >> $outfile
		echo  >> $outfile
	fi
done

