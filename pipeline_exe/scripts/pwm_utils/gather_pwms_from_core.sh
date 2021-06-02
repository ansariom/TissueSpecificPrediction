#!/bin/bash

#this script collects all PWM matrices from HUghes data and make a single PWM file for pipeline use


indir=/dfs/Megraw_Lab/data/external/pwm_databases/lab_generated/pwm/
outfile=core_pwms.mat

rm -f $outfile

while read -r line; do
	pwm=$indir/$line".mat"
	l=`wc -l $pwm | cut -d " " -f1`
	echo $l
	if [ $l -gt 4 ]; then
		#mname=$(echo $f | cut -f 1-2 -d '.')
		mname=$line
		awk -v mname=$mname '{if (NR == 1){print "> "mname} else {print "=\t"$2"\t"$3"\t"$4"\t"$5}}' $pwm >> $outfile
		echo  >> $outfile
	fi
done < other_pwms_list.txt

