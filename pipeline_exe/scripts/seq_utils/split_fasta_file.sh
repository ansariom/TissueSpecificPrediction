#!/bin/bash

inseqs=$1
nseqs_pfile=$2
outdir=$3
base=$4

rm -f -r $outdir
if [ ! -d $outdir ];then
        mkdir $outdir
fi

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $inseqs | tail -n +2 > tmp.fa

split -a 3 -d -l $(($nseqs_pfile*2)) tmp.fa $outdir/$base.

rm -f tmp.fa
