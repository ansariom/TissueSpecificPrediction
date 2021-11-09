#!/bin/bash

# $1 is the following file: ""
#        https://krishna.gs.washington.edu/content/members/vagar/Xpresso/data/datasets/hg38_promoters_cage_corrected.bed.gz 
# $2 is the "cage_promoter_3K.bed"

cat $1 | awk -F , '(NR>1){n = $11"_"$1"_"int($6)"_"$2; s = $6-3000;e = $6+3000; printf "%s\t%d\t%d\t%s\t%s\n", $1, s, e, n, $2}' > $2
