#!/bin/bash

nucs_upstream=$1
nucs_downstream=$2
prom_upstream_length=$3
prom_downstream_length=$4
input_fasta=$5

if [ $# -lt 4 ]; then
	echo "./get_tss_region.sh [nucs_upstream] [nucs_downstream] [prom_up_len] [prom_down_len] [input_fasta]"
	exit
fi


grep ">" $input_fasta | sed 's/>//g' | awk -F ',' -v pd=$prom_downstream_length -v pu=$prom_upstream_length  -v u=$nucs_upstream -v d=$nucs_downstream '{split($0, a, "_"); tssmode = a[3]; strand = a[4]; if (strand == "-"){start = tssmode - d; end = tssmode + u} else {start = tssmode - u; end = tssmode + d}; print $1"_0\t"a[2]"\t"start"\t"end"\t.\t"strand} '

