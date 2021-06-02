#!/bin/bash


cat $1 | awk -F , '(NR>1){n = $11"_"$1"_"int($6)"_"$2; s = $6-3000;e = $6+3000; printf "%s\t%d\t%d\t%s\t%s\n", $1, s, e, n, $2}' > $2
