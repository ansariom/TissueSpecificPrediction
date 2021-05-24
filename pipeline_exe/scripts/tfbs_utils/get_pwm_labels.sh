#!/bin/bash

input_pwms=$1
outfile=$2

rm -f $outfile
grep '>' $input_pwms | cut -c 3- > $outfile
