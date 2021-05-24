#!/usr/bin/perl

use warnings;
use strict;

while(<STDIN>) {
	if (/^@.{2}\t/) {
		print $_;
	} else {
		chomp;
		my @values = split("\t", $_);
		next if $values[2] eq "*";
		print "$_\n";
	}
}
