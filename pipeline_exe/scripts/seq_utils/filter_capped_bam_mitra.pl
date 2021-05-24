#!/usr/bin/perl

use warnings;
use strict;

use lib '/nfs0/BPP/Megraw_Lab/cumbiej/perl/modules';
use ReadDelim;
use BioFunc;

my $infile = shift;
my $outfile = shift;

my %genome;
BioFunc::readFasta("/nfs0/BPP/Megraw_Lab/cumbiej/dbs/tair10/genome.fas", \%genome);

my $current = 0;
my $reader = new ReadDelim('format' => 'SAM', 'save_lines' => 1);
open(SAM, "samtools view -h $infile |");
open(OUT, "| samtools view -bS - >$outfile");
while (my $line = <SAM>) {
	chomp $line;

	$current++;

	print STDERR "\r$current" if $current % 1000 == 0 || eof(SAM);
	print STDERR "\n" if eof(SAM);

	if ($line =~ m/^\@.{2}\t/) {
		print OUT "$line\n";
		next;
	}

	my %hit = $reader->parse_line($line);

	my ($base) = $hit{'qname'} =~ m/([ACGTN])\)$/;

	if ($base ne 'G') {
#		print "$base\n";
#		print OUT "$hit{'line'}\n";
		next;
	}

	my $ref = $hit{'rname'};
	my $start = $hit{'pos'};
	my $end = $start - 1;
	while ($hit{'cigar'} =~ m/([0-9]+)[NM]/g) {
		$end += $1;
	}
	my $strand = $hit{'flag'} & 16;
	$strand = ($strand)? "-" : "+";

	my $ref_base = '';
	if ($strand eq '+') {
		$ref_base = substr($genome{$ref}, $start - 2, 1);
	} else {
		$ref_base = substr($genome{$ref}, $end, 1);
		BioFunc::rc($ref_base);
	}

	if ($ref_base ne $base) {
		my @vals = split("\t", $hit{'line'});
		#$vals[0] =~ s/G\)$/G')/;
		print OUT join("\t", @vals)."\n";
#	} else {
#		print OUT "$hit{'line'}\n";
	}
}
close(OUT);
close(SAM);

