#!/usr/bin/perl

use warnings;
use strict;

use lib '/nfs0/BPP/Megraw_Lab/cumbiej/perl/modules';
use BioFunc;
use ReadDelim;

my ($tss_file, $genome_file, %genome, $header, $reader, $upstream, $downstream);

$tss_file = shift;
$upstream = shift;
$downstream = shift;

$genome_file = shift; #"/nfs0/BPP/Megraw_Lab/cumbiej/dbs/tair10/genome.fas";


BioFunc::readFasta($genome_file, \%genome);

open(DAT, $tss_file);

$header = <DAT>;
chomp $header;
$reader = new ReadDelim('no_header' => 1, 'use_headers' => $header, 'delim' => ',', 'save_lines' => 1);

while (my $line = <DAT>) {
	my %tss_peak = $reader->parse_line($line);

	my $tss = $tss_peak{'ModeLocation'};

	my $start = $tss - $upstream;
	my $end = $tss + $downstream;
	#$start = 1 if $start < 1;
	next if $start < 1;
	#$end = length($genome{$tss_peak{'Chromosome'}}) if $end > length($genome{$tss_peak{'Chromosome'}});
	next if $end > length($genome{$tss_peak{'Chromosome'}});

	next if $tss_peak{'GeneName'} =~ m/AT[CM]G/i; # ignore chloroplast and mitochondrial tss peaks

	my $pro_seq = substr($genome{$tss_peak{'Chromosome'}}, $start - 1, $end - $start + 1);
	BioFunc::rc($pro_seq) if $tss_peak{'Strand'} eq '-';
	print ">$tss_peak{'TranscriptID'}_$tss_peak{'Chromosome'}_$tss" . "_$tss_peak{'Strand'}\n";
	for (my $x = 0; $x < length($pro_seq); $x++) {
                print substr($pro_seq, $x, 1);
                print "\n" if ($x + 1) % 80 == 0 || $x == (length($pro_seq) - 1);
        }
}


close(DAT);
