#!/usr/bin/perl
################################################################################
#
#  reannot_strand.pl
#
################################################################################

use strict;
use Getopt::Std;

# Jason's overlap search module to speed up annotator script
use lib '/nfs0/BPP/Megraw_Lab/cumbiej/perl/modules';
use Overlap;
use ReadDelim;

my $overlap_obj = new Overlap();
my $gff_reader = new ReadDelim('save_lines' => 1, 'no_headers' => 1, 'use_headers' => 'seqid,source,type,start,end,score,strand,phase,attributes');

##### Define Parameters #####

our $seperator = ",";
our %dicGffOrient = (
  '+' => 1,
  '-' => -1,
  '.' => 0,
  );
our %getGffOrient = (
   1 => '+',
  -1 => '-',
   0 => '.',
  );

##### End Define Parameters #####

my %FilterItemsArray; # hash (by chromosome) of array of hashes holding annotated items such as genes
my %maxItemLength;    # length of longest annotated item on each chromosome

sub usage {
	print STDERR "USAGE: $0 [-do] [-l...] gff_file in_items > out_ft.tab\n";
    # -l  column label for the table of overlapping items, default "ITEM_OVERLAP_LIST"
    # -z  make output header same as input items header
    # -o  report overlaps only if orientations are matching
    # -G  table in gff format
    # -k  table in "knownGene.txt" format
    # -s  table in "sno.txt" format
    # -y  table in "YalePseudo.txt" format
    # -r  table in "RetroposedGenes.txt" format
    # -c  table in "CpGIslands.txt" format
    # -e  table in "EncodeRegions.txt" format
    # -w  table in "wgRna.txt" format
    # -m  table in "all_est.txt" format
    # -p  input in peak caller format (David's Peak Caller)
	exit 1;
}

# command line interface
our (%mopt);
getopts ('dhGksyrcewml:zop',\%mopt) or &usage();
if ($mopt{h} or @ARGV<2) { &usage() }
$mopt{d} and print STDERR "\n", join(' ',$0,@ARGV), "\n";
my $olap_label = $mopt{l} || 'ITEM_OVERLAP_LIST';


# Hashes to handle parsing of TAIR GFF formatted data
my %gff_gene_id_index;

my %gff_gene_types = (	'gene' => 1,
						'pseudogene' => 1,
						'transposable_element_gene' => 1);

my %gff_transcript_types = ('mRNA' => 1,
							'transcript' => 1,
							'miRNA' => 1,
							'ncRNA' => 1,
							'mRNA_TE_gene' => 1,
							'rRNA' => 1,
							'snoRNA' => 1,
							'snRNA' => 1,
							'tRNA' => 1,
							'pseudogenic_transcript' => 1 );

my %gff_part_types = (	'CDS' => 1,
						'exon' => 1,
						'five_prime_UTR' => 1,
						'pseudogenic_exon' => 1,
						'three_prime_UTR' => 1 );

# read table
open (INTABLE,$ARGV[0]) or die "failed to open input file $ARGV[0]";
while (<INTABLE>) {
	chomp;
	my @fdgff = split(/\t/);
	my $pItem;
	my $joinedItemText = join("|", @fdgff);
	my $chr;
	my $type;
	my $trans;
	my $start;
	my $end;
	my $strand;
	my $info;

	if ($mopt{G}) {
		my %feature = $gff_reader->parse_line($_);
		my %attributes;
		foreach my $pair (split(";", $feature{'attributes'})) {
			my ($key, $value) = split("=", $pair);
			$key =~ s/^\s+//;
			$key =~ s/\s+$//;
			$value =~ s/^\s+//;
			$value =~ s/\s+$//;
			$attributes{$key} = $value;
		}
		$feature{'attributes'} = \%attributes;

		if (exists($gff_gene_types{$feature{'type'}})) {
			$feature{'transcripts'} = {};
			if (!exists($FilterItemsArray{$feature{'seqid'}})) {
				$FilterItemsArray{$feature{'seqid'}} = [];
			}
			push @{$FilterItemsArray{$feature{'seqid'}}}, \%feature;
			$gff_gene_id_index{$feature{'attributes'}->{'ID'}} = scalar(@{$FilterItemsArray{$feature{'seqid'}}}) - 1;
		} elsif (exists($gff_transcript_types{$feature{'type'}})) {
			$feature{'parts'} = [];
			my $gene_id = $feature{'attributes'}->{'Parent'};
			my $transcript_id = $feature{'attributes'}->{'ID'};
			my $gene_index = $gff_gene_id_index{$gene_id};
			$FilterItemsArray{$feature{'seqid'}}->[$gene_index]->{'transcripts'}->{$transcript_id} = \%feature;
		} elsif (exists($gff_part_types{$feature{'type'}})) {
			my ($gene_id) = $feature{'attributes'}->{'Parent'} =~ m/^([^,]+)\.[0-9]+/;
			my $transcript_id = (split(",", $feature{'attributes'}->{'Parent'}))[0];
			my $gene_index = $gff_gene_id_index{$gene_id};
			push @{$FilterItemsArray{$feature{'seqid'}}->[$gene_index]->{'transcripts'}->{$transcript_id}->{'parts'}}, \%feature;
		} else {
			next;
		}
	}
	elsif ($mopt{k}) {
	    $chr = &chr_char($fdgff[1]);
	    push @{$FilterItemsArray{$chr}}, $pItem = {
		itemName => $fdgff[0],
		itemChrom => $chr,
		itemOrient => $dicGffOrient{$fdgff[2]},
		itemStart => $fdgff[3]-1,  # start position stored in computational notation
		itemEnd => $fdgff[4],
		itemScore => "",
		itemFeature => "",
	    };
	}
	elsif ($mopt{s}) {
	    $chr = &chr_char($fdgff[1]);
	    push @{$FilterItemsArray{$chr}}, $pItem = {
		itemName => $fdgff[4],
		itemChrom => $chr,
		itemOrient => $dicGffOrient{$fdgff[6]},
		itemStart => $fdgff[2]-1,  # start position stored in computational notation
		itemEnd => $fdgff[3],
		itemScore => $fdgff[5],
		itemFeature => $fdgff[9],
	    };
	}
	elsif ($mopt{y}) {
	    $chr = &chr_char($fdgff[1]);
	    push @{$FilterItemsArray{$chr}}, $pItem = {
		itemName => $fdgff[0],
		itemChrom => $chr,
		itemOrient => $dicGffOrient{$fdgff[2]},
		itemStart => $fdgff[3]-1,  # start position stored in computational notation
		itemEnd => $fdgff[4],
		itemScore => "",
		itemFeature => "",
	    };
	}
	elsif ($mopt{r}) {
	    $chr = &chr_char($fdgff[1]);
	    push @{$FilterItemsArray{$chr}}, $pItem = {
		itemName => $fdgff[4],
		itemChrom => $chr,
		itemOrient => $dicGffOrient{$fdgff[6]},
		itemStart => $fdgff[2]-1,  # start position stored in computational notation
		itemEnd => $fdgff[3],
		itemScore => $fdgff[5],
		itemFeature => "",
	    };
	}
	elsif ($mopt{c}) {
	    $chr = &chr_char($fdgff[0]);
	    push @{$FilterItemsArray{$chr}}, $pItem = {
		itemName => $fdgff[3],
		itemChrom => $chr,
		itemOrient => "",
		itemStart => $fdgff[1]-1,  # start position stored in computational notation
		itemEnd => $fdgff[2],
		itemScore => "",
		itemFeature => "",
	    };
	}
	elsif ($mopt{e}) {
	    $chr = &chr_char($fdgff[0]);
	    push @{$FilterItemsArray{$chr}}, $pItem = {
		itemName => $fdgff[3],
		itemChrom => $chr,
		itemOrient => "",
		itemStart => $fdgff[1]-1,  # start position stored in computational notation
		itemEnd => $fdgff[2],
		itemScore => "",
		itemFeature => "",
	    };
	}
	elsif ($mopt{w}) {
	    $chr = &chr_char($fdgff[1]);
	    $type = &chr_char($fdgff[9]);
	    unless ($type eq "miRna") {
		push @{$FilterItemsArray{$chr}}, $pItem = {
		    itemName => $fdgff[4],
		    itemChrom => $chr,
		    itemOrient => $dicGffOrient{$fdgff[6]},
		    itemStart => $fdgff[2]-1,  # start position stored in computational notation
		    itemEnd => $fdgff[3],
		    itemScore => "",
		    itemFeature => "",
		};
	    }
	}
	elsif ($mopt{m}) {
	    $chr = &chr_char($fdgff[14]);
	    push @{$FilterItemsArray{$chr}}, $pItem = {
		itemName => $fdgff[10],
		itemChrom => $chr,
		itemOrient => $dicGffOrient{$fdgff[9]},
		itemStart => $fdgff[16]-1,  # start position stored in computational notation
		itemEnd => $fdgff[17],
		itemScore => "",
		itemFeature => "",
	    };
	}
	else {
	    die "failed to recognize format option for $ARGV[0]";
	}

	# update global dictionary of maximum lengths (per chromosome)
	unless (exists($maxItemLength{$chr})) {
		$maxItemLength{$chr} = 0;
	}
	my $length = $pItem->{itemEnd} - $pItem->{itemStart};
	if ($length > $maxItemLength{$chr}) {
		$maxItemLength{$chr} = $length;
	}
}
close (INTABLE);
#print STDERR "$ARGV[0] loading complete\n";


# sort the GFF items, to get ready for binary search
foreach my $key (keys(%FilterItemsArray)) {
	if ($mopt{G}){
		@{$FilterItemsArray{$key}} = sort { $a->{start} <=> $b->{start} } @{$FilterItemsArray{$key}};
	} else {
		@{$FilterItemsArray{$key}} = sort { $a->{itemStart} <=> $b->{itemStart} } @{$FilterItemsArray{$key}};
	}
}
#print STDERR "GFF annotation sort complete\n";

# loop over input items, report features
open(ITEMS, $ARGV[1]) or die "failed to open input file $ARGV[1]";
if ($mopt{p}) {
	<ITEMS>; # read out header
}
my $header = 'Chromosome,Strand,Start,End,ReadCount,ModeLocation,ModeReadCount,Shape,TranscriptLocation,TranscriptID,GeneName';
$mopt{z} = 1;
if ($mopt{z}) { # start feature table header with column labels
    chomp($header);
    my $newheader = $header.",GeneType";
    print $newheader, "\n";
} else {
    print $header, "\n";
}
my $itemCnt = 0;
while (<ITEMS>) {
    my $line = $_;
    chomp($line);
    next if $line =~ /^\#/;  # skip comment lines

    #Chr1  23095 23099 Chr1.21  0  +  3.31981874  -1 -1 2  8
    my @tmpVals = split("\t", $line);

    if (!$mopt{p}){
        $line = "$tmpVals[0],$tmpVals[5],$tmpVals[1],$tmpVals[2],$tmpVals[10],$tmpVals[9],NA,NA,NA,NA,NA";
    } else {
        $line .= ",NA,NA,NA";
    }

    my $id = $line;                           # peak id
    my @itemarr = split(/\,/, $line);
    my $chr = $itemarr[0];                    # peak full chromosome name
    my $start = $itemarr[2];                  # peak start in informatics notation
    my $end = $itemarr[3];                    # peak end in informatics notation
    my $strand = $itemarr[1];                 # strand (+/-)
    my $orient = $dicGffOrient{$itemarr[1]};  # peak orientation, numerical notation
    my $mode = $start+$itemarr[5];            # peak mode

    # David Peak Caller gives absolute coordinates for mode
    if ($mopt{p}) {
        $mode = $itemarr[5];
    }

    my %priority = ('tss', 1,
		    'five_prime_UTR', 2,
		    'gene-upstream250', 3,
		    'miRNA-upstream250', 3,
		    'gene-upstream500', 4,
		    'miRNA-upstream500', 4,
		    'gene-upstream1000', 5,
		    'miRNA-upstream1000', 5,
		    'gene-upstream1500', 6,
		    'miRNA-upstream1500', 6,
		    'gene-upstream2000', 7,
		    'miRNA-upstream2000', 8,
		    'gene-upstream2500', 8,
		    'miRNA-upstream2500', 8,
		    'gene-upstream3000', 9,
		    'miRNA-upstream3000', 9,
		    'CDS', 10,
		    'three_prime_UTR', 11,
		    'exon', 12,
		    'gene', 13,
		    'intergenic', 14);
    
    my %outlabel = (1, 'tss',
		    2, '5\'utr',
		    3, '<250',
		    4, '<500',
		    5, '<1000',
		    6, '<1500',
		    7, '<2000',
		    8, '<2500',
		    9, '<3000',
		    10, 'coding',
		    11, '3\'utr',
		    12, 'exon-nonCDS',
		    13, 'genic-other',
		    14, 'intergenic');
    
    my $maxPriority = 14;

    # list of annotated items overlapping input item
    my @overlaplist;

	# 40000 in find_overlap function indicates the longest 'item' in list isn't > 40000, and 1500 ensures
	# that I pullout items within 3kb of each other (extends item list by 1500nt on each side, and the 
	# search item by 1500nt for a total of 3kb combined)
	my @overlaps = $overlap_obj->find_overlap({'start' => $start, 'end' => $end}, $FilterItemsArray{$chr}, 40000, 1500);

	# Only keep genes that are on the same strand as peak and are 'upstream' (relative to the gene) or overlap the gene
	foreach my $overlap (@overlaps) {
		my $gene = $overlap->[1];

		# Ensure on the same strand
		next if $gene->{'strand'} ne $strand;

		# Ensure no peaks are "downstream" (relative to gene orientation) of gene
		next if $gene->{'strand'} eq '+' && $start > $gene->{'end'};
		next if $gene->{'strand'} eq '-' && $end < $gene->{'start'};

		my $isupstream = 0;

		if ($gene->{'strand'} eq '+') {
			$isupstream = 1 if $end < $gene->{'start'};
		} else {
			$isupstream = 1 if $start > $gene->{'end'};
		}

		# Ignore "promoter" regions of pseudogenes and transposable_element_genes
		next if $isupstream && ($gene->{'type'} eq 'pseudogene' || $gene->{'type'} eq 'transposable_element_gene');

		if ($gene->{'start'} <= $end && $gene->{'end'} >= $start) {
			$gene->{'direct_overlap'} = 1;
		} else {
			$gene->{'direct_overlap'} = 0;
		}

		push @overlaplist, $gene;
	}

    # TAIR transcript part priority list for sort
    my %txnpartsort = ('five_prime_UTR', 6,
		    'CDS', 5,
		    'three_prime_UTR', 4,
		    'exon', 3,
		    'protein', 2,
		    'intron', 1);


    # sort overlaplist by proximity of each item to peak mode
	my @sortedoverlaplist;
	my @sortedpartlist;
	if ($strand eq '+') {
		#prioritize genes that overlap over those that don't
		@sortedoverlaplist = sort { 
									if ($a->{'direct_overlap'} == $b->{'direct_overlap'}) {
										my $a_upstream = ($end < $a->{'start'})? 1 : 0; my $b_upstream = ($end < $b->{'start'})? 1 : 0;
										if ($a_upstream == $b_upstream) {
											abs($a->{'start'} - $mode) <=> abs($b->{'start'} - $mode);
										} else {
											$b_upstream <=> $a_upstream;
										}
									} else {
										$b->{'direct_overlap'} <=> $a->{'direct_overlap'};
									}
								  } @overlaplist;
	} else {
		#prioritize genes that overlap over those that don't
		@sortedoverlaplist = sort {
									if ($a->{'direct_overlap'} == $b->{'direct_overlap'}) {
										my $a_upstream = ($start > $a->{'end'})? 1 : 0; my $b_upstream = ($start > $b->{'end'})? 1 : 0;
										if ($a_upstream == $b_upstream) {
											abs($a->{'end'} - $mode) <=> abs($b->{'end'} - $mode);
										} else {
											$b_upstream <=> $a_upstream;
										}
									} else {
										$b->{'direct_overlap'} <=> $a->{'direct_overlap'};
									}
								  } @overlaplist;
	}

	# Find how many genes directly overlap our peak
	my $total_overlap = 0;
	foreach my $gene (@sortedoverlaplist) {
		if ($strand eq '+') {
			$total_overlap++ if $end >= $gene->{'start'};
		} else {
			$total_overlap++ if $start <= $gene->{'end'};
		}
	}

	# If no overlaps, check for cases where two genes start equally far from peak
	my $total_same_dist = 1;
	if (!$total_overlap) {
		my $prev = -1;
		foreach my $gene (@sortedoverlaplist) {
			if ($strand eq '+') {
				if ($prev != -1) {
					my $curr = $gene->{'start'} - $end;
					if ($prev == $curr) {
						$total_same_dist++;
					} else {
						last;
					}
				}
				$prev = $gene->{'start'} - $end;
			} else {
				if ($prev != -1) {
					my $curr = $start - $gene->{'end'};
					if ($prev == $curr) {
						$total_same_dist++;
					} else {
						last;
					}
				}
				$prev = $start - $gene->{'end'};
			}
		}
	}

	# TODO: implement code to handle peaks that have multiple best genes...only affects <0.5% of peaks for now...
	if ($total_overlap > 1 || $total_same_dist > 1) {

	# There is *one* best gene, but probably not one best transcript
	} elsif (scalar(@sortedoverlaplist)) {
		my $total_txn = scalar(keys(%{$sortedoverlaplist[0]->{'transcripts'}}));
		if ($total_txn > 1) {
			while (my ($txn_id, $info) = each %{$sortedoverlaplist[0]->{'transcripts'}}) {
				foreach my $part (@{$info->{'parts'}}) {
					push @sortedpartlist, $part;
				}
				my @exons = grep { $_->{'type'} =~ /exon/; } @{$info->{'parts'}};
				@exons = sort { $a->{'start'} <=> $b->{'start'} } @exons;
				for (my $x = 0; $x < scalar(@exons) - 1; $x++) {
					push @sortedpartlist, {'seqid' => $exons[$x]->{'seqid'},
												  'type'  => 'intron',
												  'start' => $exons[$x]->{'end'} + 1,
												  'end' => $exons[$x + 1]->{'start'} - 1,
												  'strand' => $exons[$x]->{'strand'},
												  'attributes' => $exons[$x]->{'attributes'}};
				}
			}
			if ($strand eq '+') {
				@sortedpartlist = sort { my $a_overlap = ($a->{'end'} >= $start && $a->{'start'} <= $end)? 1 : 0;
										 		 my $b_overlap = ($b->{'end'} >= $start && $b->{'start'} <= $end)? 1 : 0;
												 $b_overlap <=> $a_overlap ||
                                    ($txnpartsort{$b->{'type'}} <=> $txnpartsort{$a->{'type'}}) ||
												(abs($a->{'start'} - $mode) <=> abs($b->{'start'} - $mode))
									} @sortedpartlist;
			} else {
				@sortedpartlist = sort { my $a_overlap = ($a->{'end'} >= $start && $a->{'start'} <= $end)? 1 : 0;
										 		 my $b_overlap = ($b->{'end'} >= $start && $b->{'start'} <= $end)? 1 : 0;
									    		 $b_overlap <=> $a_overlap ||
                                     ($txnpartsort{$b->{'type'}} <=> $txnpartsort{$a->{'type'}}) || 
												 (abs($a->{'end'} - $mode) <=> abs($b->{'end'} - $mode))
										} @sortedpartlist;
			}
		}
	}
    
    my $geneLabel = 'NA';
    my $strandLabel = 'NA';
    my $transcriptLabel = 'NA';
    my $labelPriority = $maxPriority;
	my $struct_label = '';
    my $transcriptType = 'NA';

    # Scroll down the overlap list until an item with a gene ID is found

    my $annotItem;
    my $geneid;
    my $itemname;
    my $itemStrand;
    my $geneInfo;
    my @geneInfoParts;
    my $idpart;
    
    my $foundGeneID = 0;
    my $tssOverlap = 0;
    my $type = 'intergenic';

	if (scalar(@sortedoverlaplist)) {
		$foundGeneID = 1;
		$geneid = $sortedoverlaplist[0]->{'attributes'}->{'ID'};
		$geneLabel = $geneid;
		$type = $sortedoverlaplist[0]->{'type'};
		if ($type eq 'gene') {
			my $gene = $sortedoverlaplist[0];
			my $txn_type = $gene->{'transcripts'}->{(keys(%{$gene->{'transcripts'}}))[0]}->{'type'};
			if ($txn_type ne 'mRNA') {
				$type = $txn_type;
			}
		}
		$itemStrand = $sortedoverlaplist[0]->{'strand'};
	}

	# What part of the gene is overlapped?
    if($foundGeneID) {
		my $genePart = 'NA';
		#$transcriptType = $type;
		my @txn_ids = keys(%{$sortedoverlaplist[0]->{'transcripts'}});
		$transcriptType = $sortedoverlaplist[0]->{'transcripts'}->{$txn_ids[0]}->{'type'};
		if ($transcriptType eq 'mRNA') {
			$transcriptType = 'gene';
		} else {
			if ($sortedoverlaplist[0]->{'type'} eq 'pseudogene' || $sortedoverlaplist[0]->{'type'} eq 'transposable_element_gene') {
				$transcriptType = $sortedoverlaplist[0]->{'type'};
			}
		}

		# store strandLabel
		$strandLabel = $itemStrand;

		# Does the peak overlap the annotated TSS of the item with gene ID
		if ($itemStrand eq '+') {
		    if ($start <= $sortedoverlaplist[0]->{'start'} && $end >= $sortedoverlaplist[0]->{'start'}) {
				$tssOverlap = 1;
		    }
			if (!$tssOverlap && $end < $sortedoverlaplist[0]->{'start'}) {
				my $dist = $sortedoverlaplist[0]->{'start'} - $end;
				if ($transcriptType ne 'gene' && $transcriptType ne 'miRNA') {
					$type = 'gene';
				} else {
					$type = $transcriptType;
				}
				if ($dist < 250) {
					$type .= '-upstream250';
				} elsif ($dist < 500) {
					$type .= '-upstream500';
				} elsif ($dist < 1000) {
					$type .= '-upstream1000';
				} elsif ($dist < 1500) {
					$type .= '-upstream1500';
				} elsif ($dist < 2000) {
					$type .= '-upstream2000';
				} elsif ($dist < 2500) {
					$type .= '-upstream2500';
				} elsif ($dist < 3000) {
					$type .= '-upstream3000';
				} else {
					$type = 'NA';
				}

				# label the part of promoter for non-mRNA/miRNA genes (mostly structural RNAs)
				if ($transcriptType ne 'gene' && $transcriptType ne 'miRNA' && $type ne 'NA') {
					$struct_label = $outlabel{$priority{$type}};
				}
			}
		} else {
		    if ($start <= $sortedoverlaplist[0]->{'end'} && $end >= $sortedoverlaplist[0]->{'end'}) {
				$tssOverlap = 1;
	    	}
			if (!$tssOverlap && $start > $sortedoverlaplist[0]->{'end'}) {
				my $dist = $start - $sortedoverlaplist[0]->{'end'};
				if ($transcriptType ne 'gene' && $transcriptType ne 'miRNA') {
					$type = 'gene';
				} else {
					$type = $transcriptType;
				}
				if ($dist < 250) {
					$type .= '-upstream250';
				} elsif ($dist < 500) {
					$type .= '-upstream500';
				} elsif ($dist < 1000) {
					$type .= '-upstream1000';
				} elsif ($dist < 1500) {
					$type .= '-upstream1500';
				} elsif ($dist < 2000) {
					$type .= '-upstream2000';
				} elsif ($dist < 2500) {
					$type .= '-upstream2500';
				} elsif ($dist < 3000) {
					$type .= '-upstream3000';
				} else {
					$type = 'NA';
				}

				# label the part of promoter for non-mRNA/miRNA genes (mostly structural RNAs)
				if ($transcriptType ne 'gene' && $transcriptType ne 'miRNA' && $type ne 'NA') {
					$struct_label = $outlabel{$priority{$type}};
				}
			}
		}

		# Add 'upstream' info to type names for gene or miRNA
		if ($transcriptType eq 'gene' || $transcriptType eq 'miRNA') {
			$transcriptType = $type;
		} elsif ($type eq 'pseudogene' || $type eq 'transposable_element_gene') {
			if (!$tssOverlap) {
				$genePart = 'genic-other';
			} else {
				$genePart = 'tss';
			}

		# Get transcript label for upstream peaks
		} elsif (!($type =~ m/upstream/)) {
			if ($tssOverlap) {
				$genePart = 'tss';
			}
			my @transcripts;
			while (my ($txn_id, $txn) = each %{$sortedoverlaplist[0]->{'transcripts'}}) {
				push @transcripts, $txn;
			}
			if ($strand eq '+') {
				@transcripts = sort { abs($a->{'start'} - $mode) <=> abs($b->{'start'} - $mode) } @transcripts;
			} else {
				@transcripts = sort { abs($a->{'end'} - $mode) <=> abs($b->{'end'} - $mode) } @transcripts;
			}
			$transcriptLabel = $transcripts[0]->{'attributes'}->{'ID'};
		}
	
		# If item with gene ID is a gene, look for nearest transcript and identify gene part in subsequent list items
		if ($type eq 'gene') { 
		    # What is the transcript ID associated with this gene that has the closest start to peak mode
			my $gene = $sortedoverlaplist[0];
			my $ismRNA_test = $gene->{'transcripts'}->{(keys(%{$gene->{'transcripts'}}))[0]}->{'type'};
		    my $ismRNA = ($ismRNA_test eq 'mRNA')? 1 : 0;

		    # If no TSS overlap, what is the highest priority gene part associated with this gene that overlaps the peak
		    if ($ismRNA) {
				if ($tssOverlap) {
					$genePart = 'tss';
					my @transcripts;
					while (my ($txn_id, $txn) = each %{$sortedoverlaplist[0]->{'transcripts'}}) {
						push @transcripts, $txn;
					}
					if ($strand eq '+') {
						@transcripts = sort { abs($a->{'start'} - $mode) <=> abs($b->{'start'} - $mode) } @transcripts;
					} else {
						@transcripts = sort { abs($a->{'end'} - $mode) <=> abs($b->{'end'} - $mode) } @transcripts;
					}
					$transcriptLabel = $transcripts[0]->{'attributes'}->{'ID'};
				} else {
					if (scalar(@sortedpartlist)) {
						#($transcriptLabel) = $sortedpartlist[0]->{'attributes'}->{'Parent'} =~ m/(AT[1-5CM]G[0-9]+\.[0-9]+)/;
						($transcriptLabel) = $sortedpartlist[0]->{'attributes'}->{'Parent'} =~ m/([^\.]+\.[0-9]+)/;
						$genePart = $sortedpartlist[0]->{'type'};
					} else {
						my @transcripts;
						while (my ($txn_id, $txn) = each %{$sortedoverlaplist[0]->{'transcripts'}}) {
							push @transcripts, $txn;
						}
						$transcriptLabel = $transcripts[0]->{'attributes'}->{'ID'};
						if ($strand eq '+') {
							@transcripts = sort { abs($a->{'start'} - $mode) <=> abs($b->{'start'} - $mode) } @transcripts;
						} else {
							@transcripts = sort { abs($a->{'end'} - $mode) <=> abs($b->{'end'} - $mode) } @transcripts;
						}
						my @utrs = grep { $_->{'type'} =~ m/UTR/ } @{$transcripts[0]->{'parts'}};
						my @exons = grep { $_->{'type'} =~ m/exon/ } @{$transcripts[0]->{'parts'}};
						foreach my $utr (@utrs) {
							if ($end >= $utr->{'start'} && $start <= $utr->{'end'}) {
								$genePart = $utr->{'type'};
								last;
							}
						}
						if ($genePart eq 'NA') {
							foreach my $exon (@exons) {
								if ($end >= $exon->{'start'} && $start <= $exon->{'end'}) {
									$genePart = 'CDS';
									last;
								}
							}
							if ($genePart eq 'NA') {
								$genePart = 'intron';
							}
						}
					}
				}
	    	}
		} elsif ($type =~ /miRNA/) {
			my @transcripts;
			while (my ($txn_id, $txn) = each %{$sortedoverlaplist[0]->{'transcripts'}}) {
				push @transcripts, $txn;
			}
			if ($strand eq '+') {
				@transcripts = sort { abs($a->{'start'} - $mode) <=> abs($b->{'start'} - $mode) } @transcripts;
			} else {
				@transcripts = sort { abs($a->{'end'} - $mode) <=> abs($b->{'end'} - $mode) } @transcripts;
			}
			$transcriptLabel = $transcripts[0]->{'attributes'}->{'ID'};
		}

		# Assign output labels
		my $thisLabel;
		if ($tssOverlap) {
		    $thisLabel = 'tss';
		} elsif ($genePart ne 'NA') {
		    $thisLabel = $genePart;
		} else {
		    $thisLabel = $type;
		}

		my $label;
		if ($thisLabel ne 'intron') {
			unless (exists($priority{$thisLabel})) {
			    $thisLabel = 'gene';
			}
			$labelPriority = $priority{$thisLabel};
			$label = $outlabel{$labelPriority};
		} else {
			$label = 'intron';
		}

		# Use the structural RNA label if found earlier to denote position
		$label = $struct_label if $struct_label ne '';

		# Record overlapping gene ID item
		$itemarr[1] = $strandLabel;
		$itemarr[8] = $label;
		$itemarr[9] = $transcriptLabel;
		$itemarr[10] = $geneLabel;
		$itemarr[11] = $transcriptType;

		my $outitem = join(",", @itemarr);
		print $outitem, "\n";

	# intergenic peak
    } else { 
		my $thisLabel = 'intergenic';
		$strandLabel = $strand;
		if (exists($priority{$thisLabel})) {
		    $labelPriority = $priority{$thisLabel};
		}
		my $label = $outlabel{$labelPriority};
		$itemarr[1] = $strandLabel;
		$itemarr[8] = $label;
		$itemarr[9] = $transcriptLabel;
		$itemarr[10] = $geneLabel;
		$itemarr[11] = $transcriptType;

		my $outitem = join(",", @itemarr);
		print $outitem, "\n";	
    }

	# Update user on where we are at in processing peak file...
    $itemCnt++;
    if (!($itemCnt%1000)) {
		print STDERR "\r                                                                                                                    ";
		print STDERR "\rProcessing item $itemCnt: $id";
    }
    print STDERR "\n" if eof(ITEMS);
}
close (ITEMS);

sub chr_char {
    (my $chr=$_[0]) =~ s/^chr//i;
    return $chr;
}

sub name_chars {
    my $name = shift;
    if ($name =~ /ID\=\"(.*)\"/ ) {
	return $1;
    }
    else {
	return $name;
    }
}

##### Subroutine Scan #########################################################################################################
#
# Args: item id, item chromosome name, item start, item end, item orient
# Output: list of overlap candidates

sub Scan {
	my $leftBndryIdx;
	my $rightBndryIdx;

	# Get arguments
	# note that positions are in computational notation
	my ($id,$chrom,$start,$end,$orient) = @_;

	if (abs($orient) != 1) {
		print STDERR "$0:Scan:Unrecognized case of strand orientation ($orient) for item $id\n";
	}
	if ($start > $end) {
		print STDERR "$0:Scan:start ($start) greater than end ($end) for item $id\n";
	}

	my $leftBndry = $start - $maxItemLength{$chrom};
	my $rightBndry = $end;

	# Scan annotated items for a candidate matching locations
	($leftBndryIdx, $rightBndryIdx)
	  = &binarySeek ($FilterItemsArray{$chrom}, 'itemStart', $leftBndry, $rightBndry);
	if ($rightBndryIdx < 0) {
		# no items for this chromosome in FilterItemsArray
		$leftBndryIdx = 0;
		$rightBndryIdx = -1;
	}
	elsif ($leftBndryIdx < 0) {
		# FilterItemsArray nonempty, but need to search from start
		$leftBndryIdx = 0;
	}

	# Check which candidate annotated items are indeed overlapping
	my @overlapList = ();
	for my $i ($leftBndryIdx .. $rightBndryIdx) {
		my $pItem = $FilterItemsArray{$chrom}[$i];

		if (!$mopt{o} or $orient==$pItem->{itemOrient}) {

			# does input item overlap this annotated item?
			# note that positions are in computational notation
		  if (($start>=$pItem->{itemStart} and $end<=$pItem->{itemEnd}) or
			    ($start<=$pItem->{itemStart} and $end>$pItem->{itemStart}) or
			    ($start<$pItem->{itemEnd}    and $end>=$pItem->{itemEnd}) ) {
#				push @overlapList, $pItem->{itemName};  # push only itemName (string of joined item text)
				push @overlapList, $pItem;  # push entire item (so overlapList is an array of hashes)
			}

		}
	}

	return @overlapList;
}

sub binarySeek {
# binary searches an array for a range(!), returns index of elt to the immediate
# left and index of elt to the immediate right of the range

	my $array_ref = shift;  # reference to array of hashes to be searched
	my $hashkey = shift;    # hashkey element to be searched (array should be sorted by the value of this hash element)
	my $floor = shift;      # left end of range to be searched
	my $ceiling = shift;    # right end of range to be searched

	my @answer;
	my $high;
	my $low;
	my $probe;

	# work on floor
	$high = $#{$array_ref};
	$low = -1;
	while ($high - $low > 1) {
		$probe = int(($high + $low)/2);
		if ($array_ref->[$probe]{$hashkey} < $floor) {
			$low = $probe;
		}
		else {
			$high = $probe;
		}
	}
	$answer[0] = $low;

	# work on ceiling
	$high = $#{$array_ref};
	$low = -1;
	while ($high - $low > 1) {
		$probe = int(($high + $low)/2);
		if ($array_ref->[$probe]{$hashkey} > $ceiling) {
			$high = $probe;
		}
		else {
			$low = $probe;
		}
	}
	$answer[1] = $high;

	return @answer;
}

