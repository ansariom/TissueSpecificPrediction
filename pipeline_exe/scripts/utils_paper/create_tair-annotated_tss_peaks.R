#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
inpeak = args[1]
tair_start_file = args[2]
outfile = args[3]

#main_peaks <- read.csv("../pipeline_output_roe/aligned.peaks.annotated.capped.filtered")
main_peaks = read.csv(inpeak)
#tair_start <- read.table("data/tair10/TAIR10_GFF3_genes.tair_annotated_start.txt")
tair_start = read.table(tair_start_file)

colnames(tair_start) <- c("ModeLocation", "TranscriptID")
m = merge(main_peaks[,c("TranscriptID", "Chromosome", "Strand", "GeneName","ReadCount")], tair_start, by = "TranscriptID")

write.table(file = outfile, m[,c("TranscriptID", "Chromosome", "Strand", "GeneName","ModeLocation","ReadCount")], quote = F, col.names = T, row.names=F, sep = ",")

