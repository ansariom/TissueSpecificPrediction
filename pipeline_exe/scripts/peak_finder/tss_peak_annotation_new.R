#!/usr/bin/Rscript

leaf_peaks_file <- "~/Downloads/ibdc/Aug2018//aligned.peaks.annotated_leaf.capped"
root_peaks_file <- "~/Downloads/ibdc/Aug2018//aligned.peaks.annotated_root.capped"
gff_file <- "~/Downloads/ibdc/peaks/TAIR10_GFF3_genes.gff"

args <- commandArgs(trailingOnly = T)
leaf_peaks_file = args[1]
root_peaks_file = args[2]
gff_file <- args[3]
outfile <- args[4]

min_reads <- 50

############################################################
### Assign peaks to transcripts based on their distance
############################################################
assign_peaks <- function(df) {
  #df <- others[others$GeneName == "AT1G52240",]
  df <- df[order(abs(df$distance)),]
  df <- df[abs(df$distance) < 250,]
  peaks <- unique(df$tss_name)
  trx <- unique(df$TranscriptID)
  levels <- unique(df$distance)
  selected <- c()
  all_cands <- data.frame(GeneName=c(), tss_name=c(), TranscriptID = c())
  for (i in levels) {
    cands <- df[df$distance == i,]
    cands <- cands[!cands$TranscriptID %in% selected,]
    selected <- c(selected, cands[, "TranscriptID"])
    all_cands <- rbind(all_cands, cands[,c("GeneName", "tss_name", "TranscriptID")])
  }
  all_cands
}

leaf_peaks <- read.delim(leaf_peaks_file, sep = ",", header = T)
root_peaks <- read.delim(root_peaks_file, sep = ",", header = T)
gff <- read.delim(gff_file, sep = "\t", header = T)

gff <- gff[, c(1,3,4,5,7,9)]
colnames(gff) <- c("chr", "type", "start", "end", "strand", "ID")
gff <- gff[gff$type == "mRNA",]
library(tidyr)
gff <- extract(gff, ID, c("gene_id", "transcript_id"), regex = "([^;]+);([^;]+);\\D+")
gff <- extract(gff, gene_id, "TranscriptID", regex = "ID=([^;]+)")
gff <- extract(gff, transcript_id, "GeneName", regex = "Parent=([^;]+)")


same_start <- merge(gff, gff, by = "start")
same_start <- same_start[same_start$GeneName.x != same_start$GeneName.y,]
same_start_trx <- unique(same_start$TranscriptID.x)

leaf_peaks <- leaf_peaks[leaf_peaks$ReadCount > min_reads & leaf_peaks$TranscriptLocation != "<500",]
root_peaks <- root_peaks[root_peaks$ReadCount > min_reads & root_peaks$TranscriptLocation != "<500",]

############################################################
### Annotate Peaks
############################################################
get_annotated_peaks <- function(in_peaks, gff) {
  print("assign peaks")
  print(dim(gff))
  trx_no <- aggregate(gff$GeneName, list(gff$GeneName), length)
  colnames(trx_no) <- c("GeneName", "nTranscripts")
  
  peaks_ntranscripts <- aggregate(in_peaks$GeneName, list(in_peaks$GeneName), length)
  colnames(peaks_ntranscripts) <- c("GeneName", "nTranscripts_peaks")
  
  ref_trx_no <- merge(trx_no, peaks_ntranscripts, by = "GeneName")
  
  ## peaks covering one transcripts out of total 1 transcript
  true_one <- ref_trx_no[ref_trx_no$nTranscripts == 1 & ref_trx_no$nTranscripts_peaks >= 1,]
  peaks_g1 <- in_peaks[in_peaks$GeneName %in% true_one$GeneName,]
  print(dim(peaks_g1))
  
  ## Peaks covering one out of n true transcripts (assign all others to the peak)
  peaks1_missed_others <- ref_trx_no[ref_trx_no$nTranscripts_peaks == 1 & ref_trx_no$nTranscripts > 1,]
  peaks_g2 <- in_peaks[in_peaks$GeneName %in% peaks1_missed_others$GeneName,]
  peaks_g2$TranscriptID <- NULL
  peaks_g2 <- merge(peaks_g2, gff[,c("GeneName", "TranscriptID")], by = "GeneName")
  print(dim(peaks_g2))
  
  
  g12 <- rbind(peaks_g1, peaks_g2)
  
  ## Peaks that the transcripts are not assigned according to transcripts
  others <- in_peaks[!in_peaks$GeneName %in% g12$GeneName,]
  others <- unique(merge(others[,c("GeneName", "Start", "End", "ReadCount", "Strand")], gff[,c("GeneName", "TranscriptID", "start", "end")], by = "GeneName"))
  others$distance <- ifelse(others$Strand == "+" , others$start - others$Start, others$Start - others$end)
  others$tss_name <- paste(others$GeneName, others$Start, sep = "_")
  g3 <- aggregate(others$distance, list(others$tss_name), function(x) {return(var(x))})
  eq_dist <- g3[g3$x == 0,]
  others_eqDistance <- others[others$tss_name %in% eq_dist$Group.1,]
  peaks_g3 <- in_peaks[in_peaks$GeneName %in% others_eqDistance$GeneName,]
  peaks_g3$TranscriptID <- NULL
  peaks_g3 <- unique(merge(peaks_g3, gff[,c("GeneName", "TranscriptID")], by = "GeneName"))
  print(dim(peaks_g3))
  
  g123 <- rbind(g12, peaks_g3)
  
  others <- others[!others$tss_name %in% eq_dist$Group.1,]
  all_other_peaks <- do.call(rbind, by(others, others[, "GeneName"], assign_peaks))
  all_other_peaks <- merge(all_other_peaks, others, by = c("GeneName", "tss_name", "TranscriptID"))
  
  in_peaks$TranscriptID <- NULL
  all_other_peaks <- all_other_peaks[,c("GeneName", "TranscriptID", "Start", "End")]
  all_other_peaks <- merge(all_other_peaks, in_peaks, by = c("GeneName", "Start", "End"))
  
  all_peaks <- rbind(g123, all_other_peaks)
  
}

leaf_peaks <- get_annotated_peaks(leaf_peaks, gff)
print(dim(leaf_peaks))
root_peaks <- get_annotated_peaks(root_peaks, gff)
print(dim(root_peaks))

leaf_peaks$tissue <- "leaf"
root_peaks$tissue <- "root"
peaks <- rbind(leaf_peaks, root_peaks)
write.csv(peaks, file = outfile, quote = F, row.names = F, col.names = T)

