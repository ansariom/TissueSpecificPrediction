library(tidyr)
a <- read.delim("high_scoring_TFBS_OC/scan_out/scan.FWD.locations.tbl", col.names = c("pwm", "seq", "loc"), sep="\t")
b <-  extract(a, seq, into = c("gene_id", "chr", "start", "strand"), regex = "(.+?)_(Chr\\d)_(\\d+)_(.)", remove = F)

options(scipen=99999999)

b$tss_loc <- as.numeric(b$start) + 3001
b$genomic_start <- ifelse(b$strand == "+", b$tss_loc + b$loc, b$tss_loc - b$loc)
b$genomic_end <- b$genomic_start + 8
write.table(b, file = "high_scoring_TFBS_OC/scan_out/scan.FWD.locations.map", sep = "\t", quote=F, row.names=F, col.names=T)


a <- read.delim("high_scoring_TFBS_OC/scan_out/scan.REV.locations.tbl", col.names = c("pwm", "seq", "loc"), sep="\t")
b <-  extract(a, seq, into = c("gene_id", "chr", "start", "strand"), regex = "(.+?)_(Chr\\d)_(\\d+)_(.)", remove = F)

b$tss_loc <- as.numeric(b$start) + 3001
b$genomic_start <- ifelse(b$strand == "+", b$tss_loc + b$loc, b$tss_loc - b$loc)
b$genomic_end <- b$genomic_start + 8
write.table(b, file = "high_scoring_TFBS_OC/scan_out/scan.REV.locations.map", sep = "\t", quote=F, row.names=F, col.names=T)


#----
#echo software/compute_oc_overlap_highscore_tfbs.sh high_scoring_TFBS_OC/scan_out/scan.FWD.locations.map ibdc_roe-only/oc_peaks_leaf.bed OC_P_LEAF high_scoring_TFBS_OC/scan_out/oc_overlap/scan_fwd_leaf.txt high_scoring_TFBS_OC/scan_out/oc_overlap/temp/ 100000 40 
#| SGE_Array -q *megraw* -m 300G -r high_scoring_TFBS_OC/scan_out/oc_overlap/j1_leaf