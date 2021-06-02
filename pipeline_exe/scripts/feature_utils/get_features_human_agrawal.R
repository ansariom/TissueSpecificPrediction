#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

feat_wide_file <- args[1]
diff_file <- args[2]
diffs_output_file <- args[3]




############# 
## Read and merge features
#############

all_features_wide <- read.table(feat_wide_file, header = TRUE, stringsAsFactors = FALSE, check.names = F)
all_features_wide$tss_name <- rownames(all_features_wide)

diff <- read.table(diff_file, header = TRUE, stringsAsFactors = FALSE)


all_features_wide <- extract(all_features_wide, remove = F,
                             tss_name, c("gene_id", "chr", "loc", "strand", "offset?"), 
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)")

all_features_wide = merge(all_features_wide, diff, by = "gene_id")
row.names(all_features_wide) <- all_features_wide$tss_name
write.csv(all_features_wide[, !colnames(all_features_wide) %in% c("gene_id", "chr", "loc", "strand", "offset?", "tss_name")], file = diffs_output_file, col.names = T, row.names = T, quote = F)




