#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
feat_wide_file <- args[1]
seqcont_file <- args[2]
diffs_output_file <- args[3]


load(feat_wide_file)
features <- read.delim(seqcont_file, header = TRUE, stringsAsFactors = FALSE, check.names = F, sep = ",")
colnames(features)[1] <- "tss_name"

print("merging")
all_tss_diffs_wide  = unique(merge(all_tss_diffs_wide, features, by = "tss_name"))


save(all_tss_diffs_wide, file = diffs_output_file)

