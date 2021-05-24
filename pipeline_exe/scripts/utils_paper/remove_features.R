#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

feature_rdat <- args[1]
filter_pwms_file <- args[2]
outdir <-  args[3]
outbase <- args[4]

outfile <- paste(outdir, "/all_features_diffs_wide_", outbase ,"_removal.rdat", sep = "")

load(feature_rdat)
filter_pwms <- read.table(filter_pwms_file, col.names=c("pwm"))

cols_to_exclude <- all_features_diffs_wide[, grep(paste(filter_pwms$pwm, collapse='|'), colnames(all_features_diffs_wide)),]
all_features_diffs_wide <- all_features_diffs_wide[, !(colnames(all_features_diffs_wide) %in% colnames(cols_to_exclude))]

save(all_features_diffs_wide, file = outfile)
