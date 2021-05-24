#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

feature_rdat <- args[1]
all_tss_rdat <- args[2]
filter_pwms_file <- args[3]
outfile_features <-  args[4]
outfile_tss <-  args[5]

load(feature_rdat)
filter_pwms <- read.table(filter_pwms_file, col.names=c("pwm"))
load(all_tss_rdat)

cols_to_exclude <- all_features_diffs_wide[, grep(paste(filter_pwms$pwm, collapse='|'), colnames(all_features_diffs_wide)),]
all_features_diffs_wide <- all_features_diffs_wide[, !(colnames(all_features_diffs_wide) %in% colnames(cols_to_exclude))]


cols_to_exclude <- all_tss_diffs_wide[, grep(paste(filter_pwms$pwm, collapse='|'), colnames(all_tss_diffs_wide)),]
all_tss_diffs_wide <- all_tss_diffs_wide[, !(colnames(all_tss_diffs_wide) %in% colnames(cols_to_exclude))]

save(all_tss_diffs_wide, file = outfile_tss)
save(all_features_diffs_wide, file = outfile_features)
