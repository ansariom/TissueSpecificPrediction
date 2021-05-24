#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = T)
all_diff_features <- args[1]
feature_info_file <- args[2]
tss_info <- args[3]
outfile <- args[4]

load(all_diff_features)

feature_info <- read.table(feature_info_file, header = T, check.names = F)

diffs_classes <- read.table(tss_info, header = T)
diffs_classes$model_usage <- NULL
diffs_classes <- diffs_classes[diffs_classes$tss_name %in% all_features_diffs_wide$tss_name,]
rownames(diffs_classes) <- diffs_classes$tss_name

#diffs_colnames <- c("gene_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs", 
#                    "tech_var", "sigma_sq", "smooth_sigma_sq", "final_sigma_sq", 
#                    "tss_name", "chr", "loc", "strand", "offset?")

diffs_colnames <- c("gene_id", "pval", "qval", "b", "foldChange", "baseMean",
                              "tss_name", "chr", "loc", "strand", "offset?",  "mean_leaf_norm", "mean_root_norm")

features <- all_features_diffs_wide
rownames(features) <- features$tss_name
features <- features[, !colnames(features) %in% diffs_colnames]

save(feature_info, features, diffs_classes, file = outfile)

