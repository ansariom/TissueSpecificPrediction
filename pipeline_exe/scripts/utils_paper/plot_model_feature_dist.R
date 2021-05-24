infile = "~/Downloads/ibdc/aug2019/all_features_diffs_wide_withAT.rdat"
load(infile)

diffs_colnames <- c("gene_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs",
                    "tech_var", "sigma_sq", "smooth_sigma_sq", "final_sigma_sq",
                    "tss_name", "chr", "loc", "strand", "offset?")
features <- all_features_diffs_wide[, !colnames(all_features_diffs_wide) %in% diffs_colnames]

oc <- features[, grepl("OC_", colnames(features))]
tfbs <- features[, !colnames(features) %in% colnames(oc)]

tfbs_sd <- as.data.frame(apply(tfbs, MARGIN = 2, sd))
colnames(tfbs_sd) <- c("tfbs_sd")
hist(tfbs_sd$tfbs_sd, xlab = "variance", main = "Histogram of feature value's standard deviation \n un-altered TFBS loglik scores")

oc_sd <- as.data.frame(apply(oc, MARGIN = 2, sd))
colnames(oc_sd) <- c("oc_sd")
hist(oc_sd$oc_sd, xlab = "variance", main = "Histogram of OC feature value's standard deviation \n Original OC values")
#------- Max norm
#infile = "~/Downloads/all_features_diffs_wide_TFBSmaxNorm.rdat"
#load(infile)
features <- all_features_diffs_wide[, !colnames(all_features_diffs_wide) %in% diffs_colnames]

oc <- features[, grepl("OC_", colnames(features))]
tfbs <- features[, !colnames(features) %in% colnames(oc)]

tfbs_sd <- as.data.frame(apply(tfbs, MARGIN = 2, sd))
colnames(tfbs_sd) <- c("tfbs_sd")
hist(tfbs_sd$tfbs_sd, xlab = "variance", main = "Histogram of feature value's standard deviation \n TFBS loglik scores Normalized by Observed Maximum ")

#--------- Coefs
plot_coefs <- function(infile, model_type) {
  coef_table = read.table(infile, col.names = c("feature", "Coefficient"))
  oc_coefs <- coef_table[grepl("_OC_", coef_table$feature),]
  oc_coefs$type <- "OC"
  
  tfbs_coefs <- coef_table[!coef_table$feature %in% oc_coefs$feature,]
  tfbs_coefs$type <- "TFBS"
  
  all <- rbind(oc_coefs, tfbs_coefs)
  ggplot(data=all, aes(all$Coefficient)) + geom_histogram() + facet_wrap(~type) + 
    ggtitle(paste("Histogram of Coefficients ( ", model_type, " ) ", sep = ""))
  
}
infile = "~/Downloads/model_standard_7945/roe_only_standard_coef_table.txt"
plot_coefs(infile, "all features scaled by mean and sd")

infile = "~/Downloads/model_zeroOne_75436/roe_only_0-1_coef_table.txt"
plot_coefs(infile, "all features scaled between 0 and 1")

infile = "~/Downloads/model_sd_39520/roe_only_tfbs_sd_0-1_coef_table.txt"
plot_coefs(infile, "all features scaled by sd")

infile = "~/Downloads/model_tfbsSD_ocUnchanged_73621/roe_only_tfbs_sd_oc_unchanged_coef_table.txt"
plot_coefs(infile, "tfbs features scaled by sd, OC Unchanged")

