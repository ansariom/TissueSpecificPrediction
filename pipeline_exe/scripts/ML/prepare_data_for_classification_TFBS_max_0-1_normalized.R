#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)

# adds a "class" column of either "leaf_specific" or "root_specific" based
# on the threshold, assumes the present of column "b" (fold change) and "qval" (Q value)
# this also removes all unclassed examples
add_class <- function(input, qval_thresh, fold_thresh) {
  input$class <- NA
  input$class[input$qval <= qval_thresh & input$b < -1 * fold_thresh] <- 1#"leaf_specific"
  input$class[input$qval <= qval_thresh & input$b > fold_thresh] <- 0#"root_specific"
  input <- input[!is.na(input$class), ]
  return(input)
}

####################################
##  End functions, begin script
####################################
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4) {
  print("./prepare_data_for_classification.R [all_features.rdat] [tile_only|roe_only|tile_roe] [features_outfile] [diffs_outfile]")
  quit()
}

input_features <- args[1]
model_type = args[2]
features_outfile <- args[3]
diff_info_outfile <- args[4]
roe_max_tfbs <- args[5]

# into all_features_diffs_wide
load(input_features)
roe_max <- read.table(roe_max_tfbs, header = T)

# get tfbs features + tss_info
tfbs <- all_features_diffs_wide[, (!grepl("OC_P_", colnames(all_features_diffs_wide)) & !grepl("tile", colnames(all_features_diffs_wide))) ]


# normalize by max in each window
drop <- c("gene_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs",
                    "tech_var", "sigma_sq", "smooth_sigma_sq", "final_sigma_sq",
		    "chr", "loc", "strand", "offset?")
tfbs <- tfbs[, !colnames(tfbs) %in% drop]
tfbs <- gather(tfbs, feature, value, -tss_name)
tfbs <- merge(tfbs, roe_max, by = "feature")
tfbs$normal_val <- tfbs$value/tfbs$loglik_max
tfbs <- tfbs[, c("feature", "tss_name", "normal_val")]
tfbs <- spread(tfbs, feature, normal_val)

#normalize to 0-1
rownames(tfbs) <- tfbs$tss_name
tfbs$tss_name <- NULL
tfbs <- as.data.frame(apply(tfbs, MARGIN = 2, function(x) {x/max(x)}))
tfbs$tss_name <- rownames(tfbs)

# set rownames to tss names
rownames(all_features_diffs_wide) <- all_features_diffs_wide$tss_name
oc <- all_features_diffs_wide[, !colnames(all_features_diffs_wide) %in% colnames(tfbs)]
oc$tss_name <- rownames(oc)

all_features_diffs_wide  <- merge(tfbs, oc, by = "tss_name")

# define classes, get rid of unclassed columns
classed_features_diffs_wide <- add_class(all_features_diffs_wide, qval_thresh = 0.05, fold_thres = 4)

classed_features_diffs_wide <- classed_features_diffs_wide[complete.cases(classed_features_diffs_wide), ]
print("Overall class sizes:")
print(table(classed_features_diffs_wide$class))

# strip out the differential expression stuff
diffs_colnames <- c("gene_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs", 
                    "tech_var", "sigma_sq", "smooth_sigma_sq", "final_sigma_sq", 
                    "tss_name", "chr", "loc", "strand", "offset?")

# differential expression data
classed_diffs_info <- classed_features_diffs_wide[, diffs_colnames]
# features and class only
classed_features_class <- classed_features_diffs_wide[, !colnames(classed_features_diffs_wide) %in% diffs_colnames]

write.table(classed_features_class, file = features_outfile, row.names = T, col.names = T, quote = F, sep = ",")
write.table(classed_diffs_info, file = diff_info_outfile, row.names = T, col.names = T, quote = F, sep = ",")




