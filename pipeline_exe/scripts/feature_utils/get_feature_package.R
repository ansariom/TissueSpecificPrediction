#!/usr/bin/Rscript
library(ggplot2)
library(dplyr)
library(tidyr)


# adds a "class" column of either "leaf_specific" or "root_specific" based
# on the threshold, assumes the present of column "b" (fold change) and "qval" (Q value)
# this also removes all unclassed examples if strip_unclassed is TRUE (default)
add_class <- function(input, qval_thresh, fold_thresh, strip_unclassed = TRUE) {
  input$class <- -1000
  input$class[input$qval <= qval_thresh & input$b < -1 * fold_thresh] <- 1#"leaf_specific"
  input$class[input$qval <= qval_thresh & input$b > fold_thresh] <- 0#"root_specific"
  if(strip_unclassed) {
    input <- input[!is.na(input$class), ]
  }
  return(input)
}


####################################
##  Start runninh program
####################################
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 3) {
	print("./get_featureInfo_package.R [all_features.rdat] [tile_only|roe_only|tile_roe] [outbase]")
	quit()
}

input_features <- args[1]
tss_tissue_preds <- args[2]
train_tss_ids <- args[3]
test_tss_ids <- args[4]
coef_table <- args[5]
model_type = args[6]
outbase <- args[7]
fc <- as.numeric(args[8])

load(input_features)
all_features_diffs_wide <- all_tss_diffs_wide

#if (model_type == "tile_only") {
  ### : remove all but tiled features
#  all_features_diffs_wide <- all_features_diffs_wide[, !(grepl("(FWD|REV)", colnames(all_features_diffs_wide)) & !grepl("tile", colnames(all_features_diffs_wide))) ]
#} else if (model_type == "roe_only") {
  ### or, remove tiled features
#  all_features_diffs_wide <- all_features_diffs_wide[, !grepl("tile", colnames(all_features_diffs_wide)) ]
#}

rownames(all_features_diffs_wide) <- all_features_diffs_wide$tss_name

# define classes
classed_features_diffs_wide <- add_class(all_features_diffs_wide, qval_thresh = 0.05, fold_thres = fc, strip_unclassed = FALSE)

# NAs were introduced because many TSSs have overall OC features but not others
classed_features_diffs_wide <- classed_features_diffs_wide[complete.cases(classed_features_diffs_wide), ]

print("Overall class sizes:")
print(table(classed_features_diffs_wide$class))


# strip out the differential expression stuff
#diffs_colnames <- c("gene_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs", 
#                    "tech_var", "sigma_sq", "smooth_sigma_sq", "final_sigma_sq", 
#                    "tss_name", "chr", "loc", "strand", "offset?", "class")

diffs_colnames <- c("gene_id", "pval", "qval", "b", "foldChange", "baseMean",
                              "tss_name", "chr", "loc", "strand", "offset?", "class", "mean_leaf_norm", "mean_root_norm")

# differential expression data
diffs_info <- classed_features_diffs_wide[, diffs_colnames]
# features and class only
features <- classed_features_diffs_wide[, !colnames(classed_features_diffs_wide) %in% diffs_colnames]

classes <- classed_features_diffs_wide[, "class", drop = FALSE]


# let's break down the features into various info encoded in them...
################################################
feature_names <- data.frame(feature = colnames(features))
feature_names$type <- "other"
feature_names$type[grepl("FWD|REV", feature_names$feature)] <- "SLL"
feature_names$type[grepl("(OC_P_ROOT)|(OC_P_LEAF)", feature_names$feature)] <- "OC"

oc_features <- feature_names[feature_names$type == "OC", ]
oc_features <- extract(oc_features, feature, c("pwm", "strand", "window", "tissue"), regex = "(.+?)_(FWD|REV)_(\\d+).*_OC_P_(ROOT|LEAF)", remove = FALSE)
oc_features <- oc_features[, c("feature", "type", "pwm", "strand", "window", "tissue")]

sll_features <- feature_names[feature_names$type == "SLL", ]
sll_features <- extract(sll_features, feature, c("pwm", "strand", "window"), regex = "(.+?)_(FWD|REV)_(\\d+).*", remove = FALSE)
sll_features$tissue <- NA
sll_features <- sll_features[, c("feature", "type", "pwm", "strand", "window", "tissue")]

other_features <- feature_names[feature_names$type == "other", ]
other_features$pwm <- NA
other_features$strand <- NA
other_features$window <- NA
other_features$tissue <- NA
other_features <- other_features[, c("feature", "type", "pwm", "strand", "window", "tissue")]

feature_info <- rbind(oc_features, sll_features, other_features)

coef_table <- read.table(coef_table, sep = "\t", col.names = c("feature", "coefficient"))
feature_info <- merge(feature_info, coef_table, by = "feature")

finfo_outfile <- paste(outbase, "/" , "feature_info.tsv", sep = "")
write.table(feature_info, file = finfo_outfile, sep = "\t", row.names = F, col.names = T, quote = F)
print("feature_info saved.")
################################################

merge_by_rownames <- function(df1, df2) {
  df1$asdlkfjsdf <- rownames(df1)
  df2$asdlkfjsdf <- rownames(df2)
  outdf <- merge(df1, df2, by = "asdlkfjsdf", all = T)
  rownames(outdf) <- outdf$asdlkfjsdf
  outdf$asdlkfjsdf <- NULL
  return(outdf)
}

tss_preds <- read.table(tss_tissue_preds, col.names = c("tss_name", "prob0", "prob1"), sep = "\t")
diffs_info <- merge(diffs_info, tss_preds, by = "tss_name")

train_ids <- read.table(train_tss_ids, col.names = "tss_name", sep = "\t")
train_ids$model_usage <- "train"

test_ids <- read.table(test_tss_ids, col.names = "tss_name", sep = "\t")
test_ids$model_usage <- "test"

diffs_info <- merge(diffs_info, train_ids, by = "tss_name", all.x = T)
diffs_info <- merge(diffs_info, test_ids, by = "tss_name", all.x = T)
diffs_info$model_usage.x <- ifelse(diffs_info$model_usage.y == "test", "test", "NA")
diffs_info$model_usage.y <- NULL
diffs_info$model_usage <- diffs_info$model_usage.x
diffs_info$model_usage.x <- NULL

diff_info_outfile <- paste(outbase, "/" , "tss_info.tsv", sep = "")
write.table(diffs_info, file = diff_info_outfile, sep = "\t", row.names = F, col.names = T, quote = F)

