#!/usr/bin/Rscript

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
fc <- as.numeric(args[5])
mean_expr_file <- args[6]
expr_level <- args[7]

# into all_features_diffs_wide
load(input_features)
all <- all_features_diffs_wide
low_expr <- all[all$mean_leaf_norm < 1 & all$mean_root_norm < 1,]
all_features_diffs_wide <- all[!all$tss_name %in% low_expr$tss_name,]

if (expr_level == "med_low") {
	med_low <- all_features_diffs_wide[all_features_diffs_wide$mean_leaf_norm < 30 & all_features_diffs_wide$mean_root_norm < 30,]
	all_features_diffs_wide <- all_features_diffs_wide[!all_features_diffs_wide$tss_name %in% med_low$tss_name,]
}
if (expr_level == "med_high") {
	med_high <- all_features_diffs_wide[all_features_diffs_wide$mean_leaf_norm < 300 & all_features_diffs_wide$mean_root_norm < 300,]
	all_features_diffs_wide <- all_features_diffs_wide[!all_features_diffs_wide$tss_name %in% med_high$tss_name,]
}

if (expr_level == "new_med_low") {
	all_features_diffs_wide <- all_features_diffs_wide[all_features_diffs_wide$mean_leaf_norm > 30 & all_features_diffs_wide$mean_root_norm > 30,]	
}

if (expr_level == "new_med_high") {
        all_features_diffs_wide <- all_features_diffs_wide[(all_features_diffs_wide$mean_leaf_norm > 300 & all_features_diffs_wide$mean_root_norm > 30) | (all_features_diffs_wide$mean_leaf_norm > 30 & all_features_diffs_wide$mean_root_norm > 300),]
}

all_features_diffs_wide$mean_leaf_norm <- NULL
all_features_diffs_wide$mean_root_norm <- NULL


#########
diffs_colnames <- c("gene_id", "pval", "qval", "b", "foldChange", "baseMean",
                              "tss_name", "chr", "loc", "strand", "offset?")

if (model_type == "tile_only") {
  ### : remove all but tiled features
  all_features_diffs_wide <- all_features_diffs_wide[, !(grepl("(FWD|REV)", colnames(all_features_diffs_wide)) & !grepl("tile", colnames(all_features_diffs_wide))) ]
} else if (model_type == "roe_only") {
  ### or, remove tiled features
  all_features_diffs_wide <- all_features_diffs_wide[, !grepl("tile", colnames(all_features_diffs_wide)) ]
} else if (model_type == "oc_only") {
  missing_cols <- all_features_diffs_wide[, colnames(all_features_diffs_wide) %in% diffs_colnames]
  all_features_diffs_wide <- all_features_diffs_wide[, grepl("OC_P_", colnames(all_features_diffs_wide)) ]
  all_features_diffs_wide <- cbind(missing_cols, all_features_diffs_wide)
} else if (model_type == "tfbs_only_tile") {
  all_features_diffs_wide <- all_features_diffs_wide[, !(grepl("(FWD|REV)", colnames(all_features_diffs_wide)) & !grepl("tile", colnames(all_features_diffs_wide))) ]
  all_features_diffs_wide <- all_features_diffs_wide[, !grepl("OC_P_", colnames(all_features_diffs_wide))]
} else if (model_type == "tfbs_only" ) {
  all_features_diffs_wide <- all_features_diffs_wide[, (!grepl("OC_P_", colnames(all_features_diffs_wide)) & !grepl("tile", colnames(all_features_diffs_wide))) ]
} else if (model_type == "roc_only") {
  all_features_diffs_wide <- all_features_diffs_wide[, !grepl("_OVERALL", colnames(all_features_diffs_wide)) ]
} else if ( model_type == "roc_overall") {
  all_features_diffs_wide <- all_features_diffs_wide[, !grepl("_100", colnames(all_features_diffs_wide)) ]
} else if (model_type == "tfbs_ocoverall") {
  all_features_diffs_wide <- all_features_diffs_wide[, !grepl("OC_P_(ROOT|LEAF)",  colnames(all_features_diffs_wide))]
}

# set rownames to tss names
rownames(all_features_diffs_wide) <- all_features_diffs_wide$tss_name

# define classes, get rid of unclassed columns
classed_features_diffs_wide <- add_class(all_features_diffs_wide, qval_thresh = 0.1, fold_thres = fc)
print(dim(classed_features_diffs_wide))
classed_features_diffs_wide <- unique(classed_features_diffs_wide)
print(dim(classed_features_diffs_wide))

# NAs were introduced because many TSSs have overall OC features but not others
## todo: why is this again?
classed_features_diffs_wide <- classed_features_diffs_wide[complete.cases(classed_features_diffs_wide), ]
print("Overall class sizes:")
print(table(classed_features_diffs_wide$class))

# differential expression data
classed_diffs_info <- classed_features_diffs_wide[, diffs_colnames]
# features and class only
classed_features_class <- classed_features_diffs_wide[, !colnames(classed_features_diffs_wide) %in% diffs_colnames]
classed_features_class <- unique(classed_features_class)
print(dim(classed_features_class))


write.table(classed_features_class, file = features_outfile, row.names = T, col.names = T, quote = F, sep = ",")
write.table(classed_diffs_info, file = diff_info_outfile, row.names = T, col.names = T, quote = F, sep = ",")




