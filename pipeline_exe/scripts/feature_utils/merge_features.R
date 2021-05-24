#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 7) {
	print("This script is meant to join 4 very specifically-formatted 'feature' files into a single file for modeling and analysis.")
	print("It uses maybe like... 25G of RAM.")
	print("")
	print("Usage: features.rdat tiled_features.rdat features_long_oc_p_leaf.txt features_long_oc_p_root.txt tiled_features_long_oc_p_leaf.txt tiled_features_long_oc_p_root.txt  features_long_oc_p_overall_leaf.txt features_long_oc_p_overall_root.txt  diff_exp_results.txt output.txt")
	print("")
	quit()
}

feat_wide_file <- args[1]
oc_leaf_long_file <- args[2]
oc_root_long_file <- args[3]
oc_leaf_overall_file <- args[4]
oc_root_overall_file <- args[5]
diff_file <- args[6]
diffs_output_file <- args[7]
all_TSS_diff_outfile <- args[8]




############# 
## Read and merge features
#############

features <- read.table(feat_wide_file, header = TRUE, stringsAsFactors = FALSE, check.names = F)
features$tss_name <- rownames(features)
features_long <- gather(features, feature, value, -tss_name)
mean(features_long$value)
median(features_long$value)

oc_features_leaf <- read.table(oc_leaf_long_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)
oc_features_root <- read.table(oc_root_long_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)

oc_features_overall_leaf <- read.table(oc_leaf_overall_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)
oc_features_overall_root <- read.table(oc_root_overall_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)

all_features_long <- rbind(features_long, oc_features_leaf, oc_features_root, oc_features_overall_leaf, oc_features_overall_root)

# to here: ~ 15 mins, 15 gigs (before tiled features addition...)
print("Step1")
all_features_wide <- spread(all_features_long, feature, value)
print(which(is.na(all_features_wide), arr.ind = TRUE))
all_features_wide[is.na(all_features_wide)] <- 0
print("all_features_wide is ready" )

## cleanup a bit to save RAM
rm(features)
rm(features_long)
rm(oc_features_leaf)
rm(oc_features_root)
rm(all_features_long)
gc()

#####################################
### Make the diff table nice
#####################################

diff <- read.table(diff_file, header = TRUE, stringsAsFactors = FALSE)

# A copy of complete features
diff_all <- diff

# Mitra: lots of NAs in df > omit them
diff <- na.omit(diff)
#diff <- diff[(abs(diff$b) > 1),]


# Mitra: Now we have less infomation to dig into => runtime cost will be reduced dramtically
diff$gene_id <- diff$Accession
diff$Accession <- NULL

#####################################
### merge to big-ass table
#####################################

all_features_wide <- extract(all_features_wide, remove = FALSE,
                             tss_name, c("gene_id", "chr", "loc", "strand", "offset?"), 
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)")


diff_all$gene_id <- diff_all$Accession
diff_all$Accession <- NULL
all_tss_diffs_wide <- merge(diff_all, all_features_wide, by = "gene_id")
save(all_tss_diffs_wide, file = all_TSS_diff_outfile)

## save for python use
# strip out the differential expression stuff

diff_all <- NULL
all_tss_diffs_wide <- NULL
gc()

all_features_diffs_wide <- merge(diff, all_features_wide, by = "gene_id")  # about 4500 genes had tss peaks but no diff information; lots (18056) also had diffs info but no peaks. This is an inner join.
print("Done! We want to save them now!")
rm(all_features_wide)
rm(diff)
gc()

save(all_features_diffs_wide, file = diffs_output_file)




