#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 7) {
	print("This script is meant to join 4 very specifically-formatted 'feature' files into a single file for modeling and analysis.")
	print("It uses maybe like... 25G of RAM.")
	print("")
	print("Usage: features.rdat tiled_features.rdat features_long_oc_p_leaf.txt features_long_oc_p_root.txt tiled_features_long_oc_p_leaf.txt tiled_features_long_oc_p_root.txt  features_long_oc_p_overall_leaf.txt features_long_oc_p_overall_root.txt  diff_exp_results.txt output.txt")
	print("")
	quit()
}

tile_feature_wide_file <- args[1]
oc_leaf_overall1_file <- args[2]
#oc_leaf_overall2_file <- args[3]
oc_root_overall1_file <- args[3]
#oc_root_overall2_file <- args[5]
oc_leaf_tile_file <- args[4]
oc_root_tile_file <- args[5]
diff_file <- args[6]
diffs_output_file <- args[7]
all_TSS_diff_outfile <- args[8]


############# 
## Read and merge features
#############

tiled_features <- as.data.frame(fread(tile_feature_wide_file, header = TRUE, stringsAsFactors = FALSE))
tiled_features$tss_name <- tiled_features$V1
tiled_features$V1 <- NULL
tiled_features <- tiled_features[order(tiled_features$tss_name),]
print(dim(tiled_features))

load(oc_leaf_overall1_file)
oc_features_overall1_leaf <- all_features_wide
oc_features_overall1_leaf <- oc_features_overall1_leaf[order(oc_features_overall1_leaf$tss_name),]
oc_features_overall1_leaf$tss_name <- NULL
print(dim(oc_features_overall1_leaf))

#load(oc_leaf_overall2_file)
#oc_features_overall2_leaf <- all_features_wide
#oc_features_overall2_leaf <- oc_features_overall2_leaf[order(oc_features_overall2_leaf$tss_name),]
#print(oc_features_overall2_leaf[1:3,])
#oc_features_overall2_leaf$tss_name <- NULL

load(oc_root_overall1_file)
oc_features_overall1_root <- all_features_wide
oc_features_overall1_root <- oc_features_overall1_root[order(oc_features_overall1_root$tss_name),]
oc_features_overall1_root[1:3,]
oc_features_overall1_root$tss_name <- NULL
print(dim(oc_features_overall1_root))

#load(oc_root_overall2_file)
#oc_features_overall2_root <- all_features_wide
#oc_features_overall2_root <- oc_features_overall2_root[order(oc_features_overall2_root$tss_name),]
#oc_features_overall2_root[1:3,]
#oc_features_overall2_root$tss_name <- NULL

load(oc_leaf_tile_file)
oc_features_leaf_tile <- all_features_wide
oc_features_leaf_tile <- oc_features_leaf_tile[order(oc_features_leaf_tile$tss_name),]
oc_features_leaf_tile$tss_name <- NULL
print(dim(oc_features_leaf_tile))

load(oc_root_tile_file)
oc_features_root_tile <- all_features_wide
oc_features_root_tile <- oc_features_root_tile[order(oc_features_root_tile$tss_name),]
oc_features_root_tile$tss_name <- NULL
print(dim(oc_features_root_tile))


#all_features <- cbind(tiled_features, oc_features_overall1_leaf, oc_features_overall2_leaf, oc_features_overall1_root, oc_features_overall2_root, oc_features_leaf_tile, oc_features_root_tile)
all_features <- cbind(tiled_features, oc_features_overall1_leaf, oc_features_overall1_root, oc_features_leaf_tile, oc_features_root_tile)
save(all_features, file = "all_features_temp.R")
print("all_features_wide is ready" )
print(nrow(all_features))
print(ncol(all_features))

print("NA coordinates (if any):")
print(which(is.na(all_features), arr.ind = TRUE))



## cleanup a bit to save RAM
rm(oc_features_root_tile)
rm(oc_features_leaf_tile)
rm(tiled_features)
rm(oc_features_overall1_leaf)
#rm(oc_features_overall2_leaf)
rm(oc_features_overall1_root)
#rm(oc_features_overall2_root)
gc()

#####################################
### Make the diff table nice
#####################################

diff <- read.table(diff_file, header = TRUE, stringsAsFactors = FALSE)

# Mitra: lots of NAs in df > omit them
diff <- na.omit(diff)
diff <- diff[(abs(diff$b) > 1),]


# Mitra: Now we have less infomation to dig into => runtime cost will be reduced dramtically
diff$gene_id <- diff$Accession
diff$Accession <- NULL

#####################################
### merge to big-ass table
#####################################

all_features <- extract(all_features, remove = FALSE,
                             tss_name, c("gene_id", "chr", "loc", "strand", "offset?"), 
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)")


all_features_diffs_wide <- merge(diff, all_features, by = "gene_id")  # about 4500 genes had tss peaks but no diff information; lots (18056) also had diffs info but no peaks. This is an inner join.
print("Done! We want to save them now!")
#rm(all_features)
#rm(diff)
#gc()
print("NA 2 coordinates (if any):")
print(which(is.na(all_features_diffs_wide), arr.ind = TRUE))
all_features_diffs_wide[is.na(all_features_diffs_wide)] <- 0


save(all_features_diffs_wide, file = diffs_output_file)

print("diff features are saved!")
### save all TSS info
diff_all <- read.table(diff_file, header = TRUE, stringsAsFactors = FALSE)
diff_all$gene_id <- diff_all$Accession
diff_all$Accession <- NULL
all_tss_diffs_wide <- merge(diff_all, all_features, by = "gene_id")
save(all_tss_diffs_wide, file = all_TSS_diff_outfile)




