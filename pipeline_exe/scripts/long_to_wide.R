#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

oc_leaf_long_file <- args[1]
outfile <- args[2]



############# 
## Read and merge features
#############

all_features_long <- read.table(oc_leaf_long_file, header = FALSE, col.names = c("tss_name", "feature", "value"), stringsAsFactors = FALSE)

# to here: ~ 15 mins, 15 gigs (before tiled features addition...)
all_features_wide <- spread(all_features_long, feature, value)
print("all_features_wide is ready" )

## cleanup a bit to save RAM
rm(all_features_long)
gc()

#####################################
### Make the diff table nice
#####################################

save(all_features_wide, file = outfile)




