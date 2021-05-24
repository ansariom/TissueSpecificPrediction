library(tidyr)
library(plotly)

fc = 3
top_count = 30
nmodels=20

read_files_nochange <- function(f, indir, top_count) {
  all_coefs = read.table(paste(indir, "/",  f, sep = ""), col.names = c("feature", "coef"))
  all_coefs <- all_coefs[1:top_count,]
  all_coefs$rank <- seq(1, top_count)
  return(all_coefs)
}
desktop_mode = TRUE

# input directory contains coef tables from all model runs. 
args = commandArgs(trailingOnly = T)
indir = args[1]
feature_pack_roe = args[2]
feature_pack_tile = args[3]

if (desktop_mode) {
  indir = paste("~/Downloads/ibdc/aug2020/tile_vs_roe/", sep = "")
  feature_pack_tile = paste("~/Downloads/ibdc/aug2020/featureInfo_hardCodedSoftCoded_tile.rdat")
  feature_pack_roe = paste("~/Downloads/ibdc/aug2020/featureInfo_hardCodedSoftCoded.rdat")
}

tile_dir <- paste(indir, "/tiles/", sep = "")
roe_dir = paste(indir, "/roes/", sep = "")

fs = list.files(tile_dir)
all_tile = do.call(rbind, lapply(fs, read_files_nochange, tile_dir, top_count))
load(feature_pack_tile)
final_model_topx <- feature_info[order(-abs(feature_info$coefficient)),][1:top_count,]

all_tile = all_tile[all_tile$feature %in% final_model_topx$feature,]
ggplot(all_tile, aes(x = fct_reorder(feature, rank), y = rank)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90)) + ylab("Rank") +
  xlab("Feature Name") + ggtitle("Variablity of Feature Weight among Top 30 Features in Multiple TEP-Tiled Models") +
  theme(axis.text=element_text(size=12, face = "bold"),  axis.title=element_text(size=14,face="bold"), legend.text =element_text(size=14), title = element_text(size=14,face="bold"))

fs = list.files(roe_dir)
all_roe = do.call(rbind, lapply(fs, read_files_nochange, roe_dir, top_count))
load(feature_pack_roe)
final_model_topx <- feature_info[order(-abs(feature_info$coefficient)),][1:top_count,]
all_roe <- all_roe[all_roe$feature %in% final_model_topx$feature,]
ggplot(all_roe, aes(x = fct_reorder(feature, rank), y = rank)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90)) + ylab("Rank") +
  xlab("Feature Name") + ggtitle("Variablity of Feature Weight among Top 30 Features in Multiple TEP-ROE Models")+
  theme(axis.text=element_text(size=12, face = "bold"),  axis.title=element_text(size=14,face="bold"), legend.text =element_text(size=14), title = element_text(size=14,face="bold"))

