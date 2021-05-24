library(ggplot2)
indir = "~/Downloads/ibdc/aug2019/"
load(paste(indir, "/featureInfo_hardCodedSoftCoded.rdat", sep = ""))

## seq content plots for various expression threshold
threshold_all = 1 # all
threshold_low = 30 # med_low
threshold_high = 300 #med_high

fc = c(1,3)
t = c(1, 30, 300)

p1 <- get_ATcontent_mean(threshold_low, "med_low")
p2 <- get_ATcontent_mean(threshold_high, "med_high")
p3 <- get_ATcontent_mean(threshold_all, "all")

seq_feature <- rbind(p1,p2)
ggplot(seq_feature, aes(loc, mean, col = tissue)) + geom_line() + geom_point() + xlab("Promoter Region") + ylab("Mean AT Content") + facet_wrap(~type)

get_ATcontent_mean<- function(threshold, type= "", fc = 3) {
  med_high <- diffs_classes[diffs_classes$mean_leaf_norm < threshold & diffs_classes$mean_root_norm < threshold,]
  mh_features <- diffs_classes[!diffs_classes$tss_name %in% med_high$tss_name,]
  mh_features <- mh_features[abs(mh_features$b) > fc & mh_features$qval <= 0.1,]
  print(dim(mh_features))
  
  root <- mh_features[mh_features$b > 0,]
  shoot <- mh_features[mh_features$b < 0,]
  
  seq_features_root = features[rownames(features) %in% root$tss_name,grepl("content|Content", colnames(features)),]
  seq_features_shoot = features[rownames(features) %in% shoot$tss_name,grepl("content|Content", colnames(features)),]
  
  r <- as.data.frame(apply(seq_features_root, 2, mean))
  r$feature <- rownames(r)
  colnames(r) <- c("mean", "feature")
  r$tissue <- "Root"
  head(r)
  
  l <- as.data.frame(apply(seq_features_shoot, 2, mean))
  l$feature <- rownames(l)
  colnames(l) <- c("mean", "feature")
  l$tissue <- "Shoot"
  head(l)
  
  seq_feature <- rbind(r,l)
  seq_feature <- seq_feature[grepl("TAContent", seq_feature$feature),]
  seq_feature <- tidyr::extract(seq_feature, feature, c("feature_name", "tile"), regex = "(TAContent)(\\d+)")
  seq_feature$tile <- as.numeric(seq_feature$tile)
  seq_feature$loc <- -220 + (seq_feature$tile * 20)
  head(seq_feature)
  #seq_feature$type <- type
  seq_feature$type <- paste(threshold, fc, sep = "_")
  return(seq_feature)
}
