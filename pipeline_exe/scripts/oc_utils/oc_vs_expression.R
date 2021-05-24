
library(ggplot2)
indir = "~/Downloads/ibdc/aug2020/"
load(paste(indir,"featureInfo_hardCodedSoftCoded.rdat", sep = ""))
leaf_open_promoter <- paste(indir, "oc", "leaf_promoter_open_regions_3000-3000.txt", sep = "/")
root_open_promoter <- paste(indir, "oc", "root_promoter_open_regions_3000-3000.txt", sep = "/")
col_names = c("tss_id", "oc_id", "chr", "start", "end", "rel_start", "rel_end")
leaf_oc <- read.table(leaf_open_promoter, col.names = col_names)
root_oc <- read.table(root_open_promoter, col.names = col_names)


#just take immediate open regions
#leaf_oc <- leaf_oc[leaf_oc$rel_start > -1000 & leaf_oc$rel_end < 1000,]
#root_oc <- root_oc[root_oc$rel_start > -1000 & root_oc$rel_end < 1000,]


leaf_oc$width <- abs(leaf_oc$end - leaf_oc$start)
root_oc$width <- abs(root_oc$end - root_oc$start)

leaf_oc$tss_id <- paste(leaf_oc$tss_id, "_0", sep = "")
root_oc$tss_id <- paste(root_oc$tss_id, "_0", sep = "")
#oc_root_features <- features[,grepl("_OC_P_ROOT",colnames(features))]
#oc_leaf_features <- features[,grepl("_OC_P_LEAF",colnames(features))]

root_genes <- diffs_classes[diffs_classes$class == 0,]
leaf_genes <- diffs_classes[diffs_classes$class == 1,]

expr_in_root <- root_oc[root_oc$tss_id %in% root_genes$tss_name,]
expr_in_leaf <- leaf_oc[leaf_oc$tss_id %in% leaf_genes$tss_name,]

expr_in_root <- aggregate(expr_in_root$width, list(expr_in_root$tss_id), sum)
expr_in_leaf <- aggregate(expr_in_leaf$width, list(expr_in_leaf$tss_id), sum)

colnames(expr_in_leaf) <- c("tss_name", "oc_width")
colnames(expr_in_root) <- c("tss_name", "oc_width")

el <- merge(expr_in_leaf, leaf_genes, by = "tss_name")
er <- merge(expr_in_root, root_genes, by = "tss_name")
el$Tissue <- "Shoot"
el$percent_open <- (el$oc_width * 100)/6000 
er$percent_open <- (er$oc_width * 100)/6000 
er$Tissue <- "Root"

both <- rbind(el,er)
ggplot(both, aes(mean_leaf_norm, percent_open, shape = Tissue, col = Tissue)) + geom_point(size = 2) + theme_bw() + scale_shape_manual(values=c(16, 17))+
  scale_color_manual(values=c('#E69F00', '#56B4E9')) 
ggplot(both, aes(mean_root_norm, percent_open, shape = Tissue, col = Tissue)) + geom_point(size = 2) + theme_bw() + scale_shape_manual(values=c(16, 17))+
  scale_color_manual(values=c('#E69F00', '#56B4E9')) 

ggplot(el, aes(mean_leaf_norm, percent_open)) + geom_point() + theme_bw() + xlab("Mean Normalized Expression in Shoot") + ggtitle("OC openness vs. Expression value in Shoot Promoters")
ggplot(er, aes(mean_root_norm, percent_open)) + geom_point() + theme_bw() + xlab("Mean Normalized Expression in Root") + ggtitle("OC openness vs. Expression value in Root Promoters")

ggplot(el[el$mean_leaf_norm < 1000,], aes(mean_leaf_norm, percent_open)) + geom_point() + theme_bw() + xlab("Mean Normalized Expression in Shoot") + ggtitle("OC openness vs. Expression value in Shoot Promoters")
ggplot(er[er$mean_root_norm < 1000,], aes(mean_root_norm, percent_open)) + geom_point() + theme_bw() + xlab("Mean Normalized Expression in Root") + ggtitle("OC openness vs. Expression value in Root Promoters")


library(tidyr)
df <- el[, c("mean_leaf_norm", "mean_root_norm", "percent_open", "tss_name")]
ggplot(df, aes(mean_leaf_norm, mean_root_norm, col= percent_open)) + 
  geom_point() + scale_color_gradientn(colours = c("white", "blue", "red")) + ggtitle("%OC openness and Shoot to Root Expression Ratio")
cor(el$percent_open, el$mean_leaf_norm)

df <- er[, c("mean_leaf_norm", "mean_root_norm", "percent_open", "tss_name")]
ggplot(df, aes(mean_leaf_norm, mean_root_norm, col= percent_open)) + 
  geom_point() + scale_color_gradientn(colours = c("white", "blue", "red")) + ggtitle("%OC openness and Shoot to Root Expression Ratio")




both$root_to_shoot_ratio <- ifelse(both$mean_root_norm > both$mean_leaf_norm, both$mean_root_norm / both$mean_leaf_norm, -1 * both$mean_leaf_norm / both$mean_root_norm)
g <- ggplot(both, aes(root_to_shoot_ratio, percent_open, fill=root_to_shoot_ratio, col = root_to_shoot_ratio)) + geom_point(shape=21) + 
  scale_fill_gradientn(colors = c("limegreen", "limegreen", "white", "red")) + theme_bw() +
  scale_color_gradientn(colors = c("darkgreen", "darkred")) 
ggsave(g, file = )

ggplot(both[abs(both$root_to_shoot_ratio) < 4000,], aes(root_to_shoot_ratio, percent_open, fill=root_to_shoot_ratio, col = root_to_shoot_ratio)) + geom_point(shape=21) + 
  scale_fill_gradientn(colors = c("limegreen", "white", "red")) + theme_bw() +
  scale_color_gradientn(colors = c("darkgreen", "darkred")) 
ggplot(both[abs(both$root_to_shoot_ratio) < 1000,], aes(root_to_shoot_ratio, percent_open, fill=root_to_shoot_ratio, col = root_to_shoot_ratio)) + geom_point(shape=21) + 
  scale_fill_gradientn(colors = c("limegreen", "white", "red")) + theme_bw() +
  scale_color_gradientn(colors = c("darkgreen", "darkred")) 
