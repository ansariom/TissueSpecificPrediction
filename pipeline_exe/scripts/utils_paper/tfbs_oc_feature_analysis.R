library(tidyr)
library(ggplot2)
load("~/Downloads/ibdc/aug2020/featureInfo_hardCodedSoftCoded_tile.rdat")
outdir = "~/Downloads/ibdc/aug2020/"
model_type = "Tile"

sorted_f <- feature_info[order(-abs(feature_info$coefficient)),]
sorted_f$rank <- seq(1, nrow(sorted_f))


### Relationship between rank of TFBS and OC features

tfbs <- sorted_f[sorted_f$type == "SLL",]
tfbs$rank <- seq(1, nrow(tfbs))
oc <- sorted_f[sorted_f$type == "OC",]
oc$rank <- seq(1, nrow(oc))

if (model_type == "Tile") {
  tfbs$feature <- gsub("(.+?)_tile100", x = tfbs$feature, replacement = "\\1")
  oc <- extract(oc, feature, c("pwm", "strand", "window", "tissue"), regex = "(.+?)_(FWD|REV)_(\\d+)_tile100_OC_P_(ROOT|LEAF)", remove = FALSE)
}

oc$feature1 <- paste(oc$pwm, oc$strand, oc$window, sep = "_")
m = merge(tfbs, oc, by.x = "feature", by.y = "feature1")


head(oc)
g <- ggplot(m, aes(coefficient.x, coefficient.y)) + geom_point() + xlab("TFBS Feature Weight") + ylab("OC Feature Weight") + 
  ggtitle(paste("correlation between TFBS feature weight and corresponding OC feature weight (", model_type ," Model)", sep = "")) + 
  annotate(geom="text", x=8, y=0.08, label=paste("correlation_coefficient = ",cor(m$coefficient.x, m$coefficient.y), sep = ""), color = "red" )
g
ggsave(g, file = paste(outdir, "/" , model_type, "_tfbs_vs_oc_weights_cor.png", sep = ""))

##### TFBS and OC openness
features[1:10,1:10]
oc_f <- features[,grepl("_OC_P_", colnames(features))]
oc_f$tss_name <- rownames(oc_f)
oc_f <- gather(oc_f, key = "oc_feature", value = "percent_open", -tss_name)

oc_opennes_table <- oc_f
mean_oc_f <- aggregate(oc_opennes_table$percent_open, list(oc_opennes_table$oc_feature), mean)
colnames(mean_oc_f) <- c("feature", "mean_percent_open")
mean_oc_f$f <- gsub("(.+)_(FWD|REV)_(\\d+)_OC_P_\\D+", replacement = "\\1_\\2_\\3", x = mean_oc_f$feature)


tfbs_ocopen <- merge(tfbs, mean_oc_f, by.x = "feature", by.y = "f")
ggplot(tfbs_ocopen, aes(rank, mean_percent_open)) + geom_bo() + xlab("TFBS feature rank") + ylab("Mean OC %open") + ggtitle("correlation between TFBS feature weight and corresponding OC openness") 

oc_f$feature <- gsub("(.+)_(FWD|REV)_(\\d+)_OC_P_\\D+", replacement = "\\1_\\2_\\3", x = oc_f$oc_feature)
to = merge(tfbs[tfbs$rank  < 200,], oc_f, by = "feature")

h = to[to$rank < 50,]
ggplot(h, aes(factor(rank), percent_open)) + geom_boxplot() + xlab("TFBS feature rank") + ylab("OC %open") + ggtitle("correlation between TFBS feature weight and corresponding OC openness") 

  annotate(geom="text", x=3000, y=3000, label=paste("correlation_coefficient = ",cor(m$rank.x, m$rank.y), sep = ""), color = "red" )
