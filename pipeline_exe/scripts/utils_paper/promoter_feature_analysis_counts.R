library(tidyr)

type="normByWinL"
type="normByMaxWinL"
type="rawTFBS"
type = "normByMeanSd"
type = "normBy01"

model_indir <- "~/Downloads/ibdc/aug2019/"
path = paste(model_indir, type, sep = "")
path = paste(path, "/after_reduction", sep = "")
load(paste(path, "/featureInfo_hardCodedSoftCoded.rdat", sep = ""))

#### AUROC curve
dim(diffs_classes)
a <- diffs_classes[abs(diffs_classes$b) > 3, ]
df <- a[a$mean_leaf_norm > 300 | a$mean_root_norm > 300,]
head(df)
library(PRROC)
par(mar=c(5.1,4.1,4.1,2.1))
    
fg <- df[df$class == 0, "prob0"]
bg <- df[df$class == 1, "prob0"]

# ROC Curve    
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(roc)

pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(pr)

#####

coefs <- read.table(paste(path, "/roe_only_coef_table.txt", sep = ""),header = F, col.names = c("feature_name", "coef"))

#coefs <- read.table("~/Downloads/ibdc/Oct2018/tile_vs_roe/current_oct2018/3/roe_only_coef_table.txt", header = F, col.names = c("feature_name", "coef"))
maxscores <- read.table("~/Downloads/max_pwm_scores.txt", col.names = c("pwm", "max_score"))

tissue_tfs <- read.table("~/Downloads/tissue_specific_PWMs.txt", header = T)
tissue_tfs <- tissue_tfs[tissue_tfs$mean_leaf_norm >= 30 | tissue_tfs >= 30,]

####################################
# model weight and max PWM score
####################################

cfs <- coefs[!grepl("_OC_", coefs$feature_name),]
cfs$rank <- rownames(cfs)
top_pwms <- unique(extract(cfs, feature_name, regex = "(.+)_(FWD|REV)_\\d", remove = F, into = c("pwm", "strand")))
n <- unique(top_pwms$pwm)
score_df <- data.frame(pwm=n, rank=seq(1:length(n)))
a <- merge(score_df, maxscores, by = "pwm")
c <- cor(a$max_score,a$rank)
outfile = paste(path, "/coefRank_vs_PWMic.png", sep = "")
png(outfile)
plot(a$rank, a$max_score, xlab = "rank of feature", ylab = "max PWM score", main = paste("Relationship between rank of feature vs. IC of PWM \n", type, sep = "")) 
abline(lm(a$max_score~a$rank), col = "red") 
text(90,40, paste("cor_coef = ", c, sep = ""), col = "red")
dev.off()

#### Model weigt and expression level
####################################
transfac <- na.omit(extract(top_pwms, pwm, regex = "(M0[^_]+)_[^_]+_.*\\d$", into = "p"))
other <- top_pwms[!top_pwms$pwm %in% transfac$p,]
a1 <- merge(tissue_tfs, transfac, by.x = "pwm", by.y = "p")
a2 <- merge(other, tissue_tfs, by = "pwm")
tfs <- rbind(a1,a2)
head(tfs)
outfile = paste(path, "/coef_vs_tf_expr.png", sep = "")
png(outfile)
par(mfrow=c(1,2))
plot(tfs$coef, tfs$mean_leaf_norm, xlab = "model_coefficient", ylab = "mean normalized leaf expression", col = "darkgreen", main = type)
plot(tfs$coef, tfs$mean_root_norm, xlab = "model_coefficient", ylab = "mean normalized root expression", col = "red", main = type)
dev.off()
###
### coef plots
#par(new=TRUE)
outfile = paste(path, "/coef_vs_ranks.png", sep = "")
png(outfile)
coefs$rank <- seq(1:nrow(coefs))
plot(coefs$rank, abs(coefs$coef), xlab = "rank of feature", ylab = "coefficient (abs value)", main = type)
dev.off()
plot(coefs[1:1000,]$rank, abs(coefs[1:1000,]$coef), xlab = "rank of feature", ylab = "coefficient (abs value)", main = type)


####################################
# insilico comparision
####################################
types = c( "normByMaxWinL", "normByWinL", "rawTFBS", "normBy01", "normByMeanSd")

all_models <- c()
feature_promoters <- c()
for (type in types) {
  path = paste(model_indir, type, "/after_reduction", sep = "")
  insilico_cands <- read.table(paste(path, "/candidates.txt", sep = ""),header = F, col.names = c("diff_expr", "feature_name", "pre_prob1", "post_prob1","tss", "tissue"))
  insilico_cands$model <- type
  
  print(path)
  load(paste(path, "/featureInfo_hardCodedSoftCoded.rdat", sep = ""))
  insilico_cands <- merge(insilico_cands, diffs_classes, by.x = "tss", by.y = "tss_name")
  
  all_models = rbind(all_models, insilico_cands)  
  g <- ggplot(insilico_cands) + geom_histogram(aes(diff_expr), position = "identity") + ggtitle(paste("model = ", type, sep = ""))
  outfile = paste(path, "/insilico_dist_of_expressionChange.png", sep = "")
  ggsave(g, file = outfile)
  
  npromoters <- aggregate(insilico_cands$feature_name, list(insilico_cands$feature_name), length)
  npromoters$model <- type
  feature_promoters <- rbind(feature_promoters, npromoters)
  
  b <- npromoters
  b <- b[order(-b$x),]
  b$Group.1 <- factor(b$Group.1, levels=unique(as.character(b$Group.1)) )
  g <- ggplot(b, aes((b$Group.1), b$x)) + geom_bar(fill = "#87235f", stat = "identity") + theme(axis.text.x = element_text(angle = 90, size = 6)) +
    xlab("feature Name") + ylab("No. of promoters with diff_expr > 5%") + ggtitle(paste("model = ", type, sep = ""))
  outfile = paste(path, "/insilico_features_expressionChange.png", sep = "")
  ggsave(g, file = outfile)
}

g <- ggplot(all_models) + geom_histogram(aes(diff_expr, fill = model), alpha = 0.6, position = "identity") + ylab("No. of promoters") +
  facet_wrap(~model)  + ggtitle("Number of promoters with >5% change in expression across decision boundry") + theme_bw()
g

miss_classed <- all_models
miss_classed$class <- ifelse(miss_classed$b > 1.9, 0, ifelse(miss_classed$b < -1.9, 1, -1000) )
hard_to_decide <- miss_classed[miss_classed$class == -1000,]
miss_classed <- miss_classed[miss_classed$class != -1000,]
miss_classed$miss_classed <- ifelse(miss_classed$class == 0, ifelse(miss_classed$prob0 < 0.46, 1, -1), -1)
miss_classed$miss_classed <- ifelse(miss_classed$class == 1, ifelse(miss_classed$prob1 < 0.46, 1, miss_classed$miss_classed), miss_classed$miss_classed)
miss_classed <- miss_classed[miss_classed$miss_classed == 1,]
hard_to_decide <- unique(rbind(hard_to_decide,miss_classed[miss_classed$miss_classed == -1,]))

a <- all_models[all_models$model == types[4],]

b <- aggregate(a$feature_name, list(a$feature_name), length)
b <- b[order(-b$x),]
b$Group.1 <- factor(b$Group.1, levels=unique(as.character(b$Group.1)) )
ggplot(b, aes((b$Group.1), b$x)) + geom_bar(fill = "#87235f", stat = "identity") + theme(axis.text.x = element_text(angle = 90, size = 6)) 
b
### fidn number of promoters having high weight TFBS scores in each tissue
top20 <- cfs[1500:1550,]
get_class_counts <- function(feature_name) {
  #print(feature_name)
  df <- data.frame(root_count=c(), leaf_count=c(), noclass=c(), feature_name=c(), cutoff=c())
  #feature_name <- "M1575_1.02_REV_3"
  data <- features[,feature_name]
  h <- hist(data, main = feature_name)
  d <- h$counts
  sum_count <- 0
  for (count in d) {
    #print(count)
    sum_count = sum_count + count
    p <- sum_count / sum(h$counts)
    if (p > 0.95)
      break
  }
  
  
  idx <- which(h$counts == count, arr.ind = T)[1]
  
  threshold <- h$mids[idx]
  print(threshold)
  
  promoters <- rownames(features[features[,feature_name] > threshold,])
  table(cut(diffs_classes[diffs_classes$tss_name %in% promoters, "b"], breaks = c(-10,-2, 1.9,10)))
  counts <- table(diffs_classes[diffs_classes$tss_name %in% promoters, "class"])
  df <- data.frame(root_count=counts["0"], leaf_count=counts["1"], noclass=counts["-1000"], feature_name = feature_name, cutoff=threshold)
  print(df)
  return(df)
}

#par(mfrow=c(4,3))

classes <- table(diffs_classes$class)
all <- do.call(rbind, lapply(top20$feature_name, get_class_counts))
all$diff <- abs(all$root_count - all$leaf_count)
all_1500 <- all
hist(all$diff)
all$root_percent <- 100 * all$root_count / classes["0"]
all$leaf_percent <- 100 * all$leaf_count / classes["1"]
all$noclass_percent <- 100 * all$noclass / classes["-1000"]
write.csv(all, file = "~/Downloads/ibdc/Oct2018/highscore_tfbs_classes_for_top_weighted_features.csv", col.names = T, row.names = F, quote = F)
### ------


# openness
#############
topTFBS <- feature_info[feature_info$type == "SLL",]
hist(topTFBS$coefficient, xlab = "model weight", main = "Histogram of model weights for TFBS features")
topTFBS <- topTFBS[order(-abs(topTFBS$coefficient)),]
topTFBS <- topTFBS[1:20,]
ggplot(topTFBS, aes(factor(window), coefficient)) + geom_boxplot() + xlab("ROE window") + ggtitle("model weights for Top 20 TFBS features and ROE window")
head(topTFBS)

topOC <- feature_info[feature_info$type == "OC",]
hist(topOC$coefficient,  xlab = "model weight", main = "Histogram of model weights for OC features")
topOC <- topOC[order(-abs(topOC$coefficient)),]
topOC <- topOC[1:20,]
ggplot(topOC, aes(factor(window), coefficient)) + geom_boxplot()+ xlab("ROE window") + ggtitle("model weights for Top 20 OC features and ROE window")
head(topOC)


feature_info$feature_name <- paste(feature_info$pwm, feature_info$strand, feature_info$window, sep = "_")
oc_features <- feature_info[feature_info$feature_name %in% topTFBS$feature_name & feature_info$type == "OC",]
# coefs of top weighted features
oc_features$type <- "top_tfbs"
othe_oc_features <- feature_info[!feature_info$feature_name %in% topTFBS$feature_name & feature_info$type == "OC",]
othe_oc_features$type <- "others"
all <- rbind(oc_features,othe_oc_features)
ggplot(all) + geom_histogram(aes(all$coefficient, fill = type), alpha = 0.5, position = "identity") + ggtitle("comparision between the OC weights of top weighted TFBS features and the rest of TFBS features ")
ggplot(all) + geom_density(aes(all$coefficient, fill = type), stat = "density", position = "identity", kernel = "gaussian", alpha = 0.5)+ ggtitle("comparision between the OC weights of top weighted TFBS features and the rest of TFBS features ")
oc_features[abs(oc_features$coefficient) > 0.03,]

# What is the rank of TFBS features related to top weighted OC features
coefs$rank <- seq(1,nrow(coefs))
top_oc <- coefs[grep("_OC_", coefs$feature_name),][1:15,]
hist(top_oc$rank, main = "rank of top weighted oc features", xlab = "rank")
top_oc <- extract(top_oc, feature_name, into = c("fname"), regex = "(.+)_OC_P.*", remove = F)
tfbs_features <- feature_info[feature_info$feature %in% top_oc$fname & feature_info$type == "SLL",]
tfbs_features <- merge(tfbs_features, coefs, by.x = "feature", by.y = "feature_name")
hist(tfbs_features$rank, breaks = 50, main = "rank of TFBS features related to top-weighted OC features", xlab = "rank")

## How open are the top TFBS features?
openness <- features[, colnames(features) %in% oc_features$feature]
leaf_oc <- openness[, grepl("_LEAF", colnames(openness))]
root_oc <- openness[, grepl("_ROOT", colnames(openness))]
leaf_oc <- gather(leaf_oc)
root_oc <- gather(root_oc)
leaf_oc$tissue <- "LEAF"
root_oc$tissue <- "ROOT"

all <- rbind(root_oc, leaf_oc)
all <- extract(all, key, into = "key", regex = "(.+)_OC_P.*", remove = T)
ggplot(all) + geom_boxplot(aes(key, value, fill = tissue, alpha = 0.5) ) + theme(axis.text.x = element_text(angle = 90)) + 
  xlab("Top weighted TFBS feature") + ylab("%OC opennes") + ggtitle("OC openness for ROOT and LEAF top weighted features")
ggplot(all) + geom_violin(aes(key, value, fill = tissue, alpha = 0.5) ) + theme(axis.text.x = element_text(angle = 90)) + 
  xlab("Top weighted TFBS feature") + ylab("%OC opennes") + ggtitle("OC openness for ROOT and LEAF top weighted features")
