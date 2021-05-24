indir = "~/Downloads/ibdc/aug2019/"
load(paste(indir, "all_features_diffs_wide_withAT.rdat", sep = ""))

which(is.na(all_features_diffs_wide), arr.ind = T)

diffs_colnames <- c("gene_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs",
                    "tech_var", "sigma_sq", "smooth_sigma_sq", "final_sigma_sq",
                    "tss_name", "chr", "loc", "strand", "offset?", "")
d <- all_features_diffs_wide[, !colnames(all_features_diffs_wide) %in% diffs_colnames]
mean_d <- as.data.frame(apply(d, MARGIN = 2, FUN = mean))
mean_d$feature <- rownames(mean_d)
colnames(mean_d)[1] <- "mean_value"

## closed OC
rownames(all_features_diffs_wide) <- all_features_diffs_wide$tss_name
oc_p <- all_features_diffs_wide[,grepl("OC_", colnames(all_features_diffs_wide))]
oc_root <- oc_p[,grepl("ROOT", colnames(oc_p))]
oc_leaf <- oc_p[,grepl("LEAF", colnames(oc_p))]
r = apply(oc_root, 1, sd)
l = apply(oc_leaf, 1, sd)
length(r[r == 0])/nrow(oc_root)
length(l[l == 0])/nrow(oc_leaf)

# what are the expr of closed promoters
all_tss <- paste(indir, "/aligned.peaks.annotated.capped.filtered", sep = "")
df_tss <- read.delim(all_tss, sep = ",", header = T)
head(df_tss)
df_tss$tss_id <- paste(df_tss$TranscriptID, df_tss$Chromosome, df_tss$ModeLocation, df_tss$Strand, "0",  sep = "_")

### tissue-specific
med_high <- all_features_diffs_wide[all_features_diffs_wide$mean_leaf_norm < 300 & all_features_diffs_wide$mean_root_norm < 300,]
all_ts <- all_features_diffs_wide[!all_features_diffs_wide$tss_name %in% med_high$tss_name,]
all_ts <- all_ts[abs(all_ts$b) > 3 & all_ts$qval <= 0.1,]

tss_ts <- df_tss[df_tss$tss_id %in% all_ts$tss_name, ] 
length(unique(tss_ts$TranscriptID))
aggregate(tss_ts$tissue, list(tss_ts$tissue), length)

tss_ts_trx <- aggregate(tss_ts$TranscriptID, list(tss_ts$TranscriptID), length)
nrow(tss_ts_trx[tss_ts_trx$x >= 2,])

root_tss <- df_tss[df_tss$tissue == "root",]
leaf_tss <- df_tss[df_tss$tissue == "leaf",]

length(unique(root_tss$TranscriptID))
length(unique(leaf_tss$TranscriptID))

both <- merge(root_tss, leaf_tss, by = "TranscriptID")
both_trx <- unique(both$TranscriptID)
length(both_trx)
length(unique(both_trx[both_trx %in% tss_ts$TranscriptID]))

not_in_both <- tss_ts[!tss_ts$TranscriptID %in% both_trx,]
length(unique(not_in_both[not_in_both$tissue == "root", "TranscriptID"]))
length(unique(not_in_both[not_in_both$tissue == "leaf", "TranscriptID"]))

aggregate(not_in_both$tissue, list(not_in_both$tissue), length)

oc_closed <- c(r[r != 0], l[l != 0])
closed_tss <- df_tss[df_tss$tss_id %in% names(oc_closed),]
length(unique(closed_tss$TranscriptID))
closed_both <- unique(closed_tss[closed_tss$TranscriptID %in% both_trx, "TranscriptID"])
length(unique(closed_both))
other_closed <- closed_tss[!closed_tss$TranscriptID %in% both$TranscriptID,]
length(unique(other_closed$TranscriptID))
r_closed <- other_closed[other_closed$tissue == "root",]
length(unique(r_closed$TranscriptID))
l_closed <- other_closed[other_closed$tissue == "leaf",]
length(unique(l_closed$TranscriptID))

nrow(l_closed[l_closed$tss_id %in% all_ts$tss_name,])
nrow(r_closed[r_closed$tss_id %in% all_ts$tss_name,])
length(closed_both[closed_both %in% tss_ts_trx$Group.1])

library(tidyr)
# Extract feature types
mean_d <- extract(mean_d, feature, into = c("pwm" ,"strand", "win", "type"), regex = "(.+)_(FWD|REV)_(.)(\\D*)", remove = F)
mean_d$type[mean_d$type == ""] <- "TFBS"
mean_d$type[mean_d$type == "_OC_P_LEAF"] <- "OC_LEAF"
mean_d$type[mean_d$type == "_OC_P_ROOT"] <- "OC_ROOT"

tfbs <- mean_d[mean_d$type == "TFBS",]
hist(tfbs$mean_value, main = "Histogram of TFBS scores", xlab = "mean log-lik score (normalized)")

oc_closed_leaf <- mean_d[mean_d$mean_value < 0.3 & mean_d$type == "OC_LEAF",]
oc_open_leaf <- mean_d[mean_d$mean_value > 0.5 & mean_d$type == "OC_LEAF",]
oc_closed_root <- mean_d[mean_d$mean_value < 0.3 & mean_d$type == "OC_ROOT",]

# close OC and TFBS score
m_leaf <- merge(tfbs, oc_closed_leaf, by = c("pwm", "strand", "win"))
hist(m_leaf$mean_value.x)
high_tfbs_in_closed_leaf <- m_leaf[m_leaf$mean_value.x > 0.7,]

m_root <- merge(tfbs, oc_closed_root, by = c("pwm", "strand", "win"))
hist(m_root$mean_value.x)
high_tfbs_in_closed_root <- m_root[m_root$mean_value.x > 0.7,]
write.table(high_tfbs_in_closed_leaf, file = "~/Downloads/closed_strong_tfbs_sites_leaf.txt", quote = F, row.names = F, sep = "\t")
write.table(high_tfbs_in_closed_root, file = "~/Downloads/closed_strong_tfbs_sites_root.txt", quote = F, row.names = F, sep = "\t")

ggplot(high_tfbs_in_closed_leaf, aes(feature.x, mean_value.y)) + geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Top weighted features and Avg of % openness (LEAF)") + xlab("feature") + ylab("avg %open")
ggplot(high_tfbs_in_closed_root, aes(feature.x, mean_value.y)) + geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Top weighted features and Avg of % openness (ROOT)") + xlab("feature") + ylab("avg %open")

# OC stats for top 50 features
coef_table <- read.table("~/Downloads/ibdc/coefs/top0_coefs.txt", header = F, col.names = c("feature", "coef"))
topX_coef <- coef_table[1:50,]
topX_coef <- extract(topX_coef, feature, into = c("pwm" ,"strand", "win"), regex = "(.+)_(FWD|REV)_(.)\\D*", remove = F)

oc_leaf <- mean_d[mean_d$type == "OC_LEAF",]
oc_root <- mean_d[mean_d$type == "OC_ROOT",]
top_mean_oc_leaf <- merge(topX_coef, oc_leaf, by = c("pwm", "strand", "win"))
top_mean_oc_root <- merge(topX_coef, oc_root, by = c("pwm", "strand", "win"))

hist(top_mean_oc_leaf$mean_value, main = "OC openness for top 50 features (LEAF OC)", xlab = "mean %OC Openness")
hist(top_mean_oc_root$mean_value, main = "OC openness for top 50 features (ROOT OC)", xlab = "mean %OC Openness")

closed_features_leaf <- top_mean_oc_leaf[top_mean_oc_leaf$mean_value < 0.3,]
closed_features_root <- top_mean_oc_root[top_mean_oc_root$mean_value < 0.3,]

# oc only
oc_mean_leaf <- mean_d[grepl("_P_LEAF", mean_d$feature),]
oc_mean_leaf$tissue <- "LEAF"
oc_mean_root <- mean_d[grepl("_P_ROOT", mean_d$feature),]
oc_mean_root$tissue <- "ROOT"
all_oc <- rbind(oc_mean_leaf, oc_mean_root)

library(ggplot2)
ggplot(all_oc, aes(x = mean_value)) + geom_histogram() + facet_wrap(~tissue) +
  ggtitle("Avg %openness in all diff expressed promoters")

library(tidyr)
all_oc <- all_oc[!grepl("OVERALL", all_oc$feature),]
all_oc <- extract(all_oc, feature, into = c("strand", "win"), regex = "\\D+(FWD|REV)_(.)_\\D+", remove = F)
ggplot(all_oc, aes(x = mean_value)) + geom_histogram(aes(fill = tissue)) + facet_grid(tissue ~ win) + xlab("Mean %oc_openness") + 
  ggtitle("Avg %openness in all diff expressed promoters")
