#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)
library(data.table)

add_class <- function(input, qval_thresh, fold_thresh) {
  input$class <- NA
  input$class[input$qval <= qval_thresh & input$b < -1 * fold_thresh] <- 1#"leaf_specific"
  input$class[input$qval <= qval_thresh & input$b > fold_thresh] <- 0#"root_specific"
  input <- input[!is.na(input$class), ]
  return(input)
}

args <- commandArgs(trailingOnly = TRUE)
infile = args[1]
diff_file = args[2]

features <- as.data.frame(fread(infile, sep = ",", header = T))
rownames(features) <- features$V1
features$V1 <- NULL
#features <- features[,grepl("(REV_1|REV_3|REV_5|FWD_1|FWD_3|FWD_5)", colnames(features))]
features$tss_name <- rownames(features)

all_features_wide <- features

diff <- read.table(diff_file, header = TRUE, stringsAsFactors = FALSE)
diff_all <- diff

# Mitra: lots of NAs in df > omit them
diff <- na.omit(diff)
diff <- diff[(abs(diff$b) > 1),]
diff$gene_id <- diff$Accession
diff$Accession <- NULL

all_features_wide <- extract(all_features_wide, remove = FALSE,
                             tss_name, c("gene_id", "chr", "loc", "strand"),
                             regex = "([^_]+)_([^_]+)_([^_]+)_([^_]+)")
all_features_wide$"offset?" <- 0

all_features_diffs_wide <- merge(diff, all_features_wide, by = "gene_id")
rem <- all_features_diffs_wide[all_features_diffs_wide$mean_leaf_norm < 300 & all_features_diffs_wide$mean_root_norm < 300,]
all_features_diffs_wide <- all_features_diffs_wide[!all_features_diffs_wide$tss_name %in% rem$tss_name,]

classed_features_diffs_wide <- add_class(all_features_diffs_wide, qval_thresh = 0.1, fold_thres = 3)
df <- classed_features_diffs_wide
rownames(df) <- df$tss_name

cols <- c("gene_id", "chr", "loc", "strand", "offset?", "tss_name", "class")
cols <- c(cols, colnames(diff))
r = df[df$class == 0,]
l = df[df$class == 1,]
l <- l[, !colnames(l) %in% cols]
r <- r[, !colnames(r) %in% cols]
all_r <- r
all_l <- l
all_r$type <- "root"
all_l$type <- "leaf"

r = data.frame(apply(r, 2, mean))
l = data.frame(apply(l, 2, mean))
r$feature <- rownames(r)
l$feature <- rownames(l)

r <- extract(r, feature, into = c("pwm", "loc"), remove = T, regex = "(.+)_TAContent_(.+)")
l <- extract(l, feature, into = c("pwm", "loc"), remove = T, regex = "(.+)_TAContent_(.+)")

colnames(r) <- c("mean_AT", "pwm", "loc")
colnames(l) <- c("mean_AT", "pwm", "loc")

r$type <- "root"
l$type <- "leaf"
df <- rbind(r,l)
all_df_orig <- rbind(all_r, all_l)

library(ggplot2)
df_orig <- df

df <- df_orig
#df <- df[df$strand == "FWD",]
g <- ggplot(df, aes(as.numeric(loc), mean_AT, color = type)) + geom_line(alpha = 0.6) + facet_wrap(~pwm) + ggtitle("AT content")
ggsave(g, file = "same_plots_pwm_atcontent_root_leaf.png")

df$loc <- as.numeric(df$loc)
df <- df[df$loc > -200 & df$loc < 50,]
g <- ggplot(df, aes(as.numeric(loc), mean_AT, color = type)) + geom_line(alpha = 0.6) + facet_wrap(~pwm) + ggtitle("AT Content")
ggsave(g, file = "same_plots_pwm_atcontent_root_leaf_zoomed.png")

#df <- all_df_orig
#df <- df[df$strand == "FWD",]
#df <- gather(df, key = feature, value = LLS, -type)
#df <- extract(df, feature, into = c("pwm", "strand", "loc"), remove = T, regex = "(.+)_(FWD|REV)_\\d_(.+)")
#g <- ggplot(df, aes(as.numeric(loc), LLS, color = type)) + geom_boxplot() + facet_wrap(pwm~type, ncol = 2) + ggtitle("FWD")
#ggsave(g, file = "separate_plots_pwm_lls_root_leaf_fwd_box.png")

#g <- ggplot(df, aes(as.numeric(loc), mean_LLS, color = type)) + geom_line() + facet_wrap(pwm~type, ncol = 2) + ggtitle("FWD")
#ggsave(g, file = "separate_plots_pwm_atcontent_root_leaf_fwd.png")

#df <- df_orig
#df <- df[df$strand == "REV",]
#g <- ggplot(df, aes(as.numeric(loc), mean_LLS, color = type)) + geom_line(alpha = 0.6) + facet_wrap(~pwm) + ggtitle("REV")
#ggsave(g, file = "same_plots_pwm_atcontent_root_leaf_rev.png")

#g <- ggplot(df, aes(as.numeric(loc), mean_LLS, color = type)) + geom_line() + facet_wrap(pwm~type, ncol = 2) + ggtitle("REV")
#ggsave(g, file = "separate_plots_pwm_atcontent_root_leaf_rev.png")

#df$loc <- as.numeric(df$loc)
#df <- df[df$loc > -200 & df$loc < 50,]
#g <- ggplot(df, aes(as.numeric(loc), mean_LLS, color = type)) + geom_line(alpha = 0.6) + facet_wrap(~pwm) + ggtitle("REV")
#ggsave(g, file = "same_plots_pwm_atcontent_root_leaf_rev_zoomed.png")






