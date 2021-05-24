library(tidyr)
library(ggplot2)
top_count = 30
#orig_tile_file = "~/Downloads/ibdc/aug2020/tile_only_coef_table_-1000_500.txt"
orig_tile_file = "tile_only_coef_table_-1000_500.txt"
enhancer_tile_file = "tile_only_coef_table_-2000_500.txt"
#enhancer_tile_file = "~/Downloads/ibdc/aug2020/tile_only_coef_table_-2000_500.txt"
#outdir = "~/Downloads/ibdc/aug2020/"
outdir = "."


#far_model <- read.table(paste(indir, "tile_only_coef_table_-2000_500.txt", sep = ""), col.names = c("feature", "coefficient"))[1:top_count,]

read_file_coef <- function(input_file, tss_win = 10) {
  #input_file = paste(indir, "/",  f, sep = "")
  print(input_file)
  all_coefs = read.table(input_file, col.names = c("feature", "coef"))
  head(all_coefs)
  all_coefs$feature <- gsub("_tile100", "", all_coefs$feature)
  head(all_coefs)
  
  oc_features <- all_coefs[grepl("_OC_P", all_coefs$feature),]
  oc_features$type = "OC"
  
  oc_ov <- all_coefs[grep("OVERALL", all_coefs$feature),]
  oc_ov$type <- "OC"
  oc_ov$win <- tss_win
  oc_ov$pwm <- oc_ov$feature
  oc_ov$strand <- "FWD"
  
  seq_features = all_coefs[grepl("content|Content", all_coefs$feature),]
  seq_features$type = "seq"
  seq_features <- extract(seq_features, feature, c("pwm", "win"), regex = "(.*Content|.*content)(\\d*)", remove = F)
  seq_features$strand <- "FWD"
  seq_features[seq_features$win == "",]$win <- tss_win
  
  tfbs_features <- all_coefs[(!all_coefs$feature %in% oc_features$feature) & (!all_coefs$feature %in% seq_features$feature) & (!all_coefs$feature %in% oc_ov$feature),]
  tfbs_features$type = "TFBS"
  
  a1 = rbind(tfbs_features, oc_features)
  a1 = extract(a1, feature, c("pwm", "strand", "win"), regex = "(.+?)_(FWD|REV)_(\\d+).*", remove = F)
  all_coefs <- rbind(a1, seq_features, oc_ov)
  all_coefs <- all_coefs[order(-abs(all_coefs$coef)),]
  
  top_tfs <- all_coefs[1:top_count,]
  return(top_tfs)
  
}

orig_model <- read_file_coef(orig_tile_file, tss_win = 10)
orig_model$left <- -1000 + (as.numeric(orig_model$win) - 1) * 100
orig_model$right <- orig_model$left + 100

far_model <- read_file_coef(enhancer_tile_file, tss_win = 20)
far_model$left <- -2000 + (as.numeric(far_model$win) - 1) * 100
far_model$right <- far_model$left + 100

orig_model$type <- "original"
far_model$type <- "enhancer"
both <- rbind(orig_model, far_model)
data <- both[both$strand == "FWD",]
g <- ggplot(data, aes(y = pwm)) + labs(x = "promoter region", y = "PWM")  + 
  geom_segment(aes(x = data$left, y = pwm, xend = data$right, yend = pwm, color = type), size = 1)  +
  geom_point(aes(x = data$left, color = type, shape = type), size = 3) + scale_shape_manual(values = c(21, 10)) + 
  geom_point(aes(x = data$right, color = type, shape = type), size = 3)  + 
  scale_size(guide = "none") +
  scale_alpha_continuous(guide = FALSE) +
  #scale_x_continuous(breaks=seq(-2000, 1500, 100)) +
  xlim(-2000,500) +
  scale_color_manual(values = c("red", "blue")) + 
  theme_bw() + 
  theme(axis.text=element_text(size=12, face = "bold"),  axis.title=element_text(size=14,face="bold"), legend.text =element_text(size=14), title = element_text(size=14,face="bold")) +
  theme(strip.text = element_text(size=7,face="bold")) +
  ggtitle(paste("Comparison between Top ", top_count, " Highly Weighted Features in Tile Models (original vs enhancer) (FWD)", sep = ""))
g
ggsave(g, file = paste(outdir, "/top_", top_count, "_tile_models_compare_top_features_fwd.png", sep = ""), width = 10)

data <- both[both$strand == "REV",]
g <- ggplot(data, aes(y = pwm)) + labs(x = "promoter region", y = "PWM")  + 
  geom_segment(aes(x = data$left, y = pwm, xend = data$right, yend = pwm, color = type), size = 1)  +
  geom_point(aes(x = data$left, color = type, shape = type), size = 3) + scale_shape_manual(values = c(21, 10))+
  geom_point(aes(x = data$right, color = type, shape = type), size = 3) + 
  scale_size(guide = "none") +
  scale_alpha_continuous(guide = FALSE) +
  #scale_x_continuous(breaks=seq(-2000, 1500, 100)) +
  xlim(-2000,500) +
  scale_color_manual(values = c("red", "blue")) + 
  theme_bw() + 
  theme(axis.text=element_text(size=12, face = "bold"),  axis.title=element_text(size=14,face="bold"), legend.text =element_text(size=14), title = element_text(size=14,face="bold")) +
  theme(strip.text = element_text(size=7,face="bold")) +
  ggtitle(paste("Comparison between Top ", top_count, " Highly Weighted Features in Tile Models (original vs enhancer) (REV)", sep = ""))
g
ggsave(g, file = paste(outdir, "/top_", top_count, "_tile_models_compare_top_features_rev.png", sep = ""), width = 10, dpi = "retina")
