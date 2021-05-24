
roe_coef <- read.delim("~/Downloads/ibdc/coefs/top0_coefs.txt", col.names = c("feature", "coef"), sep = "\t")
tile_coef <- read.delim("~/Downloads/oc_only_coef_table.txt", col.names = c("feature", "coef"), sep = "\t")

#roe_coef <- read.delim("~/Downloads/ibdc/coefs/tissue-spec_coef_table.txt", col.names = c("feature", "coef"), sep = "\t")
roe_win_coords = read.table("~/Downloads/ibdc/coefs/roe_win_coordinates.txt", col.names = c("feature", "left", "right")) # coordinates of each window in ROE table


library(tidyr)
library(plotly)
roe_coef <- extract(roe_coef, feature, into = c("pwm" ,"strand", "win", "type"), regex = "(.+)_(FWD|REV)_(.)(\\D*)", remove = F)


####
tile_regx <- "(.+?)_(FWD|REV)_(.+)_tile100(.*)"
tile_coef <-  extract(tile_coef, feature, c("pwm", "strand", "win", "type"), regex = tile_regx, remove = FALSE)
tile_coef$left <- -1000 + (as.numeric(tile_coef$win) - 1) * 100
tile_coef$right <- tile_coef$left + 100

roe_coef <- tile_coef
roe_coef$rank <- seq(1, nrow(roe_coef))
####
oc_coefs <- roe_coef[roe_coef$type %in% "_OC_P_LEAF" | roe_coef$type %in% "_OC_P_ROOT",]
oc_coefs <- oc_coefs[order(oc_coefs$rank),]
oc_coefs$rank <- seq(1, nrow(oc_coefs))
oc_coefs$fullPwm <- paste(oc_coefs$pwm, oc_coefs$strand, oc_coefs$win, sep = "_")
oc_coefs$oc_feature <- paste(oc_coefs$pwm, oc_coefs$strand, oc_coefs$type, sep = "_")
oc_coefs = merge(oc_coefs, roe_win_coords, by.x = "fullPwm", by.y = "feature")
oc_coefs$fullPwm <- NULL

oc_coefs <- oc_coefs[oc_coefs$rank < 600,]

data <- oc_coefs[oc_coefs$strand == "FWD",]
ggplot(data, aes(y = oc_feature)) + labs(x = "promoter region", y = "PWM")  + 
  geom_segment(aes(x = data$left, y = oc_feature, xend = data$right, yend = oc_feature, color = type, alpha = abs(coef)), size = 3)  +
  scale_size(guide = "none") +
  scale_alpha_continuous(guide = FALSE) +
  #scale_x_continuous(breaks=seq(-2000, 1500, 100)) +
  xlim(-1000,1500) +
  scale_color_manual(values = c("red", "blue")) + 
  theme_bw() +
  theme(axis.text = element_text(size=7, color="black",face="bold"), axis.text.x = element_text(angle = 45)) + 
  theme(legend.text=element_text(size=9,face="bold")) +
  theme(strip.text = element_text(size=9,face="bold")) +
  ggtitle("Top 600 Highly Weighted OC Features in Tile Model (FWD)")
+
  facet_wrap(~model_type) 


data <- oc_coefs[oc_coefs$strand == "REV",]
ggplot(data, aes(y = oc_feature)) + labs(x = "promoter region", y = "PWM")  + 
  geom_segment(aes(x = data$left, y = oc_feature, xend = data$right, yend = oc_feature, color = type, alpha = abs(coef)), size = 2)  +
  scale_size(guide = "none") +
  scale_alpha_continuous(guide = FALSE) +
  #scale_x_continuous(breaks=seq(-2000, 1500, 100)) +
  xlim(-1000,1500) +
  scale_color_manual(values = c("red", "blue")) + 
  theme_bw() +
  theme(axis.text = element_text(size=7, color="black",face="bold"), axis.text.x = element_text(angle = 45)) + 
  theme(legend.text=element_text(size=9,face="bold")) +
  theme(strip.text = element_text(size=9,face="bold")) +
  ggtitle("Top 300 Highly Weighted OC Features in Tile Model (REV)")
+
  facet_wrap(~model_type) 





tfbs_coefs <- roe_coef[!roe_coef$feature %in% oc_coefs$feature,]
tfbs_coefs <- tfbs_coefs[order(tfbs_coefs$rank),]
tfbs_coefs$rank <- seq(1, nrow(tfbs_coefs))

oc_coefs$feature <- paste(oc_coefs$pwm, oc_coefs$strand, oc_coefs$win, sep = "_")

m <- merge(oc_coefs, tfbs_coefs, by="feature", suffixes = c(".oc", ".tfbs"))

m<- m[m$rank.tfbs < 100,]
m$diff_rank <- abs(m$rank.oc - m$rank.tfbs)
m[m$diff_rank < 5, c(1,5,7,8,9,15)]

plot(m$rank.tfbs, m$rank.oc)
cor(m$rank.x, m$rank.y)
