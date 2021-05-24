library(tidyr)

tile_dir <- "~/Downloads/ibdc/Aug2018/tile_vs_roe/tiles/"
roe_dir = "~/Downloads/ibdc/Aug2018/tile_vs_roe/roes/"

all_models_roe <- data.frame(feature=c(), coef=c())
all_models_tile <- data.frame(feature=c(), coef=c())

fs = list.files(tile_dir)
for (f in fs) {
  print(f)
  all_coefs = read.table(paste(tile_dir, "/",  f, sep = ""), col.names = c("feature", "coef"))
  top200 <- all_coefs[1:200,]
  all_models_tile <- rbind(all_models_tile, top200)
}
tile_commons <- aggregate(all_models_tile, by = list(all_models_tile$feature), FUN = length)
tile_commons$coef <- NULL
colnames(tile_commons) <- c("feature", "nobservation")
tile_commons <-  extract(tile_commons, feature, c("pwm", "strand", "win"), regex = "(.+?)_(FWD|REV)_(.+)_tile100", remove = FALSE)
tile_commons$left <- -1000 + (as.numeric(tile_commons$win) - 1) * 100
tile_commons$right <- tile_commons$left + 100
tile_commons$percent_observed <- tile_commons$nobservation/10 * 100
hist(tile_commons$percent_observed, ylim = c(0, nrow(tile_commons)), ylab = "No of features showed up in all models" ,
     xlab = "%observation in 10 models", main = "Histogram of feature observations among top 200 features in 10 models")

library(ggplot2)
data = tile_commons[tile_commons$strand == "FWD", ]
ggplot(data, aes(y = pwm)) + labs(x = "promoter region", y = "PWM")  + 
  geom_segment(aes(x = data$left, y = pwm, xend = data$right, yend = pwm, color = nobservation), size = 1)  + 
  scale_colour_gradient(low = "pink", high = "blue") +
  geom_point(aes(x = data$left, color = nobservation), size = 5, shape = 'I') +
  geom_point(aes(x = data$right, color = nobservation), size = 5, shape = 'I') + 
  scale_x_continuous(breaks=seq(-2000, 1500, 100)) +
  theme_bw() +
  theme(axis.text = element_text(size=7, color="black",face="bold")) + 
  theme(legend.text=element_text(size=9,face="bold")) +
  theme(strip.text = element_text(size=9,face="bold")) +
  ggtitle("Map of Top200 Weighted Features observed in 10 runs of Tiled model (FWD)")


data = tile_commons[tile_commons$strand == "REV", ]
ggplot(data, aes(y = pwm)) + labs(x = "promoter region", y = "PWM")  + 
  geom_segment(aes(x = data$left, y = pwm, xend = data$right, yend = pwm, color = nobservation), size = 1)  + 
  scale_colour_gradient(low = "pink", high = "blue") +
  geom_point(aes(x = data$left, color = nobservation), size = 5, shape = 'I') +
  geom_point(aes(x = data$right, color = nobservation), size = 5, shape = 'I') + 
  scale_x_continuous(breaks=seq(-2000, 1500, 100)) +
  theme_bw() +
  theme(axis.text = element_text(size=7, color="black",face="bold")) + 
  theme(legend.text=element_text(size=9,face="bold")) +
  theme(strip.text = element_text(size=9,face="bold")) +
  ggtitle("Map of Top200 Weighted Features observed in 10 runs of Tiled model (REV)")


nrow(tile_commons[tile_commons$nobservation > 7,])/nrow(tile_commons)
nrow(tile_commons[tile_commons$nobservation < 3 ,])/nrow(tile_commons)
nrow(tile_commons[tile_commons$nobservation > 7,])/nrow(tile_commons)


#------
fs = list.files(roe_dir)
for (f in fs) {
  print(f)
  all_coefs = read.table(paste(roe_dir, "/",  f, sep = ""), col.names = c("feature", "coef"))
  top200 <- all_coefs[1:200,]
  all_models_roe <- rbind(all_models_roe, top200)
}
roe_commons <- aggregate(all_models_roe, by = list(all_models_roe$feature), FUN = length)
 

#----Heatmap for all runs
tile_commons_hm_fwd <- tile_commons[tile_commons$strand == "FWD", colnames(tile_commons) %in% c("pwm", "nobservation", "left")]
tile_commons_hm_fwd <- tile_commons[tile_commons$strand == "REV", colnames(tile_commons) %in% c("pwm", "nobservation", "left")]
f <- list(
  family = "Courier New, monospace",
  size = 1,
  color = "#7f7f7f"
)
a <- list(
  showticklabels = TRUE,
  tickangle = 0,
  exponentformat = "E",
  pad = 1
)
xaxis <- list(
  range = c(-1000,500)
)
m <- list(
  l = 200,
  r = 50,
  b = 50,
  t = 50,
  pad = 1
)

library(plotly)
ggplot(data = tile_commons_hm_fwd, aes(x=left, y=pwm, fill=nobservation)) + 
  geom_tile()

d_fwd <- spread(tile_commons_hm_fwd, key = left, value = nobservation, fill = 0)
mat <- d_fwd[,2:16]
mat <- as.matrix(mat, rownames.force = F)
y <- d_fwd[,1]
x <- colnames(d_fwd)[2:16]

plot_ly( z = mat, type = "heatmap", y = y, x = as.numeric(x), colors = c( "white","#FFEBEE", "#880E4F")) %>%
  layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = "Map of Top200 Weighted Features observed in 10 runs of Tiled model (REV)")

# -------------------------
# Get top50 TFs in all runs, get TFs showing up in those top50 lists, compute avg of the coefs in each tile, show map of those avg values
#-------------------------
all_models_tile <- data.frame(feature=c(), coef=c())
all_coefs_all <- data.frame(feature=c(), coef=c())

indir = tile_dir
indir = roe_dir
fs = list.files(indir)
for (f in fs) {
  print(f)
  all_coefs = read.table(paste(indir, "/",  f, sep = ""), col.names = c("feature", "coef"))
  top200 <- all_coefs[1:50,]
  all_models_tile <- rbind(all_models_tile, top200)
  all_coefs_all <- rbind(all_coefs_all, all_coefs)
}
all_models_tile <-  extract(all_models_tile, feature, c("pwm", "strand", "win", "type"), regex = "(.+?)_(FWD|REV)_(.+)_tile100(.*)", remove = FALSE)
all_coefs_all <-  extract(all_coefs_all, feature, c("pwm", "strand", "win", "type"), regex = "(.+?)_(FWD|REV)_(.+)_tile100(.*)", remove = FALSE)

all_models_roe <-  extract(all_models_tile, feature, c("pwm", "strand", "win", "type"), regex = "(.+?)_(FWD|REV)_(.)(.*)", remove = FALSE)
all_coefs_all <-  extract(all_coefs_all, feature, c("pwm", "strand", "win", "type"), regex = "(.+?)_(FWD|REV)_(.)(.*)", remove = FALSE)

#all_models_tile <- all_models_roe
all_models_tile$fullpwm <- paste(all_models_tile$pwm, all_models_tile$type, sep = "")
tfs_top50 <- unique(all_models_tile$fullpwm)
all_coefs_all$fullpwm <- paste(all_coefs_all$pwm, all_coefs_all$type, sep = "")
coefs_top50_all <- all_coefs_all[all_coefs_all$fullpwm %in% tfs_top50,]

selected_tfs <- coefs_top50_all[, colnames(coefs_top50_all) %in% c("feature", "coef")]
avg_coefs_top50 <- aggregate(selected_tfs, by = list(selected_tfs$feature), FUN = mean)
avg_coefs_top50$feature <- NULL
colnames(avg_coefs_top50) <- c("feature", "coef_avg")
avg_coefs_top50 <- extract(avg_coefs_top50, feature, c("pwm", "strand", "win", "type"), regex = "(.+?)_(FWD|REV)_(.+)_tile100(.*)", remove = FALSE)
avg_coefs_top50$left <- -1000 + (as.numeric(avg_coefs_top50$win) - 1) * 100
avg_coefs_top50$right <- avg_coefs_top50$left + 100
avg_coefs_top50$pwm <- paste(avg_coefs_top50$pwm, avg_coefs_top50$type, sep = "")
avg_coefs_top50_fwd <- avg_coefs_top50[avg_coefs_top50$strand == "FWD", colnames(avg_coefs_top50) %in% c("pwm", "coef_avg", "left")]
avg_coefs_top50_rev <- avg_coefs_top50[avg_coefs_top50$strand == "REV", colnames(avg_coefs_top50) %in% c("pwm", "coef_avg", "left")]

d_fwd <- spread(avg_coefs_top50_fwd, key = left, value = coef_avg, fill = 0)
d_rev <- spread(avg_coefs_top50_rev, key = left, value = coef_avg, fill = 0)

mat <- d_fwd[,2:16]
mat <- as.matrix(mat, rownames.force = F)
y <- d_fwd[,1]
x <- colnames(d_fwd)[2:16]

library(plotly)
plot_ly( z = mat, type = "heatmap", y = y, x = as.numeric(x), colors = c( "#641E16","white", "white", "#1D8348")) %>%
  layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = "Average model weights for top50 TFs in 50 models (FWD)")

mat <- d_rev[,2:16]
mat <- as.matrix(mat, rownames.force = F)
y <- d_rev[,1]
x <- colnames(d_rev)[2:16]
plot_ly( z = mat, type = "heatmap", y = y, x = as.numeric(x), colors = c( "#641E16","white", "white", "#1D8348")) %>%
  layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = "Average model weights for top50 TFs in 50 models (REV)")
