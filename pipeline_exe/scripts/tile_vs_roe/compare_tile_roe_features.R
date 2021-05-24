library(tidyr)
library(plotly)

## What is the enriched region's boundries based? (FWD1 and 7 + REV 1 and 7)
fwd_roe <- read.table("~/Downloads/roe_3000.FWD.table")
rev_roe <- read.table("~/Downloads/roe_3000.REV.table")

# Analyzing tile-roe Model
# ---------------------------------
roe_tile <- read.table("~/Downloads/ibdc_results/tile_roe_4125/tile_roe_coef_table.txt", col.names = c("feature", "coefficient"))
roe_win_coords = read.table("~/Downloads/roe_win_coords.txt", col.names = c("feature", "left", "right")) # coordinates of each window in ROE table

roe_tile <- extract(roe_tile, feature, c("pwm", "strand", "window"), regex = "(.+?)_(FWD|REV)_(\\d+)", remove = F)
roe_tile$rank <- as.numeric(row.names(roe_tile))
top100 <- roe_tile[roe_tile$rank < 50,]

tile <- top100[grepl("tile100", top100$feature),]
roe <- top100[!grepl("tile100", top100$feature),]

roe = merge(roe, roe_win_coords, by = "feature")
tile$left <- -1000 + (as.numeric(tile$window) - 1) * 100
tile$right <- tile$left + 100

select <- c("pwm", "strand", "left" , "right", "rank")
data1 <- tile[, colnames(tile) %in% select]
data2 <- roe[, colnames(roe) %in% select]

data1$model_type <- "tile"
data2$model_type  <- "roe"
data <- rbind(data1, data2)
data <- data[data$strand == "REV",]

ggplot(data, aes(y = pwm)) + labs(x = "promoter region", y = "PWM")  + 
  geom_segment(aes(x = data$left, y = pwm, xend = data$right, yend = pwm, color = model_type), size = 1)  +
  geom_point(aes(x = data$left, color = model_type), size = 5, shape = 'I') +
  geom_point(aes(x = data$right, color = model_type), size = 5, shape = 'I') + 
  scale_x_continuous(breaks=seq(-1000, 500, 100)) +
  theme_bw() +
  theme(axis.text = element_text(size=9, color="black",face="bold")) + 
  theme(legend.text=element_text(size=9,face="bold")) +
  theme(strip.text = element_text(size=9,face="bold")) +
  ggtitle("roe-tile model Top 50 coefs (REV)") 


# ----------------
roe_tile_wincoords = read.table("~/Downloads/roe_tile_coords.txt")
colnames(roe_tile_wincoords) = c("feature", "start_tile", "end_tile")

roe = read.table("~/Downloads/ibdc/Aug2018/tile_vs_roe/roe_only_coef_table.txt", col.names = c("feature", "coefficient"))
tile = read.table("~/Downloads/ibdc/Aug2018/tile_vs_roe/tile_only_coef_table.txt", col.names = c("feature", "coefficient"))
#tile = read.table("~/Downloads/ibdc_results/tile_only_2270/tile_only_coef_table.txt", col.names = c("feature", "coefficient"))

# get tile coord for ROE wins


roe$type <- "other"
roe$type[grepl("FWD|REV", roe$feature)] <- "SLL"
roe$type[grepl("(OC_P_ROOT)|(OC_P_LEAF)", roe$feature)] <- "OC"
roe$rank = as.numeric(rownames(roe))
roe$ftype = "roe_only"

tile$type <- "other"
tile$type[grepl("FWD|REV", tile$feature)] <- "SLL"
tile$type[grepl("(OC_P_ROOT)|(OC_P_LEAF)", tile$feature)] <- "OC"
tile$rank = as.numeric(rownames(tile))
tile$ftype = "tile_only"


tile = extract(tile, feature, c("pwm", "strand", "window"), regex = "(.+?)_(FWD|REV)_(\\d+)_tile100", remove = FALSE)
tile$left <- -1000 + (as.numeric(tile$window) - 1) * 100
tile$right <- tile$left + 100

roe = extract(roe, feature, c("pwm", "strand", "window"), regex = "(.+?)_(FWD|REV)_(.)", remove = FALSE)
roe = merge(roe, roe_win_coords, by = "feature")


# how many pwms have at least one coefficients that fall in top 100
#----------------------------
df = data.frame(threshold = c(), roe_only_npwm = c(), tile_only_npwm=c(), n_common=c())
thresholds = c(50, 100, 200, 300, 400, 500, 1000, 2000, 5000)
for (i in thresholds) {
  tileX = tile[tile$rank < i,]
  roeX = roe[roe$rank < i,]
  
  roe_npwm_topX = length(unique(roeX$pwm))
  tile_npwm_topX = length(unique(tileX$pwm))
  
  commonX <- merge(roeX, tileX, by = c("pwm", "strand", "type"))
  ncommonX = length(unique(commonX$pwm))
  
  d <- data.frame(threshold = i, roe_only_npwm = roe_npwm_topX, tile_only_npwm = tile_npwm_topX, n_common = ncommonX)
  df = rbind(df, d)
}
df$mean_both <- apply(df[,c(2,3)], 1, mean)
df$p_common <- (df$n_common / df$mean_both ) * 100
df

######
# plot a map of all features
tile100 = tile[tile$rank <100,]
roe100 = roe[roe$rank < 100,]

hist(as.numeric(tile100$window), breaks = 15, xlab = "Tile Window (-1000 to +500 bp relative to TSS)", main = "Histogram of Observation of each tile window\n in Top 100 features")
hist(as.numeric(roe100$window), breaks = 7, xlab = "ROE Sub-window ", main = "Histogram of Observation of each ROE sub-window \nin Top 100 features")

select <- c("pwm", "strand", "left" , "right", "rank")

data1 <- tile100[, colnames(tile100) %in% select]
data2 <- roe100[, colnames(roe100) %in% select]

data1$model_type <- "tile"
data2$model_type  <- "roe"

common <- merge(data1, data2, by = "pwm")
all <- rbind(data1, data2)
data <- all[all$pwm %in% common$pwm,] 
data <- data[order(data$pwm),]

data <- all[all$strand == "FWD",]

 ggplot(data, aes(y = pwm)) + labs(x = "promoter region", y = "PWM")  + 
  geom_segment(aes(x = data$left, y = pwm, xend = data$right, yend = pwm, color = model_type), size = 1)  +
  geom_point(aes(x = data$left, color = model_type), size = 5, shape = 'I') +
  geom_point(aes(x = data$right, color = model_type), size = 5, shape = 'I') + 
  scale_x_continuous(breaks=seq(-2000, 1500, 100)) +
  theme_bw() +
  theme(axis.text = element_text(size=7, color="black",face="bold")) + 
  theme(legend.text=element_text(size=9,face="bold")) +
  theme(strip.text = element_text(size=9,face="bold")) +
  ggtitle("roe vs. tile model Top 100 coefs (FWD)")

 data <- all[all$strand == "REV",]
 ggplot(data, aes(y = pwm)) + labs(x = "promoter region", y = "PWM")  + 
   geom_segment(aes(x = data$left, y = pwm, xend = data$right, yend = pwm, color = model_type), size = 1)  +
   geom_point(aes(x = data$left, color = model_type), size = 5, shape = 'I') +
   geom_point(aes(x = data$right, color = model_type), size = 5, shape = 'I') + 
   scale_x_continuous(breaks=seq(-2000, 1500, 100)) +
   theme_bw() +
   theme(axis.text = element_text(size=7, color="black",face="bold")) + 
   theme(legend.text=element_text(size=9,face="bold")) +
   theme(strip.text = element_text(size=9,face="bold")) +
   ggtitle("roe vs. tile model Top 100 coefs (REV)")

ggsave(p, file = "~/Downloads/both_mode_top100.png", height = 30, width = 30)

#### Compare tile features with enriched regions
tile_coords <- data1
fwd_roe$pwm <- rownames(fwd_roe)
fwd_roe$strand <- "FWD"
rev_roe$pwm <- rownames(rev_roe)
rev_roe$strand <- "REV"
roe_coords <- rbind(fwd_roe, rev_roe)

thresh <- 300
tile_roe_overlaps <- merge(tile_coords, roe_coords, by = c("pwm", "strand"))
tile_roe_overlaps$in_roe <- ifelse(tile_roe_overlaps$right - tile_roe_overlaps$Left < (-1 * thresh), 0, ifelse(tile_roe_overlaps$left - tile_roe_overlaps$Right > thresh, 0, 1))
not_in_roe <- tile_roe_overlaps[tile_roe_overlaps$in_roe == 0,]

nfwd <- not_in_roe[not_in_roe$strand == "FWD",]
nrow(nfwd)
length(unique(nfwd$pwm))

nrev <- not_in_roe[not_in_roe$strand == "REV",]
nrow(nrev)
length(unique(nrev$pwm))

#------ Plot differences
############

tile100 = tile[tile$rank < 100,]
roe100 = roe[roe$rank < 100,]
select <- c("pwm", "strand", "left" , "right", "ftype")

data1 <- tile100[, colnames(tile100) %in% select]
data2 <- roe100[, colnames(roe100) %in% select]
data <- rbind(data1, data2)
data <- na.omit(data)

pwms = unique(data$pwm)

outdir = "~/Downloads/ibdc_results/plots_tile_roe"
dir.create(outdir)

for (pwm in pwms) {
  data1 <- data[data$pwm == pwm,]
  data1 <- na.omit(data1)
  print(pwm)
  p = ggplot(data1, aes(y = ftype)) + labs(x = "promoter region", y = "model type")  + 
    geom_segment(aes(x = data1$left, y = ftype, xend = data1$right, yend = ftype, color = ftype), size = 1)  +
    geom_point(aes(x = data1$left, color = ftype), size = 5, shape = 'I') +
    geom_point(aes(x = data1$right, color = ftype), size = 5, shape = 'I') + 
    facet_wrap(~strand, ncol = 1) + scale_x_continuous(breaks=seq(-1000, 500, 100)) +
    theme_bw() +
    theme(axis.text = element_text(angle = 45,size=9, color="black",face="bold")) + 
    theme(legend.text=element_text(size=9,face="bold")) +
    theme(strip.text = element_text(size=9,face="bold")) + 
    ggtitle(pwm)
  
  fname <- paste(outdir, "/", pwm, ".jpg", sep = "")
  ggsave(fname, p)
}

#########
# plot a heatmap of all features
tile50 = tile[tile$rank <100,]
roe50 = roe[roe$rank < 100,]

select <- c("pwm", "strand", "left" , "right", "rank")

data1 <- tile50[, colnames(tile50) %in% select]
data2 <- roe50[, colnames(roe50) %in% select]

######## plot both in a same heatmap
data1$type <- 1
data2$type <- -1

all <- rbind(data1, data2)
all$label <- paste(all$pwm, all$type, sep = "_")
all$rank <- all$rank * all$type

d_fwd <- all[all$strand == "FWD",]
plot_ly(data = d_fwd, x = as.numeric(d_fwd$left), y = d_fwd$label, z = as.numeric(d_fwd$rank),
        type = "heatmap", colors = c("red", "gray", "white", "gray", "blue")) %>%
  #layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = "Tiled model Top100 (FWD)")
  layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = "Tiled model ALL Features (FWD)")

######## what PWMs in same region?
data1$type <- "tile"
data2$type <- "roe"


########
d_fwd <- data1[data1$strand == "FWD",]
d_rev <- data1[data1$strand == "REV",]

r_fwd <- data2[data2$strand == "FWD",]
r_rev <- data2[data2$strand == "REV",]


#######
f <- list(
  family = "Courier New, monospace",
  size = 1,
  color = "#7f7f7f"
)
a <- list(
  showticklabels = TRUE,
  tickangle = 45,
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

plot_ly(data = d_fwd, x = as.numeric(d_fwd$left), y = d_fwd$pwm, z = as.numeric(d_fwd$rank),
        type = "heatmap", colors = c("red", "gray")) %>%
  #layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = "Tiled model Top100 (FWD)")
  layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = "Tiled model ALL Features (FWD)")

#plotly_IMAGE(p, format = "png", out_file = "~/Downloads/tile_roe_comapre/tile_top100_fwd.png")

plot_ly(data = d_rev, x = as.numeric(d_rev$left), y = d_rev$pwm, z = d_rev$rank, type = "heatmap", colors = c("orange", "white")) %>%
  #layout(yaxis = a,xaxis = xaxis, showlegend = FALSE, margin = m, title = "Tiled model Top100 (REV)")
  layout(yaxis = a,xaxis = xaxis, showlegend = FALSE, margin = m, title = "Tiled model ALL Features (REV)")

plot_ly(data = r_fwd, x = as.numeric(r_fwd$left), y = r_fwd$pwm, z = r_fwd$rank, type = "heatmap", colors = c("blue", "white")) %>%
  layout(yaxis = a, xaxis = xaxis,showlegend = FALSE, margin = m, title = "ROE model ALL Features (FWD)")
  #layout(yaxis = a, xaxis = xaxis,showlegend = FALSE, margin = m, title = "ROE model Top100 (FWD)")

plot_ly(data = r_rev, x = as.numeric(r_rev$left), y = r_rev$pwm, z = r_rev$rank, type = "heatmap", colors = c("darkgreen", "white")) %>%
  layout(yaxis = a, xaxis = xaxis,showlegend = FALSE, margin = m, title = "ROE model ALL Features (REV)")
  #layout(yaxis = a, xaxis = xaxis,showlegend = FALSE, margin = m, title = "ROE model Top100 (REV)")

### Get only top X pwms and train model with it:
pwms <- unique(tile100$pwm)
pwms

names <- c()
for (p in pwms) {
  for (s in c("FWD", "REV")) {
    for (i in seq(1,15)) {
      n <- paste(p, s, i, "tile100", sep = "_")
      names <- c(names, n)
    }
    for (i in seq(1,15)) {
      for (t in c("ROOT", "LEAF")) {
        n <- paste(p, s, i, "tile100", "OC_P", t, sep = "_")
        names <- c(names, n)
      }
    }
  }
}

write.table(names, file = "~/Downloads/pwms_top100_tile.txt", quote = F, row.names = F, col.names = F)

pwms <- unique(roe100$pwm)
pwms

names <- c()
for (p in pwms) {
  for (s in c("FWD", "REV")) {
    for (i in seq(1,15)) {
      n <- paste(p, s, i, sep = "_")
      names <- c(names, n)
    }
    for (i in seq(1,15)) {
      for (t in c("ROOT", "LEAF")) {
        n <- paste(p, s, i, "OC_P", t, sep = "_")
        names <- c(names, n)
      }
    }
  }
}

write.table(names, file = "~/Downloads/pwms_top100_roe.txt", quote = F, row.names = F, col.names = F)

  ####
merged_roe_tilecoords = merge(roe, roe_tile_wincoords, by = "feature")
merged_all <- merge(merged_roe_tilecoords, tile, by = c("pwm", "strand", "type"), suffixes = c("_roe", "_tile"))

pwms = unique(merged_all$pwm)
strands = unique(merged_all$strand)

df = data.frame(pwm=c(), strand=c(), roe_coverd_by_tile = c(), tile_covered_by_roe = c())
for (p in pwms) {
  for (s in strands) {
    f = merged_all[merged_all$pwm == p & merged_all$strand == s & merged_all$type == "SLL",]
    t = sort(unique(f$start_tile))
    r = sort(unique(as.numeric(f$window_tile)))
    total = length(union(t,r)) 
    m = intersect(r,t)
    roe_by_tile = (length(r) - length(setdiff(r,t))) / length(r)
    tile_by_roe = (length(t) - length(setdiff(t,r))) / length(t)
    
    d = data.frame(p,s,roe_by_tile, tile_by_roe)
    df = rbind(df, d)
    print(paste(p, s, roe_by_tile, tile_by_roe), sep = "\t")
  }
}

write.table(df, "roe_tile_model_comparison.txt", quote = F, row.names = F, sep = "\t")
hist(df$roe_by_tile, main = "Frequency of ROEs covered by Tiled regions\n(roe_only model coefs covered by tile_only model coefs)", xlab = "%roe_covered_by_tile")
hist(df$tile_by_roe, main = "Frequency of Tiled covered by ROE regions\n(tile_only model coefs covered by roe_only model coefs)", xlab = "%tile_covered_by_roe")

