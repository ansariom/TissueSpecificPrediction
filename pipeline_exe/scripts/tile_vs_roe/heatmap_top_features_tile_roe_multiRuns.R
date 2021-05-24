library(tidyr)
library(plotly)

fc = 3
top_count = 100
nmodels=30

#indir = paste("~/Downloads/ibdc/aug2019/tile_vs_roe/", sep = "")
indir = paste("/quobyte/Megraw_Lab/project_distributions/TissueSpecificPrediction/pipeline_exe/supplementary_materials/tile_vs_roe/", sep = "")
tile_dir <- paste(indir, "/tiles/", sep = "")
#tile_coefs <- read.table(paste(indir, "/tile_only_coef_table.txt", sep = ""), col.names = c("feature", "coef"))
#roe_coefs <- read.table(paste(indir, "/roe_only_coef_table.txt", sep = ""), col.names = c("feature", "coef"))
roe_dir = paste(indir, "/roes/", sep = "")
win_coords_file = paste(indir, "/roe_win_coords.txt", sep = "")
roe_win_coords = read.table(win_coords_file, col.names = c("feature", "left", "right")) # coordinates of each window in ROE table



#df = get_common_pwm_table(tile_coefs, roe_coefs)
# how many pwms have at least one coefficients that fall in top 100
#----------------------------
get_common_pwm_table <- function(tile, roe) {
  tile$rank <- seq(1, nrow(tile), by = 1)
  roe$rank <- seq(1, nrow(roe), by = 1)
  
  tile <- extract(tile, feature, into = c("pwm", "strand", "type"), regex = "(.+)_(FWD|REV)_\\d+_tile100(.*)")
  roe <- extract(roe, feature, into = c("pwm", "strand", "type"), regex = "(.+)_(FWD|REV)_\\d+(.*)")
  df = data.frame(threshold = c(), roe_only_npwm = c(), tile_only_npwm=c(), n_common=c())
  thresholds = c(20, 50, 100, 150, 200, 300, 400, 500, 1000, 2000, 5000)
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
  df$p_common <- (df$n_common * 100/ df$mean_both )
  return(df)
  
}

heatmap_plots <- function(top_count, topX_data, model_type) {
  f <- list(
    family = "Courier New, monospace",
    size = 1,
    color = "#7f7f7f"
  )
  a <- list(
    title = "PWM",
    showticklabels = TRUE,
    tickangle = 0,
    exponentformat = "E",
    pad = 1
  )
  xaxis <- list(
    title = "Promoter region",
    range = c(-1000,500)
  )
  m <- list(
    l = 200,
    r = 50,
    b = 50,
    t = 50,
    pad = 1
  )
  p_title = paste("Average model weights for top", top_count, " TFs in ", nmodels, " ", model_type, " models (FWD)", sep = "")
  p <- prepare_plot_data(topX_data)
  #g_fwd <- plot_ly( z = p$mat_fwd, type = "heatmap", y = p$y_fwd, x = as.numeric(p$x_fwd), colors = c( "#641E16", "white", "white" ,"#1D8348")) %>%
  #  layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = p_title)
  g_fwd <- heatmaply::heatmaply( p$mat_fwd, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "#641E16", high = "#1D8348", midpoint = 0), layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = p_title)) 

  p_title = paste("Average model weights for top", top_count, " TFs in ", nmodels, " ", model_type, " models (REV)", sep = "")
  g_rev <- plot_ly( z = p$mat_rev, type = "heatmap", y = p$y_rev, x = as.numeric(p$x_rev), colors = c( "#641E16", "white", "white" ,"#1D8348")) %>%
    layout(yaxis = a, xaxis = xaxis, showlegend = FALSE, margin = m, title = p_title)
  return(list(g_fwd=g_fwd, g_rev=g_rev))
  
}
######

read_files <- function(indir, top_count, regx) {
  all_models <- data.frame(feature=c(), coef=c())
  all_coefs_all <- data.frame(feature=c(), coef=c())
  a <- c()
  fs = list.files(indir)
  count = length(fs)
  for (f in fs) {
    print(f)
    all_coefs = read.table(paste(indir, "/",  f, sep = ""), col.names = c("feature", "coef"))
    all_coefs$feature <- gsub("_tile100", "", all_coefs$feature)
    head(all_coefs)
    
    oc_features <- all_coefs[grepl("_OC_P", all_coefs$feature),]
    oc_features$type = "OC"
    
    oc_ov <- all_coefs[grep("OVERALL", all_coefs$feature),]
    oc_ov$type <- "OC"
    oc_ov$win <- 10
    oc_ov$pwm <- oc_ov$feature
    oc_ov$strand <- "FWD"
    
    seq_features = all_coefs[grepl("content|Content", all_coefs$feature),]
    seq_features$type = "seq"
    seq_features <- extract(seq_features, feature, c("pwm", "win"), regex = "(.*Content|.*content)(\\d*)", remove = F)
    seq_features$strand <- "FWD"
    seq_features[seq_features$win == "",]$win <- 10

    tfbs_features <- all_coefs[(!all_coefs$feature %in% oc_features$feature) & (!all_coefs$feature %in% seq_features$feature) & (!all_coefs$feature %in% oc_ov$feature),]
    tfbs_features$type = "TFBS"
    
    a1 = rbind(tfbs_features, oc_features)
    a1 = extract(a1, feature, c("pwm", "strand", "win"), regex = "(.+?)_(FWD|REV)_(\\d+).*", remove = F)
    all_coefs <- rbind(a1, seq_features, oc_ov)
    all_coefs <- all_coefs[order(-abs(all_coefs$coef)),]
    
    top_tfs <- all_coefs[1:top_count,]
    a <- c(a,unique(top_tfs$pwm))
    all_models <- rbind(all_models, top_tfs)
    all_coefs_all <- rbind(all_coefs_all, all_coefs)
  }
  p = data.frame(pwm=a)
  pwm_repeats = aggregate(p, by = list(p$pwm), FUN = length)
  colnames(pwm_repeats) <- c("pwm", "nobserved")
  sel_tfs <- pwm_repeats[pwm_repeats$nobserved > (count/2), "pwm"]

  all_models_tmp <- all_models[all_models$pwm %in% sel_tfs,]
  all_models <- all_models[all_models$feature %in% all_models_tmp$feature,]
  
  all_coefs_tmp <- all_coefs[all_coefs_all$pwm %in% sel_tfs,]
  all_coefs_all <- all_coefs_all[all_coefs_all$feature %in% all_coefs_tmp$feature,]
  
  all_models$fullpwm <- paste(all_models$pwm, all_models$type, sep = "_")
  all_coefs_all$fullpwm <- paste(all_coefs_all$pwm, all_coefs_all$type, sep = "_")
  
  return(list("all_models" = all_models, "all_coefs" = all_coefs_all))
}

avg_coef_computation <- function(all_models, all_coefs, regx) {
  tfs_top50 <- unique(all_models$fullpwm)
  coefs_top50_all <- all_coefs[all_coefs$fullpwm %in% tfs_top50,]
  selected_tfs <- coefs_top50_all[, colnames(coefs_top50_all) %in% c("feature", "coef")]
  avg_coefs_top50 <- aggregate(selected_tfs, by = list(selected_tfs$feature), FUN = mean)
  avg_coefs_top50$feature <- NULL
  colnames(avg_coefs_top50) <- c("feature", "coef_avg")
  
  all_coefs$coef <- NULL
  all_coefs <- unique(all_coefs)
  avg_coefs_top50 = merge(all_coefs, avg_coefs_top50, by = "feature")
  return(avg_coefs_top50)
}

prepare_plot_data <- function(avg_coefs_top50) {
  #avg_coefs_top50 <- avg_coefs_top50_roe
  avg_coefs_top50_fwd <- avg_coefs_top50[avg_coefs_top50$strand == "FWD", colnames(avg_coefs_top50) %in% c("pwm", "coef_avg", "left")]
  avg_coefs_top50_rev <- avg_coefs_top50[avg_coefs_top50$strand == "REV", colnames(avg_coefs_top50) %in% c("pwm", "coef_avg", "left")]
  
  d_fwd <- spread(avg_coefs_top50_fwd, key = left, value = coef_avg, fill = 0.0)
  d_rev <- spread(avg_coefs_top50_rev, key = left, value = coef_avg, fill = 0.0)
  
  mat_fwd <- d_fwd[,2:length(d_fwd)]
  mat_fwd <- as.matrix(mat_fwd, rownames.force = F)
  y_fwd <- d_fwd[,1]
  x_fwd <- colnames(d_fwd)[2:length(d_fwd)]
  mat_rev <- d_rev[,2:length(d_rev)]
  mat_rev <- as.matrix(mat_rev, rownames.force = F)
  y_rev <- d_rev[,1]
  x_rev <- colnames(d_rev)[2:length(d_rev)]
  return(list("x_fwd" = x_fwd, "y_fwd" = y_fwd, "mat_fwd" = mat_fwd, "x_rev" = x_rev, "y_rev" = y_rev, "mat_rev" = mat_rev))
}


#--------------------------------------------------------
#------ heatamps for ROE and Tiled 
#--------------------------------------------------------


# ----- Tiled Models ------
l = read_files(tile_dir, top_count, tile_regx)
all_models_tile <-  l$all_models
all_coefs_tile <-  l$all_coefs
avg_coefs_top50_tile = avg_coef_computation(all_models_tile, all_coefs_tile, tile_regx)
avg_coefs_top50_tile$left <- -1000 + (as.numeric(avg_coefs_top50_tile$win) - 1) * 100
avg_coefs_top50_tile$right <- avg_coefs_top50_tile$left + 100

#avg_coefs_top50_tile[grepl("TAContent", avg_coefs_top50_tile$feature),]$left = -200 + (as.numeric(ta$win) - 1) * 20
#avg_coefs_top50_tile[grepl("TAContent", avg_coefs_top50_tile$feature),]$left = avg_coefs_top50_tile[grepl("TAContent", avg_coefs_top50_tile$feature),]$left + 20

l = heatmap_plots(top_count, avg_coefs_top50_tile, "Tiled")
l$g_fwd
l$g_rev

# ----- ROE Models ----
l = read_files(roe_dir, top_count, roe_regx)
all_models_roe <-  l$all_models
all_coefs_roe <-  l$all_coefs
avg_coefs_top50_roe = avg_coef_computation(all_models_roe, all_coefs_roe, roe_regx)
avg_coefs_top50_roe = merge(avg_coefs_top50_roe, roe_win_coords, by.x = "feature")
avg_coefs_top50_roe$fullPwm <- NULL

l = heatmap_plots(top_count, avg_coefs_top50_roe, "ROE")
l$g_fwd
l$g_rev


#--------------------------------------------------------
#------ Copmare Tiled and ROE
#--------------------------------------------------------
avg_coefs_top50_roe$coef_avg <- avg_coefs_top50_roe$coef_avg/max(abs(avg_coefs_top50_roe$coef_avg))
write.table(avg_coefs_top50_tile[order(-abs(avg_coefs_top50_tile$coef_avg)),], file = "~/Downloads/tile_top_coefs_sorted.txt", sep = "\t", quote = F, row.names = F)
avg_coefs_top50_roe$model_type = "ROE"
avg_coefs_top50_tile$model_type = "Tiled"
avg_coefs_top50_tile$coef_avg <- avg_coefs_top50_tile$coef_avg / max(abs(avg_coefs_top50_tile$coef_avg))
both <- rbind(avg_coefs_top50_roe, avg_coefs_top50_tile)
both$pwm <- paste(both$pwm, both$type, sep = "_")
# ggplot option for zoomed in view
data_fwd = both[both$strand == "FWD",]
data_rev = both[both$strand == "REV",]


## plots
#--------------------------------------------------------
data <- data_fwd
ggplot(data, aes(y = pwm)) + labs(x = "promoter region", y = "PWM")  + 
  geom_segment(aes(x = data$left, y = pwm, xend = data$right, yend = pwm, color = model_type, alpha = abs(coef_avg)), size = 4)  +
#  geom_point(aes(x = data$left, color = model_type, alpha = abs(coef_avg)), size = 5, shape = 'I') +
#  geom_point(aes(x = data$right, color = model_type, alpha = abs(coef_avg)), size = 5, shape = 'I') + 
  scale_size(guide = "none") +
  scale_alpha_continuous(guide = FALSE) +
  #scale_x_continuous(breaks=seq(-2000, 1500, 100)) +
  xlim(-1000,1000) +
  scale_color_manual(values = c("red", "blue")) + 
  theme_bw() +
  theme(axis.text = element_text(size=7, color="black",face="bold"), axis.text.x = element_text(angle = 45)) + 
  theme(legend.text=element_text(size=9,face="bold")) +
  theme(strip.text = element_text(size=9,face="bold")) +
  ggtitle(paste("Comparison between Top ", top_count, " Highly Weighted Features in ROE vs. Tile Models (FWD)", sep = ""))
+  facet_wrap(~model_type) 

data <- data_rev
ggplot(data, aes(y = pwm)) + labs(x = "promoter region", y = "PWM")  + 
  geom_segment(aes(x = data$left, y = pwm, xend = data$right, yend = pwm, color = model_type, alpha = abs(coef_avg)), size = 4)  +
  scale_size(guide = "none") +
  scale_alpha_continuous(guide = FALSE) +
  #scale_x_continuous(breaks=seq(-2000, 1500, 100)) +
  xlim(-1000,1000) +
  scale_color_manual(values = c("red", "blue")) + 
  theme_bw() +
  theme(axis.text = element_text(size=7, color="black",face="bold"), axis.text.x = element_text(angle = 45)) + 
  theme(legend.text=element_text(size=9,face="bold")) +
  theme(strip.text = element_text(size=9,face="bold")) +
  ggtitle(paste("Comparison between Top ", top_count, " Highly Weighted Features in ROE vs. Tile Models (REV)", sep = ""))
+
  facet_wrap(~model_type) 

# ----- numerical stats -----
br = both[both$model_type == "ROE",]
bt = both[both$model_type == "Tiled",]
br <- br[order(-abs(br$coef_avg)),]
bt <- bt[order(-abs(bt$coef_avg)),]
br$rank <- seq(1, nrow(br))
bt$rank <- seq(1, nrow(bt))

df = data.frame(threshold = c(), roe_only_npwm = c(), tile_only_npwm=c(), n_common=c())
thresholds = c(20, 50, 100, 150, 200, 300, 400, 500, 1000, 2000, 5000)
for (i in thresholds) {
  tileX = bt[bt$rank < i,]
  roeX = br[br$rank < i,]
  
  roe_npwm_topX = length(unique(roeX$pwm))
  tile_npwm_topX = length(unique(tileX$pwm))
  
  commonX <- merge(roeX, tileX, by = c("pwm", "strand", "type"))
  ncommonX = length(unique(commonX$pwm))
  
  d <- data.frame(threshold = i, roe_only_npwm = roe_npwm_topX, tile_only_npwm = tile_npwm_topX, n_common = ncommonX)
  df = rbind(df, d)
}
df$mean_both <- apply(df[,c(2,3)], 1, mean)
df$p_common <- (df$n_common * 100/ df$mean_both )

## common vs different
length(unique(both$pwm))
common_80run <- merge(avg_coefs_top50_roe, avg_coefs_top50_tile, by = c("pwm", "strand"))
length(unique(common_80run$pwm))
thresh <- 100
common_80run$in_roe <- ifelse(common_80run$right.x - common_80run$left.y < (-1 * thresh), 0, ifelse(common_80run$left.x - common_80run$right.y > thresh, 0, 1))
not_in_roe <- common_80run[common_80run$in_roe == 0,]
