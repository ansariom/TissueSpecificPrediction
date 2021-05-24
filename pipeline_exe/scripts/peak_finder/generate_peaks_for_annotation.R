#!/usr/bin/Rscript

library(zoo)
## Peak finding
#=============================================
argmax <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

## closest index finding
#=============================================
find_closest_to_max <- function(max_idx, mod_idxes) {
  distance_to_max <- 10000 
  closest <- -1
  for (idx in mod_idxes) { 
    distance = abs(idx - max_idx) 
    if (distance < distance_to_max) { 
      distance_to_max <- distance 
      closest <- idx 
    } 
  } 
  return(closest)
}

## closes bin finding
#=============================================
find_closest_bin <- function(search_idx, bin_indxes, win = 50) {
  distance_to_max <- 10000 
  closest <- ""
  for (idx in bin_indxes) { 
    distance = abs(search_idx - idx) 
    if (distance < distance_to_max && distance < win) { 
      distance_to_max <- distance 
      closest <- idx 
    } 
  } 
  return(closest)
}


args <- commandArgs(TRUE)
peaks_covg_infile <- args[1]
outfile <- args[2]

#covg <- read.table("Downloads/ath_leaf.f10.b10.normal.covg")
covg <- read.table(peaks_covg_infile)
covg$V5 <- NULL
colnames(covg) <- c("Chromosome", "Start", "End", "name", "Strand", "pos", "mapped_count")
covg$name <- paste(covg$name, covg$End, sep = ".")

## find mod peak
w <- 2 # win length for signal mode finder
span <- 0.3 # span size [0,1] (small values = more sensitive)

get_mod_loc <- function(df) {
  l <- df[1,]$End - df[1,]$Start
  #print(l)
  if (l < 10) {
    span = 0.8
  }
  scores <- df$mapped_count
  locs <- seq(df[1,]$Start, df[1,]$End-1)
  if (length(locs) != length(scores)) {return()}
  mod_loc <- locs[ceiling(length(locs)/2)]

  if (l > 4) {
  peaks <- argmax(locs, scores, w=w, span=span)
  mod_idxes <- peaks$i 
  
  max_idx <- which(scores == max(scores))
  if (length(mod_idxes) > 1 ) {
    peak_vals <- as.integer(peaks$y.hat[mod_idxes])
    max_mod_idx <- mod_idxes[peak_vals == max(peak_vals)]
    mod_loc <- locs[max_mod_idx][1]
  
    d <- merge(max_idx, mod_idxes)
    d$distance <- abs(d$x - d$y)
    idx <- which(d$distance == min(d$distance))[1]
    mod_loc <- locs[d[idx,1]]
  } else if (length(mod_idxes == 1)) {
    mod_loc = locs[mod_idxes[1]]
  }
  }
  df$mod_loc <- mod_loc
  df$mod_count <- scores[which(locs == mod_loc)]
  df$mapped_count <- NULL
  
  #print(df[1,"name"])
  # Peak shape
  has_peak <- FALSE
  narrow_peak <- FALSE
  df$Shape <- "WeakPeak"
  q <- quantile(scores)
  peak_to_med <- as.numeric(q[4]/q[3])
  if (is.nan(peak_to_med)) {return()}
  if (peak_to_med > 1.5)
    has_peak <- TRUE
  if (l < 150)
    narrow_peak <- TRUE
  if (has_peak & narrow_peak)
    df$Shape <- "NarrowPeak"
  else if (has_peak & !narrow_peak)
    df$Shape <- "BroadWithPeak"
  unique(df)
}

cumModeCounts <- aggregate(covg$mapped_count, list(covg$name), sum)
print("cum score computed!")
peaks_modeCounts <- do.call(rbind, by(covg[,c("Chromosome","Strand", "name", "mapped_count", "Start", "End")], covg[,"name"], get_mod_loc))
print("mod loc finished")

# compute mode location and readcount at mode
# ----------------------------------------------------------------
peaks_modeCounts <- merge(peaks_modeCounts, cumModeCounts, by.x = "name", by.y = "Group.1")

colnames(peaks_modeCounts)[6] <- "ModeLocation"
colnames(peaks_modeCounts)[7] <- "ModeReadCount"
colnames(peaks_modeCounts)[9] <- "ReadCount"

peaks_modeCounts <- peaks_modeCounts[!duplicated(peaks_modeCounts[,"name"]),]
row.names(peaks_modeCounts) <- peaks_modeCounts$name
peaks_modeCounts$name <- NULL
peaks <- peaks_modeCounts[order(row.names(peaks_modeCounts)),]

final_peaks <- peaks[,c("Chromosome","Strand","Start","End","ReadCount","ModeLocation","ModeReadCount","Shape")]
sorted_final_peaks <- final_peaks[order(final_peaks[,"Chromosome"], as.numeric(final_peaks[,"Start"])),]
write.table(sorted_final_peaks, file = outfile, sep = ",", row.names = F, quote = F)

