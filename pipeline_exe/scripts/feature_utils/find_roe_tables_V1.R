#!/usr/bin/env Rscript
# input :
# 1- dist directory
# 2- pwm labels file
# 3- output plots directory
# 4- outbase name

args <- commandArgs(trailingOnly = TRUE)
input_dist_dir <- args[1]
outdir <- args[2]
pwm_labels_file <- args[3]
outbase <- args[4]

argmax <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

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

pwms <- read.table(pwm_labels_file, header = T)
strands <- c("FWD", "REV")

# peak detection params
w <- 100
span <- 0.04

buffer_left <- 50
buffer_right <- 50
peak_mode_loc <- -6000

# bi-modality distance
max_neighbor_dist <- 500


for(strand in strands) {
  strand_outdir <- paste(outdir, "/", strand, sep = "")
  dir.create(paste(strand_outdir), , showWarnings = FALSE)
 
  outTable <- matrix(nrow=length(pwms$pwm), ncol=4)
  i <- 1 
  for(pwm in pwms$pwm) {
    #dist_file <- paste(input_dist_dir, "/", pwm, ".", strand, ".dist", sep = "")
    dist_file <- paste(input_dist_dir, "/", strand, "/", pwm, ".dist", sep = "")
    varuse <- ""
    print(dist_file)
    dist_tbl <- read.table(file=dist_file, header = F)
    head(dist_tbl)
    locs <- dist_tbl[,1] 
    scores <- dist_tbl[,2] 
    varuse <- "" 
    maxdist <- length(locs)/2 
    
    ## Use the smoothing method to find peaks 
    peaks <- argmax(locs, scores, w=w, span=span)
    smooth_plotfile <- paste(strand_outdir, "/", pwm, ".smoothed.jpg", sep = "")
    jpeg(file=smooth_plotfile) 
    plot(locs, scores, cex=0.75, col="Gray", main=paste("w = ", w, ", span = ", span, sep="")) 
    lines(locs, peaks$y.hat,  lwd=2) 
    scores.min <- min(scores) 
    sapply(peaks$i, function(i) lines(c(locs[i],locs[i]), c(scores.min, peaks$y.hat[i]),col="Red", lty=2)) 
    points(locs[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, cex=1.25) 
    dev.off()
    
    ##########################
    ## process found peaks
    ##########################
    
    ## 1 - Sort peak mods and take only top3
    # All the peaks found in 50bp windows
    mod_idxes <- peaks$i 
    
    # any peaks are found?
    if (length(mod_idxes) > 0) {
      
      # keep only modes that are within 1.5 kb
      smooth_scores <- peaks$y.hat
      mod_idxes = mod_idxes[which(locs[mod_idxes] > -1500 & locs[mod_idxes] < 1500)]
      
      if (length(mod_idxes > 2)) {
        last_idx = 2
      } else {
        last_idx = length(mod_idxes)
      }
      cutoff <- sort(smooth_scores[mod_idxes], decreasing = T)[last_idx]
      above_avg_mod_idxes <- mod_idxes[which((smooth_scores[mod_idxes] >= cutoff))]
      
      #------------------------------------------------
      ## Find high scores and their relative peaks
      ## Focus on peaks that show more score enrichment
      #------------------------------------------------
      if (length(above_avg_mod_idxes) > 0 ) {
        # which region has more maximuns (above Q3)
        extreme_max_idxs <- which(scores > quantile(density(scores)$x)[[4]] )
        
        # initialize the bins for each peak mod
        bins <- data.frame(idx = c(), count = c(), sum = c())
        for (peak_candidate_idx in above_avg_mod_idxes ) {
          item <- data.frame(idx = peak_candidate_idx, count = 1, sum = 0)
          bins <- rbind(bins, item)
        }
        
        # add scores to each bin that are within 50 bp of mod
        for (j in extreme_max_idxs) {
          closest_peak_idx <- find_closest_bin(j, bins$idx)
          if (closest_peak_idx > 0) {
            bin_idx <- which(bins$idx == closest_peak_idx)
            bins[bin_idx,2] <- bins[bin_idx, "count"] + 1
            bins[bin_idx,3] <- scores[j] + bins[bin_idx, "sum"]
          }
        }
        
        # get bins that have more max in them
        #------------------------------------
        max_bins <- bins[bins$count == max(bins$count),]
        
        if (nrow(max_bins) > 1) {
          max_idxs <- max_bins[max_bins$sum == max(max_bins$sum), "idx"]
          top_idx <-  which(scores == max(scores))
          closest_to_max_idx <- find_closest_to_max(top_idx, max_idxs)
          max_idx <- closest_to_max_idx
          # if even the sum is the same then choose first one
        } else {
          max_idx <- max_bins$idx
        }
        
        peak_mod_index <- max_idx
        leftind <- peak_mod_index
        rightind <- peak_mod_index
        #------------------------------------
        peaks_scores <- smooth_scores[above_avg_mod_idxes]
        # if the next highest peak is close, then set the right index at the second peak
        if (length(above_avg_mod_idxes) == 2) {
          if (abs(above_avg_mod_idxes[1] - above_avg_mod_idxes[2]) < max_neighbor_dist) {
            if (max_idx == above_avg_mod_idxes[1]) {
              peak_mod_index <- above_avg_mod_idxes[1]
              leftind <- peak_mod_index
              rightind <- above_avg_mod_idxes[2]
            } else {
              peak_mod_index <- above_avg_mod_idxes[2]
              leftind <- above_avg_mod_idxes[1]
              rightin <- peak_mod_index
            }
          }
        } 
      } else {
        print("no peaks were found!") 
        varuse = "NIX" 
      } 
    }else {
      print("no peaks were found!") 
      varuse = "NIX" 
    }
    
    #### Use old method to find region but using smoothed values
    if (varuse != "NIX") { 
      maxbgval <- mean(density(smooth_scores[1:leftind])$x)
      # Find Left 
      nBelow <- 0 
      while (nBelow < buffer_left) { 
        if (leftind > length(scores) | leftind < 1) { 
          varuse <- "NIX" 
          break 
        } 
        if (smooth_scores[leftind] < maxbgval) { 
          nBelow <- nBelow + 1 
        } 
        leftind <- leftind - 1 
      } 
      # Find Right
      maxbgval <- mean(density(smooth_scores[rightind:length(smooth_scores)])$x)
      print(maxbgval)
      nBelow <- 0
      while (nBelow < buffer_right) {
        if (rightind > length(scores)) {
          varuse <- "NIX"
          break
        }
        if (smooth_scores[rightind] < maxbgval) { 
          nBelow <- nBelow + 1 
        } 
        rightind <- rightind + 1 
      } 
      if ((leftind == peak_mod_index) || (rightind == peak_mod_index)) { 
        varuse <- "NIX" 
      } 
    } 
    
    ## write out peaks coordinates
    ##------------------------------
    if (varuse == "NIX") {
      outTable[i,1] <- NA
      outTable[i,2] <- NA
      outTable[i,3] <- NA
      outTable[i,4] <- NA
    } else {
      left <- locs[leftind]
      right <- locs[rightind]
      peak_mode_loc <- locs[peak_mod_index]
      halfWidth <- (right - left)/2.0;
      if (halfWidth > 300) {
        halfWidth <- 300
        left <- peak_mode_loc - halfWidth 
        right <- peak_mode_loc + halfWidth 
      } 
      outTable[i,1] <- format(peak_mode_loc, digits=2) 
      outTable[i,2] <- format(halfWidth, digits=2) 
      outTable[i,3] <- format(left, digits=2) 
      outTable[i,4] <- format(right, digits=2) 
    } 
    
    i <- i+1

    # Plot ROEs
    if (varuse != "NIX") {
      subT <- paste("Interval = [", toString(format(left, digits=2)), ",", toString(format(right, digits=2)), "]")
      mainT <- paste(pwm, subT, sep="\n")
      plot_filename <- paste(strand_outdir, "/", pwm, ".jpg", sep="")      
      jpeg(file=plot_filename)
      par(mfrow=c(2,1))
      plot(locs, scores, main=mainT, xlab="locs", ylab="scores")
      points(locs[locs >= left & locs <= right], scores[locs >= left & locs <= right], col='red', pch=1)
      plot(locs[(locs >= left - 20) & (locs <= right + 20)], scores[(locs >= left - 20) & (locs <= right + 20)], main=mainT, xlab="locs", ylab="scores")
      points(locs[locs >= left & locs <= right], scores[locs >= left & locs <= right], col='red', pch=1)
      points(peak_mode_loc, maxbgval, pch=20, col='green')
      arrows(left, maxbgval, right, maxbgval, col='green', length=0.1, code=3)
      results <- dev.off()
    }
  }
      # Print table of values
    dimnames(outTable)[[2]] <- c("MaxPeakLoc", "HalfWidth", "Left", "Right")
    dimnames(outTable)[[1]] <- pwms$pwm
    write.table(outTable, file=paste(outbase, ".", strand, ".table", sep = ""), row.names=T, col.names=T, quote=F, sep="\t")
}
  







