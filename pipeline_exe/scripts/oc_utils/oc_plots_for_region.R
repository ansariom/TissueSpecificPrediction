# ref: https://stackoverflow.com/questions/21506724/how-to-plot-overlapping-ranges-with-ggplot2
plot_oc <- function(root_oc, leaf_oc, left, right, title) {
  leaf_oc_in_region <- leaf_oc[leaf_oc$rel_start > left & leaf_oc$rel_end < right,]
  root_oc_in_region <- root_oc[root_oc$rel_start > left & root_oc$rel_end < right,]
  
  root_oc_in_region$width <- root_oc_in_region$rel_end - root_oc_in_region$rel_start
  leaf_oc_in_region$width <- leaf_oc_in_region$rel_end - leaf_oc_in_region$rel_start
  
  library(IRanges)   
  ir1 <- IRanges(root_oc_in_region$rel_start, width = root_oc_in_region$width)
  ir2 <- IRanges(leaf_oc_in_region$rel_start, width = leaf_oc_in_region$width)
  
  par(mfrow=c(2,1))
  plot_coverage(ir1, tissue = "ROOT", line_color = "#821e00")
  plot_coverage(ir2, tissue = "SHOOT", line_color = "#006e00")
  
  
}

plot_coverage <- function(ir, tissue = "ROOT", xadj = -3000.5, line_color = "#006e00") {
  ir <- IRanges(start(ir) + 3000, end(ir) + 3000)
  cov = coverage(ir)
  plotRanges(x = ir, main=paste(tissue, " Open Regions", sep = ""))
  cov <- as.vector(cov)
  mat <- cbind(seq_along(cov) + xadj, cov)
  d <- diff(cov) != 0
  mat <- rbind(cbind(mat[d,1]+1, mat[d,2]), mat)
  mat <- mat[order(mat[,1]),]
  lines(mat, col=line_color, lwd=4)
  axis(2)
}

plotRanges <- function(x, main=main_title, col="#ACA2A2", sep=0.5, ...) {
  height <- 1
  x <- IRanges(start(x) - 3000, end(x) - 3000)
  xlim = x
  if (is(xlim, "IntegerRanges")) 
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  print(par("mar"))
  par(mar=c(2,2,2,1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  #rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, border = NA)
  title(main)
  new_x = start(x) - 3000
  axis(1)
}

plot_oc_mitra <- function(ir1, ir2) {
  bins1 <- disjointBins(IRanges(start(ir1), end(ir1) + 1))
  bins2 <- disjointBins(IRanges(start(ir2), end(ir2) + 1))
  dat1 <- cbind(as.data.frame(ir1), bin = bins1)
  dat2 <- cbind(as.data.frame(ir2), bin = bins2)
  dat1$tissue = "ROOT"
  dat2$tissue = "LEAF"
  
  dat <- rbind(dat1, dat2)
  
  library(ggplot2)
  ggplot(dat) + 
    geom_rect(aes(xmin = start, xmax = end,
                  ymin = bin, ymax = bin + 0.8)) + xlab(paste("Promoter Region (", left, " to ", right, ")", sep = "")) + ylab("TSSs") + 
    theme_bw() + 
    facet_wrap(~tissue) +
    ggtitle(paste("Genome-wide Open Chromatin state relative to TSS locations \n", title, sep = ""))
  #ggtitle("OC coverage for leaf/root specific genes (diff expressed promoters only)\n root OC in leaf-specific promoters vs. leaf OC in root-specific promoters")
  #ggtitle("Genome-wide Open Chromatin state relative to TSS locations ")  
  
}

### Side by side view of OC regions
plot_oc_summary <- function(root_oc, leaf_oc, left, right, title) {
  leaf_oc_in_region <- leaf_oc[leaf_oc$rel_start > left & leaf_oc$rel_end < right,]
  root_oc_in_region <- root_oc[root_oc$rel_start > left & root_oc$rel_end < right,]
  root_oc_in_region$tissue <- -40
  leaf_oc_in_region$tissue <- 40
  
  
  all <- rbind(leaf_oc_in_region, root_oc_in_region) 
  bins = transform(all, bin = cut(rel_start, 50, dig.lab = 4, ordered_result = T))
  
  a <- list(
    showticklabels = TRUE,
    tickangle = 0,
    exponentformat = "E",
    pad = 1
  )
  b <- list(
    showticklabels = TRUE,
    tickangle = 45,
    exponentformat = "E",
    pad = 1
  )
  m <- list(
    l = 50,
    r = 50,
    b = 100,
    t = 50,
    pad = 1
  )
  
  bins$id <- as.numeric(factor(bins$tss_id))
  plot_ly(data = bins, x = bins$bin, y = bins$id, z = as.numeric(bins$tissue),
          type = "heatmap", colors = c("blue", "green"), alpha = 0.5) %>%
    layout(yaxis = a, showlegend = FALSE, margin = m, title = title, xaxis = b, xlab("promoter region"))
  
}

plot_oc_both <- function(root_oc, leaf_oc, left, right, title) {
  leaf_oc_in_region <- leaf_oc[leaf_oc$rel_start > left & leaf_oc$rel_end < right,]
  root_oc_in_region <- root_oc[root_oc$rel_start > left & root_oc$rel_end < right,]
  
  root_oc_in_region$width <- root_oc_in_region$rel_end - root_oc_in_region$rel_start
  leaf_oc_in_region$width <- leaf_oc_in_region$rel_end - leaf_oc_in_region$rel_start
  
  root_oc_in_region$tissue <- "ROOT"
  leaf_oc_in_region$tissue <- "LEAF"
  
  all <- rbind(root_oc_in_region, leaf_oc_in_region)
  
  ggplot(all, aes(y = tss_id)) + labs(x = "promoter region", y = "TSS")  + 
    geom_segment(aes(x = all$rel_start, y = tss_id, xend = all$rel_end, yend = tss_id, color = tissue), size = 1)  +
    geom_point(aes(x = all$rel_start, color = tissue), size = 5, shape = 'I') +
    geom_point(aes(x = all$rel_end, color = tissue), size = 5, shape = 'I') + 
   # scale_x_continuous(breaks=seq(-1000, 500, 100)) +
    theme_bw() +
    theme(axis.text = element_text(size=9, color="black",face="bold")) + 
    theme(legend.text=element_text(size=9,face="bold")) +
    theme(strip.text = element_text(size=9,face="bold")) +
    ggtitle("roe-tile model Top 50 coefs (REV)") 
  
#    ggplot(dat) + 
#    geom_rect(aes(xmin = start, xmax = end, 
#                  ymin = bin, ymax = bin + 0.8, color = tissue)) + xlab(paste("Promoter Region (", left, " to ", right, ")", sep = "")) + ylab("TSSs") + 
#    theme_bw() + 
#    ggtitle(paste("Genome-wide Open Chromatin state relative to TSS locations \n", title, sep = ""))
  #ggtitle("OC coverage for leaf/root specific genes (diff expressed promoters only)\n root OC in leaf-specific promoters vs. leaf OC in root-specific promoters")
  #ggtitle("Genome-wide Open Chromatin state relative to TSS locations ")  
  
}

indir = "~/Downloads/ibdc/aug2020/oc/"
leaf_open_promoter <- paste(indir, "/leaf_promoter_open_regions_3000-3000.txt", sep = "")
root_open_promoter <- paste(indir, "/root_promoter_open_regions_3000-3000.txt", sep = "")
all_tss <- paste(indir, "/aligned.peaks.annotated.capped.filtered", sep = "")

col_names = c("tss_id", "oc_id", "chr", "start", "end", "rel_start", "rel_end")
leaf_oc <- read.table(leaf_open_promoter, col.names = col_names)
root_oc <- read.table(root_open_promoter, col.names = col_names)

left <- -3000
right <- 3000

plot_oc(root_oc, leaf_oc, left, right, "ROOT and LEAF OC comparison")

#### Extract the OC coverage for leaf/root specific genes (diff expressed promoters only)
diff <- read.table("~/Downloads/diff_exp_results.txt", header = T)
diff_root <- diff[diff$b > 4 & diff$qval < 0.05,]
diff_leaf <- diff[diff$b < -4 & diff$qval < 0.05,]

library(tidyr)
root_oc <- extract(root_oc, tss_id, c("transcript_id", "chr", "pstart", "strand"), regex = "(.+?)_(Chr\\d+)_(\\d+)_(.)", remove = F)
leaf_oc <- extract(leaf_oc, tss_id, c("transcript_id", "chr", "pstart", "strand"), regex = "(.+?)_(Chr\\d+)_(\\d+)_(.)", remove = F)

root_oc_root <- root_oc[root_oc$transcript_id %in% diff_root$Accession,]
leaf_oc_root <- leaf_oc[leaf_oc$transcript_id %in% diff_root$Accession,]

leaf_oc_leaf <- leaf_oc[leaf_oc$transcript_id %in% diff_leaf$Accession,]
root_oc_leaf <- root_oc[root_oc$transcript_id %in% diff_leaf$Accession,]


plot_oc(root_oc_root, leaf_oc_leaf, left, right, "root OC in root-specific promoters vs. leaf OC in leaf-specific promoters")
plot_oc(root_oc_leaf, leaf_oc_leaf, left, right, "root OC in leaf-specific promoters vs. leaf OC in leaf-specific promoters")
plot_oc(root_oc_root, leaf_oc_root, left, right, "root OC in root-specific promoters vs. leaf OC in root-specific promoters")
plot_oc(root_oc_leaf, leaf_oc_root, left, right, "root OC in leaf-specific promoters vs. leaf OC in root-specific promoters")

plot_oc_summary(root_oc_root, leaf_oc_root, left, right, "root OC in root-specific promoters vs. leaf OC in root-specific promoters")
plot_oc_summary(root_oc_leaf, leaf_oc_leaf, left, right, "root OC in leaf-specific promoters vs. leaf OC in leaf-specific promoters")


#### Additional tests
sets<-function(start, end, group, overlap=length(unique(group))) {
  dd<-rbind(data.frame(pos=start, event=1), data.frame(pos=end, event=-1))
  dd<-aggregate(event~pos, dd, sum)
  dd<-dd[order(dd$pos),]
  dd$open <- cumsum(dd$event)
  r<-rle(dd$open>=overlap)
  ex<-cumsum(r$lengths-1 + rep(1, length(r$lengths))) 
  sx<-ex-r$lengths+1
  cbind(dd$pos[sx[r$values]],dd$pos[ex[r$values]+1])
  
} 

with(root_oc_in_region, sets(rel_start, rel_end, tss_id, 1))

aggregate(root_oc_in_region, by = "rel_start", FUN = nrow)
