library(ggplot2)
indir = "."
leaf_open_promoter <- paste(indir, "/oc/leaf_promoter_open_regions_3000-3000.txt", sep = "")
root_open_promoter <- paste(indir, "/oc/root_promoter_open_regions_3000-3000.txt", sep = "")
col_names = c("tss_id", "oc_id", "chr", "start", "end", "rel_start", "rel_end")
leaf_oc <- read.table(leaf_open_promoter, col.names = col_names)
root_oc <- read.table(root_open_promoter, col.names = col_names)

all_tss <- paste(indir, "/aligned.peaks.annotated.capped.filtered", sep = "")
de_trx <- read.table(paste(indir, "ath_root_leaf_rsem_deseq_diff_expr_results_filtered.txt", sep = ""), header = T)
up_down_trx <- de_trx[abs(de_trx$foldChange) > 3 & (de_trx$mean_leaf_norm > 300 | de_trx$mean_root_norm > 300),]
#######
df_tss <- read.delim(all_tss, sep = ",", header = T)
head(df_tss)

df_tss$tss_id <- paste(df_tss$Chromosome, df_tss$TranscriptID, df_tss$Strand, df_tss$ModeLocation, df_tss$tissue, sep = "_")
aggregate(df_tss$tss_id, list(df_tss$tissue), length)
tss_per_trx <- aggregate(df_tss$tss_id, list(df_tss$tissue, df_tss$TranscriptID), length)
head(tss_per_trx)
a = as.data.frame(table(tss_per_trx[tss_per_trx$Group.1 == "root",]$x))
a$tissue <- "root"
b = as.data.frame(table(tss_per_trx[tss_per_trx$Group.1 == "leaf",]$x))
b$tissue = "shoot"
counts = rbind(a,b)
colnames(counts) <- c("num_transcripts", "tss_count", "tissue")
ggplot(counts, aes(num_transcripts, tss_count, col = tissue, shape = tissue)) + geom_point() + theme_bw() + 
  ggtitle("Number of mapped TSS peaks  and transcript counts")

mapped_loc <- aggregate(df_tss$tss_id, list(df_tss$tissue, df_tss$TranscriptLocation), length)
colnames(mapped_loc) <- c("Tissue", "Location", "No_TSSs")
ggplot(mapped_loc, aes(x = Location, y = No_TSSs, fill = Tissue)) +
  geom_col() + ylab("Number of TSS peaks") + xlab("Transcript Location") +
  theme_bw() + ggtitle("Promoter location for mapped TSS peaks relative to Transcripts")

length(unique(df_tss$TranscriptID))

aggregate(df_tss$TranscriptID, list(df_tss$tissue, df_tss$TranscriptID), length)

root_tss <- df_tss[df_tss$tissue == "root",]
leaf_tss <- df_tss[df_tss$tissue == "leaf",]

length(unique(root_tss$TranscriptID))
length(unique(leaf_tss$TranscriptID))

length(unique(root_tss$GeneName))
length(unique(leaf_tss$GeneName))
length(unique(df_tss$GeneName))



both <- merge(root_tss, leaf_tss, by = "TranscriptID")
both$left <- ifelse(both$Start.x > both$Start.y, both$Start.x, both$Start.y) 
both$right <- ifelse(both$End.x < both$End.y, both$End.x, both$End.y)
both$overlap <- ifelse(both$Start.x > both$End.y,0, ifelse(both$End.x < both$Start.y, 0, abs(both$left - both$right)))
both$mod_dist <- abs(both$ModeLocation.x - both$ModeLocation.y)
both$tss_dist_groups <- cut(both$mod_dist, breaks = c(seq(0,4, by=5), seq(5,100, by=96), seq(101,max(both$mod_dist), by=(max(both$mod_dist)-101))), 
                            include.lowest = T)
ggplot(both, aes(factor(tss_dist_groups))) + geom_bar(stat = "count", fill="steelblue") + xlab("TSS mode distance (base pair)") + 
  geom_text(aes(label=count), vjust=1.6, color="white", size=3.5) + theme_minimal()

ggplot(both, aes(factor(tss_dist_groups), fill = tss_dist_groups)) + geom_bar(stat = "count", col = "black") + xlab("TSS mode distance (base pair)") + 
  theme_minimal() + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#A69F00")) + geom_text(stat='count', aes(label=..count..), vjust=-1) +
  scale_x_discrete(labels = c("[0,5]" = "0-5","(5,101]"  ="5-100", "(101,3.42e+03]" = ">100"))

de_tss <- both[both$TranscriptID %in% up_down_trx$Accession,]
ggplot(de_tss, aes(factor(tss_dist_groups), fill = tss_dist_groups)) + geom_bar(stat = "count", col = "black") + xlab("TSS mode distance (base pair)") + 
  theme_minimal() + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#A69F00")) + geom_text(stat='count', aes(label=..count..), vjust=-1) +
  scale_x_discrete(labels = c("[0,5]" = "0-5","(5,101]"  ="5-100", "(101,3.42e+03]" = ">100"))
ggplot(de_tss, aes(mod_dist)) + geom_histogram(binwidth = 2) + xlab("TSS mode distance (base pair)") + 
  theme_minimal() 


both$tss <- "tss"
ggplot(both, aes(tss, overlap)) + geom_violin() + ylab( "Size of overlap between TSS peaks (basepair)") +
  xlab("") + ggtitle("Size of overlap between root and shoot TSS peaks in each promoter")
mean(both$overlap)
median(both$overlap)
hist(both$overlap)
both$range <- cut(both$overlap, c(seq(-1,50, by = 10)))
table(both$range)

#------------ OC stats on not-DE transcripts not having diff TSS ---------
nde_tss <- de_tss[de_tss$mod_dist < 30,]
length(unique(nde_tss$tss_id.x))

nde_tss$tss_id.x <- paste(nde_tss$TranscriptID, nde_tss$Chromosome.x, nde_tss$ModeLocation.x, nde_tss$Strand.x, sep = "_")
nde_tss$tss_id.y <- paste(nde_tss$TranscriptID, nde_tss$Chromosome.y, nde_tss$ModeLocation.y, nde_tss$Strand.y, sep = "_")
nde_tss <- nde_tss[,c("tss_id.x", "tss_id.y", "mod_dist")]
colnames(nde_tss) <- c("tss_id.root", "tss_id.leaf", "mod_dist")

orig_nde_tss <- nde_tss

nde_tss <- merge(nde_tss, root_oc[,c("tss_id", "start", "end")], by.x = "tss_id.root", by.y = "tss_id")
length(unique(nde_tss$tss_id.root))

nde_tss <- merge(nde_tss, leaf_oc[,c("tss_id", "start", "end")], by.x = "tss_id.leaf", by.y = "tss_id", suffixes = c(".root", ".leaf"))
length(unique(nde_tss$tss_id.root))

not_open_in_root <- orig_nde_tss[!orig_nde_tss$tss_id.root %in% root_oc$tss_id & orig_nde_tss$tss_id.leaf %in% leaf_oc$tss_id,]
not_open_in_leaf <- orig_nde_tss[orig_nde_tss$tss_id.root %in% root_oc$tss_id & !orig_nde_tss$tss_id.leaf %in% leaf_oc$tss_id,]
nrow(not_open_in_root)
nrow(not_open_in_root)

nrow(orig_nde_tss[!orig_nde_tss$tss_id.leaf %in% nde_tss$tss_id.leaf,])

nde_tss$left <- ifelse(nde_tss$start.root > nde_tss$start.leaf, nde_tss$start.root, nde_tss$start.leaf) 
nde_tss$right <- ifelse(nde_tss$end.root < nde_tss$end.leaf, nde_tss$end.root, nde_tss$end.leaf)
nde_tss$l <- ifelse(nde_tss$start.root > nde_tss$start.leaf, nde_tss$start.leaf, nde_tss$start.root) 
nde_tss$r <- ifelse(nde_tss$end.root < nde_tss$end.leaf, nde_tss$end.leaf, nde_tss$end.root)

nde_tss$overlap <- ifelse(nde_tss$start.root > nde_tss$end.leaf,0, ifelse(nde_tss$end.root < nde_tss$start.leaf, 0, 
                                                                          abs(nde_tss$start.root - nde_tss$end.leaf)))
nde_tss$percent_overlap <- nde_tss$overlap / (nde_tss$r - nde_tss$l)
all <- aggregate(nde_tss$percent_overlap, list(nde_tss$tss_id.root, nde_tss$tss_id.leaf, nde_tss$start.root, nde_tss$end.root), max)
head(all)

hist(all$x, xlab = "% OC overlap between root and shoot promoters", main = "Histogram of open regions comparision between DE transcripts \n having differenet TSS locations in root and shoot")

all_avg <- aggregate(all$x, list(all$Group.1, all$Group.2), mean)
hist(all_avg$x, xlab = "% OC overlap between root and shoot promoters", main = "Histogram of open regions comparision between DE transcripts \n having similar TSS locations in root and shoot")


a = all[all$x < 0.1,]
nrow(all[all$x < 0.1,])
length(unique(a$Group.2))

length(unique(nde_tss$TranscriptID))
length(unique(de_tss$TranscriptID))
nde_tss <- nde_tss[, c("tss_id.leaf", "tss_id.root", "mod_dist", "overlap")]

a = de_tss[1,]
a$tss_id.x
a$tss_id.y

#######


all_peaks <- read.delim(all_tss, sep = ",", header = T)
all_peaks$pstart <- all_peaks$Start - 3000
all_peaks$tss_name <- paste(all_peaks$TranscriptID, all_peaks$Chromosome, all_peaks$ModeLocation, all_peaks$Strand, sep = "_")

col_names = c("tss_id", "oc_id", "chr", "start", "end", "rel_start", "rel_end")
leaf_oc <- read.table(leaf_open_promoter, col.names = col_names)
root_oc <- read.table(root_open_promoter, col.names = col_names)

closed_root <- all_peaks[!all_peaks$tss_name %in% root_oc$tss_id,]
nrow(closed_root)

closed_shoot <- all_peaks[!all_peaks$tss_name %in% leaf_oc$tss_id,]
nrow(closed_shoot)

closed_both <- all_peaks[!all_peaks$tss_name %in% leaf_oc$tss_id & !all_peaks$tss_name %in% root_oc$tss_id,]
nrow(closed_both)

m = merge(all_peaks, root_oc, by.x = "tss_name", by.y = "tss_id")
close_root <- all_peaks[all_peaks$tss_name %in% root_oc$tss_id,]
length(unique(root_oc$tss_id))
length(unique(leaf_oc$tss_id))

root_oc <- leaf_oc

#all
left = -1000
right = 1000
root_oc$allL <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$allR <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)


# -100 to +100
left = -200
right = 200

root_oc$tssL <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$tssR <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)

# 100 to 500
left = 200
right = 600
root_oc$down500L <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$down500R <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)
  
# 1000 to 3000
left = 600
right = 1000
root_oc$down1000L <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$down1000R <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)

left = 1000
right = 3000
root_oc$down3000L <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$down3000R <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)

# -500 to -1000
left = -200
right = -600
root_oc$up500L <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$up500R <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)
  
# 1000 to 3000
left = -600
right = -1000
root_oc$up1000L <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$up1000R <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)

left = -1000
right = -3000
root_oc$up3000L <- ifelse(root_oc$rel_start < left, left, root_oc$rel_start)
root_oc$up3000R <- ifelse(root_oc$rel_end > right, right, root_oc$rel_end)

###########
df <- data.frame(location = c(), bp=c())
#all
left = -1000
right = 1000
root_oc$all <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$allL - root_oc$allR) ))
r = root_oc[root_oc$all > 0, "all"]
d = data.frame(location = "All(-1000_+1000)", bp = r)
df <- rbind(df, d)

left = -1000
right = -3000
root_oc$up3000 <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$up3000L - root_oc$up3000R) ))
r = root_oc[root_oc$up3000 > 0, "up3000"]
d = data.frame(location = "-1000_-3000", bp = r)
#df <- rbind(df, d)

left = -600
right = -1000
root_oc$up1000 <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$up1000L - root_oc$up1000R) ))
r = root_oc[root_oc$up1000 > 0, "up1000"]
d = data.frame(location = "-600_-1000", bp = r)
df <- rbind(df, d)

left = -200
right = -600
root_oc$up500 <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$up500L - root_oc$up500R) ))
r = root_oc[root_oc$up500 > 0, "up500"]
d = data.frame(location = "-200_-600", bp = r)
df <- rbind(df, d)

left = -200
right = 200
root_oc$TSS <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$tssR - root_oc$tssL) ))
r = root_oc[root_oc$TSS > 0, "TSS"]
d = data.frame(location = "TSS", bp = r)
df <- rbind(df, d)

left = 200
right = 600
root_oc$down500 <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$down500L - root_oc$down500R) ))
r = root_oc[root_oc$down500 > 0, "down500"]
d = data.frame(location = "200-600", bp = r)
df <- rbind(df, d)

left = 600
right = 1000
root_oc$down1000 <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$down1000L - root_oc$down1000R) ))
r = root_oc[root_oc$down1000 > 0, "down1000"]
d = data.frame(location = "600-1000", bp = r)
df <- rbind(df, d)

left = 1000
right = 3000
root_oc$down3000 <- ifelse(root_oc$rel_end < left, 0, ifelse(root_oc$rel_start > right, 0, abs(root_oc$down3000L - root_oc$down3000R) ))
r = root_oc[root_oc$down3000 > 0, "down3000"]
d = data.frame(location = "1000-3000", bp = r)
#df <- rbind(df, d)

ggplot(df, aes(x = location, y = bp, fill = location)) +
  geom_violin() + ylab("Size of Open region (basepair)") + xlab("Relative Location to TSS") +
  theme_bw() + ggtitle("Coverage Regions of Open Chromatin (Shoot)")

table(df$location)

v = root_oc[root_oc$down3000 >0,]
head(df)

aggregate(df$bp, list(df$location), mean)

