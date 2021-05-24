indir = "tss_peak_calls_root/"
indir = "eval_depth/tss_peak_calls_leaf/"
itr = seq(1, 5)
sample_no = seq(10, 90, by = 10)


read_peaks <- function(percent_read, itr) {
  infile = paste(indir, "/sample.N", itr, ".", percent_read, ".anno.csv", sep = "")
  print(infile)
  df <- read.csv(infile)
  print(head(df))
  df <- df[,c("GeneName", "TranscriptLocation")]
  df <- df[df$TranscriptLocation %in% c("tss", "5'utr","<250", "3'utr", "coding", "intergenic"),]
  df$percent_reads <- percent_read
  df$iteration = itr
  return(df)
}

all <- do.call(rbind, lapply(sample_no, read_peaks, 2 ))
m = aggregate(all$GeneName, list(all$TranscriptLocation, all$percent_reads), length)
head(m)
colnames(m) <- c("Peak_Location", "percent_sampled_reads", "count")
ggplot(m, aes(percent_sampled_reads, count, fill = Peak_Location)) + geom_col() +scale_fill_brewer(palette="Dark2") + theme_bw() + xlab("% Sampled TSS reads") + ggtitle("nanoCAGEXL Root") +
  theme(axis.text=element_text(size=12),  axis.title=element_text(size=12,face="bold"), legend.text =element_text(size=14), title = element_text(size=14,face="bold"))
ggplot(m, aes(percent_sampled_reads, count, fill = Peak_Location)) + geom_col() +scale_fill_brewer(palette="Dark2") + theme_bw() + xlab("% Sampled TSS reads") + ggtitle("nanoCAGEXL Shoot")+
  theme(axis.text=element_text(size=12),  axis.title=element_text(size=12,face="bold"), legend.text =element_text(size=14), title = element_text(size=14,face="bold"))
