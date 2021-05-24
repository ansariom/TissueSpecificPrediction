#!/usr/bin/Rscript
library(tidyr)
args = commandArgs(trailingOnly = T)

indir = args[1]
outbase = args[2]

fs <- list.files(indir)
all <- c(auROC=c(), auPRC=c(), group = c())
for (f in fs) {
  infile = paste(indir, "/", f, sep = "")
  df <- read.table(infile, header = F, col.names = c("auROC", "auPRC"))
  if (nrow(df) == 0) next
  df$group <- basename(infile)
  all <- rbind(df, all)
}
print(head(all))
all <- extract(all, group, into = "group", regex = "performances_top(\\d+).txt")
all$group <- as.factor(as.numeric(all$group))
library(ggplot2)
outfile = paste(outbase, "_performance_plot_auROC.png", sep = "")
all <- gather(all, key = "metric", value = "value", -group)
all$metric = fct_relevel(all$metric, "auROC", "auPRC")
g = ggplot(all, aes(group, value, color = metric)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ xlab("Number of Top Weight Features") + ylab("") +
  theme(axis.text=element_text(size=12),  axis.title=element_text(size=14), legend.text =element_text(size=14), title = element_text(size=14,face="bold"), text = element_text(size=14,face="bold")) +facet_wrap(~metric)
g

ggsave(g, file = outfile)
#-------
