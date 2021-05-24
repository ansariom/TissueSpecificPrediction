#!/usr/bin/Rscript
library(tidyr)
args = commandArgs(trailingOnly = T)

indir = args[1]
outbase = args[2]

fs <- list.files(indir)
all <- c(feature=c(), coef=c(), group = c())
for (f in fs) {
  infile = paste(indir, "/", f, sep = "")
  df <- read.table(infile, header = F, col.names = c("feature", "coef", "group"))
  if (nrow(df) == 0) next
  all <- rbind(df, all)
}
print(head(all))
#all <- extract(all, group, into = "group", regex = "performances_top(\\d+).txt")
#all$group <- as.factor(as.numeric(all$group))
library(ggplot2)

outfile = paste(outbase, ".png", sep = "")
g = ggplot(all, aes(feature, coef, color = group)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) + facet_wrap(~group, ncol = 1)
ggsave(g, file = outfile)


#-------
