#!/usr/bin/Rscript

library(tidyr)

args = commandArgs(trailingOnly = T)
in_coefs_file = args[1]
outdir = args[2]
from = as.numeric(args[3])
to = as.numeric(args[4])
add_by = as.numeric(args[5])

coef_table <- read.table(in_coefs_file, col.names=c("feature", "coef"))

topx_counts <- seq(from ,to, by = add_by)

coef_table  <- extract(coef_table, feature, into = c("pwm", "strand", "win"), regex = "(.+?)_(FWD|REV)_(\\d+)", remove = F)
coef_table <- as.data.frame(unique(coef_table[, "pwm"]))
for (i in topx_counts) {
  outfile = paste(outdir, "/top", i, "_pwms.txt", sep = "")
  write.table(coef_table[1:i,], file = outfile, col.names = F, row.names = F, quote = F)
}
