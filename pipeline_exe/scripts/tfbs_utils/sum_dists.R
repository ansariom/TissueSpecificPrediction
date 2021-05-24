#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

all_sums <- read.table(infile, header = F)

result <- aggregate(all_sums, by = list(all_sums$V1), sum)
colnames(result) <- c("loc", "none", "value")
result$none <- NULL

write.table(result, file = outfile, col.names=F, row.names=F, quote=F, sep = "\t")

