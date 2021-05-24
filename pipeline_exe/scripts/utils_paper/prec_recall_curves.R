#!/usr/bin/Rscript
library(PRROC)

args <- commandArgs(trailingOnly = TRUE)
#featureInfoFile = "~/Downloads/featureInfo_hardCodedSoftCoded.rdat"

featureInfoFile <- args[1]
outdir <- args[2]

load(featureInfoFile)

d <- diffs_classes[diffs_classes$class != -1000,]
fg <- d[d$class == 1, "prob1"]
bg <- d[d$class == 0, "prob1"]
roc<-roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)

outfile = paste(outdir, "/roc_curve.png", sep = "")
png(outfile)
plot(roc)
dev.off()

a <- roc$curve
colnames(a) <- c("fpr", "recall", "threshold")
outfile = paste(outdir, "/roc_curve.table.txt", sep = "")
write.table(a, file = outfile, sep = "\t", col.names = T, row.names = F, quote = F)

prc<-pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
outfile = paste(outdir, "/prc_curve.png", sep = "")
png(outfile)
plot(prc)
dev.off()

a <- prc$curve
colnames(a) <- c("recall", "precision", "threshold")
outfile = paste(outdir, "/prc_curve.table.txt", sep = "")
write.table(a, file = outfile, sep = "\t", col.names = T, row.names = F, quote = F)
