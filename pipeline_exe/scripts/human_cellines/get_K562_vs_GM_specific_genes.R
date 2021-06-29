##### Data #####
## https://krishna.gs.washington.edu/content/members/vagar/Xpresso/data/datasets/
################

indir = "."
rnaseq = read.table(paste(indir, "/57epigenomes.RPKM_.pc", sep = ""), header = T)

rpkm = rnaseq
rpkm = log10(rpkm + 0.1)
rpkm$deltaActual = rpkm$K562-rpkm$GM12878

# regenerate the original work results
a= rpkm
smoothScatter(a$GM12878, a$K562, cex.axis=2, cex.lab=2, bty="n", xlab="GM12878 expression level (log10)", ylab="K562 expression level (log10)", xlim=c(-1, 4), ylim=c(-1, 4), las=1, cex=.5)
abline(1,1, col="red")
abline(-1,1, col="red")
text(0, 4, labels = paste("Upregulated in K562:", nrow(a[a$deltaActual > 1,])), offset = 0.5, col="black")
text(3, -1, labels = paste("Upregulated in GM12878:", nrow(a[a$deltaActual < -1,])), offset = 0.5, col="black")

highlight_df <- rpkm[abs(rpkm$fc) >= 10,]

rpkm %>%
  ggplot(aes(x=log10(GM12878),y=log10(K562))) +
  geom_point(alpha = 0.2) + theme_bw() +
  geom_point(data=highlight_df,
             aes(x=log10(GM12878),y=log10(K562)),
             col = "red", size = 3)

# save the gene list for analysis
rpkm$expressed_at <- ifelse(rpkm$deltaActual < -1, "GM12878", ifelse(rpkm$deltaActual > 1, "K562", "-1000"))
table(rpkm$expressed_at)
rpkm = rpkm[rpkm$expressed_at != "-1000",2:5]
rpkm$class = ifelse(rpkm$expressed_at == "K562", 1, 0)
rpkm$gene_id = rownames(rpkm)
write.table(rpkm, file = paste(indir, "/agrawal_human_K562_GM12878_diff_expressed_gennes.txt", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T)
