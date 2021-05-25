###################################################################################################

##      data: https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/
##            https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/

####################################################################################################
library(gtools)

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

indir = "~/workspace/TissueSpecificPrediction/publication_data_and_model_distribution/human_cellines_Agrawal_paper/"
# expression data is normalized as rpkm! 
rnaseq = read.table(paste(indir, "/57epigenomes.RPKM_.pc", sep = ""), header = T)

EIds = read.table(paste(indir, "/EIds_to_celline.txt", sep = ""))
rownames(rnaseq) <- rnaseq$gene_id
head(rnaseq)
head(EIds)
cellines = c("Universal_Human_Reference", "K562", "GM12878")
eids_cells = EIds[EIds$V2 %in% cellines,]
rnaseq = rnaseq[,eids_cells$V1]
colnames(rnaseq) <- eids_cells[which(eids_cells$V1 == colnames(rnaseq), arr.ind = T), "V2"]

#asinh.rpkm = as.data.frame(quantile_normalisation(asinh(rnaseq)))
rpkm = rnaseq
rpkm = log10(rpkm + 0.1)
#rpkm$fc = foldchange(rpkm$GM12878 + 1, rpkm$K562 + 1)
#rpkm$log2fc = foldchange2logratio(rpkm$fc, base = 2)
rpkm$deltaActual = rpkm$K562-rpkm$GM12878

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
  
          


                
rpkm = rpkm[rpkm$Universal_Human_Reference > 0,]
rpkm$UK562 = rpkm$K562 / rpkm$Universal_Human_Reference
rpkm$UGM = rpkm$GM12878 / rpkm$Universal_Human_Reference
head(rpkm)
#rpkm = rpkm[rpkm$UGM > 0.5 | rpkm$UK562 > 0.5,]
#rpkm = rpkm[rpkm$GM12878 > 0 | rpkm$K562 > 0,]
rpkm$log10_GM12878 = log10(rpkm$GM12878 + 0.0001)
rpkm$log10_K562 = log10(rpkm$K562 + 0.0001)
rpkm$log_FC = log2((rpkm$GM12878 + 1)/(rpkm$K562 + 1))
rpkm$log_logFC = log2((rpkm$GM12878 + 1)/(rpkm$K562 + 1))
plot(rpkm$log10_GM12878, rpkm$log10_K562, xlab = "GM12878  log10(expression)" , ylab = "K562 log10(expression)")
rpkm$expressed_at <- ifelse(rpkm$log2fc > 2.1, "GM12878", ifelse(rpkm$log2fc < -2.1, "K562", "-1000")) 
table(rpkm$expressed_at)
a = rpkm[(rpkm$log10_GM12878 > 0 & rpkm$log10_K562 < 0) | (rpkm$log10_GM12878 < 0 & rpkm$log10_K562 > 0),]
#plot(asinh.rpkm$GM12878, asinh.rpkm$Universal_Human_Reference)
#plot(asinh.rpkm$K562, asinh.rpkm$Universal_Human_Reference)
#asinh.rpkm$GM12878 = asinh.rpkm$GM12878/asinh.rpkm$Universal_Human_Reference
#asinh.rpkm$K562 = asinh.rpkm$K562/asinh.rpkm$Universal_Human_Reference

#uhrNorm.asinh.rpkm = asinh.rpkm
#plot(uhrNorm.asinh.rpkm$GM12878, uhrNorm.asinh.rpkm$K562)

#norm.data = uhrNorm.asinh.rpkm
#norm.data = na.omit(norm.data)
#norm.data = norm.data[norm.data$Universal_Human_Reference > 0,]
#norm.data = norm.data[norm.data$GM12878 > 15 | norm.data$K562 > 15,]
#norm.data = norm.data + 0.001

norm.data$logFC_KG = log2(norm.data$K562 / norm.data$GM12878)
norm.data$expressed_at <- ifelse(norm.data$logFC_KG > 10, "K562", ifelse(norm.data$logFC_KG < -10, "GM12878", "-1000")) 
table(norm.data$expressed_at)

head(norm.data)
diff_expr = norm.data
rsem_header = c("Accession",	"baseMean",	"mean_leaf_norm",	"mean_root_norm",	"foldChange",	"b",	"pval",	"qval")
diff_expr$Accession = rownames(diff_expr)
diff_expr$baseMean = diff_expr$Universal_Human_Reference
diff_expr$Universal_Human_Reference <- NULL
write.table(diff_expr, file = paste(indir, "/agrawal_human_K562_GM12878_diff_expressed_gennes.txt", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T)

############
## TSS PEAKS

# input list is annotated tpm counts for TSS peaks annoted with nearby gene
# 1st column is peak's genome coordinate, 2nd col is the gene Enterz ID and the 6 following cols are replicate counts for each cell.
# the following command has been used to extract the data from original file:
##
# 1- get the columns that contains Encode counts for non-treated cells as follows:
# >>>>> bash: grep "#" hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt | grep -n '[cell line name]'

# 2- get the related columns along with the coordinate of the peaks and gene annotation from Original file.
# >>>>> bash: grep -v "#" hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt | grep entrezgene | cut -f1,5,177-179,1335-1337 > ~/workspace/IBDC/db/agrawal_db/CAGE_peaks_genes-only_K562_GM12878_tpm_ann.txt
############
cage = read.table(paste(indir, "CAGE_peaks_genes-only_K562_GM12878_tpm_ann.txt", sep = ""), header = T)
head(cage)
colnames(cage) <- c("peak_coord", "gene_id", "GM12878_1", "GM12878_2", "GM12878_3", "K562_1", "K562_2", "K562_3")

# mean of replicates per cell
cage$K562 <- apply(cage[,6:8], 1, mean)
cage$GM12878 <- apply(cage[,3:5], 1, mean)
cage$Chr = unlist(lapply(cage$peak_coord, function(x) return(unlist(strsplit(x, split = ":"))[1])))
cage$gene_id = unlist(lapply(cage$gene_id, function(x) return(unlist(strsplit(x, split = ":"))[2])))
cage$Strand = unlist(lapply(cage$peak_coord, function(x) return(unlist(strsplit(x, split = ","))[2])))
cage$Start = unlist(lapply(cage$peak_coord, function(x) return(unlist(strsplit(unlist(strsplit(x, split = "\\.\\."))[1], split = ":"))[2])))
cage$End = unlist(lapply(cage$peak_coord, function(x) return(unlist(strsplit(unlist(strsplit(x, split = "\\.\\."))[2], split = ","))[1])))

cage = cage[,c("Chr", "Start", "End", "Strand", "gene_id", "GM12878", "K562")]
cage = cage[cage$K562 >0 | cage$GM12878 > 0,]
plot(cage$GM12878, cage$K562)
plot(asinh(cage$GM12878), asinh(cage$K562))
df = as.data.frame(table(cage$Chr))
ggplot(df, aes(df$Var1, Freq)) + geom_col() + ylab("Number of TSS peaks") + xlab("Chromosome")
cage = gather(cage, cell_line, tpm, -Chr, -Start, -End, -Strand, -gene_id)
head(cage)
cage = cage[cage$tpm > 1,]
ggplot(cage, aes(asinh(tpm))) + geom_histogram() + facet_wrap(~cell_line)
head(cage)
cage$tpm_range = cut(cage$tpm, breaks = c(0, 5, 10, 50, 100, 10000))
table(cage$tpm_range)
table(cage$cell_line)


enterz = read.delim(paste(indir, "/ncbi_genes_human.txt", sep = ""), sep = "\t", header = T)[,c(2,4)]
ensml = read.delim(paste(indir, "/ensmbl_gene_list.txt", sep = ""), sep = "\t", col.names = c("ens.id", "ens.name"))
gene_table = merge(enterz, ensml, by.x = "Approved.symbol", by.y = "ens.name")
write.table(gene_table, file = paste(indir, "/enterz_ensembl_gene_conversion_table.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

head(cage)
head(norm.data)


## get data ready for Model

TEP_header = c("Chromosome","Strand","Start","End","ReadCount","ModeLocation","ModeReadCount","Shape","TranscriptLocation","TranscriptID","GeneName","GeneType","X..Capped","tissue")
cage = merge(cage, gene_table, by.x = "gene_id", by.y = "NCBI.gene.ID")
cage = cage[,c("Chr", "Strand", "Start", "End", "tpm", "ens.id", "cell_line")] 
colnames(cage)[1:5] <- TEP_header[1:5]
colnames(cage)[6:7] <- c("GeneName", "tissue")

cage$ModeLocation <- as.numeric(cage$Start) + abs(as.numeric(cage$End) - as.numeric(cage$Start))/2
cage$ModeReadCount = cage$ReadCount/2
cage$Shape = "None"
cage$TranscriptLocation = "tss"
cage$TranscriptID = cage$GeneName
cage$GeneType = "gene"
cage$X..Capped = 100
head(cage)
write.csv(cage[,TEP_header], file = paste(indir, "/FANTOM5_K562_GM12878_CAGE_TSS_peaks_forTEP.csv", sep = ""), quote = F, row.names = F, col.names = T)
