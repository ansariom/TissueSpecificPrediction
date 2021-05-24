###################################################################################################

##      data: https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/
##            https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/

####################################################################################################

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

asinh.rpkm = as.data.frame(quantile_normalisation(asinh(rnaseq)))
plot(asinh.rpkm$GM12878, asinh.rpkm$K562)
plot(asinh.rpkm$GM12878, asinh.rpkm$Universal_Human_Reference)
plot(asinh.rpkm$K562, asinh.rpkm$Universal_Human_Reference)
asinh.rpkm$GM12878 = asinh.rpkm$GM12878/asinh.rpkm$Universal_Human_Reference
asinh.rpkm$K562 = asinh.rpkm$K562/asinh.rpkm$Universal_Human_Reference

uhrNorm.asinh.rpkm = asinh.rpkm
plot(uhrNorm.asinh.rpkm$GM12878, uhrNorm.asinh.rpkm$K562)

norm.data = uhrNorm.asinh.rpkm
norm.data = na.omit(norm.data)
norm.data = norm.data[norm.data$Universal_Human_Reference > 0,]
norm.data = norm.data[norm.data$GM12878 > 15 | norm.data$K562 > 15,]
norm.data = norm.data + 0.001

norm.data$logFC_KG = log2(norm.data$K562 / norm.data$GM12878)
norm.data$expressed_at <- ifelse(norm.data$logFC_KG > 3, "K562", ifelse(norm.data$logFC_KG < -3, "GM12878", "-1000")) 
table(norm.data$expressed_at)

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
cage = cage[,c("Chr", "Strand", "Start", "End", "tpm", "ens.id")] 
colnames(cage)[1:5] <- TEP_header[1:5]
colnames(cage)[6] <- "GeneName"
cage$ModeLocation <- as.numeric(cage$Start) + abs(as.numeric(cage$End) - as.numeric(cage$Start))/2
cage$ModeReadCount = cage$ReadCount/2
head(cage)
cage$Shape = "None"
cage$TranscriptLocation = "tss"
cage$TranscriptID = cage$GeneName
cage$GeneType = "gene"
cage$X..Capped = 100
