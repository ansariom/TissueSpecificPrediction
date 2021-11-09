#!/usr/bin/Rscript

library(tidyr)
huges_file <- "ath_Hughes_list.txt"
trans_file <- "transfac_genes_mats_multi_rows.txt"
tair10_families = "gene_families_sep_29_09_update.txt"
#diff_file <- "ath_root_leaf_rsem_deseq_diff_expr_results_filtered.txt"
#chosen_pwm_file <- "chosen_PWMs_V3.txt"

# fc = 3
fc = 3
coef_file <- "roe_only_coef_table.txt"

#outfile = paste("model_weight_coefficient_table_with_gene_expr_fc_", fc, ".csv", sep = "")
outfile = paste("coef_table_with_TF_faimily.txt", sep = "")

huges_genes <- read.delim(huges_file, header = T, sep = "\t" )[,c("Motif_ID", "TF_Name", "Family_Name", "TF_Status")]
transcfac_genes <- read.delim(trans_file, sep = "\t", col.names = c("pwm_id", "gene_id"))
tair_families = read.delim(tair10_families, sep = "\t", header = T)[,c(1,6)]
tair_families$gene_id <- toupper(tair_families$Genomic_Locus_Tag)
tair_families$Genomic_Locus_Tag <- NULL

diff <- read.table(diff_file, header = T)
coef_table <- read.table(coef_file, col.names = c("feature", "coef"))
coef_table <- extract(coef_table, feature, c("pwm", "strand", "window"), regex = "(.+?)_(FWD|REV)_(\\d+)", remove = F)
coef_table <- extract(coef_table, pwm, c("mat_id", "suffix"), regex = "(M\\d+)_(.*)", remove = F)
coef_table$mat_id <- ifelse(coef_table$suffix == "1.02", paste(coef_table$mat_id, "_", "1.02", sep = ""), coef_table$mat_id)

a = as.data.frame(unique(coef_table[,c("mat_id")]))
colnames(a) <- "mat_id"
a1 = merge(a, huges_genes, by.x = "mat_id", by.y="Motif_ID")
a2 = merge(a, transcfac_genes, by.x = "mat_id", by.y="pwm_id")
a2 = merge(a2, tair_families, by = "gene_id")
a2$f = gsub("(.+) Transcription Factor\\D*", "\\1", a2$Gene_Family)
a2$f = gsub("basic \\D+\\((\\D+)\\)", "\\1", a2$f)

out = merge(a, unique(a1[,c(1,3)], by = "mat_id"))
a2 = unique(a2[,c(2,4)])
colnames(a2) <- c("mat_id", "Family_Name")
out2 = merge(a, a2, by = "mat_id"
out = rbind(out, out2)
final_out = merge(coef_table, out, by = "mat_id")

final_out = final_out[, c("feature", "coef", "Family_Name")]
final_out = final_out[order(-abs(final_out$coef)),]
write.table(final_out, outfile, quote = F, sep = "\t", row.names = F)

### if you want to add diff information as well!
mat_table <- rbind(transcfac_genes, huges_genes)
diff <- na.omit(diff)
diff <- extract(diff, Accession, c("gene_id"), regex= "(AT\\dG\\d+).\\d", remove = F)

mat_diff <- merge(mat_table, diff, by = "gene_id")
scols <- c("pval", "qval", "b", "mean_leaf_norm", "mean_root_norm")
avg_expr <- aggregate(mat_diff[, colnames(mat_diff) %in% scols], list(gene_id = mat_diff$gene_id), mean)
max_expr <- aggregate(mat_diff[, colnames(mat_diff) %in% scols], list(gene_id = mat_diff$gene_id), max)
min_expr <- aggregate(mat_diff[, colnames(mat_diff) %in% scols], list(gene_id = mat_diff$gene_id), min)

avg_max_exp  <- merge(avg_expr, max_expr, by = "gene_id", suffixes=c("_avg", "_max"))
avg_min_max_exp <- merge(avg_max_exp, min_expr, by = "gene_id")
colnames(avg_min_max_exp)[which(colnames(avg_min_max_exp) == "b")] <- "b_min"
colnames(avg_min_max_exp)[which(colnames(avg_min_max_exp) == "pval")] <- "pval_min"
colnames(avg_min_max_exp)[which(colnames(avg_min_max_exp) == "qval")] <- "qval_min"
colnames(avg_min_max_exp)[which(colnames(avg_min_max_exp) == "mean_leaf_norm")] <- "mean_leaf_norm_min"
colnames(avg_min_max_exp)[which(colnames(avg_min_max_exp) == "mean_root_norm")] <- "mean_root_norm_min"

sel_cols <- c("feature", "coef", "gene_id","b_avg", "b_min", "b_max" , "pval_avg", "pval_min" ,"pval_max", "qval_avg", "qval_min", "qval_max", "mean_leaf_norm_avg", "mean_leaf_norm_min", "mean_leaf_norm_max", "mean_root_norm_avg", "mean_root_norm_min", "mean_root_norm_max")
final_out <- final_out[, sel_cols]
final_out <- final_out[order(-abs(final_out$coef)),]
final_out$tissue <- ifelse(final_out$b_avg > 2, "root", ifelse(final_out$b_avg < -2, "leaf", "None"))
final_out$coef_indication_tissue <- ifelse(final_out$coef < 0, "root", "leaf")
