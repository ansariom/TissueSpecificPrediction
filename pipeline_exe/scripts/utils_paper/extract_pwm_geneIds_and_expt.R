#!/usr/bin/Rscript

library(tidyr)
huges_file <- "huges_genes_mats_multi_rows.txt"
trans_file <- "transfac_genes_mats_multi_rows.txt"
diff_file <- "../db/ath_root_leaf_rsem_deseq_diff_expr_results_filtered.txt"
chosen_pwm_file <- "chosen_PWMs_V3.txt"

# fc = 3
fc = 3
coef_file <- "../ibdc_roe-only/med_high/3/785440/roe_only_coef_table.txt"

# fc 2.5
# fc = 2.5
#coef_file <- "../ibdc_roe-only/med_high/2.5/956178/roe_only_coef_table.txt"
outfile = paste("model_weight_coefficient_table_with_gene_expr_fc_", fc, ".csv", sep = "")

huges_genes <- read.table(huges_file, header = F, col.names=c("mat_id", "gene_id"))
transcfac_genes <- read.delim(trans_file, sep = "\t", col.names = c("mat_id", "gene_id"))
diff <- read.table(diff_file, header = T)
coef_table <- read.table(coef_file, col.names = c("feature", "coef"))

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

coef_table <- extract(coef_table, feature, c("pwm", "strand", "window"), regex = "(.+?)_(FWD|REV)_(\\d+)", remove = F)
coef_table <- extract(coef_table, pwm, c("mat_id", "suffix"), regex = "(M\\d+)_(.*)", remove = F)
coef_table$mat_id <- ifelse(coef_table$suffix == "1.02", paste(coef_table$mat_id, "_", "1.02", sep = ""), coef_table$mat_id)


coef_genes <- merge(coef_table, mat_table, by = "mat_id", all.x = T)
final_out <- merge(coef_genes, avg_min_max_exp, by = "gene_id", all.x = T)

sel_cols <- c("feature", "coef", "gene_id","b_avg", "b_min", "b_max" , "pval_avg", "pval_min" ,"pval_max", "qval_avg", "qval_min", "qval_max", "mean_leaf_norm_avg", "mean_leaf_norm_min", "mean_leaf_norm_max", "mean_root_norm_avg", "mean_root_norm_min", "mean_root_norm_max")
final_out <- final_out[, sel_cols]
final_out <- final_out[order(-abs(final_out$coef)),]
final_out$tissue <- ifelse(final_out$b_avg > 2, "root", ifelse(final_out$b_avg < -2, "leaf", "None"))
final_out$coef_indication_tissue <- ifelse(final_out$coef < 0, "root", "leaf")
final_out$feature_type <- "TFBS"
final_out$feature_type[grepl("(OC_P_ROOT)|(OC_P_LEAF)", final_out$feature)] <- "OC"


write.csv(final_out, outfile, quote = F, sep = ",", row.names = F)
