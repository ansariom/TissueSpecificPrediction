diffs_colnames <- c("gene_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs",
                    "tech_var", "sigma_sq", "smooth_sigma_sq", "final_sigma_sq",
                    "chr", "loc", "strand", "offset?")

infile = ""
load(infile)

row.names(all_features_diffs_wide) <- all_features_diffs_wide$tss_name
diff_cols <- all_features_diffs_wide[, colnames(all_features_diffs_wide) %in% diffs_colnames]

a <- all_features_diffs_wide[, !colnames(all_features_diffs_wide) %in% diffs_colnames]
a <- gather(a, feature, value, -tss_name)
a <- extract(a, feature, c("pwm", "strand", "win", "type"), regex = "(.+)_(FWD|REV)_(.)(\\D*)", remove=F)

oc <- a[grep("_OC_", a$type),]
tfbs <- a[!a$feature %in% oc$feature,]
m <- merge(tfbs, oc, by = c("tss_name", "pwm", "strand", "win"))
m$value <- m$value.x * m$value.y
m <- m[,c("tss_name", "feature.y", "value")]
colnames(m) <- c("tss_name", "feature", "value")

diff_cols$tss_name <- rownames(diff_cols)
k <- merge(p, diff_cols, by = "tss_name")