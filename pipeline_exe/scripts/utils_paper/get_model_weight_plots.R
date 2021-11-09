library(ggplot2)
a = read.delim("~/Downloads/coef_table_with_TF_faimily.txt", header = T, sep = "\t")
head(a)

a$type[grepl("FWD|REV", a$feature)] <- "TFBS"
a$type[grepl("(OC_P_ROOT)|(OC_P_LEAF)", a$feature)] <- "OC"
a$coef_abs <- abs(a$coef)
a$Family_Name[a$Family_Name == "MADS box"] <- "MADS-box" 
a$Family_Name[a$Family_Name == "MADS"] <- "MADS-box"
a$Family_Name[a$Family_Name == "NAC"] <- "NAC/NAM"
a$Family_Name[a$Family_Name == "MYB"] <- "Myb/SANT"
a$Family_Name[a$Family_Name == "MYB-Related"] <- "Myb/SANT"
a$Family_Name[a$Family_Name == "E2F-DP"] <- "E2F"
a$Family_Name[a$Family_Name == "AP2-EREBP"] <- "AP2"
a$Family_Name[a$Family_Name == "BBR/BPC-family of GAGA-motif binding transcription factors"] <- "BBR/BPC"
a = a[!grepl("Superfamily", a$Family_Name),] 
a = a[!a$Family_Name %in% '`',]
a = a[a$coef != 0,]


family_sumWs = aggregate(a$coef_abs, list(a$Family_Name, a$type), sum)
head(family_sumWs)
colnames(family_sumWs) <- c("TF_family", "feature_type", "sum_coefs")
family_sumWs$percent_weight = family_sumWs$sum_coefs / sum(family_sumWs$sum_coefs)

family_sumWs$TF_family <- factor(family_sumWs$TF_family,                               
                  levels = unique(family_sumWs$TF_family[order(family_sumWs$percent_weight, decreasing = TRUE)]))

ggplot(family_sumWs, aes(TF_family, percent_weight, fill = feature_type, group=feature_type)) + geom_bar(position="stack", stat="identity") + coord_flip()

######
## considering the weight signs

head(a)
a$tissue <- ifelse(a$coef < 0, -1, 1)
family_sumWs = aggregate(a$coef_abs, list(a$Family_Name, a$type,a$tissue), sum)
head(family_sumWs)
colnames(family_sumWs) <- c("TF_family", "feature_type", "tissue", "sum_coefs")
family_sumWs$percent_weight = (family_sumWs$sum_coefs * family_sumWs$tissue)/ sum(family_sumWs$sum_coefs)

family_sumWs$TF_family <- factor(family_sumWs$TF_family,                               
                                 levels = unique(family_sumWs$TF_family[order(family_sumWs$percent_weight, decreasing = TRUE)]))

family_sumWs$tissue[family_sumWs$tissue == -1] <- "Root"
family_sumWs$tissue[family_sumWs$tissue == 1] <- "Shoot"

ggplot(family_sumWs, aes(TF_family, abs(percent_weight), fill = feature_type, group=feature_type)) + geom_bar(position="stack", stat="identity") + coord_flip() +
   theme_bw() + theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold")) +
  ylab("% Model Weight") + xlab("TF-Family of Model Features")

family_sumWs$featureType <- paste(family_sumWs$feature_type, family_sumWs$tissue, sep = "_")
ggplot(family_sumWs, aes(TF_family, percent_weight, fill = featureType, group=feature_type)) + geom_bar(position="stack", stat="identity") + coord_flip() +
  scale_fill_manual(values = c("dodgerblue", "orangered1", "sienna4", "springgreen3")) + theme_bw() + theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold")) + 
  ylab("% Model Weight") + xlab("TF-Family of Model Features")


write.table(a, file = "~/Downloads/for_zach_model_weights_with_TF_family.txt", sep = "\t", quote = F, row.names = F)



