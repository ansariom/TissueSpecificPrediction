#!/usr/bin/Rscript

tf_info_file = "/dfs/Megraw_Lab/data/public_db/HughesPWM/Homo_sapiens_2021_05_21_3_48_pm/TF_Information.txt"
df = read.delim(tf_info_file, header = T, sep="\t")

# Take experimental PWMs only -- > 2891 PWMs
df = df[df$MSource_Type %in% c("SELEX", "PBM"),]
df = df[df$TF_Status == "D",]
# take only one sample from each TF
a = do.call(rbind, by(df, df[,c("DBID")], function(x) {return(x[1,])}))

write.table(a$Motif_ID, "TF_ids_filtered_SELEX_PBM_1perTF.txt", sep = "\t", quote = F, col.names = F, row.names = F)
