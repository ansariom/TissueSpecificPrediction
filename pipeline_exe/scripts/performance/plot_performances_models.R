#!/usr/bin/Rscript

infile1 = "~/Downloads/ibdc/aug2020/main.performance_roe.txt"
infile2 = "~/Downloads/ibdc/aug2020/main.performance_tile.txt"
name1 = "TEP-ROE"
name2 = "TEP-Tiled"
model_type = "TEP"
outdir = "~/Downloads/ibdc/aug2020/"

args = commandArgs(trailingOnly = T)

infile1 = args[1]
infile2 = args[2]
name1 = args[3]
name2 = args[4]
model_type= args[5]
outdir = args[6]

outfile_auc = paste(outdir, model_type, "_performance.plots.ROE-Tile.png", sep = "")
#outfile_prc = paste(infile, ".plots.auprc.png", sep = "")

df1 <- read.table(infile1, header = F, col.names = c("auROC", "auPRC", "Expr_Level", "fold_change"))
df1$model_name <- "TEP-ROE"
df2 <- read.table(infile2, header = F, col.names = c("auROC", "auPRC", "Expr_Level", "fold_change"))
df2$model_name <- "TEP-Tiled"

both = rbind(df1,df2)
both$fold_change = NULL
both$Expr_Level = NULL

library(ggplot2)
library(tidyr)

ggplot(both, aes(factor(model_name), auROC,  col = model_name)) + geom_boxplot() + xlab("Model Name") + 
  geom_jitter(position=position_jitter(0.1))+ scale_color_brewer(palette="Dark2") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold") ) + 
  ggtitle("TEP-ROE & TEP_Tiled Model Performance Stability - auROC ") + theme_bw()+
  theme(axis.text=element_text(size=12),  axis.title=element_text(size=12,face="bold"), legend.text =element_text(size=14), title = element_text(size=14,face="bold"))

ggplot(both, aes(factor(model_name), auPRC,  col = model_name)) + geom_boxplot() + xlab("Model Name") + 
  geom_jitter(position=position_jitter(0.1))+ scale_color_brewer(palette="Dark2") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold") ) + 
  geom_jitter(position=position_jitter(0.1))+ scale_color_brewer(palette="Dark2") + 
  ggtitle("TEP-ROE & TEP_Tiled Model Performance Stability - auPRC ") + theme_bw()+
  theme(axis.text=element_text(size=12),  axis.title=element_text(size=12,face="bold"), legend.text =element_text(size=14), title = element_text(size=14,face="bold"))


df1$fold_change <- as.factor(df1$fold_change)
df1 <- gather(df1, perf_metric, value, -fold_change, -Expr_Level)
p1 = ggplot(df1, aes(fold_change, value, color = fold_change)) + geom_boxplot() + facet_grid(perf_metric ~ Expr_Level) + ggtitle(paste(name1, " Performance", sep = "")) + ylim(0.8,0.96)

df2$fold_change <- as.factor(df2$fold_change)
df2 <- gather(df2, perf_metric, value, -fold_change, -Expr_Level)
p2 = ggplot(df2, aes(fold_change, value, color = fold_change)) + geom_boxplot() + facet_grid(perf_metric ~ Expr_Level) + ggtitle(paste(name2, " Performance", sep = "")) + ylim(0.8,0.96)

library(grid)
library(gridExtra)

g = grid.arrange(p1, p2, nrow = 1)

ggsave(g, file = outfile_auc, dpi = 320, width = 12)


