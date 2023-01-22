library(openxlsx)
library(data.table)
library(dplyr)
library(circlize)
library(ggplot2)
library(ggfortify)
library(singscore)
library(ComplexHeatmap)
library(factoextra)
library(GenVisR)
library(foreach)
library(data.table)
library(ggpubr)
setwd("~/Box Sync/UCSF/Manuscript/Submissions/GitHub/")

### source functions
source("./functions.R")

load(file="./data/WCDT.5groups.210.RData")

### RNA-seq data
TPM.WCDT <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes.txt",sep="\t",header = T,data.table = F)
TPM.PROMOTE <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt",sep="\t",header = T,data.table = F)

RNA.TPM.df <-  right_join(TPM.WCDT,TPM.PROMOTE)
RNA.TPM.df <- RNA.TPM.df[!duplicated(RNA.TPM.df$gene_name),]
rownames(RNA.TPM.df) <- RNA.TPM.df$gene_name;RNA.TPM.df<- RNA.TPM.df[,-c(1,which(colnames(RNA.TPM.df)=='gene_name'))]

######
# KLF family expression in mCRPC

KLF.compare.df <- as.data.frame(t(log2(RNA.TPM.df+1)[rownames(RNA.TPM.df)%in%c("KLF3","KLF1","KLF14","KLF4","KLF6","KLF9","KLF5"),]))
KLF.compare.df$samples <- rownames(KLF.compare.df)
KLF.compare.df <- left_join(Nelson.samples.PCA[,-3],KLF.compare.df)
KLF.compare.df$Subtypes <- factor(KLF.compare.df$Subtypes,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE-','AR+/NE+'))
my_comparisons <- list(c("AR-/NE-", "AR-/NE+"),c("AR-/NE-", "ARL/NE-"),c("AR-/NE-", "AR+/NE-"),c("AR-/NE-", "AR+/NE+"))

KLF.df <- melt(KLF.compare.df)
KLF.df$variable <- factor(KLF.df$variable,levels=c('KLF1','KLF3','KLF4','KLF5','KLF6','KLF9','KLF13','KLF14'))

SFigure9_boxplot <- ggboxplot(KLF.df,x = 'Subtypes',y='value',fill='Subtypes',ylab = 'log2 (TPM+1)',xlab = 'Subtypes',
                              bxp.errorbar = T,bxp.errorbar.width = 0.05,add = 'jitter',add.params = list(alpha=0.7,size=0.6))+
    facet_wrap(~variable,nrow = 2)+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+
    stat_compare_means(size=4,label.x=3,fontface='bold',label.y=8)+
    stat_compare_means(size=6,tip.length = 0.015,label.y=8,label.y.npc = 0.5,
                       step.increase = 0.2,comparisons = my_comparisons,hide.ns = F,label = 'p.signif',p.adjust.methods='fdr',method='t.test')+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),
          axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5,size=15),
          axis.text.y = element_text(size=15),
          legend.position = 'right',strip.text = element_text(face='bold',size=12.5),
          legend.title = element_text(size=15,face = 'bold'),
          legend.key.size = unit(10,units = 'mm'),legend.text = element_text(size=13),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())+coord_cartesian(ylim=c(0,14))

