library(data.table)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(plyr)
setwd("~/Box Sync/UCSF/Manuscript/Submissions/GitHub/")
source("./functions.R")

load(file="./data/WCDT.5groups.210.RData")
TPM.WCDT <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes.txt",sep="\t",header = T,data.table = F)
TPM.PROMOTE <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt",sep="\t",header = T,data.table = F)

RNA.TPM.df <-  right_join(TPM.WCDT, TPM.PROMOTE)
RNA.TPM.df <- RNA.TPM.df[!duplicated(RNA.TPM.df$gene_name),]
rownames(RNA.TPM.df) <- RNA.TPM.df$gene_name;RNA.TPM.df<- RNA.TPM.df[,-c(1,which(colnames(RNA.TPM.df)=='gene_name'))]

CHD7.compare.df <- as.data.frame(t(log2(RNA.TPM.df+1)[rownames(RNA.TPM.df)%in%toupper(c('CHD7','SOX2')),]))
CHD7.compare.df$samples <- rownames(CHD7.compare.df)
CHD7.compare.df <- left_join(Nelson.samples.PCA[,-3],CHD7.compare.df)
CHD7.compare.df$Subtypes <- factor(CHD7.compare.df$Subtypes,levels=c('AR-/NE-','AR-/NE+','AR+/NE-','ARL/NE-','AR+/NE+'))
my_comparisons <- list(c("AR-/NE+", "AR-/NE-"),c("AR-/NE+", "AR+/NE-"),c("AR-/NE+", "ARL/NE-"),c("AR-/NE+", "AR+/NE+"))


expression.SOX2.bplot <- ggboxplot(CHD7.compare.df,x='Subtypes',y='SOX2',ylab = 'log2 (TPM+1)',xlab='Subtypes',fill='Subtypes',add='jitter',add.params = list(alpha=0.7,color='black',fill='black'),
                                   bxp.errorbar = T,bxp.errorbar.width = 0.05)+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+
    scale_color_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+
    stat_compare_means(size=5,label.x=3.3,fontface='bold',label.y=10)+
    stat_compare_means(size=8,tip.length = 0.015,label.y=10.2,label.y.npc = 0.5,step.increase = 0.1,comparisons = my_comparisons,hide.ns = F,label = 'p.signif',p.adjust.methods='fdr')+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),
          axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5,size=15),
          axis.text.y = element_text(size=15),
          legend.position = 'none',strip.text = element_text(face='bold',size=12.5),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())+coord_cartesian(ylim=c(0,15))

CHD7.SOX2.compare.df.scatterplot <- ggscatter(CHD7.compare.df,x='CHD7',y='SOX2',xlab='CHD7 (log2 TPM+1)',ylab = 'SOX2 (log2 TPM+1)',color='black',fill='Subtypes',shape = 21,size=5,conf.int = T,
                                              add = "reg.line",
                                              add.params = list(fill = "lightgray"))+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+ stat_cor(method = "pearson", label.x = 2, label.y = 9,p.accuracy = 0.001, r.accuracy = 0.01,cex=8)+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          legend.position = 'none',strip.text = element_text(face='bold',size=12.5),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())


SFigure6.bp <- plot_grid(expression.SOX2.bplot,NULL,CHD7.SOX2.compare.df.scatterplot,align = 'hv',ncol = 3,rel_widths = c(1,0.2,1))


