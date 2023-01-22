library(data.table)
library(dplyr)
library(circlize)
library(ggplot2)
library(ComplexHeatmap)
library(rtracklayer)
library(DESeq2)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(cowplot)
library(msigdbr)
library(fgsea)

setwd("~/Box Sync/UCSF/Manuscript/Submissions/Github/")
source("./functions.R")

load(file = "./data/NEPC.gene.list.RData")
load(file="./data/WCDT.5groups.210.RData")
load(file="./data/counts_and_TPMs/Localized.Normal.TPM.RData")

Localized.Normal.TPM <- Localized.Normal.TPM[Localized.Normal.TPM$gene_name %in%NEPC.gene.list$Nelson$GENE_ID,]
rownames(Localized.Normal.TPM) <- Localized.Normal.TPM$gene_name;Localized.Normal.TPM <- Localized.Normal.TPM[,-1]

col_fun = colorRamp2(c(-4,-2.5,0,2.5,4), c("darkblue","royalblue2","white","gold3","gold4"))

TPM.WCDT <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes.txt",sep="\t",header = T,data.table = F)
TPM.PROMOTE <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt",sep="\t",header = T,data.table = F)

RNA.TPM.df <-  right_join(TPM.WCDT, TPM.PROMOTE)
RNA.TPM.df <- RNA.TPM.df[!duplicated(RNA.TPM.df$gene_name),]
RNA.TPM.df <- RNA.TPM.df[RNA.TPM.df$gene_name %in%NEPC.gene.list$Nelson$GENE_ID,]
RNA.TPM.df <- RNA.TPM.df[order(match(RNA.TPM.df$gene_name,NEPC.gene.list$Nelson$GENE_ID)),]
rownames(RNA.TPM.df) <- RNA.TPM.df$gene_name;RNA.TPM.df<- RNA.TPM.df[,-c(1,which(colnames(RNA.TPM.df)=='gene_name'))]
RNA.TPM.df <- RNA.TPM.df[,colnames(RNA.TPM.df)%in%Nelson.samples.PCA$samples]

mix.TPM <- merge(RNA.TPM.df,Localized.Normal.TPM,by='row.names',sort=F)
rownames(mix.TPM) <- mix.TPM$Row.names;mix.TPM <- mix.TPM[,-1]

colData.Nelson = list(GO=c("NE1"='skyblue2',
                           "NE2"='skyblue4',"AR"='forestgreen',"SQUAM"='tomato'))

rowData.Nelson <- data.frame(GO=NEPC.gene.list$Nelson[NEPC.gene.list$Nelson$GENE_ID %in% rownames(mix.TPM),]$class)

rowAnnot.Nelson = rowAnnotation(
    GO=rowData.Nelson$GO,
    col=colData.Nelson,show_annotation_name=F,show_legend=F,annotation_name_gp= gpar(fontsize = 12))
mix.TPM <- log2(mix.TPM+1);mix.TPM <- t(scale(t(mix.TPM),scale = T,center = T));mix.TPM <- as.matrix(mix.TPM)

set.seed(200)
group = kmeans(t(mix.TPM),nstart = 1,centers =5,iter.max = 10000000,algorithm = 'Lloyd')$cluster

my_colour = list('Disease stage' = c('Normal'='tomato','Localized' = 'gray','mCRPC'='black'),
                 'Molecular subtypes'=structure(c(subtype.colors$colors,'white'),
                                                names = c(as.character(subtype.colors$subtypes),"NA")))
type.df <- ifelse(colnames(mix.TPM)%in%colnames(RNA.TPM.df),"mCRPC",
                  ifelse(colnames(mix.TPM)%like%'CPCG','Localized',
                         ifelse(colnames(mix.TPM)%like%'SRR','Normal',NA)))
subtype.df <- ifelse(colnames(mix.TPM)%in%Nelson.samples.PCA[Nelson.samples.PCA$Subtypes=='AR-/NE-',]$samples,"AR-/NE-",
                     ifelse(colnames(mix.TPM)%in%Nelson.samples.PCA[Nelson.samples.PCA$Subtypes=='AR-/NE+',]$samples,"AR-/NE+","NA"))

subtype.df <- ifelse(colnames(mix.TPM)%in%Nelson.samples.PCA[Nelson.samples.PCA$Subtypes=='AR-/NE-',]$samples,"AR-/NE-",
                     ifelse(colnames(mix.TPM)%in%Nelson.samples.PCA[Nelson.samples.PCA$Subtypes=='AR+/NE-',]$samples,"AR+/NE-",
                            ifelse(colnames(mix.TPM)%in%Nelson.samples.PCA[Nelson.samples.PCA$Subtypes=='ARL/NE-',]$samples,"ARL/NE-",
                                   ifelse(colnames(mix.TPM)%in%Nelson.samples.PCA[Nelson.samples.PCA$Subtypes=='AR-/NE+',]$samples,"AR-/NE+",
                                          ifelse(colnames(mix.TPM)%in%Nelson.samples.PCA[Nelson.samples.PCA$Subtypes=='AR+/NE+',]$samples,"AR+/NE+","NA")))))


colAnnot.df = columnAnnotation(
    'Disease stage'=type.df[-40],
    'Molecular subtypes'=subtype.df[-40],simple_anno_size = unit(0.2, "cm"),
    col=my_colour,gp=gpar(col='black',lwd=0.1),show_legend=F,
    annotation_name_gp= gpar(fontsize = 4))


#####
ht_opt$TITLE_PADDING = unit(c(5,5), "points")
set.seed(5)
SFigure2.hm <- draw(Heatmap(matrix = mix.TPM[,-40],name = 'Zscore of log2(TPM)',border = T,cluster_rows = F,heatmap_legend_param = list(direction='vertical'),
                            top_annotation = colAnnot.df,na_col = 'white',
                            column_dend_reorder = T,
                            col = col_fun,
                            rect_gp = gpar(col = "black", lwd = 0.05),show_heatmap_legend = F,
                            show_column_dend = T,
                            show_column_names = F,
                            row_split  = factor(rowData.Nelson$GO,levels=c('NE1','NE2','AR','SQUAM')),row_names_gp = gpar(fontsize=4),cluster_row_slices = T,row_dend_reorder = F,
                            column_names_gp = gpar(fontsize=5),row_title_gp = gpar(fontsize=5,fontface='bold'),clustering_distance_columns = "pearson"
))
