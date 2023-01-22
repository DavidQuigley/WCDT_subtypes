library(data.table)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(plyr)
library(foreach)
library(readxl)
library(genefilter)
library(circlize)
library(universalmotif)
library(motifStack)
library(ggVennDiagram)
library(gtools)
library(GenomicRanges)
library(stringr)
library(purrr)
library(genomation)
library(annotatr)
library(graphics)
library(methylSig)
library(plyranges)
library(ComplexHeatmap)
library(msigdbr)
library(clusterProfiler)
library(ggrepel)
library(ggthemes)

setwd("~/Box Sync/UCSF/Manuscript/Submissions/GitHub/")
source("./functions.R")

load(file="./data/WCDT.5groups.210.RData")
load(file="./data/2022_11_04_WCDT_poppy.Rdata")
load(file = "./data/NEPC.gene.list.RData")

Hallmarks_db <- msigdbr(species = "Homo sapiens",category = c('H'))
WIKI_db <- msigdbr(species = "Homo sapiens",category = c('C2'),subcategory = "CP:WIKIPATHWAYS")
Hallmarks_list = split(x = Hallmarks_db$gene_symbol, f = Hallmarks_db$gs_id)
WIKI_list = split(x = WIKI_db$gene_symbol, f = WIKI_db$gs_id)
genes.list.db <- c(Hallmarks_list)
genes.db <- rbind(Hallmarks_db)
pid2gene <- genes.db %>% dplyr::select(gs_id,human_gene_symbol) #TERM2GENE
pid2name <- genes.db %>% dplyr::select(gs_id,gs_name) #TERM2NAME

chromosomes = c(paste('chr', 1:22, sep=''), 'chrX', 'chrY')
genes_info <- c("hg38_genes_3UTRs","hg38_genes_5UTRs","hg38_genes_intergenic","hg38_genes_introns","hg38_genes_promoters","hg38_genes_exons")

genes.GR <- build_annotations(genome = 'hg38',annotations = genes_info)
genes.GR <- genes.GR[genes.GR@seqnames%in%chromosomes]

TPM.WCDT <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes.txt",sep="\t",header = T,data.table = F)
TPM.PROMOTE <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt",sep="\t",header = T,data.table = F)

RNA.TPM.df <-  right_join(TPM.WCDT, TPM.PROMOTE)
RNA.TPM.df <- RNA.TPM.df[!duplicated(RNA.TPM.df$gene_name),]
rownames(RNA.TPM.df) <- RNA.TPM.df$gene_name;RNA.TPM.df<- RNA.TPM.df[,-c(1,which(colnames(RNA.TPM.df)=='gene_name'))]

###### Figure 4A
########
TF_list = read_xlsx("./data/1-s2.0-S0092867418301065-mmc2.xlsx",sheet = 2,na = 'NA')
TF_list <- c(TF_list[TF_list$`Is TF?`=='Yes',2])
TF_list <- na.omit(TF_list$...2)

subs_ARhigh = which( Nelson.samples.PCA$Subtypes == "AR+/NE-" | Nelson.samples.PCA$Subtypes == "AR+/NE+" )
subs_ARlow = which( Nelson.samples.PCA$Subtypes != "AR+/NE-" & Nelson.samples.PCA$Subtypes != "AR+/NE+" )
subs_DN = which( Nelson.samples.PCA$Subtypes == "AR-/NE-" )
subs_adeno = which( Nelson.samples.PCA$Subtypes == "AR+/NE-" )
subs_Nepc = which( Nelson.samples.PCA$Subtypes == "AR-/NE+" )
subtype_list = c("AR-/NE-", "AR-/NE+", "ARL/NE-", "AR+/NE-", "AR+/NE+")
M = data.matrix( log2(1+RNA.TPM.df) )
M = M[ rowMedians(M)>0.5 | rowMax(M)>3,]
match.idx(Nelson.samples.PCA$samples, dimnames(M)[[2]] )

M_tf = M[match.idx( dimnames(M)[[1]], TF_list )$idx.A,]

# differential expression
low_logp = 3

fc_NEPC = -3
logp_NEPC = 8

fc_DP = -2
logp_DP = 3

fc_adeno = -1
logp_adeno = low_logp

fc_DN = -1.5
logp_DN = low_logp

fc_low = -2
logp_low = low_logp

ts_NEPC = rowttests(M_tf, factor(Nelson.samples.PCA$Subtypes=="AR-/NE+"))
ts_NEPC$logp = -1*log10(ts_NEPC$p.value)
ts_NEPC$logp_cor = -1*log10( p.adjust( ts_NEPC$p.value, method="BH"))
symbols_de_NEPC = dimnames(ts_NEPC)[[1]][ ts_NEPC$dm<(fc_NEPC) & ts_NEPC$logp_cor>logp_NEPC ]

ts_DP = rowttests(M_tf, factor(Nelson.samples.PCA$Subtypes=="AR+/NE+"))
ts_DP$logp = -1*log10(ts_DP$p.value)
ts_DP$logp_cor = -1*log10( p.adjust( ts_DP$p.value, method="BH"))
symbols_de_DP = dimnames(ts_DP)[[1]][ ts_DP$dm<(fc_DP) & ts_DP$logp_cor>logp_DP ]

ts_adeno = rowttests(M_tf, factor(Nelson.samples.PCA$Subtypes=="AR+/NE-"))
ts_adeno$logp = -1*log10(ts_adeno$p.value)
ts_adeno$logp_cor = -1*log10( p.adjust( ts_adeno$p.value, method="BH"))
symbols_de_adeno = dimnames(ts_adeno)[[1]][ ts_adeno$dm<(fc_adeno) & ts_adeno$logp_cor>logp_adeno ]

ts_DN = rowttests(M_tf, factor(Nelson.samples.PCA$Subtypes=="AR-/NE-"))
ts_DN$logp = -1*log10(ts_DN$p.value)
ts_DN$logp_cor = -1*log10( p.adjust( ts_DN$p.value, method="BH"))
symbols_de_DN = dimnames(ts_DN)[[1]][ ts_DN$dm<(fc_DN) & ts_DN$logp_cor>logp_DN ]

ts_low = rowttests(M_tf, factor(Nelson.samples.PCA$Subtypes=="ARL/NE-"))
ts_low$logp = -1*log10(ts_low$p.value)
ts_low$logp_cor = -1*log10( p.adjust( ts_low$p.value, method="BH"))
symbols_de_low = dimnames(ts_low)[[1]][ ts_low$dm<(fc_low) & ts_low$logp_cor>logp_low ]

idx = which( Nelson.samples.PCA$Subtypes=="AR-/NE-" | Nelson.samples.PCA$Subtypes=="AR-/NE+" )
ts_DN_vs_NEPC = rowttests(M[,idx], factor(Nelson.samples.PCA$Subtypes[idx]=="AR-/NE-"))
symbols_de_DN_vs_NEPC = dimnames(ts_DN_vs_NEPC)[[1]][ which( abs(ts_DN_vs_NEPC$dm)>(0.5) ) ]

st_order = rep(0, length(Nelson.samples.PCA$Subtypes))
st_order[Nelson.samples.PCA$Subtypes=="AR-/NE+"] = 1
st_order[Nelson.samples.PCA$Subtypes=="ARL/NE-"] = 2
st_order[Nelson.samples.PCA$Subtypes=="AR+/NE+"] = 4
st_order[Nelson.samples.PCA$Subtypes=="AR+/NE-"] = 3

gene_list = c(symbols_de_DN, symbols_de_NEPC, symbols_de_low, symbols_de_DP, symbols_de_adeno)
gene_source = c(rep(0, length(symbols_de_DN)),
                rep(1, length(symbols_de_NEPC)),
                rep(2, length(symbols_de_low)),
                rep(3, length(symbols_de_DP)),
                rep(4, length(symbols_de_adeno)))

ord = order(gene_list)
gene_list = gene_list[ord]
gene_source = gene_source[ord]

ord = order( st_order )
M_std = (M_tf-rowMeans(M_tf)) / rowSds(M_tf)
M_std = M_std[gene_list,ord]
col_fun = colorRamp2(c(-4,-2,0,2,4), c("darkblue","royalblue","white","tomato","tomato4"))

set.seed(5)
Figure4A_heatmap <- draw(Heatmap(matrix = M_std,name = 'TF',border = T,
                            cluster_rows = F,
                            cluster_columns = F,heatmap_legend_ = list(direction='horizontal',labels_gp=gpar(fontsize = 8),legend_height=unit(0.8,'mm'),title_gp = gpar(fontsize = 9,fontface='bold')),
                            column_title = c("AR-/NE-","AR-/NE+","ARL/NE-","AR+/NE-","AR+/NE+"),
                            column_split = st_order[ord],
                            row_split=gene_source,row_title = NULL,
                            column_dend_reorder = T,
                            col = col_fun,
                            rect_gp = gpar(col = "white", lwd = 0.2),show_heatmap_legend = T,show_parent_dend_line = F,
                            show_column_names = F,cluster_column_slices = F,
                            row_names_gp = gpar(fontsize=6),cluster_row_slices = F,row_dend_reorder = F,
                            column_title_gp = gpar(fontsize=7,fontface='bold'),column_gap = unit(1.2,'mm')
),
heatmap_legend_side='bottom',annotation_legend_side='bottom')
ht_opt(RESET = TRUE)
########
# Figure 4B
#######
### loading files

ARPC.hyper.files <- na.omit(mixedsort(list.files("./data/HOMER/hyper/ARPC/knownResults/",pattern = ".motif",full.names = T)))
ARLPC.hyper.files <- na.omit(mixedsort(list.files("./data/HOMER/hyper/ARLPC/knownResults/",pattern = ".motif",full.names = T)))
DNPC.hyper.files <- na.omit(mixedsort(list.files("./data/HOMER/hyper/DNPC/knownResults/",pattern = ".motif",full.names = T)))
SCNPC.hyper.files <- na.omit(mixedsort(list.files("./data/HOMER/hyper/SCNPC/knownResults/",pattern = ".motif",full.names = T)))
Amp.hyper.files <- na.omit(mixedsort(list.files("./data/HOMER/hyper/Amp/knownResults/",pattern = ".motif",full.names = T)))

ARPC.hypo.files <- na.omit(mixedsort(list.files("./data/HOMER/hypo/ARPC/knownResults/",pattern = ".motif",full.names = T)))
ARLPC.hypo.files <- na.omit(mixedsort(list.files("./data/HOMER/hypo/ARLPC/knownResults/",pattern = ".motif",full.names = T)))
DNPC.hypo.files <- na.omit(mixedsort(list.files("./data/HOMER/hypo/DNPC/knownResults/",pattern = ".motif",full.names = T)))
SCNPC.hypo.files <- na.omit(mixedsort(list.files("./data/HOMER/hypo/SCNPC/knownResults/",pattern = ".motif",full.names = T)))
Amp.hypo.files <- na.omit(mixedsort(list.files("./data/HOMER/hypo/Amp/knownResults/",pattern = ".motif",full.names = T)))


##### AR+/NE-
Motifs.ARPC.hyper.df <- as.data.frame(foreach(i=ARPC.hyper.files,.combine='rbind',.errorhandling = 'remove') %do% {
    n <- toupper(read_homer(file = i)@name)
    p <- as.numeric(read_homer(file = i)@pval)
    logOdds <- as.numeric(read_homer(file = i)@extrainfo)
    logP <- as.numeric(fread(i,header = F,fill = T)[1,4]*-1)
    stat <- logP*logOdds
    df <- data.frame(subtype='AR+/NE-',motif=n,p.val=p,logOdds=logOdds,logP=logP,stat=stat)
    df <- df[!duplicated(df$motif),]
    df})
Motifs.ARPC.hyper.df <- Motifs.ARPC.hyper.df[order(Motifs.ARPC.hyper.df$logP,decreasing = T),]
Motifs.ARPC.hyper.df$rank <- 1:nrow(Motifs.ARPC.hyper.df)

Motifs.ARPC.hypo.df <- as.data.frame(foreach(i=ARPC.hypo.files,.combine='rbind',.errorhandling = 'remove') %do% {
    n <- toupper(read_homer(file = i)@name)
    p <- as.numeric(read_homer(file = i)@pval)
    logOdds <- as.numeric(read_homer(file = i)@extrainfo)
    logP <- as.numeric(fread(i,header = F,fill = T)[1,4]*-1)
    stat <- logP*logOdds
    df <- data.frame(subtype='AR+/NE-',motif=n,p.val=p,logOdds=logOdds,logP=logP,stat=stat)
    df})
Motifs.ARPC.hypo.df <- Motifs.ARPC.hypo.df[order(Motifs.ARPC.hypo.df$logP,decreasing = T),]
Motifs.ARPC.hypo.df$rank <- 1:nrow(Motifs.ARPC.hypo.df)


##### ARL/NE-
Motifs.ARLPC.hyper.df <- as.data.frame(foreach(i=ARLPC.hyper.files,.combine='rbind',.errorhandling = 'remove') %do% {
    n <- toupper(read_homer(file = i)@name)
    p <- as.numeric(read_homer(file = i)@pval)
    logOdds <- as.numeric(read_homer(file = i)@extrainfo)
    logP <- as.numeric(fread(i,header = F,fill = T)[1,4]*-1)
    stat <- logP*logOdds
    df <- data.frame(subtype='ARL/NE-',motif=n,p.val=p,logOdds=logOdds,logP=logP,stat=stat)
    df <- df[!duplicated(df$motif),]
    df})
Motifs.ARLPC.hyper.df <- Motifs.ARLPC.hyper.df[order(Motifs.ARLPC.hyper.df$logP,decreasing = T),]
Motifs.ARLPC.hyper.df$rank <- 1:nrow(Motifs.ARLPC.hyper.df)

Motifs.ARLPC.hypo.df <- as.data.frame(foreach(i=ARLPC.hypo.files,.combine='rbind',.errorhandling = 'remove') %do% {
    n <- toupper(read_homer(file = i)@name)
    p <- as.numeric(read_homer(file = i)@pval)
    logOdds <- as.numeric(read_homer(file = i)@extrainfo)
    logP <- as.numeric(fread(i,header = F,fill = T)[1,4]*-1)
    stat <- logP*logOdds
    df <- data.frame(subtype='ARL/NE-',motif=n,p.val=p,logOdds=logOdds,logP=logP,stat=stat)
    df})
Motifs.ARLPC.hypo.df <- Motifs.ARLPC.hypo.df[order(Motifs.ARLPC.hypo.df$logP,decreasing = T),]
Motifs.ARLPC.hypo.df$rank <- 1:nrow(Motifs.ARLPC.hypo.df)


##### AR-/NE-
Motifs.DNPC.hyper.df <- as.data.frame(foreach(i=DNPC.hyper.files,.combine='rbind',.errorhandling = 'remove') %do% {
    n <- toupper(read_homer(file = i)@name)
    p <- as.numeric(read_homer(file = i)@pval)
    logOdds <- as.numeric(read_homer(file = i)@extrainfo)
    logP <- as.numeric(fread(i,header = F,fill = T)[1,4]*-1)
    stat <- logP*logOdds
    df <- data.frame(subtype='AR-/NE-',motif=n,p.val=p,logOdds=logOdds,logP=logP,stat=stat)
    df <- df[!duplicated(df$motif),]
    df})
Motifs.DNPC.hyper.df <- Motifs.DNPC.hyper.df[order(Motifs.DNPC.hyper.df$logP,decreasing = T),]
Motifs.DNPC.hyper.df$rank <- 1:nrow(Motifs.DNPC.hyper.df)

Motifs.DNPC.hypo.df <- as.data.frame(foreach(i=DNPC.hypo.files,.combine='rbind',.errorhandling = 'remove') %do% {
    n <- toupper(read_homer(file = i)@name)
    p <- as.numeric(read_homer(file = i)@pval)
    logOdds <- as.numeric(read_homer(file = i)@extrainfo)
    logP <- as.numeric(fread(i,header = F,fill = T)[1,4]*-1)
    stat <- logP*logOdds
    df <- data.frame(subtype='AR-/NE-',motif=n,p.val=p,logOdds=logOdds,logP=logP,stat=stat)
    df})
Motifs.DNPC.hypo.df <- Motifs.DNPC.hypo.df[order(Motifs.DNPC.hypo.df$logP,decreasing = T),]
Motifs.DNPC.hypo.df$rank <- 1:nrow(Motifs.DNPC.hypo.df)



##### AR-/NE+
Motifs.SCNPC.hyper.df <- as.data.frame(foreach(i=SCNPC.hyper.files,.combine='rbind',.errorhandling = 'remove') %do% {
    n <- toupper(read_homer(file = i)@name)
    p <- as.numeric(read_homer(file = i)@pval)
    logOdds <- as.numeric(read_homer(file = i)@extrainfo)
    logP <- as.numeric(fread(i,header = F,fill = T)[1,4]*-1)
    stat <- logP*logOdds
    df <- data.frame(subtype='AR-/NE+',motif=n,p.val=p,logOdds=logOdds,logP=logP,stat=stat)
    df <- df[!duplicated(df$motif),]
    df})
Motifs.SCNPC.hyper.df <- Motifs.SCNPC.hyper.df[order(Motifs.SCNPC.hyper.df$logP,decreasing = T),]
Motifs.SCNPC.hyper.df$rank <- 1:nrow(Motifs.SCNPC.hyper.df)

Motifs.SCNPC.hypo.df <- as.data.frame(foreach(i=SCNPC.hypo.files,.combine='rbind',.errorhandling = 'remove') %do% {
    n <- toupper(read_homer(file = i)@name)
    p <- as.numeric(read_homer(file = i)@pval)
    logOdds <- as.numeric(read_homer(file = i)@extrainfo)
    logP <- as.numeric(fread(i,header = F,fill = T)[1,4]*-1)
    stat <- logP*logOdds
    df <- data.frame(subtype='AR-/NE+',motif=n,p.val=p,logOdds=logOdds,logP=logP,stat=stat)
    df})
Motifs.SCNPC.hypo.df <- Motifs.SCNPC.hypo.df[order(Motifs.SCNPC.hypo.df$logP,decreasing = T),]
Motifs.SCNPC.hypo.df$rank <- 1:nrow(Motifs.SCNPC.hypo.df)


##### AR+/NE+
Motifs.Amp.hyper.df <- as.data.frame(foreach(i=Amp.hyper.files,.combine='rbind',.errorhandling = 'remove') %do% {
    n <- toupper(read_homer(file = i)@name)
    p <- as.numeric(read_homer(file = i)@pval)
    logOdds <- as.numeric(read_homer(file = i)@extrainfo)
    logP <- as.numeric(fread(i,header = F,fill = T)[1,4]*-1)
    stat <- logP*logOdds
    df <- data.frame(subtype='AR+/NE+',motif=n,p.val=p,logOdds=logOdds,logP=logP,stat=stat)
    df <- df[!duplicated(df$motif),]
    df})
Motifs.Amp.hyper.df <- Motifs.Amp.hyper.df[order(Motifs.Amp.hyper.df$logP,decreasing = T),]
Motifs.Amp.hyper.df$rank <- 1:nrow(Motifs.Amp.hyper.df)

Motifs.Amp.hypo.df <- as.data.frame(foreach(i=Amp.hypo.files,.combine='rbind',.errorhandling = 'remove') %do% {
    n <- toupper(read_homer(file = i)@name)
    p <- as.numeric(read_homer(file = i)@pval)
    logOdds <- as.numeric(fread(i,header = F,fill = T)[1,3])
    logP <- as.numeric(fread(i,header = F,fill = T)[1,4]*-1)
    stat <- logP*logOdds
    df <- data.frame(subtype='AR+/NE+',motif=n,p.val=p,logOdds=logOdds,logP=logP,stat=stat)
    df})
Motifs.Amp.hypo.df <- Motifs.Amp.hypo.df[order(Motifs.Amp.hypo.df$logP,decreasing = T),]
Motifs.Amp.hypo.df$rank <- 1:nrow(Motifs.Amp.hypo.df)


top30.hyper.df <- unique(na.omit(c(Motifs.ARPC.hyper.df$motif[1:30],
                                   Motifs.ARLPC.hyper.df$motif[1:30],
                                   Motifs.DNPC.hyper.df$motif[1:30],
                                   Motifs.SCNPC.hyper.df$motif[1:30],
                                   Motifs.Amp.hyper.df$motif[1:30])))

top20.hypo.df <- unique(na.omit(c(Motifs.DNPC.hypo.df$motif[1:20],
                                  Motifs.SCNPC.hypo.df$motif[1:20],
                                  Motifs.ARLPC.hypo.df$motif[1:20],
                                  Motifs.ARPC.hypo.df$motif[1:20],
                                  Motifs.Amp.hypo.df$motif[1:20])))


hypo.hits <- rbind(Motifs.ARPC.hypo.df[Motifs.ARPC.hypo.df$motif%in%top20.hypo.df,],
                   Motifs.ARLPC.hypo.df[Motifs.ARLPC.hypo.df$motif%in%top20.hypo.df,],
                   Motifs.DNPC.hypo.df[Motifs.DNPC.hypo.df$motif%in%top20.hypo.df,],
                   Motifs.SCNPC.hypo.df[Motifs.SCNPC.hypo.df$motif%in%top20.hypo.df,],
                   Motifs.Amp.hypo.df[Motifs.Amp.hypo.df$motif%in%top20.hypo.df,])

hypo.hits$subtype <- factor(hypo.hits$subtype,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE-','AR+/NE+'))
hypo.hits$motif <- factor(hypo.hits$motif,levels=top20.hypo.df)
hypo.hits <- hypo.hits[order(hypo.hits$motif),]


hypo.df <- data.frame(motif=unique(c(hypo.hits$motif)))
hypo.df <- merge(hypo.df,hypo.hits[hypo.hits$subtype=='AR-/NE-',c('motif','rank')],by='motif',all.x=T,sort=F)
hypo.df <- merge(hypo.df,hypo.hits[hypo.hits$subtype=='AR-/NE+',c('motif','rank')],by='motif',all.x=T,sort=F)
hypo.df <- merge(hypo.df,hypo.hits[hypo.hits$subtype=='ARL/NE-',c('motif','rank')],by='motif',all.x=T,sort=F)
hypo.df <- merge(hypo.df,hypo.hits[hypo.hits$subtype=='AR+/NE-',c('motif','rank')],by='motif',all.x=T,sort=F)
hypo.df <- merge(hypo.df,hypo.hits[hypo.hits$subtype=='AR+/NE+',c('motif','rank')],by='motif',all.x=T,sort=F)

colnames(hypo.df) <- c('motifs','AR-/NE-','AR-/NE+','ARL/NE-',
                       'AR+/NE-','AR+/NE+')

hypo.df <- hypo.df[!duplicated(hypo.df$motifs),]
hypo.df <- hypo.df[order(hypo.df$motifs),]
rownames(hypo.df) <- hypo.df$motifs
set.seed(100)
col_fun.hypo =colorRamp2(c(400,200,100,50,10,1), c("seashell","pink2","salmon","red3","tomato3","red4"))


Figure4B_heatmap <- draw(Heatmap(matrix = t(hypo.df[,-1]),name = 'Rank',border = T,
                             heatmap_legend_param = list(direction='horizontal',at=c(1,10,50,100,200,300,400),labels=c(1,10,50,100,200,300,400),labels_gp=gpar(fontsize = 12),title_gp = gpar(fontsize = 17,fontface='bold')),
                             cluster_columns = F,
                             col = col_fun.hypo,
                             cluster_rows = F,na_col = 'white',row_dend_gp = gpar(col = "black", lwd = 1.5),
                             rect_gp = gpar(col = "black", lwd = 0.5),show_row_names = T,column_names_rot = 45,
                             show_heatmap_legend = T,
                             show_column_names = T,row_names_side = 'left',row_title_gp = gpar(fontsize=15,fontface='bold'),
                             width = nrow(hypo.df)*unit(4.5, "mm"),column_title = '',
                             height = ncol(hypo.df)*unit(4.5, "mm"),
                             column_names_gp = gpar(fontsize=10),column_title_gp = gpar(fontsize=15,fontface='bold'),
                             row_names_gp = gpar(fontsize=12,fontface='bold')),heatmap_legend_side='top',annotation_legend_side='bottom')

########

###### Figure 4C

########
Motifs.SCNPC.hypo.df$labels <- ifelse(Motifs.SCNPC.hypo.df$motif %in% c('ASCL1','SNAIL1','NEUROD1','NEUROG2','KLF5','KLF1','KLF3','SOX3','SOX2','SOX21','SOX10','OLIG2','BCL11B'),
                                      Motifs.SCNPC.hypo.df$motif,"")
Motifs.SCNPC.hypo.df$color <- ifelse(Motifs.SCNPC.hypo.df$motif=='KLF5','black',
                                     ifelse(Motifs.SCNPC.hypo.df$rank<=15,subtype.colors[subtype.colors$subtypes=='AR-/NE+',]$colors,'gray'))
Motifs.SCNPC.hypo.df$size <- ifelse(Motifs.SCNPC.hypo.df$motif=='KLF5',3,2)
Motifs.SCNPC.hypo.df$alpha <- ifelse(Motifs.SCNPC.hypo.df$labels=="",0.2,1)

tmp.SCNPC <- Motifs.SCNPC.hypo.df[Motifs.SCNPC.hypo.df$rank<=300,]
tmp.SCNPC <- tmp.SCNPC[!duplicated(tmp.SCNPC$motif),]
tmp.SCNPC$motif <- factor(tmp.SCNPC$motif,levels = unique(tmp.SCNPC$motif))

Figure4C_NEPC <- ggplot(tmp.SCNPC,aes(x=motif,y=rank,label=labels,color=color))+
    geom_point(size=tmp.SCNPC$size/2,alpha=tmp.SCNPC$alpha)+
    geom_hline(yintercept = 20,linetype='dashed',color='red',alpha=0.6)+
    geom_text_repel(max.overlaps = 35,nudge_y = 3,nudge_x = 1,size=tmp.SCNPC$size*2,fontface='bold',segment.size=0.3)+
    scale_color_manual(values = tmp.SCNPC$color,breaks = tmp.SCNPC$color)+
    scale_y_continuous(name='Rank',breaks=c(1,20,100,200,257,300))+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    labs(x='Motifs',y='Rank',title='AR-/NE+')+
    theme(axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=c(15,15,15,15,22,15),face = c('plain','plain','plain','plain','bold','plain')),title = element_text(face = "bold"),
          legend.position = 'none',strip.text = element_text(face='bold',size=12.5),
          plot.title=element_text(size=26,
                                  hjust = 0.5),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())



Motifs.DNPC.hypo.df$labels <- ifelse(Motifs.DNPC.hypo.df$motif %in% c('KLF1','KLF5','KLF3','SOX3','SOX2','SOX21','SOX10'),
                                     Motifs.DNPC.hypo.df$motif,"")
Motifs.DNPC.hypo.df$color <- ifelse(Motifs.DNPC.hypo.df$motif=='KLF5','black',
                                    ifelse(Motifs.DNPC.hypo.df$rank<=20,subtype.colors[subtype.colors$subtypes=='AR-/NE-',]$colors,'gray'))
Motifs.DNPC.hypo.df$size <- ifelse(Motifs.DNPC.hypo.df$motif=='KLF5',3,2)
Motifs.DNPC.hypo.df$alpha <- ifelse(Motifs.DNPC.hypo.df$labels=="",0.2,1)

tmp.DNPC <- Motifs.DNPC.hypo.df[Motifs.DNPC.hypo.df$rank<=200,]
tmp.DNPC <- tmp.DNPC[!duplicated(tmp.DNPC$motif),]
tmp.DNPC$motif <- factor(tmp.DNPC$motif,levels = unique(tmp.DNPC$motif))

Figure4C_DNPC <- ggplot(tmp.DNPC,aes(x=motif,y=rank,label=labels,color=color))+
    geom_point(size=tmp.DNPC$size/2,alpha=tmp.DNPC$alpha)+
    geom_hline(yintercept = 20,linetype='dashed',color='red',alpha=0.6)+
    geom_text_repel(max.overlaps = 30,nudge_y = 3,nudge_x = 1,size=tmp.DNPC$size*2,fontface='bold',segment.size=0.3)+
    scale_color_manual(values = tmp.DNPC$color,breaks = tmp.DNPC$color)+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    labs(x='Motifs',y='Rank',title='AR-/NE-')+
    scale_y_continuous(name='Rank',breaks=c(1,20,100,200))+
    theme(axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          plot.title=element_text(size=26,
                                  hjust = 0.5),
          axis.text.y = element_text(size=c(22,15,15,15),face = c('bold','plain','plain','plain')),title = element_text(face = "bold"),
          legend.position = 'none',strip.text = element_text(face='bold',size=12.5),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())


Figure4C_rankplot <- plot_grid(Figure4C_DNPC,Figure4C_NEPC+theme(axis.title.y= element_blank()),align = 'hv',ncol = 2)


### Figure 4D

KLF5.motif <- fread("./data/HOMER/klf5_motif_annotation.txt",data.table = F)
KLF5.motif <- KLF5.motif[KLF5.motif$`PeakID (cmd=annotatePeaks.pl /data1/projects/WCDT_smallcell_2021/5mc/HOMER/input/background.DMR.hypo.bed hg38 -m /data1/home/alundberg/Softwares/homer/data/knownTFs/motifs/klf5.motif -noann -size given -cpu 20 -mbed /data1/projects/WCDT_smallcell_2021/5mc/HOMER/output/Annotation/hypo/klf5_motif_annotation.bed)`!="",]

set.seed(100)
KLF5.motif.hypo <- enricher(KLF5.motif$`Gene Name`,pAdjustMethod = "BH",TERM2GENE = pid2gene,TERM2NAME = pid2name)
KLF5.motif.hypo.res <- KLF5.motif.hypo@result[KLF5.motif.hypo@result$p.adjust <= 0.1,]

KLF5.motif.hypo.res$Description <- gsub("_"," ",KLF5.motif.hypo.res$Description)
KLF5.motif.hypo.res$Description <- gsub("[.]"," ",KLF5.motif.hypo.res$Description)
KLF5.motif.hypo.res$Description <- str_to_title(KLF5.motif.hypo.res$Description)
KLF5.motif.hypo.res$Description <- gsub("Wp ","",KLF5.motif.hypo.res$Description)
KLF5.motif.hypo.res$Description <- gsub("Hallmark ","",KLF5.motif.hypo.res$Description)

KLF5.motif.hypo.res$logFDR <- -log10(KLF5.motif.hypo.res$p.adjust)

Figure4D_GSEA <- ggbarplot(KLF5.motif.hypo.res,x = 'Description',y='logFDR',color ='black',fill='p.adjust',
                             xlab = '',ylab = "-log10 (FDR)",title='')+
    geom_hline(yintercept = 1.3,linetype='dashed',color='gray20',lwd=1)+
    rotate()+
    scale_y_reverse(breaks = c(0,1,2,3),labels = c("0",1,2,3))+
    scale_fill_continuous_tableau(palette = "Blue-Teal")+
    scale_x_discrete(position = 'top',limits=rev)+
    theme_bw(base_rect_size = 2,base_size = 15)+
    theme(axis.title.y = element_text(size=15,face = 'bold'),axis.title.x = element_text(size=17,face = 'bold'),
          axis.text.y = element_text(angle=-0,vjust = 0,size=20,face = 'bold'),title = element_text(face='bold'),
          axis.text.x = element_text(size=20),strip.text = element_text(face='bold',size=15),
          strip.background =element_rect(fill='white',color='white'),legend.position = 'none',
          panel.grid = element_blank())

########

### Figure 4E

########
KLF5.compare.df <- as.data.frame(t(log2(RNA.TPM.df+1)[rownames(RNA.TPM.df)%in%toupper(c('KLF5','KRT5','KRT8','RB1','CCNB2')),]))
KLF5.compare.df$samples <- rownames(KLF5.compare.df)
KLF5.compare.df <- left_join(Nelson.samples.PCA[,-3],KLF5.compare.df)
KLF5.compare.df$Subtypes <- factor(KLF5.compare.df$Subtypes,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE-','AR+/NE+'))
my_comparisons <- list(c("AR-/NE+", "AR-/NE-"),c("AR-/NE+", "AR+/NE-"),c("AR-/NE+", "ARL/NE-"),c("AR-/NE+", "AR+/NE+"))

melt.df <- melt(KLF5.compare.df,id.vars=c('samples','Subtypes','KLF5'))

Figure4E.plot <- ggscatter(melt.df,x='KLF5',y='value',xlab='KLF5',ylab='',color='black',fill='Subtypes',shape = 21,size=2.5,conf.int = T,
                           add = "reg.line",
                           add.params = list(fill = "lightgray",size=0.2))+facet_wrap(~Subtypes,ncol = 3)+
    facet_grid(~variable~Subtypes,switch='y')+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+ stat_cor(method = "spearman", label.x = 0.2, label.y = 9,p.accuracy = 0.001, r.accuracy = 0.01,cex=3.7)+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=14),axis.title.x = element_text(size=14,face = 'bold'),
          axis.text.x = element_text(size=12.5),
          axis.text.y = element_text(size=12.5),strip.placement = "outside",
          legend.position = 'right',strip.text = element_text(face='bold',size=14),
          legend.text = element_text(size = 15),legend.title = element_text(face = 'bold',size=17),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())

########

#### Figure 4F
########
WGS.available <- SO_WCDT_WGS[names(SO_WCDT_WGS)%in%Nelson.samples.PCA$samples]

TPM.log10.df <- log10(RNA.TPM.df+1)

RB1.KLF5.df <- as.data.frame(foreach(i=names(WGS.available),.combine = 'rbind') %do% {
    KLF5.cn <- WGS.available[[i]]$CNA_genes['KLF5',2]
    KLF5.exprs <- TPM.log10.df['KLF5',i]
    RB1.cn <- WGS.available[[i]]$CNA_genes['RB1',2]
    RB1.exprs <- TPM.log10.df['RB1',i]
    c(KLF5.cn=KLF5.cn,KLF5.exprs=KLF5.exprs,RB1.cn=RB1.cn,RB1.exprs=RB1.exprs)})
RB1.KLF5.df$samples <- names(WGS.available)
RB1.KLF5.df <- merge(RB1.KLF5.df,Nelson.samples.PCA,by='samples')


RB1.RB1.cn.sp <- ggscatter(RB1.KLF5.df[RB1.KLF5.df$Subtypes%in%c('AR-/NE-','AR+/NE-'),],x='RB1.cn',y='RB1.exprs',xlab='RB1 copies',ylab='RB1\nexpression',color='black',fill='Subtypes',shape = 21,size=2.5,conf.int = T,
                           add.params = list(fill = "lightgray",size=0.2))+facet_wrap(~Subtypes,ncol = 3)+
    facet_wrap(~Subtypes,scales = 'free')+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=14,face = 'bold'),axis.title.x = element_text(size=14,face = 'bold'),
          axis.text.x = element_text(size=12.5),
          axis.text.y = element_text(size=12.5),strip.placement = "outside",
          legend.position = 'none',strip.text = element_text(face='bold',size=17),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())

KLF5.KLF5.cn.sp <- ggscatter(RB1.KLF5.df[RB1.KLF5.df$Subtypes%in%c('AR-/NE-','AR+/NE-'),],x='KLF5.cn',y='KLF5.exprs',xlab='KLF5 copies',ylab='KLF5\nexpression',color='black',fill='Subtypes',shape = 21,size=2.5,conf.int = T,
                             add.params = list(fill = "lightgray",size=0.2))+facet_wrap(~Subtypes,ncol = 3)+
    facet_wrap(~Subtypes,scales = 'free',labeller = as_labeller(c("AR-/NE-"="","AR+/NE-"="")))+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=14,face = 'bold'),axis.title.x = element_text(size=14,face = 'bold'),
          axis.text.x = element_text(size=12.5),
          axis.text.y = element_text(size=12.5),strip.placement = "outside",
          legend.position = 'none',strip.text = element_text(face='bold',size=17),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())

KLF5.RB1.cn.sp <- ggscatter(RB1.KLF5.df[RB1.KLF5.df$Subtypes%in%c('AR-/NE-','AR+/NE-'),],x='RB1.cn',y='KLF5.exprs',xlab='RB1 copies',ylab='KLF5\nexpression',color='black',fill='Subtypes',shape = 21,size=2.5,conf.int = T,
                            add.params = list(fill = "lightgray",size=0.2))+facet_wrap(~Subtypes,ncol = 3)+
    facet_wrap(~Subtypes,scales = 'free',labeller = as_labeller(c("AR-/NE-"="","AR+/NE-"="")))+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=14,face = 'bold'),axis.title.x = element_text(size=14,face = 'bold'),
          axis.text.x = element_text(size=12.5),
          axis.text.y = element_text(size=12.5),strip.placement = "outside",
          legend.position = 'none',strip.text = element_text(face='bold',size=17),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())

Figure4F.plot <- ggarrange(RB1.RB1.cn.sp,KLF5.KLF5.cn.sp,KLF5.RB1.cn.sp,nrow = 3,align = 'hv')
########


##### Supplementary Figure 8

hyper.hits <- rbind(Motifs.ARPC.hyper.df[Motifs.ARPC.hyper.df$motif%in%top30.hyper.df,],
                    Motifs.ARLPC.hyper.df[Motifs.ARLPC.hyper.df$motif%in%top30.hyper.df,],
                    Motifs.DNPC.hyper.df[Motifs.DNPC.hyper.df$motif%in%top30.hyper.df,],
                    Motifs.SCNPC.hyper.df[Motifs.SCNPC.hyper.df$motif%in%top30.hyper.df,],
                    Motifs.Amp.hyper.df[Motifs.Amp.hyper.df$motif%in%top30.hyper.df,])

hyper.hits$subtype <- factor(hyper.hits$subtype,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE-','AR+/NE+'))
hyper.hits$motif <- factor(hyper.hits$motif,levels=top30.hyper.df)
hyper.hits <- hyper.hits[order(hyper.hits$motif),]


hyper.df <- data.frame(motif=unique(c(hyper.hits$motif)))
hyper.df <- merge(hyper.df,hyper.hits[hyper.hits$subtype=='AR-/NE-',c('motif','rank')],by='motif',all.x=T,sort=F)
hyper.df <- merge(hyper.df,hyper.hits[hyper.hits$subtype=='AR-/NE+',c('motif','rank')],by='motif',all.x=T,sort=F)
hyper.df <- merge(hyper.df,hyper.hits[hyper.hits$subtype=='ARL/NE-',c('motif','rank')],by='motif',all.x=T,sort=F)
hyper.df <- merge(hyper.df,hyper.hits[hyper.hits$subtype=='AR+/NE-',c('motif','rank')],by='motif',all.x=T,sort=F)
hyper.df <- merge(hyper.df,hyper.hits[hyper.hits$subtype=='AR+/NE+',c('motif','rank')],by='motif',all.x=T,sort=F)

colnames(hyper.df) <- c('motifs','AR-/NE-','AR-/NE+','ARL/NE-',
                        'AR+/NE-','AR+/NE+')

hyper.df <- hyper.df[!duplicated(hyper.df$motifs),]
hyper.df <- hyper.df[order(hyper.df$motifs),]
rownames(hyper.df) <- hyper.df$motifs
set.seed(100)

col_fun.hyper =colorRamp2(c(350,100,10,1),c("slategray1","blue","blue2","darkblue"))


SFigure8_heatmap <- draw(Heatmap(matrix = t(hyper.df[,-1]),name = 'Rank',border = T,
                              heatmap_legend_param = list(direction='horizontal',at=c(1,10,50,100,200,300,400),labels=c(1,10,50,100,200,300,400),labels_gp=gpar(fontsize = 12),title_gp = gpar(fontsize = 17,fontface='bold')),
                              cluster_columns = F,
                              col = col_fun.hyper,
                              cluster_rows = F,na_col = 'white',row_dend_gp = gpar(col = "black", lwd = 1.5),
                              rect_gp = gpar(col = "black", lwd = 0.5),show_row_names = T,column_names_rot = 45,
                              show_heatmap_legend = T,
                              show_column_names = T,row_names_side = 'left',row_title_gp = gpar(fontsize=15,fontface='bold'),
                              width = nrow(hyper.df)*unit(2.7, "mm"),column_title = 'Hypermethylated',
                              height = ncol(hyper.df)*unit(5, "mm"),
                              column_names_gp = gpar(fontsize=6),column_title_gp = gpar(fontsize=15,fontface='bold'),
                              row_names_gp = gpar(fontsize=10,fontface='bold')),heatmap_legend_side='bottom',annotation_legend_side='bottom')





