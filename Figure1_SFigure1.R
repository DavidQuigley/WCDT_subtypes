library(openxlsx)
library(data.table)
library(dplyr)
library(circlize)
library(ggplot2)
library(ggfortify)
library(singscore)
library(ComplexHeatmap)
library(GSVA)
library(matrixStats)
library(plot.matrix)
library(rstatix)
library(msigdbr)
library(stringr)
library(ggpubr)
setwd("~/Box Sync/UCSF/Manuscript/Submissions/GitHub/")
### source functions
source("./functions.R")

#### loading gene lists

### gene panels from labraque et al.
load(file = "./data/NEPC.gene.list.RData")

### RNA-seq data
TPM.WCDT <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes.txt",sep="\t",header = T,data.table = F)
TPM.PROMOTE <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt",sep="\t",header = T,data.table = F)

RNA.TPM.df <-  right_join(TPM.WCDT,TPM.PROMOTE)
RNA.TPM.df <- RNA.TPM.df[!duplicated(RNA.TPM.df$gene_name),]
rownames(RNA.TPM.df) <- RNA.TPM.df$gene_name;RNA.TPM.df<- RNA.TPM.df[,-c(1,which(colnames(RNA.TPM.df)=='gene_name'))]

col_fun = colorRamp2(c(-4,-1,0,1,4), c("darkblue","royalblue","white","gold3","gold4"))
my_colour = list('Binary class' = c('t-SCNC' = 'black','Adeno'='white'),
                 'Treatment' = c('TRUE'='indianred4','FALSE'='salmon','no-data'='white'),
                 'Molecular subtypes'=structure(c(subtype.colors$colors),
                                                names = c(subtype.colors$subtypes)))



### selecting t-SCNC samples using beltran et. al. gene list - binary classification
#####

# subset.tscnc.tpm.data <- RNA.TPM.df[rownames(RNA.TPM.df)%in%NEPC.gene.list$beltran$GENE_ID,colnames(RNA.TPM.df)%in%matched_samples$Match_ID]
# subset.tscnc.tpm.data <- t(scale(t(log2(subset.tscnc.tpm.data+1))))
# set.seed(5)
# beltran.heatmap <- draw(Heatmap(matrix = subset.tscnc.tpm.data,name = 'log2(TPM+1)',border = T,cluster_rows = T,heatmap_legend_param = list(direction='horizontal'),
#                                 cluster_columns = T,col=col_fun,
#                                 #top_annotation = colAnnot.bel,
#                                 row_names_gp = gpar(fontsize=8),column_title = "Nature Med. gene set",
#                                 column_names_gp = gpar(fontsize=6),column_km = 2,clustering_distance_rows  = "pearson",clustering_distance_columns  = "pearson"),heatmap_legend_side='bottom')
# small_cluster.beltran.RNAseq <- sort(colnames(subset.tscnc.tpm.data)[column_order(beltran.heatmap)[[2]]])

# "DTB-003-BL"  "DTB-032-BL"  "DTB-036-BL"  "DTB-040-BL"  "DTB-130-BL"  "DTB-135-PRO" "DTB-205-BL"  "DTB-218-BL"  "DTB-231-BL"
#####

#### Figure 1A
######
tSCNC <- c("DTB-003-BL","DTB-032-BL","DTB-036-BL","DTB-040-BL","DTB-130-BL","DTB-135-PRO",
           "DTB-205-BL","DTB-218-BL","DTB-231-BL")

Nelson.tpm.data <- RNA.TPM.df[rownames(RNA.TPM.df)%in%NEPC.gene.list$Nelson$GENE_ID,]
rownames(Nelson.tpm.data)[9] <- 'TARP' # found the alias name in the cohort, rename it to the original publication name

### order genes based on Nelson list
Nelson.tpm.data <- Nelson.tpm.data[order(match(rownames(Nelson.tpm.data),NEPC.gene.list$Nelson$original_names)),]
colData.Nelson = list(GO=c("NE1"='white',
                           "NE2"='white',"AR"='white',"SQUAM"='white'))

rowData.Nelson <- data.frame(GO=NEPC.gene.list$Nelson[NEPC.gene.list$Nelson$original_names %in% rownames(Nelson.tpm.data),]$class)

Nelson.tpm.data <- log2(Nelson.tpm.data+1);Nelson.tpm.data <- as.matrix(Nelson.tpm.data)
Nelson.tpm.data <- t(scale(t(Nelson.tpm.data)))

ht_opt$TITLE_PADDING = unit(c(15,15), "points")
set.seed(20)
group = kmeans(t(Nelson.tpm.data),nstart = 1,centers = 5,iter.max = 10000)$cluster

# tested with the heatmap
Nelsons_clusters <- list('AR-/NE-' = names(group[group==4]),
                         'AR+/NE+' = names(group[group==2]),
                         'AR+/NE-' = names(group[group==1]),
                         'AR-/NE+' = names(group[group==5]),
                         'ARL/NE-' = names(group[group==3]))

group <- factor(group,labels = c('AR+/NE-','AR+/NE+','ARL/NE-','AR-/NE-','AR-/NE+'))
group <- factor(group,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE-','AR+/NE+'))

table(group)

# group
# AR-/NE- AR-/NE+ ARL/NE- AR+/NE- AR+/NE+
#     13       7      49     132       9

##### PCAs - Supplemental Figure 1B

Nelson.tpm.PCA <- t(Nelson.tpm.data)
Nelson.samples.PCA <- data.frame(samples=rownames(Nelson.tpm.PCA),Subtypes=ifelse(rownames(Nelson.tpm.PCA)%in%Nelsons_clusters[[1]], names(Nelsons_clusters)[1],
                                                                                  ifelse(rownames(Nelson.tpm.PCA)%in%Nelsons_clusters[[2]], names(Nelsons_clusters)[2],
                                                                                         ifelse(rownames(Nelson.tpm.PCA)%in%Nelsons_clusters[[3]], names(Nelsons_clusters)[3],
                                                                                                ifelse(rownames(Nelson.tpm.PCA)%in%Nelsons_clusters[[4]], names(Nelsons_clusters)[4],
                                                                                                       ifelse(rownames(Nelson.tpm.PCA)%in%Nelsons_clusters[[5]], names(Nelsons_clusters)[5],
                                                                                                              NA))))))


Nelson.data.PCA <- merge(Nelson.tpm.PCA,Nelson.samples.PCA,by.x='row.names',by.y='samples',sort=F);rownames(Nelson.data.PCA) <- Nelson.data.PCA$Row.names;Nelson.data.PCA <- Nelson.data.PCA[,-1]
Nelson.samples.PCA$binary <- ifelse(Nelson.samples.PCA$samples %in% tSCNC,'t-SCNC','Adeno')

tmp.clusters <- Nelsons_clusters
tmp.PCA <- Nelson.data.PCA

Nelson.PCA <- prcomp(tmp.PCA[,-c(ncol(tmp.PCA))])

SFigure1B_PCA <- autoplot(Nelson.PCA,data=Nelson.data.PCA,colour='Subtypes',frame=F,size=3,label.size=2,type=23)+
    scale_fill_manual(breaks = c(subtype.colors$subtypes),values = c(subtype.colors$colors))+
    scale_color_manual(breaks = c(subtype.colors$subtypes),values = c(subtype.colors$colors))+
    theme_bw(base_size = 15,base_rect_size = 2)+theme(strip.background.x = element_blank(),panel.grid = element_blank(),
                                                      strip.text.x = element_blank())



### 3D PCA - Supplemental figure 1C

idx.PCA <- factor(Nelson.samples.PCA$Subtypes,levels=c(subtype.colors$subtypes))
PCA.colors <- c(subtype.colors$colors)[idx.PCA]
dev.variable <- Nelson.PCA$sdev^2/sum(Nelson.PCA$sdev^2)


scatterplot3d::scatterplot3d(Nelson.PCA$x[,c(1,3,2)],color =PCA.colors,pch = 16, grid=TRUE, box=FALSE,type = 'p',angle = 25,cex.symbols = 1,
                             xlab =paste0("PC1: ",round(dev.variable[1]*100,2),"%"),
                             ylab =paste0("PC3: ",round(dev.variable[3]*100,2),"%"),
                             zlab =paste0("PC2: ",round(dev.variable[2]*100,2),"%"))



#### adding colSide columns
gex_rank <- rankGenes(as.matrix(RNA.TPM.df[,colnames(Nelson.tpm.data)]))
AR.panel <- NEPC.gene.list$Nelson[NEPC.gene.list$Nelson$class=='AR',]$GENE_ID
NE.panel <- NEPC.gene.list$Nelson[NEPC.gene.list$Nelson$class%in%c('NE1','NE2'),]$GENE_ID
AR.score <- simpleScore(rankData = gex_rank, upSet =AR.panel[-8],downSet = AR.panel[8],knownDirection = TRUE, centerScore = TRUE, subSamples = NULL)[,"TotalScore",drop=F]
NEPC.score <- simpleScore(rankData = gex_rank, upSet =NE.panel,knownDirection = TRUE, centerScore = TRUE, subSamples = NULL)[,"TotalScore",drop=F]

#all(rownames(AR.score)==rownames(NEPC.score))
# TRUE
signatures.df <- merge(AR.score,NEPC.score,by='row.names');signatures.df <- signatures.df[ order(match(signatures.df$Row.names, colnames(Nelson.tpm.data))), ]
colnames(signatures.df) <- c('samples','AR.score','NEPC.score')
signatures.df$binary <- ifelse(signatures.df$samples %in% tSCNC,'t-SCNC','Adeno')
#all(signatures.df$samples==colnames(Nelson.tpm.data))
# TRUE

colAnnot.df = columnAnnotation(
    'AR signature' = anno_simple(signatures.df$AR.score, col = colorRamp2(c(-1, 1), c("white", "black")), na_col = "white", gp = gpar(col = "black",lwd=0.4)),
    'NEPC score' = anno_simple(signatures.df$NEPC.score, col = colorRamp2(c(-0.35, 0.35), c("white", "black")), na_col = "white", gp = gpar(col = "black",lwd=0.4)),
    'Binary class'=signatures.df$binary,
    col=my_colour,annotation_height=unit(c(rep(0.3,3)), "cm"),gp=gpar(col="black",lwd=0.4),show_legend=F,
    annotation_name_gp= gpar(fontsize = 6))


ht_opt$TITLE_PADDING = unit(c(5,5), "points")

set.seed(5)
Figure1A_heatmap <- draw(Heatmap(matrix = Nelson.tpm.data,name = 'zscore of log2(TPM)',border = T,cluster_rows = F,heatmap_legend_param = list(direction='horizontal'),
                                       column_split = factor(group,levels=levels(group)),
                                       column_dend_reorder = T,
                                       col = col_fun,column_dend_height = unit(10,'mm'),column_dend_side = 'bottom',
                                       rect_gp = gpar(col = "black", lwd = 0.4),show_heatmap_legend = F,show_parent_dend_line = F,
                                       top_annotation = colAnnot.df,show_column_dend = F,
                                       show_column_names = F,cluster_column_slices = F,
                                       row_split  = factor(rowData.Nelson$GO,levels=c('NE1','NE2','AR','SQUAM')),row_names_gp = gpar(fontsize=5),cluster_row_slices = F,row_dend_reorder = F,
                                       clustering_distance_columns ="pearson",row_title_gp = gpar(fontsize=6,fontface='bold'),
                                       column_title_gp = gpar(fill = c("#7570B3","#E7298A","#D95F02","#1B9E77","#66A61E"),fontsize=3.9,fontface='bold',border='black',col='white')),
                               heatmap_legend_side='bottom',annotation_legend_side='bottom')
ht_opt(RESET = TRUE)
######

#### Supplementary figure 1A

dendogram.top = columnAnnotation(
    'AR signature' = anno_simple(signatures.df$AR.score, col = colorRamp2(c(-1, 1), c("white", "black")), na_col = "white", gp = gpar(col = "black",lwd=0.4)),
    'NEPC score' = anno_simple(signatures.df$NEPC.score, col = colorRamp2(c(-0.35, 0.35), c("white", "black")), na_col = "white", gp = gpar(col = "black",lwd=0.4)),
    'Binary class'=signatures.df$binary,
    # 'Molecular subtypes'=WCDT_treatment$Subtypes,
    col=my_colour,annotation_height=unit(c(rep(0.6,3)), "cm"),gp=gpar(col="black",lwd=0.4),show_legend=F,
    annotation_name_gp= gpar(fontsize = 10))



ht_opt$TITLE_PADDING = unit(c(20,20), "points")

set.seed(5)

SFigure1A_dendogram <- draw(Heatmap(matrix = Nelson.tpm.data,name = 'zscore of log2(TPM)',border = F,cluster_rows = F,heatmap_legend_param = list(direction='horizontal'),
                             column_split = factor(group,levels=levels(group)),
                             column_dend_reorder = F,
                             col = col_fun,column_dend_height = unit(200,'mm'),column_dend_side = 'top',
                             show_heatmap_legend = F,show_parent_dend_line = F,
                             top_annotation = dendogram.top,show_column_dend = T,
                             show_column_names = F,cluster_column_slices = T,row_names_gp = gpar(fontsize=1,fontface='bold'),
                             clustering_distance_columns ="pearson",row_title_gp = gpar(fontsize=1,fontface='bold'),
                             column_title_gp = gpar(fill = c("#E7298A","#7570B3","#66A61E","#1B9E77","#D95F02"),fontsize=15,fontface='bold',border='black',col='white')),
                     heatmap_legend_side='bottom',annotation_legend_side='bottom')


### Figure 1B
######
Beltran.df <-read.table(file = "./data/nepc_wcm_2016/data_RNA_Seq_mRNA_median_all_sample_Zscores.txt",sep = "\t",header=T)
Beltran.RNA.df <- subset(Beltran.df,subset=Beltran.df$Hugo_Symbol%in%NEPC.gene.list$Nelson$GENE_ID)

rownames(Beltran.RNA.df) <- Beltran.RNA.df$Hugo_Symbol;Beltran.RNA.df <- Beltran.RNA.df[,-c(1,2)]
Beltran.RNA.df <- na.omit(Beltran.RNA.df)
Beltran.RNA.df <- Beltran.RNA.df[order(match(rownames(Beltran.RNA.df),NEPC.gene.list$Nelson$GENE_ID)),]
rownames(Beltran.RNA.df)[16] <- 'TARP' # found the alias name in the cohort, rename it to the original publication name
Beltran.clinical.df <- read.table(file = "./data/nepc_wcm_2016/data_clinical_sample.txt",sep = "\t",header=T)[c('SAMPLE_ID','DISEASE_CODE')]
Beltran.clinical.df$binary <- ifelse(Beltran.clinical.df$DISEASE_CODE%in% 'CRPC-Adeno','Adeno','t-SCNC')
Beltran.clinical.df <- subset(Beltran.clinical.df,subset=Beltran.clinical.df$SAMPLE_ID%in%colnames(Beltran.RNA.df))
Beltran.clinical.df <- Beltran.clinical.df[order(match(Beltran.clinical.df$SAMPLE_ID,colnames(Beltran.RNA.df))),]



ht_opt$TITLE_PADDING = unit(c(15,15), "points")
set.seed(300)
group = kmeans(t(Beltran.RNA.df),nstart = 1,centers = 5,iter.max = 1000)$cluster

# tested with the heatmap

Beltran_clusters <- list('ARL/NE-' = names(group[group==1]),
                         'AR+/NE+' = names(group[group==2]),
                         'AR-/NE+' = names(group[group==3]),
                         'AR-/NE-' = names(group[group==4]),
                         'AR+/NE-' = names(group[group==5]))

group <- factor(group,labels = c('ARL/NE-','AR+/NE+','AR-/NE+','AR-/NE-','AR+/NE-'))
group <- factor(group,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE-','AR+/NE+'))

table(group)
# AR-/NE- AR-/NE+ ARL/NE- AR+/NE- AR+/NE+
#     8      10      12      15       4


Beltran.classes <- data.frame(samples=unlist(Beltran_clusters),Subtypes=gsub("^[0-9]|[0-9]|[0-9]$", "", names(unlist(Beltran_clusters))))
Beltran.rank.df <- na.omit(Beltran.df)
rownames(Beltran.rank.df) <- Beltran.rank.df$Hugo_Symbol
Beltran.rank.df <- Beltran.rank.df[,-c(1:2)]
Beltran.gex_rank <- rankGenes(as.matrix(Beltran.rank.df))

Beltran.AR.score <- simpleScore(rankData = Beltran.gex_rank, upSet =rownames(Beltran.gex_rank)[rownames(Beltran.gex_rank)%in%AR.panel[-8]],downSet = rownames(Beltran.gex_rank)[rownames(Beltran.gex_rank)%in%AR.panel[8]],knownDirection = TRUE, centerScore = TRUE, subSamples = NULL)[,"TotalScore",drop=F]
Beltran.NEPC.score <- simpleScore(rankData = Beltran.gex_rank, upSet =rownames(Beltran.gex_rank)[rownames(Beltran.gex_rank)%in%NE.panel],knownDirection = TRUE, centerScore = TRUE, subSamples = NULL)[,"TotalScore",drop=F]

#all(rownames(Beltran.AR.score)==rownames(Beltran.NEPC.score))
# TRUE

Beltran.signatures.df <- merge(Beltran.AR.score,Beltran.NEPC.score,by='row.names');Beltran.signatures.df <- Beltran.signatures.df[ order(match(Beltran.signatures.df$Row.names, colnames(Beltran.RNA.df))), ]
colnames(Beltran.signatures.df) <- c('SAMPLE_ID','AR.score','NEPC.score')
#all(Beltran.signatures.df$samples==colnames(Beltran.RNA.df))
# TRUE

Beltran.signatures.df <- merge(Beltran.signatures.df,Beltran.clinical.df[,-2],by='SAMPLE_ID',sort=F)

colAnnot.Beltran = columnAnnotation(
    'AR signature' = anno_simple(Beltran.signatures.df$AR.score, col = colorRamp2(c(-1, 1), c("white", "black")), na_col = "white", gp = gpar(col = "black")),
    'NEPC score' = anno_simple(Beltran.signatures.df$NEPC.score, col = colorRamp2(c(-0.35, 0.35), c("white", "black")), na_col = "white", gp = gpar(col = "black")),
    'Binary class'=Beltran.signatures.df$binary,
    col=my_colour,annotation_height=unit(c(rep(0.3,3)), "cm"),gp=gpar(col="black",lwd=0.7),show_legend=F,
    annotation_name_gp= gpar(fontsize = 8))

rowData.Beltran <- data.frame(GO=NEPC.gene.list$Nelson[NEPC.gene.list$Nelson$original_names %in% rownames(Beltran.RNA.df),]$class)

ht_opt$TITLE_PADDING = unit(c(5,5), "points")

set.seed(100)
Figure1B_heatmap <- draw(Heatmap(matrix = Beltran.RNA.df,name = 'zscore of log2(TPM)',border = T,cluster_rows = F,heatmap_legend_param = list(direction='horizontal'),
                                column_split = factor(group,levels=levels(group)),
                                column_dend_reorder = T,
                                col = col_fun,column_dend_height = unit(10,'mm'),column_dend_side = 'bottom',
                                rect_gp = gpar(col = "black", lwd = 0.4),show_heatmap_legend = F,show_parent_dend_line = F,
                                top_annotation = colAnnot.Beltran,show_column_dend = F,
                                show_column_names = F,cluster_column_slices = F,
                                row_split  = factor(rowData.Beltran$GO,levels=c('NE1','NE2','AR','SQUAM')),row_names_gp = gpar(fontsize=5),cluster_row_slices = F,row_dend_reorder = F,
                                column_names_gp = gpar(fontsize=50),clustering_distance_columns ="pearson",row_title_gp = gpar(fontsize=6,fontface='bold'),
                                column_title_gp = gpar(fill = c("#7570B3","#E7298A","#D95F02","#1B9E77","#66A61E"),fontsize=7,fontface='bold',border='black',col='white')),
                        heatmap_legend_side='bottom',annotation_legend_side='bottom')

ht_opt(RESET = TRUE)
######

### Figure 1C
######


Hallmarks_db <- msigdbr(species = "Homo sapiens",category = c('H'))
Oncogenic_db <- msigdbr(species = "Homo sapiens",category = c('C6'))
WIKI_db <- msigdbr(species = "Homo sapiens",category = c('C2'),subcategory = "CP:WIKIPATHWAYS")

Hallmarks_list = split(x = Hallmarks_db$gene_symbol, f = Hallmarks_db$gs_name)
WIKI_list = split(x = WIKI_db$gene_symbol, f = WIKI_db$gs_name)
genes.list.db <- c(Hallmarks_list,WIKI_list)

genes.db <- rbind(Hallmarks_db,WIKI_db)
pid2gene <- genes.db %>% dplyr::select(gs_id,human_entrez_gene) #TERM2GENE
pid2name <- genes.db %>% dplyr::select(gs_id,gs_name) #TERM2NAME
Count.WCDT <- fread(file = "./data/counts_and_TPMs/Counts_Genes.txt",sep="\t",header = T,data.table = F)
Count.PROMOTE <- fread(file = "./data/counts_and_TPMs/Counts_Genes_PROMOTE.txt",sep="\t",header = T,data.table = F)

RNA.counts.df <-  right_join(Count.WCDT, Count.PROMOTE)
RNA.counts.df <- RNA.counts.df[!duplicated(RNA.counts.df$gene_name),]
rownames(RNA.counts.df) <- RNA.counts.df$gene_name;RNA.counts.df<- RNA.counts.df[,-c(1,which(colnames(RNA.counts.df)%in%c('seqid','gene_type','gene_name')))]

RNA.counts.df.subset <- RNA.counts.df[,colnames(RNA.counts.df)%in%as.character(unlist(Nelsons_clusters))]

##### takes time to re-run, already saved as RData file "GSVA.df.210.RData", loaded
# gsva.df.subset <- gsva(expr = as.matrix(RNA.counts.df.subset),gset.idx.list=genes.list.db,mx.diff=F,
#                 min.sz=5,tau=0.25,kcdf="Poisson",method="ssgsea")
#
# save(gsva.df.subset, file = "./data/GSVA.df.210.RData")
load(file = "./data/GSVA.df.210.RData")


ARPC.enrich.mean <- apply(gsva.df.subset[,match(Nelsons_clusters$`AR+/NE-`, colnames(gsva.df.subset))], 1, mean, na.rm=T)
DNPC.enrich.mean <- apply(gsva.df.subset[,match(Nelsons_clusters$`AR-/NE-`, colnames(gsva.df.subset))], 1, mean, na.rm=T)
SCNPC.enrich.mean <- apply(gsva.df.subset[,match(Nelsons_clusters$`AR-/NE+`, colnames(gsva.df.subset))], 1, mean, na.rm=T)
Amp.enrich.mean <- apply(gsva.df.subset[,match(Nelsons_clusters$`AR+/NE+`, colnames(gsva.df.subset))], 1, mean, na.rm=T)
ARLPC.enrich.mean <- apply(gsva.df.subset[,match(Nelsons_clusters$`ARL/NE-`, colnames(gsva.df.subset))], 1, mean, na.rm=T)

### rank the pathways based on variability
enrichment.df <- data.frame(pathways=unique(c(names(ARPC.enrich.mean[abs(ARPC.enrich.mean)>0.3]),
                                              names(DNPC.enrich.mean[abs(DNPC.enrich.mean)>0.3]),
                                              names(SCNPC.enrich.mean[abs(SCNPC.enrich.mean)>0.3]),
                                              names(ARLPC.enrich.mean[abs(ARLPC.enrich.mean)>0.3]),
                                              names(Amp.enrich.mean[abs(Amp.enrich.mean)>0.3]))))
enrichment.df <- merge(enrichment.df,data.frame(DNPC=DNPC.enrich.mean),by.x='pathways',by.y='row.names',all.x=T)
enrichment.df <- merge(enrichment.df,data.frame(SCNPC=SCNPC.enrich.mean),by.x='pathways',by.y='row.names',all.x=T)
enrichment.df <- merge(enrichment.df,data.frame(ARLPC=ARLPC.enrich.mean),by.x='pathways',by.y='row.names',all.x=T)
enrichment.df <- merge(enrichment.df,data.frame(ARPC=ARPC.enrich.mean),by.x='pathways',by.y='row.names',all.x=T)
enrichment.df <- merge(enrichment.df,data.frame(Amp=Amp.enrich.mean),by.x='pathways',by.y='row.names',all.x=T)
enrichment.df$variance <- rowVars(as.matrix(enrichment.df[,-1]))
enrichment.df <- enrichment.df[enrichment.df$pathways%like%'HALLMARK',]

enrichment.df <- enrichment.df[order(enrichment.df$variance,decreasing = T),][1:30,]
enrichment.df$pathways <- gsub("WP_","",enrichment.df$pathways)
enrichment.df$pathways <- gsub("_"," ",enrichment.df$pathways)
enrichment.df$pathways <- gsub("HALLMARK ","",enrichment.df$pathways)
enrichment.df$pathways <- gsub("[.]"," ",enrichment.df$pathways)
enrichment.df$pathways <- str_to_title(enrichment.df$pathways)

rownames(enrichment.df)<-enrichment.df$pathways
row.order <- enrichment.df[order(enrichment.df$variance,decreasing = T),]$pathways
enrichment.df <- enrichment.df[,-c(1,7)]
colnames(enrichment.df) <- c('AR-/NE-','AR-/NE+','ARL/NE-',
                             'AR+/NE-','AR+/NE+')
enrichment.df <- t(scale(t(enrichment.df)))



set.seed(100)
col_fun =colorRamp2(c(-2,0,2), c("royalblue","white","darkred"))

Figure1C_heatmap <- draw(Heatmap(matrix = enrichment.df,name = 'Mean enrichment (z-score)',border = T,
                               heatmap_legend_param = list(direction='horizontal'),
                               cluster_columns = F,
                               col = col_fun,column_split = factor(c('AR-','AR-','AR+','AR+','AR+'),levels=c('AR-','AR+')),
                               cluster_rows = T,na_col = 'white',row_dend_gp = gpar(col = "black", lwd = 1.5),
                               rect_gp = gpar(col = "black", lwd = 1),column_gap = unit(5,'mm'),
                               show_heatmap_legend = F,row_order = row.order,
                               show_column_names = T,row_title_gp = gpar(fontsize=4,fontface='bold'),
                               width = ncol(enrichment.df)*unit(5, "mm"),
                               height = nrow(enrichment.df)*unit(3, "mm"),
                               column_names_gp = gpar(fontsize=8,fontface='bold'),
                               column_title_gp = gpar(fontsize=10,fontface='bold'),

                               row_names_gp = gpar(fontsize=10)),
                       heatmap_legend_side='bottom',annotation_legend_side='bottom')
######



