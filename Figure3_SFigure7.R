library(Gviz)
library(data.table)
library(ggVennDiagram)
library(ggpubr)
library(ggplot2)
library(biomaRt)
library(cowplot)
library(ggthemes)
library(dplyr)
library(AnnotationDbi)
library(GenomicRanges)
library(ComplexHeatmap)
library("org.Hs.eg.db")
library(rtracklayer)

setwd("~/Box Sync/UCSF/Manuscript/Submissions/GitHub/")
source("./functions.R")

gene.loc <- fread("./data/GOI.hg38.location.bed",header = T,data.table = F) ### genebody +/- 10kb
load(file="./data/gene.H3K27ac.DNASE1.tracks.GR.RData")
load(file="./data/WCDT.5groups.210.RData")
load(file="./data/Matched.WGS.WGBS.samples.RData")
load("./data/DMR.data.RData")
load(file="./data/Matched.WGS.WGBS.samples.RData")
load("./data/CHD7.mean.methylation.RData")
load(file = "./data/NEPC.gene.list.RData")

ASCL1.loc <- fread("./data/ASCL1.bed",header = F,data.table = F) ### genebody +/- 10kb
ASCL1.loc$Cell <- c('SCLC: NCI-H889  ','SCLC: NCI-H2107  ','SCLC: NCI-H2107  ','SCLC: NCI-H128  ',
                    'SCLC: NCI-H128  ','SCLC: NCI-H2107  ','Neuroblastoma: Kelly  ','SCLC: NCI-H1184  ','NE-NSCLC: NCI-1755  ')
ASCL1.GR <- makeGRangesFromDataFrame(ASCL1.loc[,-c(2:4,6,9)],seqnames.field = "V1",start.field = "V7",end.field = "V8",keep.extra.columns = T)


#### CHD7 gene expression
TPM.WCDT <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes.txt",sep="\t",header = T,data.table = F)
TPM.PROMOTE <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt",sep="\t",header = T,data.table = F)

RNA.TPM.df <-  right_join(TPM.WCDT, TPM.PROMOTE)
RNA.TPM.df <- RNA.TPM.df[!duplicated(RNA.TPM.df$gene_name),]
rownames(RNA.TPM.df) <- RNA.TPM.df$gene_name;RNA.TPM.df<- RNA.TPM.df[,-c(1,which(colnames(RNA.TPM.df)=='gene_name'))]

genes.annotation <- fread(file = "./data/DriverGenePanel.38.tsv")
genelist1 <- read.table("./data/GENE.LIST.CELL.txt",header = T,sep = "\t")
genelist2 <- openxlsx::read.xlsx("./data/prostate_cancer.genes.xlsx",sheet = 2)
genes.of.interest <- unique(c(genelist2$Gene,genelist1$GENES))
genes.of.interest <- left_join(data.frame(gene=genes.of.interest),genes.annotation[,c('gene','likelihoodType')])
genes.of.interest <- na.omit(genes.of.interest)
genes.of.interest$id <- mapIds(org.Hs.eg.db, genes.of.interest$gene, 'ENTREZID', 'SYMBOL')

expression.df <- as.data.frame(t(log2(RNA.TPM.df+1)[c('AR','PTEN','RB1','BRCA2','CHD7'),]))
expression.df$samples <- rownames(expression.df)
expression.df <- left_join(Nelson.samples.PCA[,-3],expression.df)
expression.df$Subtypes <- factor(expression.df$Subtypes,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE-','AR+/NE+'))

my_comparisons.exprs <- list(c("AR-/NE+", "AR-/NE-"),c("AR-/NE+", "ARL/NE-"),
                             c("AR-/NE+", "AR+/NE-"),c("AR-/NE+", "AR+/NE+"))


#### Figure 3A
######
Figure3A_boxplot <- ggboxplot(melt(expression.df[c('samples','Subtypes','CHD7')]),x='Subtypes',y='value',ylab = 'log2 (TPM+1)',xlab='Subtypes',fill='Subtypes',
                          bxp.errorbar = T,bxp.errorbar.width = 0.05)+
    facet_wrap(~variable,ncol = 3)+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+
    stat_compare_means(size=4,label.x=3.5,fontface='bold',
                       label.y=7)+
    stat_compare_means(size=9,tip.length = 0.01,
                       label.y=7,
                       label.y.npc = 0.5,step.increase = 0.145,comparisons = my_comparisons.exprs,hide.ns = F,label = 'p.signif',p.adjust.methods='fdr')+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=20),axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=20),axis.ticks.x  = element_blank(),
          title = element_text(face = 'bold',size=22),
          legend.position = 'none',strip.text = element_text(face='bold',size=20),legend.key.size = unit(10,'mm'),legend.text = element_text(size=15),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())+coord_cartesian(ylim = c(0,10.5))
######


bm <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")
Nelson.subtype.methyl <- Nelson.samples.PCA[Nelson.samples.PCA$samples%in% colnames(CHD7),]

DSS.df <- na.omit(DMRs.data)
DSS.df <- makeGRangesFromDataFrame(DSS.df,keep.extra.columns = T,ignore.strand = T,seqnames.field = "chr",start.field = "start",end.field = "end")
DSS.df <- sort(DSS.df)

centromere.coords = read.delim(file="./data/grch38_centromeres.txt", sep='\t', header=TRUE,stringsAsFactors=FALSE)
colnames(centromere.coords) <- c('chrom','start','end')
telomere.coords = read.delim(file="./data/grch_telomeres.txt", sep='\t', header=TRUE,stringsAsFactors=FALSE)
telomere.coords.5prime = telomere.coords[seq(from=1,to=nrow(telomere.coords), by=2),]
telomere.coords.3prime = telomere.coords[seq(from=2,to=nrow(telomere.coords), by=2),]

blacklist.regions = read.delim(file="./data/grch38_gap.txt", sep='\t', header=TRUE,stringsAsFactors=FALSE)
blacklist.regions <- blacklist.regions[,1:3];colnames(blacklist.regions) <- c('chrom','start','end')
blacklist.regions <- rbind.data.frame(blacklist.regions, centromere.coords,telomere.coords)
blacklist.regions.gr <- makeGRangesFromDataFrame(blacklist.regions)
DSS.df <- DSS.df[-queryHits(findOverlaps(DSS.df, blacklist.regions.gr, type='any')),]


tmp.DMR <- unique(DSS.df[seqnames(DSS.df)=='chr8'&start(DSS.df)>=gene.loc[gene.loc$gene=='CHD7',]$start+10000&
                             end(DSS.df)<=gene.loc[gene.loc$gene=='CHD7',]$end-10000
                         &DSS.df$subtype%in%c('SCNPC')&abs(DSS.df$diff.Methy)>0.1])
highlights.H3k27 <- subsetByOverlaps(tmp.DMR,H3K27ac.GR,minoverlap = 50)

# Each DHS is annotated with a signal level which can be used as a confidence score,
# as well as estimates for its summit position and summit variability (as captured by the 'core' coordinates). #mean_signal
#https://www.encodeproject.org/annotations/ENCSR857UZV/

DNASE1.CHD7 <- DNASE1.GR[seqnames(DNASE1.GR)=='chr8'&start(DNASE1.GR)>=gene.loc[gene.loc$gene=='CHD7',]$start&end(DNASE1.GR)<=gene.loc[gene.loc$gene=='CHD7',]$end]
H3k27ac.CHD7 <- H3K27ac.GR[seqnames(H3K27ac.GR)=='chr8'&start(H3K27ac.GR)>=gene.loc[gene.loc$gene=='CHD7',]$start&end(H3K27ac.GR)<=gene.loc[gene.loc$gene=='CHD7',]$end]
promoter.CHD7 <- makeGRangesFromDataFrame(data.frame(chr="chr8", end=60678810,start= 60678462,strand = "*")) # EH38E3836203 from ENCODE
#### PLOTTING

idxTrack <- IdeogramTrack(genome="hg38", chromosome="chr8",fontsize=20,fontface='bold',fontcolor='black')
axTrack <- GenomeAxisTrack(fontsize=12,fontface='bold',fontcolor='black',col='black',labelPos='above',add53 = TRUE,add35=T)
biomTrack <- BiomartGeneRegionTrack(genome = "hg38", name = "refseq",filter = list(with_refseq_mrna = TRUE),col.line = 'black', col = 'royalblue',fill='royalblue',collapseTranscripts='meta',
                                    symbol = "CHD7", biomart = bm,color='black',background.title='white',fontcolor.title='black',fontsize=14,protein_coding='royalblue',pseudogene='black',utr5='royalblue',utr3='royalblue')
H3K27ac.Track <- (AnnotationTrack(H3k27ac.CHD7,
                                  genome = 'hg38',name='H3K27ac',fill='black',col='black',background.title='white',fontcolor.title='black',fontsize=14,stacking='dense'))
promoter.Track <- (AnnotationTrack(promoter.CHD7,
                                   genome = 'hg38',fill='gold',col='gold',background.title='white',fontcolor.title='black',fontsize=14,stacking='dense'))

DNASE1.Track <- (DataTrack(DNASE1.CHD7[,'mean_signal'],col.histogram='black',
                           genome = 'hg38',name='DHSs',fill.histogram = "black",background.title='white',col.axis='black',fontcolor.title='black',type=c('histogram'),fontsize=14))

DMR.Track <- (AnnotationTrack(highlights.H3k27,group = c('1','2','3','4'),feature = c('1','2','3','4'),
                              genome = 'hg38',name='DMR',fill=subtype.colors[subtype.colors$subtypes=='AR-/NE+',]$colors,
                              col=subtype.colors[subtype.colors$subtypes=='AR-/NE+',]$colors,background.title='white',
                              fontcolor.title='black',fontsize=14,stacking='dense'))

ASCL1.groups <- AnnotationTrack(start = c(ASCL1.loc$V2),
                                width = c(ASCL1.loc$V3-ASCL1.loc$V2),
                                chromosome = "chr8",fill='darkred',col='darkred',
                                strand = rep(c("*"),nrow(ASCL1.loc)),collapse=F,
                                group = ASCL1.loc$Cell,feature = ASCL1.loc$Cell,background.title='white',
                                fontcolor.title='black',fontsize=14,
                                genome = "hg38", name = "ASCL1")

promoter.ht <- OverlayTrack(trackList = list(H3K27ac.Track,promoter.Track),background.title='white',fontcolor.title='black',fontsize=14)

ht <- HighlightTrack(trackList = list(biomTrack,
                                      promoter.ht,
                                      DMR.Track,
                                      DNASE1.Track,
                                      ASCL1.groups),
start = c(
    gene.loc[gene.loc$gene=='CHD7',]$start+10000,
    gene.loc[gene.loc$gene=='CHD7',]$end-10000),
width = c(
    20,20),col=c('forestgreen','brown'),fill=c('forestgreen','brown'),alpha=c(1),lty=c(3,3),
chromosome = 8,inBackground=T)



##### Figure 3B
########

plotTracks(list(idxTrack,
                axTrack,
                ht),groupAnnotation='group',
           foo = "darkred", bar = "darkgreen",just.group="left",collapse=F,
           fontcolor.group='black',fontsize.group=12,fontface.group='bold',fontface.main ='bold',sizes = c(0.04,0.1,0.05,0.05,0.05,0.05,0.1))

########

#### matched WGS, WGBS with RNA-seq data

Nelson.samples.PCA <- Nelson.samples.PCA[Nelson.samples.PCA$samples%in%Matched.WGS.WGBS.samples$samples,]
Nelson.samples.PCA$Subtypes <- factor(Nelson.samples.PCA$Subtypes,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE-','AR+/NE+'))
Nelson.samples.PCA <- Nelson.samples.PCA[order(Nelson.samples.PCA$Subtypes),]


my_comparisons <- list(c("AR-/NE-", "AR-/NE+"),c("AR-/NE-", "AR+/NE-"),c("AR-/NE+", "AR+/NE-"))
methylation.df <- data.frame(DMR1=colMeans(CHD7[CHD7$pos>=start(highlights.H3k27)[1]&CHD7$pos<=end(highlights.H3k27)[1],-c(1:2)],na.rm = T),
                             DMR2=colMeans(CHD7[CHD7$pos>=start(highlights.H3k27)[2]&CHD7$pos<=end(highlights.H3k27)[2],-c(1:2)],na.rm = T),
                             DMR3=colMeans(CHD7[CHD7$pos>=start(highlights.H3k27)[3]&CHD7$pos<=end(highlights.H3k27)[3],-c(1:2)],na.rm = T),
                             DMR4=colMeans(CHD7[CHD7$pos>=start(highlights.H3k27)[4]&CHD7$pos<=end(highlights.H3k27)[4],-c(1:2)],na.rm = T))


methylation.df$samples <- rownames(methylation.df)
methylation.df <- left_join(Nelson.samples.PCA[,-3],methylation.df)
methylation.df <- methylation.df[methylation.df$Subtypes%out%c('AR+/NE+','ARL/NE-'),]
methylation.df$Subtypes <- factor(methylation.df$Subtypes,levels=c('AR-/NE-','AR-/NE+','AR+/NE-'))



###### Adding correlation between DMR 1-4 methylation and CHD7 gene expression in AR- groups and AR+/NE-

expression.df <- expression.df[expression.df$Subtypes%out%c('AR+/NE+','ARL/NE-'),]
expression.df$Subtypes <- factor(expression.df$Subtypes,levels=c('AR-/NE-','AR-/NE+','AR+/NE-'))

methylation.df <- merge(methylation.df,expression.df[,-2],by='samples',all.x=T)


##### methylation level - boxplots (Figure 3C-F) and Venn-diagram Figure 3G
######
Figure3C_boxplot <- ggboxplot(methylation.df,x='Subtypes',y='DMR1',ylab = 'Mean methylation (%)',xlab='',fill = 'Subtypes', title = paste0('DMR1 (',chrom(highlights.H3k27)[1],": ",start(highlights.H3k27)[1]," - ",end(highlights.H3k27)[1],")"),
                         bxp.errorbar = T,bxp.errorbar.width = 0.05,add = 'jitter',add.params = list(size=1.2))+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+
    scale_y_continuous(labels = scales::percent,expand =c(0,0.1),breaks = c(0.25,0.5,0.75,1))+

    stat_compare_means(size=6,label.x=1,fontface='bold',label.y=0.27)+
    stat_compare_means(size=12,tip.length = 0.001,label.y=0.7,label.y.npc = 0.5,step.increase = 0.02,comparisons = my_comparisons,hide.ns = F,label = 'p.signif',p.adjust.methods='fdr')+
    stat_cor(data =methylation.df,mapping = aes(x=methylation.df$DMR1,y = methylation.df$CHD7),
             method = "pearson",p.accuracy = 0.001, r.accuracy = 0.01,label.x=0.71,label.y=0.21,size=5)+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=20),axis.title.x = element_text(size=20),
          axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5,size=20),
          axis.text.y = element_text(size=20),title = element_text(face = 'bold',size=13),
          legend.position = 'none',strip.text = element_text(face='bold',size=12.5),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())

Figure3D_boxplot <- ggboxplot(methylation.df,x='Subtypes',y='DMR2',ylab = 'Mean methylation (%)',xlab='',fill = 'Subtypes', title = paste0('DMR2 (',chrom(highlights.H3k27)[2],": ",start(highlights.H3k27)[2]," - ",end(highlights.H3k27)[2],")"),
                         bxp.errorbar = T,bxp.errorbar.width = 0.05,add = 'jitter',add.params = list(size=1.2))+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+
    scale_y_continuous(labels = scales::percent,expand =c(0,0.1),breaks = c(0.25,0.5,0.75,1))+

    stat_compare_means(size=6,label.x=1,fontface='bold',label.y=0.27)+
    stat_compare_means(size=12,tip.length = 0.001,label.y=0.7,label.y.npc = 0.5,step.increase = 0.02,comparisons = my_comparisons,hide.ns = F,label = 'p.signif',p.adjust.methods='fdr')+
    stat_cor(data =methylation.df,mapping = aes(x=methylation.df$DMR2,y = methylation.df$CHD7),
             method = "pearson",p.accuracy = 0.001, r.accuracy = 0.01,label.x=0.7,label.y=0.20,size=5)+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=20),axis.title.x = element_text(size=20),
          axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5,size=20),
          axis.text.y = element_text(size=20),title = element_text(face = 'bold',size=13),
          legend.position = 'none',strip.text = element_text(face='bold',size=12.5),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())


Figure3E_boxplot <- ggboxplot(methylation.df,x='Subtypes',y='DMR3',ylab = 'Mean methylation (%)',xlab='',fill = 'Subtypes', title = paste0('DMR3 (',chrom(highlights.H3k27)[3],": ",start(highlights.H3k27)[3]," - ",end(highlights.H3k27)[3],")"),
                         bxp.errorbar = T,bxp.errorbar.width = 0.05,add = 'jitter',add.params = list(size=1.2))+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+
    scale_y_continuous(labels = scales::percent,expand =c(0,0.12),breaks = c(0.25,0.5,0.75,1))+

    stat_compare_means(size=6,label.x=1,fontface='bold',label.y=0.16)+
    stat_compare_means(size=12,tip.length = 0.001,label.y=0.7,label.y.npc = 0.5,step.increase = 0.02,comparisons = my_comparisons,hide.ns = F,label = 'p.signif',p.adjust.methods='fdr')+
    stat_cor(data =methylation.df,mapping = aes(x=methylation.df$DMR3,y = methylation.df$CHD7),
             method = "pearson",p.accuracy = 0.001, r.accuracy = 0.01,label.x=0.7,label.y=0.09,size=5)+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=20),axis.title.x = element_text(size=20),
          axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5,size=20),
          axis.text.y = element_text(size=20),title = element_text(face = 'bold',size=13),
          legend.position = 'none',strip.text = element_text(face='bold',size=12.5),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())


Figure3F_boxplot <- ggboxplot(methylation.df,x='Subtypes',y='DMR4',ylab = 'Mean methylation (%)',xlab='',fill = 'Subtypes', title = paste0('DMR4 (',chrom(highlights.H3k27)[4],": ",start(highlights.H3k27)[4]," - ",end(highlights.H3k27)[4],")"),
                         bxp.errorbar = T,bxp.errorbar.width = 0.05,add = 'jitter',add.params = list(size=1.2))+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+
    scale_y_continuous(labels = scales::percent,expand =c(0,0.1),breaks = c(0.25,0.5,0.75,1))+

    stat_compare_means(size=6,label.x=1,fontface='bold',label.y=0.33)+
    stat_compare_means(size=12,tip.length = 0.001,label.y=0.7,label.y.npc = 0.5,step.increase = 0.02,comparisons = my_comparisons,hide.ns = F,label = 'p.signif',p.adjust.methods='fdr')+
    stat_cor(data =methylation.df,mapping = aes(x=methylation.df$DMR4,y = methylation.df$CHD7),
             method = "pearson",p.accuracy = 0.001, r.accuracy = 0.01,label.x=0.68,label.y=0.27,size=5)+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=20),axis.title.x = element_text(size=20),
          axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5,size=20),
          axis.text.y = element_text(size=20),title = element_text(face = 'bold',size=13),
          legend.position = 'none',strip.text = element_text(face='bold',size=12.5),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())

Figure3C_F <- plot_grid(Figure3C_boxplot+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank()),Figure3D_boxplot+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank()),
          Figure3E_boxplot,Figure3F_boxplot+theme(axis.title.y = element_blank()),ncol = 2,align = 'hv')


##### checking the DMR regions - what motifs are connecting to them

DNPC.genes <- openxlsx::read.xlsx(xlsxFile = "./data/SData1.xlsx",sheet = "AR- NE-")
DNPC.genes <- DNPC.genes[abs(DNPC.genes$log2FoldChange)>=1&DNPC.genes$padj<0.01,]
SCNPC.genes <- openxlsx::read.xlsx(xlsxFile = "./data/SData1.xlsx",sheet = "AR- NE+")
SCNPC.genes <- SCNPC.genes[abs(SCNPC.genes$log2FoldChange)>=1&SCNPC.genes$padj<0.01,]
ARLPC.genes <- openxlsx::read.xlsx(xlsxFile = "./data/SData1.xlsx",sheet = "ARL NE-")
ARLPC.genes <- ARLPC.genes[abs(ARLPC.genes$log2FoldChange)>=1&ARLPC.genes$padj<0.01,]
Amp.genes <- openxlsx::read.xlsx(xlsxFile = "./data/SData1.xlsx",sheet = "AR+ NE+")
Amp.genes <- Amp.genes[abs(Amp.genes$log2FoldChange)>=1&Amp.genes$padj<0.01,]

DMR1 <- na.omit(fread("./data/FIMO/DMR_1/fimo.tsv")[,c('motif_alt_id','sequence_name','q-value')]);DMR1$ID <- toupper(gsub(".*[.]","",DMR1$motif_alt_id));DMR1 <- DMR1[!duplicated(DMR1$ID),]
DMR1 <- DMR1[DMR1$`q-value`<0.05,];DMR1 <- DMR1[order(DMR1$`q-value`,decreasing = F),];DMR1 <- DMR1[1:10,];DMR1 <-na.omit(DMR1)
DMR2 <- na.omit(fread("./data/FIMO/DMR_2/fimo.tsv")[,c('motif_alt_id','sequence_name','q-value')]);DMR2$ID <- toupper(gsub(".*[.]","",DMR2$motif_alt_id));DMR2 <- DMR2[!duplicated(DMR2$ID),]
DMR2 <- DMR2[DMR2$`q-value`<0.05,];DMR2 <- DMR2[order(DMR2$`q-value`,decreasing = F),];DMR2 <- DMR2[1:10,];DMR2 <-na.omit(DMR2)
DMR3 <- na.omit(fread("./data/FIMO/DMR_3/fimo.tsv")[,c('motif_alt_id','sequence_name','q-value')]);DMR3$ID <- toupper(gsub(".*[.]","",DMR3$motif_alt_id));DMR3 <- DMR3[!duplicated(DMR3$ID),]
DMR3 <- DMR3[DMR3$`q-value`<0.05,];DMR3 <- DMR3[order(DMR3$`q-value`,decreasing = F),];DMR3 <- DMR3[1:10,];DMR3 <-na.omit(DMR3)
DMR4 <- na.omit(fread("./data/FIMO/DMR_4/fimo.tsv")[,c('motif_alt_id','sequence_name','q-value')]);DMR4$ID <- toupper(gsub(".*[.]","",DMR4$motif_alt_id));DMR4 <- DMR4[!duplicated(DMR4$ID),]
DMR4 <- DMR4[DMR4$`q-value`<0.05,];DMR4 <- DMR4[order(DMR4$`q-value`,decreasing = F),];DMR4 <- DMR4[1:10,];DMR4 <-na.omit(DMR4)

venn.motif <- list(DMR1=DMR1$ID,
                   DMR2=DMR2$ID,
                   DMR3=DMR3$ID,
                   DMR4=DMR4$ID)


venn.motif.data <- process_data(Venn(venn.motif))


Figure3G_venndiagram <- ggplot() +
    geom_sf(data = venn_region(venn.motif.data),fill='white') + ggtitle('')+
    geom_sf(aes(color = name), data = venn_setedge(venn.motif.data), show.legend = TRUE, size = 2) +
    geom_sf_text(aes(label = name,fontface='bold'),nudge_x = c(0.05,0.05,-0.05,-0.05),size=7,data = venn_setlabel(venn.motif.data)) +
    geom_sf_text(aes(label = count),
                 data = venn_region(venn.motif.data),
                 size = 10) +
    scale_color_colorblind()+
    scale_fill_colorblind()+
    theme_void()+
    theme(legend.position = 'none',title = element_text(face = 'bold'))+
    annotate("segment", x = 0.2, y = 0.83, xend = 0.3, yend = 0.77,
             arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
    annotate("text", x=0.14,y=0.82, label = "ASCL1", color = "black",
             size = 5, fontface = "bold") +
    annotate("text", x=0.14,y=0.85, label = "BCL11B", color = "black",
             size = 5, fontface = "bold") +
    annotate("segment", x = 0.8, y = 0.8, xend = 0.7, yend = 0.77,
             arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
    annotate("text", x=0.86,y=0.82, label = "OLIG2", color = "black",
             size = 5, fontface = "bold")+
    annotate("text", x=0.86,y=0.85, label = "NEUROG2", color = "black",
             size = 5, fontface = "bold")

######



######### data from Beltran et al. - Suppelementary Figure 7

Beltran.df <-read.table(file = "./data/nepc_wcm_2016/data_RNA_Seq_mRNA_median_all_sample_Zscores.txt",sep = "\t",header=T)
Beltran.RNA.df <- subset(Beltran.df,subset=Beltran.df$Hugo_Symbol%in%NEPC.gene.list$Nelson$GENE_ID)
rownames(Beltran.RNA.df) <- Beltran.RNA.df$Hugo_Symbol;Beltran.RNA.df <- Beltran.RNA.df[,-c(1,2)]
Beltran.RNA.df <- na.omit(Beltran.RNA.df)
Beltran.RNA.df <- Beltran.RNA.df[order(match(rownames(Beltran.RNA.df),NEPC.gene.list$Nelson$GENE_ID)),]
Beltran.clinical.df <- read.table(file = "./data/nepc_wcm_2016/data_clinical_sample.txt",sep = "\t",header=T)[c('SAMPLE_ID','DISEASE_CODE')]
Beltran.clinical.df$binary <- ifelse(Beltran.clinical.df$DISEASE_CODE%in% 'CRPC-Adeno','Adeno','t-SCNC')
Beltran.clinical.df <- subset(Beltran.clinical.df,subset=Beltran.clinical.df$SAMPLE_ID%in%colnames(Beltran.RNA.df))
Beltran.clinical.df <- Beltran.clinical.df[order(match(Beltran.clinical.df$SAMPLE_ID,colnames(Beltran.RNA.df))),]

#######
ht_opt$TITLE_PADDING = unit(c(15,15), "points")
set.seed(300)
group = kmeans(t(Beltran.RNA.df),nstart = 1,centers = 5,iter.max = 1000)$cluster
# based on the heatmap in figure 1B
Beltran_clusters <- list('ARL/NE-' = names(group[group==1]),
                         'AR+/NE+' = names(group[group==2]),
                         'AR-/NE+' = names(group[group==3]),
                         'AR-/NE-' = names(group[group==4]),
                         'AR+/NE-' = names(group[group==5]))


Beltran.classes <- data.frame(samples=unlist(Beltran_clusters),Subtypes=gsub("^[0-9]|[0-9]|[0-9]$", "", names(unlist(Beltran_clusters))))

beltran.CHD7.df <- t(Beltran.df[Beltran.df$Hugo_Symbol%in%c('CHD7'),-c(1,2)])
colnames(beltran.CHD7.df) <- c('CHD7')

beltran.subtype <- data.frame(Subtypes=names(unlist2(Beltran_clusters)),samples=as.character(unlist2(Beltran_clusters)))
beltran.subtype <- merge(beltran.subtype,beltran.CHD7.df,by.x='samples',by.y='row.names')
beltran.subtype$Subtypes <- factor(beltran.subtype$Subtypes,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE-','AR+/NE+'))


my_comparisons <- list(c("AR-/NE+", "AR-/NE-"), c("AR-/NE+", "ARL/NE-"), c("AR-/NE+", "AR+/NE-"), c("AR-/NE+", "AR+/NE+"))

SFigure7_boxplot <- ggboxplot(beltran.subtype,x='Subtypes',y='CHD7',ylab = 'z-score',xlab='Subtypes',fill='Subtypes',
                                    bxp.errorbar = T,bxp.errorbar.width = 0.05)+
    scale_fill_manual(breaks=subtype.colors$subtypes,
                      values=subtype.colors$colors,name='Subtypes')+
    stat_compare_means(size=5,label.x=3.5,fontface='bold',label.y=2.5)+
    stat_compare_means(size=8,tip.length = 0.015,label.y=2.5,label.y.npc = 0.5,step.increase = 0.1,comparisons = my_comparisons,hide.ns = F,label = 'p.signif',p.adjust.methods='fdr',method = 't.test')+
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),
          axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5,size=15),legend.title = element_text(size=15,face = 'bold'),
          axis.text.y = element_text(size=15),legend.key.size = unit(10,units = 'mm'),legend.text = element_text(size=13),
          legend.position = 'none',strip.text = element_text(face='bold',size=12.5),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())




