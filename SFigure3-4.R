### SFigure 3 and SFigure 4 - DE genes and GSEA BP
library(data.table)
library(openxlsx)
library(DESeq2)
library(VennDiagram)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(foreach)
library(ggpubr)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(cowplot)
library(tibble)
setwd("~/Box Sync/UCSF/Manuscript/Submissions/GitHub")
source("./functions.R")

load(file="./data/2022_11_04_WCDT_poppy.Rdata")
load(file="./data/WCDT.5groups.210.RData")
load(file="./data/Matched.WGS.WGBS.samples.RData")
load(file = "./data/NEPC.gene.list.RData")

Count.WCDT <- fread(file = "./data/counts_and_TPMs/Counts_Genes.txt",sep="\t",header = T,data.table = F)
Count.PROMOTE <- fread(file = "./data/counts_and_TPMs/Counts_Genes_PROMOTE.txt",sep="\t",header = T,data.table = F)

RNA.counts.df <-  right_join(Count.WCDT, Count.PROMOTE)
RNA.counts.df <- RNA.counts.df[!duplicated(RNA.counts.df$gene_name),]
rownames(RNA.counts.df) <- RNA.counts.df$gene_name;RNA.counts.df<- RNA.counts.df[,-c(1,which(colnames(RNA.counts.df)%in%c('seqid','gene_type','gene_name')))]

gene.df <- fread("./data/ensembl2sym.txt",data.table = F)
prot.genes <- gene.df %>% dplyr::filter(gene.df$type%in%"protein_coding")
Nelson.samples.PCA <- Nelson.samples.PCA[Nelson.samples.PCA$samples%in%Matched.WGS.WGBS.samples$samples,]

Nelson.samples.PCA$Subtypes.ID <- as.factor(ifelse(Nelson.samples.PCA$Subtypes%in%"ARL/NE-",'ARLPC',
                                                   ifelse(Nelson.samples.PCA$Subtypes%in%"AR+/NE-",'ARPC',
                                                          ifelse(Nelson.samples.PCA$Subtypes%in%"AR-/NE+",'SCNPC',
                                                                 ifelse(Nelson.samples.PCA$Subtypes%in%"AR+/NE+",'Amphicrine',
                                                                        ifelse(Nelson.samples.PCA$Subtypes%in%"AR-/NE-",'DNPC',NA))))))



purity.ploidy.df <- foreach(i=names(SO_WCDT_WGS),.combine = 'rbind') %do% {
    data.frame(SO_WCDT_WGS[[i]][c('purity','ploidy')])}
purity.ploidy.df$samples <- names(SO_WCDT_WGS)
purity.ploidy.df <- purity.ploidy.df[purity.ploidy.df$samples%in%Nelson.samples.PCA$samples,]

Nelson.samples.PCA <- merge(purity.ploidy.df,Nelson.samples.PCA,by='samples')

DE.counts <- RNA.counts.df[,colnames(RNA.counts.df)%in%Nelson.samples.PCA$samples]
DE.counts <- DE.counts[,order(match(colnames(DE.counts),Nelson.samples.PCA$samples))]

# all(colnames(DE.counts)==Nelson.samples.PCA$samples)
#TRUE

sample_annot <- data.frame(condition=Nelson.samples.PCA$Subtypes.ID,
                           purity=Nelson.samples.PCA$purity,
                           ploidy=Nelson.samples.PCA$ploidy,
                           sample=colnames(DE.counts),
                           row.names = colnames(DE.counts))

dds <- DESeqDataSetFromMatrix(countData=DE.counts, colData=sample_annot, design= ~ condition+purity+ploidy)
dds <- dds[rowSums(counts(dds))>0, ]
set.seed(100)
dds <- DESeq(dds)

res.dds.ARPC.ARLPC <- lfcShrink(dds = dds,contrast = c("condition","ARLPC","ARPC"),res=results(dds,contrast = c("condition","ARLPC","ARPC")),type = 'ashr');res.dds.ARPC.ARLPC <- res.dds.ARPC.ARLPC[order(res.dds.ARPC.ARLPC$padj), ]
res.dds.ARPC.DNPC <- lfcShrink(dds = dds,contrast = c("condition","DNPC","ARPC"),res=results(dds,contrast = c("condition","DNPC","ARPC")),type = 'ashr');res.dds.ARPC.DNPC <- res.dds.ARPC.DNPC[order(res.dds.ARPC.DNPC$padj), ]
res.dds.ARPC.SCNPC <- lfcShrink(dds = dds,contrast = c("condition","SCNPC","ARPC"),res=results(dds,contrast = c("condition","SCNPC","ARPC")),type = 'ashr');res.dds.ARPC.SCNPC <- res.dds.ARPC.SCNPC[order(res.dds.ARPC.SCNPC$padj), ]
res.dds.ARPC.Amphicrine <- lfcShrink(dds = dds,contrast = c("condition","Amphicrine","ARPC"),res=results(dds,contrast = c("condition","Amphicrine","ARPC")),type = 'ashr');res.dds.ARPC.Amphicrine <- res.dds.ARPC.Amphicrine[order(res.dds.ARPC.Amphicrine$padj), ]

res.dds.all <- list('DNPC'=res.dds.ARPC.DNPC,'ARLPC'=res.dds.ARPC.ARLPC,'SCNPC'=res.dds.ARPC.SCNPC,'Amphicrine'=res.dds.ARPC.Amphicrine)


########

##### generate significant results with log2FC more than 1 FDR < 0.01

######## UPREGULATED
res.fdr01.FC2.SIG.ARPC.up.ARLPC <- subset(res.dds.ARPC.ARLPC,res.dds.ARPC.ARLPC$padj <= 0.01&(res.dds.ARPC.ARLPC$log2FoldChange)>=1)
res.fdr01.FC2.SIG.ARPC.up.DNPC <- subset(res.dds.ARPC.DNPC,res.dds.ARPC.DNPC$padj <= 0.01&(res.dds.ARPC.DNPC$log2FoldChange)>=1)
res.fdr01.FC2.SIG.ARPC.up.SCNPC <- subset(res.dds.ARPC.SCNPC,res.dds.ARPC.SCNPC$padj <= 0.01&(res.dds.ARPC.SCNPC$log2FoldChange)>=1)
res.fdr01.FC2.SIG.ARPC.up.Amphicrine <- subset(res.dds.ARPC.Amphicrine,res.dds.ARPC.Amphicrine$padj <= 0.01&(res.dds.ARPC.Amphicrine$log2FoldChange)>=1)

res.fdr01.FC2.SIG.up.genes <- list('DNPC'=rownames(res.fdr01.FC2.SIG.ARPC.up.DNPC),'ARLPC'=rownames(res.fdr01.FC2.SIG.ARPC.up.ARLPC),
                                   'SCNPC'=rownames(res.fdr01.FC2.SIG.ARPC.up.SCNPC),'Amphicrine'=rownames(res.fdr01.FC2.SIG.ARPC.up.Amphicrine))
res.fdr01.FC2.SIG.up.values <- list('DNPC'=(res.fdr01.FC2.SIG.ARPC.up.DNPC),'ARLPC'=(res.fdr01.FC2.SIG.ARPC.up.ARLPC),
                                    'SCNPC'=(res.fdr01.FC2.SIG.ARPC.up.SCNPC),'Amphicrine'=(res.fdr01.FC2.SIG.ARPC.up.Amphicrine))


##########
######## DOWNREGULATED

res.fdr01.FC2.SIG.ARPC.dn.ARLPC <- subset(res.dds.ARPC.ARLPC,res.dds.ARPC.ARLPC$padj <= 0.01&(res.dds.ARPC.ARLPC$log2FoldChange)<=-1)
res.fdr01.FC2.SIG.ARPC.dn.DNPC <- subset(res.dds.ARPC.DNPC,res.dds.ARPC.DNPC$padj <= 0.01&(res.dds.ARPC.DNPC$log2FoldChange)<=-1)
res.fdr01.FC2.SIG.ARPC.dn.SCNPC <- subset(res.dds.ARPC.SCNPC,res.dds.ARPC.SCNPC$padj <= 0.01&(res.dds.ARPC.SCNPC$log2FoldChange)<=-1)
res.fdr01.FC2.SIG.ARPC.dn.Amphicrine <- subset(res.dds.ARPC.Amphicrine,res.dds.ARPC.Amphicrine$padj <= 0.01&(res.dds.ARPC.Amphicrine$log2FoldChange)<=-1)

res.fdr01.FC2.SIG.dn.genes <- list('DNPC'=rownames(res.fdr01.FC2.SIG.ARPC.dn.DNPC),'ARLPC'=rownames(res.fdr01.FC2.SIG.ARPC.dn.ARLPC),
                                   'SCNPC'=rownames(res.fdr01.FC2.SIG.ARPC.dn.SCNPC),'Amphicrine'=rownames(res.fdr01.FC2.SIG.ARPC.dn.Amphicrine))
res.fdr01.FC2.SIG.dn.values <- list('DNPC'=(res.fdr01.FC2.SIG.ARPC.dn.DNPC),'ARLPC'=(res.fdr01.FC2.SIG.ARPC.dn.ARLPC),
                                    'SCNPC'=(res.fdr01.FC2.SIG.ARPC.dn.SCNPC),'Amphicrine'=(res.fdr01.FC2.SIG.ARPC.dn.Amphicrine))

# names(res.fdr01.FC2.SIG.dn.genes)
# [1] "DNPC"       "ARLPC"      "SCNPC"      "Amphicrine"

unique.genes <- VennDiagram::get.venn.partitions(res.fdr01.FC2.SIG.up.genes)

fdr01.FC2.unique.up.genes <- list(DNPC_only=unique.genes[unique.genes$DNPC==T&unique.genes$ARLPC==F&unique.genes$SCNPC==F&unique.genes$Amphicrine==F,]$..values..,
                                                ARlow_only=unique.genes[unique.genes$DNPC==F&unique.genes$ARLPC==T&unique.genes$SCNPC==F&unique.genes$Amphicrine==F,]$..values..,
                                                NE_only=unique.genes[unique.genes$DNPC==F&unique.genes$ARLPC==F&unique.genes$SCNPC==T&unique.genes$Amphicrine==F,]$..values..,
                                                Amphicrine_only=unique.genes[unique.genes$DNPC==F&unique.genes$ARLPC==F&unique.genes$SCNPC==F&unique.genes$Amphicrine==T,]$..values..)

unique.genes <- VennDiagram::get.venn.partitions(res.fdr01.FC2.SIG.dn.genes)

fdr01.FC2.unique.dn.genes <- list(DNPC_only=unique.genes[unique.genes$DNPC==T&unique.genes$ARLPC==F&unique.genes$SCNPC==F&unique.genes$Amphicrine==F,]$..values..,
                                                ARlow_only=unique.genes[unique.genes$DNPC==F&unique.genes$ARLPC==T&unique.genes$SCNPC==F&unique.genes$Amphicrine==F,]$..values..,
                                                NE_only=unique.genes[unique.genes$DNPC==F&unique.genes$ARLPC==F&unique.genes$SCNPC==T&unique.genes$Amphicrine==F,]$..values..,
                                                Amphicrine_only=unique.genes[unique.genes$DNPC==F&unique.genes$ARLPC==F&unique.genes$SCNPC==F&unique.genes$Amphicrine==T,]$..values..)


######## Volcano plots
#### DNPC

DNPC.genes <- unique(unlist(c(fdr01.FC2.unique.dn.genes$DNPC_only$`15`,fdr01.FC2.unique.up.genes$DNPC_only$`15`)))

vol.DNPC <- na.omit(as.data.frame(res.dds.ARPC.DNPC))
vol.DNPC$log10p <- -log10(vol.DNPC$padj)

vol.DNPC$Status = ifelse(vol.DNPC$log10p > 2 & vol.DNPC$log2FoldChange >= 1,'Up',
                         ifelse(vol.DNPC$log10p > 2 & vol.DNPC$log2FoldChange < -1,'Down',
                                'Stable'))
id.labels.DNPC <- which(vol.DNPC$log10p > 2&abs(vol.DNPC$log2FoldChange)>1&rownames(vol.DNPC)%in% c(prot.genes$name,NEPC.gene.list$Nelson$GENE_ID,NEPC.gene.list$beltran$GENE_ID))

vol.DNPC$names <- NA;vol.DNPC[id.labels.DNPC,]$names <- rownames(vol.DNPC[id.labels.DNPC,])

vol.DNPC$Specific <- factor(ifelse(vol.DNPC$log10p > 2&abs(vol.DNPC$log2FoldChange)>1&vol.DNPC$names%in%DNPC.genes,'Yes','No'),levels=c('Yes','No'))

set.seed(100)
volcanoplot.DNPC <- ggplot(data = vol.DNPC,
                           aes(x = log2FoldChange,
                               y = vol.DNPC$log10p,
                               colour=Status,label=names)) +
    geom_point(alpha=0.8, size=2) +
    geom_label_repel(max.overlaps = 10,na.rm = T,aes(fill=Specific),color='black',min.segment.length = unit(0, 'lines'),
                     segment.size  = 0.2,
                     segment.color = "grey50",
                     direction     = "x",)+
    scale_color_manual(values=c("blue", "grey","red"))+
    scale_fill_manual(values=c('gold','white'),breaks = c('Yes','No'),name='Subtype specific')+
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = 2,lty=4,col="black",lwd=0.8) +
    labs(x=expression(Log[2]~FC),
         y=expression(-Log[10]~FDR),
         title="AR-/NE-")  +
    theme_bw(base_size = 15,base_rect_size = 2)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="right")


#### SCNPC
SCNPC.genes <- unique(unlist(c(fdr01.FC2.unique.dn.genes$NE_only$`12`,fdr01.FC2.unique.up.genes$NE_only$`12`)))

vol.SCNPC <- na.omit(as.data.frame(res.dds.ARPC.SCNPC))
vol.SCNPC$log10p <- -log10(vol.SCNPC$padj)

vol.SCNPC$Status = ifelse(vol.SCNPC$log10p > 2 & vol.SCNPC$log2FoldChange >= 1,'Up',
                          ifelse(vol.SCNPC$log10p > 2 & vol.SCNPC$log2FoldChange < -1,'Down',
                                 'Stable'))
id.labels.SCNPC <- which(vol.SCNPC$log10p > 2&abs(vol.SCNPC$log2FoldChange)>1&rownames(vol.SCNPC)%in% c(prot.genes$name,NEPC.gene.list$Nelson$GENE_ID,NEPC.gene.list$beltran$GENE_ID))

vol.SCNPC$names <- NA;vol.SCNPC[id.labels.SCNPC,]$names <- rownames(vol.SCNPC[id.labels.SCNPC,])

vol.SCNPC$Specific <- factor(ifelse(vol.SCNPC$log10p > 2&abs(vol.SCNPC$log2FoldChange)>1&vol.SCNPC$names%in%SCNPC.genes,'Yes','No'),levels=c('Yes','No'))

set.seed(100)
volcanoplot.SCNPC <- ggplot(data = vol.SCNPC,
                            aes(x = log2FoldChange,
                                y = vol.SCNPC$log10p,
                                colour=Status,label=names)) +
    geom_point(alpha=0.8, size=2) +
    geom_label_repel(max.overlaps = 10,na.rm = T,aes(fill=Specific),color='black',min.segment.length = unit(0, 'lines'),
                     segment.size  = 0.2,
                     segment.color = "grey50",
                     direction     = "x",)+
    scale_color_manual(values=c("blue", "grey","red"))+
    scale_fill_manual(values=c('gold','white'),breaks = c('Yes','No'),name='Subtype specific')+
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = 2,lty=4,col="black",lwd=0.8) +
    labs(x=expression(Log[2]~FC),
         y=expression(-Log[10]~FDR),
         title="AR-/NE+")  +
    theme_bw(base_size = 15,base_rect_size = 2)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="right")



#### ARLPC

ARLPC.genes <- unique(unlist(c(fdr01.FC2.unique.dn.genes$ARlow_only$`14`,fdr01.FC2.unique.up.genes$ARlow_only$`14`)))

vol.ARLPC <- na.omit(as.data.frame(res.dds.ARPC.ARLPC))
vol.ARLPC$log10p <- -log10(vol.ARLPC$padj)

vol.ARLPC$Status = ifelse(vol.ARLPC$log10p > 2 & vol.ARLPC$log2FoldChange >= 1,'Up',
                          ifelse(vol.ARLPC$log10p > 2 & vol.ARLPC$log2FoldChange < -1,'Down',
                                 'Stable'))
id.labels.ARLPC <- which(vol.ARLPC$log10p > 2&abs(vol.ARLPC$log2FoldChange)>1&rownames(vol.ARLPC)%in% c(prot.genes$name,NEPC.gene.list$Nelson$GENE_ID,NEPC.gene.list$beltran$GENE_ID))

vol.ARLPC$names <- NA;vol.ARLPC[id.labels.ARLPC,]$names <- rownames(vol.ARLPC[id.labels.ARLPC,])

vol.ARLPC$Specific <- factor(ifelse(vol.ARLPC$log10p > 2&abs(vol.ARLPC$log2FoldChange)>1&vol.ARLPC$names%in%ARLPC.genes,'Yes','No'),levels=c('Yes','No'))

set.seed(100)
volcanoplot.ARLPC <- ggplot(data = vol.ARLPC,
                            aes(x = log2FoldChange,
                                y = vol.ARLPC$log10p,
                                colour=Status,label=names)) +
    geom_point(alpha=0.8, size=2) +
    geom_label_repel(max.overlaps = 10,na.rm = T,aes(fill=Specific),color='black',min.segment.length = unit(0, 'lines'),
                     segment.size  = 0.2,
                     segment.color = "grey50",
                     direction     = "x",)+
    scale_color_manual(values=c("blue", "grey","red"))+
    scale_fill_manual(values=c('gold','white'),breaks = c('Yes','No'),name='Subtype specific')+
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = 2,lty=4,col="black",lwd=0.8) +
    labs(x=expression(Log[2]~FC),
         y=expression(-Log[10]~FDR),
         title="ARL/NE-")  +
    theme_bw(base_size = 15,base_rect_size = 2)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="right")


### AR+/NE+
Amphicrine.genes <- unique(unlist(c(fdr01.FC2.unique.dn.genes$Amphicrine_only$`8`,fdr01.FC2.unique.up.genes$Amphicrine_only$`8`)))

vol.Amphicrine <- na.omit(as.data.frame(res.dds.ARPC.Amphicrine))
vol.Amphicrine$log10p <- -log10(vol.Amphicrine$padj)

vol.Amphicrine$Status = ifelse(vol.Amphicrine$log10p > 2 & vol.Amphicrine$log2FoldChange >= 1,'Up',
                               ifelse(vol.Amphicrine$log10p > 2 & vol.Amphicrine$log2FoldChange < -1,'Down',
                                      'Stable'))
id.labels.Amphicrine <- which(vol.Amphicrine$log10p > 2&abs(vol.Amphicrine$log2FoldChange)>1&rownames(vol.Amphicrine)%in% c(prot.genes$name,NEPC.gene.list$Nelson$GENE_ID,NEPC.gene.list$beltran$GENE_ID))

vol.Amphicrine$names <- NA;vol.Amphicrine[id.labels.Amphicrine,]$names <- rownames(vol.Amphicrine[id.labels.Amphicrine,])

vol.Amphicrine$Specific <- factor(ifelse(vol.Amphicrine$log10p > 2&abs(vol.Amphicrine$log2FoldChange)>1&vol.Amphicrine$names%in%Amphicrine.genes,'Yes','No'),levels=c('Yes','No'))

set.seed(100)
volcanoplot.Amphicrine <- ggplot(data = vol.Amphicrine,
                                 aes(x = log2FoldChange,
                                     y = vol.Amphicrine$log10p,
                                     colour=Status,label=names)) +
    geom_point(alpha=0.8, size=2) +
    geom_label_repel(max.overlaps = 10,na.rm = T,aes(fill=Specific),color='black',min.segment.length = unit(0, 'lines'),
                     segment.size  = 0.2,
                     segment.color = "grey50",
                     direction     = "x",)+
    scale_color_manual(values=c("blue", "grey","red"))+
    scale_fill_manual(values=c('gold','white'),breaks = c('Yes','No'),name='Subtype specific')+
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = 2,lty=4,col="black",lwd=0.8) +
    labs(x=expression(Log[2]~FC),
         y=expression(-Log[10]~FDR),
         title="AR+/NE+")  +
    theme_bw(base_size = 15,base_rect_size = 2)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="right")


#### number of up/down regulated genes in each subtype

bar.df.genes <- rbind(data.frame(Subtypes=c('AR-/NE-','ARL/NE-','AR-/NE+','AR+/NE+'),status='Up',numbers=unlist(lapply(res.fdr01.FC2.SIG.up.genes,length))),
                      data.frame(Subtypes=c('AR-/NE-','ARL/NE-','AR-/NE+','AR+/NE+'),status='Down',numbers=unlist(lapply(res.fdr01.FC2.SIG.dn.genes,length))))

bar.df.genes$Subtypes <- factor(bar.df.genes$Subtypes,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE+'))

DE.genes.numbers <- ggbarplot(bar.df.genes,x='Subtypes',y='numbers',fill = 'status',ylab = 'No. DE genes')+
    scale_fill_manual(breaks = c("Up","Down"),values = c("red","blue"))+
    theme_bw(base_size = 15,base_rect_size = 2)+
    theme(panel.spacing.x = unit(1, "mm"),
          plot.title = element_blank(),
          axis.title.y = element_text(size=10),
          axis.text.x = element_text(angle = 0,size=15),
          panel.grid = element_blank(),
          strip.text = element_text(size=8),
          strip.background  = element_blank(),legend.position = 'none')


SFigure3_volcano <- ggarrange(volcanoplot.DNPC+theme(panel.grid = element_blank()),
                          volcanoplot.SCNPC+theme(axis.title.y = element_blank(),panel.grid = element_blank()),
                          volcanoplot.ARLPC+theme(axis.title.y = element_blank(),panel.grid = element_blank()),
                          volcanoplot.Amphicrine+theme(panel.grid = element_blank()),
                          DE.genes.numbers,
                          nrow=2,ncol=3,common.legend = T,legend = 'none',align = 'hv')

###### Supplementary figure 4
### Biological processes GSEA

### DNPC

DNPC.genes <- rbind(res.fdr01.FC2.SIG.up.values$DNPC[fdr01.FC2.unique.up.genes$DNPC_only$`15`,],
                    res.fdr01.FC2.SIG.dn.values$DNPC[fdr01.FC2.unique.dn.genes$DNPC_only$`15`,])
DNPC.vector <- DNPC.genes[order(DNPC.genes$log2FoldChange,decreasing = T),2]
names(DNPC.vector) <- rownames(DNPC.genes[order(DNPC.genes$log2FoldChange,decreasing = T),])

set.seed(100)
gse.DNPC <- gseGO(geneList=DNPC.vector,
                  ont ="BP",
                  keyType = "SYMBOL",seed = 100,
                  nPerm = 5000,
                  minGSSize = 5,
                  maxGSSize = 800,
                  pvalueCutoff = 0.1,
                  verbose = TRUE,
                  OrgDb = "org.Hs.eg.db",pAdjustMethod = "BH")

DNPC.dotplot <- dotplot(gse.DNPC, showCategory = 10, title = "AR-/NE-" , split=".sign",color='p.adjust') +
    facet_grid(.~.sign)+scale_color_gradient(high = "azure3",low = subtype.colors[subtype.colors$subtypes=='AR-/NE-',]$colors)+
    scale_size(range = c(3,9),breaks = c(25,50,75))+
    theme_bw(base_rect_size = 1.5,base_size = 15)+
    theme(axis.title.y = element_text(size=15,face = 'bold'),axis.title.x = element_text(size=15,face = 'bold'),
          axis.text.y = element_text(angle=15,vjust = 0.5,size=10,face = 'bold'),title = element_text(face='bold'),
          axis.text.x = element_text(size=15),strip.text = element_text(face='bold',size=15),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())


### ARLPC
## no enrichment

### SCNPC

SCNPC.genes <- rbind(res.fdr01.FC2.SIG.up.values$SCNPC[fdr01.FC2.unique.up.genes$NE_only$`12`,],
                     res.fdr01.FC2.SIG.dn.values$SCNPC[fdr01.FC2.unique.dn.genes$NE_only$`12`,])
SCNPC.vector <- SCNPC.genes[order(SCNPC.genes$log2FoldChange,decreasing = T),2]
names(SCNPC.vector) <- rownames(SCNPC.genes[order(SCNPC.genes$log2FoldChange,decreasing = T),])

set.seed(100)
gse.SCNPC <- gseGO(geneList=SCNPC.vector,
                   ont ="BP",
                   keyType = "SYMBOL",
                   nPerm = 5000,
                   minGSSize = 5,
                   maxGSSize = 800,seed = 100,
                   pvalueCutoff = 0.1,
                   verbose = TRUE,
                   OrgDb = "org.Hs.eg.db",pAdjustMethod = "BH")

SCNPC.dotplot <- dotplot(gse.SCNPC, showCategory = 10, title = "AR-/NE+" , split=".sign",color='p.adjust') +
    facet_grid(.~.sign)+scale_color_gradient(high = "azure3",low = subtype.colors[subtype.colors$subtypes=='AR-/NE+',]$colors)+
    scale_size(range = c(3,9),breaks = c(25,50,75))+
    theme_bw(base_rect_size = 1.5,base_size = 15)+
    theme(axis.title.y = element_text(size=15,face = 'bold'),axis.title.x = element_text(size=15,face = 'bold'),
          axis.text.y = element_text(angle=15,vjust = 0.5,size=10,face = 'bold'),title = element_text(face='bold'),
          axis.text.x = element_text(size=15),strip.text = element_text(face='bold',size=15),
          strip.background =element_rect(fill='white',color='white'),
          panel.grid = element_blank())


### Amp
## no enrichment


SFigure4.plot <- plot_grid(DNPC.dotplot,SCNPC.dotplot,align = 'hv',ncol = 2)

