library(rowr)
library(annotatr)
library(mutSignatures)
library(kableExtra)
library(deconstructSigs)
library(plyranges)
library(tidyr)
library(rcompanion)
library(plyr)
library(ggthemes)
library(openxlsx)
library(ggpubr)
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(BSgenome.Hsapiens.UCSC.hg38)
library(CHORD)
library(foreach)
library(data.table)

setwd("~/Box Sync/UCSF/Manuscript/Submissions/GitHub/")
source("./functions.R")

chromosomes = c(paste('chr', 1:22, sep=''), 'chrX', 'chrY')
genes_info <- sort(builtin_annotations())[c(186:196)]
# [1] "hg38_genes_3UTRs"                "hg38_genes_5UTRs"                "hg38_genes_cds"                  "hg38_genes_exonintronboundaries" "hg38_genes_exons"
# [6] "hg38_genes_firstexons"           "hg38_genes_intergenic"           "hg38_genes_intronexonboundaries" "hg38_genes_introns"              "hg38_genes_promoters"
# [11] "hg38_lncrna_gencode"

genes.GR <- build_annotations(genome = 'hg38',annotations = genes_info)
genes.GR <- genes.GR[genes.GR@seqnames%in%chromosomes]

# loading data
load(file="./data/2022_11_04_WCDT_poppy.Rdata")
load(file = "./data/NEPC.gene.list.RData")
load(file="./data/WCDT.5groups.210.RData")
load(file="./data/Matched.WGS.WGBS.samples.RData")
load(file = "./data/HRD.status.RData")

Nelson.samples.PCA <- Nelson.samples.PCA[Nelson.samples.PCA$samples%in%Matched.WGS.WGBS.samples$samples,]
Nelson.samples.PCA$Subtypes <- factor(Nelson.samples.PCA$Subtypes,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE-','AR+/NE+'))
Nelson.samples.PCA <- Nelson.samples.PCA[order(Nelson.samples.PCA$Subtypes),]

##### Generated HRD status using CHORD
# Pan-cancer landscape of homologous recombination deficiency
# Luan Nguyen, John Martens, Arne Van Hoeck, Edwin Cuppen. Nat Commun 11, 5584 (2020). https://www.nature.com/articles/s41467-020-19406-4
######
# HRD.status <- foreach(i=names(WGS.available),.combine = 'rbind')%do%{
#     SNV.df <-  WGS.available[[i]]$strelka_mutect[!complete.cases(WGS.available[[i]]$strelka_mutect[c("DP_indel_ref_T","DP_indel_alt_T")]),c('chrom','pos','ref','alt')]
#     indel.df <- WGS.available[[i]]$strelka_mutect[complete.cases(WGS.available[[i]]$strelka_mutect[c("DP_indel_ref_T","DP_indel_alt_T")]),c('chrom','pos','ref','alt')]
#     SV.df <- WGS.available[[i]]$list_sv_gripss[,c('svtype','sv_length')]
#     contexts <- extractSigsChord(df.snv = SNV.df,df.indel = indel.df,df.sv = SV.df,sample.name = i,
#                                  sv.caller = 'gridss',ref.genome = BSgenome.Hsapiens.UCSC.hg38,verbose = F)
#     HRD.df <- chordPredict(contexts,verbose = F)
#     HRD.df}
#
# HRD.df <- merge(HRD.status,Nelson.samples.PCA[,-3],by.x='sample',by.y='samples')
# save(file = "./data/HRD.status.RData",HRD.df)
######

#### using gene list from cell 2018 Quigley et al., oncogenic driver of prostate cancer

genelist1 <- read.table("./data/GENE.LIST.CELL.txt",header = T,sep = "\t")
genelist2 <- read.xlsx("./data/prostate_cancer.genes.xlsx",sheet = 2)
genes.of.interest <- unique(c(genelist2$Gene,genelist1$GENES))

WGS.available <- SO_WCDT_WGS[names(SO_WCDT_WGS)%in%Nelson.samples.PCA$samples]

WGS.info <- as.data.frame(foreach(i=names(WGS.available),.combine = "rbind") %do% {
    WG_dup <- WGS.available[[i]]$wholeGenomeDuplication
    msStatus <- WGS.available[[i]]$msStatus

    ploidy <- WGS.available[[i]]$ploidy
    purity <- WGS.available[[i]]$purity*100
    tmb_MB <- WGS.available[[i]]$tmbPerMb
    tmb_MB_log2 <- log2(WGS.available[[i]]$tmbPerMb)
    microHomo <- log2(WGS.available[[i]]$nt_for_microhomology)
    dip.pur <- WGS.available[[i]]$diploidProportion
    c(WG_dup=WG_dup,msStatus=msStatus,ploidy=ploidy,purity=purity,tmb_MB=tmb_MB,tmb_MB_log2=tmb_MB_log2,microHomo=microHomo,
      dip.pur=dip.pur)});WGS.info <- tibble::add_column(.data = WGS.info,PATIENT_ID=stringr::str_extract(names(WGS.available), "[^-]*-[^-]*"),SAMPLE_ID=names(WGS.available),.before = "WG_dup")

WGS.info <- merge(WGS.info,Nelson.samples.PCA,by.x='SAMPLE_ID',by.y='samples',all.x=T,sort=F)
WGS.info$Subtypes <- factor(WGS.info$Subtypes,levels = c("AR-/NE-","AR-/NE+","ARL/NE-","AR+/NE-","AR+/NE+"))

WGS.info[,c("ploidy","purity","tmb_MB","tmb_MB_log2","microHomo","dip.pur")] = apply(WGS.info[,c("ploidy","purity","tmb_MB","tmb_MB_log2","microHomo","dip.pur")], 2, function(x) as.numeric(as.character(x)))
WGS.info$WG_dup <- factor(WGS.info$WG_dup,levels=c('true','false'),labels = c('Yes','No'));WGS.info$msStatus <- factor(WGS.info$msStatus)
WGS.info <- WGS.info[order(WGS.info$Subtypes),]

gridss.output <- foreach(i=names(WGS.available)) %do% {
    table(WGS.available[[i]]$list_sv_gripss$svtype)};names(gridss.output) <- names(WGS.available)

gridss.df <- as.data.frame(bind_rows(gridss.output));rownames(gridss.df) <- names(WGS.available)
gridss.df[is.na(gridss.df)] <- 0;gridss.df <- gridss.df[WGS.info$SAMPLE_ID,]

# "BND" "DEL" "DUP" "INS" "INV" "SGL"
# "Interchromosomal translocation","Deletion","Duplication","Insertions","Inversions","Single break end"

table(WGS.info$Subtypes)
# AR-/NE- AR-/NE+ ARL/NE- AR+/NE- AR+/NE+
#     7       6      32      77       6

##### Mutation matrix

subset.onco <- list()
for(i in 1:length(names(WGS.available))){
    sample_id <- names(WGS.available)[i]
    cat(format(Sys.time(), "%a %b %d %X %Y"), "START:", sample_id, "\n", sep=" ")
    strelka_mutect <- subset(WGS.available[[i]]$strelka_mutect,WGS.available[[i]]$strelka_mutect$gene %in% genes.of.interest)[c('gene','effect')]
    subset.onco[[sample_id]]$strelka_mutect <- strelka_mutect
    cat(format(Sys.time(), "%a %b %d %X %Y"), "PROCESSED:", sample_id, "\n", sep=" ")}

## only choose these mutational impacts, from Sequence ontology

mut.impact <- c("transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained","transcript_amplification","frameshift_variant","stop_lost","start_lost", # high impact
                "inframe_insertion","inframe_deletion","missense_variant",
                "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant","frameshift_variant&splice_region_variant",
                "missense_variant&splice_region_variant",
                "splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant",
                "splice_acceptor_variant&intron_variant",
                "splice_acceptor_variant&splice_region_variant&intron_variant&non_coding_transcript_exon_variant",
                "splice_donor_variant&intron_variant")

mat.df <- matrix(ncol = 1, nrow=length(genes.of.interest));mat.df[,1] = genes.of.interest
colnames(mat.df) <- 'gene'
foreach(i=names(subset.onco),.combine = 'cbind.fill') %do% {subset.onco[[i]]$strelka_mutect$effect[subset.onco[[i]]$strelka_mutect$effect %out% mut.impact] <- NA}

mat.df <- foreach(i=names(subset.onco)) %do% {
    onco.df <- subset.onco[[i]]$strelka_mutect[c('gene','effect')] %>% group_by(effect,gene) %>% dplyr::summarise() %>% as.data.frame()}
mat.df <- lapply(mat.df,function(x) x[!duplicated(x$gene),])

names(mat.df) <- names(subset.onco);tmp.df <- data.frame(gene=genes.of.interest)

mat.onco <- foreach(i=names(mat.df)) %do% {
    merge(tmp.df,mat.df[[i]],by='gene',all.x=T)}
names(mat.onco) <- names(subset.onco)
mat.onco <- do.call(cbind,mat.onco)

### rename row names with gene names
rownames(mat.onco) = mat.onco$`DTB-005-BL.gene`;mat.onco <-  mat.onco[,c(F,T)];colnames(mat.onco) = names(WGS.available)

### list of alterations in the WCDT samples

# names(table(unlist(mat.onco)))
# [1] "frameshift_variant"
# [2] "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant"
# [3] "frameshift_variant&splice_region_variant"
# [4] "missense_variant"
# [5] "missense_variant&splice_region_variant"
# [6] "splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant"
# [7] "splice_acceptor_variant&intron_variant"
# [8] "splice_donor_variant&intron_variant"
# [9] "stop_gained"
# [10] "stop_lost"


mat.onco[is.na(mat.onco)] <- ""
mat.onco[(mat.onco)=='frameshift_variant'] <- "Frameshift"
mat.onco[(mat.onco)=='frameshift_variant&splice_region_variant'] <- "Frameshift"
mat.onco[(mat.onco)=='frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant'] <- "Frameshift"
mat.onco[(mat.onco)=='missense_variant'] <- "Missense"
mat.onco[(mat.onco)=='missense_variant&splice_region_variant'] <- "Missense"
mat.onco[(mat.onco)=='splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant'] <- "Splice acceptor"
mat.onco[(mat.onco)=='splice_acceptor_variant&intron_variant'] <- "Splice acceptor"
mat.onco[(mat.onco)=='splice_donor_variant&intron_variant'] <- "Splice donor"
mat.onco[(mat.onco)=='stop_gained'] <- "Stop gained"
mat.onco[(mat.onco)=='stop_lost'] <- "Stop lost"


### Amplification matrix

sex.chr <- c('chrX','chrY')
Amp.matrix <- list()

purple.genes <- data.frame(genes=rownames(WGS.available$`DTB-003-BL`$CNA_genes))
purple.genes <- merge(purple.genes,unique(as.data.frame(genes.GR[genes.GR$symbol%in%genes.of.interest])[,c(1,9)]),
                      by.x='genes',by.y='symbol',sort=F)
purple.genes <- purple.genes[,c(2,1)]
colnames(purple.genes) <- c('chrom','genes')

for(i in 1:length(names(WGS.available))){
    sample_id <- names(WGS.available)[i]
    cat(format(Sys.time(), "%a %b %d %X %Y"), "START:", sample_id, "\n", sep=" ")
    Amp <- data.frame(purple.genes,WGS.available[[i]]$CNA_genes[purple.genes$genes,])
    Amp.matrix[[sample_id]] <- Amp
    Amp.matrix[[i]]$ploidy <- WGS.available[[i]]$ploidy

    Amp.matrix[[i]]$Amp.stat <-
        ifelse(Amp.matrix[[i]]$chrom %in% sex.chr & Amp.matrix[[i]]$weighted_mean > Amp.matrix[[i]]$ploidy*0.9,'AMP',
               ifelse(Amp.matrix[[i]]$chrom %out% sex.chr &Amp.matrix[[i]]$weighted_mean > Amp.matrix[[i]]$ploidy*1.95,'AMP',
                      ifelse(Amp.matrix[[i]]$chrom %in% sex.chr & Amp.matrix[[i]]$weighted_mean < 0.75,'DEL',
                             ifelse(Amp.matrix[[i]]$chrom %out% sex.chr & Amp.matrix[[i]]$weighted_mean < 1.1,'DEL',''))))
    Amp.matrix[[i]]$Biallelic <- ifelse(Amp.matrix[[i]]$weighted_mean < 0.5,'biallelic','')
    cat(format(Sys.time(), "%a %b %d %X %Y"), "PROCESSED:", sample_id, "\n", sep=" ")}

## subset by genes of interest
subset.AMP <- lapply(Amp.matrix,function(x) x[x$gene%in%genes.of.interest,])
names(subset.AMP) <- names(Amp.matrix)


mat.AMP <- matrix(ncol = 1,nrow=length(purple.genes$genes));mat.AMP[,1] = purple.genes$genes;colnames(mat.AMP) <- 'genes'
mat.AMP <- foreach(i=names(subset.AMP),.combine = 'cbind.fill') %do% {merge(mat.AMP,subset.AMP[[i]][,c(2,6:7)],by='genes')}
rownames(mat.AMP) = mat.AMP$genes;mat.AMP <-  mat.AMP[,c(F,T,T)];colnames(mat.AMP) = rep(names(subset.AMP),each=2)
gene.idx <- intersect(rownames(mat.onco),rownames(mat.AMP))

### combine Amplification and somatic alteration data

mat.AMP.united <- foreach(i=1:ncol(mat.onco),.combine = 'cbind') %do%{
    df <- cbind(mat.AMP[gene.idx,c(T,F)][i],mat.AMP[gene.idx,c(F,T)][i])
    colnames(df) <- c('AMP','BI')
    df[df == ""]<- NA
    df <- df %>% unite('test',sep = ';',remove = T,na.rm = T)
};colnames(mat.AMP.united) <- colnames(mat.onco)

mut.matrix <- foreach(i=1:ncol(mat.onco),.combine = 'cbind') %do%{
    df <- cbind(mat.onco[gene.idx,i],mat.AMP[gene.idx,c(T,F)][i],mat.AMP[gene.idx,c(F,T)][i])
    colnames(df) <- c('ONCO','AMP','BI')
    df[df == ""]<- NA
    df <- df %>% unite('test',sep = ';',remove = T,na.rm = T)
};colnames(mut.matrix) <- colnames(mat.onco)


### adding AR.enhancer to the matrix - AR enhancer region reported in cell 2018 and Nature genetics Zhao et al.
AR.enhance.GR <- GRanges(seqnames="chrX",ranges=IRanges(start = 66895158, end = 66910158),strand="*")
ARenhancer <- lapply(WGS.available,function(x) makeGRangesFromDataFrame(x$segments,keep.extra.columns = T))
ARenhancer <- endoapply(ARenhancer,subsetByOverlaps,AR.enhance.GR)

for(i in 1:length(ARenhancer)){
    ARenhancer[[i]]$ploidy <- WGS.available[[i]]$ploidy
    ARenhancer[[i]]$Amp.stat <-
        ifelse(ARenhancer[[i]]$copyNumber > ARenhancer[[i]]$ploidy*0.9,'AMP',
               ifelse(Amp.matrix[[i]]$chrom %in% sex.chr & Amp.matrix[[i]]$maxCopyNumber < 0.75 , 'DEL',''))}

### DTB-261-BL has two reports in the region two different locations within that region

ARenhancer <- unlist(sapply(ARenhancer,function(x) unique(x$Amp.stat)))[-114]
names(ARenhancer)[114] <- "DTB-261-BL"
mut.matrix <- rbind(mut.matrix,"AR.enhancer"=ARenhancer)


#annoData <- toGRanges(EnsDb.Hsapiens.v86,feature='exon')
# excluded SGLs
SV.output <- foreach(i=names(WGS.available)) %do% {
    WGS.available[[i]]$list_sv_gripss[WGS.available[[i]]$list_sv_gripss$svtype!='SGL',]};names(SV.output) <- names(WGS.available)
fusion.output <- foreach(i=names(WGS.available)) %do% {
    WGS.available[[i]]$fusions_linx};names(fusion.output) <- names(WGS.available)
fusion.output <- foreach(i=names(WGS.available)) %do% {
    fusion.output[[i]][fusion.output[[i]]$geneStart%in%c('ETV1','ETV4','ETV5','ERG')|fusion.output[[i]]$geneEnd%in%c('ETV1','ETV4','ETV5','ERG'),]}
names(fusion.output) <- names(WGS.available)


fusion.list <- foreach(i=names(WGS.available)) %do% {
    df <- fusion.output[[i]][,'name']};names(fusion.list) <- names(WGS.available)
fusion.list <- data.frame(fusion=unlist(fusion.list))

### manually check the fusions if not reported by the algorithm
### genes location hg38, Genecard and NCBI reports
#ERG chr21:38367261-38661783    TMRSS2 chr21:41464300-41531116 ETV1 chr7:13891229-13991425 ETV4 chr17:43527843-43579620 ETV5 chr3:186046314-186110318

fusion.ch3 <- lapply(SV.output,function(x) x[x$chrom_1 %in% "chr3",])
fusion.ch7 <- lapply(SV.output,function(x) x[x$chrom_1 %in% "chr7",])
fusion.chr17 <- lapply(SV.output,function(x) x[x$chrom_1 %in% "chr17",])
fusion.chr21 <- lapply(SV.output,function(x) x[x$chrom_1 %in% "chr21"&x$chrom_2 %in% "chr21",])

ETV1.fusion <- lapply(fusion.ch7,function(x) x[x$start_1>13891000&x$start_1<13991500,])
ETV4.fusion <- lapply(fusion.chr17,function(x) x[x$start_1>43527000&x$start_1<43579650,])
ETV5.fusion <- lapply(fusion.ch3,function(x) x[x$start_1>186046314&x$start_1<186110318,])
ERG.fusion <- lapply(fusion.chr21,function(x) x[x$start_1>37000000&x$start_1<39000000&x$start_2>41000000|x$start_2>37000000&x$start_2<39000000&x$start_1>41000000,])

ETV1.df <- makeGRangesFromDataFrame(do.call(rbind,ETV1.fusion)[c('chrom_2','start_2','end_2','svtype')],seqnames.field = 'chrom_2',start.field = 'start_2',end.field = 'end_2',keep.extra.columns = T)
ETV4.df <- makeGRangesFromDataFrame(do.call(rbind,ETV4.fusion)[c('chrom_2','start_2','end_2','svtype')],seqnames.field = 'chrom_2',start.field = 'start_2',end.field = 'end_2',keep.extra.columns = T)
ETV5.df <- makeGRangesFromDataFrame(do.call(rbind,ETV5.fusion)[c('chrom_2','start_2','end_2','svtype')],seqnames.field = 'chrom_2',start.field = 'start_2',end.field = 'end_2',keep.extra.columns = T)
ERG.df <- makeGRangesFromDataFrame(do.call(rbind,ERG.fusion)[c('chrom_2','start_2','end_2','svtype')],seqnames.field = 'chrom_2',start.field = 'start_2',end.field = 'end_2',keep.extra.columns = T)

#### annotate regions

ETV1.df.fusion <- join_nearest(ETV1.df,genes.GR,distance = T)
ETV1.df.fusion$SAMPLE_ID <- gsub("[.].*","",names(ETV1.df.fusion))
ETV1.df.fusion <- as.data.frame(ETV1.df.fusion)[,c('SAMPLE_ID','symbol','distance','svtype')];ETV1.df.fusion <- ETV1.df.fusion[ETV1.df.fusion$symbol%out%c('ETV1',"")&ETV1.df.fusion$distance<200,]
ETV1.df.fusion <- ETV1.df.fusion[,1:2]
ETV1.df.fusion <- ETV1.df.fusion[!is.na(ETV1.df.fusion$symbol),]
ETV1.df.fusion$fusion <- paste0('ETV1_',ETV1.df.fusion$symbol)

ETV4.df.fusion <- join_nearest(ETV4.df,genes.GR,distance = T)
ETV4.df.fusion$SAMPLE_ID <- gsub("[.].*","",names(ETV4.df.fusion))
ETV4.df.fusion <- as.data.frame(ETV4.df.fusion)[,c('SAMPLE_ID','symbol','distance','svtype')];ETV4.df.fusion <- ETV4.df.fusion[ETV4.df.fusion$symbol%out%c('ETV4',"")&ETV4.df.fusion$distance<200,]
ETV4.df.fusion <- ETV4.df.fusion[,1:2]
ETV4.df.fusion <- ETV4.df.fusion[!is.na(ETV4.df.fusion$symbol),]
ETV4.df.fusion$fusion <- paste0('ETV4_',ETV4.df.fusion$symbol)

ETV5.df.fusion <- join_nearest(ETV5.df,genes.GR,distance = T)
ETV5.df.fusion$SAMPLE_ID <- gsub("[.].*","",names(ETV5.df.fusion))
ETV5.df.fusion <- as.data.frame(ETV5.df.fusion)[,c('SAMPLE_ID','symbol','distance','svtype')];ETV5.df.fusion <- ETV5.df.fusion[ETV5.df.fusion$symbol%out%c('ETV5',"")&ETV5.df.fusion$distance<200,]
ETV5.df.fusion <- ETV5.df.fusion[,1:2]
ETV5.df.fusion <- ETV5.df.fusion[!is.na(ETV5.df.fusion$symbol),]
ETV5.df.fusion$fusion <- paste0('ETV5_',ETV5.df.fusion$symbol)

ERG.df.fusion <- join_nearest(ERG.df,genes.GR,distance = T)
ERG.df.fusion$SAMPLE_ID <- gsub("[.].*","",names(ERG.df.fusion))
ERG.df.fusion <- as.data.frame(ERG.df.fusion)[,c('SAMPLE_ID','symbol','distance','svtype')];ERG.df.fusion <- ERG.df.fusion[ERG.df.fusion$symbol%out%c('ERG',"")&ERG.df.fusion$distance<200,]
ERG.df.fusion <- ERG.df.fusion[,1:2]
ERG.df.fusion <- ERG.df.fusion[!is.na(ERG.df.fusion$symbol),]
ERG.df.fusion$fusion <- paste0('ERG_',ERG.df.fusion$symbol)


##### fusion for each subtype - Supplementary Table 3
##########

fusion.list <- tibble::add_column('SAMPLE_ID'=rownames(fusion.list),.data = fusion.list,.before='fusion')

fusion.df <- rbind(unique(ETV1.df.fusion),
                   unique(ETV4.df.fusion),
                   unique(ETV5.df.fusion),
                   unique(ERG.df.fusion))
fusion.df <- fusion.df[,c(1,3)]

fusion.data <- rbind(fusion.df,fusion.list)
rownames(fusion.data) <- 1:nrow(fusion.data)

fusion.data$SAMPLE_ID <- gsub(fusion.data$SAMPLE_ID,pattern = "^[0-9]|[0-9]$",replacement = "")
fusion.data$gene1 <- gsub(fusion.data$fusion,pattern = "_.*",replacement = "")
fusion.data$gene2 <- gsub(fusion.data$fusion,pattern = ".*_",replacement = "")
fusion.data$gene3 <- ifelse(fusion.data$gene1 %in% c('ETV1','ETV4','ETV5','ERG'),fusion.data$gene1,
                            ifelse(fusion.data$gene2 %in% c('ETV1','ETV4','ETV5','ERG'),fusion.data$gene2,NA))
fusion.data$gene4 <- ifelse(fusion.data$gene1 %out% c('ETV1','ETV4','ETV5','ERG'),fusion.data$gene1,
                            ifelse(fusion.data$gene2 %out% c('ETV1','ETV4','ETV5','ERG'),fusion.data$gene2,NA))

fusion.data$Fusion <- paste0(fusion.data$gene3,"_",fusion.data$gene4);fusion.data <- na.omit(fusion.data)
fusion.data <- fusion.data[,c('SAMPLE_ID','gene3','gene4','Fusion')]
fusion.data <- unique(fusion.data)

fusion.data.frame <- merge(fusion.data,WGS.info[c('SAMPLE_ID','Subtypes')],all.y=T)
fusion.data.frame$TMPRSS2_ERG <-fusion.data.frame$Fusion %in% 'ERG_TMPRSS2'
fusion.data.frame$ETS_fusion <- !is.na(fusion.data.frame$Fusion)

fusion.table <- fusion.data.frame[,-1] %>% group_by(Subtypes,gene3,gene4) %>% dplyr::count(Fusion) %>% drop_na %>% as.data.frame()
fusion.table$Fusion <- paste0(fusion.table$gene4,' (',fusion.table$n,')')
fusion.table <- fusion.table[c('Subtypes','gene3','Fusion')]
fusion.table$Subtypes <- factor(fusion.table$Subtypes,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE-','AR+/NE+'))
fusion.table$Gene <- factor(fusion.table$gene3,levels=c('ETV1','ETV4','ETV5','ERG'))
fusion.table <- fusion.table[c(c('Subtypes','Gene','Fusion'))]
fusion.table <- fusion.table %>% arrange(Subtypes,Gene)

# write.xlsx(x = fusion.table,file = "./tables/STable3.xlsx",asTable = T) - Suppelementary Table 3 excel file

##########

ETS.family.df <- (fusion.data.frame[!duplicated(fusion.data.frame$SAMPLE_ID),c('SAMPLE_ID','ETS_fusion')])

TMPRSS2_tmp <- distinct(fusion.data.frame[,c('SAMPLE_ID','TMPRSS2_ERG')],.keep_all = F)
### samples that have two calls - linx and manual check are not in agreement (one of them is T or F)
idx <- TMPRSS2_tmp[which(duplicated(TMPRSS2_tmp$SAMPLE_ID)),]$SAMPLE_ID
TMPRSS2_tmp[TMPRSS2_tmp$SAMPLE_ID%in%idx,]$TMPRSS2_ERG <- T
TMPRSS2_ERG.df  <- distinct(TMPRSS2_tmp[,c('SAMPLE_ID','TMPRSS2_ERG')],.keep_all = F)

fusion_group <- merge(ETS.family.df,TMPRSS2_ERG.df,by='SAMPLE_ID')
fusion_group <- merge(Nelson.samples.PCA[,-3],fusion_group,by.x='samples',by.y='SAMPLE_ID',all.x=T,sort=T)
colnames(fusion_group)[1] <- 'SAMPLE_ID'

WGS.info <- merge(WGS.info,fusion_group[,-2],by='SAMPLE_ID',sort=F,all.x=T)

##### preparing data for visualization - Figure 2A

########
colside.cols = list(
    "ETS_fusion" = c('TRUE' = 'black','FALSE' = 'white'),
    "ETV4_fusion" = c('TRUE' = 'black','FALSE' = 'white'),
    "TMRSS2_ERG" = c('TRUE' = 'black','FALSE' = 'white'),
    'Binary class' = c('t-SCNC' = 'orange','Adeno'='forestgreen'),
    'Subtypes' = c(structure(c(subtype.colors$colors),
                             names = c(subtype.colors$subtypes))),
    'Biopsy site'=c('Lymph_node'='slateblue4','Bone'='peachpuff3','Other'='green','Liver'='red','No-data'='white'),
    'SV'=c("BND"=ggsci::pal_d3()(5)[1],"DEL"=ggsci::pal_d3()(5)[2],
           "DUP"=ggsci::pal_d3()(5)[3],"INS"=ggsci::pal_d3()(5)[4],"INV"=ggsci::pal_d3()(5)[5]))

WGS.info$tmb_MB.log2 <- log2(WGS.info$tmb_MB)

top.Annot.df = columnAnnotation(
    "ETS_fusion"=WGS.info$ETS_fusion,
    "TMRSS2_ERG"=WGS.info$TMPRSS2_ERG,
    'Tumour purity' = anno_simple(WGS.info$purity, col = colorRamp2(c(0, 100), c("white", "darkslateblue")), na_col = "white", gp = gpar(col = "black")),
    'Tumour ploidy' = anno_simple(WGS.info$ploidy, col = colorRamp2(c(1, 6), c("white", "darkorange2")),na_col = "white", gp = gpar(col = "black")),
    col=colside.cols,annotation_height=unit(c(rep(0.3,4)), "cm"),gp=gpar(col="black",lwd=0.7),show_legend=F,
    annotation_name_gp= gpar(fontsize = 8,lwd=0.4),annotation_name_side = "left")

SV.annot = columnAnnotation("Number of SV"=anno_barplot(gridss.df,gp =gpar(fill=colside.cols$SV,rot=45),bar_width = 0.7, height = unit(3, "cm"),baseline = 0),
                            col=colside.cols$SV,annotation_name_gp= gpar(fontsize = 10),annotation_name_rot=-90)

colors.WGS<- structure(c("white","purple","forestgreen","orange1","sienna3","seagreen2","navajowhite4","firebrick1",
                         "royalblue","black","salmon",'black'),
                       names = c("NULL","Frameshift","Missense","Splice acceptor","Splice donor","Stop gained","Stop lost",
                                 "AMP",
                                 "DEL","biallelic","germline",'HRD'))

main.genes <- c('AR','AR.enhancer','PTEN','RB1','TP53','BRCA2','MYC','CHD7')

### adding germline data on BRCA2 - checked on IGV, matched with GATK4 germline caller

mut.matrix['BRCA2',c('DTB-097-PRO')] <- c('DEL;germline')
mut.matrix['BRCA2',c('DTB-003-BL')] <- c('germline')
HRD.subset <- HRD.df[HRD.df$p_hrd>=0.5,]$sample

mut.matrix['BRCA2',HRD.subset] <- paste0('HRD;',mut.matrix['BRCA2',HRD.subset])
oncoPrint.mat <- mut.matrix[,as.character(WGS.info$SAMPLE_ID)]

# all(WGS.info$SAMPLE_ID==colnames(oncoPrint.mat))
#[1] TRUE

ht_opt$TITLE_PADDING = unit(c(5,5), "points")

Figure2A_heatmap <- oncoPrint(mat = oncoPrint.mat[main.genes,], top_annotation = top.Annot.df,alter_fun = alter_fun,show_pct = T,pct_digits = 0,
                                       col=colors.WGS,remove_empty_rows = T,row_names_gp = gpar(fontsize=8,fontface='bold'),pct_gp = gpar(fontsize=8),
                                       pct_side = 'right',row_names_side = 'left',
                                       row_labels = c('AR','AR enhancer','PTEN','RB1','TP53','BRCA2','MYC','CHD7'),
                                       column_split = WGS.info$Subtypes,column_gap=unit(2,'mm'),
                                       show_heatmap_legend=F,row_title_gp=gpar(fontsize=8,fontface='bold'),
                                       column_order = WGS.info$SAMPLE_ID,show_column_names = F,
                                       column_title_gp = gpar(fill = c("#7570B3","#E7298A","#D95F02","#1B9E77","#66A61E"),fontsize=8,fontface='bold',border='black',col='white'))

########


### Figure 2B
########
sub.matrix <- mut.matrix[c('AR','PTEN','RB1','TP53','MYC','BRCA2','CHD7'),]
tmp.data <- as.data.frame(t(apply(sub.matrix,2,function(x){ifelse(x%in%'DEL;biallelic',"Biallelic Deletion",
                                                                  ifelse(x%in%"AMP","Amplification",
                                                                         ifelse(x%in%"DEL","Deletion",
                                                                                ifelse(x%in%c(""),NA,
                                                                                       ifelse(x%in%c("Missense;AMP"),'Multi-hit_AMP',
                                                                                              ifelse(x%in%c("DEL;germline","Frameshift;DEL","HRD;DEL","HRD;DEL;biallelic","HRD;Frameshift;DEL","HRD;germline",
                                                                                                            "HRD;Missense;DEL","HRD;Stop gained;DEL","Missense;DEL","Splice acceptor;DEL","Splice donor;DEL","Stop gained;DEL"),'Multi-hit_inactive',x))))))})))

colnames(tmp.data) <- rownames(sub.matrix)
tmp.data <- merge(Nelson.samples.PCA[,-3],tmp.data,by.x='samples',by.y='row.names')

tmp.DNPC <- tmp.data[tmp.data$Subtypes%in%'AR-/NE-',]
DNPC.df <- apply(tmp.DNPC[,-c(1:2)],2,table)
DNPC.df <- DNPC.df[lapply(DNPC.df,length)>0];DNPC.df <- ldply(DNPC.df,data.frame);DNPC.df$Freq <- DNPC.df$Freq/(nrow(tmp.DNPC))*100
DNPC.df$Subtype <- 'AR-/NE-'
colnames(DNPC.df) <- c('gene','Alteration','freq','Subtypes')


tmp.SCNPC <- tmp.data[tmp.data$Subtypes%in%'AR-/NE+',]
SCNPC.df <- apply(tmp.SCNPC[,-c(1:2)],2,table)
SCNPC.df <- SCNPC.df[lapply(SCNPC.df,length)>0];SCNPC.df <- ldply(SCNPC.df,data.frame);SCNPC.df$Freq <- SCNPC.df$Freq/(nrow(tmp.SCNPC))*100
SCNPC.df$Subtype <- 'AR-/NE+'
colnames(SCNPC.df) <- c('gene','Alteration','freq','Subtypes')


tmp.ARLPC <- tmp.data[tmp.data$Subtypes%in%'ARL/NE-',]
ARLPC.df <- apply(tmp.ARLPC[,-c(1:2)],2,table)
ARLPC.df <- ARLPC.df[lapply(ARLPC.df,length)>0];ARLPC.df <- ldply(ARLPC.df,data.frame);ARLPC.df$Freq <- ARLPC.df$Freq/(nrow(tmp.ARLPC))*100
ARLPC.df$Subtype <- 'ARL/NE-'
colnames(ARLPC.df) <- c('gene','Alteration','freq','Subtypes')

tmp.AMPC <- tmp.data[tmp.data$Subtypes%in%'AR+/NE+',]
AMPC.df <- apply(tmp.AMPC[,-c(1:2)],2,table)
AMPC.df <- AMPC.df[lapply(AMPC.df,length)>0];AMPC.df <- ldply(AMPC.df,data.frame);AMPC.df$Freq <- AMPC.df$Freq/(nrow(tmp.AMPC))*100
AMPC.df$Subtype <- 'AR+/NE+'
colnames(AMPC.df) <- c('gene','Alteration','freq','Subtypes')

tmp.ARPC <- tmp.data[tmp.data$Subtypes%in%'AR+/NE-',]
ARPC.df <- apply(tmp.ARPC[,-c(1:2)],2,table)
ARPC.df <- ARPC.df[lapply(ARPC.df,length)>0];ARPC.df <- ldply(ARPC.df,data.frame);ARPC.df$Freq <- ARPC.df$Freq/(nrow(tmp.ARPC))*100
ARPC.df$Subtype <- 'AR+/NE-'
colnames(ARPC.df) <- c('gene','Alteration','freq','Subtypes')


mut.freq.df <- rbind(DNPC.df,SCNPC.df,ARLPC.df,ARPC.df,AMPC.df)

mut.freq.df$gene <- factor(mut.freq.df$gene,levels=c('AR','PTEN','RB1','TP53','MYC','BRCA2','CHD7'))
mut.freq.df$Alteration <- factor(mut.freq.df$Alteration,levels=c('Multi-hit_AMP','Multi-hit_inactive','Biallelic Deletion','Amplification','Deletion',
                                                                 "Stop lost","Stop gained","Splice donor","Splice acceptor","Frameshift","Missense"))

Figure2B_barplot <- ggbarplot(mut.freq.df,x='Subtypes',y='freq',legend='none',fill = 'Alteration',ggtheme = theme_bw(base_size = 15,base_rect_size = 1.5),xlab = 'Subtypes',ylab='Alteration Frequency (%)')+
    scale_fill_manual(breaks=levels(mut.freq.df$Alteration),values = c("deeppink","slateblue1","black","firebrick1","royalblue",
                                                                       "navajowhite4","seagreen2","sienna3","orange1","purple","forestgreen"))+
    facet_wrap(~gene,ncol = 4)+
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5),
          strip.background =element_rect(fill='white',color='white'),strip.text = element_text(face = 'bold',size = 15),
          panel.grid = element_blank())
########


##### Supplementary figure 5.
ETV1.df <- do.call(rbind,ETV1.fusion)[c('chrom_1','start_1','end_1','chrom_2','start_2','end_2')]
ETV1.df$SAMPLE_ID <- gsub("[.].*","",rownames(ETV1.df))
ETV4.df <- do.call(rbind,ETV4.fusion)[c('chrom_1','start_1','end_1','chrom_2','start_2','end_2')]
ETV4.df$SAMPLE_ID <- gsub("[.].*","",rownames(ETV4.df))
ETV5.df <- do.call(rbind,ETV5.fusion)[c('chrom_1','start_1','end_1','chrom_2','start_2','end_2')]
ETV5.df$SAMPLE_ID <- gsub("[.].*","",rownames(ETV5.df))
ERG.df <- do.call(rbind,ERG.fusion)[c('chrom_1','start_1','end_1','chrom_2','start_2','end_2')]
ERG.df$SAMPLE_ID <- gsub("[.].*","",rownames(ERG.df))


ETS.fusion.df <- rbind(ETV1.df,ETV4.df,ETV5.df,ERG.df)
ETS.fusion.df <- merge(ETS.fusion.df,Nelson.samples.PCA[c('samples','Subtypes')],by.x='SAMPLE_ID',by.y='samples')
ETS.fusion.df$value <- as.numeric(factor(ETS.fusion.df$Subtypes,levels = c('AR-/NE-','AR-/NE+','AR+/NE-','AR+/NE+','ARL/NE-')))
ETS.fusion.df <- na.omit(ETS.fusion.df)

gene.df <- fread("./data/ensembl2sym.txt",data.table = F)
label.genes <- gene.df[gene.df$name%in%c('ETV1','ETV4','ETV5','ERG',fusion.data.frame$gene4,'CTB-175E5.4'),][c('chr','start','end','name')]
colnames(label.genes) <- c('chr','start','end','label')
label.genes <- label.genes[!duplicated(label.genes$label),]
label.subset <- label.genes[label.genes$label %in% unique(c('ETV1','ETV4','ERG','SLC30A4','TMPRSS2','CTB-175E5.4',unique(fusion.data.frame[fusion.data.frame$Subtypes%in%c('AR-/NE+','AR-/NE-'),]$gene4))),]

### Circos plot - Supplementary Figure 5A
########
circos.clear()
circos.initializeWithIdeogram(plotType =NULL,species = 'hg38')
circos.genomicLabels(label.subset, labels.column = 4, side = "outside",labels.side = "outside",connection_height = 0.08,
                     col = 'black', line_col ='black',niceFacing = T,cex = 0.9)

circos.track(ylim = c(0, 1),panel.fun = function(x, y) {
    chr = CELL_META$sector.numeric.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, col = 'gray97')
    circos.text(mean(xlim), mean(ylim),chr, cex = 1, col = "black",
                facing = "inside", niceFacing = TRUE)
}, track.height = 0.1, bg.border = NA)

circos.genomicLink(lwd = 2,region1 =ETS.fusion.df[,2:4],region2 = ETS.fusion.df[,5:7],border = ifelse(ETS.fusion.df$value%in%1, '#7570B3',
                                                                                                      ifelse(ETS.fusion.df$value%in%2, '#E7298A',
                                                                                                             ifelse(ETS.fusion.df$value%in%3, '#1B9E77',
                                                                                                                    ifelse(ETS.fusion.df$value%in%4, '#66A61E',
                                                                                                                           ifelse(ETS.fusion.df$value%in%5,'#D95F02',
                                                                                                                                  'black'))))))
########

#### Supplementary Figure 5B

#######
TPM.WCDT <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes.txt",sep="\t",header = T,data.table = F)
TPM.PROMOTE <- fread(file = "./data/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt",sep="\t",header = T,data.table = F)

RNA.TPM.df <-  right_join(TPM.WCDT,TPM.PROMOTE)
RNA.TPM.df <- RNA.TPM.df[!duplicated(RNA.TPM.df$gene_name),]
rownames(RNA.TPM.df) <- RNA.TPM.df$gene_name;RNA.TPM.df<- RNA.TPM.df[,-c(1,which(colnames(RNA.TPM.df)=='gene_name'))]

RNA.data <- RNA.TPM.df[c('ETV4','CTB-175E5.4'),]
RNA.data <- na.omit(RNA.data)
RNA.data <- log2(RNA.data+1)
CTB15E5fusion.bp <- as.data.frame(t(RNA.data[,c(Nelson.samples.PCA[Nelson.samples.PCA$Subtypes%in%c('AR-/NE-'),]$samples)]))

CTB15E5fusion.bp$ID <- rownames(CTB15E5fusion.bp)
CTB15E5fusion.bp <- merge(CTB15E5fusion.bp,Nelson.samples.PCA[,1:2],by.x='ID',by.y='samples',all.x=T)
CTB15E5fusion.bp$color <- factor(CTB15E5fusion.bp$ID %in% "DTB-191-BL",levels=c('TRUE','FALSE'))
melt.CTB15E5fusion.bp <- melt(CTB15E5fusion.bp)


SFigure5B_barplot <- ggbarplot(melt.CTB15E5fusion.bp,x = 'ID',y = 'value',fill = 'color',palette = c('red','grey'),legend.title='Fusion detected',
                                        ylab = 'log2(TPM+1)',xlab = '',ggtheme = theme_bw(base_size = 18,base_rect_size = 1.5))+
    facet_wrap(~variable,scales = 'free_y')+
    theme(axis.text.x = element_text(angle = 90),legend.position = 'bottom',axis.text.y = element_text(size=18),axis.title = element_text(size=18),
          strip.background =element_rect(fill='white',color='white'),strip.text = element_text(face = 'bold',size = 15),
          panel.grid = element_blank())
#######
