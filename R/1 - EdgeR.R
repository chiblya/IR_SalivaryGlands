
## This script was used for the analysis of human biopsies of irradiated and control salivary
## glands in Manuscript # _____________ doi: 

#title: "EdgeR Human SG RNAseq"
#author: "Alex Chibly"


#Analysis libraries
library(edgeR)
library(limma)
library(RColorBrewer)
library(DESeq2)
library(gplots)
library(dplyr)
library(EDASeq)
library(ggplot2)
library(RUVSeq)
library(Hmisc)
library(Matrix)
library(superheat)
library(reshape2)
library(pheatmap)
library(viridis)
library(gage)
library(fgsea)
library(corrplot)
library(tidyr)

## Exploratory data analysis
#We are going to perform an initial exploratory data analysis on the dataset.
#We first load the read counts and do some sanity checks.
#The Annotated.Counts file has geneIDs
#SampleInfo contains the treatment and sample type for all samples

Annotated.Counts <- read.delim("data/Annotated_Counts.txt")
SampleInfo <- read.delim("data/SampleInfo.txt")
summary(Annotated.Counts)
Annotated.Counts <- distinct(Annotated.Counts, external_gene_name, .keep_all = T)
rownames(Annotated.Counts) <- Annotated.Counts[,1]

#Artificial IERCC spikes were used - they can be removed
scrubbedRawCountFile <- Annotated.Counts[rownames(Annotated.Counts)[grep("_|ERCC-", rownames(Annotated.Counts), invert = TRUE)], ]

Annotated.Counts <- scrubbedRawCountFile[,-1]
table(colnames(Annotated.Counts)==SampleInfo$sample)
rownames(SampleInfo) <- SampleInfo$sample


#Based on preliminary PCA analysis, we drop outlier samples 
drop<-c("s16_PAR", "s15_SMG", "s03_PAR", "s03_SMG") ##remove bad samples
Annotated.Counts<-Annotated.Counts[,!(names(Annotated.Counts) %in% drop)]
SampleInfo <- SampleInfo[!(rownames(SampleInfo) %in% drop), ]

print("Number of genes and samples:")
dim(Annotated.Counts)

#We observe the there are **56528** genes quantified, strongly biased due to non-expressing genes (median = 0). Therefore we need to remove these genes from subsequent steps, as they would greatly affect the creation of dispersion models and normalization steps. Because very low expression genes  increase noise and affect the downstream model fitting steps, we set  a conservative filtering threshold: **keep genes with at least 10 reads in at least 5 samples**.
minreads = 10
minsamples = 5
filter <- apply(Annotated.Counts, 1, function(x) length(x[x>minreads])>=minsamples)
counts <- Annotated.Counts[filter,]
geneCounts <- rownames(counts)[grep("^ERCC-", rownames(counts), invert = TRUE)]
erccCounts <- rownames(counts)[grep("^ERCC-", rownames(counts))]
dim(counts)
Annotated.Counts <- counts
head(Annotated.Counts)

#We have **discarded 36623** genes with very low or no expression (and therefore non informative and very noisy) and **kept 19905** genes,  bumping the median from 0 to about 100 on average. 
#Next, we create variables for condition, sample type, and treatment which will be used for downstream statistical analysis:

SampleInfo<-SampleInfo[!(rownames(SampleInfo) %in% drop),] 
SampleInfo
##We create a EDASeq Expression set object to store the counts and the phenotypic data.
condition <- as.factor(SampleInfo$condition)
Treatment <- as.factor(SampleInfo$treatment)
Gland <- as.factor(SampleInfo$gland)

#Then, we create the expression matrix and perform initial PCA analysis on raw (non-normalized) counts
set <- newSeqExpressionSet(as.matrix(Annotated.Counts), phenoData = SampleInfo)
colors <-brewer.pal(4, "Set2")
par(las=2)
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[condition]) 
plotPCA(set, col=colors[condition], cex=1.2)
dev.off()
#Some of the basic descriptive statistics are library-size dependent.

library_size <-as.data.frame(colSums(counts))
library_size$samples <- factor(rownames(library_size), levels = rownames(library_size))
library_size$million_reads <- as.vector(library_size$`colSums(counts)`*1e-6)
library_size$treatment <- SampleInfo$treatment
library_size$gland <- SampleInfo$gland
library_size <- library_size[,2:5]
ggplot(library_size) + geom_bar(aes(x = samples, y = million_reads, fill = treatment), stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#These differences results in variation across expression profiles that need to be normalized before comparison:
counts.melt <- melt(counts)
ggplot(counts.melt) + geom_density(aes(x = value, colour = variable)) + labs(x = "read counts") + theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"),legend.position='right', legend.text=element_text(size=14), legend.title=element_text(size=14)) + scale_x_log10() + ggtitle("Gene expression profiles before normalization") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(plot.subtitle = element_text(vjust = 1), 
                                                                       plot.caption = element_text(vjust = 1), 
                                                                       plot.title = element_text(hjust = 0.5), 
                                                                       legend.text = element_text(size = 12), 
                                                                       legend.key = element_rect(fill = NA), 
                                                                       legend.background = element_rect(fill = NA), 
                                                                       legend.position = c(0.80, 0.5)) + theme(legend.position = c(0.85, 0.5))


#This can also be represented as the more familiar *box-whisker* plots and principal component analysis (PCA) plot to better understand sample relations:
colors <- brewer.pal(nlevels(Treatment), "Set2")
par(mai=c(1.5,1,1,1))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[Treatment], las=2)
par(mai=c(1,1,1,1))
print("PCA colored by treatment (CTRL vs IR):")
plotPCA(set, col=colors[Treatment], cex=1.2)
print("PCA colored by type of gland (PG vs SMG")
plotPCA(set, col=colors[Gland], cex=1.2)
dev.off()
#The PCA analysis reveals that samples tend to cluster based on treatment (PC1 separates the irradiated from the control) while there appears some variance driven by gland type (SMG over PAR) on PC2.

##Normalization
#  We are going to explore two normalization procedures to see which one produces the best results to identify differential expression (DE), and therefore to better separate the different covariates in this experiment. 
# We are going to use the facilities provided by the package RUVSeq. The first step is to perform a between sample normalization, based on applying a factor derived from each sample library size and other distributional differences.
#We will use the `upper` method, a scaling normalization that forces the median of the upper quartile of each sample to be the same. This is described in detail here:
#_J. H. Bullard, E. A. Purdom, K. D. Hansen and S. Dudoit (2010). Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments. BMC Bioinformatics Vol. 11, Article 94._

set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[SampleInfo$condition], las = 2)
plotPCA(set, col=colors[SampleInfo$condition], cex=1.2)

#This baseline normalization improves significantly the variances and median
#gene expression accross samples (therefore enabling direct comparison) and
#improves the resolution of covariates based on PCA analysis. Now this
#**"upper"** normalization is standard and can be combined with additional steps
#to further refine the normalization across samples. We will explore a
#normalization method based on genes that are not affected by the covariates. We
#can define these genes based on an **"in silico"**  method. In this case, we
#remove the top 5000 most variable genes and use the rest to normalize. We first
#run a preliminary DE pass on the standard (upper quartile) normalized data and
#then use that information to compute the empirical **in silico** negative
#controls and associated `W_1`normalization factors.

design<- model.matrix(~condition, data=pData(set))
design
y <- DGEList(counts=counts(set), group =pData(set)$condition)
y <- calcNormFactors(y, method = "RLE")
y <-estimateGLMCommonDisp(y, design)
y <-estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
res <- residuals(fit, type="deviance")

#We perform normalization based on empirical controls and W_1 factors
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]
set <-RUVg(set, empirical, k=1)
pData(set)
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[SampleInfo$treatment], las=2)
plotPCA(set, col=colors[SampleInfo$condition], cex=1.2)

#In this particular set, the standard normalization method and normalization
#with residuals produce tight distributions and good covariate resolution in
#the PCA plot. 
# the computation of empirical controls produces tighter distributions.


#Using generalized linear models we will try to extract DE genes across the
#conditions of interest, which in this case are genes that are regulated by IRR
#treatment, and the different gene expression profiles between gland types. We
#can run an initial estimation of DE and see how well these three normalization
#methods perform:

# BtwLane normalization using upper quartile:
design<- model.matrix(~0+ condition, data=pData(set))
y <- DGEList(counts=counts(set), group =pData(set)$condition)
y <- calcNormFactors(y, method = "RLE")
y <-estimateGLMCommonDisp(y, design)
y <-estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
colnames(fit)
lrt <- glmLRT(fit, contrast = c(-1,0,1,0))
topTags(lrt, n = 20) ##List of top 20 DE genes when comparing CTRL PG vs CTRL SMG.

#Normalization with empirical negative controls:
design<- model.matrix(~0+condition+W_1, data=pData(set))
y <- DGEList(counts=counts(set), group=pData(set)$condition)
y <- calcNormFactors(y, method = "RLE")
y <-estimateGLMCommonDisp(y, design)
y <-estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)

lrt <- glmLRT(fit, contrast = c(-1,0,1,0,0))
topTags(lrt, n = 20) ##List of top 20 DE genes when comparing PAR to CTRL vs IR.

#While they both perform reasonably well, empirical normalization produces 
#higher p values (less significance), resulting in a more stringent analysys. 
#Between lane normalization produces a similar gene list with extremely low
# p values. 
# we will use the empirical normalization factors to perform the DE analysis using
# the EdgeR method.


## preEdgeR analysis
### Recalculate y with desired design and normalization
design<- model.matrix(~0+condition +W_1, data=pData(set))
y <- DGEList(counts=counts(set), group =pData(set)$condition)
y <- calcNormFactors(y, method = "RLE")
y <-estimateGLMCommonDisp(y, design)
y <-estimateGLMTagwiseDisp(y, design)

## transformation
rawCounts <- y$counts
ylog2 <- cpm(y,log=TRUE,normalized.lib.sizes=TRUE,prior.count=2)
ylog2data= as.data.frame(ylog2) # prior count like avelogcpm
ylog2data$geneID <- rownames(ylog2data)
ndata <- cpm(y,log=FALSE,normalized.lib.sizes=TRUE)*1e6
ndataf = as.data.frame(ndata)
ndataf$geneID <- rownames(ndataf)

## saving data tables of normalized counts
write.table(rawCounts, file="output/GlandsCombined_edgeR_RawCounts.txt",sep="\t",col.names=NA, quote = FALSE)
write.table(ylog2data,file="output/GlandsCombined_edgeR_normalized_log2_CPMs.txt",sep="\t",col.names=NA, quote = FALSE)
write.table(ndataf,file="output/GlandsCombined_edgeR_normalized_CPMs.txt",sep="\t",col.names=NA, quote = FALSE)

#Printing normalized gene expression profiles
print("Below is a representation of the gene expression profiles after EdgeR normalization in counts per million (CPM)")
df.m <- melt(as.data.frame(ndata))
print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),legend.position='right', legend.text=element_text(size=8), legend.title=element_text(size=8)) + scale_x_log10())

## clustering / heatmap
hmcol <- colorRampPalette(brewer.pal(9, "Reds"))(50)
distylog2=dist(t(ylog2))
mat = as.matrix(distylog2)
print("Below is a representation of the sample correlations and distances (unsupervised clustering).")
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(6, 6),cexRow = 0.7, cexCol = 0.7)

## MDS plots 
tx.colors <- c("firebrick2", "palegreen2", "goldenrod2", "dodgerblue3")[SampleInfo$condition]
print("Let's look at some PCA representations...")
plotMDS(y, method = "bcv", col="black" ,pch = 21, cex = 2, bg= tx.colors, lwd=2,cex.lab=1.2, cex.axis=1.2)
plotMDS(y, method = "bcv", col=tx.colors, lwd=2,cex.lab=1.2, cex.axis=1.2)

print("BCV PCA calculates distances based on biological coefficient of variation. A set of top genes are chosen that have largest biological variation between the libraries (those with largest genewise dispersion treating all libraries as one group). Then the distance between each pair of libraries (columns) is the biological coefficient of variation (square root of the common dispersion) between those two libraries alone, using the top genes.
The number of genes (top) chosen for this exercise should roughly correspond to the number of differentially expressed genes with materially large fold-changes (550 by default)")

## plotMDS(y) default
plotMDS(y, col=tx.colors,method="logFC" , main="MDS plot logFC", cex.lab=1.2, cex.axis=1.2) 
plotMDS(y, method = "logFC", col="black" ,pch = 21, cex = 2, bg= tx.colors, lwd=2, cex.lab=1.2, cex.axis=1.2)
print("LogFC PCA analysis representing distances between samples based on log2 fold changes")
plotBCV(y, main="BCV plot", cex=0.45, cex.lab=1.5, cex.axis=1.5)

## PCA analysis
pr2=prcomp(t(ylog2))
print("PCA plot using log2 transformed CPM data")
graphX <- pr2$x[,1]
graphY <- pr2$x[,2]
graphhLabels <- colnames(ylog2)
pcadf <- data.frame(graphX,graphY,graphhLabels)

library(Cairo)
cond.colors <- c("firebrick2", "palegreen2", "goldenrod2", "dodgerblue3")[SampleInfo$condition]

pca.plot <- ggplot(data=pcadf, aes(x= graphX, y= graphY))  + 
  theme(axis.title = element_text(face = "bold", size = 14),
        title = element_text(size=16), 
        axis.text = element_text(size=14, color = "black"), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.5,linetype = "solid") ,
        plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0.5, margin = unit(c(0, 0, 0.5, 0), "cm"))) +
  labs(title = "PCA (Log2 transformed CPM)", x = "PC1", y = "PC2") + 
  scale_shape(solid = FALSE) + theme(plot.margin = unit(c(1,1,1,1), "cm"))

svg(filename = "output/PCA colored.svg", width = 5,height = 5)
pca.plot+ geom_point(shape = 21, fill = cond.colors, color = "black", size = 5)
dev.off()

svg(filename = "output/PCA empty.svg", width = 5,height = 5)
pca.plot+ geom_point(shape = 21, fill = "white", color = "black", size = 5)
dev.off()

## EdgeR Analysis 
## Create model matrix with chosen method (standard EdgeR normalization)
design<- model.matrix(~0+condition+W_1, data=pData(set))
y <- DGEList(counts=counts(set), group =pData(set)$condition)
y <- calcNormFactors(y, method = "RLE")
y <-estimateGLMCommonDisp(y, design)
y <-estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
colnames(fit)

degCTRL <- glmLRT(fit, contrast = c(-1,0,1,0,0))
n=dim(y$counts)[1]
ttCTRL=topTags(degCTRL, n=n)

## Parotid vs CTRL controls
#Below are the top 20 DE genes (by logFC) between Parotid and Submandibular glands (controls) and their normalized counts at p<=0.01 . The annotated full list will be saved as BtwLane_EdgeR_PAR_vs_CTRL_DE_genes.txt

topTags(degCTRL, n=20,adjust.method="BH", sort.by="logFC", p.value=0.01)
topCTRL<-cpm(ylog2)[rownames(topTags(degCTRL, n=50, adjust.method="BH", sort.by="logFC", p.value=0.01)),]
head(topCTRL, n=20)
matCTRL <-topCTRL[,c(1:7, 14:19)]
rownames(matCTRL)<-sapply(strsplit(rownames(matCTRL), ".", fixed = TRUE), "[", 1)
dfCTRL<- as.data.frame(pData(set)[c(c(1:7, 14:19)),c("gland","treatment")])
par(mai=c(1.5,1.5,1.5,1.5))
pheatmap(matCTRL, annotation_col=dfCTRL,fontsize=10, scale = "row", color = magma(6), breaks = c(-2:3))
resCTRL = as.data.frame(ttCTRL)
resCTRL.sig <- resCTRL[resCTRL$PValue<0.01,]
dev.off()
## Number of DEGS with p-val <0.01
print(dim(resCTRL.sig))
write.table(resCTRL,file="output/Empirical_EdgeR_PGCTRL_vs_SMGCTRL_DE_genes.txt",sep="\t",col.names=NA, quote = FALSE)

## Control vs IR (both glands combined)
#Below are the top 20 DE genes (by logFC) between Control and Irradiated glands and their normalized counts at p<=0.01 . The annotated full list will be saved as EdgeR_CONTROL_vs_IRRADIATED_DE_genes.txt

colnames(fit)
degIRR=glmLRT(fit, contrast = c(-1/2,1/2,-1/2,1/2,0))
n=dim(y$counts)[1]
ttIRR=topTags(degIRR, n=n)
topTags(degIRR, n=20,adjust.method="BH", sort.by="logFC", p.value=0.01)
topIRR<-cpm(ylog2)[rownames(topTags(degIRR, n=50, adjust.method="BH", sort.by="logFC", p.value=0.01)),]
head(topIRR, n=20)
matIRR <-topIRR
rownames(matIRR)<-sapply(strsplit(rownames(matIRR), ".", fixed = TRUE), "[", 1)
dfIRR<- as.data.frame(pData(set)[,c("gland","treatment")])
par(mai=c(1.5,1.5,1.5,1.5))
pheatmap(matIRR, annotation_col=dfIRR,fontsize=9,scale = "row", color = magma(6), breaks=c(-2:3)) 
resIRR = as.data.frame(ttIRR)
resIRR.sig <- resIRR[resIRR$PValue<0.01, ]
resIRR.sig$gene <- rownames(resIRR.sig)
dev.off()
## Number of DEGS with p-val <0.01
print(dim(resIRR.sig))
write.table(resIRR,file="output/Empirical_EdgeR_CONTROL_vs_IRRADIATED_DE_genes.txt",sep="\t",col.names=NA, quote = FALSE)

########## SMG analysis ##############
colnames(fit)
degSMG <- glmLRT(fit, contrast = c(0,0,-1,1,0))
n=dim(y$counts)[1]
ttSMG=topTags(degSMG, n=n)

topTags(degSMG, n=20,adjust.method="BH", sort.by="logFC", p.value=0.01)
topSMG<-cpm(ylog2)[rownames(topTags(degSMG, n=50, adjust.method="BH", sort.by="logFC", p.value=0.01)),]
head(topSMG, n=20)

matSMG <-topSMG[,c(1:13)]
rownames(matSMG)<-sapply(strsplit(rownames(matSMG), ".", fixed = TRUE), "[", 1)
dfSMG<- as.data.frame(pData(set)[c(1:13),c("gland","treatment")])
par(mai=c(1.5,1.5,1.5,1.5))
pheatmap(matSMG, annotation_col=dfSMG,fontsize=10, scale = "row", color = magma(6), breaks = c(-2:3))
resSMG = as.data.frame(ttSMG)
resSMG$gene <- rownames(resSMG)
write.table(resSMG,file="output/Empirical_EdgeR_SMG Ctrl vs IR_DE_genes.txt",sep="\t",col.names=NA, quote = FALSE)

SMG.DEGs <- merge(resSMG,ylog2data, by.x="gene", by.y="geneID")
SMG.DEGs <- SMG.DEGs[,1:19]
counts$gene <- rownames(counts)
SMG.DEGs.raw <- merge(resSMG,counts, by.x="gene", by.y="gene")
SMG.DEGs.raw <-SMG.DEGs.raw[,1:19]

## Number of DEGS with p-val <0.01
resSMG.sig <- resSMG[resSMG$PValue<0.05,]
resSMG.sig <- rbind(resSMG.sig[resSMG.sig$logFC< -1, ],resSMG.sig[resSMG.sig$logFC>1,] )
resSMG.sig$gene <- rownames(resSMG.sig)

#####PG ANALYSIS #####
colnames(fit)
degPG <- glmLRT(fit, contrast = c(-1,1,0,0,0))
n=dim(y$counts)[1]
ttPG=topTags(degPG, n=n)

topTags(degPG, n=20,adjust.method="BH", sort.by="logFC", p.value=0.01)
topPG<-cpm(ylog2)[rownames(topTags(degPG, n=50, adjust.method="BH", sort.by="logFC", p.value=0.01)),]
head(topPG, n=20)
matPG <-topPG[,c(14:24)]
rownames(matPG)<-sapply(strsplit(rownames(matPG), ".", fixed = TRUE), "[", 1)
dfPG<- as.data.frame(pData(set)[c(14:24),c("gland","treatment")])
par(mai=c(1.5,1.5,1.5,1.5))
pheatmap(matPG, annotation_col=dfPG,fontsize=10, scale = "row", color = magma(6), breaks = c(-2:3))
resPG = as.data.frame(ttPG)
resPG$gene <- rownames(resPG)
write.table(resPG,file="output/Empirical_EdgeR_PG Ctrl vs IR_DE_genes.txt",sep="\t",col.names=NA, quote = FALSE)

PG.DEGs <- merge(resPG,ylog2data, by.x="gene", by.y="geneID")
PG.DEGs <- PG.DEGs[,c(1:6, 20:30)]
PG.DEGs.raw <- merge(resPG,counts, by.x="gene", by.y="gene")
PG.DEGs.raw <-PG.DEGs.raw[,c(1:6, 20:30)]

## Number of DEGS with p-val <0.01
resPG.sig <- resPG[resPG$PValue<0.05,]
resPG.sig <- rbind(resPG.sig[resPG.sig$logFC< -1, ],resPG.sig[resPG.sig$logFC>1,] )
resPG.sig$gene <- rownames(resPG.sig)
dev.off()


genes.oi <- c("AQP5", "BHLHA15", "LPO", "OPRPN", "CAR6", "AMY1A", "AMY1B",
                  "XBP1", "CDH1", "BPIFA2","DCPP1", "MUC5B", "MUC7", "KRT5", "KRT14", 
                  "KRT19", "KRT8", "KRT7", "TUBB3", "VIP", "TH", "GFRA2", "NRTN", 
                  "CHRM1", "CHRM3", "RET", "GFRA1", "GFRA3", "NGFR", "ACHE", "CHAT", 
                  "VACHT", "SMGC", "MUC19", "PRB1", "PRB2", "PRB4", 
                  "HTN3", "VIM", "COL1A1", "COL3A1")
genes.oi.PGdf <- subset(resPG, rownames(resPG) %in% genes.oi)
genes.oi.PGdf$gene <- rownames(genes.oi.PGdf)
genes.oi.SMGdf <- subset(resSMG, rownames(resSMG) %in% genes.oi)
genes.oi.SMGdf$gene <- rownames(genes.oi.SMGdf)
genes.oi.both <- merge(genes.oi.PGdf, genes.oi.SMGdf, by="gene")

write.csv(genes.oi.both, file = "output/selected_genes.csv")

 ### filter out non-coding RNAs for heatmaps
rownames(SMG.DEGs) <- SMG.DEGs$gene
scrubbed.SMG.DEGs <-  SMG.DEGs[rownames(SMG.DEGs)[grep("AC0", rownames(SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("AC", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("AL0", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("AL1", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("AL2", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("AL3", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("AL4", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("AL5", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("AL6", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("AL7", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("AL8", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("AP0", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("ENSG", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("LINC", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("RN7", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("SNOR", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("SCARN", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("RNU", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("-AS1", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("RASL", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <-  scrubbed.SMG.DEGs[rownames(scrubbed.SMG.DEGs)[grep("RPL", rownames(scrubbed.SMG.DEGs), invert = TRUE) ], ]
scrubbed.SMG.DEGs <- scrubbed.SMG.DEGs[scrubbed.SMG.DEGs$PValue<0.05,]
scrubbed.SMG.DEGs <- rbind(scrubbed.SMG.DEGs[scrubbed.SMG.DEGs$logFC< -1, ],scrubbed.SMG.DEGs[scrubbed.SMG.DEGs$logFC>1,] )
scrubbed.SMG.DEGs <- scrubbed.SMG.DEGs[order(scrubbed.SMG.DEGs$logFC, decreasing = T),]
topSMG <- rbind(head(scrubbed.SMG.DEGs, 25), tail(scrubbed.SMG.DEGs, 25))
matSMG <- topSMG[,c(7:19)]
pdf(file = "output/Top50 SMG DEGs heatmap (no pseudogenes).pdf",width = 5.5, height = 9, useDingbats = F)
pheatmap(matSMG, annotation_col=dfSMG,fontsize=10, 
         scale = "row", color = magma(6), breaks = c(-2:3), 
         cluster_cols = F, cluster_rows = F, border_color = "black")
pheatmap(matSMG, annotation_col=dfSMG,fontsize=10, 
         scale = "row", color = magma(6), breaks = c(-2:3), 
         cluster_cols = F, cluster_rows = F, border_color = NA)

dev.off()

# 
rownames(PG.DEGs) <- PG.DEGs$gene
scrubbed.PG.DEGs <-  PG.DEGs[rownames(PG.DEGs)[grep("AC0", rownames(PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("AC", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("AL0", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("AL1", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("AL2", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("AL3", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("AL4", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("AL5", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("AL6", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("AL7", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("AL8", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("AP0", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("ENSG", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("LINC", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("RN7", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("AJ0", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("SNOR", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("-AS1", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]
scrubbed.PG.DEGs <-  scrubbed.PG.DEGs[rownames(scrubbed.PG.DEGs)[grep("RNU", rownames(scrubbed.PG.DEGs), invert = TRUE) ], ]

scrubbed.PG.DEGs <- scrubbed.PG.DEGs[scrubbed.PG.DEGs$PValue<0.05,]
scrubbed.PG.DEGs <- rbind(scrubbed.PG.DEGs[scrubbed.PG.DEGs$logFC< -1, ],scrubbed.PG.DEGs[scrubbed.PG.DEGs$logFC>1,] )
scrubbed.PG.DEGs <- scrubbed.PG.DEGs[order(scrubbed.PG.DEGs$logFC, decreasing = T),]
topPG <- rbind(head(scrubbed.PG.DEGs, 25), tail(scrubbed.PG.DEGs, 25))
matPG <- topPG[,c(7:17)]

pdf(file = "output/Top50 PAR DEGs heatmap (no pseudogenes).pdf",width = 5.5, height = 9, useDingbats = F)
pheatmap(matPG, annotation_col=dfPG,fontsize=10, 
         scale = "row", color = magma(6), breaks = c(-2:3), 
         cluster_cols = F, cluster_rows = F, border_color = "black")
pheatmap(matPG, annotation_col=dfPG,fontsize=10, 
         scale = "row", color = magma(6), breaks = c(-2:3), 
         cluster_cols = F, cluster_rows = F, border_color = NA)
dev.off()

### Upload results from pathway analysis
SMG.pathway.analysis <- read.delim("data/SMG pathway analysis.txt") #file contains 443observations and 5 columns
PG.pathway.analysis <- read.delim("data//PG pathway analysis.txt") #file contains 499observations and 5 columns

### Filter out non-significant pathways with a log p valua >1.3 (p<0.05)
SMG.pathway.analysis <- SMG.pathway.analysis[SMG.pathway.analysis$X.log.p.value.>1.3,] # 64 significant pathways
PG.pathway.analysis <- PG.pathway.analysis[PG.pathway.analysis$X.log.p.value.>1.3,] #101 significant pathways
names(SMG.pathway.analysis)[2] <- "pVal"
names(PG.pathway.analysis)[2] <- "pVal"

SMG.pathway.analysis <- SMG.pathway.analysis[SMG.pathway.analysis$Ingenuity.Canonical.Pathways %in% PG.pathway.analysis$Ingenuity.Canonical.Pathways, ]
SMG.pathway.analysis <- merge(SMG.pathway.analysis, PG.pathway.analysis, by="Ingenuity.Canonical.Pathways")
SMG.pathway.analysis$score <- SMG.pathway.analysis$pVal.x + SMG.pathway.analysis$pVal.y
SMG.pathway.analysis <- SMG.pathway.analysis[order(SMG.pathway.analysis$score, decreasing = T),]
write.table(SMG.pathway.analysis, file = "output/Common top pathways.txt", sep = "\t",row.names = F)

#### Make plot with overlapping genes in different pathways
library(UpSetR)


ECM.smg.genes <- unique(c(unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.x[1]), ",")), 
                          unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.x[2]), ",")),
                          unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.x[7]), ","))))
nerve.smg.genes <- unique(c(unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.x[3]), ",")), 
                            unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.x[5]), ","))))
inf.smg.genes <- unique(c(unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.x[4]), ",")), 
                          unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.x[8]), ",")),
                          unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.x[10]), ","))))

ECM.pg.genes <- unique(c(unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.y[1]), ",")), 
                          unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.y[2]), ",")),
                          unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.y[7]), ","))))
nerve.pg.genes <- unique(c(unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.y[3]), ",")), 
                            unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.y[5]), ","))))
inf.pg.genes <- unique(c(unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.y[4]), ",")), 
                          unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.y[8]), ",")),
                          unlist(strsplit(as.character(SMG.pathway.analysis$Molecules.y[10]), ","))))


mypathways <- as.character(c("Neurotrophin signaling", "Fibrosis", "Inflammation"))
common.genes <-unique(c(nerve.smg.genes[nerve.smg.genes %in% nerve.pg.genes],
                        ECM.smg.genes[ECM.smg.genes %in% ECM.pg.genes], 
                 inf.smg.genes[inf.smg.genes %in% inf.pg.genes]))
common.genes[70] <- "WISP1" #Fix name for consistency

test.data <- matrix(nrow = length(common.genes), ncol = length(mypathways))
rownames(test.data) <- common.genes
colnames(test.data) <- mypathways

test.data[]<-1
setlist <- list(nerve.smg.genes,ECM.smg.genes,inf.smg.genes)

for (j in 1:length(setlist)) {
  for (i in 1:length(common.genes)) {
    if (rownames(test.data)[i] %in% setlist[[j]] ==FALSE){
      test.data[i,j] = 0
    }
  }
}

test.data.df <- as.data.frame(test.data)
test.data.melted <- melt(test.data)

library(circlize)


test.data.melted <- test.data.melted[,c(2,1,3)]
border_mat2 = matrix("black", ncol = 3, nrow = nrow(test.data.melted))
colnames(border_mat2) = c("Neurotrophin signaling", "Fibrosis", "Inflammation")
rownames(border_mat2) = test.data.melted$Var1
border_mat2[,2:3] = "blue"

test.data.melted$color <- NA
test.data.melted[test.data.melted$Var2 %in% "Neurotrophin signaling" & test.data.melted$value>0, ]$color <- "black" 

pdf("output/chord diagram pathway genes.pdf", useDingbats = F, width = 8,height = 8)
chordDiagram(test.data.melted, link.lwd = 1, link.lty = 1,grid.col = c("Fibrosis" = "#EEDE0150", "Neurotrophin signaling" = "firebrick2", "Inflammation" = "#0148EE50"),
                       symmetric = F, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1,
                       link.arr.length = 0.1, link.arr.type = "big.arrow", link.border = test.data.melted$color,
                       link.largest.ontop = T, link.arr.lty = 1,grid.border = 1, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test.data.melted))))))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()


# generate plot
p <- ggplot(test.data.melted, aes(x =Var1, y = Var2)) 

common.mat <- ylog2data[,1:24]
balloonheatmap <- subset(common.mat[rownames(common.mat) %in% common.genes,])
balloonheatmap <- balloonheatmap[order(match(rownames(balloonheatmap), rownames(test.data))),]

pdf("output/pathway balloon plot.pdf", useDingbats = F, width = 13, height = 2.5)
p+geom_point( aes(size=value, fill=value),shape=21, colour="black")+
  theme(panel.background=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  scale_size_area(max_size=4)+
  #Add labels to axes
  labs(x="Gene", y="Pathway") + theme(axis.text.x = element_text(angle=90, colour = "black", face = "italic", vjust = 0.5,hjust = 1, size = 10), axis.text.y = element_text(color = "black", face = "bold", size = 6))
dev.off()

pdf("output/pathway heatmap plot.pdf", useDingbats = F, width = 4, height = 9.5)
pheatmap(balloonheatmap, annotation_col=dfIRR, 
         fontsize = 7, scale = "row", color = magma(6), 
         breaks = c(-2:3), cluster_rows = FALSE, border_color = "black")
dev.off()

### load upstream regulator files
SMG.Upregs <- read.delim("data/SMG pathway analysis - Upregs (p05).txt")
PG.Upregs <- read.delim("data/PG pathway analysis - Upregs (p05).txt")
test <- SMG.Upregs[SMG.Upregs$Upstream.Regulator %in% PG.Upregs$Upstream.Regulator,]

SMG.upregs.sig <- SMG.Upregs[SMG.Upregs$Upstream.Regulator %in% resSMG.sig$gene, ]
SMG.Upregs.receptors <- SMG.upregs.sig[!(SMG.upregs.sig$Molecule.Type == "transcription regulator"), ]
SMG.Upregs.receptors <- SMG.Upregs.receptors[order(SMG.Upregs.receptors$Expr.Log.Ratio), ]
SMG.Upregs.tfs <- SMG.upregs.sig[(SMG.upregs.sig$Molecule.Type == "transcription regulator"), ]
SMG.Upregs.tfs <- SMG.Upregs.tfs[order(SMG.Upregs.tfs$Expr.Log.Ratio), ]

PG.upregs.sig <- PG.Upregs[PG.Upregs$Upstream.Regulator %in% resPG.sig$gene, ]
PG.Upregs.receptors <- PG.upregs.sig[!(PG.upregs.sig$Molecule.Type == "transcription regulator"), ]
PG.Upregs.tfs <- PG.upregs.sig[(PG.upregs.sig$Molecule.Type == "transcription regulator"), ]
PG.Upregs.receptors <- PG.Upregs.receptors[order(PG.Upregs.receptors$Expr.Log.Ratio), ]
PG.Upregs.tfs <- PG.Upregs.tfs[order(PG.Upregs.tfs$Expr.Log.Ratio), ]

##common upstream regulators
common.upregs <- SMG.upregs.sig[SMG.upregs.sig$Upstream.Regulator %in% PG.upregs.sig$Upstream.Regulator, ]
common.upreg.receptors <- common.upregs[!(common.upregs$Molecule.Type == "transcription regulator"), ]
common.upreg.tfs <- common.upregs[(common.upregs$Molecule.Type == "transcription regulator"), ]
common.upreg.receptors <- common.upreg.receptors[order(common.upreg.receptors$Expr.Log.Ratio), ]
common.upreg.tfs <- common.upreg.tfs[order(common.upreg.tfs$Expr.Log.Ratio), ]

pdf(file = "output/Upstream regulator heatmap.pdf", width = 6,height = 5, useDingbats = F)
pheatmap(subset(common.mat[rownames(common.mat) %in% common.upregs$Upstream.Regulator,]), 
         annotation_col=dfIRR, fontsize = 10, scale = "row", color = magma(6), 
         breaks = c(-2:3), cluster_cols = T, cluster_rows = T, border_color = "black")
dev.off()

######## GENERATE BOX PLOTS WITH NTRK GENES
library(tidyverse)
library(hrbrthemes)

expression.matrix <- ylog2data[,1:24]
dataTEST <- subset(expression.matrix, rownames(expression.matrix) %in% c("NGF", "NTF4", "NTF3", "BDNF"))
dataTEST <- dataTEST[order(rownames(dataTEST)),]
dataTEST <- as.data.frame(t(dataTEST))
dataTEST$gland <- rownames(dataTEST)

mdataTEST <- melt(dataTEST, id=c("gland"))
dataTEST <- merge(mdataTEST, SampleInfo, by.x  = "gland", by.y = "sample")

ggplot(dataTEST, aes(x=variable, y=value, fill=condition)) + 
  geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.6) + theme(axis.text = element_text(size = 12, colour = "black"), axis.line = element_line(size = 0.5), panel.background = element_blank()) + xlab("") +ylab("Scaled expression value") + facet_wrap(~variable, scale="free", ncol = 4)

####### EXPORT FINAL TABLES ########
write.csv(PG.DEGs, "output/Human PG EdgeR analysis - DEG with log2 cpms.csv")
write.csv(SMG.DEGs, "output/Human SMG EdgeR analysis - DEG with log2 cpms.csv")
write.csv(resCTRL, "output/Human SMG vs PG control EdgeR analysis.csv")
write.csv(PG.DEGs.raw, "output/Human PG EdgeR analysis - DEG with raw counts.csv")
write.csv(SMG.DEGs.raw, "output/Human SMG EdgeR analysis - DEG with raw counts.csv")
write.csv(SMG.upregs.sig, "output/Human SMG Upstream Regulators.csv")
write.csv(PG.upregs.sig, "output/Human PG Upstream Regulators.csv")
write.csv(common.upregs, "output/Human Common Upstream Regulators.csv")
#write.csv(common.upregs2, "output/Human Common Upstream Regulators(PG).csv")
write.table(expression.matrix, "output/Expression matrix.txt", sep = "/t")
write.csv(expression.matrix, "output/Expression matrix.csv")

save(resPG.sig,resSMG.sig,resCTRL,SMG.upregs.sig,PG.upregs.sig,SampleInfo,expression.matrix, file = "output/Final Human gene lists.RData")

#### END OF SCRIPT ####


