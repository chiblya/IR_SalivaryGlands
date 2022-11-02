rm(list=ls())
source("~/CustomFunctions/pubPlots.R")
source("~/CustomFunctions/UtilFunctions.R")

library(SummarizedExperiment)
library(MultiAssayExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(Seurat)
library(SeuratDisk)
library(pheatmap)

# :::::::::::::::::::::::::::::::::::::::::::: 
# CODE TO ADDRESS REVIEWERES COMMENTS
# :::::::::::::::::::::::::::::::::::::::::::: 

cpm <- read.delim("../../output/GlandsCombined_edgeR_normalized_CPMs.txt", row.names = 1, sep = "\t")
logcpm <- read.delim("../../output/GlandsCombined_edgeR_normalized_log2_CPMs.txt", row.names = 1, sep = "\t")
counts= read.delim("../../output/GlandsCombined_edgeR_RawCounts.txt", row.names = 1, sep = "\t")
SampleInfo <- read.delim("../../data/SampleInfo.txt")
clinData = read.csv("../../data/clinData.csv")

logcpm = logcpm[,-25]
cpm = cpm[,-25]

cData = merge(SampleInfo, clinData, by.x = "sample", by.y = "Sample")
cData$gland = factor(cData$gland)
cData$condition = factor(cData$condition)
cData$treatment = factor(cData$treatment)

cData$TimeCat = factor(case_when(cData$Time <12 ~ "Acute", 
                                 cData$Time >= 12 ~ "Chronic", 
                                 is.na(cData$Time) ~ "Naive"), levels = c("Naive", "Acute", "Chronic"))


rownames(cData) = cData$sample
stopifnot(all(colnames(logcpm) %in% cData$sample))
cData = cData[match(colnames(logcpm), cData$sample),]


SE_counts = SummarizedExperiment(assays = list(counts = counts), colData = cData)
SE_logcpm = SummarizedExperiment(assays = list(logcpmed = logcpm), colData = cData)


saveRDS(SE_counts, "Raw_counts_SE.rds")
saveRDS(SE_logcpm, "logCPM_SE.rds")

# ::::::::::::::::::::::::::::::::::::::::::::::::::::: #
#           VISUALIZATION OF SELECTED GENES
# ::::::::::::::::::::::::::::::::::::::::::::::::::::: #

fulldf = cbind(cData, t(logcpm))
rownames(logcpm)[grep(pattern = "NTRK", x = rownames(logcpm))]

genes = c("NGFR", "NTRK1", "NTRK2", "NTRK3", "NGF")

geneplots = lapply(genes, function(g){
  clean_boxPlot(fulldf, xVar = "TimeCat", yVar = g, palette = RT.palette, stats = F) + bigPlot() + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab("") +
    stat_compare_means(ref.group = "Naive", label = "p.signif")
})

ggarrange(plotlist = geneplots, ncol = 5, nrow = 1)


# :::::::::::::::::::::::::::::::::::::::::::: 
# LOAD TABULA SAPIENS DATA
# :::::::::::::::::::::::::::::::::::::::::::: 


Convert("TS_Salivary_Gland.h5ad", dest = "h5seurat", overwrite = TRUE)
humanSC <- LoadH5Seurat("TS_Salivary_Gland.h5seurat", assays = "RNA")
#hfile = Connect(filename = "TS_Salivary_Gland.h5seurat")
#hfile$index()
#hfile[["assays/RNA"]]


humanSC = SplitObject(humanSC, split.by = "donor")
humanSC = lapply(X = humanSC, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = humanSC)
anchors <- FindIntegrationAnchors(object.list = humanSC, anchor.features = features)
humanSC <- IntegrateData(anchorset = anchors)


DefaultAssay(humanSC) <- "integrated"

# Run the standard workflow for visualization and clustering
humanSC <- ScaleData(humanSC, verbose = FALSE)
humanSC <- RunPCA(humanSC, npcs = 30, verbose = FALSE)
humanSC <- RunUMAP(humanSC, reduction = "pca", dims = 1:30)
humanSC <- FindNeighbors(humanSC, reduction = "pca", dims = 1:30)
humanSC <- FindClusters(humanSC, resolution = 0.5)


DimPlot(humanSC, label = T) + NoLegend()
DimPlot(humanSC, group.by = "cell_ontology_class", label = T) + NoLegend()
#humanSC = SetIdent(humanSC, value = "cell_ontology_class")
humanMarkers = FindAllMarkers(humanSC, assay = "RNA", only.pos = T, max.cells.per.ident = 2000, random.seed = 3512)
humanMarkers = humanMarkers[humanMarkers$p_val_adj<0.05, ]
topHmarkers = humanMarkers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

plotMarkers = c("EPCAM", "PRB1", "HTN3", "AQP5", "BHLHA15", "BPIFA2", "SMGC",
                "KRT14", "KRT5", "ACTA2", "CNN1",
                "GSTT1",  "ASCL3", "CFTR", "KLK1", "KRT19", "WFDC18",
                "PECAM1", "TUBB3", "NCAM1", "GFRA3", "KIT",
                "ADGRE1",  "GZMA", "VIM", "COL1A1", "TWIST1", "ALAS2")

matrix.cell.markers <- AverageExpression(humanSC, features = unique(plotMarkers), assays = "RNA")
matrix.cell.markers <- matrix.cell.markers$RNA 
library(pheatmap)
pheatmap(matrix.cell.markers, scale = "row", cluster_cols = T, cluster_rows = F)
write.csv(humanMarkers, file = "unsupervised degs tabulaSapiens.csv")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                         ANNOTATE CELLS 
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

humanSC = RenameIdents(humanSC, 
                       "0" = "Acinar",
                       "1" = "Acinar",
                       "2" = "T-cell",
                       "3" = "B-cell",
                       "4" = "ID",
                       "5" = "Fibroblasts",
                       "6" = "B-cell",
                       "7" = "B-cell",
                       "8" = "Fibroblasts",
                       "9" = "Macrophages",
                       "10" = "NK cell",
                       "11" = "Basal duct",
                       "12" = "Acinar",
                       "13" = "ID",
                       "14" = "Endothelial",
                       "15" = "Ionocytes",
                       "16" = "T-cell",
                       "17" = "Plasma cells",
                       "18" = "Monocytes",
                       "19" = "Myoepithelial",
                       "20" = "Pericytes",
                       "21" = "Neutrophils",
                       "22" = "SD",
                       "23" = "Proliferative",
                       "24" = "Lymphatic")

DimPlot(humanSC, label = T, label.size = 3.5, repel = T) + NoLegend()
humanSC$CellType = Idents(humanSC)
#saveRDS(humanSC, file = "TabulaSapiens_SG_reclustered.rds")
humanSC = readRDS("TabulaSapiens_SG_reclustered.rds")

humanMarkers = FindAllMarkers(humanSC, assay = "RNA", only.pos = T, max.cells.per.ident = 2000, random.seed = 3512)
humanMarkers = humanMarkers[humanMarkers$p_val_adj<0.05, ]
topHmarkers = humanMarkers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
saveRDS(humanMarkers, file = "../Reviews/Humansc_cellMarkers.rds")


matrix.cell.markers <- AverageExpression(humanSC, features = unique(topHmarkers$gene), assays = "RNA")
matrix.cell.markers <- matrix.cell.markers$RNA 
pheatmap(matrix.cell.markers, scale = "row", cluster_cols = F, cluster_rows = F)


DefaultAssay(humanSC) = "RNA"
subsetgenes = c("NGF","BDNF", "NTF3", "NTF4", "NTRK1", "NTRK2", "NTRK3", "NGFR")
pdf("HumanSG dot plot neurotrophin genes.pdf", useDingbats = F, width = 5.5, height = 3.5)
DotPlot(humanSC, features = subsetgenes, dot.min = 0.05, col.min = 0) + theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())
dev.off()

#genes = c("ACTA2", "CNN1", "KRT14", "BHLHA15", "AQP5", "ASCL3")
#DotPlot(humanSC, features = genes, dot.min = 0.05, col.min = 0) + theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())

# ::::::::::: Subset SMG gland ::::::::::::::::: #
smgsubset = subset(humanSC, anatomical_information == "Submandibular")
pgsubset = subset(humanSC, anatomical_information %in% c("Parotid", "parotid"))
DotPlot(pgsubset, features = subsetgenes, dot.min = 0.05, col.min = 0) + theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())


# :::::::::::::::::::::::::::::::::::::::::::::: 
# proportion of cells per gland
# :::::::::::::::::::::::::::::::::::::::::::::: 

metadata = humanSC@meta.data
metadata$gland = metadata$anatomical_information
metadata$gland[metadata$gland == "parotid"] = "Parotid"
propGraph_discVar(data = metadata, yVar = "gland", discVar = "CellType") + 
  scale_x_discrete(label = levels(metadata$CellType)) +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=10)) 


# Load human data
hPG.degs <- read.csv("../../output/Human PG EdgeR analysis - DEG with log2 cpms.csv")
hSMG.degs <- read.csv("../../output/Human SMG EdgeR analysis - DEG with log2 cpms.csv")

#generate expression matrix with all genes for visualization
expression.matrix <- merge(hPG.degs[,c(2, 8:18)], hSMG.degs[,c(2,8:20)], by = "gene")
#expression.matrix$gene <- tolower(expression.matrix$gene)
#expression.matrix$gene <- capitalize(expression.matrix$gene)
rownames(expression.matrix) <- expression.matrix$gene
expression.matrix <- as.matrix(expression.matrix[,-1])

#hPG.degs$gene <- tolower(hPG.degs$gene)
#hPG.degs$gene <- capitalize(hPG.degs$gene)
#hSMG.degs$gene <- tolower(hSMG.degs$gene)
#hSMG.degs$gene <- capitalize(hSMG.degs$gene)


# Filter non-significant genes (p>0.05)
hPG.degs <- hPG.degs[hPG.degs$PValue<0.05, ]
hSMG.degs <- hSMG.degs[hSMG.degs$PValue<0.05, ]
hPG.degs <- hPG.degs[hPG.degs$logFC< -1 | hPG.degs$logFC> 1, ]
hSMG.degs <- hSMG.degs[hSMG.degs$logFC>1 | hSMG.degs$logFC< -1, ]


allDEGs = merge(hPG.degs[, -1], hSMG.degs[,-1], by = "gene", all = T)
FCvalues.10x.SMG <- reshape2::dcast(data = humanMarkers,formula = gene~cluster,fun.aggregate = sum,value.var = "avg_logFC")
SMG.DEGs.specific.10x <- base::merge(allDEGs, FCvalues.10x.SMG, by.x = "gene", by.y = "gene") #We use SMG only as 10x data does not contain PG

# write.csv(SMG.DEGs.specific.10x, file = "Cellular localization of SMG-IR DEGs in scRNAseq.csv", row.names = T)
n = rep(0, length(unique(humanSC$CellType)))

counter<- data.frame("NoGenes" =  n, row.names = names(SMG.DEGs.specific.10x[36:53])) #we use columns 21:43 which contain cell data
for (i in 36:53) {
  counter[i-35,] <- length(which(SMG.DEGs.specific.10x[,i]>0))
}
counter$CellType <- rownames(counter)
counter <- counter[order(counter$NoGenes, decreasing = T),]


ggplot(data=counter, aes(x=reorder(CellType, -NoGenes), y=NoGenes)) +
  geom_bar(stat="identity", color="black", width = 0.85, fill="blue") +
  theme_classic() + theme(axis.text.x = element_text(colour = "black", size = 14, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14, colour = "black"), line = element_line(size = 0.5))


upregs <- read.csv("../../output/Human Common Upstream Regulators.csv")
DotPlot(humanSC, features = upregs$Upstream.Regulator,dot.min = 0.15, col.min = 0) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"), axis.title = element_blank())

genes.to.heatmap <- upregs$Upstream.Regulator #select Upstream regulator genes for visualization
genes.to.heatmap <-c("NTRK2", "NGFR", "FGF1", "FGF10", "FGFR1", "SFRP2", "ETV4", "ETV5","BMP7", "NOTCH1","JAG1" , 
  "JAG2", "SNAI2", "TFAP2C", "TP63",        "ITGB3",   "CX3CL1",  "IGF2" ,      "PGF",   "NOSTRIN")           
              

upregs.avgexpression <- AverageExpression(object = humanSC, assays = "RNA", slot = "data", features = genes.to.heatmap)
upregs.avgexpression <- upregs.avgexpression$RNA
pheatmap(upregs.avgexpression,scale = "row")

DotPlot(humanSC, features = genes.to.heatmap,dot.min = 0.15, col.min = 0) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"), axis.title = element_blank())



# write.csv(SMG.DEGs.specific.10x, file = "Cellular localization of SMG-IR DEGs in scRNAseq.csv", row.names = T)
n = rep(0, length(unique(humanSC$CellType)))
upregs.10x <- SMG.DEGs.specific.10x[SMG.DEGs.specific.10x$gene %in% genes.to.heatmap, ]

counter<- data.frame("NoGenes" =  n, row.names = names(upregs.10x[36:53])) #we use columns 21:43 which contain cell data
for (i in 36:53) {
  counter[i-35,] <- length(which(upregs.10x[,i]>0))
}
counter$CellType <- rownames(counter)
counter <- counter[order(counter$NoGenes, decreasing = T),]


ggplot(data=counter, aes(x=reorder(CellType, -NoGenes), y=NoGenes)) +
  geom_bar(stat="identity", color="black", width = 0.85, fill="blue") +
  theme_classic() + theme(axis.text.x = element_text(colour = "black", size = 14, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14, colour = "black"), line = element_line(size = 0.5))



subsetgenes = c("NGF","BDNF", "NTF3", "NTF4", "NTRK1", "NTRK2", "NTRK3", "NGFR")
ntfAvgExp = AverageExpression(humanSC,  assays = "RNA", slot = "data", features = subsetgenes)$RNA
pheatmap(ntfAvgExp, scale = "row")



# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
#     Cell-specific localization of Upstream Regulators (from IPA)
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

upregs <- read.csv("../../output/Human Common Upstream Regulators.csv")
upregs <- upregs[order(upregs$Expr.Log.Ratio, decreasing = T), ]

## calculate average scaled expression of all genes from integrated dataset
genes.to.heatmap <- upregs$Upstream.Regulator #select Upstream regulator genes for visualization
upregs.avgexpression <- AverageExpression(object = humanSC, assays = "RNA", slot = "data", features = genes.to.heatmap)
upregs.avgexpression <- upregs.avgexpression$RNA

library(pheatmap)
pheatmap(upregs.avgexpression,scale = "row")


##### Ligand-Receptor analysis with Upstream regulators ######
# install LigandReceptor v1 package from GitHub
# devtools::install_github("chiblyaa/LigandReceptor")

# Load table of ligand-receptor pairs from published manuscript; https://doi.org/10.1038/ncomms8866
library(LigandReceptor)
Ligand.Receptor.Pairs <- LigandReceptor::LRdatabase
Ligand.Receptor.Pairs$Pair <- toupper(Ligand.Receptor.Pairs$Pair)
Ligand.Receptor.Pairs$Ligand <- toupper(Ligand.Receptor.Pairs$Ligand)
Ligand.Receptor.Pairs$Receptor <- toupper(Ligand.Receptor.Pairs$Receptor)


# create table of ligand-receptor pairs:
names(humanMarkers)[2] = "avg_logFC"
upregpairstable <- LigandReceptorPairsTable(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = upregs$Upstream.Regulator)
write.csv(upregpairstable, file = "Ligand-Receptor pairs with upstream regulators (HUMAN).csv")
colors = c(my.palette, third.palette)[1:18]
names(colors) = levels(humanMarkers$cluster)


# chord plot for MEC interactions
pdf("ChordPlot Upregs all- human.pdf", useDingbats = F, height = 6,width = 9)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = upregs$Upstream.Regulator, cellcolors = colors)
legend("right",   # location of legend
       legend = names(colors), # categories or elements to render in
       # the legend
       fill = colors, bty = "n", cex = 0.8, xpd = TRUE) 
dev.off()



# chord plot for MEC interactions
pdf("ChordPlot Upregs from MECs - human.pdf", useDingbats = F, height = 6,width = 9)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = upregs$Upstream.Regulator, cellcolors = colors, from = "Myoepithelial")
legend("right",   # location of legend
       legend = names(colors), # categories or elements to render in
       # the legend
       fill = colors, bty = "n", cex = 0.8, xpd = TRUE) 
dev.off()

pdf("ChordPlot Upregs to MECs.pdf", useDingbats = F, height = 6,width = 9)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = upregs$Upstream.Regulator, cellcolors = colors, to = "Myoepithelial")
legend("right",   # location of legend
       legend = names(colors), # categories or elements to render in
       # the legend
       fill = colors, bty = "n", cex = 0.8, xpd = TRUE) 
dev.off()




PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = c("JAG1", "NOTCH1") , cellcolors = colors)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = c("FGF1", "ETV4", "FGF10", "FGFR1"), cellcolors = colors)
#PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Nostrin", cellcolors = colors)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Ngfr", cellcolors = colors)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Fgf10", cellcolors = colors)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Igf2", cellcolors = colors)
#PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Sfrp2", cellcolors = colors)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Itgb3", cellcolors = colors)
#PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Cxc3cl1", cellcolors = colors)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Fgf1", cellcolors = colors)
#PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Etv5", cellcolors = colors)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Ntrk2", cellcolors = colors)
#PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Snai2", cellcolors = colors)
#PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Trp63", cellcolors = colors)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Bmp7", cellcolors = colors)
#PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Etv4", cellcolors = colors)
#PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Jag2", cellcolors = colors)
#PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Tfap2c", cellcolors = colors)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Fgfr1", cellcolors = colors)
PairsPlot(seuratDEGS = humanMarkers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Pgf", cellcolors = colors)

