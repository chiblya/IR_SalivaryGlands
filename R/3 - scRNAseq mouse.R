library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(pheatmap)
library(superheat)
library(Hmisc)
library(reshape2)
library(circlize)

##### Load annotated seurat objects from published manuscript ####
e12smg <- readRDS("../../10X Manuscript Projects/10X revisions/E12 SMG Annotated (SEURAT v3).rds")
e14smg <- readRDS("../../10X Manuscript Projects/10X revisions/E14 SMG Annotated (SEURAT v3).rds")
e16smg <- readRDS("../../10X Manuscript Projects/10X revisions/E16 SMG Annotated (SEURAT v3).rds")
p1smg <- readRDS("../../10X Manuscript Projects/10X revisions/P1 SMG Annotated (SEURAT v3).rds")
p30smg <- readRDS("../../10X Manuscript Projects/10X revisions/P30_Male_and_female combined - annotated (split from integrated).rds")
p120smg <- readRDS("../../10X Manuscript Projects/10X revisions/Adult SMG annotated (split from Integrated).rds")

# P1<- DimPlot(e12smg, pt.size = .5,group.by = "CellType",label = F, repel = T, label.size = 6, cols = colors) +NoLegend() + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
# P2<- DimPlot(e14smg, pt.size = .5,group.by = "CellType",label = F, repel = T, label.size = 6, cols = colors) +NoLegend()+ theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
# P3<- DimPlot(e16smg, pt.size = .5,group.by = "CellType",label = F, repel = T, label.size = 6, cols = colors) +NoLegend()+ theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
# P4<- DimPlot(p1smg, pt.size = .5,group.by = "CellType",label = F, repel = T, label.size = 6, cols = colors)  +NoLegend()+ theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
# P5<- DimPlot(p120smg, pt.size = .5,group.by = "CellType",label = F, repel = T, label.size = 6, cols = colors) +NoLegend()+ theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
# plot.list <- list(P1,P2,P3,P4,P5)
# CombinePlots(plots = plot.list, ncol = 3,legend = "bottomright")

##### Set color theme and consistent order of labels based on cell type #####
colors <- c("#eb7341","#8a308c","#4ab539","#ffb558", "#006f35","#b65ea4","#b8cd42","#f3e3cb", 
            "#adf4d0",   "#565e00","#47b5d5","#ff9e96","#01908c", "#7152f6","#e10060", 
            "#345a78","#9c001b","#eee799","#e74134", "#0074c9",  "#7f4a21", "#c3c1ff","#fbe54b") #9a8e62" - mesenchymal

cell.types <- c("Seromucous acinar", "Ascl3+ duct", "Basal duct", "Serous acinar", "Bpifa2+ Proacinar", "End bud", "Endothelial",
                "Erythroid", "GCT", "Glial cells", "Gfra3+ ID", "Krt19+ duct", "Macrophages", "Mast cells",
                "Mitotic cells", "Myoepithelial", "Nerves", "NK cells", "Gstt1+ ID", "Smgc+ Proacinar",
                "Pericytes", "Striated duct", "Stromal")
names(colors) <- cell.types

stage.colors <- c(
  "#c34d54",
  "#d7cd9e",
  "#019f84",
  "#92d8d4",
  "#0076af",
  "#6b0081"
)
names(stage.colors) <- c("E12", "E14", "E16", "P1", "P30", "Adult")

cell.type.levels <- c("End bud", "Krt19+ duct", "Basal duct", 
                      "Myoepithelial", "Bpifa2+ Proacinar", 
                      "Smgc+ Proacinar", "Mitotic cells",
                      "Serous acinar", "Seromucous acinar", 
                      "Gstt1+ ID", "Gfra3+ ID", "Ascl3+ duct", 
                      "GCT", "Striated duct", "Endothelial", 
                      "Pericytes", "Nerves",  "Glial cells",   
                      "Macrophages", "Mast cells", "NK cells", 
                      "Stromal", "Erythroid")

##### E16 SMG visualization ####
e16smg <- RenameIdents(e16smg, "Mesenchyme" = "Stromal", "Smooth muscle" = "Pericytes")
Idents(e16smg) <- factor(Idents(e16smg), levels = cell.type.levels)
e16smg$CellType <- Idents(e16smg)
pdf("E16 UMAP.pdf", useDingbats = F, width = 6, height = 4.5)
DimPlot(e16smg, cols = colors)
dev.off()

subsetgenes = c("Ngf","Bdnf", "Ntf3", "Ntf5", "Ntrk1", "Ntrk2", "Ntrk3", "Ngfr")
pdf("e16 dot plot neurotrophin genes.pdf", useDingbats = F, width = 5.5, height = 3.5)
DotPlot(e16smg, features = subsetgenes, dot.min = 0.05, col.min = 0) + theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())
dev.off()

######### Integration of datasets #######
# create list of objects for analysis
smg.integrated.list <- list(e16smg, p1smg, p30smg, p120smg) # SEURAT objects are already normalized and scaled

# find anchors and perform integration to correct for batch effect. 
# Integration will allow us to compare MECs from different stages
# to evaluate differential gene expression across development
smg.integrated.anchors <- FindIntegrationAnchors(smg.integrated.list, dims = 1:30)
smg.integrated <- IntegrateData(smg.integrated.anchors, dims = 1:30)

#remove list to save memory
remove(smg.integrated.list)

# Run the standard workflow for visualization and clustering
smg.integrated <- ScaleData(smg.integrated, verbose = FALSE, vars.to.regress = "percent.mt")
smg.integrated <- RunPCA(smg.integrated, npcs = 30)
smg.integrated <- RunUMAP(smg.integrated, reduction = "pca", dims = 1:30)
# An optimal resolution for clustering was chosen based on appropriate
# discrimination between cell types, which were already annotated in the SEURAT object
smg.integrated <- FindNeighbors(smg.integrated, reduction = "pca", dims = 1:30) 
smg.integrated <- FindClusters(smg.integrated, resolution = 0.9)

smg.integrated <- SetIdent(smg.integrated, value = "CellType")
#update cell type labels
smg.integrated <- RenameIdents(smg.integrated, "Smgc+" = "Gstt1+ ID", "Acinar" = "Seromucous acinar",
                               "Bpifa2+" = "Serous acinar", "Intercalated duct" = "Gfra3+ ID", 
                               "Mesenchyme" = "Stromal", "Smooth muscle" = "Pericytes")
DimPlot(smg.integrated, cols = colors)


Idents(smg.integrated) <- factor(Idents(smg.integrated), levels = cell.type.levels)
smg.integrated[["CellType.fixed"]] <- Idents(smg.integrated)

saveRDS(smg.integrated, file = "E16-Adult SMG integrated annotated.rds")
##### UMAP of integrated dataset #####
pdf(file = "UMAP E16-Adult SMG Integrated.pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(smg.integrated, reduction = "umap", label = T, label.size = 5, repel = F, group.by = "seurat_clusters") + NoLegend()
dev.off()
pdf(file = "UMAP E16-Adult SMG Integrated (annotated).pdf", useDingbats = F, width = 9, height = 5)
DimPlot(smg.integrated, reduction = "umap", label = F, group.by = "CellType.fixed", cols = colors)
dev.off()


DefaultAssay(smg.integrated) <- "RNA"
SEROUS <- FindMarkers(smg.integrated, ident.1 = "Seromucous acinar", only.pos = T)

genes.to.heatmap1 <- c("Epcam", "Krt19", "Krt14", "Krt5", "Acta2", "Cnn1", "Aqp5", "Bhlha15", "Bpifa2", "Smgc",
                      "Mki67", "Dcpp1", "Mucl2", "Car6",
                       "Gstt1", "Gfra3", "Ascl3", "Cftr", "Klk1",
                       "Pecam1", "Rgs5", "Tubb3", "Ncam1", 
                      "Adgre1", "Kit", "Gzma", "Vim", "Col1a1", "Twist1", "Alas2")

matrix.cell.markers <- AverageExpression(smg.integrated, features = unique(genes.to.heatmap1), assays = "RNA")
matrix.cell.markers <- matrix.cell.markers$RNA 

pdf("Heatmap e16-adult SMG integrated known cell markers.pdf", useDingbats = F, width = 6,height = 7)
pheatmap(matrix.cell.markers, scale = "row", cluster_cols = F, cluster_rows = F)
dev.off()

pdf("Dot Plot e16-adult SMG integrated known cell markers.pdf", useDingbats = F, width = 7,height = 5)
DotPlot(smg.integrated, features = genes.to.heatmap1, dot.min = 0.05, dot.scale = 4,group.by = "CellType.fixed") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 9), axis.title = element_blank(), axis.text.y = element_text(size = 9))
dev.off()


###### Determine cell type markers for all cell types in mouse SMG #####

DefaultAssay(smg.integrated) <- "RNA"
smg.integrated <- SetIdent(smg.integrated, value = "CellType.fixed")
smg.cell.markers <- FindAllMarkers(smg.integrated, min.pct = 0.2, only.pos = T)
smg.cell.markers <- smg.cell.markers[smg.cell.markers$p_val_adj<0.05, ]

write.csv(smg.cell.markers, file = "Cell-defining genes E16-adult SMG integrated.csv", row.names = T)

### Identify MEC markers and Basal progenitor markers
basalvsMEC <- FindMarkers(smg.integrated, ident.1 = "Basal duct", ident.2 = "Myoepithelial", logfc.threshold = 0.25, group.by = "CellType.fixed", only.pos = F)
basalvsMEC.adult <- FindMarkers(p300smg, ident.1 = "Basal duct", ident.2 = "Myoepithelial", logfc.threshold = 0.25, group.by = "CellType", only.pos = F)
basalmarkers <- FindMarkers(smg.integrated, ident.1 = "Basal duct", logfc.threshold = 0.25, only.pos = T)

DotPlot(smg.integrated, features = c("Ly6d","Krt14", "Krt17", "Sfn","Smoc2", "Cldn4", "Krt18", "Wfdc2", "Acta2", "Myl9", "Meg3", "Igfbp5", "Tagln"), dot.min = 0.05) + theme(axis.text.x = element_text(angle = 90), axis.title = element_blank())

##############################################################################################################
##############################################################################################################

###### Combine Human RNAseq with 10X single-cell RNAseq from mouse SMG #######

# Load human data
hPG.degs <- read.csv("../Human PG EdgeR analysis - DEG with log2 cpms.csv")
hSMG.degs <- read.csv("../EdgeR Both Glands combined/Human SMG EdgeR analysis - DEG with log2 cpms.csv")

#generate expression matrix with all genes for visualization
expression.matrix <- merge(hPG.degs[,c(2, 8:18)], hSMG.degs[,c(2,8:20)], by = "gene")
expression.matrix$gene <- tolower(expression.matrix$gene)
expression.matrix$gene <- capitalize(expression.matrix$gene)
rownames(expression.matrix) <- expression.matrix$gene
expression.matrix <- as.matrix(expression.matrix[,-1])

hPG.degs$gene <- tolower(hPG.degs$gene)
hPG.degs$gene <- capitalize(hPG.degs$gene)
hSMG.degs$gene <- tolower(hSMG.degs$gene)
hSMG.degs$gene <- capitalize(hSMG.degs$gene)

# Filter non-significant genes (p>0.05)
hPG.degs <- hPG.degs[hPG.degs$PValue<0.05, ]
hSMG.degs <- hSMG.degs[hSMG.degs$PValue<0.05, ]
hPG.degs <- hPG.degs[hPG.degs$logFC< -1 | hPG.degs$logFC> 1, ]
hSMG.degs <- hSMG.degs[hSMG.degs$logFC>1 | hSMG.degs$logFC< -1, ]

### Determine cellular localization of HUMAN IR-DEGs in 10x data

### Create subset with relevant populations
relevant.cells <- c("Seromucous acinar", "Ascl3+ duct", "Basal duct", "Serous acinar", "Endothelial",
                    "Erythroid", "GCT", "Glial cells", "Gfra3+ ID", "Krt19+ duct", "Macrophages", "Mast cells",
                    "Myoepithelial", "Nerves", "NK cells", "Gstt1+ ID",
                    "Pericytes", "Striated duct", "Stromal")

### Generate table of DEGs with cell-specific localization
relevant.subset <- subset(smg.integrated, idents = relevant.cells)
relevant.genes <- FindAllMarkers(relevant.subset, only.pos = T) #determine markers again without fetal-exclusive cells
FCvalues.10x.SMG <- reshape2::dcast(data = relevant.genes,formula = gene~cluster,fun.aggregate = sum,value.var = "avg_logFC")
SMG.DEGs.specific.10x <- base::merge(hSMG.degs, FCvalues.10x.SMG, by.x = "gene", by.y = "gene") #We use SMG only as 10x data does not contain PG
write.csv(SMG.DEGs.specific.10x, file = "Cellular localization of SMG-IR DEGs in scRNAseq.csv", row.names = T)


`%notin%` <- Negate(`%in%`)
colors.subset <- colors[names(colors) %notin% c("End bud", "Bpifa2+ Proacinar", "Smgc+ Proacinar", "Mitotic cells")]

counter<- data.frame("NoGenes" =  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), row.names = names(SMG.DEGs.specific.10x[21:39])) #we use columns 21:43 which contain cell data
for (i in 21:39) {
  counter[i-20,] <- length(which(SMG.DEGs.specific.10x[,i]>0))
}
counter$CellType <- rownames(counter)
counter <- counter[order(counter$NoGenes, decreasing = T),]
write.csv(counter, file = "Number of cell-enriched HUMAN DEGs in 10x SMG.csv", row.names = T)

p<-ggplot(data=counter, aes(x=reorder(CellType, -NoGenes), y=NoGenes)) +
  geom_bar(stat="identity", color="black", width = 0.85, fill="blue")

pdf("DEGs localization.pdf", useDingbats = F, width = 7,height = 5.5)
p + theme_classic() + theme(axis.text.x = element_text(colour = "black", size = 14, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14, colour = "black"), line = element_line(size = 0.5))
dev.off()

##### Cell-specific localization of Upstream Regulators (from IPA) #####
upregs <- read.csv("../Human Glands - Final Run/EdgeR Both Glands combined/Human Common Upstream Regulators.csv")
upregs$Upstream.Regulator <- tolower(upregs$Upstream.Regulator) 
upregs$Upstream.Regulator <- capitalize(upregs$Upstream.Regulator)
upregs <- upregs[order(upregs$Expr.Log.Ratio, decreasing = T), ]
upregs$Upstream.Regulator <- gsub("Tp63", "Trp63", upregs$Upstream.Regulator)

pdf("Upregs dotplot.pdf", useDingbats = F, width = 7.5, height = 6)
DotPlot(smg.integrated, features = upregs$Upstream.Regulator, 
        idents = nonfetal,
        dot.min = 0.05) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"), axis.title = element_blank())
dev.off()

## calculate average scaled expression of all genes from integrated dataset
genes.to.heatmap <- upregs$Upstream.Regulator #select Upstream regulator genes for visualization

upregs.avgexpression <- AverageExpression(object = relevant.subset, assays = "RNA", slot = "data", features = genes.to.heatmap)
upregs.avgexpression <- upregs.avgexpression$RNA

pdf(file = "Heatmap - Upreg gene expression.pdf", width = 5, height = 5, useDingbats = F)
pheatmap(upregs.avgexpression,scale = "row")
dev.off()

##### Ligand-Receptor analysis with Upstream regulators ######
# install LigandReceptor v1 package from GitHub
#devtools::install_github("chiblyaa/LigandReceptor")

# Load table of ligand-receptor pairs from published manuscript; https://doi.org/10.1038/ncomms8866
Ligand.Receptor.Pairs <- read.delim("../data/Ligand-Receptor Pairs.txt")
Ligand.Receptor.Pairs$Pair <- tolower(Ligand.Receptor.Pairs$Pair)
Ligand.Receptor.Pairs$Pair <- capitalize(Ligand.Receptor.Pairs$Pair)
Ligand.Receptor.Pairs$Ligand <- tolower(Ligand.Receptor.Pairs$Ligand)
Ligand.Receptor.Pairs$Ligand <- capitalize(Ligand.Receptor.Pairs$Ligand)
Ligand.Receptor.Pairs$Receptor <- tolower(Ligand.Receptor.Pairs$Receptor)
Ligand.Receptor.Pairs$Receptor <- capitalize(Ligand.Receptor.Pairs$Receptor)

Ligand.Receptor.Pairs <- data.frame(lapply(Ligand.Receptor.Pairs, function(x) {
  gsub("Ntf4", "Ntf5", x) 
}))

# create table of ligand-receptor pairs:
upregpairstable <- LigandReceptorPairsTable(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = upregs$Upstream.Regulator)
write.csv(upregpairstable, file = "Ligand-Receptor pairs with upstream regulators.csv")

# chord plot for MEC interactions
pdf("ChordPlot Upregs from MECs.pdf", useDingbats = F, height = 6,width = 9)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = upregs$Upstream.Regulator, cellcolors = colors.subset, from = "Myoepithelial")
legend("right",   # location of legend
       legend = names(colors), # categories or elements to render in
       # the legend
       fill = colors, bty = "n", cex = 0.8, xpd = TRUE) 
dev.off()

pdf("ChordPlot Upregs to MECs.pdf", useDingbats = F, height = 6,width = 9)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = upregs$Upstream.Regulator, cellcolors = colors.subset, to = "Myoepithelial")
legend("right",   # location of legend
       legend = names(colors), # categories or elements to render in
       # the legend
       fill = colors, bty = "n", cex = 0.8, xpd = TRUE) 
dev.off()

pdf("ChordPlot Upregs all.pdf", useDingbats = F, height = 6,width = 9)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = upregs$Upstream.Regulator, cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
legend("right",   # location of legend
       legend = names(colors), # categories or elements to render in
       # the legend
       fill = colors, bty = "n", cex = 0.8, xpd = TRUE) 
dev.off()


pdf("ChordPlot individual upregs.pdf", useDingbats = F, height = 6,width = 9)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Jag1", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Notch1", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
#PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Nostrin", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Ngfr", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Fgf10", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Igf2", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
#PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Sfrp2", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Itgb3", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
#PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Cxc3cl1", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Fgf1", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
#PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Etv5", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Ntrk2", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
#PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Snai2", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
#PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Trp63", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Bmp7", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
#PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Etv4", cellcolors = colors, from = relevant.cells, to= relevant.cells)
#PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Jag2", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
#PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Tfap2c", cellcolors = colors, from = relevant.cells, to= relevant.cells)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Fgfr1", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = "Pgf", cellcolors = colors.subset, from = relevant.cells, to= relevant.cells)
legend("right",   # location of legend
       legend = names(colors), # categories or elements to render in
       # the legend
       fill = colors, bty = "n", cex = 0.8, xpd = TRUE) 
dev.off()

pdf("ChordPlot from MECs all-genes.pdf", useDingbats = F, height = 6,width = 9)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, cellcolors = colors.subset, from = "Myoepithelial", to= relevant.cells)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, cellcolors = colors.subset, from = relevant.cells, to= "Myoepithelial")
legend("right",   # location of legend
       legend = names(colors), # categories or elements to render in
       # the legend
       fill = colors, bty = "n", cex = 0.8, xpd = TRUE) 
dev.off()

###### count highly expressed ligands and receptors per cell type based on defining genes ####
number.of.cellypes <- length(relevant.cells)
  
ligand_list <- vector(mode = "list", length = number.of.cellypes) # make a list of cell types
names(ligand_list) <- names(colors.subset)

for (i in 1:number.of.cellypes) {
  temptable <- relevant.genes[relevant.genes$cluster %in% names(ligand_list)[i],]
  ligand_list[[i]] <- temptable[temptable$gene %in% Ligand.Receptor.Pairs$Ligand, ]
}

receptor_list <- vector(mode = "list", length = number.of.cellypes)
names(receptor_list) <- names(colors.subset)

for (i in 1:number.of.cellypes) {
  temptable <- relevant.genes[relevant.genes$cluster %in% names(receptor_list)[i],]
  receptor_list[[i]] <- temptable[temptable$gene %in% Ligand.Receptor.Pairs$Receptor, ]
}

receptor.counter = replicate(number.of.cellypes, 0)
ligand.counter = replicate(number.of.cellypes, 0)
names(receptor.counter) <- names(ligand_list)
names(ligand.counter) <- names(ligand_list)
for (i in 1:number.of.cellypes) {
  receptor.counter[i]<- dim(receptor_list[[i]])[1]
  ligand.counter[i]<- dim(ligand_list[[i]])[1]
}
receptor.counter <- data.frame(receptor.counter)
ligand.counter <- data.frame(ligand.counter)

receptor.counter$CellType <- rownames(receptor.counter)
names(receptor.counter)[1] <-"NoGenes"
ligand.counter$CellType <- rownames(ligand.counter)
names(ligand.counter)[1] <-"NoGenes"

## graph with number of enriched receptors and ligands per cell type 
pr<-ggplot(data=data.frame(receptor.counter), aes(x=CellType, y=NoGenes)) +
  geom_bar(stat="identity", color="black", width = 0.85, fill="blue")
pl<-ggplot(data=data.frame(ligand.counter), aes(CellType, y=NoGenes)) +
  geom_bar(stat="identity", color="black", width = 0.85, fill="blue")

pdf("Number of enriched receptors and ligands per cell.pdf", useDingbats = F, width = 6,height = 4)
pr + theme_classic() + theme(axis.text.x = element_text(colour = "black", size = 11, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14, colour = "black"), line = element_line(size = 0.5))
pl + theme_classic() + theme(axis.text.x = element_text(colour = "black", size = 11, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14, colour = "black"), line = element_line(size = 0.5))
dev.off()


##### Ligand-receptor analysis with neurotrophin signaling genes in integrated dataset #####
ntfpairs <- LigandReceptorPairsTable(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = c("Ngf","Bdnf", "Ntf3", "Ntf5", "Ntrk1", "Ntrk2", "Ntrk3", "Ngfr"), from = relevant.cells, to = relevant.cells)
ntfpairs <- ntfpairs[ntfpairs$value>0, ]
write.csv(ntfpairs, file = "Neurotrophin signaling pairs - integrated dataset.csv")

pdf("ChordPlot Neurotrophin genes in nonfetal cells.pdf", useDingbats = F, height = 6,width = 9)
PairsPlot(seuratDEGS = relevant.genes, LRdatabase = Ligand.Receptor.Pairs, cellcolors = colors.subset, subsetgenes = c("Ngf","Bdnf", "Ntf3", "Ntf5", "Ntrk1", "Ntrk2", "Ntrk3", "Ngfr"), from = relevant.cells, to = relevant.cells)
legend("right",   # location of legend
       legend = names(subsetcolors), # categories or elements to render in
       # the legend
       fill = subsetcolors, bty = "n", cex = 0.8, xpd = TRUE)
dev.off()

ntfgenes = c("Ngf","Bdnf", "Ntf3", "Ntf5", "Ntrk1", "Ntrk2", "Ntrk3", "Ngfr")

pdf("dot plot neurotrophin genes.pdf", useDingbats = F, width = 5.5, height = 4)
DotPlot(relevant.subset, features = ntfgenes, idents = relevant.cells, dot.min = 0.05) + theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())
dev.off()

ntfs.avgexpression <- AverageExpression(object = relevant.subset, assays = "RNA", slot = "data", features = ntfgenes)
ntfs.avgexpression <- ntfs.avgexpression$RNA

pdf(file = "Heatmap - NTFS gene expression.pdf", width = 6, height = 4, useDingbats = F)
pheatmap(ntfs.avgexpression,scale = "row")
dev.off()

###### NEUROTROPHIN SIGNALING GENES IN HUMAN SAMPLES #####
### Neurotrophin genes were predicted as central players
### Generate box plots with neurotrophin genes across human samples
library(tidyverse)
library(hrbrthemes)
library(reshape2)
library(viridis)

SampleInfo <- read.delim("../data/SampleInfo.txt")
rownames(SampleInfo) <- SampleInfo$sample
drop<-c("s16_PAR", "s15_SMG", "s03_PAR", "s03_SMG") ##remove bad samples
SampleInfo <- SampleInfo[!(rownames(SampleInfo) %in% drop), ]

ntfData <- subset(expression.matrix, rownames(expression.matrix) %in% c("Ngfr", "Ntrk1","Ntrk2", "Ntrk3"))
ntfData <- ntfData[order(rownames(ntfData)),]
ntfData <- as.data.frame(t(ntfData))
ntfData$gland <- rownames(ntfData)

mntfData <- melt(ntfData, id=c("gland"))
ntfData <- merge(mntfData, SampleInfo, by.x  = "gland", by.y = "sample")

ntrkbox <- ggplot(ntfData, aes(x=variable, y=value, fill=condition)) + 
  geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.6) + theme(axis.text = element_text(size = 12, colour = "black"), axis.line = element_line(size = 0.5), panel.background = element_blank()) + xlab("") +ylab("Scaled expression value") + facet_wrap(~variable, scale="free", ncol = 4)

pdf("Box plots NTRK genes - human DEGs.pdf", height = 2.5, width = 7, useDingbats = F)
ntrkbox
dev.off()

ntfData <- subset(expression.matrix, rownames(expression.matrix) %in% c("Ngf", "Ntf3","Ntf4", "Bdnf"))
ntfData <- ntfData[order(rownames(ntfData)),]
ntfData <- as.data.frame(t(ntfData))
ntfData$gland <- rownames(ntfData)

mntfData <- melt(ntfData, id=c("gland"))
ntfData <- merge(mntfData, SampleInfo, by.x  = "gland", by.y = "sample")

ntfbox <- ggplot(ntfData, aes(x=variable, y=value, fill=condition)) + 
  geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.6) + theme(axis.text = element_text(size = 12, colour = "black"), axis.line = element_line(size = 0.5), panel.background = element_blank()) + xlab("") +ylab("Scaled expression value") + facet_wrap(~variable, scale="free", ncol = 4)

pdf("Box plots NTF genes - human DEGs.pdf", height = 2.5, width = 7, useDingbats = F)
ntfbox
dev.off()


save(e16smg, smg.integrated, hSMG.degs, upregs, Ligand.Receptor.Pairs, file = "Files for MEC analysis.RData")
