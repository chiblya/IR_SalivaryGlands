library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(pheatmap)
library(superheat)
library(Hmisc)
library(reshape2)
library(circlize)
library (monocle)
library(dplyr)
library(cowplot)
library(htmlwidgets)
library(SeuratWrappers)
library(SingleCellExperiment)
library(loomR)
library(scater)
library(patchwork)

##### Load necessary files #####

load("Files for MEC analysis.RData")

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

###### Trajectory analysis with E16 SMG #####
DimPlot(e16smg)
e16smg <- SetIdent(e16smg, value = "CellType")
e16smg <- RenameIdents(e16smg, "Mesenchyme" = "Stromal", "Smooth muscle" = "Pericytes")
Idents(e16smg) <- factor(Idents(e16smg), levels = cell.type.levels)
e16smg$CellType <- Idents(e16smg)

epi.subset <- smg.integrated[,  smg.integrated$CellType %in% c("Myoepithelial")] #, "End bud", "Krt19+ duct", "Basal duct"
DimPlot(epi.subset)

epi.sce <- as.SingleCellExperiment(epi.subset)

plotExpression(epi.sce, features = "Acta2", x = "stage") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plotPCA(epi.sce, colour_by = "stage")
plotUMAP(epi.sce, colour_by = "stage")


class(epi.sce)
structure(epi.sce)
table(epi.sce$CellType)

# Re-order the levels of the factor storing the cell developmental stage.
epi.sce$stage <- factor(epi.sce$stage, levels = c("E16", "P1", "P30", "Adult"))
# Run PCA on epi data. Use the runPCA function from the SingleCellExperiment package.
epi.sce <- runPCA(epi.sce, ncomponents = 50)
# Use the reducedDim function to access the PCA and store the results. 
pca <- reducedDim(epi.sce, "PCA")
# Describe how the PCA is stored in a matrix. Why does it have this structure?
head(pca)

# Add PCA data to the epi.sce object.
epi.sce$PC1 <- pca[, 1]
epi.sce$PC2 <- pca[, 2]

# Plot PC biplot with cells colored by stage. 
# colData(epi.sce) accesses the cell metadata DataFrame object for epi.sce.
# Look at Figure 1A of the paper as a comparison to your PC biplot.
library(vipor)
library(ggbeeswarm)
library(ggthemes)

ggplot(as.data.frame(colData(epi.sce)), aes(x = PC1, y = PC2, color = stage)) + geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")

# PCA is a simple approach and can be good to compare to more complex algorithms 
# designed to capture differentiation processes. As a simple measure of pseudotime 
# we can use the coordinates of PC1.
# Plot PC1 vs stage. 
epi.sce$pseudotime_PC1 <- rank(epi.sce$PC1)  # rank cells by their PC1 score
ggplot(as.data.frame(colData(epi.sce)), aes(x = pseudotime_PC1, y = stage, 
                                             colour = stage)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")

#  Prepare a counts matrix with labeled rows and columns. 
epi <- logcounts(epi.sce)  # access log-transformed counts matrix
cellLabels <- epi.sce$stage
colnames(epi) <- cellLabels

# Make a diffusion map.
library(destiny)
dm <- DiffusionMap(epi.sce)
sigmas <- find_sigmas(data = t(as.matrix(epi)))
dm <- DiffusionMap(data = epi.sce, sigma = optimal_sigma(sigmas),n_pcs = 50)


# Optional: Try different sigma values when making diffusion map.
# dm <- DiffusionMap(t(epi), sigma = "local")  # use local option to set sigma
# sigmas <- find_sigmas(t(epi), verbose = FALSE)  # find optimal sigma
# dm <- DiffusionMap(t(epi), sigma = optimal_sigma(sigmas))  

# Plot diffusion component 1 vs diffusion component 2 (DC1 vs DC2). 
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  Timepoint = epi.sce$stage)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point() + scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()
# Try plotting higher diffusion components against one another.

# Next, let us use the first diffusion component (DC1) as a measure of pseudotime.
# How does the separation by cell stage look?
epi.sce$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])    # rank cells by their dpt
ggplot(as.data.frame(colData(epi.sce)), 
       aes(x = pseudotime_diffusionmap, 
           y = stage, colour = stage)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Diffusion component 1 (DC1)") + ylab("Timepoint") +
  ggtitle("Cells ordered by DC1")

plot(eigenvalues(dm), ylim = 0:1, pch = 20, xlab = 'Diffusion component (DC)', ylab = 'Eigenvalue')

# What happens if you run the diffusion map on the PCs? Why would one do this?
rownames(pca) <- cellLabels
dm <- DiffusionMap(pca)

# Diffusion pseudotime calculation. 
# Set index or tip of pseudotime calculation to be a zygotic cell (cell 268). 
dpt <- DPT(dm, tips = 1)


# Plot DC1 vs DC2 and color the cells by their inferred diffusion pseudotime.
# We can accesss diffusion pseudotime via dpt$dpt.
df <- data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2], 
                 dptval = dpt$dpt, stage = epi.sce$stage)
p1 <- ggplot(df) + geom_point(aes(x = DC1, y = DC2, color = dptval))
p2 <- ggplot(df) + geom_point(aes(x = DC1, y = DC2, color = stage))
p <- plot_grid(p1, p2)
p

# Plot diffusion pseudotime vs timepoint. 
# Which separates the data better, DC1 or diffusion pseudotime?
epi.sce$pseudotime_dpt <- rank(dpt$dpt) 
pdf("Cells ordered by Diffusion map pseudotime.pdf", useDingbats = F, width = 4.5,height = 3)
ggplot(as.data.frame(colData(epi.sce)), 
       aes(x = pseudotime_dpt, 
           y = stage, colour = stage)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Diffusion map pseudotime (dpt)") +
  ylab("Timepoint") +
  ggtitle("Cells ordered by diffusion map pseudotime")
dev.off()

saveRDS(epi.sce, file = "MECs SCE file pre-Slingshot.rds")

#### Slingshot #####
library(slingshot)
library(Seurat)

# Read the Slingshot documentation (?slingshot) and then run Slingshot below. 
# Given your understanding of the algorithm and the documentation, what is one 
# major set of parameters we omitted here when running Slingshot?
sce <- slingshot(epi.sce, reducedDim = 'PCA', clusterLabels = epi.sce$stage)  # no clusters

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(50, alpha = 1)

pdf("PCA plot slingshot.pdf", useDingbats = F, height = 5, width = 6)

plot(reducedDims(sce)$PCA, col = hcl.colors(50)[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1) 
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black')
legend('topright', title = 'Pseudotime', col = hcl.colors(5), legend=c('0','15','30', "45", "60"), pch=16)
dev.off()

pdf("MECS PCA plot colored by stage.pdf", useDingbats = F, height = 3, width = 4)
ggplot(as.data.frame(colData(sce)), aes(x = PC1, y = PC2, color = stage)) + geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")
dev.off()

# Plot Slingshot pseudotime vs cell stage. 

pdf("Cells ordered by Slingshot pseudotime.pdf", useDingbats = F, width = 4.5,height = 3)
ggplot(as.data.frame(colData(epi.sce)), aes(x = sce$slingPseudotime_1, y = stage, 
                                             colour = stage)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")
dev.off()


###### DIFFERENTIAL EXPRESSION ########
library(gam)

# Only look at the 1,000 most variable genes when identifying temporally expressesd genes.
# Identify the variable genes by ranking all genes by their variance.
Y <- log2(counts(deng_SCE) + 1)
var1K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:1000]
Y <- Y[var1K, ]  # only counts for variable genes

# Fit GAM for each gene using pseudotime as independent variable.
t <- deng_SCE$slingPseudotime_1
gam.pval <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

# Identify genes with the most significant time-dependent model fit.
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]  

# Prepare and plot a heatmap of the top genes that vary their expression over pseudotime.
require(clusterExperiment)
heatdata <- as.matrix(gcdata@data[rownames(gcdata@data) %in% topgenes, order(t, na.last = NA)])
heatclus <- gcdata@ident[order(t, na.last = NA)]
png(paste0(mydir, "heatmap_time_genes.png"), width=10, height=10, units = "in", res=200)
ce <- ClusterExperiment(heatdata, heatclus, transformation = log1p)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed', cexRow = 1.5, fontsize = 15)
dev.off()

##### Ligand-receptor analysis with neurotrophin signaling genes in E16SMG #####


e16smg.markers <- FindAllMarkers(e16smg, logfc.threshold = 0.25,only.pos = T)

E16ntfpairs <- LigandReceptorPairsTable(seuratDEGS = e16smg.markers, LRdatabase = Ligand.Receptor.Pairs, subsetgenes = c("Ngf","Bdnf", "Ntf3", "Ntf5", "Ntrk1", "Ntrk2", "Ntrk3", "Ngfr"))
E16ntfpairs <- E16ntfpairs[E16ntfpairs$value>0, ]
write.csv(E16ntfpairs, file = "E16 Neurotrophin signaling pairs.csv")

pdf("ChordPlot Neurotrophin genes in E16 SMG.pdf", useDingbats = F, height = 6,width = 9)
PairsPlot(seuratDEGS = e16smg.markers, LRdatabase = Ligand.Receptor.Pairs, cellcolors = colors[names(colors) %in% e16smg.markers$cluster], subsetgenes = c("Ngf","Bdnf", "Ntf3", "Ntf5", "Ntrk1", "Ntrk2", "Ntrk3", "Ngfr"))
legend("right",   # location of legend
       legend = names(subsetcolors), # categories or elements to render in
       # the legend
       fill = subsetcolors, bty = "n", cex = 0.8, xpd = TRUE)

dev.off()




##### MEC markers #####
## determine IR-DEGs enriched in MECs
mec.hSMG.degs <- merge(relevant.genes[relevant.genes$cluster %in% "Myoepithelial", ],  hSMG.degs, by = "gene") #61 GENES
write.csv(mec.hSMG.degs, "MEC genes affected in irradiated human SMG.csv") 
mec.hSMG.degs <- mec.hSMG.degs[order(mec.hSMG.degs$avg_logFC, decreasing = T), ]

MEC.avgexpression <- AverageExpression(object = relevant.subset, assays = "RNA", slot = "data", features = mec.hSMG.degs$gene)
MEC.avgexpression <- MEC.avgexpression$RNA

pdf(file = "Heatmap - MEC DEGS gene expression.pdf", width =5, height = 10, useDingbats = F)
pheatmap(MEC.avgexpression,scale = "row")
dev.off()

pdf("Dot Plot of 61 MEC genes affected in human IR SMG.pdf", useDingbats = F, width = 13,height = 6)
DotPlot(relevant.subset, features = mec.hSMG.degs$gene, dot.min = 0.10) + theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())
dev.off()

## cross-ref MEC-genes with upstream regulators
MEC.upregs <- merge(MECmarkers, upregs, by.x = "gene", by.y = "Upstream.Regulator")


####### Analysis of MEC development #######
#### Identify genes in MEC differentiation ####
smg.integrated$celltype.stage <- paste0(smg.integrated$stage, "_", smg.integrated$CellType.fixed)
smg.integrated <- SetIdent(smg.integrated, value = "celltype.stage")
DimPlot(smg.integrated)

MECmarkers <- relevant.genes[relevant.genes$cluster %in% "Myoepithelial", ]
MECmarkers <- MECmarkers[MECmarkers$p_val_adj<0.05, ]

DefaultAssay(smg.integrated) <- "RNA"
MEC.e16.to.p1 <- FindMarkers(smg.integrated, ident.1 = "E16_Myoepithelial", ident.2 = "P1_Myoepithelial", only.pos = F,min.pct = 0.25)
MEC.e16.to.p1 <- MEC.e16.to.p1[MEC.e16.to.p1$p_val_adj<0.05,]
MEC.e16.to.p1 <- MEC.e16.to.p1[rownames(MEC.e16.to.p1) %in% MECmarkers$gene, ]
MEC.e16.to.p1$gene <- rownames(MEC.e16.to.p1)

MEC.p1.to.adult <- FindMarkers(smg.integrated, ident.1 = "P1_Myoepithelial", ident.2 = c("P30_Myoepithelial", "Adult_Myoepithelial"), only.pos = F,min.pct = 0.25)
MEC.p1.to.adult <- MEC.p1.to.adult[MEC.p1.to.adult$p_val_adj<0.05,]
MEC.p1.to.adult <- MEC.p1.to.adult[rownames(MEC.p1.to.adult) %in% MECmarkers$gene, ]
MEC.p1.to.adult$gene <- rownames(MEC.p1.to.adult)

write.csv(MEC.e16.to.p1, file = "DEGs MECs E16 vs P1.csv")
write.csv(MEC.p1.to.adult, file = "DEGs MECs P1 vs P30_adult.csv")

###### create MEC subset ####
smg.integrated <- SetIdent(smg.integrated, value = "CellType.fixed")
DimPlot(smg.integrated)
mec.subset<-subset(smg.integrated, idents = c("Myoepithelial"))
DimPlot(mec.subset)
## Order identities for visualization
mec.subset <- SetIdent(mec.subset, value = "stage")
Idents(mec.subset) <- factor(Idents(mec.subset), levels = c("E16", "P1", "P30", "Adult"))
mec.subset$stage <- Idents(mec.subset)
DimPlot(mec.subset)


## Violin plots of selected genes from heatmap
DefaultAssay(mec.subset) <- "integrated"

pdf("Violin plots selected MEC DEGS.pdf", useDingbats = F, width = 9, height = 2.5)
VlnPlot(mec.subset, features = c("Krt15", "Smoc2", "Ucma", "Smr3a", "2010007H06Rik"), group.by = "stage", pt.size = 0, cols = stage.colors, ncol = 5)
VlnPlot(mec.subset, features = c("Col9a3", "Cpm", "Ntng1", "Wnt6", "Col9a2"), group.by = "stage", pt.size = 0, cols = stage.colors, ncol = 5)
VlnPlot(mec.subset, features = c("Ntrk3", "Spon2", "Sfrp1", "Pgf", "Il17b"), group.by = "stage", pt.size = 0, cols = stage.colors, ncol = 5)
VlnPlot(mec.subset, features = c("Emid1", "Cryab", "Col4a2", "Pdgfa", "Igfbp5"), group.by = "stage", pt.size = 0, cols = stage.colors, ncol = 5)
VlnPlot(mec.subset, features = c("Nrtn", "Ntf3", "Ngf", "Fgfr1", "Dlk2"), group.by = "stage", pt.size = 0, cols = stage.colors, ncol = 5)
VlnPlot(mec.subset, features = c("Krt14", "Acta2", "Cnn1", "Myh11"), group.by = "stage", pt.size = 0, cols = stage.colors, ncol = 5)
dev.off()






################
## Load file from transcription factors
## download list of transcription factors from Schmeier et al (2016): https://tools.sschmeier.com/tcof/home/
Transcription.factors <- read.csv("data/TranscriptionFactors.csv")
names(Transcription.factors)[1] <- "Gene"







