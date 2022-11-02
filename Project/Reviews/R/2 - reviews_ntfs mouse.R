# :::::::::::::::::::::::::::::::::::::::::::: 
# Mouse scRNAseq
# :::::::::::::::::::::::::::::::::::::::::::: 

esmg = readRDS(file = "../../data/Embryonic SMG Integrated.rds")
psmg = readRDS(file = "../../data/Postnatal SMG Integrated.rds")

splito <- SplitObject(esmg, split.by = "stage")

e12smg = splito$E12
e14smg <- splito$E14
e16smg <- splito$E16

splito <- SplitObject(psmg, split.by = "stage")

p1smg <- splito$P1
p30smg <- splito$P30
p120smg <- splito$Adult

remove(e12smg, e14smg, splito, esmg, psmg)


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
DimPlot(e16smg, cols = colors)
mousegenes = c("Ngf", "Ntf3", "Ntf5", "Ntrk1", "Ntrk2", "Ntrk3", "Ngfr")

DotPlot(e16smg, features = mousegenes, dot.min = 0.05, col.min = 0, assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())


smg.integrated = readRDS("../../output/E16-Adult SMG integrated annotated.rds")
DotPlot(smg.integrated, features = mousegenes, dot.min = 0.05, col.min = 0, assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())


DotPlot(e16smg, features = mousegenes, dot.min = 0.05, col.min = 0, assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())
DotPlot(p1smg, features = mousegenes, dot.min = 0.05, col.min = 0, assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())

DotPlot(p30smg, features = mousegenes, dot.min = 0.05, col.min = 0, assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())
