# Exploratory analysis (unpublished data)

## CORRELATION ANALYSIS ###
library(corrplot)
library(RColorBrewer)
library(tidyr)
library(psych)

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

gene<-as.numeric(expression.matrix["AQP5",])
correlations<-apply(expression.matrix,1,function(x){cor.test(gene,x)})
df <- data.frame(matrix(unlist(correlations), nrow=19888, byrow=T),stringsAsFactors=FALSE)
rownames(df) <- names(correlations)
colnames(df) <- names(correlations$`A1BG-AS1`)
names(df)[10] <- "conf.int2"
df$p.value <- as.numeric(df$p.value)
df$estimate <- as.numeric(df$estimate)
df <- df[df$p.value <0.01,]
df$gene <- rownames(df)
df <- df[abs(df$estimate) >0.7, ]

correlations2 <- as.data.frame(apply(expression.matrix,1,function(x){cor(gene,x)}))
correlations2$gene <- rownames(correlations2)
names(correlations2)[1] <- "estimate"
correlations2 <- correlations2[order(correlations2$estimate), ]
correlations2 <- correlations2[correlations2$gene %in% rownames(df),]

# matrix of the p-value of the correlation
genescorr <- rbind(tail(correlations2, 10), head(correlations2, 10))
genescorr <- genescorr[order(genescorr$estimate),]
t.expression.matrix <- subset(expression.matrix, rownames(expression.matrix) %in% genescorr$gene)

correlations3<-apply(t.expression.matrix,1,function(x){cor(gene,x)})
correlations4 <- as.data.frame(correlations3)
correlations4$gene <- rownames(correlations4)
correlations4 <- correlations2[order(correlations4$correlations3), ]

t.expression.matrix <- t(t.expression.matrix)
p.mat <- cor.mtest(t.expression.matrix)
head(p.mat[, 1:5])
test.run.subset <- cor(t.expression.matrix)
dev.off()
corrplot(test.run.subset, method="color",bg = "transparent",diag = T,outline = "black",addgrid.col = "black", tl.cex = 0.7,
         tl.col = "black", order = "FPC", 
         col = rev(brewer.pal(n=8, name="Spectral")),type = "upper",sig.level = 0.01, p.mat = p.mat, 
         insig = "pch", mar = c(0, 0, 0, 0)) 

as.data.frame(t(expression.matrix)) %>% 
  gather(key = variable, value = values, c("BHLHA15", "AQP5")) %>% 
  ggplot(aes(NGFR, values)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  scale_y_continuous(trans = "log2") + 
  scale_x_continuous(trans = "log2") + theme(axis.text = element_text(size = 12, colour = "black"), axis.line = element_line(size = 0.5), panel.background = element_blank()) + xlab("NGFR") + facet_wrap(~variable, scale="free")




## GSEA Analysis of pathway-enriched genes
Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

## Download GO annotation file for biological processes from MSigDB
GO.file = "../EdgeR Both Glands combined/c2.cp.reactome.v7.1.symbols.gmt"
go.pathways <- gmtPathways(GO.file)
head(go.pathways, 2)




######## CORRELATION ANALYSIS PAIRED WITH GSEA ##########
######## NTRK1 ##########
gene<-as.numeric(expression.matrix["NTRK1",])
correlations<-apply(expression.matrix,1,function(x){cor.test(gene,x)})
Ntrk1 <- data.frame(matrix(unlist(correlations), nrow=19888, byrow=T),stringsAsFactors=FALSE)
rownames(Ntrk1) <- names(correlations)
colnames(Ntrk1) <- names(correlations$`A1BG-AS1`)
names(Ntrk1)[10] <- "conf.int2"
Ntrk1$p.value <- as.numeric(Ntrk1$p.value)
Ntrk1$estimate <- as.numeric(Ntrk1$estimate)
Ntrk1 <- Ntrk1[Ntrk1$p.value <0.01,]
Ntrk1 <- Ntrk1[abs(Ntrk1$estimate) >0.75,]
Ntrk1$gene <- rownames(Ntrk1)
Ntrk1 <- Ntrk1[order(Ntrk1$p.value), ]
write.csv(Ntrk1, file = "NTRK1-correlated genes.csv")

### GSEA 
gene.list.GSEA <- resIRR[resIRR$gene %in% Ntrk1$gene, ] #Used CTRL vs IR stats for both glands combined
gene.list <- gene.list.GSEA$logFC
names(gene.list) <- gene.list.GSEA$gene
gene.list = sort(gene.list, decreasing = TRUE)
gene.list = gene.list[!duplicated(names(gene.list))]
GOres.Ntrk1 = fgsea(go.pathways, gene.list, minSize=15, maxSize = 500, nperm=1000)
head(GOres.Ntrk1[order(padj, -abs(NES)), ], n=10)
GOres.Ntrk1.10 <- head(GOres.Ntrk1[order(padj, -abs(NES)), ], n=10)
GOres.Ntrk1.10$leadingEdge <- as.character(GOres.Ntrk1.10$leadingEdge)
write.table(GOres.Ntrk1.10, file = "GO_Ntrk1 genes.txt", sep = "\t", row.names = F)
plotEnrichment(go.pathways[[GOres.Ntrk1.10[5]$pathway]], gene.list) 

fgseaRes <- GOres.Ntrk1
fgseaResTidy <- GOres.Ntrk1 %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score",
       title="GO processes (GSEA)") + 
  theme(axis.text = element_text(color = "black"), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_fill_brewer(palette = "Set2")

######## Ntrk2 ##########
gene<-as.numeric(expression.matrix["NTRK2",])
correlations<-apply(expression.matrix,1,function(x){cor.test(gene,x)})
Ntrk2 <- data.frame(matrix(unlist(correlations), nrow=19888, byrow=T),stringsAsFactors=FALSE)
rownames(Ntrk2) <- names(correlations)
colnames(Ntrk2) <- names(correlations$`A1BG-AS1`)
names(Ntrk2)[10] <- "conf.int2"
Ntrk2$p.value <- as.numeric(Ntrk2$p.value)
Ntrk2$estimate <- as.numeric(Ntrk2$estimate)
Ntrk2 <- Ntrk2[Ntrk2$p.value <0.01,]
Ntrk2 <- Ntrk2[abs(Ntrk2$estimate) >0.75,]
Ntrk2$gene <- rownames(Ntrk2)
Ntrk2 <- Ntrk2[order(Ntrk2$p.value), ]

### GSEA 
gene.list.GSEA <- resIRR[resIRR$gene %in% Ntrk2$gene, ] #Used CTRL vs IR stats for both glands combined
gene.list <- gene.list.GSEA$logFC
names(gene.list) <- gene.list.GSEA$gene
gene.list = sort(gene.list, decreasing = TRUE)
gene.list = gene.list[!duplicated(names(gene.list))]
GOres.Ntrk2 = fgsea(go.pathways, gene.list, minSize=15, maxSize = 500, nperm=1000)
head(GOres.Ntrk2[order(padj, -abs(NES)), ], n=10)
GOres.Ntrk2.10 <- head(GOres.Ntrk2[order(padj, -abs(NES)), ], n=10)
GOres.Ntrk2.10$leadingEdge <- as.character(GOres.Ntrk2.10$leadingEdge)
write.table(GOres.Ntrk2.10, file = "GO_Ntrk2 genes.txt", sep = "\t", row.names = F)
plotEnrichment(go.pathways[[GOres.Ntrk2.10[1]$pathway]], gene.list) 

fgseaRes <- GOres.Ntrk2
fgseaResTidy <- GOres.Ntrk2 %>%
  as_tibble() %>%
  arrange(desc(NES))

ggplot(fgseaResTidy[c(1:10),], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score",
       title="GO processes (GSEA)") + 
  theme(axis.text = element_text(color = "black"), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_fill_brewer(palette = "Set2")

######## Ntrk3 ##########
gene<-as.numeric(expression.matrix["NTRK3",])
correlations<-apply(expression.matrix,1,function(x){cor.test(gene,x)})
Ntrk3 <- data.frame(matrix(unlist(correlations), nrow=19888, byrow=T),stringsAsFactors=FALSE)
rownames(Ntrk3) <- names(correlations)
colnames(Ntrk3) <- names(correlations$`A1BG-AS1`)
names(Ntrk3)[10] <- "conf.int2"
Ntrk3$p.value <- as.numeric(Ntrk3$p.value)
Ntrk3$estimate <- as.numeric(Ntrk3$estimate)
Ntrk3 <- Ntrk3[Ntrk3$p.value <0.01,]
Ntrk3 <- Ntrk3[abs(Ntrk3$estimate) >0.75,]
Ntrk3$gene <- rownames(Ntrk3)
Ntrk3 <- Ntrk3[order(Ntrk3$p.value), ]

### GSEA 
gene.list.GSEA <- resIRR[resIRR$gene %in% Ntrk3$gene, ] #Used CTRL vs IR stats for both glands combined
gene.list <- gene.list.GSEA$logFC
names(gene.list) <- gene.list.GSEA$gene
gene.list = sort(gene.list, decreasing = TRUE)
gene.list = gene.list[!duplicated(names(gene.list))]
GOres.Ntrk3 = fgsea(go.pathways, gene.list, minSize=15, maxSize = 500, nperm=1000)
head(GOres.Ntrk3[order(padj, -abs(NES)), ], n=10)
GOres.Ntrk3.10 <- head(GOres.Ntrk3[order(padj, -abs(NES)), ], n=10)
GOres.Ntrk3.10$leadingEdge <- as.character(GOres.Ntrk3.10$leadingEdge)
write.table(GOres.Ntrk3.10, file = "GO_Ntrk3 genes.txt", sep = "\t", row.names = F)
plotEnrichment(go.pathways[[GOres.Ntrk3.10[1]$pathway]], gene.list) 

fgseaRes <- GOres.Ntrk3
fgseaResTidy <- GOres.Ntrk3 %>%
  as_tibble() %>%
  arrange(desc(NES))

ggplot(fgseaResTidy[c(1:10),], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score",
       title="GO processes (GSEA)") + 
  theme(axis.text = element_text(color = "black"), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_fill_brewer(palette = "Set2")

######## NGFR ##########
gene<-as.numeric(expression.matrix["NGFR",])
correlations<-apply(expression.matrix,1,function(x){cor.test(gene,x)})
NGFR <- data.frame(matrix(unlist(correlations), nrow=19888, byrow=T),stringsAsFactors=FALSE)
rownames(NGFR) <- names(correlations)
colnames(NGFR) <- names(correlations$`A1BG-AS1`)
names(NGFR)[10] <- "conf.int2"
NGFR$p.value <- as.numeric(NGFR$p.value)
NGFR$estimate <- as.numeric(NGFR$estimate)
NGFR <- NGFR[NGFR$p.value <0.01,]
NGFR <- NGFR[abs(NGFR$estimate) >0.75, ]
NGFR$gene <- rownames(NGFR)
NGFR <- NGFR[order(NGFR$p.value), ]

### GSEA 
gene.list.GSEA <- resIRR[resIRR$gene %in% NGFR$gene, ] #Used CTRL vs IR stats for both glands combined
gene.list <- gene.list.GSEA$logFC
names(gene.list) <- gene.list.GSEA$gene
gene.list = sort(gene.list, decreasing = TRUE)
gene.list = gene.list[!duplicated(names(gene.list))]
GOres.NGFR = fgsea(go.pathways, gene.list, minSize=15, maxSize = 500, nperm=1000)

head(GOres.NGFR[order(padj, -abs(NES)), ], n=10)
GOres.NGFR.10 <- head(GOres.NGFR[order(padj, -abs(NES)), ], n=10)
GOres.NGFR.10$leadingEdge <- as.character(GOres.NGFR.10$leadingEdge)
write.table(GOres.NGFR.10, file = "GO_NGFR genes.txt", sep = "\t", row.names = F)
plotEnrichment(go.pathways[[GOres.NGFR.10[1]$pathway]], gene.list) 

fgseaRes <- GOres.NGFR
fgseaResTidy <- GOres.NGFR %>%
  as_tibble() %>%
  arrange(desc(NES))

ggplot(fgseaResTidy[c(1:10),], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score",
       title="GO processes (GSEA)") + 
  theme(axis.text = element_text(color = "black"), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_fill_brewer(palette = "Set2")


######
#####
#####
common.corr.genes <- Ntrk1[Ntrk1$gene %in% Ntrk2$gene, ]
common.corr.genes <- common.corr.genes[common.corr.genes$gene %in% Ntrk3$gene, ]
common.corr.genes <- common.corr.genes[common.corr.genes$gene %in% NGFR$gene, ]




######## AQP5 ##########
gene<-as.numeric(expression.matrix["GFRA2",])
correlations<-apply(expression.matrix,1,function(x){cor.test(gene,x)})
AQP5 <- data.frame(matrix(unlist(correlations), nrow=19888, byrow=T),stringsAsFactors=FALSE)
rownames(AQP5) <- names(correlations)
colnames(AQP5) <- names(correlations$`A1BG-AS1`)
names(AQP5)[10] <- "conf.int2"
AQP5$p.value <- as.numeric(AQP5$p.value)
AQP5$estimate <- as.numeric(AQP5$estimate)
AQP5 <- AQP5[AQP5$p.value <0.01,]
AQP5 <- AQP5[abs(AQP5$estimate) >0.75, ]
AQP5$gene <- rownames(AQP5)

### GSEA 
gene.list.GSEA <- resIRR[resIRR$gene %in% AQP5$gene, ] #Used CTRL vs IR stats for both glands combined
gene.list <- gene.list.GSEA$logFC
names(gene.list) <- gene.list.GSEA$gene
gene.list = sort(gene.list, decreasing = TRUE)
gene.list = gene.list[!duplicated(names(gene.list))]
GOres.AQP5 = fgsea(pathways = go.pathways, gene.list, minSize=15, maxSize = 500, nperm=1000)
head(GOres.AQP5[order(padj, -abs(NES)), ], n=10)
GOres.AQP5.10 <- head(GOres.AQP5[order(padj, -abs(NES)), ], n=10)
GOres.AQP5.10$leadingEdge <- as.character(GOres.AQP5.10$leadingEdge)
write.table(GOres.AQP5.10, file = "GO_AQP5 genes.txt", sep = "\t", row.names = F)
plotEnrichment(go.pathways[[GOres.AQP5.10[3]$pathway]], gene.list) 

fgseaRes <- GOres.AQP5
fgseaResTidy <- GOres.AQP5 %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy[c(1:10),], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score",
       title="GO processes (GSEA)") + 
  theme(axis.text = element_text(color = "black"), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_fill_brewer(palette = "Set2")


