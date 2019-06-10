rm(list=ls())
library(ggplot2)
library(Seurat)
library(dplyr)
library(scater)
library(SingleCellExperiment)
library(pcaMethods)
library(pcaReduce)
library(monocle)
library(swne)

deng_SCE <- readRDS("deng-reads.rds")
deng_SCE$cell_type2 <- factor(
  deng_SCE$cell_type2,
  levels = c("zy", "early2cell", "mid2cell", "late2cell",
             "4cell", "8cell", "16cell", "earlyblast",
             "midblast", "lateblast")
)
cellLabels <- deng_SCE$cell_type2
deng <- counts(deng_SCE)
colnames(deng) <- cellLabels

d<- deng_SCE
d <- deng_SCE[!duplicated(rownames(d)), ]
colnames(d) <- 1:ncol(d)
geneNames <- rownames(d)
rownames(d) <- 1:nrow(d)

pd <- data.frame(timepoint = cellLabels)
pd <- new("AnnotatedDataFrame", data=pd)
fd <- data.frame(gene_short_name = rownames(deng_SCE))
fd <- new("AnnotatedDataFrame", data=fd)

deng_SCE <- newCellDataSet(as.matrix(logcounts(d)),phenoData = pd,featureData = fd,lowerDetectionLimit=0.5,expressionFamily=negbinomial.size())
deng_SCE <- detectGenes(deng_SCE, min_expr = 0.1)
print(head(fData(deng_SCE)))
expressed_genes <- row.names(subset(fData(deng_SCE), num_cells_expressed >= 2)) #genes expressed in at least 2 cells
print(length(expressed_genes))

deng_SCE <- estimateSizeFactors(deng_SCE)
deng_SCE <- estimateDispersions(deng_SCE) 
disp_table <- dispersionTable(deng_SCE)
ordering_genes <- subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1.2 * dispersion_fit)$gene_id#simluation 1.5
length(ordering_genes)# HEE: 0.5,2

d <- deng_SCE[which(rownames(deng_SCE) %in% ordering_genes),]
#d <- deng_SCE[!duplicated(rownames(d)), ]
colnames(d) <- 1:ncol(d)
geneNames <- rownames(d)
rownames(d) <- 1:nrow(d)
dim(d)

deng_SCE <- d
deng_SCE$cell_type2 <- cellLabels
#PCA
deng_PCA <- runPCA(deng_SCE)
plotPCA(deng_PCA, colour_by = "cell_type2")
#TSNE
deng_TSNE <- runTSNE(deng_SCE, rand_seed = 1)
plotTSNE(deng_TSNE, colour_by = "cell_type2")
#swne
swne.embedding <- RunSWNE(exprs(deng_SCE), k = 16)
PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = deng_SCE$cell_type2,
         do.label = T, label.size = 3.5, pt.size = 1.5, show.legend = F,
         seed = 42)
