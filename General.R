# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(AnnotationHub)
library(ensembldb)
library(stringr)
library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)
library(data.table)
library(gridExtra)
library("sscClust")
library(SingleR)
library(RColorBrewer)

tung_counts <- read.table("GSE99254_NSCLC.TCell.S12346.count.txt", sep = "\t", header = TRUE, fill=TRUE)
# There are some rows that don't have a gene symbol, which gives error if we want to use this column as row names.
tung_counts$symbol[is.na(tung_counts$symbol)] <- tung_counts$geneID[is.na(tung_counts$symbol)]
tung_counts2 <- tung_counts[,-1]
tung_counts3 <- tung_counts2[,-1]
rownames(tung_counts3) <- tung_counts2[,1]
tung_counts4 <- tung_counts3[ , order(names(tung_counts3))] # Both files have to be in the same order

tung_annotation <- read.table("CellInfo.txt", sep = "\t", header = TRUE)
tung_annotation2 <- tung_annotation[order(tung_annotation$sample), ] 
# Difference in naming: one file uses . instead of -
tung_annotation2$sample <- str_replace(tung_annotation2$sample,"-", ".")
tung_annotation2$sample <- str_replace(tung_annotation2$sample,"-", ".")
tung_annotation2$sample <- str_replace(tung_annotation2$sample,"-", ".")
# Both files have to be in the same order
tung_annotation3 <- tung_annotation2[,-1]
rownames(tung_annotation3) <- tung_annotation2[,1]

# note that the data passed to the assay slot has to be a matrix!
tung <- SingleCellExperiment(
  assays = list(counts = as.matrix(tung_counts4)),
  colData = tung_annotation3
)

sce_Seu <- as.Seurat(tung, data= NULL)

#######################################
# SELECT THE CD8 CELLS
#######################################
sce_CD8 <- subset(x = sce_Seu, subset = type == "CD8")
sce_CD8$Sample <- Idents(sce_CD8)

# Identification of highly variable features
sce_CD8 <- FindVariableFeatures(object= sce_CD8, nfeatures = 1500)
# Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(sce_CD8), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sce_CD8)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2

# Shifts the expression of each gene, so that the mean expression across cells is 0
sce_CD8 <- ScaleData(object = sce_CD8)
#Perform linear dimensional reduction. We perform PCA on the scaled data. 
sce_CD8 <- RunPCA(object = sce_CD8, features = VariableFeatures(object = sce_CD8))
print(sce_CD8[["pca"]], dims = 1:9, nfeatures = 30)
sce_CD8 <- RunTSNE(object = sce_CD8, tsne.method = "Rtsne")

sce_CD8 <- FindNeighbors(object = sce_CD8, dims = 1:50)

sce_CD8 <- FindClusters(object= sce_CD8)
sce_CD8 <- RunUMAP(object = sce_CD8, dims = 1:20)

DimPlot(sce_CD8, reduction = "umap", label = TRUE)
DimPlot(sce_CD8, reduction = "tsne", label = TRUE)

#Create plot
DimPlot(sce_CD8, reduction = "tsne",group.by = "sampleType", cols = c("red", "orange", "blue"))
