library(monocle)
library(SingleCellExperiment)
library(Seurat)
library(SeuratWrappers)
library(DESeq2)

sce <- as.SingleCellExperiment(sce_CD8)

# Monocle
cds <- as.CellDataSet(sce_CD8, assay = "originalexp")
cds

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

mk <- FindAllMarkers(sce_CD8,assay = "originalexp", min.pct = 0.3, min.diff.pct = 0.2, only.pos = TRUE)
mk[1:3,]

deg <- mk

deg <- deg[which(deg$cluster %in% unique(sce_CD8@active.ident)), ]
sel.gene <- unique(deg$gene)
cds <- monocle::setOrderingFilter(cds, sel.gene)

cds <- monocle::reduceDimension(cds, method = 'DDRTree')
cds <- monocle::orderCells(cds)

monocle::plot_cell_trajectory(cds, color_by = "sce_CD8@meta.data$Patient")+
  scale_color_manual(values = c("P0616A"= "orange3",
                                "P0616P"= "blue",
                                "P0617"= "green",
                                "P0619"= "yellow",
                                "P0913"= "red",
                                "P0729"= "cyan",
                                "P0706" = "darkmagenta",
                                "P1010" = "deepskyblue",
                                "P1011"= "deeppink2",
                                "P1118" = "darkolivegreen",
                                "P1120" = "darkgoldenrod2",
                                "P1202"= "darkturquoise",
                                "P1208" = "orange",
                                "P1219" = "black"))


monocle::plot_cell_trajectory(cds, color_by = "sce_CD8@active.ident") +
  scale_color_manual(values = c("Texh" = "hotpink3", 
                                "GZMK+ Tcms" = "blue4",
                                "Teff" = "khaki3",
                                "TLdys" = "springgreen4",
                                "CX3CR1+ Tems"="lightskyblue1",
                                "LEF1+ Tstem" = "indianred4",
                                "GZMH+ Trms"="lavenderblush3",
                                "XCL1+ Trms" = "lightsalmon",
                                "KLF2+ Tcms" = "mediumorchid3"))


monocle::plot_cell_trajectory(cds, color_by = "State")
monocle::plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~State)
monocle::plot_cell_trajectory(cds, color_by = "Pseudotime")

############################
# Each patient separately
############################

"sce_CD8_P" <- sce_CD8  # Change idents for patient figure
Idents(sce_CD8_P) <- "Patient"

"P1208" 

sce_P1208 <- subset(sce_CD8_P, idents ="P1208")
Obj <- sce_P1208
cds <- as.CellDataSet(Obj, assay = "originalexp")
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
mk <- FindAllMarkers(Obj,assay = "originalexp", min.pct = 0.3, min.diff.pct = 0.2, only.pos = TRUE)
deg <- mk
deg <- deg[which(deg$cluster %in% unique(Obj@meta.data$CellType)), ]
sel.gene <- unique(deg$gene)
cds <- monocle::setOrderingFilter(cds, sel.gene)
cds <- monocle::reduceDimension(cds, method = 'DDRTree')
cds <- monocle::orderCells(cds)
monocle::plot_cell_trajectory(cds, color_by = "Obj@meta.data$CellType") +
  scale_color_manual(values = c("Texh" = "hotpink3", 
                                "GZMK+ Tcms" = "blue4",
                                "Teff" = "khaki3",
                                "TLdys" = "springgreen4",
                                "CX3CR1+ Tems"="lightskyblue1",
                                "LEF1+ Tstem" = "indianred4",
                                "GZMH+ Trms"="lavenderblush3",
                                "XCL1+ Trms" = "lightsalmon",
                                "KLF2+ Tcms" = "mediumorchid3"))
