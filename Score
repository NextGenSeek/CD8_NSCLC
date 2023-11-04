Idents(sce_CD8) <- sce_CD8@meta.data$majorCluster
ExhVsOther <- FindMarkers(sce_CD8, ident.1 = "CD8_C6-LAYN")
ExhVsOther$genes <- rownames(ExhVsOther)

Exh <- c("CXCL13", "HAVCR2", "CCL3", "TNFRSF9","ENTPD1")
Naive <- c("CCR7", "TCF7", "LEF1","SELL")
Cytotox <- c("PRF1", "IFNG", "GNLY", "NKG7","GZMB", "GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7")

# Select genes of interest 
gene.set <- Exh
# Get mean expression of genes of interest per cell
mean.exp <- colMeans(x = sce_CD8@assays$originalexp[gene.set, ], na.rm = TRUE)
# Add mean expression values in 'object@meta.data$gene.set.score'
if (all(names(x = mean.exp) == rownames(x = sce_CD8@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  sce_CD8@meta.data$Exh.score <- mean.exp
}
# Plot mean expression using Seurat::FeaturePlot()
FeaturePlot(object = sce_CD8, 
            reduction = "tsne",
            features = "Exh.score")
            
####### Cytotox/naiveness/exhaustion score for regulon expression
filename <- as.data.frame(t(sce_CD8@assays[["originalexp"]]@data))
filename_df <- data.frame(colData(sce)[, names(colData(sce)), drop=FALSE])
idx3 <- match(rownames(filename),rownames(filename_df))
ScoreGene <- filename_df
ScoreGene$PRDM1 <- filename$PRDM1
ScoreGene$TBX21 <- filename$TBX21
ScoreGene$ZEB2 <- filename$ZEB2
ScoreGene$TRPS1 <- filename$TRPS1
ScoreGene$E2F1 <- filename$E2F1
ScoreGene$ETV1 <- filename$ETV1
ScoreGene$ARNTL2 <- filename$ARNTL2
ScoreGene$JUN <- filename$JUN

ggplot(ScoreGene, aes(x= PRDM1, y=Naive.score, color=CellType)) +
  geom_point()  +
  geom_smooth(method= "loess", se=FALSE, aes(colour = NA), colour = "black") +
  stat_cor(method = "pearson", label.x = 0, label.y = 3) +
  scale_color_manual(values = c("Texh" = "hotpink3", 
                                "GZMK+ Tcms" = "blue4",
                                "Teff" = "khaki3",
                                "TLdys" = "springgreen4",
                                "CX3CR1+ Tems"="lightskyblue1",
                                "LEF1+ Tstem" = "indianred4",
                                "GZMH+ Trms"="lavenderblush3",
                                "XCL1+ Trms" = "lightsalmon",
                                "KLF2+ Tcms" = "mediumorchid3"))
