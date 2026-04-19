library(RTN)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(snow)
library(RedeR)
library(tibble)
library(Rmisc)
library(ggpubr)

# Download dorothea_hs_pancancer

dorothea_hs_pancancer$source <- dorothea_hs_pancancer$tf
DoroPan <- dorothea_hs_pancancer[,-1]

TFS_Doro <- as.data.frame(DoroPan$source)
regulatoryElements <- TFS_Doro$`DoroPan$source`

#################################################

df <- sce_CD8[["originalexp"]]@data
raw.mat <- as.matrix(df)
filt.mat <- QCTransform(raw.mat) # if needed
cpm.mat <- CPMTransform(filt.mat)

# Run the TNI constructor
rtni_all <- tni.constructor(expData = cpm.mat, 
                             regulatoryElements = regulatoryElements)

# Compute the reference regulatory network by permutation and bootstrap analyses.
options(cluster=snow::makeCluster(spec=4, "SOCK"))
rtni_all <- tni.permutation(rtni_all, nPermutations = 100, pValueCutoff = 1e-3)
rtni_all <- tni.bootstrap(rtni_all)
stopCluster(getOption("cluster"))

# Compute the DPI-filtered regulatory network: ARACNe algorithm
rtni_all <- tni.dpi.filter(rtni_all, eps = 0)

# summary
tni.regulon.summary(rtni_all)
regulons_all <- tni.get(rtni_all, what = "regulons.and.mode")

############
ExprMt <- GetAssayData(sce_CD8, slot = "data")
ExprDf <- as.data.frame(ExprMt)

Test4 <- unlist(ExprDf[,1])

names(Test4) <- rownames(raw.mat)

cluster1.markers <- FindMarkers(sce_CD8, 
                                ident.1= c("LEF1+ Tstem"),
                                ident.2= "Teff",
                                logfc.threshold = 0.25,
                                min.cells.feature = 3,
                                min.cells.group = 3,
                                min.pct = 0.1) # Settings to give all genes that were detected in one of the groups

markersTest <- cluster1.markers$avg_log2FC
names(markersTest) <- rownames(cluster1.markers)

rtna <- tni2tna.preprocess(object = rtni, 
                           phenotype = Test4, 
                           hits = regulatoryElements,
                           phenoIDs = NULL,
                           verbose = TRUE)

rtna <- tna.mra(rtna, pValueCutoff = 0.05, minRegulonSize = 10, tfs = regulatoryElements)

mra <- tna.get(rtna, what="mra", ntop = -1)
head(mra, n=10)

mra$Pvalue2 <- log(mra$Pvalue, base=exp(1))

ggplot(mra, aes(x = reorder(Regulon, Pvalue2), y = Pvalue2)) + 
  geom_bar(stat = "identity") +
  scale_y_continuous(trans='log2') +
  scale_y_reverse() +
  ylab("-log(Pvalue)") +
  theme(axis.text.x = element_text(size = 12, angle = 90)) +
  xlab("")

# The absolute value of a weight represents the MI value, while the sign (+/-) indicates the predicted mode of action based on the Pearson’s correlation between the regulator and its targets.
g <- tni.graph(rtni_cyto, regulatoryElements = c("gene1","gene2","gene3"))

# Plot most important regulons
rdp <- RedPort()
calld(rdp)
addGraph(rdp, g, layout=NULL)
addLegend.color(rdp, g, type="edge")
addLegend.shape(rdp, g)
relax(rdp, ps = TRUE)

#################################################

#Comparison of the top TFs

mra_dys$Name <- "Dysfunctional"
mra_exh$Name <- "Exhausted"
mra_cyto$Name <- "Cytotoxic"

Combined <-rbind(mra_dys,mra_exh,mra_cyto)
Combined$Pvalue2 <- log(Combined$Pvalue, base=exp(1))

# Select only regulons that were >P0.5 in one of the groups
ListRegulons <- c("gene1","gene2","gene3")
Combined_short <- dplyr::filter(Combined, Regulon %in% ListRegulons)
Combined_short <- transform(Combined_short, Regulon= factor(Regulon, levels=ListRegulons))

ggplot(Combined_short, aes(x = Regulon, y = Pvalue2, fill = Name)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_y_continuous(trans='log2') +
  scale_y_reverse() +
  ylab("-log(Pvalue)") +
  theme(axis.text.x = element_text(size = 12, angle = 90)) +
  xlab("") +
  scale_fill_brewer(palette="Accent")

#################################################
# Collect gene information data regulons

sce_combined <- as.SingleCellExperiment(sce_CD8_combined)
SCE_combined <- logcounts(sce_combined)  # access log-transformed counts matrix
SCEmt_combined <- as.matrix(SCE_combined)
cellLabels_combined <- sce_CD8_combined@active.ident
Cell_Labels <- as.data.frame(cellLabels_combined)

TPM_reg_combined <- as.data.frame(SCEmt_combined[rownames(SCEmt_combined) %in% ListRegulons, ])
TPM_reg_combined$Gene <- rownames(TPM_reg_combined)
TPM <- as.data.frame(t(TPM_reg_combined))
# Add CellType information
idx <- match(rownames(TPM), rownames(Cell_Labels))
TPM$CellType <- Cell_Labels$cellLabels_combined[idx]

TPM_df <- as.data.frame(TPM)
TPM_df$CellType <- as.character(TPM_df$CellType)
TPM <- TPM_df$CellType %>%
  dplyr::mutate(Branch = case_when(
    endsWith(CellType, "CX3CR1+ Tems") ~ "cytotoxic",
    startsWith(CellType, "Teff") ~ "cytotoxic",
    startsWith(CellType, "Texh") ~ "exhausted",
    #  startsWith(CellType, "GZMH+ Trms") ~ "exhausted",
    startsWith(CellType, "TLdys") ~ "dysfunctional"
  ))

TPM2 <- TPM[-425,]
TPM_T <- pivot_longer(TPM2,c(1:17))
TPM_T$value <- as.numeric(TPM_T$value)

ForDots <- summarySE(TPM_T, measurevar="value", groupvars=c("CellType","name"))
DF <-as.data.frame(ForDots)

ggplot(DF) +
  geom_bar( aes(x=name, y= value, fill=CellType), position="dodge", stat="identity") +
  geom_errorbar( aes(x=name, ymin= value-sd, ymax=value+sd), width=0.4, colour="grey", alpha=0.9, linewidth=0.7, position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
