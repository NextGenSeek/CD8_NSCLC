https://github.com/califano-lab/single-cell-pipeline

source('/PATH/TO/process-utils.R')
source('/PATH/TO/cluster-functions.R')
source('/PATH/TO/viper-utils.R')
library(ggplot2)
library(ggpubr)
library(viper)
library(pheatmap)
library(RColorBrewer)
library(MUDAN)
library(umap)
library(Matrix)
library(SeuratObject)
library(viper)
library(BisqueRNA)

df <- sce_CD8[["originalexp"]]@data
raw.mat <- as.matrix(df)
filt.mat <- QCTransform(raw.mat) # if needed
cpm.mat <- CPMTransform(filt.mat)

# run these if data was not normalized or transformed
cpm.mat <- CPMTransform(filt.mat)
rank.mat <- RankTransform(cpm.mat)
ARACNeTable(rank.mat, "/PATH/TO/CD8.tsv")

Run ARACNE
bash ./run_aracne.sh

#Load expression matrix used to make network.txt file in ARACNe-AP
count.mt <- read.table(file = "CD8.tsv", sep = "\t", header = T, quote = "")
#set gene symbol column as row names
rownames(count.mt) = make.names(count.mt$gene, unique = TRUE)
count.mt <- as.matrix(count.mt[,-1])
#Load network.txt file generated from ARACNe-AP
net <- read.table(file = "network.txt", sep = "\t", header = T, quote = "")
#Remove pvalue column
net2 <- net %>% dplyr::select(! pvalue)
#save with no row names or col names, insert end of line character
write_delim(net2, file = "network_all.txt", na = "NA", append = F, col_names = F, delim = "\t", quote = "none", escape = "none", eol = "\n")
#Set connection to pre-processed network file
net3 <- "/PATH/TO/network_all.txt"
#Make regulon object
regulons <- aracne2regulon(afile = net3, eset = count.mt, format = "3col")

cell_names_CD8 <- sce_CD8@meta.data$Sample
sce_CD8_exp <- SeuratToExpressionSet(sce_CD8)

SeuratObject.es = ExpressionSet(assayData = as.matrix(GetAssayData(sce_CD8)), 
                                phenoData =  new("AnnotatedDataFrame",
                                sce_CD8@meta.data))

signature <- rowTtest(SeuratObject.es, "CellType",c("Teff","Texh","TLdys"),"LEF1+ Tstem")$statistic
dnull <- ttestNull(SeuratObject.es, "CellType",c("Teff","Texh","TLdys"),"LEF1+ Tstem", repos=TRUE, per=100)
dim(dnull)

mrs <- msviper(signature, regulons, dnull, minsize = 25, adaptive.size = TRUE)
summary(mrs)
mrs_df <- as.data.frame(summary(mrs, mrs = 1053))

library(stringi)
mrs2 <- ledge(mrs)
summary(mrs2)
