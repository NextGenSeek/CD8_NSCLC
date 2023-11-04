LUAD <- read.delim("/PATH/TO/Human__TCGA_LUAD__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct")
GTEX <- read.delim("/PATH/TO/bulk-gex-v8-rna-seq-counts-by-tissue-gene_reads_2017-06-05_v8_lung.gct", sep = '')

LUAD_TPM = apply(LUAD2, 2, function(x) x / sum(as.numeric(x)) * 10^6) %>% as.data.frame()

############################################################
#Average expression master regulons
############################################################

LUAD_short <- LUAD[LUAD$attrib_name %in% ListRegulons,]
GTEX_short <- GTEX[GTEX$Description %in% ListRegulons,]

# Calculate mean and SD
LUAD_short$Mean <-apply(LUAD_short[2:516],1,mean)
LUAD_short$SD <-apply(LUAD_short[2:516],1,sd)
GTEX_short$Mean <-apply(GTEX_short[4:581],1,mean)
GTEX_short$SD <- apply(GTEX_short[4:581], 1, sd)

# Make figure
ggplot(LUAD_short) +
  geom_bar( aes(x=attrib_name, y= Mean), stat="identity") +
  geom_errorbar( aes(x=attrib_name, ymin= Mean-SD, ymax=Mean+SD), width=0.4, colour="grey", alpha=0.9, linewidth=0.7) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(trans='log10')

ggplot(GTEX_short) +
  geom_bar( aes(x=Description, y= Mean), stat="identity") +
  geom_errorbar( aes(x=Description, ymin= Mean-SD, ymax=Mean+SD), width=0.4, colour="grey", alpha=0.9, linewidth=0.7) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(trans='log10')

############################################################
#Average expression for the signature genes

# Extract target genes
ListRegulons <- c("gene","gene","gene")
Reg_gene <- subset(DoroPan, source=="gene", select = target)
LUAD_gene <- LUAD[LUAD$attrib_name %in% Reg_gene$target,]
GTEX_gene <- GTEX[GTEX$Description %in% Reg_gene$target,]

LUAD_gene$Mean <-apply(LUAD_gene[2:516],1,mean)
LUAD_gene$SD <-apply(LUAD_gene[2:516],1,sd)
LUAD_gene$Regulon <- "gene"

LUAD_combined_targets <- rbind(LUAD_gene1,LUAD_gene2,LUAD_gene3)

GTEX_gene$Mean <-apply(GTEX_gene[4:581],1,mean)
GTEX_gene$SD <-apply(GTEX_gene[4:581],1,sd)
GTEX_gene$Regulon <- "gene"

GTEX_combined_targets <- rbind(GTEX_gene1,GTEX_gene2,GTEX_gene3)

ggplot(LUAD_combined_targets, aes(x= Regulon, y= Mean) ,stat="identity") + 
  geom_boxplot(outlier.colour="black", outlier.shape = 20) +
  scale_y_continuous(trans='log10')

ggplot(GTEX_combined_targets, aes(x= Regulon, y= Mean) ,stat="identity") + 
  geom_boxplot(outlier.colour="black", outlier.shape = 20) +
  scale_y_continuous(trans='log10')

######################################################
#Survival analysis:
######################################################

library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condsurv)

ezfun::set_ccf_palette("contrast")
# Download TCGA-LUAD.survival.tsv
SurvivalFile <- read.delim("/PATH/TO/TCGA-LUAD.survival.tsv")

head(SurvivalFile[,c("OS.time","X_PATIENT","OS")], n=20)

#Survival times
Surv(SurvivalFile$OS.time,SurvivalFile$OS)[1:10]

s1 <- survfit(Surv(OS.time, OS) ~ 1, data = SurvivalFile)
str(s1)

##### Kaplan-Meier plots ####
# overall survival
survfit2(Surv(OS.time, OS) ~ 1, data = SurvivalFile) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + 
  add_confidence_interval()

# But first we need to match gene expression levels (master regulon and signature genes) to the data

# We're only interested in the CD8+ T-cell fractions
Cibersort_CD8T <- as.data.frame(res_ciber[,4])
LUAD_T <- setNames(data.frame(t(LUAD[,-1])), LUAD[,1])
                   
idx2 <- match(rownames(Cibersort_CD8T),row.names(LUAD_T))
LUAD_T$Cibersort <- Cibersort_CD8T$`res_ciber[, 4]` [idx2]

# Row means
LUAD_T$Mean <- rowMeans(LUAD_T)

# multiply with CD8 (enkel van de signature genes per regulon)
LUAD_T$MeanCib <- LUAD_T$Mean * LUAD_T$Cibersort

LUAD_short_T <- setNames(data.frame(t(LUAD_short[,-1])), LUAD_short[,1])
SurvivalFile$patient = gsub("-",".", SurvivalFile$X_PATIENT)

# Survival analysis for signature genes

#### gene
LUAD_gene_T <- LUAD_gene[,-c(517:519)]
LUAD_gene_T <- setNames(data.frame(t(LUAD_gene_T[,-1])), LUAD_gene[,1])
LUAD_gene_T$Mean <- rowMeans(LUAD_gene_T)

idx2 <- match(rownames(Cibersort_CD8T),row.names(LUAD_gene_T))
#LUAD_gene_T$Cibersort <- Cibersort_CD8T$`res_ciber[, 4]` [idx2]
LUAD_gene_T$Cibersort <- Cibersort_CD8T$ranged_gene [idx2]
LUAD_gene_T$MeanCib <- LUAD_gene_T$Mean * LUAD_gene_T$Cibersort
LUAD_gene_T <- LUAD_gene_T[LUAD_gene_T$Cibersort !=0,]

LUAD_gene_T_top_bottom <- LUAD_gene_T[order(LUAD_gene_T$MeanCib, decreasing = TRUE), ]

LUAD_gene_T_top <- as.vector(head(rownames(LUAD_gene_T_top_bottom), n=100))
LUAD_gene_T_bottom <- as.vector(tail(rownames(LUAD_gene_T_top_bottom), n=100))

LUAD_gene_T_top_bottom$level <- ifelse(rownames(LUAD_gene_T_top_bottom) %in% LUAD_gene_T_top, "high",
                                       ifelse(rownames(LUAD_gene_T_top_bottom) %in% LUAD_gene_T_bottom, "low",
                                              "mid"))

idx1 <- match(SurvivalFile$patient, rownames(LUAD_gene_T_top_bottom)) 
SurvivalFile$gene_S <- LUAD_gene_T_top_bottom$level [idx1]

# comparing survival times between groups
survdiff(Surv(OS.time, OS) ~ gene_S, data = SurvivalFile)

coxph(Surv(OS.time, OS) ~ gene_S, data = SurvivalFile)

coxph(
  Surv(OS.time, OS) ~ gene_S, 
  #  subset = gene_S == c("high","low"), 
  data = SurvivalFile
) %>% 
  tbl_regression(exp = TRUE)

survfit2(Surv(OS.time, OS) ~ gene_S, data = SurvivalFile) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + 
  add_confidence_interval()
