#read gene expression matrix 
input <- "/PATH/TO/Luad.txt"

sig_matrix <- "/PATH/TO/LM22.txt"

#Run CIBERSORT abs 
#The number of permutation
cibersort_perm = 100
#Quantile normalization of input mixture, default = FALSE for RNA-Seq data
cibersort_qn = FALSE
#whether to apply absolute mode in cibersort
cibersort_abs = TRUE
#sig.score = for each mixture sample, define S as the median expression,level of all genes in the signature matrix divided by the median expression level of all genes in the mixture. Multiple cell subset fractions by S.
cibersort_abs_method = "sig.score"
res_ciber <- CIBERSORT(sig_matrix, 
                       input, 
                       perm = cibersort_perm, 
                       QN = cibersort_qn,
                       saveLocation = "/Users/lindadmin/Desktop/Impactys/test_cybersort")

head(res_ciber,3)

# We're only interested in the CD8+ T-cell fractions
Cibersort_CD8T <- as.data.frame(res_ciber[,4])
Cibersort_CD8T$scaled <- scale(Cibersort_CD8T$`res_ciber[, 4]`, center = TRUE, scale = TRUE)

LUAD_T <- setNames(data.frame(t(LUAD[,-1])), LUAD[,1])

idx2 <- match(rownames(Cibersort_CD8T),row.names(LUAD_T))
LUAD_T$Cibersort <- Cibersort_CD8T$`res_ciber[, 4]` [idx2]

LUAD_Cibersort <- LUAD_T %>% across(funs(.*Cibersort))
LUAD_Cibersort_short <- LUAD_Cibersort[colnames(LUAD_Cibersort) %in% ListRegulons]
write.csv(LUAD_Cibersort_short,"LUAD_after_cibersort.csv")
