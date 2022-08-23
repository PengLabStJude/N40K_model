# load package
library(edgeR)
library(limma)
library(tidyverse)
library(readxl)

mouse_expr <- read.csv("raw_data/mouse_raw_count_v1.csv")
mouse_splicing <- read_excel("/raw_data/Comprehensive_Splicing_Deficiency_Score_mouse_v2.xlsx",sheet = 4) 

# select genes shared by both RNA expression dataset and RNA splicing dataset
mouse_expr_geneid <- mouse_expr %>% select(geneId, transcript)
mouse_splicing_geneid <- mouse_splicing %>% select(GeneID, GeneName)
geneid_common <- inner_join(mouse_expr_geneid, mouse_splicing_geneid, by=c("transcript"="GeneID")) %>% select(transcript)
geneid_genename_common <- inner_join(mouse_expr_geneid, mouse_splicing_geneid, by=c("transcript"="GeneID")) %>% select(transcript, GeneName)

mouse_expr_fi <- left_join(geneid_common, mouse_expr, by= "transcript")
mouse_splicing_fi <- left_join(geneid_common, mouse_splicing, by=c("transcript"="GeneID"))

mouse_data <- mouse_expr_fi[,c(3:50)]
row.names(mouse_data) <- mouse_expr_fi[,1]


# define group
groupLabels <- c(rep("WT_3m", 3), rep("N40K_3m", 3),rep("WT_6m", 3), rep("N40K_6m", 3), rep("N40K_12m", 3), rep("WT_12m", 3), rep("HPC_WT",3), rep("HPC_N40K",3), rep("HPC_FAD",3), rep("HPC_dTg",3), rep("CTX_WT",3),rep("CTX_N40K",3), rep("CTX_FAD",3), rep("CTX_dTg",3), rep("vector", 3), rep("N40K", 3))
TS <- factor(groupLabels, levels = c("WT_3m", "N40K_3m","WT_6m", "N40K_6m","N40K_12m","WT_12m", "HPC_WT", "HPC_N40K","HPC_FAD","HPC_dTg","CTX_WT","CTX_N40K","CTX_FAD","CTX_dTg","vector","N40K"))

design <- model.matrix(~ 0 + TS )
colnames(design) <- levels(TS)
dge <- DGEList(counts = mouse_data, group = groupLabels)


# filter out low expressed genes
cpm_cutoff=10
cutoff <- as.vector(cpm(cpm_cutoff,mean(dge$samples$lib.size) ) )
keep <- rowSums(cpm(dge) > cutoff) >= min(as.numeric(table(groupLabels)))
dge <- dge[keep, keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

v <- voom(dge, design, plot = F)
fit <- lmFit(v, design)
cont.matrix <- makeContrasts(N40K_3mvsWT_3m = (N40K_3m - WT_3m), N40K_6mvsWT_6m = (N40K_6m - WT_6m), N40K_12mvsWT_12m = (N40K_12m - WT_12m), HPC_N40KvsHPC_WT = (HPC_N40K - HPC_WT),HPC_FADvsHPC_WT = (HPC_FAD - HPC_WT), HPC_dTgvsHPC_WT = (HPC_dTg - HPC_WT),CTX_N40KvsCTX_WT = (CTX_N40K - CTX_WT), CTX_FADvsCTX_WT = (CTX_FAD - CTX_WT), CTX_dTgvsCTX_WT = (CTX_dTg - CTX_WT), N40Kvsvector = (N40K - vector), levels = design)
fitcon <- contrasts.fit(fit, cont.matrix)
fitcon <- eBayes(fitcon)

results1 <- topTable(fitcon, n = Inf, sort.by= "P", coef="N40K_3mvsWT_3m") %>% tibble::rownames_to_column(., "geneId") %>% rename(logFC_3m=logFC, AveExpr_3m=AveExpr, t_3m=t, P.Value_3m=P.Value, adj.P.Val_3m = adj.P.Val, B_3m=B)

#outFile1 <- paste0("N40K_3mvsWT_3m", "_diff.txt")
#write.table(results1, file = outFile1, sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)

results2 <- topTable(fitcon, n = Inf, sort.by="P", coef="N40K_6mvsWT_6m") %>% tibble::rownames_to_column(., "geneId") %>% rename(logFC_6m=logFC, AveExpr_6m=AveExpr, t_6m=t, P.Value_6m=P.Value, adj.P.Val_6m = adj.P.Val, B_6m=B)


results3 <- topTable(fitcon, n = Inf, sort.by="P", coef="N40K_12mvsWT_12m") %>% tibble::rownames_to_column(., "geneId") %>% rename(logFC_12m=logFC, AveExpr_12m=AveExpr, t_12m=t, P.Value_12m=P.Value, adj.P.Val_12m = adj.P.Val, B_12m=B)

results4 <- topTable(fitcon, n = Inf, sort.by="P", coef="HPC_N40KvsHPC_WT") %>% tibble::rownames_to_column(., "geneId") %>% rename(logFC_HPC_N40K=logFC, AveExpr_HPC_N40K=AveExpr, t_HPC_N40K=t, P.Value_HPC_N40K=P.Value, adj.P.Val_HPC_N40K = adj.P.Val, B_HPC_N40K=B)

results5 <- topTable(fitcon, n = Inf, sort.by="P", coef="HPC_FADvsHPC_WT") %>% tibble::rownames_to_column(., "geneId") %>% rename(logFC_HPC_FAD=logFC, AveExpr_HPC_FAD=AveExpr, t_HPC_FAD=t, P.Value_HPC_FAD=P.Value, adj.P.Val_HPC_FAD = adj.P.Val, B_HPC_FAD=B)

results6 <- topTable(fitcon, n = Inf, sort.by="P", coef="HPC_dTgvsHPC_WT") %>% tibble::rownames_to_column(., "geneId") %>% rename(logFC_HPC_dTg=logFC, AveExpr_HPC_dTg=AveExpr, t_HPC_dTg=t, P.Value_HPC_dTg=P.Value, adj.P.Val_HPC_dTg = adj.P.Val, B_HPC_dTg=B)

results7 <- topTable(fitcon, n = Inf, sort.by="P", coef="CTX_N40KvsCTX_WT") %>% tibble::rownames_to_column(., "geneId") %>% rename(logFC_CTX_N40K=logFC, AveExpr_CTX_N40K=AveExpr, t_CTX_N40K=t, P.Value_CTX_N40K=P.Value, adj.P.Val_CTX_N40K = adj.P.Val, B_CTX_N40K=B)

results8 <- topTable(fitcon, n = Inf, sort.by="P", coef="CTX_FADvsCTX_WT") %>% tibble::rownames_to_column(., "geneId") %>% rename(logFC_CTX_FAD=logFC, AveExpr_CTX_FAD=AveExpr, t_CTX_FAD=t, P.Value_CTX_FAD=P.Value, adj.P.Val_CTX_FAD = adj.P.Val, B_CTX_FAD=B)

results9 <- topTable(fitcon, n = Inf, sort.by="P", coef="CTX_dTgvsCTX_WT") %>% tibble::rownames_to_column(., "geneId") %>% rename(logFC_CTX_dTg=logFC, AveExpr_CTX_dTg=AveExpr, t_CTX_dTg=t, P.Value_CTX_dTg=P.Value, adj.P.Val_CTX_dTg = adj.P.Val, B_CTX_dTg=B)

results10 <- topTable(fitcon, n = Inf, sort.by="P", coef="N40Kvsvector") %>% tibble::rownames_to_column(., "geneId") %>% rename(logFC_vec=logFC, AveExpr_vec=AveExpr, t_vec=t, P.Value_vec=P.Value, adj.P.Val_vec = adj.P.Val, B_vec=B)

mouse_data1 <- data.frame(v$E) %>% rownames_to_column(.,var = "geneId")

mouse_results_all0 <- inner_join(mouse_data1, results1, by="geneId") %>% inner_join(., results2, by="geneId") %>% inner_join(., results3, by="geneId") %>% inner_join(., results4, by="geneId") %>% inner_join(., results5, by="geneId") %>% inner_join(., results6, by="geneId") %>% inner_join(., results7, by="geneId") %>% inner_join(., results8, by="geneId") %>% inner_join(., results9, by="geneId") %>% inner_join(., results10, by="geneId") 

mouse_results_all <- inner_join(geneid_genename_common, mouse_results_all0, by=c("transcript"="geneId"))

N40K_3m <- grep("3m", colnames(mouse_results_all))
N40K_12m <- grep("12m", colnames(mouse_results_all))
HPC <- grep("HPC", colnames(mouse_results_all))
CTX <- grep("CTX", colnames(mouse_results_all))
N40K_3m_12m <- mouse_results_all[,c(1,2,N40K_3m,N40K_12m)]
dTg_HPC <- mouse_results_all[,c(1,2,HPC)]
dTg_CTX <- mouse_results_all[,c(1,2,CTX)]
vector <- mouse_results_all[,c(1,2,44:49,104:109)]

library(openxlsx)
wb1 <- createWorkbook()
addWorksheet(wb1, "N40K_3m_12m")
addWorksheet(wb1, "dTg_HPC")
addWorksheet(wb1, "dTg_CTX")
addWorksheet(wb1, "vector")

writeDataTable(wb1,"N40K_3m_12m", x=N40K_3m_12m)
writeDataTable(wb1,"dTg_HPC", x=dTg_HPC )
writeDataTable(wb1,"dTg_CTX", x=dTg_CTX)
writeDataTable(wb1,"vector", x=vector)
saveWorkbook(wb1, "analysis_results/Mouse_RNA_expression_analysis_v1.xlsx", overwrite = TRUE)

## use original data
mouse_results_all1 <- inner_join(mouse_expr, results1, by= c("transcript"="geneId")) %>% inner_join(., results2, by=c("transcript"="geneId")) %>% inner_join(., results3, by=c("transcript"="geneId")) %>% inner_join(., results4, by=c("transcript"="geneId")) %>% inner_join(., results5, by=c("transcript"="geneId")) %>% inner_join(., results6, by=c("transcript"="geneId")) %>% inner_join(., results7, by=c("transcript"="geneId")) %>% inner_join(., results8, by=c("transcript"="geneId")) %>% inner_join(., results9, by=c("transcript"="geneId")) %>% inner_join(., results10, by=c("transcript"="geneId")) 
