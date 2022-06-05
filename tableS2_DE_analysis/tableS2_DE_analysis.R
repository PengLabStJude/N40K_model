library(dplyr)
library(preprocessCore)
library(Biobase)
library(limma)
library(MASS)


###################### read protein quantification data ########################
# read protein quantification data from file 'table_s2.csv'
table_s2_quan_df <- read.csv("inputs/table_s2_original.csv", header = T)

# extract the protein information and intensity matrix 
table_s2_prot_info <- table_s2_quan_df[, 1:8]
table_s2_quan_matrix <- table_s2_quan_df[, 9:28] %>% as.matrix()

# use protein accession No. as row names of the matrix
rownames(table_s2_quan_matrix) <- table_s2_quan_df$Protein.Accession..

# log2 transformation of protein intensity
table_s2_quan_matrix_log <- log2(table_s2_quan_matrix + 1)

# quantile normalization
table_s2_quan_matrix_log_norm <- normalize.quantiles(table_s2_quan_matrix_log, copy = T)
rownames(table_s2_quan_matrix_log_norm) <- rownames(table_s2_quan_matrix_log)
colnames(table_s2_quan_matrix_log_norm) <- colnames(table_s2_quan_matrix_log)


############################## log2FC calculation ##############################
# calculate mean values of Control and AD groups, respectively
group_mean_ctrl <- rowMeans(table_s2_quan_matrix_log_norm[, 1:10])
group_mean_ad <- rowMeans(table_s2_quan_matrix_log_norm[, 11:20])

# calculate log2FC for each protein
log2FC <- data.frame(log2FC = group_mean_ad - group_mean_ctrl)
# combined with protein information
log2FC_anno <- cbind(table_s2_prot_info, log2FC)


############################## limma DE analysis ###############################
# create ExpressionSet object used for limma analysis
table_s2_quan_eSet <- ExpressionSet(assayData = table_s2_quan_matrix_log_norm)

# create design matrix (allocate samples into groups)
design <- model.matrix(~ 0 + factor(c(rep("Ctrl",10), rep("AD",10)), levels = c('Ctrl', 'AD'))) 
colnames(design) <- c("Ctrl", "AD")

# fit linear model for each protein
fit <- lmFit(table_s2_quan_eSet, design)

# construct the contrast matrix (AD vs. Ctrl)
cont.matrix <- makeContrasts(ADvsCtrl=AD-Ctrl, levels=design)

# compute contrasts from linear model fit
fit2 <- contrasts.fit(fit, cont.matrix)

# empirical Bayes statistics for differential expression
fit2 <- eBayes(fit2)

# extract differential expression results
table_s2_DE_results <- topTable(fit2, adjust="BH", number = dim(table_s2_quan_eSet)[1])

# add protein information to the DE results table
matched_prot_info <- table_s2_prot_info[match(rownames(table_s2_DE_results), table_s2_prot_info$Protein.Accession..), ]
table_s2_DE_results_anno <- cbind(matched_prot_info, table_s2_DE_results)


########################### calculate SD cutoff ################################
################### calculate SD for intra-group fold change ###################
lfcMean = 0
lfcSD = 0
n = 0
m = 0

# 1. For Ctrl group, calculate fold change of each pair of samples
colCtrl = grep("Ctrl", colnames(table_s2_quan_matrix_log_norm))
# generate possible combination of samples 
combMatrix = combn(length(colCtrl), 2)
for (i in 1:ncol(combMatrix)) {
  
  # calculate fold change between two samples for each protein
  lfc = table_s2_quan_matrix_log_norm[, colCtrl[combMatrix[1, i]]] - table_s2_quan_matrix_log_norm[, colCtrl[combMatrix[2, i]]]
  # remove NA and tails
  lfc = lfc[!is.na(lfc)]
  lfc = lfc[lfc > quantile(lfc, 0.1) & lfc < quantile(lfc, 0.9)]
  # fit a normal distribution for the log2FC
  fit = fitdistr(lfc, "normal")
  # add up mean and sd of each log2FC distribution
  lfcMean = lfcMean + fit$estimate[1]
  lfcSD = lfcSD + fit$estimate[2]
  
  # draw density plots for each log2FC distribution
  if (i == 1) {
    plot(density(lfc), xlim = c(-2, 2), ylim = c(0, 2.5))
  } else {
    lines(density(lfc), xlim = c(-2, 2), ylim = c(0, 2.5))
  }
  
  # number of total proteins used for fitting
  n = n + fit$n
  # number of times of comparison
  m = m + 1
}

# 2. For AD group, calculate fold change of each pair of samples
colAD = grep("AD", colnames(table_s2_quan_matrix_log_norm))
combMatrix = combn(length(colAD), 2)
for (i in 1:ncol(combMatrix)) {
  
  # calculate fold change between two samples for each protein
  lfc = table_s2_quan_matrix_log_norm[, colAD[combMatrix[1, i]]] - table_s2_quan_matrix_log_norm[, colAD[combMatrix[2, i]]]
  # remove NA and tails
  lfc = lfc[!is.na(lfc)]
  lfc = lfc[lfc > quantile(lfc, 0.1) & lfc < quantile(lfc, 0.9)]
  # fit a normal distribution for the log2FC
  fit = fitdistr(lfc, "normal")
  # add up mean and sd of each log2FC distribution
  lfcMean = lfcMean + fit$estimate[1]
  lfcSD = lfcSD + fit$estimate[2]
  
  # draw density plots for each log2FC distribution
  lines(density(lfc), xlim = c(-2, 2), ylim = c(0, 1.2), col = "red")
  
  # number of total proteins used for fitting
  n = n + fit$n
  # number of times of comparison
  m = m + 1

}

# calculate mean values for lfcMean and lfcSD
lfcMean = lfcMean / m   # 0.00527
lfcSD = lfcSD / m       # 0.292

# generate a theoretical normal distribution based on lfcMean and lfcSD
lines(density(rnorm(round(n / ncol(table_s2_quan_matrix_log)), lfcMean, lfcSD)), col = "blue", lwd = 5)

################# calculate SD for inter-group fold change #####################
# remove tails
lfc_AD_Ctrl = log2FC$log2FC
lfc_AD_Ctrl_rm_tails = lfc_AD_Ctrl[lfc_AD_Ctrl > quantile(lfc_AD_Ctrl, 0.1) & lfc_AD_Ctrl < quantile(lfc_AD_Ctrl, 0.9)]
# fit a normal distribution for the log2FC
fit = fitdistr(lfc_AD_Ctrl_rm_tails, "normal")
# add up mean and sd of each log2FC distribution
lfc_AD_Ctrl_Mean = fit$estimate[1]  # -0.00306
lfc_AD_Ctrl_SD = fit$estimate[2]    # 0.156


################## remove isoforms, contaminants and keratins ##################
# some proteins have multiple isoforms, keep only one of them based on the "Protein Group" information
proteinGroup_split <- stringr::str_split_fixed(table_s2_DE_results_anno$Protein.Group., "\\.", 2)
# extract the last 3 digits of Protein Group No.
table_s2_DE_results_anno$proteinGroup_2 <- proteinGroup_split[,2] %>% as.numeric()
# for each unique Gene Name, keep only one protein
## remove proteins with no Gene names
table_s2_DE_results_anno_filter <- table_s2_DE_results_anno[!is.na(table_s2_DE_results_anno$GN),]
table_s2_DE_results_anno_filter <- table_s2_DE_results_anno_filter %>% group_by(GN) %>% 
  slice_min(order_by = proteinGroup_2, n=1, with_ties = F)

# remove contaminants and Keratins
## remove contaminants
con_index <- stringr::str_detect(table_s2_DE_results_anno_filter$Protein.Accession.., "co\\|")
table_s2_DE_results_anno_filter <- table_s2_DE_results_anno_filter[!con_index, ]
## remove keratins
keratin_index <- stringr::str_detect(table_s2_DE_results_anno_filter$Protein.Description, "Keratin\\,")
table_s2_DE_results_anno_filter <- table_s2_DE_results_anno_filter[!keratin_index, ]


######### define significant DE proteins based on log2FC and FDR cutoff ########
targetDE <- table_s2_DE_results_anno_filter[table_s2_DE_results_anno_filter$logFC > 0.3 & 
                                              table_s2_DE_results_anno_filter$adj.P.Val < 0.05, ]


######################### output the final DE results ##########################
write.csv(table_s2_DE_results_anno_filter[,-15], file = "outputs/table_s2_DE_results_anno_filter.csv", row.names = F)
