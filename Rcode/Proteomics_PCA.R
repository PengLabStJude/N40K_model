library(tidyverse)
library(readxl)
library(MKmisc)
library(openxlsx)
library(scales)

df <- read_excel("./results/N40K_b2_b3_final_v4.xlsx", sheet = 1) 
df_log2 <- log2(df[,c(15:33)])
df2_t <- data.frame(t(df_log2))

pca <- prcomp(df2_t, scale = TRUE)
pca_sum <- summary(pca)
scores = as.data.frame(pca$x)

xlimr = (round(max(scores[1])/10,0) + 1)*10
xliml = (round(min(scores[1])/10,0) - 1)*10
ylimu = (round(max(scores[2])/10,0) + 1)*10
ylimd = (round(min(scores[2])/10,0) - 1)*10
xlab = percent(pca_sum$importance[2], accuracy=0.01)
ylab = percent(pca_sum$importance[5], accuracy=0.01)

pdf("./results/PCA_all_protein(pt).pdf")
plot(scores[,1:2], col=c(rep("#0000FF",3), rep("#A00000",3), rep("#008000",2), rep("#000000",5),rep("#008000",2), rep("#A00000",2), rep("#0000FF",2)), cex=0.8, pch=19, xlim = c(xliml,xlimr), ylim=c(ylimd,ylimu), main="PCA", xlab = c("PC1", xlab), ylab = c("PC2", ylab))
#text(x=scores[,1], y=scores[,2], labels=c(rownames(scores)), cex= 0.7, pos=3)
dev.off()
