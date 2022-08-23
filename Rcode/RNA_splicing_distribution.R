library(readxl)
library(limma)
library(tidyverse)
N40K_neuron_splic <- read_excel("./Mouse_RNA_splicing_analysis_noNA.xlsx",sheet = 4) %>% filter(logFC_vec_N40K > 2*0.24 & adj.P.Val_vec_N40K < 0.2) %>% mutate_at(., c(3:8), log10)

N40K_neuron_splic$N40K <- rowMeans(N40K_neuron_splic[,6:8])
N40K_neuron_splic$WT <- rowMeans(N40K_neuron_splic[,3:5])

df_m0 = data.frame(cbind(N40K_neuron_splic$N40K, N40K_neuron_splic$WT))
colnames(df_m0) <- c("N40K","WT")
df_m1 <- df_m0 %>% pivot_longer(c("N40K","WT"), names_to = "id", values_to = "value")

ggplot(df_m1, aes(x = value, color = id))+
  geom_density(kernel="gaussian", size = 0.5)+ 
  scale_color_manual(values = c("#A00000","black"),name = NULL)+
  theme_classic()+
  theme(legend.position = "none")+
  xlim(-4,0)+
  #scale_x_continuous(breaks = c(-2,-1,0))+
  labs(title=NULL,x =NULL, y = NULL)
ggsave("./distribution/N40K_neuron_splicing_DE_dis_v1.pdf")
