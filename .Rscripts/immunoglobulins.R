library(ggplot2)
library(tidyverse)
library(viridis)
library(viridisLite)
library(readxl)
library(ggpubr)
Immunoglobulins<- read_excel("Immunoglobulins.xlsx")
head(Immunoglobulins)
Immunoglobulins$replicate <- as.factor(Immunoglobulins$replicate)
head(Immunoglobulins)
my_comparisons <- list(c("control", "UGAM"),c("Normal","Disease"))
Immunoglobulins$treatment <- factor(Immunoglobulins$treatment, 
                             levels=c("control", "UGAM", "Normal", "Disease"))
p_antibody <- ggplot(Immunoglobulins, aes(x=treatment, y=concentrate, fill=treatment)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size =1) +
  geom_jitter(shape=16, size=1, alpha =0.8,position = position_jitter(0.18)) +
  stat_compare_means(comparisons = my_comparisons, size=2, exact = FALSE) + # add pairwise comparison
  #     stat_compare_means()  + # add global p-value +
  labs(
    #         title = expression(paste("Antibodies influenced by AMUP and DSS ",italic("Astragalus membranaceus")," ultrafine particles")),
    #         subtitle = "(limited to characters with more than 100 appearances)",
    y = "immunoglobulins £¨¦Ìg/ml£©", x = "") +
  theme(axis.text.x = element_text(angle=45))+
  facet_wrap(~items, scales= "free", nrow=1)
p_Immunoglobulins