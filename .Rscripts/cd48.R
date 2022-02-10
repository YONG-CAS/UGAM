library(ggplot2)
library(tidyverse)
library(viridis)
library(viridisLite)
library(readxl)
library(ggpubr)
cd48 <- read_excel("data.xlsx", sheet = "FCM")
head(cd48)
cd48$treatment
my_comparisons <- list(c("control", "UGAM"),c("Normal","Disease"))
cd48$treatment <- factor(cd48$treatment, 
                         levels=c("control", "UGAM", "Normal", "Disease"))
p_cd <- ggplot(cd48,aes(x=treatment, y=concentrate,fill=treatment)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size =1) +
  geom_jitter(shape=16, size=1.5, position = position_jitter(0.18), alpha =0.8) +
  stat_compare_means(comparisons = my_comparisons, size=4, exact = FALSE) + 
  labs(title = "CD4/CD8 T cells",
       y = "CD4+/CD8+ T cell ratios", x = "") +
  theme(axis.text.x = element_text(angle=45))
p_cd
res <- 124
svg("CD4CD8 T cells.svg",width = 560/res, height = 1080/res)
p_cd
dev.off()