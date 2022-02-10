library(ggplot2)
library(tidyverse)
library(viridis)
library(viridisLite)
library(readxl)
library(ggpubr)
cytokines <- read_excel("data.xlsx", sheet = "inflammatory cytokines")
head(cytokines)
cytokines$replicate <- as.factor(cytokines$replicate)
head(cytokines)
my_comparisons <- list(c("control", "UGAM"),c("Normal","Disease"))
cytokines$treatment <- factor(cytokines$treatment, 
                              levels=c("control", "UGAM", "Normal", "Disease"))
p_cytokines <- ggplot(cytokines, aes(x=treatment, y=concentrate, fill=treatment)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size =1) +
  geom_jitter(shape=16, size=1, position = position_jitter(0.18), alpha = 0.8) +
  stat_compare_means(comparisons = my_comparisons, size=2, exact = FALSE) + # add pairwise comparison
  #     stat_compare_means()  + # add global p-value
  labs(title = expression(paste("Inflammatory cytokines by ",italic("Astragalus membranaceus"), " ultrafine particles")),
       #         subtitle = "(limited to characters with more than 100 appearances)",
       y = "Inflammatory cytokines£¨pg/ml£©", x = "") +
  theme(axis.text.x = element_text(angle=45)) +
  facet_wrap(~items, scales="free", nrow=2)
p_cytokines
svg("cytokines.svg")
p_cytokines
dev.off()