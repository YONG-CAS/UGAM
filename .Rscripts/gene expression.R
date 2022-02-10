library(ggplot2)
library(tidyverse)
library(viridis)
library(viridisLite)
library(readxl)
library(ggpubr)
gene_expression <- read_excel("data.xlsx", sheet = "gene")
head(gene_expression)
gene_expression$treatment <- factor(gene_expression$treatment, 
                                    levels=c("control", "UGAM", "Normal", "Disease"))

my_comparisons <- list( c("control", "UGAM"),c("Normal","Disease"))
p_gene_expression <- ggplot(gene_expression, aes(x=treatment, y=concentrate, fill=treatment)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size =1) +
  geom_jitter(shape=16, size=1.0, alpha = 0.8, position = position_jitter(0.18)) +
  stat_compare_means(comparisons = my_comparisons, size=2, exact = FALSE) + # add pairwise comparison
  #     stat_compare_means()  + # add global p-value +
  labs(title = "Gene expression of complement proteins",
       #         subtitle = "(limited to characters with more than 100 appearances)",
       y = "Gene expressions ", x = "")+ 
  theme(axis.text.x = element_text(angle=45))+
  facet_wrap(~index, scales="free", nrow=2)
p_gene_expression
svg("gene_expression.svg")
p_gene_expression
dev.off()