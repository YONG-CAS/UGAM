library(heatmaply)
library(ggplot2)
library(tidyverse)
library(plotly)
library(viridis)
library(viridisLite)
library(readxl)
library(ggpubr)
library(dplyr)
library(Cairo)
library(hrbrthemes)
library(rlist)
library(orca)
library(gridExtra)
library(ggtree)
library(magrittr)
library(tidyr)
library(d3heatmap)
# library(aplot)
fulldata <- read_excel("fulldataset.xls", sheet =  "periphery")
head(fulldata)
dataset <- fulldata[, 6:21]
head(dataset)
cor_data <- round(cor(dataset, use = "complete.obs", method = "pearson"), 3)
svg("cor_map2.svg", height= 8, width= 8)
heatmap(as.matrix(cor_data))
dev.off()
library(corrplot)
svg("cor_map.svg", height= 8, width= 8)
corrplot(cor(cor_data),
         method ="number",
         type = "upper")
# ggpairs(cor_data)
# ggsave("cor_map2.svg", width = 8, height = 8, units = "cm")
dev.off()
createCorrelationPlots <- function(path, filename, df) {
  corrPlot <- ggpairs(
    df, diag=list(continuous="density"), axisLabels='show')
  png(file.path(path, filename), height=1000, width=1000)
  print(corrPlot)
  dev.off()
}
library(ComplexHeatmap)
library(circlize)
res_list <- readRDS("meth.rds")
summary(res_list)
type <- res_list$type
mat_meth <- res_list$mat_meth
mat_expr <- res_list$mat_expr
direction <- res_list$direction
cor_pvalue <- res_list$cor_pvalue
gene_type <- res_list$gene_type
anno_gene <- res_list$anno_gene
dist <- res_list$dis
anno_enhancer <- res_list$anno_enhancer
##首先计算甲基化矩阵的列聚类，以便可以将表达矩阵中的列调整为具有与甲基化矩阵中相同的列顺序。
column_tree <- hclust(dist(t(mat_meth)))
column_order <- column_tree$order
library(RColorBrewer)
#定义甲基化表达水平颜色，从0/blue-0.5/white-1/red渐变
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
# definition the direction of color
direction_col <- c("hyper" = "red", "hypo"= "blue")
# definition gene expression level
expr_col_fun <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
pvalue_col_fun <- colorRamp2(c(0,2,4), c("white", "white", "red"))
gene_type_col = structure(brewer.pal(length(unique(gene_type)), "Set3"), 
                          names = unique(gene_type))
anno_gene_col = structure(brewer.pal(length(unique(anno_gene)), "Set1"), 
                          names = unique(anno_gene))
#定义距离颜色
dist_col_fun = colorRamp2(c(0, 10000), c("black", "white"))
#定义增强子相关颜色
enhancer_col_fun = colorRamp2(c(0, 1), c("white", "orange"))
names(ht_global_opt())   
ht_global_opt(
  heatmap_legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
  heatmap_legend_labels_gp = gpar(fontsize = 8), 
  heatmap_column_names_gp = gpar(fontsize = 8),
  heatmap_column_title_gp = gpar(fontsize = 10),
  heatmap_row_title_gp = gpar(fontsize = 8),ADD=TRUE
)
ha = HeatmapAnnotation(type = type, 
                       col = list(type = c("Tumor" = "pink", "Control" = "royalblue")),
                       annotation_name_side = "left")
ha2 = HeatmapAnnotation(type = type, 
                        col = list(type = c("Tumor" = "pink", "Control" = "royalblue")), 
                        show_legend = FALSE)
ht_list = Heatmap(mat_meth, name = "methylation", col = meth_col_fun,
                  column_order= column_order,
                  top_annotation = ha, column_title = "Methylation") +
  Heatmap(direction, name = "direction", col = direction_col) +
  Heatmap(mat_expr[, column_tree$order], name = "expression", 
          col = expr_col_fun, 
          column_order = column_order, 
          top_annotation = ha2, column_title = "Expression") +
  Heatmap(cor_pvalue, name = "-log10(cor_p)", col = pvalue_col_fun)+ 
  Heatmap(gene_type, name = "gene type", col = gene_type_col) +
  Heatmap(anno_gene, name = "anno_gene", col = anno_gene_col) +
  Heatmap(dist, name = "dist_tss", col = dist_col_fun) +
  Heatmap(anno_enhancer, name = "anno_enhancer", col = enhancer_col_fun, 
          cluster_columns = FALSE, column_title = "Enhancer")
draw(ht_list, km = 2, split = direction,
     column_title = "Comprehensive correspondence between methylation,
expression and other genomic features", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
     merge_legends = TRUE, heatmap_legend_side = "bottom")
ht_global_opt(RESET = TRUE)