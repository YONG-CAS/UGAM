library(ggplot2)
library(tidyverse)
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
library(pheatmap)
complement_full <- read_excel("data.xlsx",sheet="complement ")
head(complement_full)
my_comparisons <- list( c("control", "UGAM"),c("Normal","Disease"))
complement_full$treatment <- factor(complement_full$treatment, 
                                    levels=c("control", "UGAM", "Normal", "Disease"))
_complement <- ggplot(complement_full, aes(x=treatment, y=concentrate, fill=treatment)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1.2) +
  geom_jitter(shape=16, size=1.0, position = position_jitter(0.25), alpha =0.8) +
  stat_compare_means(comparisons = my_comparisons, size=2, exact = TRUE) + # add pairwise comparison
  labs(y = "Concentrate£¨mg/g£©", x = "Complement structural proteins and regulators")+
  theme(axis.text.x = element_text(angle=45))+
  facet_wrap(~index, scales="free", nrow=2)
svg("complement.svg")
p_complement
dev.off()

p_complement <- ggplot(complement_full, aes(x=treatment, y=concentrate, fill=treatment)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1.2) +
  geom_jitter(shape=16, size=1.5, position = position_jitter(0.25), alpha =0.58) +
  stat_compare_means(comparisons = my_comparisons, size=2, exact = TRUE) + # add pairwise comparison
  labs(y = "Concentrate£¨mg/g£©", x = "Complement structural proteins and regulators")+
  theme(axis.text.x = element_text(angle=45))+
  facet_wrap(~index,  scales="free", nrow=2)
svg("complement111.svg")
p_complement
dev.off()

complement <- read_excel("complement.xlsx")
head(complement)
complement$replicate <- as.factor(complement$replicate)
complement$daily_age <- as.factor(complement$daily_age)
head(complement)
complement$treatment <- factor(complement$treatment, 
                               levels=c("control", "UGAM", "Normal", "Disease"))
complement_01 <- transform(complement, treatment=NULL)
head(complement_01)
library(heatmaply)
complement_index <- complement_01[,c(6:12)]
head(complement_index)
pheatmap_complement <- pheatmap(complement_index
                                , scale="row"
                                , show_rownames = F
                                , cluster_rows = T
                                #          ,border_color = "grey60"
                                , cluster_distance_rows = "correlation" 
                                , clustering_distance_cols = "manhattan"
                                ,color=colorRampPalette(c("green", "black","red"))(50)
                                #                                 , col = RColorBrewer::brewer.pal(9, "RdBu"),
                                , column_names_gp = grid::gpar(fontsize=60)
)
print(pheatmap_complement)
save_pheatmap_svg <- function(x, filename, width=8, height=8){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  svg(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_svg(pheatmap_complement,"pheatmap4444.svg")



library(grid)
library(lattice)
ggarrange(p_complement, save_pheatmap_svg, nrow=1)
library(ComplexHeatmap)
library(grid)
pheatmap(complement_index, 
         , show_rownames = F
         , cluster_distance_rows = "correlation" 
         #         , clustering_distance_cols = "manhattan"
         ,row_names_gp = grid::gpar(fontsize = 20, angle=135)
)
complement001 <- complement_01[,c(6, 7, 8, 3, 5)]
head(complement001)
ggpairs(complement001)
p1 <- ggplot(data = complement, mapping = aes(x=C1q, y=C1r), colour = treatment) +  # ` ` symbols is used to collect the column name with space
  geom_point(data = transform(complement, treatment=NULL), mapping = aes(x=C1q, y = C1r), colour = "grey65") + 
  geom_point(colour="blue") +    
  stat_smooth(method = "loess") +
  facet_wrap(vars(treatment))
p1
p1 <- ggplot(data = complement, mapping = aes(x=C4bp, y=`Factor H`), colour = treatment) +  # ` ` symbols is used to collect the column name with space
  geom_point(data = transform(complement, treatment=NULL), mapping = aes(x=C4bp, y = `Factor.H`), colour = "grey65") + 
  geom_point(colour="blue") +    
  stat_smooth(method = "lm") +
  facet_wrap(vars(treatment))
svg("correlation_c4bp_factorH.svg")
p1
dev.off()


library(ggcorrplot)
head(tips)
library(GGally)
data(tips, package = "reshape")
ggpairs(
  tips[,c(1, 3, 4, 2)],
  upper = list(continous = "identity", combo = "box_no_facet"),
  lower = list(continous = "points", combo = "dot_no_facet")
)
head(iris)
ggpairs(iris)
p1 + geom_point(data = transform(complement, treatment=NULL), mapping = aes(x=C4bp, y = `Factor.H`), colour = "grey65")
a <- c(1,3,6)
b <- c(2, 4, 7)
frame <- data.frame(a, b)
frame
frame$sum <- frame$a + frame$b
frame
frame <- transform(frame)
frame
frame <- transform(frame, b=c(4, 29, 0))
frame
frame <- transform(frame, sum = a+b, mean=(a+b)/2)
frame
base <- ggplot(mpg, aes(hwy, class)) +
  geom_count() + 
  scale_x_binned(n.breaks = 20)
base
library(ggpmisc)
library(ggpp)
my.formula <- y ~ x
complement$type <- factor(complement$type, levels=c("immune_enhancing", "disease_model"))
p1 <- ggplot(data = complement, mapping = aes(x=C4bp, y=`Factor H`, colour=treatment)) +  # ` ` symbols is used to collect the column name with space
  geom_point() + 
  #     stat_smooth(method = "loess", formula = my.formula) +
  facet_wrap(vars(complement$type))
p1
p1 <- ggplot(data = complement, mapping = aes(x=C3, y=C4)) +
  geom_point() + 
  stat_smooth(method = lm, formula = y ~ x, size = 1)
p1
p1 <- ggplot(data = complement, mapping = aes(x=C1q, y=C1s)) +
  geom_point() + 
  stat_smooth(method = lm)
p1
p1 <- ggplot(data = complement, mapping = aes(x=C1s, y=C1q, colour = C1r)) +
  geom_point() + 
  stat_smooth(method = lm) +
  #     scale_fill_viridis_c(guide="legend") 
  scale_size_continuous(range=c(1, 15))
p1
# visualize : Specify the comparisons you want
my_comparisons <- list( c("C", "D"), c("C", "B"), c("C", "A"),c("B", "A"),c("B","D"),c("A","D"))
p_s <- ggplot(complement, aes(x=group, y=C1s, fill = group)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) +
  geom_jitter(shape=16, size = 2.2, position=position_jitter(0.2)) +
  stat_compare_means(comparisons = my_comparisons, label.y=c(500,530,560,590,620,650)) + # add pairwise comparison
  stat_compare_means(label.y=45)  + # add global p-value
  labs(title = expression(paste("C1s value with ", italic("Astragalus membranaceus"), "ultrafine particle"))
       , size =24)
svg("c1s.svg")
print(p_s)
dev.off()

# Basic box plot for C3 different components
my_comparisons <- list(c("B", "A"),c("C", "B"), c("C", "D"), c("C", "A"), c("A","D"), c("B","D"))
p_c1q <- ggplot(complement, aes(x=group, y=C1q, fill = group)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) + 
  geom_jitter(shape=16, size = 2.5, position=position_jitter(0.2)) +
  stat_compare_means(comparisons = my_comparisons, label.y=c(730,760,790,820,850,880)) + # add pairwise comparison
  stat_compare_means(label.y=220) +  # add global p-value
  labs(title = expression(paste("C1q value with ", italic("Astragalus membranaceus"), "ultrafine particle"))
       , size =24)

svg("c1q.svg")
print(p_c1q)
dev.off()

# Basic box plot for C1r different components
my_comparisons <- list( c("B", "A"),c("C", "B"),c("C", "D"), c("C", "A"), c("B","D"), c("A", "D"))
p_r <- ggplot(complement, aes(x=group, y=C1r, fill = group)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) +
  geom_jitter(shape=16, size = 2.5,
              position=position_jitter(0.2)) +
  stat_compare_means(comparisons = my_comparisons, label.y=c(320,340,360,380,400,420)) + # add pairwise comparison
  stat_compare_means(label.y=45) +  # add global p-value
  labs(title = expression(paste("C1r value with ", italic("Astragalus membranaceus"), "ultrafine particle"))
       , size =24)

svg("c1r.svg")
print(p_r)
dev.off()

my_comparisons <- list( c("B", "A"),c("B","C"),c("C", "D"), c("A","C"), c("B", "D"), c("D", "A"))
p_c3 <- ggplot(complement, aes(x=group, y=C3, fill = group)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  geom_jitter(shape=16, size = 2.5, position=position_jitter(0.2)) +
  stat_compare_means(comparisons = my_comparisons, label.y=c(58,60,62,64,67,69)) + # add pairwise comparison
  stat_compare_means(label.y=20) +  # add global p-value
  labs(title = expression(paste("C3 within ", italic("Astragalus membranaceus"), "ultrafine particle"))
       , size =24)
svg("c3.svg")
p_c3
dev.off()


ggplot(mpg, aes(displ, hwy)) +
  geom_point(data = transform(mpg, class = NULL), colour = "grey85") +
  geom_point() +
  facet_wrap(vars(class))

ggplot(mpg, aes(displ, hwy)) +
  geom_point(data = transform(mpg, class = NULL), colour = "grey85") +
  geom_point() +
  facet_wrap(vars(class))

ggplot(economics_long, aes(date, value)) +
  geom_line() +
  facet_wrap(vars(variable), scales = "free_y", nrow = 2, strip.position = "top") +
  theme(strip.background = element_blank(), strip.placement = "outside")
# }