library(ggplot2)
library("lattice")
library(viridis)
library(dplyr)


setwd("/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/data/")
data <- read.table("out.geno.ld", header = TRUE, sep = "\t")

data$distance <- data$POS2 - data$POS1

data1 <- data[1:10000,]


p <- ggplot(data1, aes(x=distance, y=R.2)) + 
  geom_point(alpha = 0.02) +
  geom_smooth(method = 'loess', se = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(aspect.ratio=1,
        legend.position="none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),)


#ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/ld.pdf', plot = p, width = 15, height = 12, unit = 'cm')
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/ld_loess.png', plot = p, width = 15, height = 12, unit = 'cm')


#===scatterSmooth===#

pdf(file="ld_smooth.pdf",width=15,height=12)
p <- smoothScatter(data$distance, data$R.2)
dev.off()

#===stat_binhex===#

#http://auguga.blogspot.com/2016/09/r-log-transformed-color-density-scale.html
#myColor <- RColorBrewer::brewer.pal(9, "YlGnBu")
myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
#myColor <- RColorBrewer::brewer.pal(11, "BrBG")

myColor_scale_fill <- scale_fill_gradientn(colours = myColor, trans='log10')


p <- ggplot(data, aes(x=distance, y=R.2)) + 
  myColor_scale_fill + 
  stat_binhex(bins=100) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(aspect.ratio=1,
        legend.position = "none", 
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),)

#ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/ld_stat_binhex5.png', plot = p, width = 15, height = 12, unit = 'cm')
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/ld_stat_binhex5_classic.pdf', plot = p, width = 15, height = 12, unit = 'cm')

