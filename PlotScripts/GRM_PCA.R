library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(ggfortify)
library(plyr)
library(dplyr)
library(ggpubr)


#https://slowkow.com/notes/pheatmap-tutorial/
setwd('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/data/')

data = read.csv("GRM_Clover_Fullfiltering_20200728.csv", sep = ',')

data1 <- as.matrix(data)
GRM <- data1/mean(diag(data1))

#==================================================================================================================#
#                                                                                                                  #
#                                                       GRM                                                        #
#                                                                                                                  #
#==================================================================================================================#

my_palette = rev(brewer.pal(9, "BuPu"))

pheatmap(GRM,
         cluster_rows = F, cluster_cols = F,
         show_colnames     = T,
         show_rownames     = T,
         border_color = NA,
         col = my_palette,
         width = 17, 
         height = 17,
         filename = "/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/GRM.pdf")

#LEGEND
pheatmap(GRM,
         cluster_rows = F, cluster_cols = F,
         legend = TRUE,
         show_colnames     = T,
         show_rownames     = T,
         border_color = NA,
         col = my_palette,
         width = 17, 
         height = 17,
         filename = "/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/GRM_legend.pdf")



# BREAKS
quantile_breaks <- function(xs, n = 15) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(GRM, n = 15)

pheatmap(GRM,
         cluster_rows = F, cluster_cols = F,
         show_colnames     = T,
         show_rownames     = T,
         border_color = NA,
         col = magma(length(mat_breaks) -1),
         breaks = mat_breaks,
         width = 17, 
         height = 17,
         filename = "/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/GRM2.pdf")

#==================================================================================================================#
#                                                                                                                  #
#                                                       PCA                                                        #
#                                                                                                                  #
#==================================================================================================================#


d.pca <- prcomp(data,
                center = TRUE,
                scale. = TRUE)

percentVar <- d.pca$sdev^2 / sum( d.pca$sdev^2 )
percentVar <- round(100 * percentVar)

# Add clover variety column 
clover.names <- colnames(GRM)
p2.data <- data.frame(do.call(rbind, strsplit(as.character(clover.names), "")))
clover.names <-  paste(p2.data$X1, p2.data$X2, p2.data$X3, p2.data$X4, p2.data$X5, sep='') 

d.pca.rot <- as.data.frame(d.pca$rotation)

########################
######### PC12 #########
########################
grouping.colours <- c("#fcb831", "#ffdc38", "#b01533", "#42192a", "#426467", 
                      "#fa2853", "#f98988", "#f6c29d", "#bcbd99", "#71a189", 
                      "#7a0177", "#238b45", "#a2acc4", "#5a0019", "#883860",
                      "#5ac9e1", "#98d4cf", "#01665e", "#ee7125", "#d53e4f")

(pcaplot <- ggplot(d.pca, aes(PC1, PC2, color = clover.names , fill = clover.names)) +
   coord_fixed() +
   geom_point(size=3, aes(color = clover.names)) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    scale_fill_manual(values = grouping.colours) +
    scale_colour_manual(values = grouping.colours) +
   stat_chull(geom = "polygon", alpha = 0.2) +
   theme_classic() +
   theme(aspect.ratio=1,
         legend.title = element_blank(),
         axis.text.x = element_text(size = 20),
         axis.text.y = element_text(size = 20),
         axis.title=element_text(size = 22),
         legend.text = element_text(size = 20),
         legend.background = element_rect(color = NA)))


ggsave('GRM_PC12.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')









