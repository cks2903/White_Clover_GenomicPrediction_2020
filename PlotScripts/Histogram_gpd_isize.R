library(ggplot2)
library("lattice")


setwd("/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/data/")
#old data: data_rhiz <- read.table("gpd_isize_hist.csv", header = TRUE, sep = ";")
data_rhiz <- read.table("gpd_ResCorfulld6.csv", header = TRUE, sep = ";")
head(data_rhiz)


#================================================#
#                     Hist GPD                   #
#================================================#
# Histograms for gpd for individual round replicates 
{
 
  data11 <- subset(data_rhiz, round ==1 & replicate ==1)
  data12 <- subset(data_rhiz, round ==1 & replicate ==2)
  data21 <- subset(data_rhiz, round ==2 & replicate ==1)
  data22 <- subset(data_rhiz, round ==2 & replicate ==2)
  
  # X axis limits 
  x <- max(data_rhiz$growth_per_day, na.rm = T)
  
  p11 <- ggplot(data11, aes(growth_per_day))+
    geom_histogram(fill="gray44") +
    theme_classic() +
    scale_x_continuous(limits = c(0, x), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  p11
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/gpd11.pdf', plot = p11, width = 15, height = 12, unit = 'cm')
  
  p12 <- ggplot(data12, aes(growth_per_day))+
    geom_histogram(fill="gray44") +
    theme_classic() +
    scale_x_continuous(limits = c(0, x), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  p12
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/gpd12.pdf', plot = p12,  width = 15, height = 12, unit = 'cm')
  
  p21 <- ggplot(data21, aes(growth_per_day))+
    geom_histogram(fill="gray44") +
    theme_classic() +
    scale_x_continuous(limits = c(0, x), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  p21
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/gpd21.pdf', plot = p21,  width = 15, height = 12, unit = 'cm')
  
  p22 <- ggplot(data22, aes(growth_per_day))+
    geom_histogram(fill="gray44") +
    theme_classic() +
    scale_x_continuous(limits = c(0, x), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  p22
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/gpd22.pdf', plot = p22,  width = 15, height = 12, unit = 'cm')
  
}




#================================================#
#                   Hist iSize                   #
#================================================#
# Histograms for iSize for individual round replicates 
{
  
  data11 <- subset(data_rhiz, round ==1 & replicate ==1)
  data12 <- subset(data_rhiz, round ==1 & replicate ==2)
  data21 <- subset(data_rhiz, round ==2 & replicate ==1)
  data22 <- subset(data_rhiz, round ==2 & replicate ==2)
  
  # X axis limits 
  x <- max(data_rhiz$InitialSize, na.rm = T)
  
  p11 <- ggplot(data11, aes(InitialSize))+
    geom_histogram(fill="gray77") +
    theme_classic() +
    scale_x_continuous(limits = c(0, x), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 200), expand = c(0, 0)) +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  p11
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/iSize11.pdf', plot = p11, width = 15, height = 12, unit = 'cm')
  
  p12 <- ggplot(data12, aes(InitialSize))+
    geom_histogram(fill="gray77") +
    theme_classic() +
    scale_x_continuous(limits = c(0, x), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 200), expand = c(0, 0)) +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  p12
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/iSize12.pdf', plot = p12,  width = 15, height = 12, unit = 'cm')
  
  p21 <- ggplot(data21, aes(InitialSize))+
    geom_histogram(fill="gray77") +
    theme_classic() +
    scale_x_continuous(limits = c(0, x), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 200), expand = c(0, 0)) +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  p21
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/iSize21.pdf', plot = p21,  width = 15, height = 12, unit = 'cm')
  
  p22 <- ggplot(data22, aes(InitialSize))+
    geom_histogram(fill="gray77") +
    theme_classic() +
    scale_x_continuous(limits = c(0, x), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 200), expand = c(0, 0)) +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  p22
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/iSize22.pdf', plot = p22,  width = 15, height = 12, unit = 'cm')
  
}