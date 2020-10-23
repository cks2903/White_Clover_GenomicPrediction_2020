library(ggplot2)
library("lattice")
library(agricolae)


setwd("/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/data/")
data <- read.csv2("20201021_correlation.csv", header = TRUE, sep = ";")
head(data)

#================================================#
#                 Correlations                   #
#================================================#
correlation(data$InitialSize, data$gpd_NoCor, method = "pearson") # rho = 0.7244938  , p-value = 0
correlation(data$InitialSize, data$gpd_ResCor, method = "pearson") # rho = 0.1863677   , p-value = 0
correlation(data$InitialSize, data$gpd_FixCor, method = "pearson") # rho = 0.03526336    , p-value = 0.08465404 
correlation(data$gpd_NoCor, data$gpd_ResCor, method = "pearson") # rho =  0.773945     , p-value = 0
correlation(data$gpd_NoCor, data$gpd_FixCor, method = "pearson") # rho =  0.7144007    , p-value = 0
correlation(data$gpd_ResCor, data$gpd_FixCor, method = "pearson") # rho =   0.932936     , p-value = 0


#================================================#
#          Initial size correlations             #
#================================================#
 
{
  
  p11 <- ggplot(data, aes(InitialSize, gpd_NoCor))+
    geom_point() +
    theme_classic() +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  p11
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/cor_InitialSize_gpdnoCor.pdf', 
         plot = p11, width = 15, height = 12, unit = 'cm')
 
  p12 <- ggplot(data, aes(InitialSize, gpd_ResCor))+
    geom_point() +
    theme_classic() +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  p12
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/cor_InitialSize_gpdResCor.pdf', 
         plot = p12, width = 15, height = 12, unit = 'cm') 
  
  p13 <- ggplot(data, aes(InitialSize, gpd_FixCor))+
    geom_point() +
    theme_classic() +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  p13
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/cor_InitialSize_gpdFixedCor.pdf', 
         plot = p13, width = 15, height = 12, unit = 'cm') 
}


#================================================#
#            GPD noCor correlations              #
#================================================#

{
  
  p21 <- ggplot(data, aes(gpd_NoCor, gpd_ResCor))+
    geom_point() +
    theme_classic() +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  p21
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/cor_gpdnoCor_ResCor.pdf', 
         plot = p21, width = 15, height = 12, unit = 'cm')
  
  p22 <- ggplot(data, aes(gpd_NoCor, gpd_FixCor))+
    geom_point() +
    theme_classic() +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  p22
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/cor_gpdnoCor_FixedCor.pdf', 
         plot = p22, width = 15, height = 12, unit = 'cm')
}


#================================================#
#            GPD ResCor correlations              #
#================================================#

{
  
  p31 <- ggplot(data, aes(gpd_ResCor, gpd_FixCor))+
    geom_point() +
    theme_classic() +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  p31
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/cor_ResCor_fixedCor.pdf', 
         plot = p31, width = 15, height = 12, unit = 'cm')
  
}

#================================================#
#                   HISTOGRAMS                   #
#================================================#
# Histograms for individual traits
{
  

  p4 <- ggplot(data, aes(InitialSize))+
    geom_histogram(fill="black") +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  p4
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/cor_hist_InitialSize.pdf', 
         plot = p4, width = 15, height = 12, unit = 'cm')
  
  p5 <- ggplot(data, aes(gpd_NoCor))+
    geom_histogram(fill="black") +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  p5
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/cor_hist_gpd_NoCor.pdf', 
         plot = p5, width = 15, height = 12, unit = 'cm')
  
  p6 <- ggplot(data, aes(gpd_ResCor))+
    geom_histogram(fill="black") +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  p6
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/cor_hist_gpd_ResCor.pdf', 
         plot = p6, width = 15, height = 12, unit = 'cm')
  
  p7 <- ggplot(data, aes(gpd_FixCor))+
    geom_histogram(fill="black") +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(aspect.ratio=1,
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  p7
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/cor_hist_gpd_FixCor.pdf', 
         plot = p7, width = 15, height = 12, unit = 'cm')
}

