

library(ggplot2)
library(agricolae)

setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_F1_20210222")

# gpd

gpdF1 = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_F1_20210222/gpd/Correlations_GBLUP_GPD.txt",sep="\t",head=T)
dim(gpdF1)
p = ggplot(gpdF1,aes(Mean.Parental.predicted.GEBV,F1.pop.mean.dryweight))+
  geom_point(data=gpdF1, aes(y= F1.pop.mean.dryweight, x= Mean.Parental.predicted.GEBV, color=F1.population, size=1.5)) +
  scale_color_manual(values=c("#242c5e","#426468","#443659","#580b1d","#98d3ce","#71a188","#b01f35","#f6c29d","#f48988"))+
  theme_classic() +
  geom_smooth(method = "lm", se=F, colour="black") +
  annotate(geom = "text", x = 0.01, y = 3.0, label="italic(R) == 0.95  (p < 0.001)", parse=TRUE,hjust = "left")



correlation(gpdF1$Mean.Parental.predicted.GEBV,gpdF1$F1.pop.mean.dryweight)
correlation(gpdF1$Mean.Parental.predicted.GEBV,gpdF1$F1.pop.mean.freshweight)

ggsave(paste('F1_gpd_',Sys.Date(),'.pdf',sep=""), plot = p, width = 23, height = 15, unit = 'cm')



# isize

iSizeF1 = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_F1_20210222/iSize/Correlations_GBLUP_iSize.txt",sep="\t",head=T)
dim(iSizeF1)
p = ggplot(iSizeF1,aes(Mean.Parental.predicted.GEBV,F1.pop.mean.dryweight))+
  geom_point(data=iSizeF1, aes(y= F1.pop.mean.dryweight, x= Mean.Parental.predicted.GEBV, color=F1.population, size=1.5)) +
  scale_color_manual(values=c("#242c5e","#426468","#443659","#580b1d","#98d3ce","#71a188","#b01f35","#f6c29d","#f48988"))+
  theme_classic() +
  geom_smooth(method = "lm", se=F, colour="black") +
  annotate(geom = "text", x = 800, y = 3.0, label="italic(R) == 0.94  (p < 0.001)", parse=TRUE,hjust = "left")



correlation(iSizeF1$Mean.Parental.predicted.GEBV,iSizeF1$F1.pop.mean.dryweight)
correlation(iSizeF1$Mean.Parental.predicted.GEBV,iSizeF1$F1.pop.mean.freshweight)

ggsave(paste('F1_iSize_',Sys.Date(),'.pdf',sep=""), plot = p, width = 23, height = 15, unit = 'cm')



# phenotypic correlations, gpd

phenotypesParents = read.table("PhenotypicSelection.csv", head =T, sep=",")
phenotypesParents=as.data.frame(phenotypesParents)

p = ggplot(phenotypesParents,aes(Mean.Parental.gpd,F1.pop.mean.dryweight))+
  geom_point(data=gpdF1, aes(y= phenotypesParents$F1.pop.mean.dryweight, x= phenotypesParents$Mean.Parental.gpd, color=F1.population, size=1.5)) +
  scale_color_manual(values=c("#242c5e","#426468","#443659","#580b1d","#98d3ce","#71a188","#b01f35","#f6c29d","#f48988"))+
  theme_classic() +
  geom_smooth(method = "lm", se=F, colour="black") +
  annotate(geom = "text", x = 0.53, y = 3.0, label="italic(R) == 0.92  (p < 0.001)", parse=TRUE,hjust = "left")



correlation(phenotypesParents$Mean.Parental.gpd,gpdF1$F1.pop.mean.dryweight)
correlation(phenotypesParents$Mean.Parental.gpd,gpdF1$F1.pop.mean.freshweight)

ggsave(paste('F1_gpd_phenotypicsel',Sys.Date(),'.pdf',sep=""), plot = p, width = 23, height = 15, unit = 'cm')



# phenotypic correlations, iSize

p = ggplot(phenotypesParents,aes(Mean.Parental.iSize,F1.pop.mean.dryweight))+
  geom_point(data=gpdF1, aes(y= phenotypesParents$F1.pop.mean.dryweight, x= phenotypesParents$Mean.Parental.iSize, color=F1.population, size=1.5)) +
  scale_color_manual(values=c("#242c5e","#426468","#443659","#580b1d","#98d3ce","#71a188","#b01f35","#f6c29d","#f48988"))+
  theme_classic() +
  geom_smooth(method = "lm", se=F, colour="black") +
  annotate(geom = "text", x = 20000, y = 3.0, label="italic(R) == 0.78  (p < 0.05)", parse=TRUE,hjust = "left")



correlation(phenotypesParents$Mean.Parental.iSize,gpdF1$F1.pop.mean.dryweight)
correlation(phenotypesParents$Mean.Parental.iSize,gpdF1$F1.pop.mean.freshweight)

ggsave(paste('F1_iSize_phenotypicsel',Sys.Date(),'.pdf',sep=""), plot = p, width = 23, height = 15, unit = 'cm')


