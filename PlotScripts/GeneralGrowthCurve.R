library(ggplot2)

{
  setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/New_fixation_Trait_20210106")
  library("lme4")
  #library("BGLR")
  library("BayzR")
  library("parallel")
  d=read.table("greenhouse_data_normalized_24092019.csv",header=T,sep=";",stringsAsFactors = FALSE)
  dim(d)
}

d=d[-which(d$Rhizobium=="NO"),]

plot <- ggplot(d, aes(NormTime,Size)) + geom_point(alpha=0.01) +
  scale_x_continuous(name="dpi", seq(0,70,10))+
  theme_classic()
p=plot + geom_smooth(method = "gam",se = FALSE, color="green")

ggsave(paste('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GrowthPlot',Sys.Date(),'.png',sep="_"), plot = p, width = 20, height = 10, unit = 'cm')
ggsave(paste('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GrowthPlot',Sys.Date(),'.pdf',sep="_"), plot = p, width = 20, height = 10, unit = 'cm')


# then for only unioculated plants
{
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/New_fixation_Trait_20210106")
library("lme4")
#library("BGLR")
library("BayzR")
library("parallel")
d=read.table("greenhouse_data_normalized_24092019.csv",header=T,sep=";",stringsAsFactors = FALSE)
dim(d)
}

d=d[which(d$Rhizobium=="NO"),]

plot <- ggplot(d, aes(NormTime,Size)) + geom_point(alpha=0.08) +
  scale_x_continuous(name="dpi", seq(0,70,10)) +
theme_classic() 
p=plot + geom_smooth(method = "gam",se = FALSE, color="green")

ggsave(paste('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GrowthPlot_UninoculatedPlants',Sys.Date(),'.png',sep="_"), plot = p, width = 20, height = 10, unit = 'cm')
ggsave(paste('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GrowthPlot_UninoculatedPlants',Sys.Date(),'.pdf',sep="_"), plot = p, width = 20, height = 10, unit = 'cm')

