# load libraries
library(ggplot2)

# Load in data
{
  setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/New_fixation_Trait_20210106")
  library("lme4")
  #library("BGLR")
  library("BayzR")
  library("parallel")
  d=read.table("greenhouse_data_normalized_24092019.csv",header=T,sep=";",stringsAsFactors = FALSE)
  dim(d)
}

# Apply the same filtering and corrections as in the Genomic Prediction of yield
{
  before=dim(d)[1] # number of observations before filtering
  d=d[-which(d$Rhizobium=="SM73"),]
  d=d[-which(d$Rhizobium=="NO"),]
  d=d[-which(d$n.stolons<4),]
  
  d$Clover[which(d$Clover=="AAran_0104")]="Aaran_0104"
  d$Clover[which(d$Clover=="AAran_0206")]="Aaran_0206"
  
  d$roundRep <- paste(d$Round, d$Replicate, sep='_')
  
  after=dim(d)[1]
  before-after
  print(paste("removed",before-after,"observations,",after, "remains"),sep="")
}

# subset for the barcodes of interest
d2942 = subset(d,Barcode=="2942")
d3081 = subset(d,Barcode=="3081")
d1911 = subset(d,Barcode=="1911")

p1 = ggplot(d2942) + geom_line(aes(x = NormTime, y = Size)) + theme_classic() + ylim(0,200000) + xlim(1,50)
p2 = ggplot(d3081) + geom_line(aes(x = NormTime, y = Size)) + theme_classic()  + ylim(0,200000) + xlim(1,50)
p3 = ggplot(d1911) + geom_line(aes(x = NormTime, y = Size)) + theme_classic()  + ylim(0,200000) + xlim(1,50)

ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/growthcurves_2942.pdf",plot= p1,width=6,height=5,unit="cm")
ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/growthcurves_3081.pdf",plot= p2,width=6,height=5,unit="cm")
ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/growthcurves_1911.pdf",plot= p3,width=6,height=5,unit="cm")


