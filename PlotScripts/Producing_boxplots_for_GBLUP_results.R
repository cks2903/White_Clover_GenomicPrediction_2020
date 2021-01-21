library("ggplot2")
library("wesanderson")

setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/")

# accuracies and force zero
# Boxplots to compare prediction of different iSize corrected traits
Prediction_correlations=read.table("Correlations_gpd_20210121.csv",head=T,sep=",",stringsAsFactors = F)
Prediction_correlations=as.data.frame(Prediction_correlations)
head(Prediction_correlations)
Prediction_correlations$Trait[which(Prediction_correlations$Trait=="Day11to25")]="gpi"
Prediction_correlations$Trait[which(Prediction_correlations$Trait=="Day11to25Cor")]="gpiCor"

Prediction_correlations_true=Prediction_correlations

# force all negative values to be zero
forcetobezero=which(Prediction_correlations$Correlation<0)
Prediction_correlations$Correlation[forcetobezero]=as.numeric(0)

# Figure 4 in manuscript
#png(file="PredictionAccuraciesOfDifferent_gpdTraits_20200911.png",width=20,height=10,units="cm",res=300)
ggplot(Prediction_correlations, aes(x=Trait, y=Correlation,fill=Trait)) + 
  geom_boxplot(aes(fill = Trait))+
  geom_point(aes(fill = Trait), size = 2, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values=c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[2],"brown4","cadetblue",wes_palette("Rushmore1")[4])) +
  ylim(c(min(Prediction_correlations$Correlation),0.615)) +
  theme_classic() +theme(axis.text.x = element_text(size = 10, angle = 90,hjust=1),axis.text.y = element_text(size = 10), 
                         axis.title=element_text(size=15))+
  scale_y_continuous(breaks = round(seq(-0.25, 0.60, by = 0.05),1))

ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/PredictionAccuraciesOfDifferent_gpdTraits_20201118.pdf", width =30, height = 20, units = "cm",useDingbats=FALSE)
#dev.off()


# Figure 4 in manuscript
#png(file="PredictionAccuraciesOfDifferent_gpdTraits_20200911_realvalues.png",width=20,height=10,units="cm",res=300)
ggplot(Prediction_correlations_true, aes(x=Trait, y=Correlation,fill=Trait)) + 
  geom_boxplot(aes(fill = Trait))+
  geom_point(aes(fill = Trait), size = 2, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values=c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[2],"brown4","cadetblue",wes_palette("Rushmore1")[4])) +
  ylim(c(min(Prediction_correlations_true$Correlation),0.615)) +
  theme_classic() +theme(axis.text.x = element_text(size = 10, angle = 90,hjust=1),axis.text.y = element_text(size = 10), 
                         axis.title=element_text(size=15)) +
scale_y_continuous(breaks = round(seq(-0.25, 0.60, by = 0.05),1))

ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/PredictionCorOfDifferent_gpdTraits_20210121_realvalues.pdf", width =30, height = 20, units = "cm",useDingbats=FALSE)
#dev.off()

# Figure 4 in manuscript
#png(file="PredictionAccuraciesOfDifferent_gpdTraits_20200911_realvalues.png",width=20,height=10,units="cm",res=300)
ggplot(Prediction_correlations_true, aes(x=Trait, y=Accuracy,fill=Trait)) + 
  geom_boxplot(aes(fill = Trait))+
  geom_point(aes(fill = Trait), size = 2, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values=c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[2],"brown4","cadetblue",wes_palette("Rushmore1")[4])) +
  ylim(c(min(Prediction_correlations_true$Accuracy),0.615)) +
  theme_classic() +theme(axis.text.x = element_text(size = 10, angle = 90,hjust=1),axis.text.y = element_text(size = 10), 
                         axis.title=element_text(size=15)) +
  scale_y_continuous(breaks = round(seq(-0.25, 0.60, by = 0.05),1))

ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/PredictionAccOfDifferent_gpdTraits_20210121_realvalues.pdf", width =30, height = 20, units = "cm",useDingbats=FALSE)
#dev.off()
