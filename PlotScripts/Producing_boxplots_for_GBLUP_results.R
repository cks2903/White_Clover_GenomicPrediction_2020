llibrary("ggplot2")
library("wesanderson")

setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd")

# accuracies and force zero
# Boxplots to compare prediction of different iSize corrected traits
Prediction_accuracies=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/ResultsOverview_differentTraits_new.csv",head=T,sep=",")
Prediction_accuracies=as.data.frame(Prediction_accuracies)
head(Prediction_accuracies)
Prediction_accuracies_true=Prediction_accuracies

# force all negative values to be zero
forcetobezero=which(Prediction_accuracies$Accuracy<0)
Prediction_accuracies$Accuracy[forcetobezero]=as.numeric(0)

# Figure 4 in manuscript
#png(file="PredictionAccuraciesOfDifferent_gpdTraits_20200911.png",width=20,height=10,units="cm",res=300)
ggplot(Prediction_accuracies, aes(x=Trait, y=Accuracy,fill=Trait)) + 
  geom_boxplot(aes(fill = Trait))+
  geom_point(aes(fill = Trait), size = 2, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) +
  ylim(c(min(Prediction_accuracies$Accuracy),0.615)) +
  theme_classic() +theme(axis.text.x = element_text(size = 10, angle = 90,hjust=1),axis.text.y = element_text(size = 10), 
                         axis.title=element_text(size=15))
ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/PredictionAccuraciesOfDifferent_gpdTraits_20201118.pdf", width =30, height = 20, units = "cm",useDingbats=FALSE)
#dev.off()


# Figure 4 in manuscript
#png(file="PredictionAccuraciesOfDifferent_gpdTraits_20200911_realvalues.png",width=20,height=10,units="cm",res=300)
ggplot(Prediction_accuracies_true, aes(x=Trait, y=Accuracy,fill=Trait)) + 
  geom_boxplot(aes(fill = Trait))+
  geom_point(aes(fill = Trait), size = 2, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) +
  ylim(c(min(Prediction_accuracies_true$Accuracy),0.615)) +
  theme_classic() +theme(axis.text.x = element_text(size = 10, angle = 90,hjust=1),axis.text.y = element_text(size = 10), 
                         axis.title=element_text(size=15))
ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/PredictionAccuraciesOfDifferent_gpdTraits_20201118_realvalues.pdf", width =30, height = 20, units = "cm",useDingbats=FALSE)
#dev.off()

# Figure 4 in manuscript
#png(file="PredictionCorrelationsOfDifferent_gpdTraits_20200911.png",width=20,height=10,units="cm",res=300)

ggplot(Prediction_accuracies_true, aes(x=Trait, y=Correlation,fill=Trait)) + 
  geom_boxplot(aes(fill = Trait))+
  geom_point(aes(fill = Trait), size = 2, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) +
  ylim(c(min(Prediction_accuracies_true$Correlation),max(Prediction_accuracies_true$Correlation))) +
  theme_classic() +theme(axis.text.x = element_text(size = 10, angle = 90,hjust=1),axis.text.y = element_text(size = 10), 
                         axis.title=element_text(size=15))
ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/PredictionCorrelationsOfDifferent_gpdTraits_20201118.pdf", width =30, height = 20, units = "cm",useDingbats=FALSE)
