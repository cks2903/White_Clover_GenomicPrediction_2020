
library("ggplot2")

setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831")

NoCorData=read.table("replicate_reduction_gpdNoCor.txt",head=T)
dim(NoCorData)

ResCorData=read.table("replicate_reduction_gpdResCor.txt",head=T)
dim(ResCorData)


# Combined

ResCorData_means_corrected=aggregate(ResCorData$cor_CorrectedPheno_GEBV, list(ResCorData$Replicates), mean)
colnames(ResCorData_means_corrected)=c("Replicates","cor_CorrectedPheno_GEBV")

NoCorData_means_corrected=aggregate(NoCorData$cor_CorrectedPheno_GEBV, list(NoCorData$Replicates), mean)
colnames(NoCorData_means_corrected)=c("Replicates","cor_CorrectedPheno_GEBV")


# plot together
p1 = ggplot() + 
  geom_line(data = ResCorData_means_corrected, aes(x = Replicates, y = cor_CorrectedPheno_GEBV), color = "#0B775E") +
  geom_line(data = NoCorData_means_corrected, aes(x = Replicates, y = cor_CorrectedPheno_GEBV), color = "#35274A") +
  scale_x_reverse() +
  xlab('Replicates') +
  ylab('Correlation') +
  theme_bw()

p1

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/CorrectedPhenoCorrelation_gpdNoCorANDgpdResCor.pdf', plot = p1, width = 15, height = 12, unit = 'cm')

# with standard deviations
sd_rescor=aggregate(ResCorData$cor_gpdResCor_GEBV, list(ResCorData$Replicates), sd)
ResCorData_means_corrected$sd=sd_rescor[,2]
ResCorData_means_corrected$upper=ResCorData_means_corrected$cor+ResCorData_means_corrected$sd
ResCorData_means_corrected$lower=ResCorData_means_corrected$cor-ResCorData_means_corrected$sd

sd_nocor=aggregate(NoCorData$cor_gpdNoCor_GEBV, list(NoCorData$Replicates), sd)
NoCorData_means_corrected$sd=sd_nocor[,2]
NoCorData_means_corrected$upper=NoCorData_means_corrected$cor+NoCorData_means_corrected$sd
NoCorData_means_corrected$lower=NoCorData_means_corrected$cor-NoCorData_means_corrected$sd

together=rbind(NoCorData_means_corrected,ResCorData_means_corrected)
together$trait=c(rep("gpd_NoCor",9),rep("gpd_ResCor",9))




p1.5 = ggplot(together, aes(Replicates, cor_CorrectedPheno_GEBV, colour= trait)) + 
  geom_line(aes(group = trait)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  scale_x_reverse() +
  xlab('Replicates') +
  ylab('Correlation') +
  scale_color_manual(values = c("#35274A","#0B775E"))+
  theme_bw()
p1.5

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/Correlation_gpdNoCorANDgpdResCor_pred_SDs.pdf', plot = p1.5, width = 15, height = 12, unit = 'cm')





# Combined

ResCorData_means_gpdResCor=aggregate(ResCorData$cor_gpdResCor_GEBV, list(ResCorData$Replicates), mean)
colnames(ResCorData_means_gpdResCor)=c("Replicates","cor_gpd_ResCor_GEBV")

NoCorData_means_gpdNoCor=aggregate(NoCorData$cor_gpdNoCor_GEBV, list(NoCorData$Replicates), mean)
colnames(NoCorData_means_gpdNoCor)=c("Replicates","cor_gpd_NoCor_GEBV")


# plot together
p2 = ggplot() + 
  geom_line(data = ResCorData_means_gpdResCor, aes(x = Replicates, y = cor_gpd_ResCor_GEBV), color = "#0B775E") +
  geom_line(data = NoCorData_means_gpdNoCor, aes(x = Replicates, y = cor_gpd_NoCor_GEBV), color = "#35274A") +
  scale_x_reverse() +
  xlab('Replicates') +
  ylab('Correlation') +
  theme_bw()

p2

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/Correlation_gpdNoCorANDgpdResCor.pdf', plot = p2, width = 15, height = 12, unit = 'cm')

t.test(ResCorData$cor_CorrectedPheno_GEBV[which(ResCorData$Replicates==2)], NoCorData$cor_CorrectedPheno_GEBV[which(NoCorData$Replicates==2)], paired = FALSE, alternative = "greater") # N.S
t.test(ResCorData$cor_CorrectedPheno_GEBV[which(ResCorData$Replicates==3)], NoCorData$cor_CorrectedPheno_GEBV[which(NoCorData$Replicates==3)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_CorrectedPheno_GEBV[which(ResCorData$Replicates==4)], NoCorData$cor_CorrectedPheno_GEBV[which(NoCorData$Replicates==4)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_CorrectedPheno_GEBV[which(ResCorData$Replicates==5)], NoCorData$cor_CorrectedPheno_GEBV[which(NoCorData$Replicates==5)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_CorrectedPheno_GEBV[which(ResCorData$Replicates==6)], NoCorData$cor_CorrectedPheno_GEBV[which(NoCorData$Replicates==6)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_CorrectedPheno_GEBV[which(ResCorData$Replicates==7)], NoCorData$cor_CorrectedPheno_GEBV[which(NoCorData$Replicates==7)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_CorrectedPheno_GEBV[which(ResCorData$Replicates==8)], NoCorData$cor_CorrectedPheno_GEBV[which(NoCorData$Replicates==8)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_CorrectedPheno_GEBV[which(ResCorData$Replicates==9)], NoCorData$cor_CorrectedPheno_GEBV[which(NoCorData$Replicates==9)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_CorrectedPheno_GEBV[which(ResCorData$Replicates==10)], NoCorData$cor_CorrectedPheno_GEBV[which(NoCorData$Replicates==10)], paired = FALSE, alternative = "greater") # ***

t.test(ResCorData$cor_gpdResCor_GEBV[which(ResCorData$Replicates==2)], NoCorData$cor_gpdNoCor_GEBV[which(NoCorData$Replicates==2)], paired = FALSE, alternative = "greater") # N.S
t.test(ResCorData$cor_gpdResCor_GEBV[which(ResCorData$Replicates==3)], NoCorData$cor_gpdNoCor_GEBV[which(NoCorData$Replicates==3)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_gpdResCor_GEBV[which(ResCorData$Replicates==4)], NoCorData$cor_gpdNoCor_GEBV[which(NoCorData$Replicates==4)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_gpdResCor_GEBV[which(ResCorData$Replicates==5)], NoCorData$cor_gpdNoCor_GEBV[which(NoCorData$Replicates==5)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_gpdResCor_GEBV[which(ResCorData$Replicates==6)], NoCorData$cor_gpdNoCor_GEBV[which(NoCorData$Replicates==6)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_gpdResCor_GEBV[which(ResCorData$Replicates==7)], NoCorData$cor_gpdNoCor_GEBV[which(NoCorData$Replicates==7)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_gpdResCor_GEBV[which(ResCorData$Replicates==8)], NoCorData$cor_gpdNoCor_GEBV[which(NoCorData$Replicates==8)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_gpdResCor_GEBV[which(ResCorData$Replicates==9)], NoCorData$cor_gpdNoCor_GEBV[which(NoCorData$Replicates==9)], paired = FALSE, alternative = "greater") # ***
t.test(ResCorData$cor_gpdResCor_GEBV[which(ResCorData$Replicates==10)], NoCorData$cor_gpdNoCor_GEBV[which(NoCorData$Replicates==10)], paired = FALSE, alternative = "greater") # ***



