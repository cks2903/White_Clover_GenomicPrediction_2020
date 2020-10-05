

library("ggplot2")
library(psych)

setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831")

# in F0
NoCorData=read.table("replicate_reduction_gpdNoCor.txt",head=T)
dim(NoCorData)

# in F0
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
  theme_classic()

p1

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/CorrectedPhenoCorrelation_gpdNoCorANDgpdResCor.pdf', plot = p1, width = 15, height = 15, unit = 'cm')

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
  theme_classic()
p1.5

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/Correlation_gpdNoCorANDgpdResCor_pred_SDs.pdf', plot = p1.5, width = 15, height = 15, unit = 'cm')





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
  theme_classic()

p2

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/Correlation_gpdNoCorANDgpdResCor.pdf', plot = p2, width = 15, height = 15, unit = 'cm')

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














############################
####    f1 population
############################

# in F0
NoCorData_F1=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/F1_pred_gpd_NoCor/replicate_reduction_F1summary.txt",head=T)
dim(NoCorData_F1)

# in F0
ResCorData_F1=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/F1_pred_gpd_ResCor/replicate_reduction_F1summary.txt",head=T)
dim(NoCorData_F1)

# Combined
ResCorData_means_corrected_F1=aggregate(ResCorData_F1$correlation, list(ResCorData_F1$Replicates), mean)
colnames(ResCorData_means_corrected_F1)=c("Replicates","cor_CorrectedPheno_GEBV")

NoCorData_means_corrected_F1=aggregate(NoCorData_F1$correlation, list(NoCorData_F1$Replicates), mean)
colnames(NoCorData_means_corrected_F1)=c("Replicates","cor_CorrectedPheno_GEBV")

# plot together
p11 = ggplot() + 
  geom_line(data = ResCorData_means_corrected_F1, aes(x = Replicates, y = cor_CorrectedPheno_GEBV), color = "#0B775E") +
  geom_line(data = NoCorData_means_corrected_F1, aes(x = Replicates, y = cor_CorrectedPheno_GEBV), color = "#35274A") +
  scale_x_reverse() +
  xlab('Replicates') +
  ylab('Correlation') +
  theme_classic()

p11

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/CorrectedPhenoCorrelation_gpdNoCorANDgpdResCor_F1.pdf', plot = p11, width = 15, height = 15, unit = 'cm')




# test if they differ significantly from each other in F1.
# remember they are actually dependend as (GEBV1,F1_Avg), (GEBV2,F1_Avg)
# Therefore Hotelling-Williams-test

rawresults_rescorgpd_F1=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/F1_pred_gpd_ResCor/replicate_reduction_F1Results.txt",head=T,sep="\t")
head(rawresults_rescorgpd_F1)
dim(rawresults_rescorgpd_F1)

rawresults_nocorgpd_F1=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/F1_pred_gpd_NoCor/replicate_reduction_F1Results.txt",head=T,sep="\t")
head(rawresults_nocorgpd_F1)
dim(rawresults_nocorgpd_F1)

# 10 reps
Ten_reps_F1_gpdResCor=rawresults_rescorgpd_F1[which(rawresults_rescorgpd_F1$Replicates==10),]
Ten_reps_F1_gpdResCor_means=aggregate(Ten_reps_F1_gpdResCor$Mean.Parental.predicted.GEBV, list(Ten_reps_F1_gpdResCor$Round), mean)

Ten_reps_F1_gpdNoCor=rawresults_nocorgpd_F1[which(rawresults_nocorgpd_F1$Replicates==10),]
Ten_reps_F1_gpdNoCor_means=aggregate(Ten_reps_F1_gpdNoCor$Mean.Parental.predicted.GEBV, list(Ten_reps_F1_gpdNoCor$Round), mean)

correlation_ofGEBVs=cor(Ten_reps_F1_gpdResCor_means[,2],Ten_reps_F1_gpdNoCor_means[,2])
r.test(n=900, r12=ResCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[9], r13=NoCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[9],r23=correlation_ofGEBVs) # N.S.

# 9 reps
Nine_reps_F1_gpdResCor=rawresults_rescorgpd_F1[which(rawresults_rescorgpd_F1$Replicates==9),]
Nine_reps_F1_gpdNoCor=rawresults_nocorgpd_F1[which(rawresults_nocorgpd_F1$Replicates==9),]
r23=cor(Nine_reps_F1_gpdResCor$Mean.Parental.predicted.GEBV,Nine_reps_F1_gpdNoCor$Mean.Parental.predicted.GEBV)
r.test(n=900, r12=ResCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[8], r13=NoCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[8],r23=r23) # N.S

# 8 reps
Eight_reps_F1_gpdResCor=rawresults_rescorgpd_F1[which(rawresults_rescorgpd_F1$Replicates==8),]
Eight_reps_F1_gpdNoCor=rawresults_nocorgpd_F1[which(rawresults_nocorgpd_F1$Replicates==8),]
r23=cor(Eight_reps_F1_gpdResCor$Mean.Parental.predicted.GEBV,Eight_reps_F1_gpdNoCor$Mean.Parental.predicted.GEBV)
r.test(n=900, r12=ResCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[7], r13=NoCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[7],r23=r23) # N.S

# 7 reps
Seven_reps_F1_gpdResCor=rawresults_rescorgpd_F1[which(rawresults_rescorgpd_F1$Replicates==7),]
Seven_reps_F1_gpdNoCor=rawresults_nocorgpd_F1[which(rawresults_nocorgpd_F1$Replicates==7),]
r23=cor(Seven_reps_F1_gpdResCor$Mean.Parental.predicted.GEBV,Seven_reps_F1_gpdNoCor$Mean.Parental.predicted.GEBV)
r.test(n=900, r12=ResCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[6], r13=NoCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[6],r23=r23) # N.S

# 6 reps
Six_reps_F1_gpdResCor=rawresults_rescorgpd_F1[which(rawresults_rescorgpd_F1$Replicates==6),]
Six_reps_F1_gpdNoCor=rawresults_nocorgpd_F1[which(rawresults_nocorgpd_F1$Replicates==6),]
r23=cor(Six_reps_F1_gpdResCor$Mean.Parental.predicted.GEBV,Six_reps_F1_gpdNoCor$Mean.Parental.predicted.GEBV)
r.test(n=900, r12=ResCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[5], r13=NoCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[5],r23=r23) # N.S


# 5 reps
five_reps_F1_gpdResCor=rawresults_rescorgpd_F1[which(rawresults_rescorgpd_F1$Replicates==5),]
five_reps_F1_gpdNoCor=rawresults_nocorgpd_F1[which(rawresults_nocorgpd_F1$Replicates==5),]
r23=cor(five_reps_F1_gpdResCor$Mean.Parental.predicted.GEBV,five_reps_F1_gpdNoCor$Mean.Parental.predicted.GEBV)
r.test(n=900, r12=ResCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[4], r13=NoCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[4],r23=r23) # N.S

# 4 reps
four_reps_F1_gpdResCor=rawresults_rescorgpd_F1[which(rawresults_rescorgpd_F1$Replicates==4),]
four_reps_F1_gpdNoCor=rawresults_nocorgpd_F1[which(rawresults_nocorgpd_F1$Replicates==4),]
r23=cor(four_reps_F1_gpdResCor$Mean.Parental.predicted.GEBV,four_reps_F1_gpdNoCor$Mean.Parental.predicted.GEBV)
r.test(n=900, r12=ResCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[3], r13=NoCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[3],r23=r23) # N.S

# 3 reps
three_reps_F1_gpdResCor=rawresults_rescorgpd_F1[which(rawresults_rescorgpd_F1$Replicates==3),]
three_reps_F1_gpdNoCor=rawresults_nocorgpd_F1[which(rawresults_nocorgpd_F1$Replicates==3),]
r23=cor(three_reps_F1_gpdResCor$Mean.Parental.predicted.GEBV,three_reps_F1_gpdNoCor$Mean.Parental.predicted.GEBV)
r.test(n=900, r12=ResCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[2], r13=NoCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[2],r23=r23) # N.S

# 2 reps
two_reps_F1_gpdResCor=rawresults_rescorgpd_F1[which(rawresults_rescorgpd_F1$Replicates==2),]
two_reps_F1_gpdNoCor=rawresults_nocorgpd_F1[which(rawresults_nocorgpd_F1$Replicates==2),]
r23=cor(two_reps_F1_gpdResCor$Mean.Parental.predicted.GEBV,two_reps_F1_gpdNoCor$Mean.Parental.predicted.GEBV)
r.test(n=900, r12=ResCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[1], r13=NoCorData_means_corrected_F1$cor_CorrectedPheno_GEBV[1],r23=r23) # N.S
