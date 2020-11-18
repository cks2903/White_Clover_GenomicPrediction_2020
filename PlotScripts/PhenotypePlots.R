# Load libraries
{
  library("ggplot2")
  library("tidyverse")
  library("wesanderson")
  library("lme4")
}

# set working directory
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/")

# Load data
{
  data=read.table("gpdtraits_fulld6_20201118.txt",sep="\t",header=T)
  head(data)
  
}



# Calculate means and order according to means
growth_per_day_means=aggregate(as.numeric(data$growth_per_day), list(data$Clovershort), mean)
colnames(growth_per_day_means)=c("Clover","growth_per_day")

growth_per_day_means_ordered <- growth_per_day_means[order(growth_per_day_means$growth_per_day),]
individual_order=growth_per_day_means_ordered$Clover

keyDF <- data.frame(key=individual_order,weight=1:length(individual_order))
growth_per_day_means_ordered <- merge(growth_per_day_means_ordered,keyDF,by.x='Clover',by.y='key',all.x=T,all.y=F)
growth_per_day_means_reordered <- growth_per_day_means_ordered[order(growth_per_day_means_ordered$weight),]
growth_per_day_means_reordered$Clover <- factor(growth_per_day_means_reordered$Clover,levels = individual_order)

# put data in order
data_ready <- merge(data,keyDF,by.x='Clovershort',by.y='key',all.x=T,all.y=F)
data_ordered <- data_ready[order(data_ready$weight),]
data_ordered$Clovershort <- factor(data_ordered$Clovershort,levels = individual_order)



# Make scatter plot, growth per day
#sp<-ggplot(data_ordered, aes(x= Clovershort, y=growth_per_day)) + geom_point(alpha=0.1,size=5)+
  #scale_color_manual(values=wes_palette("Rushmore1")[2]) + 
  #theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

#sp1 = sp + geom_point(data=growth_per_day_means_reordered,aes(Clover,growth_per_day),size=7.5,alpha=0.9)
#sp1

# Make boxplot plot, growth per day, clover x-axis
sp1<-ggplot(data_ordered, aes(x= reorder(Clovershort, growth_per_day, FUN =median), y=growth_per_day)) +
  geom_boxplot(fill=wes_palette("Rushmore1")[4],alpha=0.7) +
  ylim(-1.1,1.6) +
  xlab("Clover genotype") +
  ylab("growth per day") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp1
ggsave("growthperday_clover.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)


# plots BLUPs instead
Model1 = lmer(growth_per_day ~ factor(NS) + factor(EW) + factor(Rhizobium) + inoculation_date + (1|Clovershort), data=data_ordered) 
CloverBLUPs=ranef(Model1)$Clover
CloverBLUPs=as.data.frame(CloverBLUPs)
CloverBLUPs$Clovershort=rownames(CloverBLUPs)
colnames(CloverBLUPs)=c("BLUPsClover","Clovershort")
data_ordered_with_BLUPS1=merge(data_ordered,CloverBLUPs,by="Clovershort")


sp1.5<-ggplot(data_ordered_with_BLUPS1, aes(x= reorder(Clovershort, BLUPsClover, FUN =median), y=BLUPsClover)) +
  geom_point(color=wes_palette("Rushmore1")[4],alpha=0.7) +
  ylim(-0.3,0.3) +
  xlab("Clover genotype") +
  ylab("Clover BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp1.5
ggsave("growthperdayBLUPs_clover.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)


# Make boxplot plot, growth per day, rhiz x-axis
sp2<-ggplot(data_ordered, aes(x= reorder(Rhizobium, growth_per_day, FUN =median), y=growth_per_day)) +
  geom_boxplot(fill=wes_palette("Rushmore1")[4],alpha=0.7) +
  ylim(-1.1,1.6) +
  xlab("Clover rhizobium") +
  ylab("growth per day") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp2
ggsave("growthperday_rhiz.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)

# plots BLUPs instead
Model2 = lmer(growth_per_day ~ factor(NS) + factor(EW) + factor(Clovershort) + inoculation_date + (1|Rhizobium), data=data_ordered) 
RhizobiumBLUPs=ranef(Model2)$Rhizobium
RhizobiumBLUPs=as.data.frame(RhizobiumBLUPs)
RhizobiumBLUPs$Rhizobium=rownames(RhizobiumBLUPs)
colnames(RhizobiumBLUPs)=c("BLUPsRhiz","Rhizobium")
data_ordered_with_BLUPS2=merge(data_ordered,RhizobiumBLUPs,by="Rhizobium")


sp2.5<-ggplot(data_ordered_with_BLUPS2, aes(x= reorder(Rhizobium, BLUPsRhiz, FUN =median), y=BLUPsRhiz)) +
  geom_point(color=wes_palette("Rushmore1")[4],alpha=0.7) +
  ylim(-0.3,0.3) +
  xlab("Rhizobium genotype") +
  ylab("Rhizobium BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp2.5
ggsave("growthperdayBLUPs_rhiz.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)


#### for rescor

# Calculate means and order according to means
resCor_means=aggregate(as.numeric(data$gpdResCor), list(data$Clovershort), mean)
colnames(resCor_means)=c("Clover","gpd_rescor")

resCor_means_ordered <- resCor_means[order(resCor_means$gpd_rescor),]
individual_order=resCor_means_ordered$Clover

keyDF <- data.frame(key=individual_order,weight=1:length(individual_order))
resCor_means_ordered <- merge(resCor_means_ordered,keyDF,by.x='Clover',by.y='key',all.x=T,all.y=F)
resCor_means_ordered_reordered <- resCor_means_ordered[order(resCor_means_ordered$weight),]
resCor_means_ordered_reordered$Clover <- factor(resCor_means_ordered_reordered$Clover,levels = individual_order)

# put data in order
data_ready <- merge(data,keyDF,by.x='Clovershort',by.y='key',all.x=T,all.y=F)
data_ordered1 <- data_ready[order(data_ready$weight),]
data_ordered1$Clovershort <- factor(data_ordered1$Clovershort,levels = individual_order)





# Make boxplot plot, gpdRescor, clover x-axis
sp3<-ggplot(data_ordered1, aes(x= reorder(Clovershort, gpdResCor, FUN =median), y=gpdResCor)) +
  geom_boxplot(fill=wes_palette("Rushmore1")[3],alpha=0.7) +
  ylim(-1.1,1.6) +
  xlab("Clover genotype") +
  ylab("gpd_ResCor") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp3
ggsave("growthperdayBLUPs_clover.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)


# plots BLUPs instead
Model3 = lmer(gpdResCor ~ factor(NS) + factor(EW) + factor(Rhizobium) + inoculation_date + (1|Clovershort), data=data_ordered1) 
CloverBLUPs=ranef(Model3)$Clover
CloverBLUPs=as.data.frame(CloverBLUPs)
CloverBLUPs$Clovershort=rownames(CloverBLUPs)
colnames(CloverBLUPs)=c("BLUPsClover","Clovershort")
data_ordered_with_BLUPS3=merge(data_ordered1,CloverBLUPs,by="Clovershort")


sp3.5<-ggplot(data_ordered_with_BLUPS3, aes(x= reorder(Clovershort, BLUPsClover, FUN =median), y=BLUPsClover)) +
  geom_point(color=wes_palette("Rushmore1")[3],alpha=0.7) +
  ylim(-0.3,0.3) +
  xlab("Clover genotype") +
  ylab("Clover BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp3.5
ggsave("growthperdayBLUPs_clover(gpdres).pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)


# Make boxplot plot, gpdRescor, rhiz x-axis
sp4<-ggplot(data_ordered1, aes(x= reorder(Rhizobium, gpdResCor, FUN =median), y=gpdResCor)) +
  geom_boxplot(fill=wes_palette("Rushmore1")[3],alpha=0.7) +
  ylim(-1.1,1.6) +
  xlab("Rhizobium genotype") +
  ylab("gpd_ResCor") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp4
ggsave("gpd_rescor_rhiz.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)


# plots BLUPs instead
Model4 = lmer(gpdResCor ~ factor(NS) + factor(EW) + factor(Clovershort) + inoculation_date + (1|Rhizobium), data=data_ordered1) 
RhizobiumBLUPs=ranef(Model4)$Rhizobium
RhizobiumBLUPs=as.data.frame(RhizobiumBLUPs)
RhizobiumBLUPs$Rhizobium=rownames(RhizobiumBLUPs)
colnames(RhizobiumBLUPs)=c("BLUPsRhiz","Rhizobium")
data_ordered_with_BLUPS4=merge(data_ordered1,RhizobiumBLUPs,by="Rhizobium")


sp4.5<-ggplot(data_ordered_with_BLUPS4, aes(x= reorder(Rhizobium, BLUPsRhiz, FUN =median), y=BLUPsRhiz)) +
  geom_point(color=wes_palette("Rushmore1")[3],alpha=0.7) +
  ylim(-0.3,0.3) +
  xlab("Rhizobium genotype") +
  ylab("Rhizobium BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp4.5
ggsave("growthperdayBLUPs_rhiz(gpdres).pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)



#### for fixcor

# Calculate means and order according to means
FixCor_means=aggregate(as.numeric(data$gpdFixCor), list(data$Clovershort), mean)
colnames(FixCor_means)=c("Clover","gpd_fixcor")

FixCor_means_ordered <- FixCor_means[order(FixCor_means$gpd_fixcor),]
individual_order=FixCor_means_ordered$Clover

keyDF <- data.frame(key=individual_order,weight=1:length(individual_order))
FixCor_means_ordered <- merge(FixCor_means_ordered,keyDF,by.x='Clover',by.y='key',all.x=T,all.y=F)
FixCor_means_ordered_reordered <- FixCor_means_ordered[order(FixCor_means_ordered$weight),]
FixCor_means_ordered_reordered$Clover <- factor(FixCor_means_ordered_reordered$Clover,levels = individual_order)

# put data in order
data_ready <- merge(data,keyDF,by.x='Clovershort',by.y='key',all.x=T,all.y=F)
data_ordered10 <- data_ready[order(data_ready$weight),]
data_ordered10$Clovershort <- factor(data_ordered10$Clovershort,levels = individual_order)



# Make boxplot plot, gpdRescor, clover x-axis
sp11<-ggplot(data_ordered10, aes(x= reorder(Clovershort, gpdFixCor, FUN =median), y=gpdFixCor)) +
  geom_boxplot(fill=wes_palette("Rushmore1")[2],alpha=0.7) +
  ylim(-1.2,1.6) +
  xlab("Clover genotype") +
  ylab("gpd_FixCor") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp11
ggsave("gpd_fixcor_clover.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)


# plots BLUPs instead
Model10 = lmer(gpdFixCor ~ factor(NS) + factor(EW) + factor(Rhizobium) + inoculation_date + (1|Clovershort), data=data_ordered10) 
CloverBLUPs=ranef(Model10)$Clover
CloverBLUPs=as.data.frame(CloverBLUPs)
CloverBLUPs$Clovershort=rownames(CloverBLUPs)
colnames(CloverBLUPs)=c("BLUPsClover","Clovershort")
data_ordered_with_BLUPS10=merge(data_ordered10,CloverBLUPs,by="Clovershort")


sp12<-ggplot(data_ordered_with_BLUPS10, aes(x= reorder(Clovershort, BLUPsClover, FUN =median), y=BLUPsClover)) +
  geom_point(color=wes_palette("Rushmore1")[2],alpha=0.7) +
  ylim(-0.3,0.3) +
  xlab("Clover genotype") +
  ylab("Clover BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp12
ggsave("growthperdayBLUPs_clover(gpdfixcor).pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)


# Make boxplot plot, gpdRescor, rhiz x-axis
sp13<-ggplot(data_ordered_with_BLUPS10, aes(x= reorder(Rhizobium, gpdFixCor, FUN =median), y=gpdFixCor)) +
  geom_boxplot(fill=wes_palette("Rushmore1")[2],alpha=0.7) +
  ylim(-1.2,1.6) +
  xlab("Rhizobium genotype") +
  ylab("gpd_FixCor") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp13
ggsave("gpd_fixcor_rhiz.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)


# plots BLUPs instead
Model14 = lmer(gpdFixCor ~ factor(NS) + factor(EW) + factor(Clovershort) + inoculation_date + (1|Rhizobium), data=data_ordered10) 
RhizobiumBLUPs_fixed=ranef(Model14)$Rhizobium
RhizobiumBLUPs_fixed=as.data.frame(RhizobiumBLUPs_fixed)
RhizobiumBLUPs_fixed$Rhizobium=rownames(RhizobiumBLUPs_fixed)
colnames(RhizobiumBLUPs_fixed)=c("BLUPsRhiz","Rhizobium")
data_ordered_with_BLUPS14=merge(data_ordered10,RhizobiumBLUPs_fixed,by="Rhizobium")


sp14<-ggplot(data_ordered_with_BLUPS14, aes(x= reorder(Rhizobium, BLUPsRhiz, FUN =median), y=BLUPsRhiz)) +
  geom_point(color=wes_palette("Rushmore1")[2],alpha=0.7) +
  ylim(-0.3,0.3) +
  xlab("Rhizobium genotype") +
  ylab("Rhizobium BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp14
ggsave("growthperdayBLUPs_rhiz(gpdfix).pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)




# plot together for comparison, gpd values, clover
#x1 = mean(data_ordered$growth_per_day)
#x2 = mean(data_ordered1$gpd_dryweight_cor)

#sp5<-ggplot(data_ordered, aes(x= reorder(Clovershort, growth_per_day, FUN =median), y=growth_per_day)) +
  #geom_boxplot(fill=wes_palette("Rushmore1")[4],alpha=0.7) +
  #ylim(c(-0.7,1.7)) +
  #xlab("Clover genotype") +
  #ylab("growth per day") +
  #theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
#sp5


#sp6<-ggplot(data_ordered1, aes(x= reorder(Clovershort, gpd_dryweight_cor, FUN =median), y=gpd_dryweight_cor+(x1+(-x2))))  +
  #geom_boxplot(fill=wes_palette("Rushmore1")[3],alpha=0.7) +
  #ylim(c(-0.7,1.7)) +
  #xlab("Clover genotype") +
  #ylab("gpd_ResCor") +
  #theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
#sp6

require(gridExtra)
#plot=grid.arrange(sp5, sp6, ncol=2)
#ggsave("CloverPhenoComparison.pdf",plot, width =60, height = 20, units = "cm",useDingbats=FALSE)


# plot together for comparison, BLUPS clover
#sp5<-ggplot(data_ordered_with_BLUPS1, aes(x= reorder(Clovershort, BLUPsClover, FUN =median), y=BLUPsClover)) +
  #geom_point(color=wes_palette("Rushmore1")[4],alpha=0.7) +
  #ylim(c(-0.3,0.3)) +
  #xlab("Clover genotype") +
  #ylab("BLUPs Clover") +
  #theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
#sp5


#sp6<-ggplot(data_ordered_with_BLUPS3, aes(x= reorder(Clovershort, BLUPsClover, FUN =median), y=BLUPsClover))  +
  #geom_point(color=wes_palette("Rushmore1")[3],alpha=0.7) +
  #ylim(c(-0.3,0.3)) +
  #xlab("Clover genotype") +
  #ylab("BLUPs Clover (gpdResCor)") +
  #theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
#sp6

#require(gridExtra)
#plot=grid.arrange(sp5, sp6, ncol=2)
#ggsave("CloverBLUPComparison.pdf",plot, width =60, height = 20, units = "cm",useDingbats=FALSE)

######## clover in same plot

sp6.5=ggplot() +
  geom_point(data = data_ordered_with_BLUPS3, aes(x= reorder(Clovershort, BLUPsClover, FUN =median), y=BLUPsClover), color = "#0B775E",cex=5,alpha=0.0) +
  geom_point(data = data_ordered_with_BLUPS1, aes(x = Clovershort, y = BLUPsClover), color = "#35274A",cex=5,alpha=0.08) +
  geom_point(data = data_ordered_with_BLUPS10, aes(x = Clovershort, y = BLUPsClover), color = wes_palette("Rushmore1")[2],cex=5,alpha=0.08) +
  geom_point(data = data_ordered_with_BLUPS3, aes(x= reorder(Clovershort, BLUPsClover, FUN =median), y=BLUPsClover), color = "#0B775E",cex=5,alpha=0.08) +
  expand_limits(x = 150) +
  ylim(c(-0.3,0.3)) +
  xlab("Clover genotype") +
  ylab("BLUPs Clover") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp6.5
ggsave("CloverBLUPComparison_v2.pdf",sp6.5, width =30, height = 20, units = "cm",useDingbats=FALSE)




# plot together for comparison
#sp7<-ggplot(data_ordered, aes(x= reorder(Rhizobium, growth_per_day, FUN =median), y=growth_per_day)) +
  #geom_boxplot(fill=wes_palette("Rushmore1")[4],alpha=0.7) +
  #ylim(c(-0.7,1.7)) +
  #xlab("Rhizobium genotype") +
  #ylab("growth per day") +
  #theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
#sp7


#sp8<-ggplot(data_ordered1, aes(x= reorder(Rhizobium, gpd_dryweight_cor, FUN =median), y=gpd_dryweight_cor+(x1+(-x2))))  +
  #geom_boxplot(fill=wes_palette("Rushmore1")[3],alpha=0.7) +
  #ylim(c(-0.7,1.7)) +
  #xlab("Rhizobium genotype") +
  #ylab("gpd_ResCor") +
  #theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
#sp8

#require(gridExtra)
#plot=grid.arrange(sp7, sp8, ncol=2)
#ggsave("RhizobiumPhenoComparison.pdf",plot, width =60, height = 20, units = "cm",useDingbats=FALSE)


#sp5<-ggplot(data_ordered_with_BLUPS2, aes(x= reorder(Rhizobium, BLUPsRhiz, FUN =median), y=BLUPsRhiz)) +
  #geom_point(color=wes_palette("Rushmore1")[4],alpha=0.7) +
  #ylim(c(-0.3,0.3)) +
  #xlab("Rhizobium genotype") +
  #ylab("BLUPs Rhizobium") +
  #theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
#sp5


#sp6<-ggplot(data_ordered_with_BLUPS4, aes(x= reorder(Rhizobium, BLUPsRhiz, FUN =median), y=BLUPsRhiz))  +
  #geom_point(color=wes_palette("Rushmore1")[3],alpha=0.7) +
  #ylim(c(-0.3,0.3)) +
  #xlab("Rhizobium genotype") +
  #ylab("BLUPs Rhizobium (gpdResCor)") +
  #theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
#sp6

#require(gridExtra)
#plot=grid.arrange(sp5, sp6, ncol=2)
#plot
#ggsave("RhizobiumBLUPComparison.pdf",plot, width =60, height = 20, units = "cm",useDingbats=FALSE)

# rhiz in same plot
sp7=ggplot() +
  geom_point(data = data_ordered_with_BLUPS4, aes(x = reorder(Rhizobium, BLUPsRhiz, FUN =median), y = BLUPsRhiz), color = "#0B775E",cex=5,alpha=0.00) +
  geom_point(data = data_ordered_with_BLUPS2, aes(x = Rhizobium, y = BLUPsRhiz), color = "#35274A",cex=5,alpha=0.08) +
  geom_point(data = data_ordered_with_BLUPS14, aes(x = Rhizobium, y = BLUPsRhiz), color = wes_palette("Rushmore1")[2],cex=5,alpha=0.08) +
  geom_point(data = data_ordered_with_BLUPS4, aes(x = reorder(Rhizobium, BLUPsRhiz, FUN =median), y = BLUPsRhiz), color = "#0B775E",cex=5,alpha=0.08) +
  expand_limits(x = 175) +
  ylim(c(-0.3,0.3)) +
  xlab("Rhizobium genotype") +
  ylab("BLUPs Rhizobium") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp7
ggsave("RhizobiumBLUPComparison_v2.pdf",sp7, width =30, height = 20, units = "cm",useDingbats=FALSE)
