# Load libraries
{
  library("ggplot2")
  library("tidyverse")
  library("wesanderson")
  
}

# set working directory
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/")

# Load data
{
  data=read.table("gpd_ResCorfulld6.txt",sep="\t",header=T)
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


# Make boxplot plot, growth per day, rhiz x-axis
sp2<-ggplot(data_ordered, aes(x= reorder(Rhizobium, growth_per_day, FUN =median), y=growth_per_day)) +
  geom_boxplot(fill=wes_palette("Rushmore1")[4],alpha=0.7) +
  ylim(-1.1,1.6) +
  xlab("Clover rhizobium") +
  ylab("growth per day") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp2
ggsave("growthperday_rhiz.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)






#### for rescor

# Calculate means and order according to means
resCor_means=aggregate(as.numeric(data$gpd_dryweight_cor), list(data$Clovershort), mean)
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
sp3<-ggplot(data_ordered1, aes(x= reorder(Clovershort, gpd_dryweight_cor, FUN =median), y=gpd_dryweight_cor)) +
  geom_boxplot(fill=wes_palette("Rushmore1")[3],alpha=0.7) +
  ylim(-1.1,1.6) +
  xlab("Clover genotype") +
  ylab("gpd_ResCor") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp3
ggsave("gpd_rescor_clover.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)


# Make boxplot plot, gpdRescor, rhiz x-axis
sp4<-ggplot(data_ordered1, aes(x= reorder(Rhizobium, gpd_dryweight_cor, FUN =median), y=gpd_dryweight_cor)) +
  geom_boxplot(fill=wes_palette("Rushmore1")[3],alpha=0.7) +
  ylim(-1.1,1.6) +
  xlab("Rhizobium genotype") +
  ylab("gpd_ResCor") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp4
ggsave("gpd_rescor_rhiz.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)





# plot together for comparison
x1 = mean(data_ordered$growth_per_day)
x2 = mean(data_ordered1$gpd_dryweight_cor)

sp5<-ggplot(data_ordered, aes(x= reorder(Clovershort, growth_per_day, FUN =median), y=growth_per_day)) +
  geom_boxplot(fill=wes_palette("Rushmore1")[4],alpha=0.7) +
  ylim(c(-0.7,1.7)) +
  xlab("Clover genotype") +
  ylab("growth per day") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp5


sp6<-ggplot(data_ordered1, aes(x= reorder(Clovershort, gpd_dryweight_cor, FUN =median), y=gpd_dryweight_cor+(x1+(-x2))))  +
  geom_boxplot(fill=wes_palette("Rushmore1")[3],alpha=0.7) +
  ylim(c(-0.7,1.7)) +
  xlab("Clover genotype") +
  ylab("gpd_ResCor") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp6

require(gridExtra)
plot=grid.arrange(sp5, sp6, ncol=2)
ggsave("CloverPhenoComparison.pdf",plot, width =60, height = 20, units = "cm",useDingbats=FALSE)



# plot together for comparison


sp7<-ggplot(data_ordered, aes(x= reorder(Rhizobium, growth_per_day, FUN =median), y=growth_per_day)) +
  geom_boxplot(fill=wes_palette("Rushmore1")[4],alpha=0.7) +
  ylim(c(-0.7,1.7)) +
  xlab("Rhizobium genotype") +
  ylab("growth per day") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp7


sp8<-ggplot(data_ordered1, aes(x= reorder(Rhizobium, gpd_dryweight_cor, FUN =median), y=gpd_dryweight_cor+(x1+(-x2))))  +
  geom_boxplot(fill=wes_palette("Rushmore1")[3],alpha=0.7) +
  ylim(c(-0.7,1.7)) +
  xlab("Rhizobium genotype") +
  ylab("gpd_ResCor") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp8

require(gridExtra)
plot=grid.arrange(sp7, sp8, ncol=2)
ggsave("RhizobiumPhenoComparison.pdf",plot, width =60, height = 20, units = "cm",useDingbats=FALSE)

