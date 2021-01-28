
###############################################################
###############################################################
#       A script for investigation the GEBVs for each trait and 
#       correlation for each round
###############################################################
###############################################################


###############################################################
###############################################################
#
#         based on averages, so one value pr. genotype
###############################################################
###############################################################

library(corrgram)
library(ggplot2)
library("wesanderson")

# iSize

# Load predictions from iSize
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/iSize_20210120/")
# Load files

list.files()

file_list <- list.files(pattern="Predictions_")

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset1")){
    dataset1 <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset1")){
    pred_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset1<-rbind(dataset1, pred_dataset)
    rm(pred_dataset)
  }
  
}

nrow(dataset1)
dataset1=as.data.frame(dataset1)
head(dataset1)

# Calculate mean for each individual for the 100 times its GEBV was estimated
meaniSize=aggregate(as.numeric(dataset1$Observed), list(dataset1$ID), mean)[,2]
meangebv=aggregate(as.numeric(dataset1$GEBV), list(dataset1$ID), mean)
df_iSize=cbind(meangebv,meaniSize)
colnames(df_iSize)=c("Individual","iSize_GEBV","iSize_obs")
head(df_iSize)


dataframeWithResults=rep(0,2)
dataframeWithResults=as.data.frame(dataframeWithResults)
dataframeWithResults=t(dataframeWithResults)
colnames(dataframeWithResults)=c("iSize","Round")

# Load data from GBLUP gpd_ResCor
AllCorrelationResults=list.files(pattern="^Correlation")


for (i in seq(1:length(AllCorrelationResults))){
  filename=AllCorrelationResults[i]
  round=str_sub(filename, 22, 23)
  if (str_sub(round,2,2)=="."){
    round=str_sub(round,1,1)
  }else if (str_sub(filename,22,24)=="100"){
    round=str_sub(filename,22,24)
  }
  table=read.table(AllCorrelationResults[i],stringsAsFactors = F)
  table_df=as.data.frame(table)
  table_df$Round=as.numeric(round)
  colnames(table_df)[1]=c("iSize")
  dataframeWithResults=rbind(dataframeWithResults,table_df)
}

dataframeWithResults_iSize=dataframeWithResults[2:nrow(dataframeWithResults),]
dim(dataframeWithResults_iSize)





# gpdResCor
# Load predictions from gpdResCor
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/gpd_20210120/")
# Load files

list.files()

file_list <- list.files(pattern="Predictions_")

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset2")){
    dataset2 <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset2")){
    pred_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset2<-rbind(dataset2, pred_dataset)
    rm(pred_dataset)
  }
  
}

nrow(dataset2)
dataset2=as.data.frame(dataset2)
head(dataset2)

# Calculate mean for each individual for the 100 times its GEBV was estimated
meangpdrescor=aggregate(as.numeric(dataset2$Observed), list(dataset2$ID), mean)[,2]
meangebv=aggregate(as.numeric(dataset2$GEBV), list(dataset2$ID), mean)
df_gpdRescor=cbind(meangebv,meangpdrescor)
colnames(df_gpdRescor)=c("Individual","gpd_GEBV","gpd_obs")
head(df_gpdRescor)




dataframeWithResults=rep(0,2)
dataframeWithResults=as.data.frame(dataframeWithResults)
dataframeWithResults=t(dataframeWithResults)
colnames(dataframeWithResults)=c("gpd","Round")

# Load data from GBLUP gpd_ResCor
AllCorrelationResults=list.files(pattern="^Correlation")


for (i in seq(1:length(AllCorrelationResults))){
  filename=AllCorrelationResults[i]
  round=str_sub(filename, 22, 23)
  if (str_sub(round,2,2)=="."){
    round=str_sub(round,1,1)
  }else if (str_sub(filename,22,24)=="100"){
    round=str_sub(filename,22,24)
  }
  table=read.table(AllCorrelationResults[i],stringsAsFactors = F)
  table_df=as.data.frame(table)
  table_df$Round=as.numeric(round)
  colnames(table_df)[1]=c("gpd")
  dataframeWithResults=rbind(dataframeWithResults,table_df)
}

dataframeWithResults_gpd=dataframeWithResults[2:nrow(dataframeWithResults),]
dim(dataframeWithResults_gpd)




# gpdFxCor
# Load predictions from gpdFixCor
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/gpdFixCor_afteraveraging_20210121/")
# Load files

list.files()

file_list <- list.files(pattern="Predictions_")

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset3")){
    dataset3 <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset3")){
    pred_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset3<-rbind(dataset3, pred_dataset)
    rm(pred_dataset)
  }
  
}

nrow(dataset3)
dataset3=as.data.frame(dataset3)
head(dataset3)

# Calculate mean for each individual for the 100 times its GEBV was estimated
meangpdFixCor=aggregate(as.numeric(dataset3$Observed), list(dataset3$ID), mean)[,2]
meangebv=aggregate(as.numeric(dataset3$GEBV), list(dataset3$ID), mean)
df_gpdFixCor=cbind(meangebv,meangpdFixCor)
colnames(df_gpdFixCor)=c("Individual","gpdFixCor_GEBV","gpdFixCor_obs")
head(df_gpdFixCor)


dataframeWithResults=rep(0,2)
dataframeWithResults=as.data.frame(dataframeWithResults)
dataframeWithResults=t(dataframeWithResults)
colnames(dataframeWithResults)=c("gpdFixCor","Round")

# Load data from GBLUP gpd_ResCor
AllCorrelationResults=list.files(pattern="^Correlation")


for (i in seq(1:length(AllCorrelationResults))){
  filename=AllCorrelationResults[i]
  round=str_sub(filename, 22, 23)
  if (str_sub(round,2,2)=="."){
    round=str_sub(round,1,1)
  }else if (str_sub(filename,22,24)=="100"){
    round=str_sub(filename,22,24)
  }
  table=read.table(AllCorrelationResults[i],stringsAsFactors = F)
  table_df=as.data.frame(table)
  table_df$Round=as.numeric(round)
  colnames(table_df)[1]=c("gpdFixCor")
  dataframeWithResults=rbind(dataframeWithResults,table_df)
}

dataframeWithResults_gpdFixCor=dataframeWithResults[2:nrow(dataframeWithResults),]
dim(dataframeWithResults_gpdFixCor)





# gpi
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/gpdDay11to25/")
# Load files

list.files()

file_list <- list.files(pattern="Predictions_")

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset4")){
    dataset4 <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset4")){
    pred_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset4<-rbind(dataset4, pred_dataset)
    rm(pred_dataset)
  }
  
}

nrow(dataset4)
dataset4=as.data.frame(dataset4)
head(dataset4)

# Calculate mean for each individual for the 100 times its GEBV was estimated
meangpi=aggregate(as.numeric(dataset4$Observed), list(dataset4$ID), mean)[,2]
meangebv=aggregate(as.numeric(dataset4$GEBV), list(dataset4$ID), mean)
df_meangpi=cbind(meangebv,meangpi)
colnames(df_meangpi)=c("Individual","gpi_GEBV","gpi_obs")
head(df_meangpi)


dataframeWithResults=rep(0,2)
dataframeWithResults=as.data.frame(dataframeWithResults)
dataframeWithResults=t(dataframeWithResults)
colnames(dataframeWithResults)=c("gpi","Round")

# Load data from GBLUP gpd_ResCor
AllCorrelationResults=list.files(pattern="^Correlation")


for (i in seq(1:length(AllCorrelationResults))){
  filename=AllCorrelationResults[i]
  round=str_sub(filename, 22, 23)
  if (str_sub(round,2,2)=="."){
    round=str_sub(round,1,1)
  }else if (str_sub(filename,22,24)=="100"){
    round=str_sub(filename,22,24)
  }
  table=read.table(AllCorrelationResults[i],stringsAsFactors = F)
  table_df=as.data.frame(table)
  table_df$Round=as.numeric(round)
  colnames(table_df)[1]=c("gpi")
  dataframeWithResults=rbind(dataframeWithResults,table_df)
}

dataframeWithResults_gpi=dataframeWithResults[2:nrow(dataframeWithResults),]
dim(dataframeWithResults_gpi)



# gpiCor
# Load predictions from gpdday11to25
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/gpdDay11to25_correctedForiSize/")
# Load files

list.files()

file_list <- list.files(pattern="Predictions_")

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset5")){
    dataset5 <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset5")){
    pred_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset5<-rbind(dataset5, pred_dataset)
    rm(pred_dataset)
  }
  
}

nrow(dataset5)
dataset5=as.data.frame(dataset5)
head(dataset5)

# Calculate mean for each individual for the 100 times its GEBV was estimated
meangpiCor=aggregate(as.numeric(dataset5$Observed), list(dataset5$ID), mean)[,2]
meangebv=aggregate(as.numeric(dataset5$GEBV), list(dataset5$ID), mean)
df_gpiCor=cbind(meangebv,meangpiCor)
colnames(df_gpiCor)=c("Individual","gpiCor_GEBV","gpiCor_obs")
head(df_gpiCor)




dataframeWithResults=rep(0,2)
dataframeWithResults=as.data.frame(dataframeWithResults)
dataframeWithResults=t(dataframeWithResults)
colnames(dataframeWithResults)=c("gpiCor","Round")

# Load data from GBLUP gpd_ResCor
AllCorrelationResults=list.files(pattern="^Correlation")


for (i in seq(1:length(AllCorrelationResults))){
  filename=AllCorrelationResults[i]
  round=str_sub(filename, 22, 23)
  if (str_sub(round,2,2)=="."){
    round=str_sub(round,1,1)
  }else if (str_sub(filename,22,24)=="100"){
    round=str_sub(filename,22,24)
  }
  table=read.table(AllCorrelationResults[i],stringsAsFactors = F)
  table_df=as.data.frame(table)
  table_df$Round=as.numeric(round)
  colnames(table_df)[1]=c("gpiCor")
  dataframeWithResults=rbind(dataframeWithResults,table_df)
}

dataframeWithResults_gpiCor=dataframeWithResults[2:nrow(dataframeWithResults),]
dim(dataframeWithResults_gpiCor)





####### plot result for each round
# Figure 4C

head(dataframeWithResults_gpiCor)


p3 = ggplot() + 
  geom_line(data = dataframeWithResults_iSize, aes(x = Round, y = iSize), color = wes_palette("Rushmore1")[4]) +
  geom_line(data = dataframeWithResults_gpd, aes(x = Round, y = gpd), color = "#0B775E") +
  geom_line(data = dataframeWithResults_gpdFixCor, aes(x = Round, y = gpdFixCor), color =  wes_palette("Rushmore1")[2]) +
  geom_line(data = dataframeWithResults_gpi, aes(x = Round, y = gpi), color = "brown4") +
  geom_line(data = dataframeWithResults_gpiCor, aes(x = Round, y = gpiCor), color = "cadetblue") +
  
  xlab('Round') +
  ylab('Correlation') +
  theme_classic()

p3
ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/Correlations_pr_round_AllTraits.pdf', plot = p3, width = 30, height = 20, unit = 'cm')



####### plot average GEBV for each traits
# Figure 4E

head(dataset1)
head(dataset2)
head(dataset3)
head(dataset4)
head(dataset5)

Merged=rep(0,nrow(dataset5)*5)
Merged=as.data.frame(Merged)
Merged$Ind=dataset1$ID
head(Merged)
Merged$trait=rep(0,nrow(dataset5)*5)
Merged$GEBV=rep(0,nrow(dataset5)*5)
Merged=Merged[,2:ncol(Merged)]

Merged$trait[1:14645]="iSize"
Merged$GEBV[1:14645]=dataset1$GEBV

Merged$trait[14646:29290]="gpd"
Merged$GEBV[14646:29290]=dataset2$GEBV

Merged$trait[29291:43935]="gpdFixCor"
Merged$GEBV[29291:43935]=dataset3$GEBV

Merged$trait[43936:58580]="gpi"
Merged$GEBV[43936:58580]=dataset4$GEBV

Merged$trait[58581:73225]="gpiCor"
Merged$GEBV[58581:73225]=dataset5$GEBV

sp<-ggplot(Merged, aes(x= Ind, y=GEBV, color=trait)) + geom_point(alpha=0.01,size=5)+
  scale_color_manual(values= c(wes_palette("Rushmore1")[4],wes_palette("Rushmore1")[2],wes_palette("Rushmore1")[3],"brown4", "cadetblue")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

sp

df_ordered <- df_gpdRescor[order(df_gpdRescor$gpd_GEBV),] # sort after GEBVs for gpd
individual_order=df_ordered$Individual

keyDF <- data.frame(key=individual_order,weight=1:length(individual_order))
colnames(keyDF)[1]="Ind"
combined_ready1 <- merge(Merged,keyDF,by='Ind',all.x=T,all.y=F)
combined_ordered1 <- combined_ready1[order(combined_ready1$weight),c('Ind','trait','GEBV', 'weight')]
combined_ordered1$Ind <- factor(combined_ordered1$Ind,levels = individual_order)


sp1 = sp + geom_point(data=combined_ordered1,aes(Ind,GEBV,color=trait),size=7.5,alpha=0.9)
sp1


# plot individually

# isize
iSize_newdf = combined_ordered1[which(combined_ordered1$trait=="iSize"),]
iSize_newdf_means = aggregate(as.numeric(iSize_newdf$GEBV),list(iSize_newdf$Ind), mean)
colnames(iSize_newdf_means)=c("Ind","GEBV")

sp<-ggplot(iSize_newdf, aes(x= Ind, y=GEBV)) + geom_point(alpha=0.05,size=2,col=wes_palette("Rushmore1")[4])+
  #scale_color_manual(values= c(wes_palette("Rushmore1")[4],wes_palette("Rushmore1")[2],wes_palette("Rushmore1")[3],"brown4", "cadetblue")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))+  expand_limits(x = 175)
sp1 = sp + geom_point(data=iSize_newdf_means,aes(Ind,GEBV),size=3,alpha=0.9,col=wes_palette("Rushmore1")[4]) +expand_limits(x = 175)
sp1
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/GEBVs_iSize_sortedAfterAvggpd_",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)


# gpd
df_gpdRescor
gpd_newdf = combined_ordered1[which(combined_ordered1$trait=="gpd"),]
gpd_newdf_means = aggregate(as.numeric(gpd_newdf$GEBV),list(gpd_newdf$Ind), mean)
colnames(gpd_newdf_means)=c("Ind","GEBV")

sp<-ggplot(gpd_newdf, aes(x= Ind, y=GEBV)) + geom_point(alpha=0.05,size=2,col=wes_palette("Rushmore1")[3])+
  #scale_color_manual(values= c(wes_palette("Rushmore1")[4],wes_palette("Rushmore1")[2],wes_palette("Rushmore1")[3],"brown4", "cadetblue")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))+  expand_limits(x = 175)
sp1 = sp + geom_point(data=gpd_newdf_means,aes(Ind,GEBV),size=3,alpha=0.9,col=wes_palette("Rushmore1")[3]) +expand_limits(x = 175)
sp1
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/GEBVs_gpd_sortedAfterAvggpd_",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)

# gpdFixcor
gpdFix_newdf = combined_ordered1[which(combined_ordered1$trait=="gpdFixCor"),]
gpdFix_newdf_means = aggregate(as.numeric(gpdFix_newdf$GEBV),list(gpdFix_newdf$Ind), mean)
colnames(gpdFix_newdf_means)=c("Ind","GEBV")

sp<-ggplot(gpdFix_newdf, aes(x= Ind, y=GEBV)) + geom_point(alpha=0.05,size=2,col=wes_palette("Rushmore1")[2])+
  #scale_color_manual(values= c(wes_palette("Rushmore1")[4],wes_palette("Rushmore1")[2],wes_palette("Rushmore1")[3],"brown4", "cadetblue")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))+  expand_limits(x = 175)
sp1 = sp + geom_point(data=gpdFix_newdf_means,aes(Ind,GEBV),size=3,alpha=0.9,col=wes_palette("Rushmore1")[2]) +expand_limits(x = 175)
sp1
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/GEBVs_gpdFixCor_sortedAfterAvggpd_",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)


# gpi
gpi_newdf = combined_ordered1[which(combined_ordered1$trait=="gpi"),]
gpi_newdf_means = aggregate(as.numeric(gpi_newdf$GEBV),list(gpi_newdf$Ind), mean)
colnames(gpi_newdf_means)=c("Ind","GEBV")

sp<-ggplot(gpi_newdf, aes(x= Ind, y=GEBV)) + geom_point(alpha=0.05,size=2,col="brown4")+
  #scale_color_manual(values= c(wes_palette("Rushmore1")[4],wes_palette("Rushmore1")[2],wes_palette("Rushmore1")[3],"brown4", "cadetblue")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))+  expand_limits(x = 175)
sp1 = sp + geom_point(data=gpi_newdf_means,aes(Ind,GEBV),size=3,alpha=0.9,col="brown4") +expand_limits(x = 175)
sp1
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/GEBVs_gpi_sortedAfterAvggpi_",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)

# gpiCor
gpiCor_newdf = combined_ordered1[which(combined_ordered1$trait=="gpiCor"),]
gpiCor_newdf_means = aggregate(as.numeric(gpiCor_newdf$GEBV),list(gpiCor_newdf$Ind), mean)
colnames(gpiCor_newdf_means)=c("Ind","GEBV")

sp<-ggplot(gpiCor_newdf, aes(x= Ind, y=GEBV)) + geom_point(alpha=0.05,size=2,col="cadetblue")+
  #scale_color_manual(values= c(wes_palette("Rushmore1")[4],wes_palette("Rushmore1")[2],wes_palette("Rushmore1")[3],"brown4", "cadetblue")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))+  expand_limits(x = 175)
sp1 = sp + geom_point(data=gpiCor_newdf_means,aes(Ind,GEBV),size=3,alpha=0.9,col="cadetblue") +expand_limits(x = 175)
sp1
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/GEBVs_gpiCor_sortedAfterAvggpi_",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)


