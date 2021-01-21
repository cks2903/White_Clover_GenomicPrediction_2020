
###############################################################
###############################################################
#       A script for investigation the relationships
#       between different traits and their predicted
#       GEBVs
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



# gpdNoCor
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
meangpdNoCor=aggregate(as.numeric(dataset3$Observed), list(dataset3$ID), mean)[,2]
meangebv=aggregate(as.numeric(dataset3$GEBV), list(dataset3$ID), mean)
df_gpdNoCor=cbind(meangebv,meangpdNoCor)
colnames(df_gpdNoCor)=c("Individual","gpdFixCor_GEBV","gpdFixCor_obs")
head(df_gpdNoCor)



# gpi
# Load predictions from gpdFixCor
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
meangpdFixCor=aggregate(as.numeric(dataset4$Observed), list(dataset4$ID), mean)[,2]
meangebv=aggregate(as.numeric(dataset4$GEBV), list(dataset4$ID), mean)
df_gpdFixCor=cbind(meangebv,meangpdFixCor)
colnames(df_gpdFixCor)=c("Individual","gpi_GEBV","gpi_obs")
head(df_gpdFixCor)




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
meangpdDay11to25=aggregate(as.numeric(dataset5$Observed), list(dataset5$ID), mean)[,2]
meangebv=aggregate(as.numeric(dataset5$GEBV), list(dataset5$ID), mean)
df_gpdDay11to25=cbind(meangebv,meangpdDay11to25)
colnames(df_gpdDay11to25)=c("Individual","gpiCor_GEBV","gpiCor_obs")
head(df_gpdDay11to25)



# merge all
Large_df= cbind(df_iSize,df_gpdRescor[,2:3],df_gpdNoCor[,2:3],df_gpdFixCor[,2:3],df_gpdDay11to25[,2:3])
head(Large_df)

Correlations_numeric_matrix = Large_df[,2:ncol(Large_df)]
# correlations



path="/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/"
filename = paste(path,"correlations_",Sys.Date(),".pdf",sep="")

corrgram(Correlations_numeric_matrix,
         main="Correlations between different yield traits and their prediction",
         lower.panel=panel.pts, upper.panel=panel.conf,
         diag.panel=panel.density)




# see if you can predict one trait from another??

# predict gpdResCor from isize GEBVs and gpdDay11to25 GEBVS, partitipation into two growth phases
fit<- lm(gpdResCor_obs ~ iSize_GEBV + gpdDay11to25_GEBV, data = Large_df)
summary(fit)
TotalGEBV=model.matrix( ~ iSize_GEBV + gpdDay11to25_GEBV, data=Large_df) %*% fit$coefficients
plot(Large_df$gpdResCor_obs ~ TotalGEBV)
cor(data.matrix(Large_df$gpdResCor_obs), data.matrix(TotalGEBV)) # 0.40

# predict gpdResCor from isize alone
fit2<- lm(gpdResCor_obs ~iSize_GEBV, data = Large_df)
summary(fit2)
TotalGEBV2=model.matrix( ~ iSize_GEBV, data=Large_df) %*% fit2$coefficients
plot(Large_df$gpdResCor_obs ~ TotalGEBV2)
cor(data.matrix(Large_df$gpdResCor_obs), data.matrix(TotalGEBV2)) # 0.39

# predict gpdResCor from gpdDay11to25 alone
fit3<- lm(gpdResCor_obs ~gpdDay11to25_GEBV, data = Large_df)
summary(fit3)
TotalGEBV3=model.matrix( ~ gpdDay11to25_GEBV, data=Large_df) %*% fit3$coefficients
plot(Large_df$gpdResCor_obs ~ TotalGEBV3)
cor(data.matrix(Large_df$gpdResCor_obs), data.matrix(TotalGEBV3)) # 0.23






