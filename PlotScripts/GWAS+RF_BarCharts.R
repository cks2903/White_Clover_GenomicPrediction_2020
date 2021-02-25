#########################################
#########################################
##  This is a script to display the 
## prediction accuracies achieved
## when combining gwas on a trn pop
## with RF using the top 25 or top 200
## most significant SNPs
## or random SNPs
## traits used for GWAS are iSize and gpd
## traits predicted are iSize or gpd
## there will be 4 combinations
#########################################
#########################################

library(ggplot2)


#############################################
#############################################
##          Prediction of iSize           ##
#############################################
#############################################


setwd("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_iSizeFromiSize")

# Load files with correlations from RF of gpd following GWAS of gpd
list.files()


# top 25 most significant SNPS when predicting from isize
file_list <- list.files(pattern="Correlations_iSize_RF_top25SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationstop25_isizefromisize")){
    correlationstop25_isizefromisize <- read.table(file, header=TRUE, sep=",")

  
  # if the merged dataset does exist, append to it
  }else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationstop25_isizefromisize<-rbind(correlationstop25_isizefromisize, pred_dataset)
    rm(pred_dataset)
  }
}

correlationstop25_isizefromisize=as.data.frame(correlationstop25_isizefromisize)
colnames(correlationstop25_isizefromisize)=c("cor")
nrow(correlationstop25_isizefromisize)==100 #check


# top 200 most significant SNPS when predicting from iSize
file_list <- list.files(pattern="Correlations_iSize_RF_top200SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationstop200_isizefromisize")){
    correlationstop200_isizefromisize <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationstop200_isizefromisize<-rbind(correlationstop200_isizefromisize, pred_dataset)
    rm(pred_dataset)
  }
}

correlationstop200_isizefromisize=as.data.frame(correlationstop200_isizefromisize)
colnames(correlationstop200_isizefromisize)=c("cor")
nrow(correlationstop200_isizefromisize)==100 #check



# top 500 most significant SNPS when predicting from iSize
file_list <- list.files(pattern="Correlations_iSize_RF_top500SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationstop500_isizefromisize")){
    correlationstop500_isizefromisize <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationstop500_isizefromisize<-rbind(correlationstop500_isizefromisize, pred_dataset)
    rm(pred_dataset)
  }
}

correlationstop500_isizefromisize=as.data.frame(correlationstop500_isizefromisize)
colnames(correlationstop500_isizefromisize)=c("cor")
nrow(correlationstop500_isizefromisize)==10 #check


# 25 random SNPS when predicting from iSize
file_list <- list.files(pattern="Correlations_iSize_RF_random25SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationsrandom25_isizefromisize")){
    correlationsrandom25_isizefromisize <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
 else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationsrandom25_isizefromisize<-rbind(correlationsrandom25_isizefromisize, pred_dataset)
    rm(pred_dataset)
  }
}
correlationsrandom25_isizefromisize=as.data.frame(correlationsrandom25_isizefromisize)
colnames(correlationsrandom25_isizefromisize)=c("cor")
nrow(correlationsrandom25_isizefromisize)==10 #check

# 200 random SNPS when predicting from iSize
file_list <- list.files(pattern="Correlations_iSize_RF_random200SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationsrandom200_isizefromisize")){
    correlationsrandom200_isizefromisize <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationsrandom200_isizefromisize<-rbind(correlationsrandom200_isizefromisize, pred_dataset)
    rm(pred_dataset)
  }
}

correlationsrandom200_isizefromisize=as.data.frame(correlationsrandom200_isizefromisize)
colnames(correlationsrandom200_isizefromisize)=c("cor")
nrow(correlationsrandom200_isizefromisize)==10 #check

# 500 random SNPS  when predicting from iSize
file_list <- list.files(pattern="Correlations_iSize_RF_random500SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationsrandom500_isizefromisize")){
    correlationsrandom500_isizefromisize <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationsrandom500_isizefromisize<-rbind(correlationsrandom500_isizefromisize, pred_dataset)
    rm(pred_dataset)
  }
}

correlationsrandom500_isizefromisize=as.data.frame(correlationsrandom500_isizefromisize)
colnames(correlationsrandom500_isizefromisize)=c("cor")
nrow(correlationsrandom500_isizefromisize)==10 #check




setwd("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_iSizeFromgpd")


# top 25 most significant SNPS when predicting from gpd
file_list <- list.files(pattern="Correlations_iSize_RF_top25SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationstop25_isizefromgpd")){
    correlationstop25_isizefromgpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationstop25_isizefromgpd<-rbind(correlationstop25_isizefromgpd, pred_dataset)
    rm(pred_dataset)
  
  }
}
correlationstop25_isizefromgpd=as.data.frame(correlationstop25_isizefromgpd)
colnames(correlationstop25_isizefromgpd)=c("cor")
nrow(correlationstop25_isizefromgpd)==10 #check


# top 200 most significant SNPS when predicting from iSize
file_list <- list.files(pattern="Correlations_iSize_RF_top200SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationstop200isize_isizefromgpd")){
    correlationstop200isize_isizefromgpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationstop200isize_isizefromgpd<-rbind(correlationstop200isize_isizefromgpd, pred_dataset)
    rm(pred_dataset)
  }
}

correlationstop200isize_isizefromgpd=as.data.frame(correlationstop200isize_isizefromgpd)
colnames(correlationstop200isize_isizefromgpd)=c("cor")
nrow(correlationstop200isize_isizefromgpd)==10 #check


# top 500 most significant SNPS when predicting from iSize
file_list <- list.files(pattern="Correlations_iSize_RF_top500SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationstop500_isizefromgpd")){
    correlationstop500_isizefromgpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationstop500_isizefromgpd<-rbind(correlationstop500_isizefromgpd, pred_dataset)
    rm(pred_dataset)
  }
}

correlationstop500_isizefromgpd=as.data.frame(correlationstop500_isizefromgpd)
colnames(correlationstop500_isizefromgpd)=c("cor")
nrow(correlationstop500_isizefromgpd)==10 #check


# 25 random SNPS when predicting from gpd
file_list <- list.files(pattern="Correlations_iSize_RF_random25SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationsrandom25isize_isizefromgpd")){
    correlationsrandom25isize_isizefromgpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationsrandom25isize_isizefromgpd<-rbind(correlationsrandom25isize_isizefromgpd, pred_dataset)
    rm(pred_dataset)
  }
}

correlationsrandom25isize_isizefromgpd=as.data.frame(correlationsrandom25isize_isizefromgpd)
colnames(correlationsrandom25isize_isizefromgpd)=c("cor")
nrow(correlationsrandom25isize_isizefromgpd)==10 #check


# 200 random SNPS when predicting from isize
file_list <- list.files(pattern="Correlations_iSize_RF_random200SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationsrandom200isize_isizefromgpd")){
    correlationsrandom200isize_isizefromgpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationsrandom200isize_isizefromgpd<-rbind(correlationsrandom200isize_isizefromgpd, pred_dataset)
    rm(pred_dataset)
  }
}

correlationsrandom200isize_isizefromgpd=as.data.frame(correlationsrandom200isize_isizefromgpd)
colnames(correlationsrandom200isize_isizefromgpd)=c("cor")
nrow(correlationsrandom200isize_isizefromgpd)==10 #check

# 500 random SNPS  when predicting from gpd
file_list <- list.files(pattern="Correlations_iSize_RF_random500SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationsrandom500_isizefromgpd")){
    correlationsrandom500_isizefromgpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationsrandom500_isizefromgpd<-rbind(correlationsrandom500_isizefromgpd, pred_dataset)
    rm(pred_dataset)
  }
}

correlationsrandom500_isizefromgpd=as.data.frame(correlationsrandom500_isizefromgpd)
colnames(correlationsrandom500_isizefromgpd)=c("cor")
nrow(correlationsrandom500_isizefromgpd)==10 #check








# Collect data in one dataframe for plotting

correlationstop25_isizefromisize_mean=mean(correlationstop25_isizefromisize$cor)
correlationstop25_isizefromisize_sd=sd(correlationstop25_isizefromisize$cor)
correlationstop25_isizefromisize_se=correlationstop25_isizefromisize_sd/(10^0.5)

correlationstop25_isizefromisize_avg=cbind(correlationstop25_isizefromisize_mean,correlationstop25_isizefromisize_sd,correlationstop25_isizefromisize_se)
colnames(correlationstop25_isizefromisize_avg)=c("Cor","SD","SE")

correlationstop200_isizefromisize_mean=mean(correlationstop200_isizefromisize$cor)
correlationstop200_isizefromisize_sd=sd(correlationstop200_isizefromisize$cor)
correlationstop200_isizefromisize_se=correlationstop200_isizefromisize_sd/(10^0.5)

correlationstop200_isizefromisize_avg=cbind(correlationstop200_isizefromisize_mean,correlationstop200_isizefromisize_sd,correlationstop200_isizefromisize_se)
colnames(correlationstop200_isizefromisize_avg)=c("Cor","SD","SE")

correlationstop500_isizefromisize_mean=mean(correlationstop500_isizefromisize$cor)
correlationstop500_isizefromisize_sd=sd(correlationstop500_isizefromisize$cor)
correlationstop500_isizefromisize_se=correlationstop500_isizefromisize_sd/(10^0.5)

correlationstop500_isizefromisize_avg=cbind(correlationstop500_isizefromisize_mean,correlationstop500_isizefromisize_sd,correlationstop500_isizefromisize_se)
colnames(correlationstop500_isizefromisize_avg)=c("Cor","SD","SE")

correlationsrandom25_isizefromisize_mean=mean(correlationsrandom25_isizefromisize$cor)
correlationsrandom25_isizefromisize_sd=sd(correlationsrandom25_isizefromisize$cor)
correlationsrandom25_isizefromisize_se=correlationsrandom25_isizefromisize_sd/(10^0.5)

correlationsrandom25_isizefromisize_avg=cbind(correlationsrandom25_isizefromisize_mean,correlationsrandom25_isizefromisize_sd,correlationsrandom25_isizefromisize_se)
colnames(correlationsrandom25_isizefromisize_avg)=c("Cor","SD","SE")

correlationsrandom200_isizefromisize_mean=mean(correlationsrandom200_isizefromisize$cor)
correlationsrandom200_isizefromisize_sd=sd(correlationsrandom200_isizefromisize$cor)
correlationsrandom200_isizefromisize_se=correlationsrandom200_isizefromisize_sd/(10^0.5)

correlationsrandom200_isizefromisize_avg=cbind(correlationsrandom200_isizefromisize_mean,correlationsrandom200_isizefromisize_sd,correlationsrandom200_isizefromisize_se)
colnames(correlationsrandom200_isizefromisize_avg)=c("Cor","SD","SE")

correlationsrandom500_isizefromisize_mean=mean(correlationsrandom500_isizefromisize$cor)
correlationsrandom500_isizefromisize_sd=sd(correlationsrandom500_isizefromisize$cor)
correlationsrandom500_isizefromisize_se=correlationsrandom500_isizefromisize_sd/(10^0.5)

correlationsrandom500_isizefromisize_avg=cbind(correlationsrandom500_isizefromisize_mean,correlationsrandom500_isizefromisize_sd,correlationsrandom500_isizefromisize_se)
colnames(correlationsrandom500_isizefromisize_avg)=c("Cor","SD","SE")


correlationstop25_isizefromgpd_mean=mean(correlationstop25_isizefromgpd$cor)
correlationstop25_isizefromgpd_sd=sd(correlationstop25_isizefromgpd$cor)
correlationstop25_isizefromgpd_se=correlationstop25_isizefromgpd_sd/(10^0.5)

correlationstop25_isizefromgpd_avg=cbind(correlationstop25_isizefromgpd_mean,correlationstop25_isizefromgpd_sd,correlationstop25_isizefromgpd_se)
colnames(correlationstop25_isizefromgpd_avg)=c("Cor","SD","SE")

correlationstop200_isizefromgpd_mean=mean(correlationstop200isize_isizefromgpd$cor)
correlationstop200_isizefromgpd_sd=sd(correlationstop200isize_isizefromgpd$cor)
correlationstop200_isizefromgpd_se=correlationstop200_isizefromgpd_sd/(10^0.5)

correlationstop200_isizefromgpd_avg=cbind(correlationstop200_isizefromgpd_mean,correlationstop200_isizefromgpd_sd,correlationstop200_isizefromgpd_se)
colnames(correlationstop200_isizefromgpd_avg)=c("Cor","SD","SE")

correlationstop500_isizefromgpd_mean=mean(correlationstop500_isizefromgpd$cor)
correlationstop500_isizefromgpd_sd=sd(correlationstop500_isizefromgpd$cor)
correlationstop500_isizefromgpd_se=correlationstop500_isizefromgpd_sd/(10^0.5)

correlationstop500_isizefromgpd_avg=cbind(correlationstop500_isizefromgpd_mean,correlationstop500_isizefromgpd_sd,correlationstop500_isizefromgpd_se)
colnames(correlationstop500_isizefromgpd_avg)=c("Cor","SD","SE")

correlationsrandom25_isizefromgpd_mean=mean(correlationsrandom25isize_isizefromgpd$cor)
correlationsrandom25_isizefromgpd_sd=sd(correlationsrandom25isize_isizefromgpd$cor)
correlationsrandom25_isizefromgpd_se=correlationsrandom25_isizefromgpd_sd/(10^0.5)

correlationsrandom25_isizefromgpd_avg=cbind(correlationsrandom25_isizefromgpd_mean,correlationsrandom25_isizefromgpd_sd,correlationsrandom25_isizefromgpd_se)
colnames(correlationsrandom25_isizefromgpd_avg)=c("Cor","SD","SE")

correlationsrandom200_isizefromgpd_mean=mean(correlationsrandom200isize_isizefromgpd$cor)
correlationsrandom200_isizefromgpd_sd=sd(correlationsrandom200isize_isizefromgpd$cor)
correlationsrandom200_isizefromgpd_se=correlationsrandom200_isizefromgpd_sd/(10^0.5)
correlationsrandom200_isizefromgpd_avg=cbind(correlationsrandom200_isizefromgpd_mean,correlationsrandom200_isizefromgpd_sd,correlationsrandom200_isizefromgpd_se)
colnames(correlationsrandom200_isizefromgpd_avg)=c("Cor","SD","SE")

correlationsrandom500_isizefromgpd_mean=mean(correlationsrandom500_isizefromgpd$cor)
correlationsrandom500_isizefromgpd_sd=sd(correlationsrandom500_isizefromgpd$cor)
correlationsrandom500_isizefromgpd_se=correlationsrandom500_isizefromgpd_sd/(10^0.5)

correlationsrandom500_isizefromgpd_avg=cbind(correlationsrandom500_isizefromgpd_mean,correlationsrandom500_isizefromgpd_sd,correlationsrandom500_isizefromgpd_se)
colnames(correlationsrandom500_isizefromgpd_avg)=c("Cor","SD","SE")

d=rbind(correlationstop25_isizefromisize_avg,correlationstop200_isizefromisize_avg,correlationstop500_isizefromisize_avg,correlationsrandom25_isizefromisize_avg,correlationsrandom200_isizefromisize_avg,correlationsrandom500_isizefromisize_avg,correlationstop25_isizefromgpd_avg,correlationstop200_isizefromgpd_avg,correlationstop500_isizefromgpd_avg,correlationsrandom25_isizefromgpd_avg,correlationsrandom200_isizefromgpd_avg,correlationsrandom500_isizefromgpd_avg)
d=as.data.frame(d)
nrow(d)==10*2*1 #check

d$method= rep("method",nrow(d))
d$method[1]= "25 SNPs (iSize)"
d$method[2]= "200 SNPs (iSize)"
d$method[3]= "500 SNPs (iSize)"
d$method[4]= "25 SNPs (iSize)"
d$method[5]= "200 SNPs (iSize)"
d$method[6]= "500 SNPs (iSize)"
d$method[7]= "25 SNPs (gpd)"
d$method[8]= "200 SNPs (gpd)"
d$method[9]= "500 SNPs (gpd)"
d$method[10]= "25 SNPs (gpd)"
d$method[11]= "200 SNPs (gpd)"
d$method[12]= "500 SNPs (gpd)"

d$snps= rep("snps",nrow(d))
d$snps[1]= "Top"
d$snps[2]= "Top"
d$snps[3]= "Top"
d$snps[4]= "Random"
d$snps[5]= "Random"
d$snps[6]= "Random"
d$snps[7]= "Top"
d$snps[8]= "Top"
d$snps[9]= "Top"
d$snps[10]= "Random"
d$snps[11]= "Random"
d$snps[12]= "Random"

head(d)
str(d)



# load GBLUP results as well
#GBLUPcor_iSize= read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/iSize_20210120/AllCorrelations.txt",sep="\t")
setwd("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/PlottingScripts_20210215")
GBLUPcor_iSize= read.table("AllCorrelationsGBLUPiSize.txt",sep="\t")
head(GBLUPcor_iSize)
GBLUP_iSize_avgCor=mean(GBLUPcor_iSize[,1])
GBLUP_iSize_SD=sd(GBLUPcor_iSize[,1])
GBLUP_iSize_SE=GBLUP_iSize_SD/(100^0.5)



# plotting
vector= c(GBLUP_iSize_avgCor,GBLUP_iSize_SD,GBLUP_iSize_SE,"GBLUP","All")
d1=rbind(d,vector)


d1$snps <- factor(d1$snps,levels = c("Top","Random","All"))
d1$method <- factor(d1$method,levels = c("25 SNPs (iSize)","200 SNPs (iSize)","500 SNPs (iSize)","25 SNPs (gpd)","200 SNPs (gpd)","500 SNPs (gpd)","GBLUP"))
d1$Cor=as.numeric(d1$Cor)
d1$SE=as.numeric(d1$SE)


p <- ggplot(data=d1, aes(x=method, y=Cor, fill=snps,width=c(rep(0.9,12),0.45))) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=Cor-SE, ymax=Cor+SE), width=.2,
                position=position_dodge(.9)) 

p1= p + scale_fill_manual(values=c('#CA626C','#ABABAB',"#CA626C")) + ggtitle("Prediction of iSize") +
   ylab("Correlation") +xlab('') +ylim(c(-0.01,0.5))+theme( 
     #legend.position="none",
     #panel.border = element_blank(),
     panel.background = element_blank(),
     panel.grid.major.x = element_blank(),
     panel.grid.minor.x = element_blank(),
     panel.grid.minor.y = element_blank(),
     axis.line = element_line(colour = "black"),
     axis.text.x = element_text(size = 15),
     axis.text.y = element_text(size = 15),
     axis.title = element_text(size = 22),
     title = element_text(size = 20))

ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/Correlations_predictionIsize',Sys.Date(),'.pdf',sep="_"), plot = p1, width = 45, height = 15, unit = 'cm')

dsubset=d1[-c(3,6,9,12),]

p <- ggplot(data=dsubset, aes(x=method, y=Cor, fill=snps,width=c(rep(0.9,8),0.45))) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=Cor-SE, ymax=Cor+SE), width=.2,
                position=position_dodge(.9)) 

p1= p + scale_fill_manual(values=c('#CA626C','#ABABAB','#CA626C')) + ggtitle("Prediction of iSize") +
  ylab("Correlation") +xlab('') +ylim(c(-0.01,0.5)) +theme( 
    #legend.position="none",
    #panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 22),
    title = element_text(size = 20))

ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/Correlations_predictionIsizev2',Sys.Date(),'.pdf',sep="_"), plot = p1, width = 30, height = 15, unit = 'cm')




#############################################
#############################################
##          Prediction of gpd              ##
#############################################
#############################################


setwd("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_gpdFromiSize")

# Load files with correlations from RF of gpd following GWAS of gpd
list.files()


# top 25 most significant SNPS when predicting from isize
file_list <- list.files(pattern="Correlations_gpd_RF_top25SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationstop25_gpdfromisize")){
    correlationstop25_gpdfromisize <- read.table(file, header=TRUE, sep=",")
    
    
    # if the merged dataset does exist, append to it
  }else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationstop25_gpdfromisize<-rbind(correlationstop25_gpdfromisize, pred_dataset)
    rm(pred_dataset)
  }
}

correlationstop25_gpdfromisize=as.data.frame(correlationstop25_gpdfromisize)
colnames(correlationstop25_gpdfromisize)=c("cor")
nrow(correlationstop25_gpdfromisize)==10

# top 200 most significant SNPS when predicting from iSize
file_list <- list.files(pattern="Correlations_gpd_RF_top200SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationstop200_gpdfromisize")){
    correlationstop200_gpdfromisize <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationstop200_gpdfromisize<-rbind(correlationstop200_gpdfromisize, pred_dataset)
    rm(pred_dataset)
  }
}

correlationstop200_gpdfromisize=as.data.frame(correlationstop200_gpdfromisize)
colnames(correlationstop200_gpdfromisize)=c("cor")
nrow(correlationstop200_gpdfromisize)==10


# top 500 most significant SNPS when predicting from iSize
file_list <- list.files(pattern="Correlations_gpd_RF_top500SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationstop500gpdfromisize")){
    correlationstop500gpdfromisize <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationstop500gpdfromisize<-rbind(correlationstop500gpdfromisize, pred_dataset)
    rm(pred_dataset)
  }
}

correlationstop500gpdfromisize=as.data.frame(correlationstop500gpdfromisize)
colnames(correlationstop500gpdfromisize)=c("cor")
nrow(correlationstop500gpdfromisize)==10


# 25 random SNPS when predicting from iSize
file_list <- list.files(pattern="Correlations_gpd_RF_random25SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationsrandom25_gpdfromisize")){
    correlationsrandom25_gpdfromisize <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationsrandom25_gpdfromisize<-rbind(correlationsrandom25_gpdfromisize, pred_dataset)
    rm(pred_dataset)
  }
}
correlationsrandom25_gpdfromisize=as.data.frame(correlationsrandom25_gpdfromisize)
colnames(correlationsrandom25_gpdfromisize)=c("cor")
nrow(correlationsrandom25_gpdfromisize)==10


# 200 random SNPS when predicting from iSize
file_list <- list.files(pattern="Correlations_gpd_RF_random200SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationsrandom200_gpdfromisize")){
    correlationsrandom200_gpdfromisize <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationsrandom200_gpdfromisize<-rbind(correlationsrandom200_gpdfromisize, pred_dataset)
    rm(pred_dataset)
  }
}

correlationsrandom200_gpdfromisize=as.data.frame(correlationsrandom200_gpdfromisize)
colnames(correlationsrandom200_gpdfromisize)=c("cor")
nrow(correlationsrandom200_gpdfromisize)==10


# 500 random SNPS  when predicting from iSize
file_list <- list.files(pattern="Correlations_gpd_RF_random500SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationsrandom500_gpdfromisize")){
    correlationsrandom500_gpdfromisize <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationsrandom500_gpdfromisize<-rbind(correlationsrandom500_gpdfromisize, pred_dataset)
    rm(pred_dataset)
  }
}

correlationsrandom500_gpdfromisize=as.data.frame(correlationsrandom500_gpdfromisize)
colnames(correlationsrandom500_gpdfromisize)=c("cor")
nrow(correlationsrandom500_gpdfromisize)==10





setwd("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_gpdFromgpd")


# top 25 most significant SNPS when predicting from gpd
file_list <- list.files(pattern="Correlations_gpd_RF_top25SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationstop25_gpdfromgpd")){
    correlationstop25_gpdfromgpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationstop25_gpdfromgpd<-rbind(correlationstop25_gpdfromgpd, pred_dataset)
    rm(pred_dataset)
    
  }
}
correlationstop25_gpdfromgpd=as.data.frame(correlationstop25_gpdfromgpd)
colnames(correlationstop25_gpdfromgpd)=c("cor")
nrow(correlationstop25_gpdfromgpd)==10


# top 200 most significant SNPS when predicting from iSize
file_list <- list.files(pattern="Correlations_gpd_RF_top200SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationstop200isize_gpdfromgpd")){
    correlationstop200isize_gpdfromgpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationstop200isize_gpdfromgpd<-rbind(correlationstop200isize_gpdfromgpd, pred_dataset)
    rm(pred_dataset)
  }
}

correlationstop200isize_gpdfromgpd=as.data.frame(correlationstop200isize_gpdfromgpd)
colnames(correlationstop200isize_gpdfromgpd)=c("cor")
nrow(correlationstop200isize_gpdfromgpd)==10


# top 500 most significant SNPS when predicting from iSize
file_list <- list.files(pattern="Correlations_gpd_RF_top500SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationstop500_gpdfromgpd")){
    correlationstop500_gpdfromgpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationstop500_gpdfromgpd<-rbind(correlationstop500_gpdfromgpd, pred_dataset)
    rm(pred_dataset)
  }
}

correlationstop500_gpdfromgpd=as.data.frame(correlationstop500_gpdfromgpd)
colnames(correlationstop500_gpdfromgpd)=c("cor")
nrow(correlationstop500_gpdfromgpd)==10


# 25 random SNPS when predicting from gpd
file_list <- list.files(pattern="Correlations_gpd_RF_random25SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationsrandom25isize_gpdfromgpd")){
    correlationsrandom25isize_gpdfromgpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationsrandom25isize_gpdfromgpd<-rbind(correlationsrandom25isize_gpdfromgpd, pred_dataset)
    rm(pred_dataset)
  }
}

correlationsrandom25isize_gpdfromgpd=as.data.frame(correlationsrandom25isize_gpdfromgpd)
colnames(correlationsrandom25isize_gpdfromgpd)=c("cor")
nrow(correlationsrandom25isize_gpdfromgpd)==10


# 200 random SNPS when predicting from isize
file_list <- list.files(pattern="Correlations_gpd_RF_random200SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationsrandom200gpdfromgpd")){
    correlationsrandom200gpdfromgpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationsrandom200gpdfromgpd<-rbind(correlationsrandom200gpdfromgpd, pred_dataset)
    rm(pred_dataset)
  }
}

correlationsrandom200gpdfromgpd=as.data.frame(correlationsrandom200gpdfromgpd)
colnames(correlationsrandom200gpdfromgpd)=c("cor")
nrow(correlationsrandom200gpdfromgpd)==10


# 500 random SNPS  when predicting from gpd
file_list <- list.files(pattern="Correlations_gpd_RF_random500SNPs_grouping")
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("correlationsrandom500_gpdfromgpd")){
    correlationsrandom500_gpdfromgpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  else{
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    correlationsrandom500_gpdfromgpd<-rbind(correlationsrandom500_gpdfromgpd, pred_dataset)
    rm(pred_dataset)
  }
}

correlationsrandom500_gpdfromgpd=as.data.frame(correlationsrandom500_gpdfromgpd)
colnames(correlationsrandom500_gpdfromgpd)=c("cor")
nrow(correlationsrandom500_gpdfromgpd)==10








# Collect data in one dataframe for plotting

correlationstop25_gpdfromisize_mean=mean(correlationstop25_gpdfromisize$cor)
correlationstop25_gpdfromisize_sd=sd(correlationstop25_gpdfromisize$cor)
correlationstop25_gpdfromisize_se=correlationstop25_gpdfromisize_sd/(10^0.5)
correlationstop25_gpdfromisize_avg=cbind(correlationstop25_gpdfromisize_mean,correlationstop25_gpdfromisize_sd,correlationstop25_gpdfromisize_se)
colnames(correlationstop25_gpdfromisize_avg)=c("Cor","SD","SE")

correlationstop200_gpdfromisize_mean=mean(correlationstop200_gpdfromisize$cor)
correlationstop200_gpdfromisize_sd=sd(correlationstop200_gpdfromisize$cor)
correlationstop200_gpdfromisize_se=correlationstop200_gpdfromisize_sd/(10^0.5)
correlationstop200_gpdfromisize_avg=cbind(correlationstop200_gpdfromisize_mean,correlationstop200_gpdfromisize_sd,correlationstop200_gpdfromisize_se)
colnames(correlationstop200_gpdfromisize_avg)=c("Cor","SD","SE")

correlationstop500_gpdfromisize_mean=mean(correlationstop500gpdfromisize$cor)
correlationstop500_gpdfromisize_sd=sd(correlationstop500gpdfromisize$cor)
correlationstop200_gpdfromisize_se=correlationstop500_gpdfromisize_sd/(10^0.5)
correlationstop500_gpdfromisize_avg=cbind(correlationstop500_gpdfromisize_mean,correlationstop500_gpdfromisize_sd,correlationstop200_gpdfromisize_se)
colnames(correlationstop500_gpdfromisize_avg)=c("Cor","SD","SE")

correlationsrandom25_gpdfromisize_mean=mean(correlationsrandom25_gpdfromisize$cor)
correlationsrandom25_gpdfromisize_sd=sd(correlationsrandom25_gpdfromisize$cor)
correlationsrandom25_gpdfromisize_se=correlationsrandom25_gpdfromisize_sd/(10^0.5)
correlationsrandom25_gpdfromisize_avg=cbind(correlationsrandom25_gpdfromisize_mean,correlationsrandom25_gpdfromisize_sd,correlationsrandom25_gpdfromisize_se)
colnames(correlationsrandom25_gpdfromisize_avg)=c("Cor","SD","SE")

correlationsrandom200_gpdfromisize_mean=mean(correlationsrandom200_gpdfromisize$cor)
correlationsrandom200_gpdfromisize_sd=sd(correlationsrandom200_gpdfromisize$cor)
correlationsrandom200_gpdfromisize_se=correlationsrandom200_gpdfromisize_sd/(10^0.5)
correlationsrandom200_gpdfromisize_avg=cbind(correlationsrandom200_gpdfromisize_mean,correlationsrandom200_gpdfromisize_sd,correlationsrandom200_gpdfromisize_se)
colnames(correlationsrandom200_gpdfromisize_avg)=c("Cor","SD","SE")

correlationsrandom500_gpdfromisize_mean=mean(correlationsrandom500_gpdfromisize$cor)
correlationsrandom500_gpdfromisize_sd=sd(correlationsrandom500_gpdfromisize$cor)
correlationsrandom500_gpdfromisize_se=correlationsrandom500_gpdfromisize_sd/(10^0.5)
correlationsrandom500_gpdfromisize_avg=cbind(correlationsrandom500_gpdfromisize_mean,correlationsrandom500_gpdfromisize_sd,correlationsrandom500_gpdfromisize_se)
colnames(correlationsrandom500_gpdfromisize_avg)=c("Cor","SD","SE")


correlationstop25_gpdfromgpd_mean=mean(correlationstop25_gpdfromgpd$cor)
correlationstop25_gpdfromgpd_sd=sd(correlationstop25_gpdfromgpd$cor)
correlationstop25_gpdfromgpd_se=correlationstop25_gpdfromgpd_sd/(10^0.5)
correlationstop25_gpdfromgpd_avg=cbind(correlationstop25_gpdfromgpd_mean,correlationstop25_gpdfromgpd_sd,correlationstop25_gpdfromgpd_se)
colnames(correlationstop25_gpdfromgpd_avg)=c("Cor","SD","SE")

correlationstop200_gpdfromgpd_mean=mean(correlationstop200isize_gpdfromgpd$cor)
correlationstop200_gpdfromgpd_sd=sd(correlationstop200isize_gpdfromgpd$cor)
correlationstop200_gpdfromgpd_se=correlationstop200_gpdfromgpd_sd/(10^0.5)
correlationstop200_gpdfromgpd_avg=cbind(correlationstop200_gpdfromgpd_mean,correlationstop200_gpdfromgpd_sd,correlationstop200_gpdfromgpd_se)
colnames(correlationstop200_gpdfromgpd_avg)=c("Cor","SD","SE")

correlationstop500_gpdfromgpd_mean=mean(correlationstop500_gpdfromgpd$cor)
correlationstop500_gpdfromgpd_sd=sd(correlationstop500_gpdfromgpd$cor)
correlationstop500_gpdfromgpd_se=correlationstop500_gpdfromgpd_sd/(10^0.5)
correlationstop500_gpdfromgpd_avg=cbind(correlationstop500_gpdfromgpd_mean,correlationstop500_gpdfromgpd_sd,correlationstop500_gpdfromgpd_se)
colnames(correlationstop500_gpdfromgpd_avg)=c("Cor","SD","SE")

correlationsrandom25_gpdfromgpd_mean=mean(correlationsrandom25isize_gpdfromgpd$cor)
correlationsrandom25_gpdfromgpd_sd=sd(correlationsrandom25isize_gpdfromgpd$cor)
correlationsrandom25_gpdfromgpd_se=correlationsrandom25_gpdfromgpd_sd/(10^0.5)
correlationsrandom25_gpdfromgpd_avg=cbind(correlationsrandom25_gpdfromgpd_mean,correlationsrandom25_gpdfromgpd_sd,correlationsrandom25_gpdfromgpd_se)
colnames(correlationsrandom25_gpdfromgpd_avg)=c("Cor","SD","SE")

correlationsrandom200_gpdfromgpd_mean=mean(correlationsrandom200gpdfromgpd$cor)
correlationsrandom200_gpdfromgpd_sd=sd(correlationsrandom200gpdfromgpd$cor)
correlationsrandom200_gpdfromgpd_se=correlationsrandom200_gpdfromgpd_sd/(10^0.5)
correlationsrandom200_gpdfromgpd_avg=cbind(correlationsrandom200_gpdfromgpd_mean,correlationsrandom200_gpdfromgpd_sd,correlationsrandom200_gpdfromgpd_se)
colnames(correlationsrandom200_gpdfromgpd_avg)=c("Cor","SD","SE")

correlationsrandom500_gpdfromgpd_mean=mean(correlationsrandom500_gpdfromgpd$cor)
correlationsrandom500_gpdfromgpd_sd=sd(correlationsrandom500_gpdfromgpd$cor)
correlationsrandom500_gpdfromgpd_se=correlationsrandom500_gpdfromgpd_sd/(10^0.5)
correlationsrandom500_gpdfromgpd_avg=cbind(correlationsrandom500_gpdfromgpd_mean,correlationsrandom500_gpdfromgpd_sd,correlationsrandom500_gpdfromgpd_se)
colnames(correlationsrandom500_gpdfromgpd_avg)=c("Cor","SD","SE")
d=rbind(correlationstop25_gpdfromisize_avg,correlationstop200_gpdfromisize_avg,correlationstop500_gpdfromisize_avg,correlationsrandom25_gpdfromisize_avg,correlationsrandom200_gpdfromisize_avg,correlationsrandom500_gpdfromisize_avg,correlationstop25_gpdfromgpd_avg,correlationstop200_gpdfromgpd_avg,correlationstop500_gpdfromgpd_avg,correlationsrandom25_gpdfromgpd_avg,correlationsrandom200_gpdfromgpd_avg,correlationsrandom500_gpdfromgpd_avg)
d=as.data.frame(d)
nrow(d)==10*2*1 #check

d$method= rep("method",nrow(d))
d$method[1]= "25 SNPs (iSize)"
d$method[2]= "200 SNPs (iSize)"
d$method[3]= "500 SNPs (iSize)"
d$method[4]= "25 SNPs (iSize)"
d$method[5]= "200 SNPs (iSize)"
d$method[6]= "500 SNPs (iSize)"
d$method[7]= "25 SNPs (gpd)"
d$method[8]= "200 SNPs (gpd)"
d$method[9]= "500 SNPs (gpd)"
d$method[10]= "25 SNPs (gpd)"
d$method[11]= "200 SNPs (gpd)"
d$method[12]= "500 SNPs (gpd)"

d$snps= rep("snps",nrow(d))
d$snps[1]= "Top"
d$snps[2]= "Top"
d$snps[3]= "Top"
d$snps[4]= "Random"
d$snps[5]= "Random"
d$snps[6]= "Random"
d$snps[7]= "Top"
d$snps[8]= "Top"
d$snps[9]= "Top"
d$snps[10]= "Random"
d$snps[11]= "Random"
d$snps[12]= "Random"

head(d)
str(d)


# load GBLUP results as well
#GBLUPcor_gpd= read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/gpd_20210120/AllCorrelations.txt",sep="\t")
setwd("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/PlottingScripts_20210215")
GBLUPcor_gpd= read.table("AllCorrelationsGBLUPgpd.txt",sep="\t")
head(GBLUPcor_gpd)
GBLUP_gpd_avgCor=mean(GBLUPcor_gpd[,1])
GBLUP_gpd_SD=sd(GBLUPcor_gpd[,1])
GBLUP_gpd_SE=GBLUP_iSize_SD/(100^0.5)



# plotting
vector= c(GBLUP_gpd_avgCor,GBLUP_gpd_SD,GBLUP_gpd_SE,"GBLUP","All")
d1=rbind(d,vector)

# plotting

d1$snps <- factor(d1$snps,levels = c("Top","Random","All"))
d1$method <- factor(d1$method,levels = c("25 SNPs (iSize)","200 SNPs (iSize)","500 SNPs (iSize)","25 SNPs (gpd)","200 SNPs (gpd)","500 SNPs (gpd)","GBLUP"))
d1$Cor=as.numeric(d1$Cor)
d1$SE=as.numeric(d1$SE)

p <- ggplot(data=d1, aes(x=method, y=Cor, fill=snps,width=c(rep(0.9,12),0.45))) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=Cor-SE, ymax=Cor+SE), width=.2,
                position=position_dodge(.9)) 

p1= p + scale_fill_manual(values=c('#CA626C','#ABABAB','#CA626C')) + ggtitle("Prediction of gpd") +
  ylab("Correlation") +xlab('') +ylim(c(-0.01,0.5))+theme( 
    #legend.position="none",
    #panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 22),
    title = element_text(size = 20))

ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/Correlations_predictiongpd',Sys.Date(),'.pdf',sep="_"), plot = p1, width = 30, height = 15, unit = 'cm')


dsubset=d1[-c(3,6,9,12),]

p <- ggplot(data=dsubset, aes(x=method, y=Cor, fill=snps,width=c(rep(0.9,8),0.45))) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=Cor-SE, ymax=Cor+SE), width=.2,
                position=position_dodge(.9)) 

p1= p + scale_fill_manual(values=c('#CA626C','#ABABAB','#CA626C')) + ggtitle("Prediction of gpd") +
  ylab("Correlation") +xlab('') +ylim(c(-0.01,0.5)) +theme( 
    #legend.position="none",
    #panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 22),
    title = element_text(size = 20))

ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/Correlations_predictiongpdv2',Sys.Date(),'.pdf',sep="_"), plot = p1, width = 30, height = 15, unit = 'cm')
