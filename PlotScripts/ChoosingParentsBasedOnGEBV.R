library(ggplot2)
library("dplyr")
library(tidyverse)
library(wesanderson)
library(gghighlight)
library(gridExtra)
library(agricolae)
library(corrgram)
library(ggrepel)
library(psych)

# Load predictions_gpd
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/gpd_20210120")

# Load files

list.files()


file_list <- list.files(pattern="Predictions_")

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    pred_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset<-rbind(dataset, pred_dataset)
    rm(pred_dataset)
  }
  
}

nrow(dataset)
dataset=as.data.frame(dataset)
head(dataset)

# Calculate mean for each individual
meangpd=aggregate(as.numeric(dataset$Observed), list(dataset$ID), mean)[,2]
meangebv=aggregate(as.numeric(dataset$GEBV), list(dataset$ID), mean)[,2]
ind=aggregate(as.numeric(dataset$GEBV), list(dataset$ID), mean)[,1]
df_gpd=cbind(as.character(ind),meangpd,meangebv)
colnames(df_gpd)=c("Individual","Observed","GEBVs")
head(df_gpd)




# Load predictions_gpd_FixCor
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/gpdFixCor_afteraveraging_20210121")
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

# Calculate mean for each individual
meangpdFix=aggregate(as.numeric(dataset2$Observed), list(dataset2$ID), mean)[,2]
meangebv=aggregate(as.numeric(dataset2$GEBV), list(dataset2$ID), mean)[,2]
ind=aggregate(as.numeric(dataset2$GEBV), list(dataset2$ID), mean)[,1]

df_gpdFixCor=cbind(as.character(ind),meangpdFix,meangebv)
colnames(df_gpdFixCor)=c("Individual","Observed","GEBVs")
head(df_gpdFixCor)



# Load predictions_gpi
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/gpdDay11to25/")
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

# Calculate mean for each individual
meangpi=aggregate(as.numeric(dataset3$Observed), list(dataset3$ID), mean)[,2]
ind=aggregate(as.numeric(dataset3$Observed), list(dataset3$ID), mean)[,1]

meangebv=aggregate(as.numeric(dataset3$GEBV), list(dataset3$ID), mean)[,2]
df_gpi=cbind(as.character(ind),meangpi,meangebv)
colnames(df_gpi)=c("Individual","Observed","GEBVs")
head(df_gpi)




# Load predictions_gpiCor
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/gpdDay11to25_correctedForiSize/")
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

# Calculate mean for each individual
meangpicor=aggregate(as.numeric(dataset4$Observed), list(dataset4$ID), mean)[,2]
ind=aggregate(as.numeric(dataset4$Observed), list(dataset4$ID), mean)[,1]
meangebv=aggregate(as.numeric(dataset4$GEBV), list(dataset4$ID), mean)[,2]
df_gpicor=cbind(as.character(ind),meangpicor,meangebv)
colnames(df_gpicor)=c("Individual","Observed","GEBVs")
head(df_gpicor)





# Load predictions iSize
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_averages_20201120/iSize_20210120/")
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

# Calculate mean for each individual
meanisize=aggregate(as.numeric(dataset5$Observed), list(dataset5$ID), mean)[,2]
ind=aggregate(as.numeric(dataset5$Observed), list(dataset5$ID), mean)[,1]
meangebv=aggregate(as.numeric(dataset5$GEBV), list(dataset5$ID), mean)[,2]
df_isize=cbind(as.character(ind),meanisize,meangebv)
colnames(df_isize)=c("Individual","Observed","GEBVs")
head(df_isize)




# plotting


setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter")


# combine in one plot
combined=rbind(dataset,dataset2,dataset3,dataset4,dataset5)
combined$gpd_measurement=rep("NA",nrow(combined))
combined$gpd_measurement[1:nrow(dataset)]="gpd"
combined$gpd_measurement[(nrow(dataset)+1):(nrow(dataset)+nrow(dataset2))]="gpd_FixCor"
combined$gpd_measurement[(nrow(dataset)+nrow(dataset2)+1):(nrow(dataset)+nrow(dataset2)+nrow(dataset3))]="gpi"
combined$gpd_measurement[(nrow(dataset)+nrow(dataset2)+nrow(dataset3)+1):(nrow(dataset2)+nrow(dataset3)+nrow(dataset4))]="gpiCor"
combined$gpd_measurement[(nrow(dataset2)+nrow(dataset3)+nrow(dataset4)+1):nrow(combined)]="iSize"



# Not used in manuscript
png(file=paste("Comparison_of_GEBVs_from_different_gpdMeasurements",Sys.Date(),".png",sep="_"),width=150,height=60,units="cm",res=300, pointsize=6)

e <- ggplot(combined, aes(x = ID, y = GEBV))

e2 <- e + geom_boxplot(
  aes(fill = gpd_measurement),
  position = position_dodge(0.9) 
) +
  scale_fill_manual(values = c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[2],"brown4", "cadetblue",wes_palette("Rushmore1")[4]))+
  theme(axis.text.x = element_text(angle = 90)) 
e2

dev.off() # So GEBVs of gpd_rescor is not more precise



# Not used in manuscript
png(file=paste("Comparison_of_GEBVs_from_different_gpdMeasurements_2",Sys.Date(),".png",sep="_"),width=80,height=40,units="cm",res=300, pointsize=6)

sp<-ggplot(combined, aes(x=ID, y=GEBV, color=gpd_measurement,alpha=0.1)) + geom_point()+theme(axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values = c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[2],"brown4", "cadetblue",wes_palette("Rushmore1")[4])) +
  theme(axis.text.x = element_text(angle = 90)) 

  
sp
dev.off()


mean_merged=rbind(df_gpd,df_gpdFixCor,df_gpi,df_gpicor,df_isize)
mean_merged= as.data.frame(mean_merged)
mean_merged$gpd_measurement=rep("NA",nrow(mean_merged))
mean_merged$gpd_measurement[1:145]="gpd"
mean_merged$gpd_measurement[146:290]="gpd_FixCor"
mean_merged$gpd_measurement[291:435]="gpi"
mean_merged$gpd_measurement[436:580]="gpiCor"
mean_merged$gpd_measurement[581:725]="iSize"

# Not used in manuscript
mean_merged$GEBVs=as.numeric(as.character(mean_merged$GEBVs))
mean_merged$Observed=as.numeric(as.character(mean_merged$Observed))

png(file=paste("Comparison_of_GEBVs_from_different_gpdMeasurements_3",Sys.Date(),".png",sep="_"),width=80,height=40,units="cm",res=300, pointsize=10)
sp<-ggplot(mean_merged, aes(x=Individual, y=GEBVs, color=gpd_measurement)) + geom_point() +
  scale_color_manual(values = c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[2],"brown4", "cadetblue",wes_palette("Rushmore1")[4])) 
sp
dev.off()

# Not used in manuscript
png(file=paste("Comparison_of_GEBVs_from_different_gpdMeasurements_4",Sys.Date(),".png",sep="_"),width=80,height=40,units="cm",res=300, pointsize=6)
sp2<-ggplot(mean_merged, aes(x=Individual, y=GEBVs, color=gpd_measurement)) + geom_point() + geom_line(aes(group=gpd_measurement)) +theme(axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values = c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[2],"brown4", "cadetblue",wes_palette("Rushmore1")[4])) 

sp2
dev.off()



# Sorted after avg GEBVs of gpd, Not used in manuscript
png(file=paste("Comparison_of_GEBVs_from_different_gpdMeasurements_sortedAfterAvg",Sys.Date(),".png",sep="_"),width=70,height=60,units="cm",res=300)

df_gpd=as.data.frame(df_gpd)
df_gpd$Observed=as.numeric(as.character(df_gpd$Observed))
df_gpd$GEBVs=as.numeric(as.character(df_gpd$GEBVs))

df_ordered <- df_gpd[order(df_gpd$GEBVs),]
individual_order=df_ordered$Individual

keyDF <- data.frame(key=individual_order,weight=1:length(individual_order))
colnames(combined)[1]="Individual"
combined_ready <- merge(combined,keyDF,by.x='Individual',by.y='key',all.x=T,all.y=F)
combined_ordered <- combined_ready[order(combined_ready$weight),c('Individual','Observed','GEBV','gpd_measurement', 'weight')]

combined_ordered$Individual <- factor(combined_ordered$Individual,levels = individual_order)

sp<-ggplot(combined_ordered, aes(x= Individual, y=GEBV, color=gpd_measurement)) + geom_point(alpha=0.7,size=5)+
  scale_color_manual(values = c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[2],"brown4", "cadetblue",wes_palette("Rushmore1")[4]))
sp
dev.off()


# Same as above but made on averages, Not used in manuscript
png(file=paste("Comparison_of_GEBVs_from_different_gpdMeasurements_sortedAfterAvg_Averages",Sys.Date(),".png",sep="_"),width=50,height=35,units="cm",res=300)

keyDF <- data.frame(key=individual_order,weight=1:length(individual_order))
combined_ready1 <- merge(mean_merged,keyDF,by.x='Individual',by.y='key',all.x=T,all.y=F)
combined_ordered1 <- combined_ready1[order(combined_ready1$weight),c('Individual','Observed','GEBVs','gpd_measurement', 'weight')]

combined_ordered1$Individual <- factor(combined_ordered1$Individual,levels = individual_order)

sp<-ggplot(combined_ordered1, aes(x= Individual, y=GEBVs, color=gpd_measurement)) + geom_point(alpha=0.8,size=4.5)+
  scale_color_manual(values = c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[2],"brown4", "cadetblue",wes_palette("Rushmore1")[4]))+
theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp

dev.off()


# Combine the two above plots, Figure 4 in manuscript
#png(file="Comparison_of_GEBVs_from_different_gpdMeasurements_sortedAfterAvg_Both.png",width=50,height=35,units="cm",res=300)
sp<-ggplot(combined_ordered, aes(x= Individual, y=GEBV, color=gpd_measurement)) + geom_point(alpha=0.1,size=5)+
  scale_color_manual(values = c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[2],"brown4", "cadetblue",wes_palette("Rushmore1")[4]))+
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

sp1 = sp + geom_point(data=combined_ordered1,aes(Individual,GEBVs,color=gpd_measurement),size=7.5,alpha=0.9)
sp1

ggsave("Comparison_of_GEBVs_from_different_gpdMeasurements_sortedAfterAvg_Both.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)
#dev.off()


# Figure 2 but sorted after avg GEBVs of ResCor and parents marked 
parentpop1names=c("Aoost_01","Aoost_08","Aoost_09","Banna_02","Banna_03","Banna_07")
parentpop6names=c("Aearl_08","Ccyma_03","Llanc_06","Aaran_08")
parentpop7names=c("Aearl_05","Clfin_02","Ctain_05","Mrida_04")
parentpop8names=c("Clfin_03","Ctain_05","Volin_01","Aaran_04")
parentpop9names=c("Aoost_02","Ilona_09","Llanc_09","Sster_01")
parentpop10names=c("Ilona_05","Kdike_09","Llanc_09","Aalon_03")
parentpop11names=c("Ancor_10","Borek_06","Ctain_09","Rbani_02")
parentpop12names=c("Ancor_04","Aoost_10","Clfin_08","Kdike_08")
parentpop13names=c("Aoost_01","Aoost_08","Banna_02","Rbani_02","Sster_01","Sster_06")

list1=append(parentpop1names, parentpop6names)
list2=append(list1, parentpop7names)
list3=append(list2, parentpop8names)
list4=append(list3, parentpop9names)
list5=append(list4, parentpop10names)
list6=append(list5, parentpop11names)
list7=append(list6, parentpop12names)
list8=append(list7, parentpop13names)

All_parent_names=list8



# third color
# Colour palette used in manuscript

# Marking parents in different colours depending on the pop. they belong to
# y-values are GEBVs from gpd_rescor
#gpd_rescor_data=combined_ordered[which(combined_ordered$gpd_measurement=="gpd_ResCor"),]
df_ordered$group=as.character(rep("Non-parent",nrow(df_ordered)))

df_ordered$group[which(df_ordered$Individual=="Aoost_01")]="DLF1+SUA"
df_ordered$group[which(df_ordered$Individual=="Aoost_08")]="DLF1+SUA"
df_ordered$group[which(df_ordered$Individual=="Banna_02")]="DLF1+SUA"
df_ordered$group[which(df_ordered$Individual=="Ctain_05")]="LJ_L7+LJ_L9"
df_ordered$group[which(df_ordered$Individual=="Llanc_09")]="LJ_H1+LJ_H2"
df_ordered$group[which(df_ordered$Individual=="Sster_01")]="LJ_H1+SUA"
df_ordered$group[which(df_ordered$Individual=="Rbani_02")]="LJ_H3+SUA"

df_ordered$group[which(df_ordered$Individual=="Aoost_09")]="DLF1"
df_ordered$group[which(df_ordered$Individual=="Banna_03")]="DLF1"
df_ordered$group[which(df_ordered$Individual=="Banna_07")]="DLF1"

df_ordered$group[which(df_ordered$Individual=="Aearl_08")]="LJ_L6"
df_ordered$group[which(df_ordered$Individual=="Ccyma_03")]="LJ_L6"
df_ordered$group[which(df_ordered$Individual=="Llanc_06")]="LJ_L6"
df_ordered$group[which(df_ordered$Individual=="Aaran_08")]="LJ_L6"

df_ordered$group[which(df_ordered$Individual=="Aearl_05")]="LJ_L7"
df_ordered$group[which(df_ordered$Individual=="Clfin_02")]="LJ_L7"
df_ordered$group[which(df_ordered$Individual=="Mrida_04")]="LJ_L7"

df_ordered$group[which(df_ordered$Individual=="Clfin_03")]="LJ_L9"
df_ordered$group[which(df_ordered$Individual=="Volin_01")]="LJ_L9"
df_ordered$group[which(df_ordered$Individual=="Aaran_04")]="LJ_L9"

df_ordered$group[which(df_ordered$Individual=="Aoost_02")]="LJ_H1"
df_ordered$group[which(df_ordered$Individual=="Ilona_09")]="LJ_H1"

df_ordered$group[which(df_ordered$Individual=="Ilona_05")]="LJ_H2"
df_ordered$group[which(df_ordered$Individual=="Kdike_09")]="LJ_H2"
df_ordered$group[which(df_ordered$Individual=="Llanc_09")]="LJ_H2"
df_ordered$group[which(df_ordered$Individual=="Aalon_03")]="LJ_H2"

df_ordered$group[which(df_ordered$Individual=="Ancor_10")]="LJ_H3"
df_ordered$group[which(df_ordered$Individual=="Borek_06")]="LJ_H3"
df_ordered$group[which(df_ordered$Individual=="Ctain_09")]="LJ_H3"

df_ordered$group[which(df_ordered$Individual=="Ancor_04")]="LJ_H5"
df_ordered$group[which(df_ordered$Individual=="Aoost_10")]="LJ_H5"
df_ordered$group[which(df_ordered$Individual=="Clfin_08")]="LJ_H5"
df_ordered$group[which(df_ordered$Individual=="Kdike_08")]="LJ_H5"

df_ordered$group[which(df_ordered$Individual=="Sster_06")]="SUA"


dataset$group=as.character(rep("Non-parent",nrow(dataset)))

dataset$group[which(dataset$ID=="Aoost_01")]="DLF1+SUA"
dataset$group[which(dataset$ID=="Aoost_08")]="DLF1+SUA"
dataset$group[which(dataset$ID=="Banna_02")]="DLF1+SUA"
dataset$group[which(dataset$ID=="Ctain_05")]="LJ_L7+LJ_L9"
dataset$group[which(dataset$ID=="Llanc_09")]="LJ_H1+LJ_H2"
dataset$group[which(dataset$ID=="Sster_01")]="LJ_H1+SUA"
dataset$group[which(dataset$ID=="Rbani_02")]="LJ_H3+SUA"

dataset$group[which(dataset$ID=="Aoost_09")]="DLF1"
dataset$group[which(dataset$ID=="Banna_03")]="DLF1"
dataset$group[which(dataset$ID=="Banna_07")]="DLF1"

dataset$group[which(dataset$ID=="Aearl_08")]="LJ_L6"
dataset$group[which(dataset$ID=="Ccyma_03")]="LJ_L6"
dataset$group[which(dataset$ID=="Llanc_06")]="LJ_L6"
dataset$group[which(dataset$ID=="Aaran_08")]="LJ_L6"

dataset$group[which(dataset$ID=="Aearl_05")]="LJ_L7"
dataset$group[which(dataset$ID=="Clfin_02")]="LJ_L7"
dataset$group[which(dataset$ID=="Mrida_04")]="LJ_L7"

dataset$group[which(dataset$ID=="Clfin_03")]="LJ_L9"
dataset$group[which(dataset$ID=="Volin_01")]="LJ_L9"
dataset$group[which(dataset$ID=="Aaran_04")]="LJ_L9"

dataset$group[which(dataset$ID=="Aoost_02")]="LJ_H1"
dataset$group[which(dataset$ID=="Ilona_09")]="LJ_H1"

dataset$group[which(dataset$ID=="Ilona_05")]="LJ_H2"
dataset$group[which(dataset$ID=="Kdike_09")]="LJ_H2"
dataset$group[which(dataset$ID=="Llanc_09")]="LJ_H2"
dataset$group[which(dataset$ID=="Aalon_03")]="LJ_H2"

dataset$group[which(dataset$ID=="Ancor_10")]="LJ_H3"
dataset$group[which(dataset$ID=="Borek_06")]="LJ_H3"
dataset$group[which(dataset$ID=="Ctain_09")]="LJ_H3"

dataset$group[which(dataset$ID=="Ancor_04")]="LJ_H5"
dataset$group[which(dataset$ID=="Aoost_10")]="LJ_H5"
dataset$group[which(dataset$ID=="Clfin_08")]="LJ_H5"
dataset$group[which(dataset$ID=="Kdike_08")]="LJ_H5"

dataset$group[which(dataset$ID=="Sster_06")]="SUA"


cairo_pdf(paste("ParentSelection_GEBV_gpd",Sys.Date(),".pdf",sep="_"), family="Arial Unicode MS",width=60,height=20)

Datasetorder0=dataset
colnames(Datasetorder0)[1]="Individual"
Datasetorder1= merge(Datasetorder0,keyDF,by.x='Individual',by.y='key',all.x=T,all.y=F)
Datasetorder2 = Datasetorder1[order(Datasetorder1$weight),c('Individual','Observed','GEBV','group', 'weight')]
individual_order=unique(Datasetorder2$Individual)
Datasetorder2$Individual <- factor(Datasetorder2$Individual,levels = individual_order)


sp<-ggplot(Datasetorder2, aes(x= Individual, y=GEBV, color=group)) + geom_point(alpha=0.1,size=10)+
  scale_color_manual(values=c("#242c5e","black","#443659","black","#580b1d","#98d3ce","black","#71a188","#b01f35","#f6c29d","black","#f48988","grey","#426468")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) +  coord_cartesian(clip = "off")

df_ordered$Observed=as.numeric(as.character(df_ordered$Observed))
df_ordered$GEBVs=as.numeric(as.character(df_ordered$GEBVs))

df_layer_1 <- df_ordered[df_ordered$group=="Non-parent",]
df_layer_2 <- df_ordered[df_ordered$group=="DLF1",]
df_layer_3 <- df_ordered[df_ordered$group=="LJ_L6",]
df_layer_4 <- df_ordered[df_ordered$group=="LJ_L7",]
df_layer_5 <- df_ordered[df_ordered$group=="LJ_L9",]
df_layer_6 <- df_ordered[df_ordered$group=="LJ_H1",]
df_layer_7 <- df_ordered[df_ordered$group=="LJ_H2",]
df_layer_8 <- df_ordered[df_ordered$group=="LJ_H3",]
df_layer_9 <- df_ordered[df_ordered$group=="LJ_H5",]
df_layer_10 <- df_ordered[df_ordered$group=="SUA",]
df_layer_11 <- df_ordered[df_ordered$group=="DLF1+SUA",]
df_layer_12 <- df_ordered[df_ordered$group=="LJ_H1+SUA",]
df_layer_13 <- df_ordered[df_ordered$group=="LJ_H3+SUA",]
df_layer_14 <- df_ordered[df_ordered$group=="LJ_L7+LJ_L9",]


sp1 = sp + geom_point(data=df_layer_1,aes(Individual,GEBVs),size=15,alpha=1,color="grey") +
  geom_point(data=df_layer_2,aes(Individual,GEBVs),size=15,alpha=1,color="#242c5e") + 
  geom_point(data=df_layer_3,aes(Individual,GEBVs),size=15,alpha=1,color="#b01f35") +
  geom_point(data=df_layer_4,aes(Individual,GEBVs),size=15,alpha=1,color="#f6c29d") +
  geom_point(data=df_layer_5,aes(Individual,GEBVs),size=15,alpha=1,color="#f48988") +
  geom_point(data=df_layer_6,aes(Individual,GEBVs),size=15,alpha=1,color="#443659") +
  geom_point(data=df_layer_7,aes(Individual,GEBVs),size=15,alpha=1,color="#580b1d") +
  geom_point(data=df_layer_8,aes(Individual,GEBVs),size=15,alpha=1,color="#98d3ce") +
  geom_point(data=df_layer_9,aes(Individual,GEBVs),size=15,alpha=1,color="#71a188") +
  geom_point(data=df_layer_10,aes(Individual,GEBVs),size=15,alpha=1,color="#426468") +
  #geom_point(data=df_layer_11,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_11,aes(Individual,GEBVs), shape="\u25D6", colour="#242c5e", size=30,alpha=1) +
  geom_point(data=df_layer_11,aes(Individual,GEBVs), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_12,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_12,aes(Individual,GEBVs), shape="\u25D6", colour="#443659", size=30,alpha=1) +
  geom_point(data=df_layer_12,aes(Individual,GEBVs), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_13,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") 
  geom_point(data=df_layer_13,aes(Individual,GEBVs), shape="\u25D6", colour="#98d3ce", size=30,alpha=1) +
  geom_point(data=df_layer_13,aes(Individual,GEBVs), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_14,aes(Individual,GEBVs),size=7.5,alpha=1,color="black") 
  geom_point(data=df_layer_14,aes(Individual,GEBVs), shape="\u25D6", colour="#f6c29d", size=30,alpha=1) +
  geom_point(data=df_layer_14,aes(Individual,GEBVs), shape="\u25D7", colour="#f48988", size=30,alpha=1) 

sp1

dev.off()






# Colour palette used in manuscript
# GEBVs iSize

df_isize=as.data.frame(df_isize)
df_isize$Observed=as.numeric(as.character(df_isize$Observed))
df_isize$GEBVs=as.numeric(as.character(df_isize$GEBVs))

df_isize_ordered <- df_isize[order(df_isize$GEBVs),]
individual_order=df_isize_ordered$Individual

df_isize_ordered$group=as.character(rep("Non-parent",nrow(df_isize_ordered)))

df_isize_ordered$group[which(df_isize_ordered$Individual=="Aoost_01")]="DLF1+SUA"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Aoost_08")]="DLF1+SUA"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Banna_02")]="DLF1+SUA"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Ctain_05")]="LJ_L7+LJ_L9"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Llanc_09")]="LJ_H1+LJ_H2"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Sster_01")]="LJ_H1+SUA"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Rbani_02")]="LJ_H3+SUA"

df_isize_ordered$group[which(df_isize_ordered$Individual=="Aoost_09")]="DLF1"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Banna_03")]="DLF1"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Banna_07")]="DLF1"

df_isize_ordered$group[which(df_isize_ordered$Individual=="Aearl_08")]="LJ_L6"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Ccyma_03")]="LJ_L6"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Llanc_06")]="LJ_L6"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Aaran_08")]="LJ_L6"

df_isize_ordered$group[which(df_isize_ordered$Individual=="Aearl_05")]="LJ_L7"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Clfin_02")]="LJ_L7"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Mrida_04")]="LJ_L7"

df_isize_ordered$group[which(df_isize_ordered$Individual=="Clfin_03")]="LJ_L9"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Volin_01")]="LJ_L9"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Aaran_04")]="LJ_L9"

df_isize_ordered$group[which(df_isize_ordered$Individual=="Aoost_02")]="LJ_H1"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Ilona_09")]="LJ_H1"

df_isize_ordered$group[which(df_isize_ordered$Individual=="Ilona_05")]="LJ_H2"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Kdike_09")]="LJ_H2"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Llanc_09")]="LJ_H2"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Aalon_03")]="LJ_H2"

df_isize_ordered$group[which(df_isize_ordered$Individual=="Ancor_10")]="LJ_H3"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Borek_06")]="LJ_H3"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Ctain_09")]="LJ_H3"

df_isize_ordered$group[which(df_isize_ordered$Individual=="Ancor_04")]="LJ_H5"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Aoost_10")]="LJ_H5"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Clfin_08")]="LJ_H5"
df_isize_ordered$group[which(df_isize_ordered$Individual=="Kdike_08")]="LJ_H5"

df_isize_ordered$group[which(df_isize_ordered$Individual=="Sster_06")]="SUA"


dataset5$group=as.character(rep("Non-parent",nrow(dataset5)))

dataset5$group[which(dataset5$ID=="Aoost_01")]="DLF1+SUA"
dataset5$group[which(dataset5$ID=="Aoost_08")]="DLF1+SUA"
dataset5$group[which(dataset5$ID=="Banna_02")]="DLF1+SUA"
dataset5$group[which(dataset5$ID=="Ctain_05")]="LJ_L7+LJ_L9"
dataset5$group[which(dataset5$ID=="Llanc_09")]="LJ_H1+LJ_H2"
dataset5$group[which(dataset5$ID=="Sster_01")]="LJ_H1+SUA"
dataset5$group[which(dataset5$ID=="Rbani_02")]="LJ_H3+SUA"

dataset5$group[which(dataset5$ID=="Aoost_09")]="DLF1"
dataset5$group[which(dataset5$ID=="Banna_03")]="DLF1"
dataset5$group[which(dataset5$ID=="Banna_07")]="DLF1"

dataset5$group[which(dataset5$ID=="Aearl_08")]="LJ_L6"
dataset5$group[which(dataset5$ID=="Ccyma_03")]="LJ_L6"
dataset5$group[which(dataset5$ID=="Llanc_06")]="LJ_L6"
dataset5$group[which(dataset5$ID=="Aaran_08")]="LJ_L6"

dataset5$group[which(dataset5$ID=="Aearl_05")]="LJ_L7"
dataset5$group[which(dataset5$ID=="Clfin_02")]="LJ_L7"
dataset5$group[which(dataset5$ID=="Mrida_04")]="LJ_L7"

dataset5$group[which(dataset5$ID=="Clfin_03")]="LJ_L9"
dataset5$group[which(dataset5$ID=="Volin_01")]="LJ_L9"
dataset5$group[which(dataset5$ID=="Aaran_04")]="LJ_L9"

dataset5$group[which(dataset5$ID=="Aoost_02")]="LJ_H1"
dataset5$group[which(dataset5$ID=="Ilona_09")]="LJ_H1"

dataset5$group[which(dataset5$ID=="Ilona_05")]="LJ_H2"
dataset5$group[which(dataset5$ID=="Kdike_09")]="LJ_H2"
dataset5$group[which(dataset5$ID=="Llanc_09")]="LJ_H2"
dataset5$group[which(dataset5$ID=="Aalon_03")]="LJ_H2"

dataset5$group[which(dataset5$ID=="Ancor_10")]="LJ_H3"
dataset5$group[which(dataset5$ID=="Borek_06")]="LJ_H3"
dataset5$group[which(dataset5$ID=="Ctain_09")]="LJ_H3"

dataset5$group[which(dataset5$ID=="Ancor_04")]="LJ_H5"
dataset5$group[which(dataset5$ID=="Aoost_10")]="LJ_H5"
dataset5$group[which(dataset5$ID=="Clfin_08")]="LJ_H5"
dataset5$group[which(dataset5$ID=="Kdike_08")]="LJ_H5"

dataset5$group[which(dataset5$ID=="Sster_06")]="SUA"



cairo_pdf(paste("ParentSelection_GEBV_iSize",Sys.Date(),".pdf",sep="_"), family="Arial Unicode MS",width=60,height=20)

Datasetorder50=dataset5
colnames(Datasetorder50)[1]="Individual"
individual_order= df_isize_ordered$Individual
keyDF <- data.frame(key=individual_order,weight=1:length(individual_order))

Datasetorder51= merge(Datasetorder50,keyDF,by.x='Individual',by.y='key',all.x=T,all.y=F)
Datasetorder52 = Datasetorder51[order(Datasetorder51$weight),c('Individual','Observed','GEBV','group', 'weight')]

Datasetorder52$Individual <- factor(Datasetorder52$Individual,levels = individual_order)

sp<-ggplot(Datasetorder52, aes(x= Individual, y=GEBV, color=group)) + geom_point(alpha=0.1,size=10)+
  scale_color_manual(values=c("#242c5e","black","#443659","black","#580b1d","#98d3ce","black","#71a188","#b01f35","#f6c29d","black","#f48988","grey","#426468")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) +  coord_cartesian(clip = "off")



df_isize_ordered$Observed=as.numeric(as.character(df_isize_ordered$Observed))
df_isize_ordered$GEBVs=as.numeric(as.character(df_isize_ordered$GEBVs))

df_layer_1 <- df_isize_ordered[df_isize_ordered$group=="Non-parent",]
df_layer_2 <- df_isize_ordered[df_isize_ordered$group=="DLF1",]
df_layer_3 <- df_isize_ordered[df_isize_ordered$group=="LJ_L6",]
df_layer_4 <- df_isize_ordered[df_isize_ordered$group=="LJ_L7",]
df_layer_5 <- df_isize_ordered[df_isize_ordered$group=="LJ_L9",]
df_layer_6 <- df_isize_ordered[df_isize_ordered$group=="LJ_H1",]
df_layer_7 <- df_isize_ordered[df_isize_ordered$group=="LJ_H2",]
df_layer_8 <- df_isize_ordered[df_isize_ordered$group=="LJ_H3",]
df_layer_9 <- df_isize_ordered[df_isize_ordered$group=="LJ_H5",]
df_layer_10 <- df_isize_ordered[df_isize_ordered$group=="SUA",]
df_layer_11 <- df_isize_ordered[df_isize_ordered$group=="DLF1+SUA",]
df_layer_12 <- df_isize_ordered[df_isize_ordered$group=="LJ_H1+SUA",]
df_layer_13 <- df_isize_ordered[df_isize_ordered$group=="LJ_H3+SUA",]
df_layer_14 <- df_isize_ordered[df_isize_ordered$group=="LJ_L7+LJ_L9",]


sp1 = sp + geom_point(data=df_layer_1,aes(Individual,GEBVs),size=15,alpha=1,color="grey") +
  geom_point(data=df_layer_2,aes(Individual,GEBVs),size=15,alpha=1,color="#242c5e") + 
  geom_point(data=df_layer_3,aes(Individual,GEBVs),size=15,alpha=1,color="#b01f35") +
  geom_point(data=df_layer_4,aes(Individual,GEBVs),size=15,alpha=1,color="#f6c29d") +
  geom_point(data=df_layer_5,aes(Individual,GEBVs),size=15,alpha=1,color="#f48988") +
  geom_point(data=df_layer_6,aes(Individual,GEBVs),size=15,alpha=1,color="#443659") +
  geom_point(data=df_layer_7,aes(Individual,GEBVs),size=15,alpha=1,color="#580b1d") +
  geom_point(data=df_layer_8,aes(Individual,GEBVs),size=15,alpha=1,color="#98d3ce") +
  geom_point(data=df_layer_9,aes(Individual,GEBVs),size=15,alpha=1,color="#71a188") +
  geom_point(data=df_layer_10,aes(Individual,GEBVs),size=15,alpha=1,color="#426468") +
  #geom_point(data=df_layer_11,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_11,aes(Individual,GEBVs), shape="\u25D6", colour="#242c5e", size=30,alpha=1) +
  geom_point(data=df_layer_11,aes(Individual,GEBVs), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_12,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_12,aes(Individual,GEBVs), shape="\u25D6", colour="#443659", size=30,alpha=1) +
  geom_point(data=df_layer_12,aes(Individual,GEBVs), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_13,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") 
  geom_point(data=df_layer_13,aes(Individual,GEBVs), shape="\u25D6", colour="#98d3ce", size=30,alpha=1) +
  geom_point(data=df_layer_13,aes(Individual,GEBVs), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_14,aes(Individual,GEBVs),size=7.5,alpha=1,color="black") 
  geom_point(data=df_layer_14,aes(Individual,GEBVs), shape="\u25D6", colour="#f6c29d", size=30,alpha=1) +
  geom_point(data=df_layer_14,aes(Individual,GEBVs), shape="\u25D7", colour="#f48988", size=30,alpha=1) 

sp1

dev.off()












