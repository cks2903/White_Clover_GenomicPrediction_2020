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

# Load predictions_gpd_ResCor
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_gpd_rescor_6foldCV_20200828")

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
dataset$Difference=dataset$Observed_correctedForFixed-dataset$GEBVs
dataset$Num_Difference=abs(dataset$Difference)

# Calculate mean for each individual
df=aggregate(as.numeric(dataset$Num_Difference), list(dataset$Individual), mean)
df=as.data.frame(df)
meangpd=aggregate(as.numeric(dataset$Observed_correctedForFixed), list(dataset$Individual), mean)[,2]
meangebv=aggregate(as.numeric(dataset$GEBVs), list(dataset$Individual), mean)[,2]
meandifference=aggregate(as.numeric(dataset$Difference), list(dataset$Individual), mean)[,2]
df=cbind(df,meangpd,meangebv,meandifference)
colnames(df)=c("Individual","Num_Difference","Observed_correctedForFixed","GEBVs","Difference")
head(df)




# Load predictions_gpd_FixCor
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_gpd_Fixcor_6foldCV_20200828")
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
dataset2$Difference=dataset2$Observed_correctedForFixed-dataset2$GEBVs
dataset2$Num_Difference=abs(dataset2$Difference)

# Calculate mean for each individual
df_fix=aggregate(as.numeric(dataset2$Num_Difference), list(dataset2$Individual), mean)
df_fix=as.data.frame(df_fix)
meangpd=aggregate(as.numeric(dataset2$Observed_correctedForFixed), list(dataset2$Individual), mean)[,2]
meangebv=aggregate(as.numeric(dataset2$GEBVs), list(dataset2$Individual), mean)[,2]
meandifference=aggregate(as.numeric(dataset2$Difference), list(dataset2$Individual), mean)[,2]
df_fix=cbind(df_fix,meangpd,meangebv,meandifference)
colnames(df_fix)=c("Individual","Num_Difference","Observed_correctedForFixed","GEBVs","Difference")
head(df_fix)


# Load predictions_gpd_NoCor
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_gpd_Nocor_6foldCV_20200811")
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
dataset3$Difference=dataset3$Observed_correctedForFixed-dataset3$GEBVs
dataset3$Num_Difference=abs(dataset3$Difference)

# Calculate mean for each individual
df_No=aggregate(as.numeric(dataset3$Num_Difference), list(dataset3$Individual), mean)
df_No=as.data.frame(df_No)
meangpd=aggregate(as.numeric(dataset3$Observed_correctedForFixed), list(dataset3$Individual), mean)[,2]
meangebv=aggregate(as.numeric(dataset3$GEBVs), list(dataset3$Individual), mean)[,2]
meandifference=aggregate(as.numeric(dataset3$Difference), list(dataset3$Individual), mean)[,2]
df_No=cbind(df_No,meangpd,meangebv,meandifference)
colnames(df_No)=c("Individual","Num_Difference","Observed_correctedForFixed","GEBVs","Difference")
head(df_No)



# Check if they agree on best individuals
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter")


cor(df_fix$GEBVs,df$GEBVs) #0.75
cor(df_No$GEBVs,df$GEBVs) #0.99
df_merged <-merge(df_fix,df_No, by=("Individual"), all.x=TRUE, all.y=TRUE); df 
df_merged1 <-merge(df_merged,df, by=("Individual"), all.x=TRUE, all.y=TRUE); df 

df_merged1=df_merged1[,c(1,3,4,7,8,11,12)]
colnames(df_merged1)=c("Individual","CorrectedPheno_gpd_FixCor","GEBV_gpd_FixCor","CorrectedPheno_gpd_NoCor","GEBV_gpd_NoCor","CorrectedPheno_gpd_ResCor","GEBV_gpd_ResCor")


# Not used in manuscript
png(file="GEBVS_gpd_NoCor_and_gpd_ResCor.png",width=35,height=20,units="cm",res=300, pointsize=6)

ggplot(data=df_merged1, aes(
  x = GEBV_gpd_NoCor,
  y = GEBV_gpd_ResCor))+
  #geom_point(alpha = 0.7) +
  geom_point(alpha=0.7) +
  #scale_colour_manual(values=setNames(color.codes, parents)) +
  theme_classic()  +
  labs(title="Predicted GEBVs",x="GEBV_gpd_NoCor", y = "GEBV_gpd_ResCor") +
  theme(axis.text.x = element_text(angle = 90)) 

dev.off()

# Not used in manuscript
png(file="GEBVs_gpd_ResCor.png",width=35,height=20,units="cm",res=300, pointsize=6)

ggplot(data=dataset, aes(
  x = Individual,
  y = GEBVs))+
  #geom_point(alpha = 0.7) +
  geom_point(alpha=0.7) +
  #scale_colour_manual(values=setNames(color.codes, parents)) +
  theme_classic()  +
  labs(title="Predicted GEBVs gpd_ResCor",x="Individual", y = "GEBVs") +
  theme(axis.text.x = element_text(angle = 90)) 

dev.off()


# Not used in manuscript
png(file="GEBVs_gpd_NoCor.png",width=35,height=20,units="cm",res=300, pointsize=6)

ggplot(data=dataset3, aes(
  x = Individual,
  y = GEBVs))+
  #geom_point(alpha = 0.7) +
  geom_point(alpha=0.7) +
  #scale_colour_manual(values=setNames(color.codes, parents)) +
  theme_classic()  +
  labs(title="Predicted GEBVs gpd_NoCor",x="Individual", y = "GEBVs") +
  theme(axis.text.x = element_text(angle = 90)) 

dev.off()



# Not used in manuscript
png(file="GEBVs_gpd_FixCor.png",width=35,height=20,units="cm",res=300, pointsize=6)

ggplot(data=dataset2, aes(
  x = Individual,
  y = GEBVs))+
  #geom_point(alpha = 0.7) +
  geom_point(alpha=0.7) +
  #scale_colour_manual(values=setNames(color.codes, parents)) +
  theme_classic()  +
  labs(title="Predicted GEBVs gpd_FixCor",x="Individual", y = "GEBVs") +
  theme(axis.text.x = element_text(angle = 90)) 

dev.off()





# combine in one plot
combined=rbind(dataset,dataset2,dataset3)
combined$gpd_measurement=rep("NA",nrow(combined))
combined$gpd_measurement[1:nrow(dataset)]="gpd_ResCor"
combined$gpd_measurement[(nrow(dataset)+1):(nrow(dataset)+nrow(dataset2))]="gpd_FixCor"
combined$gpd_measurement[(nrow(dataset)+nrow(dataset2)+1):(nrow(dataset)+nrow(dataset2)+nrow(dataset3))]="gpd_NoCor"

# Not used in manuscript
png(file="Comparison_of_GEBVs_from_different_gpdMeasurements.png",width=150,height=60,units="cm",res=300, pointsize=6)

e <- ggplot(combined, aes(x = Individual, y = GEBVs))

e2 <- e + geom_boxplot(
  aes(fill = gpd_measurement),
  position = position_dodge(0.9) 
) +
  scale_fill_manual(values = c("#253494", "#41b6c4","#c7e9b4"))+
  theme(axis.text.x = element_text(angle = 90)) 
e2

dev.off() # So GEBVs of gpd_rescor is not more precise



# Check if differences between GEBVs get better
# Not used in manuscript
png(file="Comparison_of_GEBVs_from_different_gpdMeasurements_2.png",width=80,height=40,units="cm",res=300, pointsize=6)

sp<-ggplot(combined, aes(x=Individual, y=GEBVs, color=gpd_measurement)) + geom_point()+theme(axis.text.x = element_text(angle = 90))
sp
dev.off()

mean_merged=rbind(df,df_fix,df_No)
mean_merged$gpd_measurement=rep("NA",nrow(mean_merged))
mean_merged$gpd_measurement[1:145]="gpd_ResCor"
mean_merged$gpd_measurement[146:290]="gpd_FixCor"
mean_merged$gpd_measurement[291:435]="gpd_NoCor"

# Not used in manuscript
png(file="Comparison_of_GEBVs_from_different_gpdMeasurements_3.png",width=80,height=40,units="cm",res=300, pointsize=10)
sp<-ggplot(mean_merged, aes(x=Individual, y=GEBVs, color=gpd_measurement)) + geom_point()
sp
dev.off()

# Not used in manuscript
png(file="Comparison_of_GEBVs_from_different_gpdMeasurements_4.png",width=80,height=40,units="cm",res=300, pointsize=6)
sp2<-ggplot(mean_merged, aes(x=Individual, y=GEBVs, color=gpd_measurement)) + geom_point() + geom_line(aes(group=gpd_measurement)) +theme(axis.text.x = element_text(angle = 90)) 
sp2
dev.off()



# Sorted after avg GEBVs of ResCor, Not used in manuscript
png(file="Comparison_of_GEBVs_from_different_gpdMeasurements_sortedAfterAvg.png",width=70,height=60,units="cm",res=300)

df_ordered <- df[order(df$GEBVs),]
individual_order=df_ordered$Individual

keyDF <- data.frame(key=individual_order,weight=1:length(individual_order))
combined_ready <- merge(combined,keyDF,by.x='Individual',by.y='key',all.x=T,all.y=F)
combined_ordered <- combined_ready[order(combined_ready$weight),c('Individual','Observed_correctedForFixed','GEBVs','Difference','Num_Difference','gpd_measurement', 'weight')]

combined_ordered$Individual <- factor(combined_ordered$Individual,levels = individual_order)

sp<-ggplot(combined_ordered, aes(x= Individual, y=GEBVs, color=gpd_measurement)) + geom_point(alpha=0.7,size=5)+
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp

dev.off()


# Same as above but made on averages, Not used in manuscript
png(file="Comparison_of_GEBVs_from_different_gpdMeasurements_sortedAfterAvg_Averages.png",width=50,height=35,units="cm",res=300)

keyDF <- data.frame(key=individual_order,weight=1:length(individual_order))
combined_ready1 <- merge(mean_merged,keyDF,by.x='Individual',by.y='key',all.x=T,all.y=F)
combined_ordered1 <- combined_ready1[order(combined_ready1$weight),c('Individual','Observed_correctedForFixed','GEBVs','Difference','Num_Difference','gpd_measurement', 'weight')]

combined_ordered1$Individual <- factor(combined_ordered1$Individual,levels = individual_order)

sp<-ggplot(combined_ordered1, aes(x= Individual, y=GEBVs, color=gpd_measurement,shape=gpd_measurement)) + geom_point(alpha=0.8,size=4.5)+
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + scale_shape_manual(values=c(15,16,18)) +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp

dev.off()


# Combine the two above plots, Figure 4 in manuscript
#png(file="Comparison_of_GEBVs_from_different_gpdMeasurements_sortedAfterAvg_Both.png",width=50,height=35,units="cm",res=300)
sp<-ggplot(combined_ordered, aes(x= Individual, y=GEBVs, color=gpd_measurement, shape = gpd_measurement)) + geom_point(alpha=0.1,size=5)+
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + 
  scale_shape_manual(values=c(16,16,16)) +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

sp1 = sp + geom_point(data=combined_ordered1,aes(Individual,GEBVs,color=gpd_measurement,shape=gpd_measurement),size=7.5,alpha=0.9)
sp1

ggsave("Comparison_of_GEBVs_from_different_gpdMeasurements_sortedAfterAvg_Both.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)
#dev.off()


# Above plot but sorted by gpd_rescor phenotype
d6data=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/JustTheWhole_d6Table.txt",sep="\t",head=T)
gpd_ResCorMeans=aggregate(as.numeric(d6data$gpd_dryweight_cor), list(d6data$Clovershort), mean)
gpd_ResCorMeans_ordered <- gpd_ResCorMeans[order(gpd_ResCorMeans[,2]),]

individual_order=gpd_ResCorMeans_ordered[,1]

keyDF <- data.frame(key=individual_order,weight=1:length(individual_order))
combined_ready2 <- merge(combined,keyDF,by.x='Individual',by.y='key',all.x=T,all.y=F)
combined_ordered2 <- combined_ready2[order(combined_ready2$weight),c('Individual','Observed_correctedForFixed','GEBVs','Difference','Num_Difference','gpd_measurement', 'weight')]

combined_ordered2$Individual <- factor(combined_ordered2$Individual,levels = individual_order)

# Same as above but ordered according to gpd_rescor phenotypes. Not used in manuscript.
#png(file="Comparison_of_GEBVs_from_different_gpdMeasurements_sortedAftergpdResCorAvg_Both.png",width=50,height=35,units="cm",res=300)
sp<-ggplot(combined_ordered2, aes(x= Individual, y=GEBVs, color=gpd_measurement, shape = gpd_measurement)) + geom_point(alpha=0.1,size=5)+
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + 
  scale_shape_manual(values=c(16,16,16)) +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

sp1 = sp + geom_point(data=combined_ordered1,aes(Individual,GEBVs,color=gpd_measurement,shape=gpd_measurement),size=7.5,alpha=0.9)

sp1
ggsave("Comparison_of_GEBVs_from_different_gpdMeasurements_sortedAftergpdResCorAvg_Both.pdf", width =60, height = 30, units = "cm",useDingbats=FALSE)
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

# Not used in manuscript
png(file="AllParentsMarked.png",width=70,height=60,units="cm",res=300)

ggplot(combined_ordered, aes(x= Individual, y=GEBVs, color=gpd_measurement)) + geom_point(alpha=0.4,size=5)+
  gghighlight(Individual %in% All_parent_names) +ggtitle("All parents of F1 populations") +
scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

dev.off()


# Marking parents in different colours depending on the pop. they belong to
# y-values are GEBVs from gpd_rescor
#gpd_rescor_data=combined_ordered[which(combined_ordered$gpd_measurement=="gpd_ResCor"),]
combined_ordered$group=as.character(rep("Non-parent",nrow(combined_ordered)))

combined_ordered$group[which(combined_ordered$Individual=="Aoost_01")]="DLF1+SUA"
combined_ordered$group[which(combined_ordered$Individual=="Aoost_08")]="DLF1+SUA"
combined_ordered$group[which(combined_ordered$Individual=="Banna_02")]="DLF1+SUA"
combined_ordered$group[which(combined_ordered$Individual=="Ctain_05")]="LJ_L7+LJ_L9"
combined_ordered$group[which(combined_ordered$Individual=="Llanc_09")]="LJ_H1+LJ_H2"
combined_ordered$group[which(combined_ordered$Individual=="Sster_01")]="LJ_H1+SUA"
combined_ordered$group[which(combined_ordered$Individual=="Rbani_02")]="LJ_H3+SUA"

combined_ordered$group[which(combined_ordered$Individual=="Aoost_09")]="DLF1"
combined_ordered$group[which(combined_ordered$Individual=="Banna_03")]="DLF1"
combined_ordered$group[which(combined_ordered$Individual=="Banna_07")]="DLF1"

combined_ordered$group[which(combined_ordered$Individual=="Aearl_08")]="LJ_L6"
combined_ordered$group[which(combined_ordered$Individual=="Ccyma_03")]="LJ_L6"
combined_ordered$group[which(combined_ordered$Individual=="Llanc_06")]="LJ_L6"
combined_ordered$group[which(combined_ordered$Individual=="Aaran_08")]="LJ_L6"

combined_ordered$group[which(combined_ordered$Individual=="Aearl_05")]="LJ_L7"
combined_ordered$group[which(combined_ordered$Individual=="Clfin_02")]="LJ_L7"
combined_ordered$group[which(combined_ordered$Individual=="Mrida_04")]="LJ_L7"

combined_ordered$group[which(combined_ordered$Individual=="Clfin_03")]="LJ_L9"
combined_ordered$group[which(combined_ordered$Individual=="Volin_01")]="LJ_L9"
combined_ordered$group[which(combined_ordered$Individual=="Aaran_04")]="LJ_L9"

combined_ordered$group[which(combined_ordered$Individual=="Aoost_02")]="LJ_H1"
combined_ordered$group[which(combined_ordered$Individual=="Ilona_09")]="LJ_H1"

combined_ordered$group[which(combined_ordered$Individual=="Ilona_05")]="LJ_H2"
combined_ordered$group[which(combined_ordered$Individual=="Kdike_09")]="LJ_H2"
combined_ordered$group[which(combined_ordered$Individual=="Llanc_09")]="LJ_H2"
combined_ordered$group[which(combined_ordered$Individual=="Aalon_03")]="LJ_H2"

combined_ordered$group[which(combined_ordered$Individual=="Ancor_10")]="LJ_H3"
combined_ordered$group[which(combined_ordered$Individual=="Borek_06")]="LJ_H3"
combined_ordered$group[which(combined_ordered$Individual=="Ctain_09")]="LJ_H3"

combined_ordered$group[which(combined_ordered$Individual=="Ancor_04")]="LJ_H5"
combined_ordered$group[which(combined_ordered$Individual=="Aoost_10")]="LJ_H5"
combined_ordered$group[which(combined_ordered$Individual=="Clfin_08")]="LJ_H5"
combined_ordered$group[which(combined_ordered$Individual=="Kdike_08")]="LJ_H5"

combined_ordered$group[which(combined_ordered$Individual=="Sster_06")]="SUA"

gpd_resCor_data=combined_ordered[which(combined_ordered$gpd_measurement=="gpd_ResCor"),]

head(gpd_resCor_data)

means_to_plot_gpdRescor=aggregate(as.numeric(as.character(gpd_resCor_data$GEBVs)), list(gpd_resCor_data$Individual), mean)
colnames(means_to_plot_gpdRescor)=c("Individual","GEBVs")
means_to_plot_gpdRescor$group=as.character(rep("Non-parent",nrow(means_to_plot_gpdRescor)))

for (i in seq(1:nrow(means_to_plot_gpdRescor))){
  line=means_to_plot_gpdRescor$Individual[i]
  group=gpd_resCor_data$group[which(gpd_resCor_data$Individual==line)][1]
  means_to_plot_gpdRescor$group[i]=group
  }
                                  
# Colour palette not used in manuscript
cairo_pdf("ParentSelection_GEBV_gpdResCor.pdf", family="Arial Unicode MS",width=27,height=20)

sp<-ggplot(gpd_resCor_data, aes(x= Individual, y=GEBVs, color=group)) + geom_point(alpha=0.1,size=10)+
  scale_color_manual(values=c("#7fcdbb","black","#41b6c4","black","#1d91c0","#225ea8","41b6c4","#253494","#ffffd9","#edf8b1","41b6c4","#c7e9b4","grey","#081d58")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) +  coord_cartesian(clip = "off")


df_layer_1 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="Non-parent",]
df_layer_2 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="DLF1",]
df_layer_3 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L6",]
df_layer_4 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L7",]
df_layer_5 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L9",]
df_layer_6 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H1",]
df_layer_7 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H2",]
df_layer_8 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H3",]
df_layer_9 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H5",]
df_layer_10 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="SUA",]
df_layer_11 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="DLF1+SUA",]
df_layer_12 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H1+SUA",]
df_layer_13 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H3+SUA",]
df_layer_14 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L7+LJ_L9",]


sp1 = sp + geom_point(data=df_layer_1,aes(Individual,GEBVs),size=15,alpha=1,color="grey") +
  geom_point(data=df_layer_2,aes(Individual,GEBVs),size=15,alpha=1,color="#7fcdbb") + 
  geom_point(data=df_layer_3,aes(Individual,GEBVs),size=15,alpha=1,color="#ffffd9") +
  geom_point(data=df_layer_4,aes(Individual,GEBVs),size=15,alpha=1,color="#edf8b1") +
  geom_point(data=df_layer_5,aes(Individual,GEBVs),size=15,alpha=1,color="#c7e9b4") +
  geom_point(data=df_layer_6,aes(Individual,GEBVs),size=15,alpha=1,color="#41b6c4") +
  geom_point(data=df_layer_7,aes(Individual,GEBVs),size=15,alpha=1,color="#1d91c0") +
  geom_point(data=df_layer_8,aes(Individual,GEBVs),size=15,alpha=1,color="#225ea8") +
  geom_point(data=df_layer_9,aes(Individual,GEBVs),size=15,alpha=1,color="#253494") +
  geom_point(data=df_layer_10,aes(Individual,GEBVs),size=15,alpha=1,color="#081d58") +
  #geom_point(data=df_layer_11,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_11,aes(Individual,GEBVs), shape="\u25D6", colour="#7fcdbb", size=30,alpha=1) +
  geom_point(data=df_layer_11,aes(Individual,GEBVs), shape="\u25D7", colour="#081d58", size=30,alpha=1) +
  #geom_point(data=df_layer_12,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_12,aes(Individual,GEBVs), shape="\u25D6", colour="#41b6c4", size=30,alpha=1) +
  geom_point(data=df_layer_12,aes(Individual,GEBVs), shape="\u25D7", colour="#081d58", size=30,alpha=1) +
  #geom_point(data=df_layer_13,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") 
  geom_point(data=df_layer_13,aes(Individual,GEBVs), shape="\u25D6", colour="#225ea8", size=30,alpha=1) +
  geom_point(data=df_layer_13,aes(Individual,GEBVs), shape="\u25D7", colour="#081d58", size=30,alpha=1) +
  #geom_point(data=df_layer_14,aes(Individual,GEBVs),size=7.5,alpha=1,color="black") 
  geom_point(data=df_layer_14,aes(Individual,GEBVs), shape="\u25D6", colour="#edf8b1", size=30,alpha=1) +
  geom_point(data=df_layer_14,aes(Individual,GEBVs), shape="\u25D7", colour="#c7e9b4", size=30,alpha=1) 
  
  
sp1

dev.off()


# some other colours
# Colour palette not used in manuscript
cairo_pdf("ParentSelection_GEBV_gpdResCor_v2.pdf", family="Arial Unicode MS",width=27,height=20)

sp<-ggplot(gpd_resCor_data, aes(x= Individual, y=GEBVs, color=group)) + geom_point(alpha=0.1,size=10)+
  scale_color_manual(values=c("#580b1d","black","#f6c29d","black","#bcbd99","#98d3ce","black","#71a188","#b01f35","#ee3055","black","#f48988","grey","#426468")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) +  coord_cartesian(clip = "off")


df_layer_1 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="Non-parent",]
df_layer_2 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="DLF1",]
df_layer_3 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L6",]
df_layer_4 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L7",]
df_layer_5 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L9",]
df_layer_6 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H1",]
df_layer_7 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H2",]
df_layer_8 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H3",]
df_layer_9 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H5",]
df_layer_10 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="SUA",]
df_layer_11 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="DLF1+SUA",]
df_layer_12 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H1+SUA",]
df_layer_13 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H3+SUA",]
df_layer_14 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L7+LJ_L9",]


sp1 = sp + geom_point(data=df_layer_1,aes(Individual,GEBVs),size=15,alpha=1,color="grey") +
  geom_point(data=df_layer_2,aes(Individual,GEBVs),size=15,alpha=1,color="#580b1d") + 
  geom_point(data=df_layer_3,aes(Individual,GEBVs),size=15,alpha=1,color="#b01f35") +
  geom_point(data=df_layer_4,aes(Individual,GEBVs),size=15,alpha=1,color="#ee3055") +
  geom_point(data=df_layer_5,aes(Individual,GEBVs),size=15,alpha=1,color="#f48988") +
  geom_point(data=df_layer_6,aes(Individual,GEBVs),size=15,alpha=1,color="#f6c29d") +
  geom_point(data=df_layer_7,aes(Individual,GEBVs),size=15,alpha=1,color="#bcbd99") +
  geom_point(data=df_layer_8,aes(Individual,GEBVs),size=15,alpha=1,color="#98d3ce") +
  geom_point(data=df_layer_9,aes(Individual,GEBVs),size=15,alpha=1,color="#71a188") +
  geom_point(data=df_layer_10,aes(Individual,GEBVs),size=15,alpha=1,color="#426468") +
  #geom_point(data=df_layer_11,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_11,aes(Individual,GEBVs), shape="\u25D6", colour="#580b1d", size=30,alpha=1) +
  geom_point(data=df_layer_11,aes(Individual,GEBVs), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_12,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_12,aes(Individual,GEBVs), shape="\u25D6", colour="#f6c29d", size=30,alpha=1) +
  geom_point(data=df_layer_12,aes(Individual,GEBVs), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_13,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") 
  geom_point(data=df_layer_13,aes(Individual,GEBVs), shape="\u25D6", colour="#98d3ce", size=30,alpha=1) +
  geom_point(data=df_layer_13,aes(Individual,GEBVs), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_14,aes(Individual,GEBVs),size=7.5,alpha=1,color="black") 
  geom_point(data=df_layer_14,aes(Individual,GEBVs), shape="\u25D6", colour="#ee3055", size=30,alpha=1) +
  geom_point(data=df_layer_14,aes(Individual,GEBVs), shape="\u25D7", colour="#f48988", size=30,alpha=1) 

sp1

dev.off()

# third color
# Colour palette used in manuscript
cairo_pdf("ParentSelection_GEBV_gpdResCor_v3.pdf", family="Arial Unicode MS",width=60,height=20)


sp<-ggplot(gpd_resCor_data, aes(x= Individual, y=GEBVs, color=group)) + geom_point(alpha=0.1,size=10)+
  scale_color_manual(values=c("#242c5e","black","#443659","black","#580b1d","#98d3ce","black","#71a188","#b01f35","#f6c29d","black","#f48988","grey","#426468")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) +  coord_cartesian(clip = "off")


df_layer_1 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="Non-parent",]
df_layer_2 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="DLF1",]
df_layer_3 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L6",]
df_layer_4 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L7",]
df_layer_5 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L9",]
df_layer_6 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H1",]
df_layer_7 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H2",]
df_layer_8 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H3",]
df_layer_9 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H5",]
df_layer_10 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="SUA",]
df_layer_11 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="DLF1+SUA",]
df_layer_12 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H1+SUA",]
df_layer_13 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H3+SUA",]
df_layer_14 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L7+LJ_L9",]


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

# fourth color
# Colour palette not used in manuscript
cairo_pdf("ParentSelection_GEBV_gpdResCor_v4.pdf", family="Arial Unicode MS",width=27,height=20)


sp<-ggplot(gpd_resCor_data, aes(x= Individual, y=GEBVs, color=group)) + geom_point(alpha=0.1,size=10)+
  scale_color_manual(values=c("#8FB794","black","#285238","black","#1C2541","#E2BAB6","black","#F4C095","#645577","#008F75","black","#826251","grey","#185A77")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) +  coord_cartesian(clip = "off")


df_layer_1 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="Non-parent",]
df_layer_2 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="DLF1",]
df_layer_3 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L6",]
df_layer_4 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L7",]
df_layer_5 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L9",]
df_layer_6 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H1",]
df_layer_7 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H2",]
df_layer_8 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H3",]
df_layer_9 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H5",]
df_layer_10 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="SUA",]
df_layer_11 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="DLF1+SUA",]
df_layer_12 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H1+SUA",]
df_layer_13 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_H3+SUA",]
df_layer_14 <- means_to_plot_gpdRescor[means_to_plot_gpdRescor$group=="LJ_L7+LJ_L9",]


sp1 = sp + geom_point(data=df_layer_1,aes(Individual,GEBVs),size=15,alpha=1,color="grey") +
  geom_point(data=df_layer_2,aes(Individual,GEBVs),size=15,alpha=1,color="#8FB794") + 
  geom_point(data=df_layer_3,aes(Individual,GEBVs),size=15,alpha=1,color="#645577") +
  geom_point(data=df_layer_4,aes(Individual,GEBVs),size=15,alpha=1,color="#008F75") +
  geom_point(data=df_layer_5,aes(Individual,GEBVs),size=15,alpha=1,color="#826251") +
  geom_point(data=df_layer_6,aes(Individual,GEBVs),size=15,alpha=1,color="#285238") +
  geom_point(data=df_layer_7,aes(Individual,GEBVs),size=15,alpha=1,color="#1C2541") +
  geom_point(data=df_layer_8,aes(Individual,GEBVs),size=15,alpha=1,color="#E2BAB6") +
  geom_point(data=df_layer_9,aes(Individual,GEBVs),size=15,alpha=1,color="#F4C095") +
  geom_point(data=df_layer_10,aes(Individual,GEBVs),size=15,alpha=1,color="#185A77") +
  #geom_point(data=df_layer_11,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_11,aes(Individual,GEBVs), shape="\u25D6", colour="#8FB794", size=30,alpha=1) +
  geom_point(data=df_layer_11,aes(Individual,GEBVs), shape="\u25D7", colour="#185A77", size=30,alpha=1) +
  #geom_point(data=df_layer_12,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_12,aes(Individual,GEBVs), shape="\u25D6", colour="#285238", size=30,alpha=1) +
  geom_point(data=df_layer_12,aes(Individual,GEBVs), shape="\u25D7", colour="#185A77", size=30,alpha=1) +
  #geom_point(data=df_layer_13,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") 
  geom_point(data=df_layer_13,aes(Individual,GEBVs), shape="\u25D6", colour="#E2BAB6", size=30,alpha=1) +
  geom_point(data=df_layer_13,aes(Individual,GEBVs), shape="\u25D7", colour="#185A77", size=30,alpha=1) +
  #geom_point(data=df_layer_14,aes(Individual,GEBVs),size=7.5,alpha=1,color="black") 
  geom_point(data=df_layer_14,aes(Individual,GEBVs), shape="\u25D6", colour="#008F75", size=30,alpha=1) +
  geom_point(data=df_layer_14,aes(Individual,GEBVs), shape="\u25D7", colour="#826251", size=30,alpha=1) 

sp1

dev.off()


# Figure 4 in manuscript
# Marking parents in different colours depending on the pop. they belong to
# y-values are GEBVs from gpd_Nocor
df_No_ordered <- df_No[order(df_No$GEBVs),]
individual_order_no=df_No_ordered$Individual

keyDF <- data.frame(key=individual_order_no,weight=1:length(individual_order_no))
combined_ready_gpdNo <- merge(combined,keyDF,by.x='Individual',by.y='key',all.x=T,all.y=F)
combined_ordered_No <- combined_ready_gpdNo[order(combined_ready_gpdNo$weight),c('Individual','Observed_correctedForFixed','GEBVs','Difference','Num_Difference','gpd_measurement', 'weight')]

combined_ordered_No$Individual <- factor(combined_ordered_No$Individual,levels = individual_order_no)


combined_ordered_No$group=as.character(rep("Non-parent",nrow(combined_ordered_No)))

combined_ordered_No$group[which(combined_ordered_No$Individual=="Aoost_01")]="DLF1+SUA"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Aoost_08")]="DLF1+SUA"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Banna_02")]="DLF1+SUA"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Ctain_05")]="LJ_L7+LJ_L9"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Llanc_09")]="LJ_H1+LJ_H2"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Sster_01")]="LJ_H1+SUA"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Rbani_02")]="LJ_H3+SUA"

combined_ordered_No$group[which(combined_ordered_No$Individual=="Aoost_09")]="DLF1"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Banna_03")]="DLF1"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Banna_07")]="DLF1"

combined_ordered_No$group[which(combined_ordered_No$Individual=="Aearl_08")]="LJ_L6"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Ccyma_03")]="LJ_L6"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Llanc_06")]="LJ_L6"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Aaran_08")]="LJ_L6"

combined_ordered_No$group[which(combined_ordered_No$Individual=="Aearl_05")]="LJ_L7"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Clfin_02")]="LJ_L7"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Mrida_04")]="LJ_L7"

combined_ordered_No$group[which(combined_ordered_No$Individual=="Clfin_03")]="LJ_L9"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Volin_01")]="LJ_L9"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Aaran_04")]="LJ_L9"

combined_ordered_No$group[which(combined_ordered_No$Individual=="Aoost_02")]="LJ_H1"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Ilona_09")]="LJ_H1"

combined_ordered_No$group[which(combined_ordered_No$Individual=="Ilona_05")]="LJ_H2"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Kdike_09")]="LJ_H2"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Llanc_09")]="LJ_H2"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Aalon_03")]="LJ_H2"

combined_ordered_No$group[which(combined_ordered_No$Individual=="Ancor_10")]="LJ_H3"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Borek_06")]="LJ_H3"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Ctain_09")]="LJ_H3"

combined_ordered_No$group[which(combined_ordered_No$Individual=="Ancor_04")]="LJ_H5"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Aoost_10")]="LJ_H5"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Clfin_08")]="LJ_H5"
combined_ordered_No$group[which(combined_ordered_No$Individual=="Kdike_08")]="LJ_H5"

combined_ordered_No$group[which(combined_ordered_No$Individual=="Sster_06")]="SUA"


gpd_NoCor_data=combined_ordered_No[which(combined_ordered_No$gpd_measurement=="gpd_NoCor"),]
head(gpd_NoCor_data)

means_to_plot_gpdNocor=aggregate(as.numeric(as.character(gpd_NoCor_data$GEBVs)), list(gpd_NoCor_data$Individual), mean)
colnames(means_to_plot_gpdNocor)=c("Individual","GEBVs")
means_to_plot_gpdNocor$group=as.character(rep("Non-parent",nrow(means_to_plot_gpdNocor)))

for (i in seq(1:nrow(means_to_plot_gpdNocor))){
  line=means_to_plot_gpdNocor$Individual[i]
  group=gpd_NoCor_data$group[which(gpd_NoCor_data$Individual==line)][1]
  means_to_plot_gpdNocor$group[i]=group
}


cairo_pdf("ParentSelection_GEBV_gpdNoCor_v3.pdf", family="Arial Unicode MS",width=60,height=20)

sp<-ggplot(gpd_NoCor_data, aes(x= Individual, y=GEBVs, color=group)) + geom_point(alpha=0.1,size=10)+
  scale_color_manual(values=c("#242c5e","black","#443659","black","#580b1d","#98d3ce","black","#71a188","#b01f35","#f6c29d","black","#f48988","grey","#426468")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) +  coord_cartesian(clip = "off")


df_layer_1 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="Non-parent",]
df_layer_2 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="DLF1",]
df_layer_3 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="LJ_L6",]
df_layer_4 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="LJ_L7",]
df_layer_5 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="LJ_L9",]
df_layer_6 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="LJ_H1",]
df_layer_7 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="LJ_H2",]
df_layer_8 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="LJ_H3",]
df_layer_9 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="LJ_H5",]
df_layer_10 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="SUA",]
df_layer_11 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="DLF1+SUA",]
df_layer_12 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="LJ_H1+SUA",]
df_layer_13 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="LJ_H3+SUA",]
df_layer_14 <- means_to_plot_gpdNocor[means_to_plot_gpdNocor$group=="LJ_L7+LJ_L9",]


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



# Figure 4
# y-values are phenotypes (gpd_NoCor)
phenotypeperf_gpdNoCor=cbind(as.character(d6data$Clovershort),as.numeric(as.character(d6data$growth_per_day)))
colnames(phenotypeperf_gpdNoCor)=c("Individual","gpd_NoCor")
phenotypeperf_gpdNoCor=as.data.frame(phenotypeperf_gpdNoCor)


gpd_NoCorMeans=aggregate(as.numeric(as.character((phenotypeperf_gpdNoCor$gpd_NoCor))), list(phenotypeperf_gpdNoCor$Individual), mean)
gpd_NoCorMeans_ordered <- gpd_NoCorMeans[order(gpd_NoCorMeans[,2]),]
individual_order_gpdNoCorPheno=gpd_NoCorMeans_ordered[,1]

keyDF <- data.frame(key=individual_order_gpdNoCorPheno,weight=1:length(individual_order_gpdNoCorPheno))
combined_ready3 <- merge(phenotypeperf_gpdNoCor,keyDF,by.x='Individual',by.y='key',all.x=T,all.y=F)
combined_ordered3 <- combined_ready3[order(combined_ready3$weight),c('Individual','gpd_NoCor','weight')]

combined_ordered3$Individual <- factor(combined_ordered3$Individual,levels = individual_order_gpdNoCorPheno)

combined_ordered3$group=as.character(rep("Non-parent",nrow(combined_ordered3)))
combined_ordered3$group[which(combined_ordered3$Individual=="Aoost_01")]="DLF1+SUA"
combined_ordered3$group[which(combined_ordered3$Individual=="Aoost_08")]="DLF1+SUA"
combined_ordered3$group[which(combined_ordered3$Individual=="Banna_02")]="DLF1+SUA"
combined_ordered3$group[which(combined_ordered3$Individual=="Ctain_05")]="LJ_L7+LJ_L9"
combined_ordered3$group[which(combined_ordered3$Individual=="Llanc_09")]="LJ_H1+LJ_H2"
combined_ordered3$group[which(combined_ordered3$Individual=="Sster_01")]="LJ_H1+SUA"
combined_ordered3$group[which(combined_ordered3$Individual=="Rbani_02")]="LJ_H3+SUA"

combined_ordered3$group[which(combined_ordered3$Individual=="Aoost_09")]="DLF1"
combined_ordered3$group[which(combined_ordered3$Individual=="Banna_03")]="DLF1"
combined_ordered3$group[which(combined_ordered3$Individual=="Banna_07")]="DLF1"

combined_ordered3$group[which(combined_ordered3$Individual=="Aearl_08")]="LJ_L6"
combined_ordered3$group[which(combined_ordered3$Individual=="Ccyma_03")]="LJ_L6"
combined_ordered3$group[which(combined_ordered3$Individual=="Llanc_06")]="LJ_L6"
combined_ordered3$group[which(combined_ordered3$Individual=="Aaran_08")]="LJ_L6"

combined_ordered3$group[which(combined_ordered3$Individual=="Aearl_05")]="LJ_L7"
combined_ordered3$group[which(combined_ordered3$Individual=="Clfin_02")]="LJ_L7"
combined_ordered3$group[which(combined_ordered3$Individual=="Mrida_04")]="LJ_L7"

combined_ordered3$group[which(combined_ordered3$Individual=="Clfin_03")]="LJ_L9"
combined_ordered3$group[which(combined_ordered3$Individual=="Volin_01")]="LJ_L9"
combined_ordered3$group[which(combined_ordered3$Individual=="Aaran_04")]="LJ_L9"

combined_ordered3$group[which(combined_ordered3$Individual=="Aoost_02")]="LJ_H1"
combined_ordered3$group[which(combined_ordered3$Individual=="Ilona_09")]="LJ_H1"

combined_ordered3$group[which(combined_ordered3$Individual=="Ilona_05")]="LJ_H2"
combined_ordered3$group[which(combined_ordered3$Individual=="Kdike_09")]="LJ_H2"
combined_ordered3$group[which(combined_ordered3$Individual=="Llanc_09")]="LJ_H2"
combined_ordered3$group[which(combined_ordered3$Individual=="Aalon_03")]="LJ_H2"

combined_ordered3$group[which(combined_ordered3$Individual=="Ancor_10")]="LJ_H3"
combined_ordered3$group[which(combined_ordered3$Individual=="Borek_06")]="LJ_H3"
combined_ordered3$group[which(combined_ordered3$Individual=="Ctain_09")]="LJ_H3"

combined_ordered3$group[which(combined_ordered3$Individual=="Ancor_04")]="LJ_H5"
combined_ordered3$group[which(combined_ordered3$Individual=="Aoost_10")]="LJ_H5"
combined_ordered3$group[which(combined_ordered3$Individual=="Clfin_08")]="LJ_H5"
combined_ordered3$group[which(combined_ordered3$Individual=="Kdike_08")]="LJ_H5"

combined_ordered3$group[which(combined_ordered3$Individual=="Sster_06")]="SUA"


gpd_NoCorMeans_ordered$group=rep("not defined yet",nrow(gpd_NoCorMeans_ordered))
colnames(gpd_NoCorMeans_ordered)=c("Individual","gpd_NoCor","group")


for (i in seq(1:nrow(gpd_NoCorMeans_ordered))){
  line=gpd_NoCorMeans_ordered$Individual[i]
  group=combined_ordered3$group[which(combined_ordered3$Individual==line)][1]
  gpd_NoCorMeans_ordered$group[i]=group
}

combined_ordered3$gpd_NoCor=as.numeric(as.character(combined_ordered3$gpd_NoCor))
cairo_pdf("ParentSelection_pheno_all_values_gpdNoCor_v3.pdf", family="Arial Unicode MS",width=60,height=20)

sp<-ggplot(combined_ordered3, aes(x= Individual, y=(gpd_NoCor), color=group)) + geom_point(alpha=0.3,size=10)+
  scale_color_manual(values=c("#242c5e","black","#443659","black","#580b1d","#98d3ce","black","#71a188","#b01f35","#f6c29d","black","#f48988","grey","#426468")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) +  coord_cartesian(clip = "off") 

df_layer_1 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="Non-parent",]
df_layer_2 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="DLF1",]
df_layer_3 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_L6",]
df_layer_4 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_L7",]
df_layer_5 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_L9",]
df_layer_6 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_H1",]
df_layer_7 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_H2",]
df_layer_8 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_H3",]
df_layer_9 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_H5",]
df_layer_10 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="SUA",]
df_layer_11 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="DLF1+SUA",]
df_layer_12 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_H1+SUA",]
df_layer_13 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_H3+SUA",]
df_layer_14 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_L7+LJ_L9",]


sp1 = sp + geom_point(data=df_layer_1,aes(Individual,as.numeric(as.character(gpd_NoCor))),size=15,alpha=1,color="grey") +
  geom_point(data=df_layer_2,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#242c5e") + 
  geom_point(data=df_layer_3,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#b01f35") +
  geom_point(data=df_layer_4,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#f6c29d") +
  geom_point(data=df_layer_5,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#f48988") +
  geom_point(data=df_layer_6,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#443659") +
  geom_point(data=df_layer_7,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#580b1d") +
  geom_point(data=df_layer_8,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#98d3ce") +
  geom_point(data=df_layer_9,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#71a188") +
  geom_point(data=df_layer_10,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#426468") +
  #geom_point(data=df_layer_11,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_11,aes(Individual,gpd_NoCor), shape="\u25D6", colour="#242c5e", size=30,alpha=1) +
  geom_point(data=df_layer_11,aes(Individual,gpd_NoCor), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_12,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_12,aes(Individual,gpd_NoCor), shape="\u25D6", colour="#443659", size=30,alpha=1) +
  geom_point(data=df_layer_12,aes(Individual,gpd_NoCor), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_13,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") 
  geom_point(data=df_layer_13,aes(Individual,gpd_NoCor), shape="\u25D6", colour="#98d3ce", size=30,alpha=1) +
  geom_point(data=df_layer_13,aes(Individual,gpd_NoCor), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_14,aes(Individual,GEBVs),size=7.5,alpha=1,color="black") 
  geom_point(data=df_layer_14,aes(Individual,gpd_NoCor), shape="\u25D6", colour="#f6c29d", size=30,alpha=1) +
  geom_point(data=df_layer_14,aes(Individual,gpd_NoCor), shape="\u25D7", colour="#f48988", size=30,alpha=1) 

sp1

dev.off()

cairo_pdf("ParentSelection_pheno_means_gpdNoCor_v3.pdf", family="Arial Unicode MS",width=60,height=20)


sp<-ggplot(combined_ordered3, aes(x= Individual, y=as.numeric(as.character(gpd_NoCor)), color=group)) + geom_point(alpha=0,size=10)+
  scale_color_manual(values=c("#242c5e","black","#443659","black","#580b1d","#98d3ce","black","#71a188","#b01f35","#f6c29d","black","#f48988","grey","#426468")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) +  coord_cartesian(clip = "off") 

df_layer_1 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="Non-parent",]
df_layer_2 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="DLF1",]
df_layer_3 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_L6",]
df_layer_4 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_L7",]
df_layer_5 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_L9",]
df_layer_6 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_H1",]
df_layer_7 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_H2",]
df_layer_8 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_H3",]
df_layer_9 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_H5",]
df_layer_10 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="SUA",]
df_layer_11 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="DLF1+SUA",]
df_layer_12 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_H1+SUA",]
df_layer_13 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_H3+SUA",]
df_layer_14 <- gpd_NoCorMeans_ordered[gpd_NoCorMeans_ordered$group=="LJ_L7+LJ_L9",]


sp1 = sp + geom_point(data=df_layer_1,aes(Individual,as.numeric(as.character(gpd_NoCor))),size=15,alpha=1,color="grey") +
  geom_point(data=df_layer_2,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#242c5e") + 
  geom_point(data=df_layer_3,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#b01f35") +
  geom_point(data=df_layer_4,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#f6c29d") +
  geom_point(data=df_layer_5,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#f48988") +
  geom_point(data=df_layer_6,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#443659") +
  geom_point(data=df_layer_7,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#580b1d") +
  geom_point(data=df_layer_8,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#98d3ce") +
  geom_point(data=df_layer_9,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#71a188") +
  geom_point(data=df_layer_10,aes(Individual,gpd_NoCor),size=15,alpha=1,color="#426468") +
  #geom_point(data=df_layer_11,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_11,aes(Individual,gpd_NoCor), shape="\u25D6", colour="#242c5e", size=30,alpha=1) +
  geom_point(data=df_layer_11,aes(Individual,gpd_NoCor), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_12,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_12,aes(Individual,gpd_NoCor), shape="\u25D6", colour="#443659", size=30,alpha=1) +
  geom_point(data=df_layer_12,aes(Individual,gpd_NoCor), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_13,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") 
  geom_point(data=df_layer_13,aes(Individual,gpd_NoCor), shape="\u25D6", colour="#98d3ce", size=30,alpha=1) +
  geom_point(data=df_layer_13,aes(Individual,gpd_NoCor), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_14,aes(Individual,GEBVs),size=7.5,alpha=1,color="black") 
  geom_point(data=df_layer_14,aes(Individual,gpd_NoCor), shape="\u25D6", colour="#f6c29d", size=30,alpha=1) +
  geom_point(data=df_layer_14,aes(Individual,gpd_NoCor), shape="\u25D7", colour="#f48988", size=30,alpha=1) 

sp1

dev.off()





# Figure 4
# y-values are phenotypes (gpd_ResCor)
phenotypeperf_gpdResCor=cbind(as.character(d6data$Clovershort),as.numeric(as.character(d6data$gpd_dryweight_cor)))
colnames(phenotypeperf_gpdResCor)=c("Individual","gpd_ResCor")
phenotypeperf_gpdResCor=as.data.frame(phenotypeperf_gpdResCor)


gpd_ResCorMeans=aggregate(as.numeric(as.character((phenotypeperf_gpdResCor$gpd_ResCor))), list(phenotypeperf_gpdResCor$Individual), mean)
gpd_ResCorMeans_ordered <- gpd_ResCorMeans[order(gpd_ResCorMeans[,2]),]
individual_order_gpdResCorPheno=gpd_ResCorMeans_ordered[,1]

keyDF <- data.frame(key=individual_order_gpdResCorPheno,weight=1:length(individual_order_gpdResCorPheno))
combined_ready4 <- merge(phenotypeperf_gpdResCor,keyDF,by.x='Individual',by.y='key',all.x=T,all.y=F)
combined_ordered4 <- combined_ready4[order(combined_ready4$weight),c('Individual','gpd_ResCor','weight')]

combined_ordered4$Individual <- factor(combined_ordered4$Individual,levels = individual_order_gpdResCorPheno)

combined_ordered4$group=as.character(rep("Non-parent",nrow(combined_ordered4)))
combined_ordered4$group[which(combined_ordered4$Individual=="Aoost_01")]="DLF1+SUA"
combined_ordered4$group[which(combined_ordered4$Individual=="Aoost_08")]="DLF1+SUA"
combined_ordered4$group[which(combined_ordered4$Individual=="Banna_02")]="DLF1+SUA"
combined_ordered4$group[which(combined_ordered4$Individual=="Ctain_05")]="LJ_L7+LJ_L9"
combined_ordered4$group[which(combined_ordered4$Individual=="Llanc_09")]="LJ_H1+LJ_H2"
combined_ordered4$group[which(combined_ordered4$Individual=="Sster_01")]="LJ_H1+SUA"
combined_ordered4$group[which(combined_ordered4$Individual=="Rbani_02")]="LJ_H3+SUA"

combined_ordered4$group[which(combined_ordered4$Individual=="Aoost_09")]="DLF1"
combined_ordered4$group[which(combined_ordered4$Individual=="Banna_03")]="DLF1"
combined_ordered4$group[which(combined_ordered4$Individual=="Banna_07")]="DLF1"

combined_ordered4$group[which(combined_ordered4$Individual=="Aearl_08")]="LJ_L6"
combined_ordered4$group[which(combined_ordered4$Individual=="Ccyma_03")]="LJ_L6"
combined_ordered4$group[which(combined_ordered4$Individual=="Llanc_06")]="LJ_L6"
combined_ordered4$group[which(combined_ordered4$Individual=="Aaran_08")]="LJ_L6"

combined_ordered4$group[which(combined_ordered4$Individual=="Aearl_05")]="LJ_L7"
combined_ordered4$group[which(combined_ordered4$Individual=="Clfin_02")]="LJ_L7"
combined_ordered4$group[which(combined_ordered4$Individual=="Mrida_04")]="LJ_L7"

combined_ordered4$group[which(combined_ordered4$Individual=="Clfin_03")]="LJ_L9"
combined_ordered4$group[which(combined_ordered4$Individual=="Volin_01")]="LJ_L9"
combined_ordered4$group[which(combined_ordered4$Individual=="Aaran_04")]="LJ_L9"

combined_ordered4$group[which(combined_ordered4$Individual=="Aoost_02")]="LJ_H1"
combined_ordered4$group[which(combined_ordered4$Individual=="Ilona_09")]="LJ_H1"

combined_ordered4$group[which(combined_ordered4$Individual=="Ilona_05")]="LJ_H2"
combined_ordered4$group[which(combined_ordered4$Individual=="Kdike_09")]="LJ_H2"
combined_ordered4$group[which(combined_ordered4$Individual=="Llanc_09")]="LJ_H2"
combined_ordered4$group[which(combined_ordered4$Individual=="Aalon_03")]="LJ_H2"

combined_ordered4$group[which(combined_ordered4$Individual=="Ancor_10")]="LJ_H3"
combined_ordered4$group[which(combined_ordered4$Individual=="Borek_06")]="LJ_H3"
combined_ordered4$group[which(combined_ordered4$Individual=="Ctain_09")]="LJ_H3"

combined_ordered4$group[which(combined_ordered4$Individual=="Ancor_04")]="LJ_H5"
combined_ordered4$group[which(combined_ordered4$Individual=="Aoost_10")]="LJ_H5"
combined_ordered4$group[which(combined_ordered4$Individual=="Clfin_08")]="LJ_H5"
combined_ordered4$group[which(combined_ordered4$Individual=="Kdike_08")]="LJ_H5"

combined_ordered4$group[which(combined_ordered4$Individual=="Sster_06")]="SUA"


gpd_ResCorMeans_ordered$group=rep("not defined yet",nrow(gpd_ResCorMeans_ordered))
colnames(gpd_ResCorMeans_ordered)=c("Individual","gpd_ResCor","group")


for (i in seq(1:nrow(gpd_ResCorMeans_ordered))){
  line=gpd_ResCorMeans_ordered$Individual[i]
  group=combined_ordered4$group[which(combined_ordered4$Individual==line)][1]
  gpd_ResCorMeans_ordered$group[i]=group
}


cairo_pdf("ParentSelection_pheno_all_values_gpdResCor_v3.pdf", family="Arial Unicode MS",width=60,height=20)

sp<-ggplot(combined_ordered4, aes(x= Individual, y=as.numeric(as.character(gpd_ResCor)), color=group)) + geom_point(alpha=0.3,size=10)+
  scale_color_manual(values=c("#242c5e","black","#443659","black","#580b1d","#98d3ce","black","#71a188","#b01f35","#f6c29d","black","#f48988","grey","#426468")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) +  coord_cartesian(clip = "off") 

df_layer_1 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="Non-parent",]
df_layer_2 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="DLF1",]
df_layer_3 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_L6",]
df_layer_4 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_L7",]
df_layer_5 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_L9",]
df_layer_6 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_H1",]
df_layer_7 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_H2",]
df_layer_8 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_H3",]
df_layer_9 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_H5",]
df_layer_10 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="SUA",]
df_layer_11 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="DLF1+SUA",]
df_layer_12 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_H1+SUA",]
df_layer_13 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_H3+SUA",]
df_layer_14 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_L7+LJ_L9",]


sp1 = sp + geom_point(data=df_layer_1,aes(Individual,gpd_ResCor),size=15,alpha=1,color="grey") +
  geom_point(data=df_layer_2,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#242c5e") + 
  geom_point(data=df_layer_3,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#b01f35") +
  geom_point(data=df_layer_4,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#f6c29d") +
  geom_point(data=df_layer_5,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#f48988") +
  geom_point(data=df_layer_6,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#443659") +
  geom_point(data=df_layer_7,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#580b1d") +
  geom_point(data=df_layer_8,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#98d3ce") +
  geom_point(data=df_layer_9,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#71a188") +
  geom_point(data=df_layer_10,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#426468") +
  #geom_point(data=df_layer_11,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_11,aes(Individual,gpd_ResCor), shape="\u25D6", colour="#242c5e", size=30,alpha=1) +
  geom_point(data=df_layer_11,aes(Individual,gpd_ResCor), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_12,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_12,aes(Individual,gpd_ResCor), shape="\u25D6", colour="#443659", size=30,alpha=1) +
  geom_point(data=df_layer_12,aes(Individual,gpd_ResCor), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_13,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") 
  geom_point(data=df_layer_13,aes(Individual,gpd_ResCor), shape="\u25D6", colour="#98d3ce", size=30,alpha=1) +
  geom_point(data=df_layer_13,aes(Individual,gpd_ResCor), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_14,aes(Individual,GEBVs),size=7.5,alpha=1,color="black") 
  geom_point(data=df_layer_14,aes(Individual,gpd_ResCor), shape="\u25D6", colour="#f6c29d", size=30,alpha=1) +
  geom_point(data=df_layer_14,aes(Individual,gpd_ResCor), shape="\u25D7", colour="#f48988", size=30,alpha=1) 

sp1

dev.off()

cairo_pdf("ParentSelection_pheno_means_gpdResCor_v3.pdf", family="Arial Unicode MS",width=60,height=20)


sp<-ggplot(combined_ordered4, aes(x= Individual, y=as.numeric(as.character(gpd_ResCor)), color=group)) + geom_point(alpha=0,size=10)+
  scale_color_manual(values=c("#242c5e","black","#443659","black","#580b1d","#98d3ce","black","#71a188","#b01f35","#f6c29d","black","#f48988","grey","#426468")) + 
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) +  coord_cartesian(clip = "off") 

df_layer_1 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="Non-parent",]
df_layer_2 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="DLF1",]
df_layer_3 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_L6",]
df_layer_4 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_L7",]
df_layer_5 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_L9",]
df_layer_6 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_H1",]
df_layer_7 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_H2",]
df_layer_8 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_H3",]
df_layer_9 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_H5",]
df_layer_10 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="SUA",]
df_layer_11 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="DLF1+SUA",]
df_layer_12 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_H1+SUA",]
df_layer_13 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_H3+SUA",]
df_layer_14 <- gpd_ResCorMeans_ordered[gpd_ResCorMeans_ordered$group=="LJ_L7+LJ_L9",]


sp1 = sp + geom_point(data=df_layer_1,aes(Individual,as.numeric(as.character(gpd_ResCor))),size=15,alpha=1,color="grey") +
  geom_point(data=df_layer_2,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#242c5e") + 
  geom_point(data=df_layer_3,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#b01f35") +
  geom_point(data=df_layer_4,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#f6c29d") +
  geom_point(data=df_layer_5,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#f48988") +
  geom_point(data=df_layer_6,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#443659") +
  geom_point(data=df_layer_7,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#580b1d") +
  geom_point(data=df_layer_8,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#98d3ce") +
  geom_point(data=df_layer_9,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#71a188") +
  geom_point(data=df_layer_10,aes(Individual,gpd_ResCor),size=15,alpha=1,color="#426468") +
  #geom_point(data=df_layer_11,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_11,aes(Individual,gpd_ResCor), shape="\u25D6", colour="#242c5e", size=30,alpha=1) +
  geom_point(data=df_layer_11,aes(Individual,gpd_ResCor), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_12,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") +
  geom_point(data=df_layer_12,aes(Individual,gpd_ResCor), shape="\u25D6", colour="#443659", size=30,alpha=1) +
  geom_point(data=df_layer_12,aes(Individual,gpd_ResCor), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_13,aes(Individual,GEBVs),size=7.5,alpha=0.9,color="black") 
  geom_point(data=df_layer_13,aes(Individual,gpd_ResCor), shape="\u25D6", colour="#98d3ce", size=30,alpha=1) +
  geom_point(data=df_layer_13,aes(Individual,gpd_ResCor), shape="\u25D7", colour="#426468", size=30,alpha=1) +
  #geom_point(data=df_layer_14,aes(Individual,GEBVs),size=7.5,alpha=1,color="black") 
  geom_point(data=df_layer_14,aes(Individual,gpd_ResCor), shape="\u25D6", colour="#f6c29d", size=30,alpha=1) +
  geom_point(data=df_layer_14,aes(Individual,gpd_ResCor), shape="\u25D7", colour="#f48988", size=30,alpha=1) 

sp1

dev.off()


# make a highlighted plot but one for each parental group
# Not used in manuscript
png(file="AllParentsMarkedDiscretelyGroupsCanBeSeenInd.png",width=210,height=180,units="cm",res=300)

p1=ggplot(combined_ordered, aes(x= Individual, y=GEBVs, color=gpd_measurement)) + geom_point(alpha=0.4,size=5)+
  gghighlight(Individual %in% parentpop1names) +ggtitle("DLF1") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p2=ggplot(combined_ordered, aes(x= Individual, y=GEBVs, color=gpd_measurement)) + geom_point(alpha=0.4,size=5)+
  gghighlight(Individual %in% parentpop6names) +ggtitle("LJ L6") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p3=ggplot(combined_ordered, aes(x= Individual, y=GEBVs, color=gpd_measurement)) + geom_point(alpha=0.4,size=5)+
  gghighlight(Individual %in% parentpop7names) +ggtitle("LJ L7") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p4=ggplot(combined_ordered, aes(x= Individual, y=GEBVs, color=gpd_measurement)) + geom_point(alpha=0.4,size=5)+
  gghighlight(Individual %in% parentpop8names) +ggtitle("LJ L9") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p5=ggplot(combined_ordered, aes(x= Individual, y=GEBVs, color=gpd_measurement)) + geom_point(alpha=0.4,size=5)+
  gghighlight(Individual %in% parentpop9names) +ggtitle("LJ H1") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p6=ggplot(combined_ordered, aes(x= Individual, y=GEBVs, color=gpd_measurement)) + geom_point(alpha=0.4,size=5)+
  gghighlight(Individual %in% parentpop10names) +ggtitle("LJ H2") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p7=ggplot(combined_ordered, aes(x= Individual, y=GEBVs, color=gpd_measurement)) + geom_point(alpha=0.4,size=5)+
  gghighlight(Individual %in% parentpop11names) +ggtitle("LJ H3") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p8=ggplot(combined_ordered, aes(x= Individual, y=GEBVs, color=gpd_measurement)) + geom_point(alpha=0.4,size=5)+
  gghighlight(Individual %in% parentpop12names) +ggtitle("LJ H5") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p9=ggplot(combined_ordered, aes(x= Individual, y=GEBVs, color=gpd_measurement)) + geom_point(alpha=0.4,size=5)+
  gghighlight(Individual %in% parentpop13names) +ggtitle("SUA") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))


grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)
dev.off()


# Same as above but with averages
# Not used in manuscript
png(file="AllParentsMarkedDiscretelyGroupsCanBeSeenInd_averages.png",width=200,height=140,units="cm",res=300)

p1=ggplot(combined_ordered1, aes(x= Individual, y=GEBVs, color=gpd_measurement,shape=gpd_measurement)) + geom_point(alpha=0.7,size=15)+
  gghighlight(Individual %in% parentpop1names, use_direct_label = FALSE,use_group_by = FALSE) +
  ggtitle("DLF1") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + scale_shape_manual(values=c(15,16,18)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 6),axis.text.y = element_text(angle = 90, size = 12), axis.title=element_text(size=22),
          plot.title = element_text(size = 35, face = "bold")) +
  coord_cartesian(clip = "off") # so it doesm't cut off parts of points to the right

p2=ggplot(combined_ordered1, aes(x= Individual, y=GEBVs, color=gpd_measurement,shape=gpd_measurement,shape=gpd_measurement)) + geom_point(alpha=0.7,size=15)+
  gghighlight(Individual %in% parentpop6names, use_direct_label = FALSE,use_group_by = FALSE) +ggtitle("LJ L6") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + scale_shape_manual(values=c(15,16,18)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 6),axis.text.y = element_text(angle = 90, size = 12), axis.title=element_text(size=22),
          plot.title = element_text(size = 35, face = "bold")) +
  coord_cartesian(clip = "off")

p3=ggplot(combined_ordered1, aes(x= Individual, y=GEBVs, color=gpd_measurement,shape=gpd_measurement)) + geom_point(alpha=0.7,size=15)+
  gghighlight(Individual %in% parentpop7names, use_direct_label = FALSE,use_group_by = FALSE) +ggtitle("LJ L7") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) +  scale_shape_manual(values=c(15,16,18)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 6),axis.text.y = element_text(angle = 90, size = 12), axis.title=element_text(size=22),
          plot.title = element_text(size = 35, face = "bold")) +
  coord_cartesian(clip = "off")

p4=ggplot(combined_ordered1, aes(x= Individual, y=GEBVs, color=gpd_measurement,shape=gpd_measurement)) + geom_point(alpha=0.7,size=15)+
  gghighlight(Individual %in% parentpop8names, use_direct_label = FALSE,use_group_by = FALSE) +ggtitle("LJ L9") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + scale_shape_manual(values=c(15,16,18)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 6),axis.text.y = element_text(angle = 90, size = 12), axis.title=element_text(size=22),
          plot.title = element_text(size = 35, face = "bold")) +
  coord_cartesian(clip = "off")

p5=ggplot(combined_ordered1, aes(x= Individual, y=GEBVs, color=gpd_measurement,shape=gpd_measurement)) + geom_point(alpha=0.7,size=15)+
  gghighlight(Individual %in% parentpop9names, use_direct_label = FALSE,use_group_by = FALSE) +ggtitle("LJ H1") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) +  scale_shape_manual(values=c(15,16,18)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 6),axis.text.y = element_text(angle = 90, size = 12), axis.title=element_text(size=22),
          plot.title = element_text(size = 35, face = "bold")) +
  coord_cartesian(clip = "off")

p6=ggplot(combined_ordered1, aes(x= Individual, y=GEBVs, color=gpd_measurement,shape=gpd_measurement)) + geom_point(alpha=0.7,size=15)+
  gghighlight(Individual %in% parentpop10names, use_direct_label = FALSE,use_group_by = FALSE) +ggtitle("LJ H2") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) +  scale_shape_manual(values=c(15,16,18)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 6),axis.text.y = element_text(angle = 90, size = 12), axis.title=element_text(size=22),
          plot.title = element_text(size = 35, face = "bold")) +
  coord_cartesian(clip = "off")

p7=ggplot(combined_ordered1, aes(x= Individual, y=GEBVs, color=gpd_measurement,shape=gpd_measurement)) + geom_point(alpha=0.7,size=15)+
  gghighlight(Individual %in% parentpop11names, use_direct_label = FALSE,use_group_by = FALSE) +ggtitle("LJ H3") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) +  scale_shape_manual(values=c(15,16,18)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 6),axis.text.y = element_text(angle = 90, size = 12), axis.title=element_text(size=22),
          plot.title = element_text(size = 35, face = "bold")) +
  coord_cartesian(clip = "off")

p8=ggplot(combined_ordered1, aes(x= Individual, y=GEBVs, color=gpd_measurement,shape=gpd_measurement)) + geom_point(alpha=0.7,size=15)+
  gghighlight(Individual %in% parentpop12names, use_direct_label = FALSE,use_group_by = FALSE) +ggtitle("LJ H5") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) +  scale_shape_manual(values=c(15,16,18)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 6),axis.text.y = element_text(angle = 90, size = 12), axis.title=element_text(size=22),
          plot.title = element_text(size = 35, face = "bold")) +
  coord_cartesian(clip = "off")

p9=ggplot(combined_ordered1, aes(x= Individual, y=GEBVs, color=gpd_measurement,shape=gpd_measurement)) + geom_point(alpha=0.7,size=15)+
  gghighlight(Individual %in% parentpop13names, use_direct_label = FALSE,use_group_by = FALSE) +ggtitle("SUA") +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,4,3)]) + scale_shape_manual(values=c(15,16,18)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 6),axis.text.y = element_text(angle = 90, size = 12), axis.title=element_text(size=22),
          plot.title = element_text(size = 35, face = "bold")) +
  coord_cartesian(clip = "off")


grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)
dev.off()


# gpd_ResCor versus different GEBVs
#in mean_merged every GEBV corresponds to the trait it is measured for
# The gpd_rescor corrected for fixed effects are the one that corresponds to the gpd_ResCor
# Not used in manuscript
# can be deleted if corrgram do not depend on it
mean_merged_only_with_gpdResCor_phenotypes=mean_merged

for (names in unique(mean_merged_only_with_gpdResCor_phenotypes$Individual)){
  indexes=which(mean_merged_only_with_gpdResCor_phenotypes$Individual==names)
  gpd_ResCor=phenotype=mean_merged_only_with_gpdResCor_phenotypes$Observed_correctedForFixed[indexes[1]]
  mean_merged_only_with_gpdResCor_phenotypes$Observed_correctedForFixed[indexes]=gpd_ResCor
}
mean_merged_only_with_gpdResCor_phenotypes_ordered <- mean_merged_only_with_gpdResCor_phenotypes[order(mean_merged_only_with_gpdResCor_phenotypes$Observed_correctedForFixed),]

png(file="gpd_ResCor_vs_differentGEBVs",width=70,height=60,units="cm",res=300)

sp<-ggplot(mean_merged_only_with_gpdResCor_phenotypes_ordered, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,3,4)])  +
  labs(x="GEBV", y = "gpd_ResCor_Corrected") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) 
dev.off()


# gpd_ResCor versus different GEBVs, 3 plots
# can be deleted if corrgram do not depend on it
# Not used in manuscript
Res_Subset_phenogpdrescor=subset(mean_merged_only_with_gpdResCor_phenotypes_ordered,gpd_measurement=="gpd_ResCor")
No_Subset_phenogpdrescor=subset(mean_merged_only_with_gpdResCor_phenotypes_ordered,gpd_measurement=="gpd_NoCor")
Fix_Subset_phenogpdrescor=subset(mean_merged_only_with_gpdResCor_phenotypes_ordered,gpd_measurement=="gpd_FixCor")

correlation(Res_Subset_phenogpdrescor$GEBVs,Res_Subset_phenogpdrescor$Observed_correctedForFixed)
correlation(No_Subset_phenogpdrescor$GEBVs,No_Subset_phenogpdrescor$Observed_correctedForFixed)
correlation(Fix_Subset_phenogpdrescor$GEBVs,Fix_Subset_phenogpdrescor$Observed_correctedForFixed)

png(file="gpd_ResCor_vs_IndividualGEBVs.png",width=70,height=100,units="cm",res=300)


p1= ggplot(Res_Subset_phenogpdrescor, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[3])  +
  labs(title="GEBVs fromgpd_ResCor (cor=0.29, p < 0.001)",x="GEBV", y = "gpd_ResCor_Corrected") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20),plot.title = element_text(color = "black", size = 40, face = "bold")) 


p2= ggplot(No_Subset_phenogpdrescor, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[4])  +
  labs(title="GEBVs from gpd_NoCor (cor=0.28, p < 0.001)",x="GEBV", y = "gpd_ResCor_Corrected") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20),plot.title = element_text(color = "black", size = 40, face = "bold")) 

p3= ggplot(Fix_Subset_phenogpdrescor, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[2])  +
  labs(title="GEBVs from gpd_FixCor (cor=0.04, p > 0.05)",x="GEBV", y = "gpd_ResCor_Corrected") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20),plot.title = element_text(color = "black", size = 40, face = "bold")) 

grid.arrange(p1, p2, p3, nrow = 3)
dev.off()



# gpd_NoCor versus different GEBVs, 3 plots
# Not used in manuscript

setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd")

mean_merged_only_with_gpdNoCor_phenotypes=mean_merged

for (names in unique(mean_merged_only_with_gpdNoCor_phenotypes$Individual)){
  indexes=which(mean_merged_only_with_gpdNoCor_phenotypes$Individual==names)
  gpd_ResCor=phenotype=mean_merged_only_with_gpdNoCor_phenotypes$Observed_correctedForFixed[indexes[3]]
  mean_merged_only_with_gpdNoCor_phenotypes$Observed_correctedForFixed[indexes]=gpd_ResCor
}

mean_merged_only_with_gpdNoCor_phenotypes_ordered <- mean_merged_only_with_gpdNoCor_phenotypes[order(mean_merged_only_with_gpdNoCor_phenotypes$Observed_correctedForFixed),]


png(file="gpd_NoCor_vs_differentGEBVs.png",width=70,height=60,units="cm",res=300)

sp<-ggplot(mean_merged_only_with_gpdNoCor_phenotypes_ordered, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,3,4)])  +
  labs(x="GEBV", y = "gpd_NoCor_Corrected") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) 

sp
dev.off()

# gpd_noCor versus different GEBVs, 3 plots
# Not used in manuscript

Res_Subset_phenogpdNocor=subset(mean_merged_only_with_gpdNoCor_phenotypes_ordered,gpd_measurement=="gpd_ResCor")
No_Subset_phenogpdNocor=subset(mean_merged_only_with_gpdNoCor_phenotypes_ordered,gpd_measurement=="gpd_NoCor")
Fix_Subset_phenogpdNocor=subset(mean_merged_only_with_gpdNoCor_phenotypes_ordered,gpd_measurement=="gpd_FixCor")

correlation(Res_Subset_phenogpdNocor$GEBVs,Res_Subset_phenogpdNocor$Observed_correctedForFixed)
correlation(No_Subset_phenogpdNocor$GEBVs,No_Subset_phenogpdNocor$Observed_correctedForFixed)
correlation(Fix_Subset_phenogpdNocor$GEBVs,Fix_Subset_phenogpdNocor$Observed_correctedForFixed)

png(file="gpd_NoCor_vs_IndividualGEBVs.png",width=70,height=100,units="cm",res=300)


p1= ggplot(Res_Subset_phenogpdNocor, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[3])  +
  labs(title="GEBVs fromgpd_ResCor (cor=0.23, p < 0.01)",x="GEBV", y = "gpd_NoCor_Corrected") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20),plot.title = element_text(color = "black", size = 40, face = "bold")) 


p2= ggplot(No_Subset_phenogpdNocor, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[4])  +
  labs(title="GEBVs from gpd_NoCor (cor=0.22, p < 0.05)",x="GEBV", y = "gpd_NoCor_Corrected") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20),plot.title = element_text(color = "black", size = 40, face = "bold")) 

p3= ggplot(Fix_Subset_phenogpdNocor, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[2])  +
  labs(title="GEBVs from gpd_FixCor (cor=-0.02, p > 0.05)",x="GEBV", y = "gpd_NoCor_Corrected") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20),plot.title = element_text(color = "black", size = 40, face = "bold")) 

grid.arrange(p1, p2, p3, nrow = 3)
dev.off()


# gpd_FixCor versus different GEBVs, 3 plots
# Not used in manuscript

mean_merged_only_with_gpdFixCor_phenotypes=mean_merged
for (names in unique(mean_merged_only_with_gpdFixCor_phenotypes$Individual)){
  indexes=which(mean_merged_only_with_gpdFixCor_phenotypes$Individual==names)
  gpd_ResCor=phenotype=mean_merged_only_with_gpdFixCor_phenotypes$Observed_correctedForFixed[indexes[2]]
  mean_merged_only_with_gpdFixCor_phenotypes$Observed_correctedForFixed[indexes]=gpd_ResCor
}

mean_merged_only_with_gpdFixCor_phenotypes_ordered <- mean_merged_only_with_gpdFixCor_phenotypes[order(mean_merged_only_with_gpdFixCor_phenotypes$Observed_correctedForFixed),]


png(file="gpd_FixCor_vs_differentGEBVs.png",width=70,height=60,units="cm",res=300)

sp<-ggplot(mean_merged_only_with_gpdFixCor_phenotypes_ordered, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[c(2,3,4)])  +
  labs(x="GEBV", y = "gpd_FixCor_Corrected") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20)) 

sp
dev.off()

# gpd_noCor versus different GEBVs, 3 plots
# Not used in manuscript

Res_Subset_phenogpdFixcor=subset(mean_merged_only_with_gpdFixCor_phenotypes_ordered,gpd_measurement=="gpd_ResCor")
No_Subset_phenogpdFixcor=subset(mean_merged_only_with_gpdFixCor_phenotypes_ordered,gpd_measurement=="gpd_NoCor")
Fix_Subset_phenogpdFixcor=subset(mean_merged_only_with_gpdFixCor_phenotypes_ordered,gpd_measurement=="gpd_FixCor")

correlation(Res_Subset_phenogpdFixcor$GEBVs,Res_Subset_phenogpdFixcor$Observed_correctedForFixed)
correlation(No_Subset_phenogpdFixcor$GEBVs,No_Subset_phenogpdFixcor$Observed_correctedForFixed)
correlation(Fix_Subset_phenogpdFixcor$GEBVs,Fix_Subset_phenogpdFixcor$Observed_correctedForFixed)

png(file="gpd_FixCor_vs_IndividualGEBVs.png",width=70,height=100,units="cm",res=300)


p1= ggplot(Res_Subset_phenogpdFixcor, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[3])  +
  labs(title="GEBVs fromgpd_ResCor (cor=0.05, p > 0.05)",x="GEBV", y = "gpd_FixCor_Corrected") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20),plot.title = element_text(color = "black", size = 40, face = "bold")) 


p2= ggplot(No_Subset_phenogpdFixcor, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[4])  +
  labs(title="GEBVs from gpd_NoCor (cor=0.02, p > 0.05)",x="GEBV", y = "gpd_FixCor_Corrected") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20),plot.title = element_text(color = "black", size = 40, face = "bold")) 

p3= ggplot(Fix_Subset_phenogpdFixcor, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[2])  +
  labs(title="GEBVs from gpd_FixCor (cor=-0.15, p > 0.05)",x="GEBV", y = "gpd_FixCor_Corrected") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20),plot.title = element_text(color = "black", size = 40, face = "bold")) 

grid.arrange(p1, p2, p3, nrow = 3)
dev.off()




# Phenotype averages versus different GEBVs, 3 plots
# Not used in manuscript

full_data_table=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/JustTheWhole_d6Table.txt",head=T)
colnames(full_data_table)

Phenotypic_averages=aggregate(as.numeric(full_data_table$growth_per_day), list(full_data_table$Clovershort), mean)
colnames(Phenotypic_averages)=c("Individual","gpd_avg")
  
PhenoAvg_GEBVsFromRescor=merge(Res_Subset_phenogpdFixcor,Phenotypic_averages,by="Individual")
PhenoAvg_GEBVsFromFixcor=merge(Fix_Subset_phenogpdFixcor,Phenotypic_averages,by="Individual")
PhenoAvg_GEBVsFromNocor=merge(No_Subset_phenogpdNocor,Phenotypic_averages,by="Individual")

correlation(PhenoAvg_GEBVsFromRescor$GEBVs,PhenoAvg_GEBVsFromRescor$gpd_avg)
correlation(PhenoAvg_GEBVsFromFixcor$GEBVs,PhenoAvg_GEBVsFromFixcor$gpd_avg)
correlation(PhenoAvg_GEBVsFromNocor$GEBVs,PhenoAvg_GEBVsFromNocor$gpd_avg)

png(file="Average_phenotypes_vs_IndividualGEBVs.png",width=70,height=100,units="cm",res=300)


p1= ggplot(PhenoAvg_GEBVsFromRescor, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[3])  +
  labs(title="GEBVs fromgpd_ResCor (cor=0.32, p < 0.001)",x="GEBV", y = "Avg. gpd of genotype") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20),plot.title = element_text(color = "black", size = 40, face = "bold")) 


p2= ggplot(PhenoAvg_GEBVsFromNocor, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[4])  +
  labs(title="GEBVs from gpd_NoCor (cor=0.07, p > 0.05)",x="GEBV", y = "Avg. gpd of genotype") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20),plot.title = element_text(color = "black", size = 40, face = "bold")) 

p3= ggplot(PhenoAvg_GEBVsFromFixcor, aes(x= GEBVs, y=Observed_correctedForFixed, color=gpd_measurement)) + geom_point(alpha=0.7,size=5) +
  scale_color_manual(values=wes_palette("Rushmore1")[2])  +
  labs(title="GEBVs from gpd_FixCor (cor=0.30, p < 0.001)",x="GEBV", y = "Avg. gpd of genotype") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20),plot.title = element_text(color = "black", size = 40, face = "bold")) 

grid.arrange(p1, p2, p3, nrow = 3)
dev.off()




############ Make F1 prediction plot, gpd_rescor

F1predictions=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/F1_prediction_20200828/gpd_rescor/AvgParentalGEBV.txt",head=T)
head(F1predictions)
dim(F1predictions)
F1predictions_t=t(F1predictions)
F1predictions_t=as.data.frame(F1predictions_t)
F1predictions_t_=F1predictions_t[2:nrow(F1predictions_t),]
F1predictions_t_$Population=as.character(row.names(F1predictions_t_))
rownames(F1predictions_t_)=NULL
colnames(F1predictions_t_)=c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","Population")
head(F1predictions_t_)


F1observations=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/F1_prediction/With Lflex/F1observations.txt",head=T)
head(F1observations)
 
mergedF1info=merge(F1predictions_t_,F1observations,by="Population")
mergedF1info_onlyPred=mergedF1info[,1:(ncol(mergedF1info)-2)]
mergedF1info$Means=rowMeans(mergedF1info_onlyPred[sapply(mergedF1info_onlyPred, is.numeric)])


correlation(as.numeric(as.character(mergedF1info$Means)),as.numeric(as.character(mergedF1info$Avg.Dryweight.F1)))

mergedF1info
mergedF1info$SD=apply(mergedF1info_onlyPred[,-1],1,sd)
  
# Not used in manuscript
#png(file="F1Prediction_avg_gpd_ResCor.png",width=20,height=14,units="cm",res=300, useDingbats=FALSE)

sp<-ggplot(mergedF1info, aes(x= Means, y=Avg.Dryweight.F1)) + geom_point(alpha=0.8,size=5,color=wes_palette("Rushmore1")[3]) +
labs(x="Avg. parental GEBV", y = "Avg. dry weight of population") +
  xlim(c(min(mergedF1info$Means),max(mergedF1info$Means)))+
  geom_smooth(method='lm',se = FALSE, color="black") +
  geom_text(x=0.012, y=3, label="R",fontface='italic') +
  geom_text(x=0.018, y=3, label="= 0.93 (p < 0.001)") +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp

ggsave("F1Prediction_avg_gpd_ResCor.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)

#dev.off()

# With labels, not used in manuscript
#png(file="F1Prediction_avg_withlabels_gpdResCor.png",width=20,height=14,units="cm",res=300)
sp<-ggplot(mergedF1info, aes(x= Means, y=Avg.Dryweight.F1)) + geom_point(alpha=0.8,size=5,color=wes_palette("Rushmore1")[3]) +
  labs(x="Avg. parental GEBV", y = "Avg. dry weight of population") +
  xlim(c(min(mergedF1info$Means),max(mergedF1info$Means)))+
  geom_text(x=0.012, y=3, label="R",fontface='italic') +
  geom_text(x=0.018, y=3, label="= 0.93 (p < 0.001)") +
  geom_smooth(method='lm',se = FALSE, color="black") +
  geom_label_repel(aes(label = mergedF1info$Population),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp
ggsave("F1Prediction_avg_withlabels_gpdResCor.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)
#dev.off()

# Figure 6 in manuscript
colors=c("#242c5e","#443659","#580b1d","#98d3ce","#71a188","#b01f35","#f6c29d","#f48988","#426468")
sp<-ggplot(mergedF1info, aes(x= Means, y=Avg.Dryweight.F1, colour=(Population))) + geom_point(alpha=1.0,size=5) +
  scale_color_manual(values=colors) +
  labs(x="Avg. parental GEBV", y = "Avg. dry weight of population") +
  xlim(c(min(mergedF1info$Means),max(mergedF1info$Means)))+
  #geom_text(x=0.012, y=3, label="R",fontface='italic') +
  #geom_text(x=0.016, y=3, label="= 0.93 (p < 0.001)") +
  geom_smooth(method='lm',se = FALSE, color="black") +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp
ggsave("F1Prediction_avg_colByPop_gpdResCor.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)

############ Make F1 prediction plot, gpd_Nocor

F1predictions=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/F1_prediction_20200828/gpd_nocor/AvgParentalGEBV.txt",head=T)
head(F1predictions)
dim(F1predictions)
F1predictions_t=t(F1predictions)
F1predictions_t=as.data.frame(F1predictions_t)
F1predictions_t_=F1predictions_t[2:nrow(F1predictions_t),]
F1predictions_t_$Population=as.character(row.names(F1predictions_t_))
rownames(F1predictions_t_)=NULL
colnames(F1predictions_t_)=c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","Population")
head(F1predictions_t_)


F1observations=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/F1_prediction/With Lflex/F1observations.txt",head=T)
head(F1observations)

mergedF1info=merge(F1predictions_t_,F1observations,by="Population")
mergedF1info_onlyPred=mergedF1info[,1:(ncol(mergedF1info)-2)]
mergedF1info$Means=rowMeans(mergedF1info_onlyPred[sapply(mergedF1info_onlyPred, is.numeric)])

mergedF1info_No=mergedF1info
correlation(as.numeric(as.character(mergedF1info$Means)),as.numeric(as.character(mergedF1info$Avg.Dryweight.F1))) #0.91

mergedF1info
mergedF1info$SD=apply(mergedF1info_onlyPred[,-1],1,sd)

# Not used in manuscript
#png(file="F1Prediction_avg_gpd_NoCor.png",width=20,height=14,units="cm",res=300)

sp<-ggplot(mergedF1info, aes(x= Means, y=Avg.Dryweight.F1)) + geom_point(alpha=0.8,size=5,color=wes_palette("Rushmore1")[4]) +
  labs(x="Avg. parental GEBV", y = "Avg. dry weight of population") +
  xlim(c(min(mergedF1info$Means),max(mergedF1info$Means)))+
  geom_text(x=0.01, y=3, label="R",fontface='italic') +
  geom_text(x=0.015, y=3, label="= 0.91 (p < 0.001)") +
  geom_smooth(method='lm',se = FALSE, color="black") +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp
ggsave("F1Prediction_avg_gpd_NoCor.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)

#dev.off()

# With labels, not used in manuscript
#png(file="F1Prediction_avg_withlabels_gpdNoCor.png",width=20,height=14,units="cm",res=300)

sp<-ggplot(mergedF1info, aes(x= Means, y=Avg.Dryweight.F1)) + geom_point(alpha=0.8,size=5,color=wes_palette("Rushmore1")[4]) +
  labs(x="Avg. parental GEBV", y = "Avg. dry weight of population") +
  xlim(c(min(mergedF1info$Means),max(mergedF1info$Means)))+
  geom_text(x=0.01, y=3, label="R",fontface='italic') +
  geom_text(x=0.015, y=3, label="= 0.91 (p < 0.001)") +
  geom_label_repel(aes(label = mergedF1info$Population),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  geom_smooth(method='lm',se = FALSE, color="black") +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp
ggsave("F1Prediction_avg_withlabels_gpdNoCor.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)

#dev.off()

# Figure 6 in manuscript
#png(file="F1Prediction_avg_withlabels_gpdNoCor.png",width=20,height=14,units="cm",res=300)

sp<-ggplot(mergedF1info, aes(x= Means, y=Avg.Dryweight.F1,colour=(Population))) + geom_point(alpha=1.0,size=5) +
  scale_color_manual(values=colors) +
  labs(x="Avg. parental GEBV", y = "Avg. dry weight of population") +
  xlim(c(min(mergedF1info$Means),max(mergedF1info$Means)))+
  #geom_text(x=0.01, y=3, label="R",fontface='italic') +
  #geom_text(x=0.015, y=3, label="= 0.91 (p < 0.001)") +
  geom_smooth(method='lm',se = FALSE, color="black") +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp

ggsave("F1Prediction_avg_colByPop_gpdNoCor.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)

#dev.off()


############ Make F1 prediction plot, gpd_Fixcor
F1predictions=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/F1_prediction_20200828/gpd_fixcor/AvgParentalGEBV.txt",head=T)
head(F1predictions)
dim(F1predictions)
F1predictions_t=t(F1predictions)
F1predictions_t=as.data.frame(F1predictions_t)
F1predictions_t_=F1predictions_t[2:nrow(F1predictions_t),]
F1predictions_t_$Population=as.character(row.names(F1predictions_t_))
rownames(F1predictions_t_)=NULL
colnames(F1predictions_t_)=c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","Population")
head(F1predictions_t_)


F1observations=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/F1_prediction/With Lflex/F1observations.txt",head=T)
head(F1observations)

mergedF1info=merge(F1predictions_t_,F1observations,by="Population")
mergedF1info_onlyPred=mergedF1info[,1:(ncol(mergedF1info)-2)]
mergedF1info$Means=rowMeans(mergedF1info_onlyPred[sapply(mergedF1info_onlyPred, is.numeric)])


correlation(as.numeric(as.character(mergedF1info$Means)),as.numeric(as.character(mergedF1info$Avg.Dryweight.F1))) #0.44 (N.S)

mergedF1info
mergedF1info$SD=apply(mergedF1info_onlyPred[,-1],1,sd)

# Not used in manuscript
#png(file="F1Prediction_avg_gpd_FixCor.png",width=20,height=14,units="cm",res=300)

sp<-ggplot(mergedF1info, aes(x= Means, y=Avg.Dryweight.F1)) + geom_point(alpha=0.8,size=5,color=wes_palette("Rushmore1")[2]) +
  labs(x="Avg. parental GEBV", y = "Avg. dry weight of population") +
  xlim(c(min(mergedF1info$Means),max(mergedF1info$Means)))+
  geom_text(x=0.0010, y=3, label="R",fontface='italic') +
  geom_text(x=0.0015, y=3, label="= 0.44 (p > 0.05)") +
  geom_smooth(method='lm',se = FALSE, color="black") +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp

#dev.off()
ggsave("F1Prediction_avg_gpdFixCor.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)


# Not used in manuscript

sp<-ggplot(mergedF1info, aes(x= Means, y=Avg.Dryweight.F1)) + geom_point(alpha=0.8,size=5,color=wes_palette("Rushmore1")[2]) +
  labs(x="Avg. parental GEBV", y = "Avg. dry weight of population") +
  xlim(c(min(mergedF1info$Means),max(mergedF1info$Means)))+
  geom_text(x=0.0010, y=3, label="R",fontface='italic') +
  geom_text(x=0.0015, y=3, label="= 0.44 (p > 0.05)") +
  geom_smooth(method='lm',se = FALSE, color="black") +
  geom_label_repel(aes(label = mergedF1info$Population),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp
ggsave("F1Prediction_avg_withlabels_gpdFixCor.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)
#dev.off()


# Supplementary file
sp<-ggplot(mergedF1info, aes(x= Means, y=Avg.Dryweight.F1,colour=(Population))) + geom_point(alpha=1.0,size=5) +
  scale_color_manual(values=colors) +
  labs(x="Avg. parental GEBV", y = "Avg. dry weight of population") +
  xlim(c(min(mergedF1info$Means),max(mergedF1info$Means)))+
  #geom_text(x=0.01, y=3, label="R",fontface='italic') +
  #geom_text(x=0.015, y=3, label="= 0.91 (p < 0.001)") +
  geom_smooth(method='lm',se = FALSE, color="black") +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp

ggsave("F1Prediction_avg_colByPop_gpdFixCor.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)






# F1 prediction from average phenotypes, not genomic prediction

# gpd_NoCor_phenotypicSelection 
d6table=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/JustTheWhole_d6Table.txt",head=T)
head(d6table)

gpd_phenotypes=d6table[,c(27,29,13)]

F1pop1=gpd_phenotypes[which(gpd_phenotypes$Clovershort %in% parentpop1names),]
F1pop2=gpd_phenotypes[which(gpd_phenotypes$Clovershort %in% parentpop6names),]
F1pop3=gpd_phenotypes[which(gpd_phenotypes$Clovershort %in% parentpop7names),]
F1pop4=gpd_phenotypes[which(gpd_phenotypes$Clovershort %in% parentpop8names),]
F1pop5=gpd_phenotypes[which(gpd_phenotypes$Clovershort %in% parentpop9names),]
F1pop6=gpd_phenotypes[which(gpd_phenotypes$Clovershort %in% parentpop10names),]
F1pop7=gpd_phenotypes[which(gpd_phenotypes$Clovershort %in% parentpop11names),]
F1pop8=gpd_phenotypes[which(gpd_phenotypes$Clovershort %in% parentpop12names),]
F1pop9=gpd_phenotypes[which(gpd_phenotypes$Clovershort %in% parentpop13names),]

Parental_means_gpdNoCor1=aggregate(as.numeric(F1pop1$growth_per_day), list(F1pop1$Clovershort), mean)
F1pop1_parental_avg_GEBV=mean(Parental_means_gpdNoCor1[,2])
Parental_means_gpdNoCor2=aggregate(as.numeric(F1pop2$growth_per_day), list(F1pop2$Clovershort), mean)
F1pop2_parental_avg_GEBV=mean(Parental_means_gpdNoCor2[,2])
Parental_means_gpdNoCor3=aggregate(as.numeric(F1pop3$growth_per_day), list(F1pop3$Clovershort), mean)
F1pop3_parental_avg_GEBV=mean(Parental_means_gpdNoCor3[,2])
Parental_means_gpdNoCor4=aggregate(as.numeric(F1pop4$growth_per_day), list(F1pop4$Clovershort), mean)
F1pop4_parental_avg_GEBV=mean(Parental_means_gpdNoCor4[,2])
Parental_means_gpdNoCor5=aggregate(as.numeric(F1pop5$growth_per_day), list(F1pop5$Clovershort), mean)
F1pop5_parental_avg_GEBV=mean(Parental_means_gpdNoCor5[,2])
Parental_means_gpdNoCor6=aggregate(as.numeric(F1pop6$growth_per_day), list(F1pop6$Clovershort), mean)
F1pop6_parental_avg_GEBV=mean(Parental_means_gpdNoCor6[,2])
Parental_means_gpdNoCor7=aggregate(as.numeric(F1pop7$growth_per_day), list(F1pop7$Clovershort), mean)
F1pop7_parental_avg_GEBV=mean(Parental_means_gpdNoCor7[,2])
Parental_means_gpdNoCor8=aggregate(as.numeric(F1pop8$growth_per_day), list(F1pop8$Clovershort), mean)
F1pop8_parental_avg_GEBV=mean(Parental_means_gpdNoCor8[,2])
Parental_means_gpdNoCor9=aggregate(as.numeric(F1pop9$growth_per_day), list(F1pop9$Clovershort), mean)
F1pop9_parental_avg_GEBV=mean(Parental_means_gpdNoCor9[,2])

F1_gpdNoCor_pheno=cbind(mergedF1info$Population,c(F1pop1_parental_avg_GEBV,F1pop5_parental_avg_GEBV,F1pop6_parental_avg_GEBV,F1pop7_parental_avg_GEBV,F1pop8_parental_avg_GEBV,F1pop2_parental_avg_GEBV,F1pop3_parental_avg_GEBV,F1pop4_parental_avg_GEBV,F1pop9_parental_avg_GEBV))
F1_gpdNoCor_pheno=as.data.frame(F1_gpdNoCor_pheno)
colnames(F1_gpdNoCor_pheno)=c("Population","Pheno_gpd_NoCor")
F1_gpdNoCor_pheno=cbind(F1_gpdNoCor_pheno,mergedF1info[,13])
colnames(F1_gpdNoCor_pheno)[3]="dryweight"

correlation(as.numeric(as.character(F1_gpdNoCor_pheno$Pheno_gpd_NoCor)),as.numeric(as.character(F1_gpdNoCor_pheno$dryweight))) #0.92, ***

# figure 6 in manuscript
#png(file="F1_gpdNocor_avg_withlabels.png",width=20,height=14,units="cm",res=300)

sp<-ggplot(F1_gpdNoCor_pheno, aes(x= as.numeric(as.character(Pheno_gpd_NoCor)), y=dryweight)) + geom_point(alpha=0.8,size=5,color=wes_palette("Rushmore1")[4]) +
  labs(x="Avg. parental gpd_NoCor", y = "Avg. dry weight of population") +
  geom_smooth(method='lm',se = FALSE, color="black") +
  geom_text(x=0.56, y=3, label="R",fontface='italic') +
  geom_text(x=0.60, y=3, label="= 0.92 (p < 0.001)") +
  geom_label_repel(aes(label = F1_gpdNoCor_pheno$Population),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp
ggsave("F1_gpdNocor_avg_withlabels.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)
#dev.off()

# figure 6 in manuscript
#png(file="F1_gpdNocor_avg_withlabels.png",width=20,height=14,units="cm",res=300)

sp<-ggplot(F1_gpdNoCor_pheno, aes(x= as.numeric(as.character(Pheno_gpd_NoCor)), y=dryweight,colour=(Population))) + geom_point(alpha=1.0,size=5) +
  scale_color_manual(values=colors) +
  labs(x="Avg. parental gpd_NoCor", y = "Avg. dry weight of population") +
  geom_smooth(method='lm',se = FALSE, color="black") +
  #geom_text(x=0.56, y=3, label="R",fontface='italic') +
  #geom_text(x=0.60, y=3, label="= 0.92 (p < 0.001)") +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp

ggsave("F1_gpdNocor_avg_ColByGroup.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)
#dev.off()

# gpd_ResCor_phenotypicSelection 
Parental_means_gpdResCor1=aggregate(as.numeric(F1pop1$gpd_dryweight_cor), list(F1pop1$Clovershort), mean)
F1pop1_parental_avg_GEBV=mean(Parental_means_gpdResCor1[,2])
Parental_means_gpdResCor2=aggregate(as.numeric(F1pop2$gpd_dryweight_cor), list(F1pop2$Clovershort), mean)
F1pop2_parental_avg_GEBV=mean(Parental_means_gpdResCor2[,2])
Parental_means_gpdResCor3=aggregate(as.numeric(F1pop3$gpd_dryweight_cor), list(F1pop3$Clovershort), mean)
F1pop3_parental_avg_GEBV=mean(Parental_means_gpdResCor3[,2])
Parental_means_gpdResCor4=aggregate(as.numeric(F1pop4$gpd_dryweight_cor), list(F1pop4$Clovershort), mean)
F1pop4_parental_avg_GEBV=mean(Parental_means_gpdResCor4[,2])
Parental_means_gpdResCor5=aggregate(as.numeric(F1pop5$gpd_dryweight_cor), list(F1pop5$Clovershort), mean)
F1pop5_parental_avg_GEBV=mean(Parental_means_gpdResCor5[,2])
Parental_means_gpdResCor6=aggregate(as.numeric(F1pop6$gpd_dryweight_cor), list(F1pop6$Clovershort), mean)
F1pop6_parental_avg_GEBV=mean(Parental_means_gpdResCor6[,2])
Parental_means_gpdResCor7=aggregate(as.numeric(F1pop7$gpd_dryweight_cor), list(F1pop7$Clovershort), mean)
F1pop7_parental_avg_GEBV=mean(Parental_means_gpdResCor7[,2])
Parental_means_gpdResCor8=aggregate(as.numeric(F1pop8$gpd_dryweight_cor), list(F1pop8$Clovershort), mean)
F1pop8_parental_avg_GEBV=mean(Parental_means_gpdResCor8[,2])
Parental_means_gpdResCor9=aggregate(as.numeric(F1pop9$gpd_dryweight_cor), list(F1pop9$Clovershort), mean)
F1pop9_parental_avg_GEBV=mean(Parental_means_gpdResCor9[,2])



F1_gpdResCor_pheno=cbind(mergedF1info$Population,c(F1pop1_parental_avg_GEBV,F1pop5_parental_avg_GEBV,F1pop6_parental_avg_GEBV,F1pop7_parental_avg_GEBV,F1pop8_parental_avg_GEBV,F1pop2_parental_avg_GEBV,F1pop3_parental_avg_GEBV,F1pop4_parental_avg_GEBV,F1pop9_parental_avg_GEBV))
F1_gpdResCor_pheno=as.data.frame(F1_gpdResCor_pheno)
colnames(F1_gpdResCor_pheno)=c("Population","Pheno_gpd_ResCor")
F1_gpdResCor_pheno=cbind(F1_gpdResCor_pheno,mergedF1info[,13])
colnames(F1_gpdResCor_pheno)[3]="dryweight"

correlation(as.numeric(as.character(F1_gpdResCor_pheno$Pheno_gpd_ResCor)),as.numeric(as.character(F1_gpdResCor_pheno$dryweight))) #0.92, ***

# figure 6 in manuscript
#png(file="F1_gpdRescor_avg_withlabels.png",width=20,height=14,units="cm",res=300)

sp<-ggplot(F1_gpdResCor_pheno, aes(x= as.numeric(as.character(Pheno_gpd_ResCor)), y=dryweight)) + geom_point(alpha=0.8,size=5,color=wes_palette("Rushmore1")[3]) +
  labs(x="Avg. parental gpd_ResCor", y = "Avg. dry weight of population") +
  geom_text(x=0.13, y=3, label="R",fontface='italic') +
  geom_text(x=0.18, y=3, label="= 0.92 (p < 0.001)") +
  geom_smooth(method='lm',se = FALSE, color="black") +
  geom_label_repel(aes(label = F1_gpdResCor_pheno$Population),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp
ggsave("F1_gpdRescor_avg_withlabels.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)
#dev.off()

sp<-ggplot(F1_gpdResCor_pheno, aes(x= as.numeric(as.character(Pheno_gpd_ResCor)), y=dryweight,colour=(Population))) + geom_point(alpha=1.0,size=5) +
  scale_color_manual(values=colors) +
  labs(x="Avg. parental gpd_ResCor", y = "Avg. dry weight of population") +
  geom_smooth(method='lm',se = FALSE, color="black") +
  #geom_text(x=0.56, y=3, label="R",fontface='italic') +
  #geom_text(x=0.60, y=3, label="= 0.92 (p < 0.001)") +
  theme_classic() +theme(axis.text.x = element_text(size = 6, angle = 90,hjust=1),axis.text.y = element_text(size = 18), 
                         axis.title=element_text(size=22))
sp

ggsave("F1_gpdRescor_avg_ColByGroup.pdf", width =20, height = 14, units = "cm",useDingbats=FALSE)

# Hoteling-Williams tests to see if correlations differ
correlation(as.numeric(as.character(F1_gpdResCor_pheno$Pheno_gpd_ResCor)),as.numeric(as.character(F1_gpdResCor_pheno$dryweight))) #0.92, ***
correlation(as.numeric(as.character(F1_gpdNoCor_pheno$Pheno_gpd_NoCor)),as.numeric(as.character(F1_gpdNoCor_pheno$dryweight))) #0.92, ***
correlation(as.numeric(as.character(F1_gpdNoCor_pheno$Pheno_gpd_NoCor)),as.numeric(as.character(F1_gpdResCor_pheno$Pheno_gpd_ResCor))) #1, **'


r.test(n=9, r12=.92, r13=.92,r23=0.99) # significant (*) # That is phenotypic selection from gpd_NoCor is significantly better than from gpd_rescor

# Is genomic selection for gpd_Rescor worse than phenotypic selection?
F1predictions=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/F1_prediction_20200828/gpd_rescor/AvgParentalGEBV.txt",head=T)
head(F1predictions)
dim(F1predictions)
F1predictions_t=t(F1predictions)
F1predictions_t=as.data.frame(F1predictions_t)
F1predictions_t_=F1predictions_t[2:nrow(F1predictions_t),]
F1predictions_t_$Population=as.character(row.names(F1predictions_t_))
rownames(F1predictions_t_)=NULL
colnames(F1predictions_t_)=c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","Population")
head(F1predictions_t_)


F1observations=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/F1_prediction/With Lflex/F1observations.txt",head=T)
head(F1observations)

mergedF1info=merge(F1predictions_t_,F1observations,by="Population")
mergedF1info_onlyPred=mergedF1info[,1:(ncol(mergedF1info)-2)]
mergedF1info$Means=rowMeans(mergedF1info_onlyPred[sapply(mergedF1info_onlyPred, is.numeric)])


correlation(as.numeric(as.character(mergedF1info$Means)),as.numeric(as.character(mergedF1info$Avg.Dryweight.F1))) #0.93 ***
correlation(as.numeric(as.character(F1_gpdResCor_pheno$Pheno_gpd_ResCor)),as.numeric(as.character(F1_gpdResCor_pheno$dryweight))) #0.92, ***
correlation(as.numeric(as.character(mergedF1info$Means)),as.numeric(as.character(F1_gpdResCor_pheno$Pheno_gpd_ResCor))) #0.96 **

r.test(n=9, r12=.93, r13=.92,r23=.96) # significant (*) # That is phenotypic selection from gpd_NoCor is significantly better than from gpd_rescor


# is gp for ResCor (0.93) better than for NoCor (0.91)
correlation(as.numeric(as.character(mergedF1info_No$Means)),as.numeric(as.character(mergedF1info_No$Avg.Dryweight.F1))) #0.93 ***
correlation(as.numeric(as.character(mergedF1info$Means)),as.numeric(as.character(mergedF1info$Avg.Dryweight.F1))) #0.93 ***
correlation(as.numeric(as.character(mergedF1info_No$Means)),as.numeric(as.character(mergedF1info$Means)))

r.test(n=9, r12=.91, r13=.93,r23=.99) # N.S.

# is gp for ResCor (0.93) better than phenotypic selection (0.92)
correlation(as.numeric(as.character(mergedF1info$Means)),as.numeric(as.character(F1_gpdResCor_pheno$Pheno_gpd_ResCor))) #0.96 ***

r.test(n=9, r12=.92, r13=.93,r23=.96) # N.S.

# is gp for NoCor (0.91) worse than phenotypic selection (0.92)
correlation(as.numeric(as.character(mergedF1info_No$Means)),as.numeric(as.character(F1_gpdNoCor_pheno$Pheno_gpd_NoCor))) #0.96 ***

r.test(n=9, r12=.92, r13=.91,r23=.96) # N.S.


#### make correlation plot, supplementary figure?
mean_merged_rescor=subset(mean_merged,gpd_measurement=="gpd_ResCor")
mean_merged_NoCor=subset(mean_merged,gpd_measurement=="gpd_NoCor")
mean_merged_FixCor=subset(mean_merged,gpd_measurement=="gpd_FixCor")

Correlations=cbind(as.character(mean_merged_rescor$Individual),mean_merged_rescor$Observed_correctedForFixed,mean_merged_rescor$GEBVs,mean_merged_NoCor$Observed_correctedForFixed,mean_merged_NoCor$GEBVs,mean_merged_FixCor$Observed_correctedForFixed,mean_merged_FixCor$GEBVs,na.omit(Phenotypic_averages[,2]))
colnames(Correlations)=c("Genotype","gpd_ResCor_Correctedpheno","gpd_ResCor_GEBV","gpd_NoCor_Correctedpheno","gpd_NoCor_GEBV","gpd_FixCor_Correctedpheno","gpd_FixCor_GEBV","Phenotypic average")
Correlations_numeric=Correlations[,2:ncol(Correlations)]
Correlations_numeric=as.data.frame(Correlations_numeric)

write.table(Correlations,"Corrgramdata_20200924.txt",sep="\t",col.names=T, row.names=F, quote=F)

Correlations_numeric[,1]=as.numeric(as.character(Correlations_numeric[,1]))
Correlations_numeric[,2]=as.numeric(as.character(Correlations_numeric[,2]))
Correlations_numeric[,3]=as.numeric(as.character(Correlations_numeric[,3]))
Correlations_numeric[,4]=as.numeric(as.character(Correlations_numeric[,4]))
Correlations_numeric[,5]=as.numeric(as.character(Correlations_numeric[,5]))
Correlations_numeric[,6]=as.numeric(as.character(Correlations_numeric[,6]))
Correlations_numeric[,7]=as.numeric(as.character(Correlations_numeric[,7]))
Correlations_numeric_matrix=data.matrix(Correlations_numeric)

#png(file="correlations.png",width=24,height=24,units="cm",res=300)

corrgram(Correlations_numeric_matrix,
         main="Correlations between different gpd corrections and their prediction",
         lower.panel=panel.pts, upper.panel=panel.conf,
         diag.panel=panel.density)

#dev.off()
ggsave("correlations.pdf", width =24, height = 24, units = "cm")







# trait investigation
# should not be in this graphing script but the d6 table was here, sooo..
# But should be moved 

cor(d6table$InitialSize,d6table$gpd_dryweight_cor) # 0.18
cor(d6table$InitialSize,d6table$growth_per_day) # 0.72

correlation(d6table$InitialSize,as.numeric(d6table$inoculation_date)) # -0.35
correlation(d6table$growth_per_day,as.numeric(d6table$inoculation_date)) # -0.51 
correlation(d6table$gpd_dryweight_cor,as.numeric(d6table$inoculation_date)) # -0.36

fit1=lmer(InitialSize  ~ (1|Clover) +as.factor(EW) +as.factor(NS) +factor(Round) +factor(Rhizobium) +factor(inoculation_date), data=d6table)
summary(fit1) #0.11 broad sense her

fit1_alt1=lmer(InitialSize  ~ (1|Clover) +as.factor(EW) +as.factor(NS) +factor(Round) +factor(Rhizobium), data=d6table)
summary(fit1_alt1) #0.11 broad sense her

fit1_alt2=lmer(InitialSize  ~ (1|Clover) +as.factor(EW) +as.factor(NS) +factor(Round), data=d6table)
summary(fit1_alt2) #0.15 broad sense her

fit2=lmer(growth_per_day ~ (1|Clover) +as.factor(EW) +as.factor(NS) +factor(Round) +factor(Rhizobium) +factor(inoculation_date), data=d6table)
summary(fit2)#0.29 broad sense her

fit3=lmer(gpd_dryweight_cor ~ (1|Clover) +as.factor(EW) +as.factor(NS) +factor(Round) +factor(Rhizobium) +factor(inoculation_date), data=d6table) 
summary(fit3)#0.45 broad sense her






# Load predictions from iSize
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_iSize_20200925/")
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
meaniSize=aggregate(as.numeric(dataset4$Observed_correctedForFixed), list(dataset4$Individual), mean)[,2]
meangebv=aggregate(as.numeric(dataset4$GEBVs), list(dataset4$Individual), mean)
df_iSize=cbind(meangebv,meaniSize)
colnames(df_iSize)=c("Individual","GEBVs","Observed_correctedForFixed")
head(df_iSize)

Correlations_numeric_matrix=as.data.frame(Correlations_numeric_matrix)
Correlations_numeric_matrix$iSize_GEBV=as.numeric(as.character(meangebv[,2]))
Correlations_numeric_matrix$iSize_pheno=as.numeric(as.character(meaniSize))

corrgram(Correlations_numeric_matrix,
         main="Correlations between different gpd corrections and their prediction",
         lower.panel=panel.pts, upper.panel=panel.conf,
         diag.panel=panel.density)

ggsave("correlations_withiSize.pdf", width =24, height = 24, units = "cm")




# Check how genotypes rank when looking at their average gpd_ResCor
Fulldata=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/JustTheWhole_d6Table.txt",head=T,sep="\t")
dim(Fulldata)

Mean_gpdResCor=aggregate(as.numeric(Fulldata$gpd_dryweight_cor), list(Fulldata$Clovershort), mean)
colnames(Mean_gpdResCor)=c("Individual","gpd_ResCor")
#order from lowest to highest
Mean_gpdResCor_ordered <- Mean_gpdResCor[order(Mean_gpdResCor$gpd_ResCor),]
individual_order_10=Mean_gpdResCor_ordered$Individual

keyDF <- data.frame(key=individual_order_10,weight=1:length(individual_order_10))
Mean_gpdResCor_ordered_ready <- merge(Mean_gpdResCor_ordered,keyDF,by.x='Individual',by.y='key',all.x=T,all.y=F)
Mean_gpdResCor_ordered_ready_ready=Mean_gpdResCor_ordered_ready[order(Mean_gpdResCor_ordered_ready$weight),]

Mean_gpdResCor_ordered_ready_ready$Individual <- factor(Mean_gpdResCor_ordered_ready_ready$Individual , levels = Mean_gpdResCor_ordered_ready_ready$Individual)


# Not used in manuscript
png(file="Gpd_resCor_versusIndividuals.png",width=70,height=60,units="cm",res=300)

sp<-ggplot(Mean_gpdResCor_ordered_ready_ready, aes(x= Individual, y=gpd_ResCor)) + geom_point(alpha=0.7,size=5,color=wes_palette("Rushmore1")[3])+
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp

dev.off()



# not used in manuscript
# sorted after avg GEBVs of ResCor and individual parents marked 
png(file="AllParentsGPDRes.png",width=210,height=180,units="cm",res=300)

p1=ggplot(Mean_gpdResCor_ordered_ready_ready, aes(x= Individual, y=gpd_ResCor)) + geom_point(alpha=0.8,size=18,color=wes_palette("Rushmore1")[3])+
  gghighlight(Individual %in% parentpop1names) +ggtitle("DLF1") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p2=ggplot(Mean_gpdResCor_ordered_ready_ready, aes(x= Individual, y=gpd_ResCor)) + geom_point(alpha=0.8,size=18,color=wes_palette("Rushmore1")[3])+
  gghighlight(Individual %in% parentpop6names) +ggtitle("LJ L6") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p3=ggplot(Mean_gpdResCor_ordered_ready_ready, aes(x= Individual, y=gpd_ResCor)) + geom_point(alpha=0.8,size=18,color=wes_palette("Rushmore1")[3])+
  gghighlight(Individual %in% parentpop7names) +ggtitle("LJ L7") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p4=ggplot(Mean_gpdResCor_ordered_ready_ready, aes(x= Individual, y=gpd_ResCor)) + geom_point(alpha=0.8,size=18,color=wes_palette("Rushmore1")[3])+
  gghighlight(Individual %in% parentpop8names) +ggtitle("LJ L9") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p5=ggplot(Mean_gpdResCor_ordered_ready_ready, aes(x= Individual, y=gpd_ResCor)) + geom_point(alpha=0.8,size=18,color=wes_palette("Rushmore1")[3])+
  gghighlight(Individual %in% parentpop9names) +ggtitle("LJ H1") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p6=ggplot(Mean_gpdResCor_ordered_ready_ready, aes(x= Individual, y=gpd_ResCor)) + geom_point(alpha=0.8,size=18,color=wes_palette("Rushmore1")[3])+
  gghighlight(Individual %in% parentpop10names) +ggtitle("LJ H2") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p7=ggplot(Mean_gpdResCor_ordered_ready_ready, aes(x= Individual, y=gpd_ResCor)) + geom_point(alpha=0.8,size=18,color=wes_palette("Rushmore1")[3])+
  gghighlight(Individual %in% parentpop11names) +ggtitle("LJ H3") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p8=ggplot(Mean_gpdResCor_ordered_ready_ready, aes(x= Individual, y=gpd_ResCor)) + geom_point(alpha=0.8,size=18,color=wes_palette("Rushmore1")[3])+
  gghighlight(Individual %in% parentpop12names) +ggtitle("LJ H5") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

p9=ggplot(Mean_gpdResCor_ordered_ready_ready, aes(x= Individual, y=gpd_ResCor)) + geom_point(alpha=0.8,size=18,color=wes_palette("Rushmore1")[3])+
  gghighlight(Individual %in% parentpop13names) +ggtitle("SUA") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)
dev.off()

