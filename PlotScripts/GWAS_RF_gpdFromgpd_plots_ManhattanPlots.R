#########################################
#########################################
##  This is a script to display the 
## importance of given SNPs for
## iSize and gpd
#########################################
#########################################

library(ggplot2)
library(dplyr)
library(viridis)


# GWAS gpd results
setwd("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/pvals_GWAS_gpd")




# upload all gwas results from 10 rounds of 6 fold CV on gpd (real files to be updated)
{
list.files()

file_list <- list.files(pattern="_pid1_gpd")
length(file_list)==60 #check
print("p-values from this many")
print(length(file_list))

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset_pvals_gpd")){
    dataset_pvals_gpd <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset_pvals_gpd")){
    pred_dataset <-read.table(file, header=TRUE, sep=",")
    dataset_pvals_gpd<-rbind(dataset_pvals_gpd, pred_dataset)
    rm(pred_dataset)
  }
}

dataset_pvals_gpd=as.data.frame(dataset_pvals_gpd)
head(dataset_pvals_gpd)
nrow(dataset_pvals_gpd)%%600==0 # check if divisable by 600. It should be
}


# Calculate average p-values
{
dataset_pvals_gpd$identifier = paste(dataset_pvals_gpd$chromosomes,dataset_pvals_gpd$positions, sep='-')
head(dataset_pvals_gpd)
Pvals=aggregate(dataset_pvals_gpd$scores, list(dataset_pvals_gpd$identifier), mean)
chromosomes=aggregate(dataset_pvals_gpd$chromosomes, list(dataset_pvals_gpd$identifier), mean)
positions=aggregate(dataset_pvals_gpd$positions, list(dataset_pvals_gpd$identifier), mean)
mafs=aggregate(dataset_pvals_gpd$mafs, list(dataset_pvals_gpd$identifier), mean)

fulldata=cbind(chromosomes[,2],positions[,2],Pvals[,2],mafs[,2])
fulldata=as.data.frame(fulldata)
colnames(fulldata)=c("chromosomes","positions","scores","mafs")
head(fulldata)

fulldata=fulldata[-which(fulldata$mafs<0.05),]
}


# Plotting
{
  #prepare for plotting
  datafiltered <- subset(fulldata, chromosomes >0)
  datafiltered$p.bh <- p.adjust(datafiltered$scores, "BH")
  p.bonferroni <- 1/nrow(datafiltered)*0.05
  
  
  don <- datafiltered %>% 
    
    # Compute chromosome size
    group_by(chromosomes) %>% 
    summarise(chr_len=max(positions)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(datafiltered, ., by=c("chromosomes"="chromosomes")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chromosomes, positions) %>%
    mutate( BPcum=positions+tot) %>%
    
    mutate( is_highlight=ifelse(p.bh<0.05, "yes", "no")) 
  
  
  # Prepare X axis 
  axisdf = don %>% group_by(chromosomes) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  
  
  # plotting
  p <- ggplot(don, aes(x=BPcum, y=-log10(scores))) +
    
    # Show all points
    geom_point( aes(color=as.factor(chromosomes)), alpha=0.8, size=1.0) +
    scale_color_manual(values = rep(c("gray33", "gray67"), 16 )) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$chromosomes, breaks= axisdf$center, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,max(-log10(don$scores))+0.3), expand = expand_scale(mult = c(0, .1))) +     # remove space between plot area and x axis
    # Add highlighted points
    geom_point(data=subset(don, is_highlight=="yes"), color="darkgoldenrod1", size=1.6) +
    
    # Add thresholds
    geom_hline(yintercept=-log10(p.bonferroni), linetype="dashed", color = "#f04546") + 
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      #panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 22),
      title = element_text(size = 20)
    ) +
    labs(y = ("-log10(p-value)"), x = ("Chromosome"), title = ("gpd"))
  
  p
}


# Upload top 200 SNPs
{
  setwd("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/top200gpd")
  
  list.files()
  
  file_list <- list.files(pattern="Round")
  
  for (file in file_list){
    
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset_top200")){
      dataset_top200 <- read.table(file, header=FALSE, sep=",")
    }
    
    # if the merged dataset does exist, append to it
    else{
      pred_dataset <-read.table(file, header=FALSE, sep=",")
      dataset_top200<-rbind(dataset_top200, pred_dataset)
      rm(pred_dataset)
    }
  }
  
  nrow(dataset_top200)==200*60
  head(dataset_top200)
  colnames(dataset_top200)[1:3]=c("chromosomes","positions","scores")
}

# Calculate number of occurences for each snp in top 200 
{
  dataset_top200$identifier = paste(dataset_top200$chromosomes,dataset_top200$positions, sep='-')

  # count occurence
  count_occurence = dataset_top200 %>% group_by(identifier) %>% mutate(count = n())
  write.table(count_occurence,"occurence_top200SNPs_gpd.txt")
  occurence_fraction = round(count_occurence$count/length(file_list),2)
  
  occurence_table= cbind(count_occurence$identifier,occurence_fraction)
  occurence_table=as.data.frame(occurence_table)
  occurence_table_nodup=occurence_table %>% distinct()
  dim(occurence_table_nodup)
  colnames(occurence_table_nodup)=c("identifier","Occurence")
  occurence_table_nodup$identifier=as.character(occurence_table_nodup$identifier)
  
  # now merge full dataframe with occurences
  tempfulldata= don
  tempfulldata$identifier = paste(tempfulldata$chromosomes,tempfulldata$positions, sep='-')
  
  merged_with_occurence = merge(tempfulldata,occurence_table_nodup,by="identifier",all.x=T)
  head(merged_with_occurence)
  merged_with_occurence$Occurence=as.numeric(as.character(merged_with_occurence$Occurence))
  #merged_with_occurence$Occurence[is.na(merged_with_occurence$Occurence)]=0
}

# plot an extra layer that assigns dot size after occurence
{
  
  extralayer1 = merged_with_occurence
  str(extralayer1)
  
  extralayer_naomit= extralayer1[-which(is.na(extralayer1$Occurence)),]
  
  p1 = p  +  geom_point(data=extralayer1,aes(x=BPcum, y=-log10(scores),size=Occurence,color=as.factor(chromosomes))) +  scale_size(range = c(1.2,6)) +     scale_color_manual(values = rep(c("gray33", "gray67"), 16 ))

  p1
  ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/Manhattanplot_gpdfromgpd',Sys.Date(),'.pdf',sep="_"), plot = p1, width = 25, height = 15, unit = 'cm')
}

# plot an extra layer that assigns dot colour after occurence
{
  
  # plotting
  p_nocolor <- ggplot(don, aes(x=BPcum, y=-log10(scores))) +
    
    # Show all points
    geom_point( aes(),color="grey20",size=2) +

    # custom X axis:
    scale_x_continuous(label = axisdf$chromosomes, breaks= axisdf$center, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,max(-log10(don$scores))+0.3), expand = expand_scale(mult = c(0, .1))) +     # remove space between plot area and x axis

    # Add thresholds
    geom_hline(yintercept=-log10(p.bonferroni), linetype="dashed", color = "#f04546") + 
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      #panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 22),
      title = element_text(size = 20)
    ) +
    labs(y = ("-log10(p-value)"), x = ("Chromosome"), title = ("gpd"))
  
  p_nocolor
  
  p_coloured <- p_nocolor + geom_point(data=extralayer_naomit, aes(x=BPcum, y=-log10(scores), color=Occurence,size=0.1)) + scale_color_viridis() 
  p_coloured
  
  ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/gpdfromgpd_Manhattan_ColourafterOccurencetop200',Sys.Date(),'.pdf',sep="_"), plot = p_coloured, width = 25, height = 15, unit = 'cm')
}


# plot an extra layer that assigns dot colour after occurence, do only color snps occuring at least in 10 % of the cases

extralayer_only90=extralayer_naomit[-which(extralayer_naomit$Occurence<0.10),]

p_coloured_only90 <- p_nocolor + geom_point(data=extralayer_only90, aes(x=BPcum, y=-log10(scores), color=Occurence,size=Occurence)) + scale_color_viridis() +  scale_size(range = c(1.5,6))
p_coloured_only90

ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/gpdfromgpd_Manhattan_ColourafterOccurenceOnly90perctop200',Sys.Date(),'.pdf',sep="_"), plot = p_coloured_only90, width = 35, height = 15, unit = 'cm')



# Load in weights when predicting gpd
{

setwd("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/weights_gpdfromgpd")
  
list.files()
file_list <- list.files(pattern="Weights_gpd_RF_top200SNPs")

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset_weight_top200")){
    dataset_weight_top200 <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset_weight_top200")){
    pred_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset_weight_top200<-rbind(dataset_weight_top200, pred_dataset)
    rm(pred_dataset)
  }
}

dataset_weight_top200=as.data.frame(dataset_weight_top200)
colnames(dataset_weight_top200)=c("Overall.importance","identifier")
head(dataset_weight_top200)
}

# average weight
{
  AvgWeight_top200=aggregate(dataset_weight_top200$Overall.importance, list(dataset_weight_top200$identifier), mean)
  head(AvgWeight_top200)
  colnames(AvgWeight_top200)=c("identifier","avgweight")
}

{
extralayer_only90_w_weightTop200= merge(extralayer_only90,AvgWeight_top200,by="identifier",all.x=T)
head(extralayer_only90_w_weightTop200)
}

p_colouredScaled_only90 <- p_nocolor + geom_point(data=extralayer_only90_w_weightTop200, aes(x=BPcum, y=-log10(scores), color=Occurence,size=avgweight)) + scale_color_viridis() +  scale_size(range = c(2,8))
p_colouredScaled_only90

ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/gpdfromisize_Manhattan_ColourafterOccurencegpdAfterWeightOnly90perc_top200',Sys.Date(),'.pdf',sep="_"), plot = p_colouredScaled_only90, width = 45, height = 15, unit = 'cm')




# then check if you can also colour black and grey according to chr 

# plotting
# plotting
p <- ggplot(don, aes(x=BPcum, y=-log10(scores))) +
  
  # Show all points
  geom_point( aes(), alpha=0.8, size=1.0) +
  #scale_color_manual(values = rep(c("gray33", "gray67"), 16 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$chromosomes, breaks= axisdf$center, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,max(-log10(don$scores))+0.3), expand = expand_scale(mult = c(0, .1))) +     # remove space between plot area and x axis
  # Add highlighted points
  
  geom_point(data=subset(don, chromosomes=="1"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="2"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="3"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="4"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="5"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="6"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="7"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="8"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="9"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="10"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="11"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="12"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="13"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="14"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="15"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="16"), color="gray67", size=1.6) +
  
  # Add thresholds
  geom_hline(yintercept=-log10(p.bonferroni), linetype="dashed", color = "#f04546") + 
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    #panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 22),
    title = element_text(size = 20)
  ) +
  labs(y = ("-log10(p-value)"), x = ("Chromosome"), title = ("gpd"))

p

p_colouredScaled_only90_ <- p  + geom_point(data=extralayer_only90_w_weightTop200, aes(x=BPcum, y=-log10(scores), color=Occurence,size=avgweight)) + scale_color_viridis(option = "plasma") +  scale_color_viridis(option = "plasma") + 
  scale_size_continuous(limits=c(0,100), breaks=c(0,10,20,30,40,50,60,70,80,90,100),labels=c("0","10","20","30","40","50","60","70","80","90","100"), name = "Weight",guide="legend",range = c(0, 8) )  + 
  theme(legend.position = "right") 
p_colouredScaled_only90_

ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/gpdfromgpd_Manhattan_ColourafterOccurenceSizeAfterWeightOnly90perc_backgroundcol_top200',Sys.Date(),'.pdf',sep="_"), plot = p_colouredScaled_only90_, width = 75, height = 15, unit = 'cm')
ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/gpdfromgpd_Manhattan_ColourafterOccurenceSizeAfterWeightOnly90perc_backgroundcol_top200',Sys.Date(),'.png',sep="_"), plot = p_colouredScaled_only90_, width = 75, height = 15, unit = 'cm')










# Upload top 25 SNPs
{
  setwd("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/top25gpd")
  
  list.files()
  
  file_list <- list.files(pattern="Round")
  
  for (file in file_list){
    
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset_top25")){
      dataset_top25 <- read.table(file, header=FALSE, sep=",")
    }
    
    # if the merged dataset does exist, append to it
    else{
      pred_dataset <-read.table(file, header=FALSE, sep=",")
      dataset_top25<-rbind(dataset_top25, pred_dataset)
      rm(pred_dataset)
    }
  }
  
  nrow(dataset_top25)==25*60
  head(dataset_top25)
  colnames(dataset_top25)[1:3]=c("chromosomes","positions","scores")
}

# Calculate number of occurences for each snp in top 200 
{
  dataset_top25$identifier = paste(dataset_top25$chromosomes,dataset_top25$positions, sep='-')
  
  # count occurence
  count_occurence = dataset_top25 %>% group_by(identifier) %>% mutate(count = n())
  write.table(count_occurence,"occurence_top25SNPs_gpd.txt")
 
   occurence_fraction = round(count_occurence$count/length(file_list),2)
  
  occurence_table= cbind(count_occurence$identifier,occurence_fraction)
  occurence_table=as.data.frame(occurence_table)
  occurence_table_nodup=occurence_table %>% distinct()
  dim(occurence_table_nodup)
  colnames(occurence_table_nodup)=c("identifier","Occurence")
  occurence_table_nodup$identifier=as.character(occurence_table_nodup$identifier)
  
  # now merge full dataframe with occurences
  tempfulldata= don
  tempfulldata$identifier = paste(tempfulldata$chromosomes,tempfulldata$positions, sep='-')
  
  merged_with_occurence = merge(tempfulldata,occurence_table_nodup,by="identifier",all.x=T)
  head(merged_with_occurence)
  merged_with_occurence$Occurence=as.numeric(as.character(merged_with_occurence$Occurence))
  #merged_with_occurence$Occurence[is.na(merged_with_occurence$Occurence)]=0
}

# plot an extra layer that assigns dot size after occurence
{
  
  extralayer1 = merged_with_occurence
  str(extralayer1)
  
  extralayer_naomit= extralayer1[-which(is.na(extralayer1$Occurence)),]
  
  p1 = p  +  geom_point(data=extralayer1,aes(x=BPcum, y=-log10(scores),size=Occurence,color=as.factor(chromosomes))) +  scale_size(range = c(1.2,6)) +     scale_color_manual(values = rep(c("gray33", "gray67"), 16 ))
  
  p1
  ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/Manhattanplot_gpdfromgpd_top25',Sys.Date(),'.pdf',sep="_"), plot = p1, width = 25, height = 15, unit = 'cm')
}

# plot an extra layer that assigns dot colour after occurence
{
  
  # plotting
  p_nocolor <- ggplot(don, aes(x=BPcum, y=-log10(scores))) +
    
    # Show all points
    geom_point( aes(),color="grey20",size=2) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$chromosomes, breaks= axisdf$center, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,max(-log10(don$scores))+0.3), expand = expand_scale(mult = c(0, .1))) +     # remove space between plot area and x axis
    
    # Add thresholds
    geom_hline(yintercept=-log10(p.bonferroni), linetype="dashed", color = "#f04546") + 
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      #panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 22),
      title = element_text(size = 20)
    ) +
    labs(y = ("-log10(p-value)"), x = ("Chromosome"), title = ("gpd"))
  
  p_nocolor
  
  p_coloured <- p_nocolor + geom_point(data=extralayer_naomit, aes(x=BPcum, y=-log10(scores), color=Occurence,size=0.1)) + scale_color_viridis() 
  p_coloured
  
  ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/gpdfromgpd_Manhattan_ColourafterOccurencetop25',Sys.Date(),'.pdf',sep="_"), plot = p_coloured, width = 25, height = 15, unit = 'cm')
}


# plot an extra layer that assigns dot colour after occurence, do only color snps occuring at least in 10 % of the cases

extralayer_only90=extralayer_naomit[-which(extralayer_naomit$Occurence<0.10),]

p_coloured_only90 <- p_nocolor + geom_point(data=extralayer_only90, aes(x=BPcum, y=-log10(scores), color=Occurence,size=Occurence)) + scale_color_viridis() +  scale_size(range = c(1.5,6))
p_coloured_only90

#ggsave(paste('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GWAS_iSizeByRF_20210201/iSizepredictioncorrelations/iSizefromiSize_Manhattan_ColourafterOccurenceOnly90perctop25',Sys.Date(),'.pdf',sep="_"), plot = p_coloured_only90, width = 35, height = 15, unit = 'cm')



# Load in weights when predicting gpd
{
  
  setwd("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/weights_gpdfromgpd")
  
  list.files()
  file_list <- list.files(pattern="Weights_gpd_RF_top25SNPs")
  
  for (file in file_list){
    
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset_weight_top25")){
      dataset_weight_top25 <- read.table(file, header=TRUE, sep="\t")
    }
    
    # if the merged dataset does exist, append to it
    else{
      pred_dataset <-read.table(file, header=TRUE, sep="\t")
      dataset_weight_top25<-rbind(dataset_weight_top25, pred_dataset)
      rm(pred_dataset)
    }
  }
  
  dataset_weight_top25=as.data.frame(dataset_weight_top25)
  colnames(dataset_weight_top25)=c("Overall.importance","identifier")
  head(dataset_weight_top25)
}

# average weight
{
  AvgWeight_top25=aggregate(dataset_weight_top25$Overall.importance, list(dataset_weight_top25$identifier), mean)
  head(AvgWeight_top25)
  colnames(AvgWeight_top25)=c("identifier","avgweight")
}

{
  extralayer_only90_w_weightTop25= merge(extralayer_only90,AvgWeight_top25,by="identifier",all.x=T)
  head(extralayer_only90_w_weightTop25)
}

p_colouredScaled_only90 <- p_nocolor + geom_point(data=extralayer_only90_w_weightTop25, aes(x=BPcum, y=-log10(scores), color=Occurence,size=avgweight)) + scale_color_viridis() +  scale_size(range = c(2,8))
p_colouredScaled_only90

ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/gpdfromgpd_Manhattan_ColourafterOccurenceSizeAfterWeightOnly90perc_top25',Sys.Date(),'.pdf',sep="_"), plot = p_colouredScaled_only90, width = 45, height = 15, unit = 'cm')




# then check if you can also colour black and grey according to chr 

# plotting
# plotting
p <- ggplot(don, aes(x=BPcum, y=-log10(scores))) +
  
  # Show all points
  geom_point( aes(), alpha=0.8, size=1.0) +
  #scale_color_manual(values = rep(c("gray33", "gray67"), 16 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$chromosomes, breaks= axisdf$center, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,max(-log10(don$scores))+0.3), expand = expand_scale(mult = c(0, .1))) +     # remove space between plot area and x axis
  # Add highlighted points
  
  geom_point(data=subset(don, chromosomes=="1"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="2"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="3"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="4"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="5"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="6"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="7"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="8"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="9"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="10"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="11"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="12"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="13"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="14"), color="gray67", size=1.6) +
  geom_point(data=subset(don, chromosomes=="15"), color="gray33", size=1.6) +
  geom_point(data=subset(don, chromosomes=="16"), color="gray67", size=1.6) +
  
  # Add thresholds
  geom_hline(yintercept=-log10(p.bonferroni), linetype="dashed", color = "#f04546") + 
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    #panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 22),
    title = element_text(size = 20)
  ) +
  labs(y = ("-log10(p-value)"), x = ("Chromosome"), title = ("gpd"))

p

p_colouredScaled_only90_ <- p  + geom_point(data=extralayer_only90_w_weightTop25, aes(x=BPcum, y=-log10(scores), color=Occurence,size=avgweight)) + 
  scale_size_continuous(limits=c(0,100),
                        breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                        labels=c("0","10","20","30","40","50","60","70","80","90","100"),
                        name = "Weight",
                        guide="legend",
                        range = c(0, 8)) +
  scale_color_viridis(option = "plasma") + theme(legend.position = "right") 
p_colouredScaled_only90_

ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/gpdfromgpd_Manhattan_ColourafterOccurencgpdAfterWeightOnly90perc_backgroundcol_top25',Sys.Date(),'.pdf',sep="_"), plot = p_colouredScaled_only90_, width = 75, height = 15, unit = 'cm')
ggsave(paste('/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/Figures/gpdfromgpd_Manhattan_ColourafterOccurencgpdAfterWeightOnly90perc_backgroundcol_top25',Sys.Date(),'.png',sep="_"), plot = p_colouredScaled_only90_, width = 75, height = 15, unit = 'cm')




