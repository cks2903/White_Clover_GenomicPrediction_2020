##################################################################
#       Display replicate reduction                              #
##################################################################

setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/ReplicateReduction_20210222/F0/GBLUP/Results")

{
  library(ggplot2)
  library(stringr)
}

####################

read_files = function(n,trait){
  #A function reading in all 
  # correlation files with 
  # n replicates for 
  # a given trait
  
  FilesOfInterest=list.files(pattern=paste("^Correlation_GBLUP_",trait,"_",n,"Replicates",sep=""))
  
  for (file in FilesOfInterest){

    
    if (!exists("dataframe")){
      currentfile= read.table(file,head=F,sep="\t")
      round_= strsplit(file,"Replicates")[[1]][2]
      round=strsplit(round_,".txt")[[1]]  
      dataframe = cbind(round,currentfile)
      colnames(dataframe) =c("Round","Cor")
      
      }
  
    # if the merged dataset does exist, append to it
    else{
      currentfile <- read.table(file,head=F,sep="\t")
      round_= strsplit(file,"Replicates")[[1]][2]
      round=strsplit(round_,".txt")[[1]]  
      pred_dataframe = cbind(round,currentfile)
      colnames(pred_dataframe) =c("Round","Cor")
    
      dataframe<-rbind(dataframe, pred_dataframe)
      rm(pred_dataframe)
    }
  }
    
    return(dataframe)
}



read_predfiles = function(n,trait){
  #A function reading in all 
  # correlation files with 
  # n replicates for 
  # a given trait
  
  FilesOfInterest=list.files(pattern=paste("^Predictions_GBLUP_",trait,"_",n,"Replicates",sep=""))
  
  for (file in FilesOfInterest){
    
    if (!exists("dataframe")){
      currentfile= read.table(file,head=T,sep="\t")
      round_= strsplit(file,"Replicates")[[1]][2]
      round=strsplit(round_,".txt")[[1]]  
      dataframe = cbind(round,currentfile)
      
    }
    
    # if the merged dataset does exist, append to it
    else{
      currentfile <- read.table(file,head=T,sep="\t")
      round_= strsplit(file,"Replicates")[[1]][2]
      round=strsplit(round_,".txt")[[1]]  
      pred_dataframe = cbind(round,currentfile)
      
      dataframe<-rbind(dataframe, pred_dataframe)
      rm(pred_dataframe)
    }
  }
  
  return(dataframe)
}



# gpd, GBLUP, F0

dataframeWithResults=data.frame(matrix(NA, nrow = 10, ncol = 3))
colnames(dataframeWithResults)=c("Reps","Cor","SE")

for (i in seq(1:10)){
  ReadData = read_files(i,"gpd")
  dataframeWithResults$Reps[i]= as.numeric(i)
  dataframeWithResults$Cor[i]= mean(na.omit(ReadData[,2]))
  dataframeWithResults$SE[i]= (sd(na.omit(ReadData[,2]))/(nrow(na.omit(ReadData))^0.5))
  
}

dataframeWithResults



p1=ggplot(dataframeWithResults, aes(x=Reps, y=Cor)) + 
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
  geom_line() +
  geom_point() +
  ylim(min(dataframeWithResults$Cor)-0.02,max(dataframeWithResults$Cor)+0.02) +
  geom_errorbar(aes(ymin=Cor-SE, ymax=Cor+SE), width=.2,
                position=position_dodge(.9))  +
  scale_x_reverse(breaks=c(seq(1,10))) +
  xlab("Replicates") +
  ylab("Correlation") +
  labs(title="F0, GBLUP, gpd")+
  theme_classic()

p1

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/ReplicateReduction_20210222/F0/RepRed_F0_GBLUP_gpd.pdf', plot = p1, width = 15, height = 15, unit = 'cm')

# iSize, GBLUP, F0

dataframeWithResults=data.frame(matrix(NA, nrow = 10, ncol = 3))
colnames(dataframeWithResults)=c("Reps","Cor","SE")

for (i in seq(1:10)){
  ReadData = read_files(i,"iSize")
  dataframeWithResults$Reps[i]= as.numeric(i)
  dataframeWithResults$Cor[i]= mean(na.omit(ReadData[,2]))
  dataframeWithResults$SE[i]= (sd(na.omit(ReadData[,2]))/(nrow(na.omit(ReadData))^0.5))
  
}

dataframeWithResults



p1=ggplot(dataframeWithResults, aes(x=Reps, y=Cor)) + 
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
  geom_line() +
  geom_point() +
  ylim(min(dataframeWithResults$Cor)-0.02,max(dataframeWithResults$Cor)+0.02) +
  geom_errorbar(aes(ymin=Cor-SE, ymax=Cor+SE), width=.2,
                position=position_dodge(.9))  +
  scale_x_reverse(breaks=c(seq(1,10))) +
  xlab("Replicates") +
  ylab("Correlation") +
  labs(title="F0, GBLUP, iSize")+
  theme_classic()

p1

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/ReplicateReduction_20210222/F0/RepRed_F0_GBLUP_iSize.pdf', plot = p1, width = 15, height = 15, unit = 'cm')





# t-test to see when correlation is significantly lower

ttest_function = function(a,b,trait){
  # testing if a number of replicates give a higher 
  #correlation than b number of replicates
  # using a t-test
  ReadData1 = read_files(a,trait)
  ReadData2 = read_files(b,trait)
  
  test_results = t.test(ReadData1$Cor,ReadData2$Cor,paired = TRUE, alternative = "greater")
  pvalue=test_results$p.value
  return(pvalue)
}

ttest_function(10,9,"iSize") # ***
ttest_function(10,8,"iSize") # ***
ttest_function(10,7,"iSize") # ***
ttest_function(10,6,"iSize") # ***
ttest_function(10,5,"iSize") # ***
ttest_function(10,4,"iSize") # ***
ttest_function(10,3,"iSize") # ***
ttest_function(10,2,"iSize") # ***
ttest_function(10,1,"iSize") # ***



FDR_function = function(a,b,trait){
  # testing if a number of replicates give a higher 
  #correlation than b number of replicates
  # using an FDR approach
  ReadData1 = read_files(a,trait)
  ReadData2 = read_files(b,trait)
  
  
  if (all(ReadData1$Round==ReadData2$Round)){
    count=0
    for (i in seq(1:100)){
      if(is.na(ReadData1$Cor[i]<ReadData2$Cor[i])) {
        check = FALSE
      }else{
        check = ReadData1$Cor[i]<ReadData2$Cor[i]
      }
        
      if (check ==TRUE){
        count=count+1
      }
    }
  }
  fdr=count/100
  return(fdr)
}

FDR_function(10,9,"iSize") # 0.23 (means that the correlation from 9reps is better than from 10 reps in 23/100 cases)
FDR_function(10,8,"iSize") # 0.14
FDR_function(10,7,"iSize") # 0.06
FDR_function(10,6,"iSize") # 0.05
FDR_function(10,5,"iSize") # 0.02
FDR_function(10,4,"iSize") # 0.01
FDR_function(10,3,"iSize") # 0.01
FDR_function(10,2,"iSize") # 0.00
FDR_function(10,1,"iSize") # 0.00



FDR_function(10,9,"gpd") # 0.37
FDR_function(10,8,"gpd") # 0.29
FDR_function(10,7,"gpd") # 0.25
FDR_function(10,6,"gpd") # 0.22
FDR_function(10,5,"gpd") # 0.13
FDR_function(10,4,"gpd") # 0.14
FDR_function(10,3,"gpd") # 0.13
FDR_function(10,2,"gpd") # 0.07
FDR_function(10,1,"gpd") # 0.04

ttest_function(10,9,"gpd") # ***
ttest_function(10,8,"gpd") # ***
ttest_function(10,7,"gpd") # ***
ttest_function(10,6,"gpd") # ***
ttest_function(10,5,"gpd") # ***
ttest_function(10,4,"gpd") # ***
ttest_function(10,3,"gpd") # ***
ttest_function(10,2,"gpd") # ***
ttest_function(10,1,"gpd") # ***




