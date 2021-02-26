##################################################################
#       Display replicate reduction                              #
##################################################################

setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/ReplicateReduction_20210222/F1/Results")

{
  library(ggplot2)
  library(stringr)
  library(psych)
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
  ylim(min(dataframeWithResults$Cor)-0.03,max(dataframeWithResults$Cor)+0.03) +
  geom_errorbar(aes(ymin=Cor-SE, ymax=Cor+SE), width=.2,
                position=position_dodge(.9))  +
  scale_x_reverse(breaks=c(seq(1,10))) +
  xlab("Replicates") +
  ylab("Correlation") +
  labs(title="F1, GBLUP, gpd")+
  theme_classic()

p1

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/ReplicateReduction_20210222/F1/RepRed_F1_GBLUP_gpd.pdf', plot = p1, width = 15, height = 15, unit = 'cm')



# iSize, GBLUP, F1

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
  ylim(min(dataframeWithResults$Cor)-0.03,max(dataframeWithResults$Cor)+0.03) +
  geom_errorbar(aes(ymin=Cor-SE, ymax=Cor+SE), width=.2,
                position=position_dodge(.9))  +
  scale_x_reverse(breaks=c(seq(1,10))) +
  xlab("Replicates") +
  ylab("Correlation") +
  labs(title="F1, GBLUP, iSize")+
  theme_classic()

p1

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/ReplicateReduction_20210222/F1/RepRed_F1_GBLUP_iSize.pdf', plot = p1, width = 15, height = 15, unit = 'cm')





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



HW_function = function(a,b,trait){
  # testing if a number of replicates give a higher 
  #correlation than b number of replicates
  # using a Hoteling-Williams test
  ReadData1 = read_predfiles(a,trait)
  ReadData2 = read_predfiles(b,trait)
  
  r12 = cor(ReadData1$X20200226_dryweight,ReadData1$GEBV)
  r13 = cor(ReadData1$X20200226_dryweight,ReadData2$GEBV)
  r23 = cor(ReadData1$GEBV,ReadData2$GEBV)
  
  n = 100*9
  testresults = r.test(n, r12, r13, r23,twotailed = FALSE) 
  pvalue=testresults$p
  return(pvalue)
}


HW_function(10,9,"iSize") # *
HW_function(10,8,"iSize") # ***
HW_function(10,7,"iSize") # ***
HW_function(10,6,"iSize") # ***
HW_function(10,5,"iSize") # ***
HW_function(10,4,"iSize") # ***
HW_function(10,3,"iSize") # ***
HW_function(10,2,"iSize") # ***
HW_function(10,1,"iSize") # ***



FDR_function = function(a,b,trait){
  # testing if a number of replicates give a higher 
  #correlation than b number of replicates
  # using an FDR approach
  ReadData1 = read_files(a,trait)
  ReadData2 = read_files(b,trait)
  
  
  if (all(ReadData1$Round==ReadData2$Round)){
    count=0
    for (i in seq(1:100)){
      if(ReadData1$Cor[i]<ReadData2$Cor[i]){
        count=count+1
      }
    }
  }
  fdr=count/100
  return(fdr)
}

FDR_function(10,9,"iSize") # 0.42 (means that the correlation from 9reps is better than from 10 reps in 42/100 cases)
FDR_function(10,8,"iSize") # 0.43
FDR_function(10,7,"iSize") # 0.34
FDR_function(10,6,"iSize") # 0.29
FDR_function(10,5,"iSize") # 0.27
FDR_function(10,4,"iSize") # 0.26
FDR_function(10,3,"iSize") # 0.13
FDR_function(10,2,"iSize") # 0.14
FDR_function(10,1,"iSize") # 0.06



FDR_function(10,9,"gpd") # 0.39
FDR_function(10,8,"gpd") # 0.30
FDR_function(10,7,"gpd") # 0.20
FDR_function(10,6,"gpd") # 0.16
FDR_function(10,5,"gpd") # 0.22
FDR_function(10,4,"gpd") # 0.15
FDR_function(10,3,"gpd") # 0.11
FDR_function(10,2,"gpd") # 0.06
FDR_function(10,1,"gpd") # 0.01


HW_function(10,9,"gpd") # *
HW_function(10,8,"gpd") # 
HW_function(10,7,"gpd") # 
HW_function(10,6,"gpd") # 
HW_function(10,5,"gpd") # 
HW_function(10,4,"gpd") # 
HW_function(10,3,"gpd") # 
HW_function(10,2,"gpd") # 
HW_function(10,1,"gpd") # 

