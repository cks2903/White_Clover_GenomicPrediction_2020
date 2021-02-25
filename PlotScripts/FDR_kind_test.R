# FDR kind of test

library(gtools)




FDR <- function(path1,prefix1,path2,prefix2){
  # iSize from iSize, top 25
  setwd(path1)
  file_list <- list.files(pattern=prefix1)
  file_list1=mixedsort(sort(file_list))

  for (i in seq(1:length(file_list1))){
    
    file=file_list1[i]

    # if the merged dataset doesn't exist, create it
    if (i==1){
      df <- read.table(file_list1[i], header=TRUE, sep=",")
    
      # if the merged dataset does exist, append to it
    }else{
      pred_dataset <-read.table(file_list1[i], header=TRUE, sep=",")
      df<-rbind(df, pred_dataset)
      rm(pred_dataset)
    }
  }

    df=as.data.frame(df)
    colnames(df)=c("Correlation")
    nrow(df)==100 #check


    # iSize from iSize, random 25
    setwd(path2)

  file_list <- list.files(pattern=prefix2)
  file_list1=mixedsort(sort(file_list))

  for (i in seq(1:length(file_list1))){
  
    # if the merged dataset doesn't exist, create it
    if (i==1){
      df1 <- read.table(file_list1[i], header=TRUE, sep=",")
    
    
      # if the merged dataset does exist, append to it
    }else{
      pred_dataset <-read.table(file_list1[i], header=TRUE, sep=",")
      df1<-rbind(df1, pred_dataset)
      rm(pred_dataset)
    }
  }

  df1=as.data.frame(df1)
  colnames(df1)=c("Correlation")
  nrow(df1)==100 #check

  count=0
  for (i in seq(1:nrow(df1))){
  
    if (df1$Cor[i] > df$Cor[i]){
      count=count+1
      print(i)
    }
    else{
      next
    }
  }

  FDR = count/100

  return(FDR)
}


FDR_compGBLUP <- function(path1,prefix1,path2,prefix2){
  # iSize from iSize, top 25
  setwd(path1)
  file_list <- list.files(pattern=prefix1)
  file_list1=mixedsort(sort(file_list))
  
  for (i in seq(1:length(file_list1))){
    
    file=file_list1[i]
    
    # if the merged dataset doesn't exist, create it
    if (i==1){
      df <- read.table(file_list1[i], header=TRUE, sep=",")
      
      # if the merged dataset does exist, append to it
    }else{
      pred_dataset <-read.table(file_list1[i], header=TRUE, sep=",")
      df<-rbind(df, pred_dataset)
      rm(pred_dataset)
    }
  }
  
  df=as.data.frame(df)
  colnames(df)=c("Correlation")
  nrow(df)==100 #check
  
  
  # iSize from iSize, random 25
  setwd(path2)
  
  file_list <- list.files(pattern=prefix2)
  file_list1=mixedsort(sort(file_list))
  
  for (i in seq(1:length(file_list1))){
    
    # if the merged dataset doesn't exist, create it
    if (i==1){
      df1 <- read.table(file_list1[i], header=FALSE, sep="\t",stringsAsFactors=FALSE)
      
      
      # if the merged dataset does exist, append to it
    }else{
      pred_dataset <-read.table(file_list1[i], header=FALSE, sep="\t",stringsAsFactors=FALSE)
      df1<-rbind(df1, pred_dataset)
      rm(pred_dataset)
    }
  }
  
  df1=as.data.frame(df1)
  colnames(df1)=c("Correlation")
  nrow(df1)==100 #check
  
  count=0
  for (i in seq(1:nrow(df1))){
    
    if (df1$Cor[i] > df$Cor[i]){
      count=count+1
      print(i)
    }
    else{
      next
    }
  }
  
  FDR = count/100
  
  return(FDR)
}

# testing how often the last measurement (random) performs better than the first (top) 

# FDR for top25 vs. random25 predicting iSizefromiSize
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_iSizeFromiSize"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_iSizeFromiSize"
prefix1="Correlations_iSize_RF_top25SNP"
prefix2="Correlations_iSize_RF_random25SNP"

FDR(path1,prefix1,path2,prefix2) #0.0


# FDR for top200 vs. random200 predicting iSizefromiSize
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_iSizeFromiSize"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_iSizeFromiSize"
prefix1="Correlations_iSize_RF_top200SNP"
prefix2="Correlations_iSize_RF_random200SNP"

FDR(path1,prefix1,path2,prefix2) #0.03 

# FDR for top25 vs. random25 predicting gpdfromiSize
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_gpdFromiSize"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_gpdFromiSize"
prefix1="Correlations_gpd_RF_top25SNP"
prefix2="Correlations_gpd_RF_random25SNP"

FDR(path1,prefix1,path2,prefix2) #0

# FDR for top200 vs. random200 predicting gpdfromiSize
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_gpdFromiSize"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_gpdFromiSize"
prefix1="Correlations_gpd_RF_top200SNP"
prefix2="Correlations_gpd_RF_random200SNP"

FDR(path1,prefix1,path2,prefix2) #0


# FDR for top25 vs. random25 predicting gpdfromgpd
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_gpdFromgpd"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_gpdFromgpd"
prefix1="Correlations_gpd_RF_top25SNP"
prefix2="Correlations_gpd_RF_random25SNP"

FDR(path1,prefix1,path2,prefix2) #0.07

# FDR for top200 vs. random200 predicting gpdfromgpd
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_gpdFromgpd"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_gpdFromgpd"
prefix1="Correlations_gpd_RF_top200SNP"
prefix2="Correlations_gpd_RF_random200SNP"

FDR(path1,prefix1,path2,prefix2) #0.01


# FDR for top25 vs. random25 predicting iSizefromgpd
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_iSizeFromgpd"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_iSizeFromgpd"
prefix1="Correlations_iSize_RF_top25SNP"
prefix2="Correlations_iSize_RF_random25SNP"

FDR(path1,prefix1,path2,prefix2)  #0.0

# FDR for top200 vs. random200 predicting isizefromgpd
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_iSizeFromgpd"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_iSizeFromgpd"
prefix1="Correlations_iSize_RF_top200SNP"
prefix2="Correlations_iSize_RF_random200SNP"

FDR(path1,prefix1,path2,prefix2) #0.01








# now compare between methods



# FDR for top200 and top 25 isizefromgpd 
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_iSizeFromgpd"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_iSizeFromgpd"
prefix1="Correlations_iSize_RF_top200SNP"
prefix2="Correlations_iSize_RF_top25SNP"

FDR(path1,prefix1,path2,prefix2) #0.06



# FDR for top200 and top 25 isizefromisize
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_iSizeFromiSize"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_iSizeFromiSize"
prefix1="Correlations_iSize_RF_top200SNP"
prefix2="Correlations_iSize_RF_top25SNP"

FDR(path1,prefix1,path2,prefix2) #0.62 = 0.38


# FDR for top200 fprom isizefromgpd  top200 fprom isizefromisize
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_iSizeFromgpd"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_iSizeFromiSize"
prefix1="Correlations_iSize_RF_top200SNP"
prefix2="Correlations_iSize_RF_top200SNP"

FDR(path1,prefix1,path2,prefix2) #0.15


# FDR for top200 fprom isizefromgpd  top200 fprom isizefromisize
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_iSizeFromgpd"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_iSizeFromiSize"
prefix1="Correlations_iSize_RF_top200SNP"
prefix2="Correlations_iSize_RF_top200SNP"

FDR(path1,prefix1,path2,prefix2) #0.15



# FDR for top200 from isizefromgpd  vs. gblup
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_iSizeFromgpd"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GBLUP_iSizeAvg_20210120"
prefix1="Correlations_iSize_RF_top200SNP"
prefix2="Correlation_GBLUP"

FDR_compGBLUP(path1,prefix1,path2,prefix2) #0.88 so GLUP performs better in 88/100 rounds

# FDR for top200 from isizefromisize  vs. gblup
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_iSizeFromiSize"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GBLUP_iSizeAvg_20210120"
prefix1="Correlations_iSize_RF_top200SNP"
prefix2="Correlation_GBLUP"

FDR_compGBLUP(path1,prefix1,path2,prefix2) #098 so GLUP performs better in 98/100 rounds


# FDR for top25 from isizefromgpd  vs. gblup
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_iSizeFromgpd"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GBLUP_iSizeAvg_20210120"
prefix1="Correlations_iSize_RF_top25SNP"
prefix2="Correlation_GBLUP"

FDR_compGBLUP(path1,prefix1,path2,prefix2) #0.99 so GLUP performs better in 99/100 rounds


# FDR for top25 from isizefromisize  vs. gblup
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_iSizeFromiSize"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GBLUP_iSizeAvg_20210120"
prefix1="Correlations_iSize_RF_top25SNP"
prefix2="Correlation_GBLUP"

FDR_compGBLUP(path1,prefix1,path2,prefix2) #0.87 so GLUP performs better in 87/100 rounds








# FDR for top25 from gpd from iSize  vs. gblup
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_gpdFromiSize"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GBLUP_gpdAvg_20210120"
prefix1="Correlations_gpd_RF_top25SNP"
prefix2="Correlation_GBLUP"
FDR_compGBLUP(path1,prefix1,path2,prefix2) #0.33 so GBLUP performs better in 33/100 rounds

# FDR for top200 from gpd from iSize  vs. gblup
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Correlations_gpdFromiSize"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GBLUP_gpdAvg_20210120"
prefix1="Correlations_gpd_RF_top200SNP"
prefix2="Correlation_GBLUP"
FDR_compGBLUP(path1,prefix1,path2,prefix2) #0.57 so GBLUP performs better in 57/100 rounds


# FDR for top25 from gpd from gpd  vs. gblup
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_gpdFromgpd"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GBLUP_gpdAvg_20210120"
prefix1="Correlations_gpd_RF_top25SNP"
prefix2="Correlation_GBLUP"
FDR_compGBLUP(path1,prefix1,path2,prefix2) #0.98 so GBLUP performs better in 98/100 rounds



# FDR for top200 from gpd from gpd  vs. gblup
path1="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_gpd_followedByRF_20210201/Correlations_gpdFromgpd"
path2="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GBLUP_gpdAvg_20210120"
prefix1="Correlations_gpd_RF_top200SNP"
prefix2="Correlation_GBLUP"
FDR_compGBLUP(path1,prefix1,path2,prefix2) #0.89 so GBLUP performs better in 89/100 rounds
