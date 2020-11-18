# Check how often gpd_rescor produces a higher correlation than gpd_fixcor
rm(list = ls())

# Load packages
{
  library("boot")
  library("stringr")
  library("ggplot2")
}

# Set working directories and load files
{
  setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_gpd_rescor_6foldCV_20200828")
  
  dataframeWithResults=rep(0,2)
  dataframeWithResults=as.data.frame(dataframeWithResults)
  dataframeWithResults=t(dataframeWithResults)
  colnames(dataframeWithResults)=c("gpd_ResCorCorrelation","Round")
  
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
    table=read.table(AllCorrelationResults[i])
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    colnames(table_df)[1]=c("gpd_ResCorCorrelation")
    dataframeWithResults=rbind(dataframeWithResults,table_df)
  }
  
  dataframeWithResults=dataframeWithResults[2:nrow(dataframeWithResults),]
  
  
  setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_gpd_Nocor_6foldCV_20200811")
  dataframeWithResults2=rep(0,2)
  dataframeWithResults2=as.data.frame(dataframeWithResults2)
  dataframeWithResults2=t(dataframeWithResults2)
  colnames(dataframeWithResults2)=c("gpd_NoCorCorrelation","Round")
  
  # Load data from GBLUP gpd_ResCor
  AllCorrelationResults=list.files(pattern="^Correlation_")
  
  for (i in seq(1:length(AllCorrelationResults))){
    filename=AllCorrelationResults[i]
    round=str_sub(filename, 22, 23)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }else if (str_sub(filename,22,24)=="100"){
      round=str_sub(filename,22,24)
    }
    table=read.table(AllCorrelationResults[i])
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    colnames(table_df)[1]=c("gpd_NoCorCorrelation")
    dataframeWithResults2=rbind(dataframeWithResults2,table_df)
  }
  
  dataframeWithResults2=dataframeWithResults2[2:nrow(dataframeWithResults2),]
  
  
  # merge two dataframes
  MergedDf=merge(dataframeWithResults,dataframeWithResults2,by="Round")
  head(MergedDf)
  
  length(which(MergedDf$gpd_ResCorCorrelation>=MergedDf$gpd_NoCorCorrelation)) #100/100 rounds the correlation of gpd_rescor is higher than gpd_nocor
  length(which(MergedDf$gpd_ResCorCorrelation<=MergedDf$gpd_NoCorCorrelation))
  
  
  p2 = ggplot() + 
    geom_line(data = MergedDf, aes(x = Round, y = gpd_ResCorCorrelation), color = "#0B775E") +
    geom_line(data = MergedDf, aes(x = Round, y = gpd_NoCorCorrelation), color = "#35274A") +
    xlab('Round') +
    ylab('Correlation') +
    theme_classic()
  
  p2
  
  ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/CorrelationGBLUP_gpdResCor_vs._gpdNoCor.pdf', plot = p2, width = 30, height = 20, unit = 'cm')
  
  
  # try adding gpd_fixcor
  
  setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GBLUP_gpd_Fixcor_6foldCV_20201113")
  dataframeWithResults3=rep(0,2)
  dataframeWithResults3=as.data.frame(dataframeWithResults3)
  dataframeWithResults3=t(dataframeWithResults3)
  colnames(dataframeWithResults3)=c("gpd_FixCorCorrelation","Round")
  
  # Load data from GBLUP gpd_ResCor
  AllCorrelationResults=list.files(pattern="^Correlation_")
  
  for (i in seq(1:length(AllCorrelationResults))){
    filename=AllCorrelationResults[i]
    round=str_sub(filename, 22, 23)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }else if (str_sub(filename,22,24)=="100"){
      round=str_sub(filename,22,24)
    }
    table=read.table(AllCorrelationResults[i])
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    colnames(table_df)[1]=c("gpd_FixCorCorrelation")
    dataframeWithResults3=rbind(dataframeWithResults3,table_df)
  }
  
  dataframeWithResults3=dataframeWithResults3[2:nrow(dataframeWithResults3),]
  
  MergedDf2=merge(MergedDf,dataframeWithResults3,by="Round")
  
  p3 = ggplot() + 
    geom_line(data = MergedDf2, aes(x = Round, y = gpd_ResCorCorrelation), color = "#0B775E") +
    geom_line(data = MergedDf2, aes(x = Round, y = gpd_NoCorCorrelation), color = "#35274A") +
    geom_line(data = MergedDf2, aes(x = Round, y = gpd_FixCorCorrelation), color = "#EABE94") +
    xlab('Round') +
    ylab('Correlation') +
    theme_classic()
  
  p3
  ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/CorrelationGBLUP_gpdResCor_vs._gpdNoCor_vs.gpdFixCor_20201118.pdf', plot = p3, width = 30, height = 20, unit = 'cm')
}
