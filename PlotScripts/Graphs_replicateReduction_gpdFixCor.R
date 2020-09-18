
# Load libraries
{
  library(ggplot2)
  library(stringr)
}

# set working directory
{
  setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/gpd_FixCor")
}

# Load files 10 reps
{
  AllPredictionResults_10reps=list.files(pattern="^Predictions_GBLUP_GPD_10")
  
  dataframeWithResults_10reps=rep(0,7)
  dataframeWithResults_10reps=as.data.frame(dataframeWithResults_10reps)
  dataframeWithResults_10reps=t(dataframeWithResults_10reps)
  colnames(dataframeWithResults_10reps)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno","Round")

  for (i in seq(1:length(AllPredictionResults_10reps))){
    filename=AllPredictionResults_10reps[i]
    round=str_sub(filename, 35, 36)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,35,37)=="100"){
      round=str_sub(filename,35,37)
    }
    table=read.table(AllPredictionResults_10reps[i])
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    
    dataframeWithResults_10reps=rbind(dataframeWithResults_10reps,table_df)
  }
  
  dataframeWithResults_10reps=dataframeWithResults_10reps[-1,]
  dataframeWithResults_10reps=dataframeWithResults_10reps[order(dataframeWithResults_10reps$Round),]
  
}

# Load files 9 reps
{
  AllPredictionResults_9reps=list.files(pattern="^Predictions_GBLUP_GPD_9")
  
  dataframeWithResults_9reps=rep(0,7)
  dataframeWithResults_9reps=as.data.frame(dataframeWithResults_9reps)
  dataframeWithResults_9reps=t(dataframeWithResults_9reps)
  colnames(dataframeWithResults_9reps)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno","Round")
  
  for (i in seq(1:length(AllPredictionResults_9reps))){
    filename=AllPredictionResults_9reps[i]
    round=str_sub(filename, 34, 35)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,34,36)=="100"){
      round=str_sub(filename,34,36)
    }
    table=read.table(AllPredictionResults_9reps[i])
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    
    dataframeWithResults_9reps=rbind(dataframeWithResults_9reps,table_df)
  }
  
  dataframeWithResults_9reps=dataframeWithResults_9reps[-1,]
  dataframeWithResults_9reps=dataframeWithResults_9reps[order(dataframeWithResults_9reps$Round),]
  
}

# Load files 8 reps
{
  AllPredictionResults_8reps=list.files(pattern="^Predictions_GBLUP_GPD_8")
  
  dataframeWithResults_8reps=rep(0,7)
  dataframeWithResults_8reps=as.data.frame(dataframeWithResults_8reps)
  dataframeWithResults_8reps=t(dataframeWithResults_8reps)
  colnames(dataframeWithResults_8reps)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno","Round")
  
  for (i in seq(1:length(AllPredictionResults_8reps))){
    filename=AllPredictionResults_8reps[i]
    round=str_sub(filename, 34, 35)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,34,36)=="100"){
      round=str_sub(filename,34,36)
    }
    table=read.table(AllPredictionResults_8reps[i])
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    
    dataframeWithResults_8reps=rbind(dataframeWithResults_8reps,table_df)
  }
  
  dataframeWithResults_8reps=dataframeWithResults_8reps[-1,]
  dataframeWithResults_8reps=dataframeWithResults_8reps[order(dataframeWithResults_8reps$Round),]
}

# Load files 7 reps
{
  AllPredictionResults_7reps=list.files(pattern="^Predictions_GBLUP_GPD_7")
  
  dataframeWithResults_7reps=rep(0,7)
  dataframeWithResults_7reps=as.data.frame(dataframeWithResults_7reps)
  dataframeWithResults_7reps=t(dataframeWithResults_7reps)
  colnames(dataframeWithResults_7reps)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno","Round")
  
  for (i in seq(1:length(AllPredictionResults_7reps))){
    filename=AllPredictionResults_7reps[i]
    round=str_sub(filename, 34, 35)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,34,36)=="100"){
      round=str_sub(filename,34,36)
    }
    table=read.table(AllPredictionResults_7reps[i])
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    
    dataframeWithResults_7reps=rbind(dataframeWithResults_7reps,table_df)
  }
  
  dataframeWithResults_7reps=dataframeWithResults_7reps[-1,]
  dataframeWithResults_7reps=dataframeWithResults_7reps[order(dataframeWithResults_7reps$Round),]
}

# Load files 6 reps
{
  AllPredictionResults_6reps=list.files(pattern="^Predictions_GBLUP_GPD_6")
  
  dataframeWithResults_6reps=rep(0,7)
  dataframeWithResults_6reps=as.data.frame(dataframeWithResults_6reps)
  dataframeWithResults_6reps=t(dataframeWithResults_6reps)
  colnames(dataframeWithResults_6reps)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno","Round")
  
  for (i in seq(1:length(AllPredictionResults_6reps))){
    filename=AllPredictionResults_6reps[i]
    round=str_sub(filename, 34, 35)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,34,36)=="100"){
      round=str_sub(filename,34,36)
    }
    table=read.table(AllPredictionResults_6reps[i])
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    
    dataframeWithResults_6reps=rbind(dataframeWithResults_6reps,table_df)
  }
  
  dataframeWithResults_6reps=dataframeWithResults_6reps[-1,]
  dataframeWithResults_6reps=dataframeWithResults_6reps[order(dataframeWithResults_6reps$Round),]
  
}

# Load files 5 reps
{
  AllPredictionResults_5reps=list.files(pattern="^Predictions_GBLUP_GPD_5")
  
  dataframeWithResults_5reps=rep(0,7)
  dataframeWithResults_5reps=as.data.frame(dataframeWithResults_5reps)
  dataframeWithResults_5reps=t(dataframeWithResults_5reps)
  colnames(dataframeWithResults_5reps)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno","Round")
  
  for (i in seq(1:length(AllPredictionResults_5reps))){
    filename=AllPredictionResults_5reps[i]
    round=str_sub(filename, 34, 35)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,34,36)=="100"){
      round=str_sub(filename,34,36)
    }
    table=read.table(AllPredictionResults_5reps[i])
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    
    dataframeWithResults_5reps=rbind(dataframeWithResults_5reps,table_df)
  }
  
  dataframeWithResults_5reps=dataframeWithResults_5reps[-1,]
  dataframeWithResults_5reps=dataframeWithResults_5reps[order(dataframeWithResults_5reps$Round),]
  
  
}

# Load files 4 reps
{
  AllPredictionResults_4reps=list.files(pattern="^Predictions_GBLUP_GPD_4")
  
  dataframeWithResults_4reps=rep(0,7)
  dataframeWithResults_4reps=as.data.frame(dataframeWithResults_4reps)
  dataframeWithResults_4reps=t(dataframeWithResults_4reps)
  colnames(dataframeWithResults_4reps)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno","Round")
  
  for (i in seq(1:length(AllPredictionResults_4reps))){
    filename=AllPredictionResults_4reps[i]
    round=str_sub(filename, 34, 35)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,34,36)=="100"){
      round=str_sub(filename,34,36)
    }
    table=read.table(AllPredictionResults_4reps[i])
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    
    dataframeWithResults_4reps=rbind(dataframeWithResults_4reps,table_df)
  }
  
  dataframeWithResults_4reps=dataframeWithResults_4reps[-1,]
  dataframeWithResults_4reps=dataframeWithResults_4reps[order(dataframeWithResults_4reps$Round),]
  
}

# Load files 3 reps
{
  AllPredictionResults_3reps=list.files(pattern="^Predictions_GBLUP_GPD_3")
  
  dataframeWithResults_3reps=rep(0,7)
  dataframeWithResults_3reps=as.data.frame(dataframeWithResults_3reps)
  dataframeWithResults_3reps=t(dataframeWithResults_3reps)
  colnames(dataframeWithResults_3reps)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno","Round")
  
  for (i in seq(1:length(AllPredictionResults_3reps))){
    filename=AllPredictionResults_3reps[i]
    round=str_sub(filename, 34, 35)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,34,36)=="100"){
      round=str_sub(filename,34,36)
    }
    table=read.table(AllPredictionResults_3reps[i])
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    
    dataframeWithResults_3reps=rbind(dataframeWithResults_3reps,table_df)
  }
  
  dataframeWithResults_3reps=dataframeWithResults_3reps[-1,]
  dataframeWithResults_3reps=dataframeWithResults_3reps[order(dataframeWithResults_3reps$Round),]
  
}

# Load files 2 reps
{
  AllPredictionResults_2reps=list.files(pattern="^Predictions_GBLUP_GPD_2")
  
  dataframeWithResults_2reps=rep(0,7)
  dataframeWithResults_2reps=as.data.frame(dataframeWithResults_2reps)
  dataframeWithResults_2reps=t(dataframeWithResults_2reps)
  colnames(dataframeWithResults_2reps)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno","Round")
  
  for (i in seq(1:length(AllPredictionResults_2reps))){
    filename=AllPredictionResults_2reps[i]
    round=str_sub(filename, 34, 35)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,34,36)=="100"){
      round=str_sub(filename,34,36)
    }
    table=read.table(AllPredictionResults_2reps[i])
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    
    dataframeWithResults_2reps=rbind(dataframeWithResults_2reps,table_df)
  }
  
  dataframeWithResults_2reps=dataframeWithResults_2reps[-1,]
  dataframeWithResults_2reps=dataframeWithResults_2reps[order(dataframeWithResults_2reps$Round),]
  
}


# Make data for plotting

numberofrounds=nrow(dataframeWithResults_2reps)/142 


ForPlotting_10REPs = data.frame(matrix(vector(), numberofrounds, 7,
                       dimnames=list(c(), c("cor_gpdFixCor_GEBV","cor_CorrectedPheno_GEBV","Average_SD_CorrectedPheno","Average_SD_gpdFixCor","Average_SE_CorrectedPheno","Average_SE_gpdFixCor","Round"))),
                stringsAsFactors=F)



for (j in unique(dataframeWithResults_10reps$Round)){
  print(j)
  subtable=dataframeWithResults_10reps[which(dataframeWithResults_10reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_10REPs$cor_gpdFixCor_GEBV[j] = cor(subtable$gpd_rescor,subtable$GEBVs)
  ForPlotting_10REPs$cor_CorrectedPheno_GEBV[j] =  cor(subtable$Observed_correctedForFixed,subtable$GEBVs)
  ForPlotting_10REPs$Average_SD_CorrectedPheno[j] = mean(subtable$SD_CorrectedPheno)
  ForPlotting_10REPs$Average_SE_CorrectedPheno[j] = ForPlotting_10REPs$Average_SD_CorrectedPheno[j]/sqrt(10)
  ForPlotting_10REPs$Average_SD_gpdFixCor[j] = mean(subtable$SD_gpdResCor)
  ForPlotting_10REPs$Average_SE_gpdFixCor[j] = ForPlotting_10REPs$Average_SD_gpdFixCor[j]/sqrt(10)
  ForPlotting_10REPs$Round[j] = as.numeric(subtable$Round[1])
}

ForPlotting_10REPs$Replicates=rep(10,numberofrounds)



# Make data for plotting

ForPlotting_9REPs = data.frame(matrix(vector(), numberofrounds, 7,
                                       dimnames=list(c(), c("cor_gpdFixCor_GEBV","cor_CorrectedPheno_GEBV","Average_SD_CorrectedPheno","Average_SD_gpdFixCor","Average_SE_CorrectedPheno","Average_SE_gpdFixCor","Round"))),
                                stringsAsFactors=F)



for (j in unique(dataframeWithResults_9reps$Round)){
  print(j)
  subtable=dataframeWithResults_9reps[which(dataframeWithResults_9reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_9REPs$cor_gpdFixCor_GEBV[j] = cor(subtable$gpd_rescor,subtable$GEBVs)
  ForPlotting_9REPs$cor_CorrectedPheno_GEBV[j] =  cor(subtable$Observed_correctedForFixed,subtable$GEBVs)
  ForPlotting_9REPs$Average_SD_CorrectedPheno[j] = mean(subtable$SD_CorrectedPheno)
  ForPlotting_9REPs$Average_SE_CorrectedPheno[j] = ForPlotting_9REPs$Average_SD_CorrectedPheno[j]/sqrt(9)
  ForPlotting_9REPs$Average_SD_gpdFixCor[j] = mean(subtable$SD_gpdResCor)
  ForPlotting_9REPs$Average_SE_gpdFixCor[j] = ForPlotting_9REPs$Average_SD_gpdFixCor[j]/sqrt(9)
  ForPlotting_9REPs$Round[j] = as.numeric(subtable$Round[1])
}

ForPlotting_9REPs$Replicates=rep(9,numberofrounds)


# Make data for plotting
ForPlotting_8REPs = data.frame(matrix(vector(), numberofrounds, 7,
                                      dimnames=list(c(), c("cor_gpdFixCor_GEBV","cor_CorrectedPheno_GEBV","Average_SD_CorrectedPheno","Average_SD_gpdFixCor","Average_SE_CorrectedPheno","Average_SE_gpdFixCor","Round"))),
                               stringsAsFactors=F)



for (j in unique(dataframeWithResults_8reps$Round)){
  print(j)
  subtable=dataframeWithResults_8reps[which(dataframeWithResults_8reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_8REPs$cor_gpdFixCor_GEBV[j] = cor(subtable$gpd_rescor,subtable$GEBVs)
  ForPlotting_8REPs$cor_CorrectedPheno_GEBV[j] =  cor(subtable$Observed_correctedForFixed,subtable$GEBVs)
  ForPlotting_8REPs$Average_SD_CorrectedPheno[j] = mean(subtable$SD_CorrectedPheno)
  ForPlotting_8REPs$Average_SE_CorrectedPheno[j] = ForPlotting_8REPs$Average_SD_CorrectedPheno[j]/sqrt(8)
  ForPlotting_8REPs$Average_SD_gpdFixCor[j] = mean(subtable$SD_gpdResCor)
  ForPlotting_8REPs$Average_SE_gpdFixCor[j] = ForPlotting_8REPs$Average_SD_gpdFixCor[j]/sqrt(8)
  ForPlotting_8REPs$Round[j] = as.numeric(subtable$Round[1])
}

ForPlotting_8REPs$Replicates=rep(8,numberofrounds)




# Make data for plotting
ForPlotting_7REPs = data.frame(matrix(vector(), numberofrounds, 7,
                                      dimnames=list(c(), c("cor_gpdFixCor_GEBV","cor_CorrectedPheno_GEBV","Average_SD_CorrectedPheno","Average_SD_gpdFixCor","Average_SE_CorrectedPheno","Average_SE_gpdFixCor","Round"))),
                               stringsAsFactors=F)



for (j in unique(dataframeWithResults_7reps$Round)){
  print(j)
  subtable=dataframeWithResults_7reps[which(dataframeWithResults_7reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_7REPs$cor_gpdFixCor_GEBV[j] = cor(subtable$gpd_rescor,subtable$GEBVs)
  ForPlotting_7REPs$cor_CorrectedPheno_GEBV[j] =  cor(subtable$Observed_correctedForFixed,subtable$GEBVs)
  ForPlotting_7REPs$Average_SD_CorrectedPheno[j] = mean(subtable$SD_CorrectedPheno)
  ForPlotting_7REPs$Average_SE_CorrectedPheno[j] = ForPlotting_7REPs$Average_SD_CorrectedPheno[j]/sqrt(7)
  ForPlotting_7REPs$Average_SD_gpdFixCor[j] = mean(subtable$SD_gpdResCor)
  ForPlotting_7REPs$Average_SE_gpdFixCor[j] = ForPlotting_7REPs$Average_SD_gpdFixCor[j]/sqrt(7)
  ForPlotting_7REPs$Round[j] = as.numeric(subtable$Round[1])
}

ForPlotting_7REPs$Replicates=rep(7,numberofrounds)




# Make data for plotting
ForPlotting_6REPs = data.frame(matrix(vector(), numberofrounds, 7,
                                      dimnames=list(c(), c("cor_gpdFixCor_GEBV","cor_CorrectedPheno_GEBV","Average_SD_CorrectedPheno","Average_SD_gpdFixCor","Average_SE_CorrectedPheno","Average_SE_gpdFixCor","Round"))),
                               stringsAsFactors=F)



for (j in unique(dataframeWithResults_6reps$Round)){
  print(j)
  subtable=dataframeWithResults_6reps[which(dataframeWithResults_6reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_6REPs$cor_gpdFixCor_GEBV[j] = cor(subtable$gpd_rescor,subtable$GEBVs)
  ForPlotting_6REPs$cor_CorrectedPheno_GEBV[j] =  cor(subtable$Observed_correctedForFixed,subtable$GEBVs)
  ForPlotting_6REPs$Average_SD_CorrectedPheno[j] = mean(subtable$SD_CorrectedPheno)
  ForPlotting_6REPs$Average_SE_CorrectedPheno[j] = ForPlotting_6REPs$Average_SD_CorrectedPheno[j]/sqrt(6)
  ForPlotting_6REPs$Average_SD_gpdFixCor[j] = mean(subtable$SD_gpdResCor)
  ForPlotting_6REPs$Average_SE_gpdFixCor[j] = ForPlotting_6REPs$Average_SD_gpdFixCor[j]/sqrt(6)
  ForPlotting_6REPs$Round[j] = as.numeric(subtable$Round[1])
}

ForPlotting_6REPs$Replicates=rep(6,numberofrounds)



# Make data for plotting
ForPlotting_5REPs = data.frame(matrix(vector(), numberofrounds, 7,
                                      dimnames=list(c(), c("cor_gpdFixCor_GEBV","cor_CorrectedPheno_GEBV","Average_SD_CorrectedPheno","Average_SD_gpdFixCor","Average_SE_CorrectedPheno","Average_SE_gpdFixCor","Round"))),
                               stringsAsFactors=F)



for (j in unique(dataframeWithResults_5reps$Round)){
  print(j)
  subtable=dataframeWithResults_5reps[which(dataframeWithResults_5reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_5REPs$cor_gpdFixCor_GEBV[j] = cor(subtable$gpd_rescor,subtable$GEBVs)
  ForPlotting_5REPs$cor_CorrectedPheno_GEBV[j] =  cor(subtable$Observed_correctedForFixed,subtable$GEBVs)
  ForPlotting_5REPs$Average_SD_CorrectedPheno[j] = mean(subtable$SD_CorrectedPheno)
  ForPlotting_5REPs$Average_SE_CorrectedPheno[j] = ForPlotting_5REPs$Average_SD_CorrectedPheno[j]/sqrt(5)
  ForPlotting_5REPs$Average_SD_gpdFixCor[j] = mean(subtable$SD_gpdResCor)
  ForPlotting_5REPs$Average_SE_gpdFixCor[j] = ForPlotting_5REPs$Average_SD_gpdFixCor[j]/sqrt(5)
  ForPlotting_5REPs$Round[j] = as.numeric(subtable$Round[1])
}

ForPlotting_5REPs$Replicates=rep(5,numberofrounds)



# Make data for plotting
ForPlotting_4REPs = data.frame(matrix(vector(), numberofrounds, 7,
                                      dimnames=list(c(), c("cor_gpdFixCor_GEBV","cor_CorrectedPheno_GEBV","Average_SD_CorrectedPheno","Average_SD_gpdFixCor","Average_SE_CorrectedPheno","Average_SE_gpdFixCor","Round"))),
                               stringsAsFactors=F)



for (j in unique(dataframeWithResults_4reps$Round)){
  print(j)
  subtable=dataframeWithResults_4reps[which(dataframeWithResults_4reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_4REPs$cor_gpdFixCor_GEBV[j] = cor(subtable$gpd_rescor,subtable$GEBVs)
  ForPlotting_4REPs$cor_CorrectedPheno_GEBV[j] =  cor(subtable$Observed_correctedForFixed,subtable$GEBVs)
  ForPlotting_4REPs$Average_SD_CorrectedPheno[j] = mean(subtable$SD_CorrectedPheno)
  ForPlotting_4REPs$Average_SE_CorrectedPheno[j] = ForPlotting_4REPs$Average_SD_CorrectedPheno[j]/sqrt(4)
  ForPlotting_4REPs$Average_SD_gpdFixCor[j] = mean(subtable$SD_gpdResCor)
  ForPlotting_4REPs$Average_SE_gpdFixCor[j] = ForPlotting_4REPs$Average_SD_gpdFixCor[j]/sqrt(4)
  ForPlotting_4REPs$Round[j] = as.numeric(subtable$Round[1])
}

ForPlotting_4REPs$Replicates=rep(4,numberofrounds)



# Make data for plotting
ForPlotting_3REPs = data.frame(matrix(vector(), numberofrounds, 7,
                                      dimnames=list(c(), c("cor_gpdFixCor_GEBV","cor_CorrectedPheno_GEBV","Average_SD_CorrectedPheno","Average_SD_gpdFixCor","Average_SE_CorrectedPheno","Average_SE_gpdFixCor","Round"))),
                               stringsAsFactors=F)



for (j in unique(dataframeWithResults_3reps$Round)){
  print(j)
  subtable=dataframeWithResults_3reps[which(dataframeWithResults_3reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_3REPs$cor_gpdFixCor_GEBV[j] = cor(subtable$gpd_rescor,subtable$GEBVs)
  ForPlotting_3REPs$cor_CorrectedPheno_GEBV[j] =  cor(subtable$Observed_correctedForFixed,subtable$GEBVs)
  ForPlotting_3REPs$Average_SD_CorrectedPheno[j] = mean(subtable$SD_CorrectedPheno)
  ForPlotting_3REPs$Average_SE_CorrectedPheno[j] = ForPlotting_3REPs$Average_SD_CorrectedPheno[j]/sqrt(3)
  ForPlotting_3REPs$Average_SD_gpdFixCor[j] = mean(subtable$SD_gpdResCor)
  ForPlotting_3REPs$Average_SE_gpdFixCor[j] = ForPlotting_3REPs$Average_SD_gpdFixCor[j]/sqrt(3)
  ForPlotting_3REPs$Round[j] = as.numeric(subtable$Round[1])
}

ForPlotting_3REPs$Replicates=rep(3,numberofrounds)

# Make data for plotting
ForPlotting_2REPs = data.frame(matrix(vector(), numberofrounds, 7,
                                      dimnames=list(c(), c("cor_gpdFixCor_GEBV","cor_CorrectedPheno_GEBV","Average_SD_CorrectedPheno","Average_SD_gpdFixCor","Average_SE_CorrectedPheno","Average_SE_gpdFixCor","Round"))),
                               stringsAsFactors=F)



for (j in unique(dataframeWithResults_2reps$Round)){
  print(j)
  subtable=dataframeWithResults_2reps[which(dataframeWithResults_2reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_2REPs$cor_gpdFixCor_GEBV[j] = cor(subtable$gpd_rescor,subtable$GEBVs)
  ForPlotting_2REPs$cor_CorrectedPheno_GEBV[j] =  cor(subtable$Observed_correctedForFixed,subtable$GEBVs)
  ForPlotting_2REPs$Average_SD_CorrectedPheno[j] = mean(subtable$SD_CorrectedPheno)
  ForPlotting_2REPs$Average_SE_CorrectedPheno[j] = ForPlotting_2REPs$Average_SD_CorrectedPheno[j]/sqrt(2)
  ForPlotting_2REPs$Average_SD_gpdFixCor[j] = mean(subtable$SD_gpdResCor)
  ForPlotting_2REPs$Average_SE_gpdFixCor[j] = ForPlotting_2REPs$Average_SD_gpdFixCor[j]/sqrt(2)
  ForPlotting_2REPs$Round[j] = as.numeric(subtable$Round[1])
}

ForPlotting_2REPs$Replicates=rep(2,numberofrounds)


# Bind all data 

AllDataForplotting=rbind(ForPlotting_2REPs,ForPlotting_3REPs,ForPlotting_4REPs,ForPlotting_5REPs,ForPlotting_6REPs,ForPlotting_7REPs,ForPlotting_8REPs,ForPlotting_9REPs,ForPlotting_10REPs)

# check
nrow(AllDataForplotting)== numberofrounds*9







# plotting



#colnames(data)=c("replicates","cor","sterrr","lower","upper")

pd <- position_dodge(0.1) # move them .05 to the left and right

p1=ggplot(AllDataForplotting, aes(x=Replicates, y=cor_gpdFixCor_GEBV)) + 
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
  #geom_line(position=pd) +
  geom_point(position=pd) +
  ylim(min(AllDataForplotting$cor_gpdFixCor_GEBV),max(AllDataForplotting$cor_gpdFixCor_GEBV)) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Correlation") +
  theme_bw()

p1

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/cor_of_gpdFixCorAndGEBV.png', plot = p1, width = 15, height = 12, unit = 'cm')


pd <- position_dodge(0.1) # move them .05 to the left and right

p2=ggplot(AllDataForplotting, aes(x=Replicates, y=cor_CorrectedPheno_GEBV)) + 
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
  #geom_line(position=pd) +
  geom_point(position=pd) +
  ylim(min(AllDataForplotting$cor_CorrectedPheno_GEBV),max(AllDataForplotting$cor_CorrectedPheno_GEBV)) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Correlation") +
  theme_bw()

p2

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/cor_of_CorrectedPheno(gpd_FixCor)AndGEBV.png', plot = p2, width = 15, height = 12, unit = 'cm')

# Correlation between gpd_NoCor and GEBVs, means with standard deviations
Means_cor_gpdFixCor_GEBV=aggregate(AllDataForplotting$cor_gpdFixCor_GEBV, list(AllDataForplotting$Replicates), mean)
Means_cor_gpdFixCor_GEBV=as.data.frame(Means_cor_gpdFixCor_GEBV)
colnames(Means_cor_gpdFixCor_GEBV)=c("Replicates","Means_cor_gpdFixCor_GEBV")
STDEV_cor_gpdFixCor_GEBV=aggregate(AllDataForplotting$cor_gpdFixCor_GEBV, list(AllDataForplotting$Replicates), sd)
Means_cor_gpdFixCor_GEBV$lower=Means_cor_gpdFixCor_GEBV$Means_cor_gpdFixCor_GEBV-STDEV_cor_gpdFixCor_GEBV[,2]
Means_cor_gpdFixCor_GEBV$upper=Means_cor_gpdFixCor_GEBV$Means_cor_gpdFixCor_GEBV+STDEV_cor_gpdFixCor_GEBV[,2]
colnames(Means_cor_gpdFixCor_GEBV)[3:4]=c("lower","upper")

p3=ggplot(Means_cor_gpdFixCor_GEBV, aes(x=Replicates, y=Means_cor_gpdFixCor_GEBV)) + 
  geom_line(position=pd) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.05, position=pd) +
  #ylim(0,0.4) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Correlation") +
  theme_bw()

p3

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/cor_of_gpdFixCorAndGEBV_averages_withSD.png', plot = p3, width = 15, height = 12, unit = 'cm')


# Correlation between gpd_NoCor and GEBVs, means with standard errors
Means_cor_gpdFixCor_GEBV_SE=aggregate(AllDataForplotting$cor_gpdFixCor_GEBV, list(AllDataForplotting$Replicates), mean)
Means_cor_gpdFixCor_GEBV_SE=as.data.frame(Means_cor_gpdFixCor_GEBV_SE)
colnames(Means_cor_gpdFixCor_GEBV_SE)=c("Replicates","Means_cor_gpdFixCor_GEBV")
STDEV_cor_gpdFixCor_GEBV=aggregate(AllDataForplotting$cor_gpdFixCor_GEBV, list(AllDataForplotting$Replicates), sd)
STERR_cor_gpdFixCor_GEBV=STDEV_cor_gpdFixCor_GEBV/sqrt(100)

Means_cor_gpdFixCor_GEBV_SE$lower=Means_cor_gpdFixCor_GEBV_SE$Means_cor_gpdFixCor_GEBV-STERR_cor_gpdFixCor_GEBV[,2]
Means_cor_gpdFixCor_GEBV_SE$upper=Means_cor_gpdFixCor_GEBV_SE$Means_cor_gpdFixCor_GEBV+STERR_cor_gpdFixCor_GEBV[,2]
colnames(Means_cor_gpdFixCor_GEBV_SE)[3:4]=c("lower","upper")

p4=ggplot(Means_cor_gpdFixCor_GEBV_SE, aes(x=Replicates, y=Means_cor_gpdFixCor_GEBV)) + 
  geom_line(position=pd) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.05, position=pd) +
  #ylim(0,0.4) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Correlation") +
  theme_bw()

p4
ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/cor_of_gpdFixCorAndGEBV_averages_withSE.png', plot = p4, width = 15, height = 12, unit = 'cm')


# Correlation between Corrected Pheno and GEBVs, means with standard deviations
Means_cor_CorrectedPheno_GEBV=aggregate(AllDataForplotting$cor_CorrectedPheno_GEBV, list(AllDataForplotting$Replicates), mean)
Means_cor_CorrectedPheno_GEBV=as.data.frame(Means_cor_CorrectedPheno_GEBV)
colnames(Means_cor_CorrectedPheno_GEBV)=c("Replicates","Means_cor_correctedpheno_GEBV")
STDEV_cor_CorrectedPheno_GEBV=aggregate(AllDataForplotting$cor_CorrectedPheno_GEBV, list(AllDataForplotting$Replicates), sd)
Means_cor_CorrectedPheno_GEBV$lower=Means_cor_CorrectedPheno_GEBV$Means_cor_correctedpheno_GEBV-STDEV_cor_CorrectedPheno_GEBV[,2]
Means_cor_CorrectedPheno_GEBV$upper=Means_cor_CorrectedPheno_GEBV$Means_cor_correctedpheno_GEBV+STDEV_cor_CorrectedPheno_GEBV[,2]
colnames(Means_cor_CorrectedPheno_GEBV)[3:4]=c("lower","upper")

p5=ggplot(Means_cor_CorrectedPheno_GEBV, aes(x=Replicates, y=Means_cor_correctedpheno_GEBV)) + 
  geom_line(position=pd) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.05, position=pd) +
  #ylim(0,0.4) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Correlation") +
  theme_bw()

p5

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/cor_of_CorrectedPheno(gpd_NoCor)AndGEBV_averages_withSD.png', plot = p5, width = 15, height = 12, unit = 'cm')


# Correlation between gpd_ResCor and GEBVs, means with standard errors
Means_cor_CorrectedPheno_GEBV_SE=aggregate(AllDataForplotting$cor_CorrectedPheno_GEBV, list(AllDataForplotting$Replicates), mean)
Means_cor_CorrectedPheno_GEBV_SE=as.data.frame(Means_cor_CorrectedPheno_GEBV_SE)
colnames(Means_cor_CorrectedPheno_GEBV_SE)=c("Replicates","Means_cor_correctedpheno_GEBV")
STDEV_cor_CorrectedPheno_GEBV=aggregate(AllDataForplotting$cor_CorrectedPheno_GEBV, list(AllDataForplotting$Replicates), sd)
STERR_cor_gpdResCor_GEBV=STDEV_cor_CorrectedPheno_GEBV/sqrt(100)

Means_cor_CorrectedPheno_GEBV_SE$lower=Means_cor_CorrectedPheno_GEBV_SE$Means_cor_correctedpheno_GEBV-STERR_cor_gpdResCor_GEBV[,2]
Means_cor_CorrectedPheno_GEBV_SE$upper=Means_cor_CorrectedPheno_GEBV_SE$Means_cor_correctedpheno_GEBV+STERR_cor_gpdResCor_GEBV[,2]
colnames(Means_cor_CorrectedPheno_GEBV_SE)[3:4]=c("lower","upper")

p6=ggplot(Means_cor_CorrectedPheno_GEBV_SE, aes(x=Replicates, y=Means_cor_correctedpheno_GEBV)) + 
  geom_line(position=pd) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.05, position=pd) +
  #ylim(0,0.4) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Correlation") +
  theme_bw()

p6

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/cor_of_CorrectedPheno(gpd_NoCor)AndGEBV_averages_withSE.png', plot = p6, width = 15, height = 12, unit = 'cm')




# plot SD of average gpd_NoCor
Mean_avg_SD_gpdFixCor=aggregate(AllDataForplotting$Average_SD_gpdFixCor, list(AllDataForplotting$Replicates), mean)
Mean_avg_SD_gpdFixCor
Mean_avg_SD_gpdFixCor=as.data.frame(Mean_avg_SD_gpdFixCor)
colnames(Mean_avg_SD_gpdFixCor)=c("Replicates","Mean_avg_SD_gpdFixCor")

p7=ggplot(Mean_avg_SD_gpdFixCor, aes(x=Replicates, y=Mean_avg_SD_gpdFixCor)) + 
  geom_line(position=pd) +
  geom_point() +
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.05, position=pd) +
  #ylim(0,0.4) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Avg. Standard deviation of gpd_FixCor estimates") +
  theme_bw()

p7

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/Mean_avg_SD_gpdFixCor.png', plot = p7, width = 15, height = 12, unit = 'cm')


# plot SD of average corrected phenotypes
Mean_avg_SD_CorrectedPheno=aggregate(AllDataForplotting$Average_SD_CorrectedPheno, list(AllDataForplotting$Replicates), mean)
Mean_avg_SD_CorrectedPheno
Mean_avg_SD_CorrectedPheno=as.data.frame(Mean_avg_SD_CorrectedPheno)
colnames(Mean_avg_SD_CorrectedPheno)=c("Replicates","Mean_avg_SD_CorrectedPheno")

p8=ggplot(Mean_avg_SD_CorrectedPheno, aes(x=Replicates, y=Mean_avg_SD_CorrectedPheno)) + 
  geom_line(position=pd) +
  geom_point() +
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.05, position=pd) +
  #ylim(0,0.4) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Avg. Standard deviation of Corrected Phenotypes estimates") +
  theme_bw()

p8

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/Mean_avg_SD_CorrectedPheno(gpd_FixCor).png', plot = p8, width = 15, height = 12, unit = 'cm')


# plot SE of average corrected phenotypes
Mean_avg_SE_CorrectedPheno=aggregate(AllDataForplotting$Average_SE_CorrectedPheno, list(AllDataForplotting$Replicates), mean)
colnames(Mean_avg_SE_CorrectedPheno)=c("Replicates","Mean_avg_SE_CorrectedPheno")


p9=ggplot(Mean_avg_SE_CorrectedPheno, aes(x=Replicates, y=Mean_avg_SE_CorrectedPheno)) + 
  geom_line(position=pd) +
  geom_point() +
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.05, position=pd) +
  #ylim(0,0.4) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Avg. Standard error of Corrected Phenotypes estimates") +
  theme_bw()

p9

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/Mean_avg_SE_CorrectedPheno(gpd_FixCor).png', plot = p9, width = 15, height = 12, unit = 'cm')


# plot SE of average gpd_resCor 
Mean_avg_SE_gpdFixCor=aggregate(AllDataForplotting$Average_SE_gpdFixCor, list(AllDataForplotting$Replicates), mean)
colnames(Mean_avg_SE_gpdFixCor)=c("Replicates","Mean_avg_SE_gpdFixcor")

p10=ggplot(Mean_avg_SE_gpdFixCor, aes(x=Replicates, y=Mean_avg_SE_gpdFixcor)) + 
  geom_line(position=pd) +
  geom_point() +
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.05, position=pd) +
  #ylim(0,0.4) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Avg. Standard error of gpd_FixCor estimates") +
  theme_bw()

p10

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/Mean_avg_SE_gpdFixCor.png', plot = p10, width = 15, height = 12, unit = 'cm')




#When is the correlation on average significantly smaller

# Check when the average correlation is significantly smaller than what is observed at 
# number of replicates equal 10.
# I think a paired-samples t test is appropriate as the 10 groupings into trn and tst populations
# are the same for each round where a replicate is removed
# thus grouping 1 cor from rep10 should be compared to grouping 1 from rep9 etc.
# Populations compared are NOT independent

# paired-samples t test is prepared as follows:
# count the difference between each pair (d)
# COmpute the mean (m) and the sd of d
# Compare the average difference to 0

 # for cor_CorrectedPheno_GEBV
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),2], AllDataForplotting[which(AllDataForplotting$Replicates==9),2], paired = FALSE, alternative = "greater") # N.S.
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),2], AllDataForplotting[which(AllDataForplotting$Replicates==8),2], paired = FALSE, alternative = "greater") # N.S
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),2], AllDataForplotting[which(AllDataForplotting$Replicates==7),2], paired = FALSE, alternative = "greater") # N.S.
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),2], AllDataForplotting[which(AllDataForplotting$Replicates==6),2], paired = FALSE, alternative = "greater") # N.S.
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),2], AllDataForplotting[which(AllDataForplotting$Replicates==5),2], paired = FALSE, alternative = "greater") # N.S.
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),2], AllDataForplotting[which(AllDataForplotting$Replicates==4),2], paired = FALSE, alternative = "greater") # N.S.
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),2], AllDataForplotting[which(AllDataForplotting$Replicates==3),2], paired = FALSE, alternative = "greater") # N.S.
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),2], AllDataForplotting[which(AllDataForplotting$Replicates==2),2], paired = FALSE, alternative = "greater") # N.S.


write.table(AllDataForplotting,file="/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/replicate_reduction_gpdFixCor.txt",sep="\t",quote=F,col.names = T, row.names = F)
