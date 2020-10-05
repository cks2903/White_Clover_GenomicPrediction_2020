# Load libraries
{
  library(ggplot2)
  library(stringr)
}

# set working directory
{
  setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/F1_pred_gpd_ResCor")
}

# load functions
{
  FDRfunc <- function(repA,repB){
    count=0
    for (n_round in seq(1:100)){
      Fixedround=AllDataForplotting[which(AllDataForplotting$Round==n_round),]
      statement=Fixedround[which(Fixedround$Replicates==repA),1]<=Fixedround[which(Fixedround$Replicates==repB),1]
      if (statement==T){
        count=count+1
      }
    }
    FDR=count/100
    return(format(round(FDR,3),nsmall = 2))
  }
  
}

# Load files 10 reps
{
  AllPredictionResults_10reps=list.files(pattern="^Correlation_Correctedphenotype_GBLUP_GPD_10Rep")
  
  dataframeWithResults_10reps=rep(0,6)
  dataframeWithResults_10reps=as.data.frame(dataframeWithResults_10reps)
  dataframeWithResults_10reps=t(dataframeWithResults_10reps)
  colnames(dataframeWithResults_10reps)=c("F1.population","Mean.Parental.predicted.GEBV","F1.pop.mean.dryweight","F1.pop.mean.freshweight","Round","Replicates")
  
  for (i in seq(1:length(AllPredictionResults_10reps))){
    filename=AllPredictionResults_10reps[i]
    round=str_sub(filename, 54, 55)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,54,56)=="100"){
      round=str_sub(filename,54,56)
    }
    
    table=read.table(AllPredictionResults_10reps[i],sep="\t",head=T)
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    table_df$Replicates=as.numeric(10)
    
    notused=c("DLF2 NR high yield","DLF3 NR high yield","DLF4 NR high yield","DLF5 NR high yield","YellowTip1 SM yellow tip")
    rmv=which(table_df$F1.population %in% notused)  # rows to delete
    table_df=table_df[-rmv,]
    
    dataframeWithResults_10reps=rbind(dataframeWithResults_10reps,table_df)
  }
  
  dataframeWithResults_10reps=dataframeWithResults_10reps[-1,]
  dataframeWithResults_10reps

}

# Load files 9 reps
{
  AllPredictionResults_9reps=list.files(pattern="^Correlation_Correctedphenotype_GBLUP_GPD_9Rep")
  
  dataframeWithResults_9reps=rep(0,6)
  dataframeWithResults_9reps=as.data.frame(dataframeWithResults_9reps)
  dataframeWithResults_9reps=t(dataframeWithResults_9reps)
  colnames(dataframeWithResults_9reps)=c("F1.population","Mean.Parental.predicted.GEBV","F1.pop.mean.dryweight","F1.pop.mean.freshweight","Round","Replicates")
  
  for (i in seq(1:length(AllPredictionResults_9reps))){
    filename=AllPredictionResults_9reps[i]
    round=str_sub(filename, 53, 54)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,53,55)=="100"){
      round=str_sub(filename,53,55)
    }
    
    table=read.table(AllPredictionResults_9reps[i],sep="\t",head=T)
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    table_df$Replicates=as.numeric(9)
    
    notused=c("DLF2 NR high yield","DLF3 NR high yield","DLF4 NR high yield","DLF5 NR high yield","YellowTip1 SM yellow tip")
    rmv=which(table_df$F1.population %in% notused)  # rows to delete
    table_df=table_df[-rmv,]
    
    dataframeWithResults_9reps=rbind(dataframeWithResults_9reps,table_df)
  }
  
  dataframeWithResults_9reps=dataframeWithResults_9reps[-1,]
  dataframeWithResults_9reps
  
}

# Load files 8 reps
{
  AllPredictionResults_8reps=list.files(pattern="^Correlation_Correctedphenotype_GBLUP_GPD_8Rep")
  
  dataframeWithResults_8reps=rep(0,6)
  dataframeWithResults_8reps=as.data.frame(dataframeWithResults_8reps)
  dataframeWithResults_8reps=t(dataframeWithResults_8reps)
  colnames(dataframeWithResults_8reps)=c("F1.population","Mean.Parental.predicted.GEBV","F1.pop.mean.dryweight","F1.pop.mean.freshweight","Round","Replicates")
  
  for (i in seq(1:length(AllPredictionResults_8reps))){
    filename=AllPredictionResults_8reps[i]
    round=str_sub(filename, 53, 54)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,53,55)=="100"){
      round=str_sub(filename,53,55)
    }
    
    table=read.table(AllPredictionResults_8reps[i],sep="\t",head=T)
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    table_df$Replicates=as.numeric(8)
    
    notused=c("DLF2 NR high yield","DLF3 NR high yield","DLF4 NR high yield","DLF5 NR high yield","YellowTip1 SM yellow tip")
    rmv=which(table_df$F1.population %in% notused)  # rows to delete
    table_df=table_df[-rmv,]
    
    dataframeWithResults_8reps=rbind(dataframeWithResults_8reps,table_df)
  }
  
  dataframeWithResults_8reps=dataframeWithResults_8reps[-1,]
  dataframeWithResults_8reps
  
}

# Load files 7 reps
{
  AllPredictionResults_7reps=list.files(pattern="^Correlation_Correctedphenotype_GBLUP_GPD_7Rep")
  
  dataframeWithResults_7reps=rep(0,6)
  dataframeWithResults_7reps=as.data.frame(dataframeWithResults_7reps)
  dataframeWithResults_7reps=t(dataframeWithResults_7reps)
  colnames(dataframeWithResults_7reps)=c("F1.population","Mean.Parental.predicted.GEBV","F1.pop.mean.dryweight","F1.pop.mean.freshweight","Round","Replicates")
  
  for (i in seq(1:length(AllPredictionResults_7reps))){
    filename=AllPredictionResults_7reps[i]
    round=str_sub(filename, 53, 54)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,53,55)=="100"){
      round=str_sub(filename,53,55)
    }
    table=read.table(AllPredictionResults_7reps[i],sep="\t",head=T)
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    table_df$Replicates=as.numeric(7)
    
    notused=c("DLF2 NR high yield","DLF3 NR high yield","DLF4 NR high yield","DLF5 NR high yield","YellowTip1 SM yellow tip")
    rmv=which(table_df$F1.population %in% notused)  # rows to delete
    table_df=table_df[-rmv,]
    
    dataframeWithResults_7reps=rbind(dataframeWithResults_7reps,table_df)
  }
  
  dataframeWithResults_7reps=dataframeWithResults_7reps[-1,]
  dataframeWithResults_7reps
  
}

# Load files 6 reps
{
  AllPredictionResults_6reps=list.files(pattern="^Correlation_Correctedphenotype_GBLUP_GPD_6Rep")
  
  dataframeWithResults_6reps=rep(0,6)
  dataframeWithResults_6reps=as.data.frame(dataframeWithResults_6reps)
  dataframeWithResults_6reps=t(dataframeWithResults_6reps)
  colnames(dataframeWithResults_6reps)=c("F1.population","Mean.Parental.predicted.GEBV","F1.pop.mean.dryweight","F1.pop.mean.freshweight","Round","Replicates")
  
  for (i in seq(1:length(AllPredictionResults_6reps))){
    filename=AllPredictionResults_6reps[i]
    round=str_sub(filename, 53, 54)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,53,55)=="100"){
      round=str_sub(filename,53,55)
    }
    
    table=read.table(AllPredictionResults_6reps[i],sep="\t",head=T)
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    table_df$Replicates=as.numeric(6)
    
    notused=c("DLF2 NR high yield","DLF3 NR high yield","DLF4 NR high yield","DLF5 NR high yield","YellowTip1 SM yellow tip")
    rmv=which(table_df$F1.population %in% notused)  # rows to delete
    table_df=table_df[-rmv,]
    
    dataframeWithResults_6reps=rbind(dataframeWithResults_6reps,table_df)
  }
  
  dataframeWithResults_6reps=dataframeWithResults_6reps[-1,]
  dataframeWithResults_6reps
  
}

# Load files 5 reps
{
  AllPredictionResults_5reps=list.files(pattern="^Correlation_Correctedphenotype_GBLUP_GPD_5Rep")
  
  dataframeWithResults_5reps=rep(0,6)
  dataframeWithResults_5reps=as.data.frame(dataframeWithResults_5reps)
  dataframeWithResults_5reps=t(dataframeWithResults_5reps)
  colnames(dataframeWithResults_5reps)=c("F1.population","Mean.Parental.predicted.GEBV","F1.pop.mean.dryweight","F1.pop.mean.freshweight","Round","Replicates")
  
  for (i in seq(1:length(AllPredictionResults_5reps))){
    filename=AllPredictionResults_5reps[i]
    round=str_sub(filename, 53, 54)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,53,55)=="100"){
      round=str_sub(filename,53,55)
    }
    
    table=read.table(AllPredictionResults_5reps[i],sep="\t",head=T)
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    table_df$Replicates=as.numeric(5)
    
    notused=c("DLF2 NR high yield","DLF3 NR high yield","DLF4 NR high yield","DLF5 NR high yield","YellowTip1 SM yellow tip")
    rmv=which(table_df$F1.population %in% notused)  # rows to delete
    table_df=table_df[-rmv,]
    
    dataframeWithResults_5reps=rbind(dataframeWithResults_5reps,table_df)
  }
  
  dataframeWithResults_5reps=dataframeWithResults_5reps[-1,]
  dataframeWithResults_5reps
  
}

# Load files 4 reps
{
  AllPredictionResults_4reps=list.files(pattern="^Correlation_Correctedphenotype_GBLUP_GPD_4Rep")
  
  dataframeWithResults_4reps=rep(0,6)
  dataframeWithResults_4reps=as.data.frame(dataframeWithResults_4reps)
  dataframeWithResults_4reps=t(dataframeWithResults_4reps)
  colnames(dataframeWithResults_4reps)=c("F1.population","Mean.Parental.predicted.GEBV","F1.pop.mean.dryweight","F1.pop.mean.freshweight","Round","Replicates")
  
  for (i in seq(1:length(AllPredictionResults_4reps))){
    filename=AllPredictionResults_4reps[i]
    round=str_sub(filename, 53, 54)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,53,55)=="100"){
      round=str_sub(filename,53,55)
    }
    
    table=read.table(AllPredictionResults_4reps[i],sep="\t",head=T)
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    table_df$Replicates=as.numeric(4)
    
    notused=c("DLF2 NR high yield","DLF3 NR high yield","DLF4 NR high yield","DLF5 NR high yield","YellowTip1 SM yellow tip")
    rmv=which(table_df$F1.population %in% notused)  # rows to delete
    table_df=table_df[-rmv,]
    
    dataframeWithResults_4reps=rbind(dataframeWithResults_4reps,table_df)
  }
  
  dataframeWithResults_4reps=dataframeWithResults_4reps[-1,]
  dataframeWithResults_4reps
  
}


# Load files 3 reps
{
  AllPredictionResults_3reps=list.files(pattern="^Correlation_Correctedphenotype_GBLUP_GPD_3Rep")
  
  dataframeWithResults_3reps=rep(0,6)
  dataframeWithResults_3reps=as.data.frame(dataframeWithResults_3reps)
  dataframeWithResults_3reps=t(dataframeWithResults_3reps)
  colnames(dataframeWithResults_3reps)=c("F1.population","Mean.Parental.predicted.GEBV","F1.pop.mean.dryweight","F1.pop.mean.freshweight","Round","Replicates")
  
  for (i in seq(1:length(AllPredictionResults_3reps))){
    filename=AllPredictionResults_3reps[i]
    round=str_sub(filename, 53, 54)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,53,55)=="100"){
      round=str_sub(filename,53,55)
    }
    table=read.table(AllPredictionResults_3reps[i],sep="\t",head=T)
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    table_df$Replicates=as.numeric(3)
    
    notused=c("DLF2 NR high yield","DLF3 NR high yield","DLF4 NR high yield","DLF5 NR high yield","YellowTip1 SM yellow tip")
    rmv=which(table_df$F1.population %in% notused)  # rows to delete
    table_df=table_df[-rmv,]
    
    dataframeWithResults_3reps=rbind(dataframeWithResults_3reps,table_df)
  }
  
  dataframeWithResults_3reps=dataframeWithResults_3reps[-1,]
  dataframeWithResults_3reps
  
}

# Load files 2 reps
{
  AllPredictionResults_2reps=list.files(pattern="^Correlation_Correctedphenotype_GBLUP_GPD_2Rep")
  
  dataframeWithResults_2reps=rep(0,6)
  dataframeWithResults_2reps=as.data.frame(dataframeWithResults_2reps)
  dataframeWithResults_2reps=t(dataframeWithResults_2reps)
  colnames(dataframeWithResults_2reps)=c("F1.population","Mean.Parental.predicted.GEBV","F1.pop.mean.dryweight","F1.pop.mean.freshweight","Round","Replicates")
  
  for (i in seq(1:length(AllPredictionResults_2reps))){
    filename=AllPredictionResults_2reps[i]
    round=str_sub(filename, 53, 54)
    if (str_sub(round,2,2)=="."){
      round=str_sub(round,1,1)
    }
    else if (str_sub(filename,53,55)=="100"){
      round=str_sub(filename,53,55)
    }
    table=read.table(AllPredictionResults_2reps[i],sep="\t",head=T)
    table_df=as.data.frame(table)
    table_df$Round=as.numeric(round)
    table_df$Replicates=as.numeric(2)
    
    notused=c("DLF2 NR high yield","DLF3 NR high yield","DLF4 NR high yield","DLF5 NR high yield","YellowTip1 SM yellow tip")
    rmv=which(table_df$F1.population %in% notused)  # rows to delete
    table_df=table_df[-rmv,]
    
    dataframeWithResults_2reps=rbind(dataframeWithResults_2reps,table_df)
  }
  
  dataframeWithResults_2reps=dataframeWithResults_2reps[-1,]
  dataframeWithResults_2reps
}
  
# save results in output file    
{
MergedData=rbind(dataframeWithResults_2reps,dataframeWithResults_3reps,dataframeWithResults_4reps,dataframeWithResults_5reps,dataframeWithResults_6reps,dataframeWithResults_7reps,dataframeWithResults_8reps,dataframeWithResults_9reps,dataframeWithResults_10reps)

# check
nrow(MergedData)== 100*9*9

write.table(MergedData,file="/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/F1_pred_gpd_ResCor/replicate_reduction_F1Results.txt",sep="\t",quote=F,col.names = T, row.names = F)
}




# Make data for plotting

numberofrounds=nrow(dataframeWithResults_2reps)/9


ForPlotting_10REPs = data.frame(matrix(vector(), numberofrounds, 2,
                                       dimnames=list(c(), c("correlation","Round"))),
                                stringsAsFactors=F)

for (j in unique(dataframeWithResults_10reps$Round)){
  print(j)
  subtable=dataframeWithResults_10reps[which(dataframeWithResults_10reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_10REPs$correlation[j] = cor(subtable$Mean.Parental.predicted.GEBV,subtable$F1.pop.mean.dryweight)
  ForPlotting_10REPs$Round[j] = subtable$Round[1]
}

ForPlotting_10REPs$Replicates=rep(10,numberofrounds)



# Make data for plotting
ForPlotting_9REPs = data.frame(matrix(vector(), numberofrounds, 2,
                                       dimnames=list(c(), c("correlation","Round"))),
                                stringsAsFactors=F)

for (j in unique(dataframeWithResults_9reps$Round)){
  print(j)
  subtable=dataframeWithResults_9reps[which(dataframeWithResults_9reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_9REPs$correlation[j] = cor(subtable$Mean.Parental.predicted.GEBV,subtable$F1.pop.mean.dryweight)
  ForPlotting_9REPs$Round[j] = subtable$Round[1]
}

ForPlotting_9REPs$Replicates=rep(9,numberofrounds)

# Make data for plotting
ForPlotting_8REPs = data.frame(matrix(vector(), numberofrounds, 2,
                                      dimnames=list(c(), c("correlation","Round"))),
                               stringsAsFactors=F)

for (j in unique(dataframeWithResults_8reps$Round)){
  print(j)
  subtable=dataframeWithResults_8reps[which(dataframeWithResults_8reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_8REPs$correlation[j] = cor(subtable$Mean.Parental.predicted.GEBV,subtable$F1.pop.mean.dryweight)
  ForPlotting_8REPs$Round[j] = subtable$Round[1]
}

ForPlotting_8REPs$Replicates=rep(8,numberofrounds)



# Make data for plotting
ForPlotting_7REPs = data.frame(matrix(vector(), numberofrounds, 2,
                                      dimnames=list(c(), c("correlation","Round"))),
                               stringsAsFactors=F)

for (j in unique(dataframeWithResults_7reps$Round)){
  print(j)
  subtable=dataframeWithResults_7reps[which(dataframeWithResults_7reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_7REPs$correlation[j] = cor(subtable$Mean.Parental.predicted.GEBV,subtable$F1.pop.mean.dryweight)
  ForPlotting_7REPs$Round[j] = subtable$Round[1]
}

ForPlotting_7REPs$Replicates=rep(7,numberofrounds)




# Make data for plotting
ForPlotting_6REPs = data.frame(matrix(vector(), numberofrounds, 2,
                                      dimnames=list(c(), c("correlation","Round"))),
                               stringsAsFactors=F)

for (j in unique(dataframeWithResults_6reps$Round)){
  print(j)
  subtable=dataframeWithResults_6reps[which(dataframeWithResults_6reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_6REPs$correlation[j] = cor(subtable$Mean.Parental.predicted.GEBV,subtable$F1.pop.mean.dryweight)
  ForPlotting_6REPs$Round[j] = subtable$Round[1]
}

ForPlotting_6REPs$Replicates=rep(6,numberofrounds)


# Make data for plotting
ForPlotting_5REPs = data.frame(matrix(vector(), numberofrounds, 2,
                                      dimnames=list(c(), c("correlation","Round"))),
                               stringsAsFactors=F)

for (j in unique(dataframeWithResults_5reps$Round)){
  print(j)
  subtable=dataframeWithResults_5reps[which(dataframeWithResults_5reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_5REPs$correlation[j] = cor(subtable$Mean.Parental.predicted.GEBV,subtable$F1.pop.mean.dryweight)
  ForPlotting_5REPs$Round[j] = subtable$Round[1]
}

ForPlotting_5REPs$Replicates=rep(5,numberofrounds)


# Make data for plotting
ForPlotting_4REPs = data.frame(matrix(vector(), numberofrounds, 2,
                                      dimnames=list(c(), c("correlation","Round"))),
                               stringsAsFactors=F)

for (j in unique(dataframeWithResults_4reps$Round)){
  print(j)
  subtable=dataframeWithResults_4reps[which(dataframeWithResults_4reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_4REPs$correlation[j] = cor(subtable$Mean.Parental.predicted.GEBV,subtable$F1.pop.mean.dryweight)
  ForPlotting_4REPs$Round[j] = subtable$Round[1]
}

ForPlotting_4REPs$Replicates=rep(4,numberofrounds)

# Make data for plotting
ForPlotting_3REPs = data.frame(matrix(vector(), numberofrounds, 2,
                                      dimnames=list(c(), c("correlation","Round"))),
                               stringsAsFactors=F)

for (j in unique(dataframeWithResults_3reps$Round)){
  print(j)
  subtable=dataframeWithResults_3reps[which(dataframeWithResults_3reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_3REPs$correlation[j] = cor(subtable$Mean.Parental.predicted.GEBV,subtable$F1.pop.mean.dryweight)
  ForPlotting_3REPs$Round[j] = subtable$Round[1]
}

ForPlotting_3REPs$Replicates=rep(3,numberofrounds)

# Make data for plotting
ForPlotting_2REPs = data.frame(matrix(vector(), numberofrounds, 2,
                                      dimnames=list(c(), c("correlation","Round"))),
                               stringsAsFactors=F)

for (j in unique(dataframeWithResults_2reps$Round)){
  print(j)
  subtable=dataframeWithResults_2reps[which(dataframeWithResults_2reps$Round==j),]
  j=as.numeric(j)
  ForPlotting_2REPs$correlation[j] = cor(subtable$Mean.Parental.predicted.GEBV,subtable$F1.pop.mean.dryweight)
  ForPlotting_2REPs$Round[j] = subtable$Round[1]
}

ForPlotting_2REPs$Replicates=rep(2,numberofrounds)


# Bind all data 

AllDataForplotting=rbind(ForPlotting_2REPs,ForPlotting_3REPs,ForPlotting_4REPs,ForPlotting_5REPs,ForPlotting_6REPs,ForPlotting_7REPs,ForPlotting_8REPs,ForPlotting_9REPs,ForPlotting_10REPs)

# check
nrow(AllDataForplotting)== numberofrounds*9


write.table(AllDataForplotting,file="/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/F1_pred_gpd_ResCor/replicate_reduction_F1summary.txt",sep="\t",quote=F,col.names = T, row.names = F)





# plotting



#colnames(data)=c("replicates","cor","sterrr","lower","upper")

pd <- position_dodge(0.1) # move them .05 to the left and right

p1=ggplot(AllDataForplotting, aes(x=Replicates, y=correlation)) + 
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
  #geom_line(position=pd) +
  geom_point(position=pd) +
  ylim(min(AllDataForplotting$correlation),max(AllDataForplotting$correlation)) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Correlation") +
  theme_classic()

p1

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/cor_of_gpdResCorAndGEBV_F1prediction.pdf', plot = p1, width = 15, height = 15, unit = 'cm')


# Correlation between gpd_ResCor and GEBVs, means with standard deviations
Means_cor_gpdResCor_GEBV=aggregate(AllDataForplotting$correlation, list(AllDataForplotting$Replicates), mean)
Means_cor_gpdResCor_GEBV=as.data.frame(Means_cor_gpdResCor_GEBV)
colnames(Means_cor_gpdResCor_GEBV)=c("Replicates","Means_cor_gpdResCor_GEBV")
STDEV_cor_gpdResCor_GEBV=aggregate(AllDataForplotting$correlation, list(AllDataForplotting$Replicates), sd)
Means_cor_gpdResCor_GEBV$lower=Means_cor_gpdResCor_GEBV$Means_cor_gpdResCor_GEBV-STDEV_cor_gpdResCor_GEBV[,2]
Means_cor_gpdResCor_GEBV$upper=Means_cor_gpdResCor_GEBV$Means_cor_gpdResCor_GEBV+STDEV_cor_gpdResCor_GEBV[,2]
colnames(Means_cor_gpdResCor_GEBV)[3:4]=c("lower","upper")

p3=ggplot(Means_cor_gpdResCor_GEBV, aes(x=Replicates, y=Means_cor_gpdResCor_GEBV)) + 
  geom_line(position=pd) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.05, position=pd) +
  #ylim(0,0.4) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Correlation") +
  theme_classic()

p3

ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/cor_of_gpdResCorAndGEBV_F1_averages_withSD.pdf', plot = p3, width = 15, height = 15, unit = 'cm')


# Correlation between gpd_ResCor and GEBVs, means with standard errors
Means_cor_gpdResCor_GEBV_SE=aggregate(AllDataForplotting$correlation, list(AllDataForplotting$Replicates), mean)
Means_cor_gpdResCor_GEBV_SE=as.data.frame(Means_cor_gpdResCor_GEBV_SE)
colnames(Means_cor_gpdResCor_GEBV_SE)=c("Replicates","Means_cor_gpdResCor_GEBV")
STDEV_cor_gpdResCor_GEBV=aggregate(AllDataForplotting$correlation, list(AllDataForplotting$Replicates), sd)
STERR_cor_gpdResCor_GEBV=STDEV_cor_gpdResCor_GEBV/sqrt(10)

Means_cor_gpdResCor_GEBV_SE$lower=Means_cor_gpdResCor_GEBV_SE$Means_cor_gpdResCor_GEBV-STERR_cor_gpdResCor_GEBV[,2]
Means_cor_gpdResCor_GEBV_SE$upper=Means_cor_gpdResCor_GEBV_SE$Means_cor_gpdResCor_GEBV+STERR_cor_gpdResCor_GEBV[,2]
colnames(Means_cor_gpdResCor_GEBV_SE)[3:4]=c("lower","upper")

p4=ggplot(Means_cor_gpdResCor_GEBV_SE, aes(x=Replicates, y=Means_cor_gpdResCor_GEBV)) + 
  geom_line(position=pd) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.05, position=pd) +
  #ylim(0,0.4) +
  scale_x_reverse() +
  xlab("Replicates") +
  ylab("Correlation") +
  theme_classic()

p4
ggsave('/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Replicate_Reduction_20200831/Figures/cor_of_gpdResCorAndGEBV_averagesF1_withSE.pdf', plot = p4, width = 15, height = 15, unit = 'cm')




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
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),1], AllDataForplotting[which(AllDataForplotting$Replicates==9),1], paired = FALSE, alternative = "greater") # N.S
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),1], AllDataForplotting[which(AllDataForplotting$Replicates==8),1], paired = FALSE, alternative = "greater") # *
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),1], AllDataForplotting[which(AllDataForplotting$Replicates==7),1], paired = FALSE, alternative = "greater") #***
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),1], AllDataForplotting[which(AllDataForplotting$Replicates==6),1], paired = FALSE, alternative = "greater") #***
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),1], AllDataForplotting[which(AllDataForplotting$Replicates==5),1], paired = FALSE, alternative = "greater") #***
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),1], AllDataForplotting[which(AllDataForplotting$Replicates==4),1], paired = FALSE, alternative = "greater") #***
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),1], AllDataForplotting[which(AllDataForplotting$Replicates==3),1], paired = FALSE, alternative = "greater") #***
t.test(AllDataForplotting[which(AllDataForplotting$Replicates==10),1], AllDataForplotting[which(AllDataForplotting$Replicates==2),1], paired = FALSE, alternative = "greater") #***


# In how many cases do the correlation at 9 replicates get lower than at 10 replicates
FDRfunc(10,9)  
FDRfunc(10,8) 
FDRfunc(10,7) 
FDRfunc(10,6) 
FDRfunc(10,5) 
FDRfunc(10,4) 
FDRfunc(10,3) 
FDRfunc(10,2) 
