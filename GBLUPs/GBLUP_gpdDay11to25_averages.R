###################################################################################################################
###################################################################################################################
### this is a script to run GBLUP on replicate data with predefined division into training and testing groups   ###
###################################################################################################################
###################################################################################################################

# Load libraries
{
  library(lme4)
  library(BGLR)
  library("parallel")
  library("methods")
  library("Matrix")
  library("MASS")
}

# Define some variables
{
  args=commandArgs(trailingOnly = TRUE)
  print(args)
  round=args[2]
}

# Define functions for GBLUP
{
GP_GBLUP<-function(testpop){
  
  ################  ################  ################  ################
  ##start by estimating GEBVs for training population individuals 
  ################  ################  ################  ################
  
  
  gpdmeans_training=gpdmeans[-testpop,] # limit the dataframe to only the individuals allowed for training the model
  gpdmeans_training_ready=na.omit(gpdmeans_training, cols = c("gpdDay11to25")) # remember that gpd na inidividuals should be removed whether or not they are in the training pop or not
  
  ind_not_in_train=gpdmeans$Individual[testpop]
  IndividualsToRemoveGRM=which(colnames(GRM1) %in% ind_not_in_train)
  
  GRM_trn = GRM1[-IndividualsToRemoveGRM,-IndividualsToRemoveGRM]
  
  # Run the GBLUP model on full training population to extract GEBVs
  yNA=gpdmeans_training_ready$gpdDay11to25
  ETA=list(list(K=GRM_trn,model="RKHS"))
  GBLUP=BGLR(y=yNA,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("GBLUP",round))
  matrix=cbind(as.character(gpdmeans_training_ready$Individual),as.numeric(gpdmeans_training_ready$gpdDay11to25),as.numeric(GBLUP$ETA[[1]]$u))
  colnames(matrix)=c("ID", "Observed", "GEBV")
  GEBV_contribution1data=as.numeric(as.character(matrix[,3]))
  

  ################  ################ 
  ## Now predict testing population 
  ################  ################  
  
  GRMforpred_test = GRM1[testpop,testpop] # GRM for individuals that will be predicted
  GRMforpred_covar = GRM1[testpop,-testpop] # Covariance between training and testing pop.
  
  #GEBVpred_contr1 = GcloverReps_covar%*%solve(GcloverReps_trn) %*% GEBV_contribution1data 
  GEBVpred = GRMforpred_covar%*%ginv(GRM_trn) %*% GEBV_contribution1data 
  #GEBVpred_contr1 = GcloverReps_covar%*%solve(GcloverReps_trn + diag(0.01, 1661, 1661)) %*% GEBV_contribution1data 
  
  # Output matrix with prediction results
  matrix1=cbind(as.character(gpdmeans$Individual[testpop]),as.numeric(gpdmeans$gpdDay11to25[testpop]),as.numeric(as.character(GEBVpred)))
  colnames(matrix1)=c("ID", "Observed", "GEBV")
  return(matrix1)
}
}
  
# Load data
{
  d <- read.table("GPD_day11to25.csv",head = T)
}


# Remove plants that has been found to not show growth between day 10 and 20 dpi and drop in growth from day 10 dpi to the last day it was measured
{
  remove = read.table("Barcodes_removed_based_on_single_Observations_2021-01-06.txt")
  #remove = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/New_fixation_Trait_20210106/Barcodes_removed_based_on_single_Observations_2021-01-06.txt")
  
  removeidx = which(d$Barcode %in% remove)
  
  if (length(removeidx)==0){
    print("barcodes that showed weird behaviour have already been removed.")
    d005=d
    
  } else{
  d005 = d[-removeidx,]
  }
}


{
  d0=na.omit(d005)
  d0$Clover=as.character(d0$Clover)
  d0$Clover=as.factor(d0$Clover)
  d2=d0
}


# Load genomic relationship matrix and make sure clover genotypes match data
{
  GRM=read.table(args[1],sep=",",header=T)
  #GRM=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Heritabilities_20201217/GRM_Clover_Fullfiltering_20200728.csv",sep=",",header=T)
  
  dim(GRM)
  d2$Clovershort <- strtrim(d2$Clover,8)
  d3=d2[order(d2$Clovershort,decreasing=F),]
  length(unique(d3$Clovershort)) #149
  d4=d3[which(d3$Clovershort %in% colnames(GRM)),]
  length(unique(d4$Clovershort)) #147 unique genotypes with GPD data
  remove=GRM[which(colnames(GRM) %in% d4$Clovershort),which(colnames(GRM) %in% d4$Clovershort)]
  print(remove)
  GRM1=GRM[which(colnames(GRM) %in% d4$Clovershort),which(colnames(GRM) %in% d4$Clovershort)]
  dim(GRM1)
  nrow(GRM1)==length(unique(d4$Clovershort))
  GRM1=data.matrix(GRM1)
  length(colnames(GRM1)==unique(d4$Clovershort))==nrow(GRM1) #check
}

# Aberpearl_07 contribute 700 of the datapoints and thus influence the variance a lot. Cut down Aberpearl_07 data so that we have only 6 different Rhizobia left like the other clovers
{
  Aberpearl_07=which(d4$Clovershort=="Aearl_07")
  Inocolums=unique(d4$Rhizobium[Aberpearl_07])
  set.seed(15)
  sample=sample(Inocolums,6)
  #sample=c("MIX","SM22","SM25","SM149C","SM31","SM155A")
  print(sample)

  which(d4$Rhizobium[Aberpearl_07]==sample[1]) #4
  which(d4$Rhizobium[Aberpearl_07]==sample[2]) #4
  which(d4$Rhizobium[Aberpearl_07]==sample[3]) #4
  which(d4$Rhizobium[Aberpearl_07]==sample[4]) #4
  which(d4$Rhizobium[Aberpearl_07]==sample[5]) #4
  which(d4$Rhizobium[Aberpearl_07]==sample[6]) #4

  remove=which((d4$Rhizobium[Aberpearl_07] %in% sample)==FALSE)
  d4=d4[-Aberpearl_07[remove],]
  nrow(d4)
}

{
  d4$roundRep <- paste(d4$Round, d4$Replicate, sep='_')
  d6=d4
  nrow(d6)
  length(which(d6$Clovershort=="Aearl_07"))
}

# Clean up
{
  d6$Rhizobium=droplevels(d6$Rhizobium) # removing levels not used in actual data
  d6$Clover=droplevels(d6$Clover) # removing levels not used in actual data
  d6=d6[order(d6$Clovershort),] # make sure it is in alphabetic order like the GRM
}


# calculate means of genotypes
{
  gpdmeans=aggregate(d6$GPD_in_interval, list(d6$Clovershort), mean)
  colnames(gpdmeans)=c("Individual","gpdDay11to25")
}


# Setting up 6-fold CV system
{
# if groups are not provided 
# Make six groups for 6-fold CV

if (is.na(args[3])){
  set.seed(NULL)
  tst=sample(1:length(unique(gpdmeans$Individual)),size=length(unique(gpdmeans$Individual)),replace=FALSE) 
  k=6
  testing_pop=split(tst, sort(tst%%k))
  
  tst1=testing_pop[1]$'0'
  tst2=testing_pop[2]$'1'
  tst3=testing_pop[3]$'2'
  tst4=testing_pop[4]$'3'
  tst5=testing_pop[5]$'4'
  tst6=testing_pop[6]$'5'
  
  testpop1=unique(gpdmeans$Individual)[tst1]
  testpop2=unique(gpdmeans$Individual)[tst2]
  testpop3=unique(gpdmeans$Individual)[tst3]
  testpop4=unique(gpdmeans$Individual)[tst4]
  testpop5=unique(gpdmeans$Individual)[tst5]
  testpop6=unique(gpdmeans$Individual)[tst6]
  
  testpop1_idx=which(gpdmeans$Individual %in% testpop1)
  testpop2_idx=which(gpdmeans$Individual %in% testpop2)
  testpop3_idx=which(gpdmeans$Individual %in% testpop3)
  testpop4_idx=which(gpdmeans$Individual %in% testpop4)
  testpop5_idx=which(gpdmeans$Individual %in% testpop5)
  testpop6_idx=which(gpdmeans$Individual %in% testpop6)
  
  grouping=list(testpop1,testpop2,testpop3,testpop4,testpop5,testpop6)
  sink(paste("grouping",round,".txt",sep=""))
  print(grouping)
  sink()
  
  tests=list(testpop1_idx,testpop2_idx,testpop3_idx,testpop4_idx,testpop5_idx,testpop6_idx)
  
} else {
  set.seed(NULL)
  #Load groups from file
  f=read.table(args[3],fill = TRUE)
  linewherenewgroupcomes=vector()
  
  for (i in 1:nrow(f)){
    beginning=as.character(f[i,1])
    secondchr=strsplit(beginning,"")[[1]][2]
    
    if (secondchr=="["){
      linewherenewgroupcomes=append(linewherenewgroupcomes,i)
    }
  }
  
  group1=f[(linewherenewgroupcomes[1]+1):(linewherenewgroupcomes[1+1]-1),2:ncol(f)]
  testpop1=c(unique(as.character(group1[,1])),unique(as.character(group1[,2])),unique(as.character(group1[,3])),unique(as.character(group1[,4])),unique(as.character(group1[,5])),unique(as.character(group1[,6])))
  
  group2=f[(linewherenewgroupcomes[2]+1):(linewherenewgroupcomes[2+1]-1),2:ncol(f)]
  testpop2=c(unique(as.character(group2[,1])),unique(as.character(group2[,2])),unique(as.character(group2[,3])),unique(as.character(group2[,4])),unique(as.character(group2[,5])),unique(as.character(group2[,6])))
  
  group3=f[(linewherenewgroupcomes[3]+1):(linewherenewgroupcomes[3+1]-1),2:ncol(f)]
  testpop3=c(unique(as.character(group3[,1])),unique(as.character(group3[,2])),unique(as.character(group3[,3])),unique(as.character(group3[,4])),unique(as.character(group3[,5])),unique(as.character(group3[,6])))
  
  group4=f[(linewherenewgroupcomes[4]+1):(linewherenewgroupcomes[4+1]-1),2:ncol(f)]
  testpop4=c(unique(as.character(group4[,1])),unique(as.character(group4[,2])),unique(as.character(group4[,3])),unique(as.character(group4[,4])),unique(as.character(group4[,5])),unique(as.character(group4[,6])))
  
  group5=f[(linewherenewgroupcomes[5]+1):(linewherenewgroupcomes[5+1]-1),2:ncol(f)]
  testpop5=c(unique(as.character(group5[,1])),unique(as.character(group5[,2])),unique(as.character(group5[,3])),unique(as.character(group5[,4])),unique(as.character(group5[,5])),unique(as.character(group5[,6])))
  
  group6=f[(linewherenewgroupcomes[6]+1):nrow(f),2:ncol(f)]
  testpop6=c(unique(as.character(group6[,1])),unique(as.character(group6[,2])),unique(as.character(group6[,3])),unique(as.character(group6[,4])),unique(as.character(group6[,5])),unique(as.character(group6[,6])))
  
  testpop1_idx=which(gpdmeans$Individual %in% testpop1)
  testpop2_idx=which(gpdmeans$Individual %in% testpop2)
  testpop3_idx=which(gpdmeans$Individual %in% testpop3)
  testpop4_idx=which(gpdmeans$Individual %in% testpop4)
  testpop5_idx=which(gpdmeans$Individual %in% testpop5)
  testpop6_idx=which(gpdmeans$Individual %in% testpop6)
  
  tests=list(testpop1_idx,testpop2_idx,testpop3_idx,testpop4_idx,testpop5_idx,testpop6_idx)
}
}

#  Calculate heritability
{
  print("Calculating heritability")
  y=gpdmeans[,2]
  
  ETA=list(list(K=GRM1,model="RKHS")) 
  GBLUP=BGLR(y=y,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("heritabilityestimation",round))
  
  print(GBLUP$ETA[[1]]$varU)
  print(GBLUP$varE)
  h2=(GBLUP$ETA[[1]]$varU)/(GBLUP$ETA[[1]]$varU+GBLUP$varE)
  write.table(h2,paste("h2",round,".txt",sep=""),sep="\t")
  
  
  #Cloverfit2 <- bayz(gpdDay11to25 ~ ranf(Individual)  + ranf(Individual,V=GRM1),
   #                  data = gpdmeans, chain=c(20000, 5000, 10))
  
  #summary(Cloverfit2) 
  #head(Cloverfit2$Samples)
  #plot(Cloverfit2$Samples[,3]) # var.ranf.Individual.GRM1
  #plot(Cloverfit2$Samples[,2]) # var.ranf.Individual
  #plot(Cloverfit2$Samples[,1]) # var.resid
  #plot(Cloverfit2$Samples[,4]) # mean
  
  
  
}  

# Prediction
{
  print("Starting GBLUP prediction")
  results=mclapply(tests,GP_GBLUP)  
}  

# Summarize and make files with results
{
  first=results[[1]]
  second=results[[2]]
  third=results[[3]]
  fourth=results[[4]]
  fifth=results[[5]]
  sixth=results[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)
  correlation=cor(as.numeric(as.character(All[,2])),as.numeric(as.character(All[,3]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_GPD",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_GPD",round,".txt",sep="")
  write.table(All,filename1,sep="\t",quote=F,row.names=F)
}  