###################################################################################################################
###################################################################################################################
### this is a script to run GBLUP on replicate data removing one replicate pr. genotype each round              ###
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

# Load data
{
  d <- read.csv("/home/cks/NChain/faststorage/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/greenhouse_area.csv", header = TRUE, sep = ",")
  f=read.csv("/home/cks/NChain/faststorage/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/2018_weight.csv",header=T,sep=";")
  colnames(f)[1]="Barcode"
  df=merge(d,f,by="Barcode")
  d=df
}

# Calculate growth pr. day
{
  d$days_of_growth <- as.Date(d$harvest_date, format = "%d/%m/%y") - as.Date(d$inoculation_date, format = "%d/%m/%y")
  d$growth_per_day <- d$weight.y/as.numeric(d$days_of_growth)
}

# Plants that has been found to not show growth between day 10 and 20 dpi according to "Heritability_of_growth_periods_20191113.R" script
{
  DeadplantBarcodes = c("1870","1976","2007","2057","2235","2277","2371","2449","2534","2716","2747","3325","3424","3465","3504","3869","3880","4704","4706","4902")
  DeadplantBarcodesidx = which(d$Barcode %in% DeadplantBarcodes)
  d05 = d[-DeadplantBarcodesidx,]
}

# Plants that has been found to drop in growth from day 10 dpi to the last day it was measured according to "Heritability_of_growth_periods_20191113.R" script
{
  WeirdMeasurementBarcodes = c("1821","3663","5135")
  WeirdMeasurementBarcodesBarcodesidx = which(d05$Barcode %in% WeirdMeasurementBarcodes)
  d005 = d05[-WeirdMeasurementBarcodesBarcodesidx,]
}

# Sort out plants that were inoculated with no rhizobium, SM73 or had n_stolon=0. These are errors and don't grow
{
  
  d0=na.omit(d005)
  d0$Clover=as.character(d0$Clover)
  d0$Clover[which(d0$Clover=="AAran_0104")]="Aaran_0104"
  d0$Clover[which(d0$Clover=="AAran_0206")]="Aaran_0206"
  d0$Clover=as.factor(d0$Clover)
  
  d0=d0[-which(d0$rhizobium=="SM73"),]
  d0=d0[-which(d0$rhizobium=="NO"),]
  d0=d0[-which(d0$n_stolons==0),]
  d0=d0[-which(d0$n_stolons==1),]
  d0=d0[-which(d0$n_stolons==2),]
  d0=d0[-which(d0$n_stolons==3),]
  
  d2=d0
}


# Load genomic relationship matrix and make sure clover genotypes match data
{
  GRM=read.table(args[1],sep=",",header=T)
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

# Make a column dedicated to RoundRep
{
  d4$roundRep <- paste(d4$Round, d4$Replicate, sep='_')
  d6=d4
  nrow(d6)
  length(which(d6$Clovershort=="Aearl_07"))
}

# Remove all genotypes that has <10 replicates
Genotypes=(unique(d6$Clovershort))

for (genotype in Genotypes){
  idx=which(d6$Clovershort==genotype)
  if (length(idx)<10){
    d6=d6[-idx,]
    print(paste(genotype,"removed",sep=" "))
  }
}

# Remove genotypes not included from GRM
{
  GRM1=GRM1[which(colnames(GRM1) %in% d6$Clovershort),which(colnames(GRM1) %in% d6$Clovershort)]
  nrow(GRM1)==length(unique(d6$Clovershort))
  GRM1=data.matrix(GRM1)
}

# Clean up
{
  d6$Rhizobium=droplevels(d6$Rhizobium) # removing levels not used in actual data
  d6$Clover=droplevels(d6$Clover) # removing levels not used in actual data
  d6=d6[order(d6$Clovershort),] # make sure it is in alphabetic order like the GRM
}

# Divide into 6 populations for GP
if (is.na(args[3])){
  set.seed(NULL)
  tst=sample(1:length(unique(d6$Clovershort)),size=length(unique(d6$Clovershort)),replace=FALSE) 
  k=6
  testing_pop=split(tst, sort(tst%%k))
  
  tst1=testing_pop[1]$'0'
  tst2=testing_pop[2]$'1'
  tst3=testing_pop[3]$'2'
  tst4=testing_pop[4]$'3'
  tst5=testing_pop[5]$'4'
  tst6=testing_pop[6]$'5'
  
  testpop1=unique(d6$Clovershort)[tst1]
  testpop2=unique(d6$Clovershort)[tst2]
  testpop3=unique(d6$Clovershort)[tst3]
  testpop4=unique(d6$Clovershort)[tst4]
  testpop5=unique(d6$Clovershort)[tst5]
  testpop6=unique(d6$Clovershort)[tst6]
  
  testpop1_idx=which(d6$Clovershort %in% testpop1)
  testpop2_idx=which(d6$Clovershort %in% testpop2)
  testpop3_idx=which(d6$Clovershort %in% testpop3)
  testpop4_idx=which(d6$Clovershort %in% testpop4)
  testpop5_idx=which(d6$Clovershort %in% testpop5)
  testpop6_idx=which(d6$Clovershort %in% testpop6)
  
  grouping=list(testpop1,testpop2,testpop3,testpop4,testpop5,testpop6)
  sink(paste("grouping",round,".txt",sep=""))
  print(grouping)
  sink()
  
  tests=list(testpop1_idx,testpop2_idx,testpop3_idx,testpop4_idx,testpop5_idx,testpop6_idx)
  
} else {
  #Load groups from file
  f=read.table(args[3],fill = TRUE)
  linewherenewgroupcomes=vector()
  set.seed(NULL)
  
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
  
}


############################################################

# Now remove replicates so each genotype has a maximum of of desired number (maxreplicates) 

removereplicates <- function(maxreplicates,dataframe){
  for (genotype in Genotypes){
    replicateidx=which(dataframe$Clovershort==genotype)
    if (length(replicateidx)>maxreplicates){
      numbertoremove=length(replicateidx)-maxreplicates
      remove=sample(replicateidx,numbertoremove)
      dataframe=dataframe[-remove,]
    }
  }
  return(dataframe)
}



Correcterfunction<-function(dataframe){
  
  # Function that takes in a dataframe and calculates the gpd_rescor and corrected phenotype
  # Correct GPD for the full effect of initial size 
  # Correction should only be based on the data that are left when reducing to the desired number
  # of replicates
  {
    lm.fit <- lmer(growth_per_day ~ InitialSize + (1|Clover), data=dataframe) 
    ycorr <- dataframe$growth_per_day - model.matrix( ~ InitialSize, data=dataframe) %*% fixef(lm.fit)
    dataframe$gpd_dryweight_cor <- ycorr #this is the new corrected dry weight
    print(paste("this is the correlation coefficient between gpd_rescor and iSize:",round(cor(dataframe$gpd_dryweight_cor,dataframe$InitialSize),2)),sep=" ")
    print(paste("this is the correlation coefficient between gpd_NoCor and iSize:",round(cor(dataframe$growth_per_day,dataframe$InitialSize),2)), sep=" ")
  }
  
  # Calculate a phenotype corrected for all fixed effects
  {
    Correctedforallfixed <- lmer(gpd_dryweight_cor ~ factor(Round) + factor(NS) + factor(EW) + factor(Rhizobium) + inoculation_date + (1|Clover), data=dataframe) 
    summary(Correctedforallfixed)
    matrixOfeffects=model.matrix( ~ factor(Round) + factor(NS) + factor(EW)  + factor(Rhizobium) + inoculation_date, data=dataframe)
    column_to_remove=which(colnames(matrixOfeffects) %in% names(fixef(Correctedforallfixed))==F)
    matrixOfeffects_new=matrixOfeffects[,-column_to_remove]
    
    ycorr <- dataframe$gpd_dryweight_cor - matrixOfeffects_new %*% fixef(Correctedforallfixed)
    dataframe$CorrectedPheno <- ycorr 
  }  
  
  return(dataframe)
}

Make_GRM_matrices<-function(dataframe){
    CloverDesign <- model.matrix(~0+dataframe$Clovershort)
    GcloverReps <- CloverDesign %*% GRM1 %*% t(CloverDesign) 
    CloverIndep <- CloverDesign %*% t(CloverDesign)
    
    out <- list(GcloverReps, CloverIndep)    
    
    return(out)
}

testpop_generator<-function(dataframe){
  #Find indexes for test population  
  testpop1_idx=which(dataframe$Clovershort %in% testpop1)
  testpop2_idx=which(dataframe$Clovershort %in% testpop2)
  testpop3_idx=which(dataframe$Clovershort %in% testpop3)
  testpop4_idx=which(dataframe$Clovershort %in% testpop4)
  testpop5_idx=which(dataframe$Clovershort %in% testpop5)
  testpop6_idx=which(dataframe$Clovershort %in% testpop6)
  
  tests=list(testpop1_idx,testpop2_idx,testpop3_idx,testpop4_idx,testpop5_idx,testpop6_idx)
  
  return(tests)
}
  
GP_GBLUP<-function(testpop){
    
  ################  ################  ################  ################
  ##start by estimating GEBVs for training population individuals 
  ################  ################  ################  ################
    
  d6_training=dataframe[-testpop,] # limit the dataframe to only the individuals allowed for training the model
  d6_training_ready=na.omit(d6_training, cols = c("gpd_dryweight_cor")) # remember that gpd na inidividuals should be removed whether or not they are in the training pop or not
    
  ind_not_in_train=dataframe$Clovershort[testpop]
  IndividualsToRemoveGRM=which(colnames(GRM1) %in% ind_not_in_train)
    
  GRM_trn = GRM1[-IndividualsToRemoveGRM,-IndividualsToRemoveGRM]
    
  GcloverReps_trn = GcloverReps[-testpop,-testpop] # a GRM for the individuals used to predict
  CloverIndep_trn= CloverIndep[-testpop,-testpop] # another GRM for the individuals used to predict
    
  # Run the GBLUP model on full training population to extract GEBVs
  yNA=d6_training_ready$gpd_dryweight_cor
    
  fixedmod=model.matrix(~factor(d6_training_ready$Round)+factor(d6_training_ready$NS)+factor(d6_training_ready$EW)+factor(d6_training_ready$Rhizobium)+factor(d6_training_ready$inoculation_date))
  ETA=list(list(K=GcloverReps_trn,model="RKHS"),list(K=CloverIndep_trn,model="RKHS"),list(X=fixedmod,model="FIXED"))
  GBLUP=BGLR(y=yNA,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("GBLUP",round))
  matrix=cbind(as.character(d6_training_ready$Clovershort),as.numeric(as.character(d6_training_ready$gpd_dryweight_cor)),as.numeric(d6_training_ready$CorrectedPheno),as.numeric(GBLUP$ETA[[1]]$u),as.numeric(GBLUP$ETA[[2]]$u),(as.numeric(GBLUP$ETA[[1]]$u)+as.numeric(GBLUP$ETA[[2]]$u)))
  colnames(matrix)=c("ID", "gpd_ResCor","Observed", "1st GEBV contribution","2nd GEBV contribution (clover independent)","totalGEBV")
  GEBV_contribution1data=as.numeric(as.character(matrix[,4]))
  GEBV_contribution1data=as.matrix(GEBV_contribution1data)
  GEBV_contribution2data=as.numeric(as.character(matrix[,5]))
    
  ################  ################ 
  ## Now predict testing population 
  ################  ################  
    
  GcloverReps_test = GcloverReps[testpop,testpop] # GRM for individuals that will be predicted
  GcloverReps_covar = GcloverReps[testpop,-testpop] # Covariance between training and testing pop.
    
  #GEBVpred_contr1 = GcloverReps_covar%*%solve(GcloverReps_trn) %*% GEBV_contribution1data 
  GEBVpred = GcloverReps_covar%*%ginv(GcloverReps_trn) %*% GEBV_contribution1data 
  #GEBVpred_contr1 = GcloverReps_covar%*%solve(GcloverReps_trn + diag(0.01, 1661, 1661)) %*% GEBV_contribution1data 
    
  # Output matrix with prediction results
  matrix1=cbind(as.character(dataframe$Clovershort[testpop]),as.numeric(dataframe$gpd_dryweight_cor[testpop]),as.numeric(dataframe$CorrectedPheno[testpop]),as.numeric(as.character(GEBVpred)))
  colnames(matrix1)=c("ID", "gpd_ResCor","Observed", "GEBV")
  return(matrix1)
  }


#Apply so maximum of replicates is 10
{
  Only10reps=removereplicates(10,d6)
  Only10reps_ready=Correcterfunction(Only10reps)
  if (nrow(Only10reps_ready)==length(unique(d6$Clovershort))*10){
    print("Number of replicates pr. genotype has been reduced to 10")
  GcloverReps=Make_GRM_matrices(Only10reps_ready)[[1]]
  CloverIndep=Make_GRM_matrices(Only10reps_ready)[[2]]
  tests=testpop_generator(Only10reps_ready)
  dataframe = Only10reps_ready
  print("Starting GBLUP prediction")
  results10=mclapply(tests, GP_GBLUP)
  }
}

# Summarize and make files with results
{
  first=results10[[1]]
  second=results10[[2]]
  third=results10[[3]]
  fourth=results10[[4]]
  fifth=results10[[5]]
  sixth=results10[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)
  
  All1=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  All2=aggregate(as.numeric(All[, 4]), list(All[,1]), mean)
  StandardDeviation_gpdRes=aggregate(as.numeric(All[, 2]), list(All[,1]), sd)
  StandardDeviation_corrected=aggregate(as.numeric(All[, 3]), list(All[,1]), sd)
  Gpdmean=aggregate(as.numeric(All[, 2]), list(All[,1]), mean)

  
  correlation1=cor(All1[,2],All2[,2]) #means of replicates
  correlation1

  correlation2=cor(Gpdmean[,2],All2[,2])
  correlation2
  
  Combined=cbind(All1[,1],Gpdmean[,2],All1[,2],All2[,2],StandardDeviation_gpdRes[,2],StandardDeviation_corrected[,2])
  colnames(Combined)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno")
  
  filename=paste("Correlation_Correctedphenotype_GBLUP_GPD_10Replicates",round,".txt",sep="")
  write.table(correlation1,filename,sep="\t",quote=F)
  
  filename=paste("Correlation_gpdResCor_GBLUP_GPD_10Replicates",round,".txt",sep="")
  write.table(correlation2,filename,sep="\t",quote=F)

  filename1=paste("Predictions_GBLUP_GPD_10Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
}  


#Apply so maximum of replicates is 9
{
  Only9reps=removereplicates(9,Only10reps)
  Only9reps_ready=Correcterfunction(Only9reps)
  if (nrow(Only9reps_ready)==length(unique(d6$Clovershort))*9){
    print("Number of replicates pr. genotype has been reduced to 9")
    GcloverReps=Make_GRM_matrices(Only9reps_ready)[[1]]
    CloverIndep=Make_GRM_matrices(Only9reps_ready)[[2]]
    tests=testpop_generator(Only9reps_ready)
    dataframe = Only9reps_ready
    print("Starting GBLUP prediction")
    results9=mclapply(tests, GP_GBLUP)
  }
}

# Summarize and make files with results
{
  first=results9[[1]]
  second=results9[[2]]
  third=results9[[3]]
  fourth=results9[[4]]
  fifth=results9[[5]]
  sixth=results9[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)
  
   
  All1=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  All2=aggregate(as.numeric(All[, 4]), list(All[,1]), mean)
  StandardDeviation_gpdRes=aggregate(as.numeric(All[, 2]), list(All[,1]), sd)
  StandardDeviation_corrected=aggregate(as.numeric(All[, 3]), list(All[,1]), sd)
  Gpdmean=aggregate(as.numeric(All[, 2]), list(All[,1]), mean)

  
  correlation1=cor(All1[,2],All2[,2]) #means of replicates
  correlation1

  correlation2=cor(Gpdmean[,2],All2[,2])
  correlation2
  
  Combined=cbind(All1[,1],Gpdmean[,2],All1[,2],All2[,2],StandardDeviation_gpdRes[,2],StandardDeviation_corrected[,2])
  colnames(Combined)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno")
  
  filename=paste("Correlation_Correctedphenotype_GBLUP_GPD_9Replicates",round,".txt",sep="")
  write.table(correlation1,filename,sep="\t",quote=F)
  
  filename=paste("Correlation_gpdResCor_GBLUP_GPD_9Replicates",round,".txt",sep="")
  write.table(correlation2,filename,sep="\t",quote=F)

  filename1=paste("Predictions_GBLUP_GPD_9Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
}  


#Apply so maximum of replicates is 8
{
  Only8reps=removereplicates(8,Only9reps)
  Only8reps_ready=Correcterfunction(Only8reps)
  if (nrow(Only8reps_ready)==length(unique(d6$Clovershort))*8){
    print("Number of replicates pr. genotype has been reduced to 8")
    GcloverReps=Make_GRM_matrices(Only8reps_ready)[[1]]
    CloverIndep=Make_GRM_matrices(Only8reps_ready)[[2]]
    tests=testpop_generator(Only8reps_ready)
    dataframe = Only8reps_ready
    print("Starting GBLUP prediction")
    results8=mclapply(tests, GP_GBLUP)
  }
}

# Summarize and make files with results
{
  first=results8[[1]]
  second=results8[[2]]
  third=results8[[3]]
  fourth=results8[[4]]
  fifth=results8[[5]]
  sixth=results8[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)
  
  All1=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  All2=aggregate(as.numeric(All[, 4]), list(All[,1]), mean)
  StandardDeviation_gpdRes=aggregate(as.numeric(All[, 2]), list(All[,1]), sd)
  StandardDeviation_corrected=aggregate(as.numeric(All[, 3]), list(All[,1]), sd)
  Gpdmean=aggregate(as.numeric(All[, 2]), list(All[,1]), mean)

  
  correlation1=cor(All1[,2],All2[,2]) #means of replicates
  correlation1

  correlation2=cor(Gpdmean[,2],All2[,2])
  correlation2
  
  Combined=cbind(All1[,1],Gpdmean[,2],All1[,2],All2[,2],StandardDeviation_gpdRes[,2],StandardDeviation_corrected[,2])
  colnames(Combined)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno")
  
  filename=paste("Correlation_Correctedphenotype_GBLUP_GPD_8Replicates",round,".txt",sep="")
  write.table(correlation1,filename,sep="\t",quote=F)
  
  filename=paste("Correlation_gpdResCor_GBLUP_GPD_8Replicates",round,".txt",sep="")
  write.table(correlation2,filename,sep="\t",quote=F)

  filename1=paste("Predictions_GBLUP_GPD_8Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
}  

#Apply so maximum of replicates is 7
{
  Only7reps=removereplicates(7,Only8reps)
  Only7reps_ready=Correcterfunction(Only7reps)
  if (nrow(Only7reps_ready)==length(unique(d6$Clovershort))*7){
    print("Number of replicates pr. genotype has been reduced to 7")
    GcloverReps=Make_GRM_matrices(Only7reps_ready)[[1]]
    CloverIndep=Make_GRM_matrices(Only7reps_ready)[[2]]
    tests=testpop_generator(Only7reps_ready)
    dataframe = Only7reps_ready
    print("Starting GBLUP prediction")
    results7=mclapply(tests, GP_GBLUP)
  }
}

# Summarize and make files with results
{
  first=results7[[1]]
  second=results7[[2]]
  third=results7[[3]]
  fourth=results7[[4]]
  fifth=results7[[5]]
  sixth=results7[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)
  
  All1=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  All2=aggregate(as.numeric(All[, 4]), list(All[,1]), mean)
  StandardDeviation_gpdRes=aggregate(as.numeric(All[, 2]), list(All[,1]), sd)
  StandardDeviation_corrected=aggregate(as.numeric(All[, 3]), list(All[,1]), sd)
  Gpdmean=aggregate(as.numeric(All[, 2]), list(All[,1]), mean)

  
  correlation1=cor(All1[,2],All2[,2]) #means of replicates
  correlation1

  correlation2=cor(Gpdmean[,2],All2[,2])
  correlation2
  
  Combined=cbind(All1[,1],Gpdmean[,2],All1[,2],All2[,2],StandardDeviation_gpdRes[,2],StandardDeviation_corrected[,2])
  colnames(Combined)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno")
  
  filename=paste("Correlation_Correctedphenotype_GBLUP_GPD_7Replicates",round,".txt",sep="")
  write.table(correlation1,filename,sep="\t",quote=F)
  
  filename=paste("Correlation_gpdResCor_GBLUP_GPD_7Replicates",round,".txt",sep="")
  write.table(correlation2,filename,sep="\t",quote=F)

  filename1=paste("Predictions_GBLUP_GPD_7Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
}  


#Apply so maximum of replicates is 6
{
  Only6reps=removereplicates(6,Only7reps)
  Only6reps_ready=Correcterfunction(Only6reps)
  if (nrow(Only6reps_ready)==length(unique(d6$Clovershort))*6){
    print("Number of replicates pr. genotype has been reduced to 6")
    GcloverReps=Make_GRM_matrices(Only6reps_ready)[[1]]
    CloverIndep=Make_GRM_matrices(Only6reps_ready)[[2]]
    tests=testpop_generator(Only6reps_ready)
    dataframe = Only6reps_ready
    print("Starting GBLUP prediction")
    results6=mclapply(tests, GP_GBLUP)
  }
}

# Summarize and make files with results
{
  first=results6[[1]]
  second=results6[[2]]
  third=results6[[3]]
  fourth=results6[[4]]
  fifth=results6[[5]]
  sixth=results6[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)
  
  All1=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  All2=aggregate(as.numeric(All[, 4]), list(All[,1]), mean)
  StandardDeviation_gpdRes=aggregate(as.numeric(All[, 2]), list(All[,1]), sd)
  StandardDeviation_corrected=aggregate(as.numeric(All[, 3]), list(All[,1]), sd)
  Gpdmean=aggregate(as.numeric(All[, 2]), list(All[,1]), mean)

  
  correlation1=cor(All1[,2],All2[,2]) #means of replicates
  correlation1

  correlation2=cor(Gpdmean[,2],All2[,2])
  correlation2
  
  Combined=cbind(All1[,1],Gpdmean[,2],All1[,2],All2[,2],StandardDeviation_gpdRes[,2],StandardDeviation_corrected[,2])
  colnames(Combined)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno")
  
  filename=paste("Correlation_Correctedphenotype_GBLUP_GPD_6Replicates",round,".txt",sep="")
  write.table(correlation1,filename,sep="\t",quote=F)
  
  filename=paste("Correlation_gpdResCor_GBLUP_GPD_6Replicates",round,".txt",sep="")
  write.table(correlation2,filename,sep="\t",quote=F)

  filename1=paste("Predictions_GBLUP_GPD_6Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
}  


#Apply so maximum of replicates is 5
{
  Only5reps=removereplicates(5,Only6reps)
  Only5reps_ready=Correcterfunction(Only5reps)
  if (nrow(Only5reps_ready)==length(unique(d6$Clovershort))*5){
    print("Number of replicates pr. genotype has been reduced to 5")
    GcloverReps=Make_GRM_matrices(Only5reps_ready)[[1]]
    CloverIndep=Make_GRM_matrices(Only5reps_ready)[[2]]
    tests=testpop_generator(Only5reps_ready)
    dataframe = Only5reps_ready
    print("Starting GBLUP prediction")
    results5=mclapply(tests, GP_GBLUP)
  }
}

# Summarize and make files with results
{
  first=results5[[1]]
  second=results5[[2]]
  third=results5[[3]]
  fourth=results5[[4]]
  fifth=results5[[5]]
  sixth=results5[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)
    
  All1=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  All2=aggregate(as.numeric(All[, 4]), list(All[,1]), mean)
  StandardDeviation_gpdRes=aggregate(as.numeric(All[, 2]), list(All[,1]), sd)
  StandardDeviation_corrected=aggregate(as.numeric(All[, 3]), list(All[,1]), sd)
  Gpdmean=aggregate(as.numeric(All[, 2]), list(All[,1]), mean)

  
  correlation1=cor(All1[,2],All2[,2]) #means of replicates
  correlation1

  correlation2=cor(Gpdmean[,2],All2[,2])
  correlation2
  
  Combined=cbind(All1[,1],Gpdmean[,2],All1[,2],All2[,2],StandardDeviation_gpdRes[,2],StandardDeviation_corrected[,2])
  colnames(Combined)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno")
  
  filename=paste("Correlation_Correctedphenotype_GBLUP_GPD_5Replicates",round,".txt",sep="")
  write.table(correlation1,filename,sep="\t",quote=F)
  
  filename=paste("Correlation_gpdResCor_GBLUP_GPD_5Replicates",round,".txt",sep="")
  write.table(correlation2,filename,sep="\t",quote=F)

  filename1=paste("Predictions_GBLUP_GPD_5Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
}  

#Apply so maximum of replicates is 4
{
  Only4reps=removereplicates(4,Only5reps)
  Only4reps_ready=Correcterfunction(Only4reps)
  if (nrow(Only4reps_ready)==length(unique(d6$Clovershort))*4){
    print("Number of replicates pr. genotype has been reduced to 4")
    GcloverReps=Make_GRM_matrices(Only4reps_ready)[[1]]
    CloverIndep=Make_GRM_matrices(Only4reps_ready)[[2]]
    tests=testpop_generator(Only4reps_ready)
    dataframe = Only4reps_ready
    print("Starting GBLUP prediction")
    results4=mclapply(tests, GP_GBLUP)
  }
}


# Summarize and make files with results
{
  first=results4[[1]]
  second=results4[[2]]
  third=results4[[3]]
  fourth=results4[[4]]
  fifth=results4[[5]]
  sixth=results4[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)
    
  All1=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  All2=aggregate(as.numeric(All[, 4]), list(All[,1]), mean)
  StandardDeviation_gpdRes=aggregate(as.numeric(All[, 2]), list(All[,1]), sd)
  StandardDeviation_corrected=aggregate(as.numeric(All[, 3]), list(All[,1]), sd)
  Gpdmean=aggregate(as.numeric(All[, 2]), list(All[,1]), mean)
  
  correlation1=cor(All1[,2],All2[,2]) #means of replicates
  correlation1

  correlation2=cor(Gpdmean[,2],All2[,2])
  correlation2
  
  Combined=cbind(All1[,1],Gpdmean[,2],All1[,2],All2[,2],StandardDeviation_gpdRes[,2],StandardDeviation_corrected[,2])
  colnames(Combined)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno")
  
  filename=paste("Correlation_Correctedphenotype_GBLUP_GPD_4Replicates",round,".txt",sep="")
  write.table(correlation1,filename,sep="\t",quote=F)
  
  filename=paste("Correlation_gpdResCor_GBLUP_GPD_4Replicates",round,".txt",sep="")
  write.table(correlation2,filename,sep="\t",quote=F)

  filename1=paste("Predictions_GBLUP_GPD_4Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
}  




#Apply so maximum of replicates is 3
{
  Only3reps=removereplicates(3,Only4reps)
  Only3reps_ready=Correcterfunction(Only3reps)
  if (nrow(Only3reps_ready)==length(unique(d6$Clovershort))*3){
    print("Number of replicates pr. genotype has been reduced to 3")
    GcloverReps=Make_GRM_matrices(Only3reps_ready)[[1]]
    CloverIndep=Make_GRM_matrices(Only3reps_ready)[[2]]
    tests=testpop_generator(Only3reps_ready)
    dataframe = Only3reps_ready
    print("Starting GBLUP prediction")
    results3=mclapply(tests, GP_GBLUP)
  }
}


# Summarize and make files with results
{
  first=results3[[1]]
  second=results3[[2]]
  third=results3[[3]]
  fourth=results3[[4]]
  fifth=results3[[5]]
  sixth=results3[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)
  
    
  All1=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  All2=aggregate(as.numeric(All[, 4]), list(All[,1]), mean)
  StandardDeviation_gpdRes=aggregate(as.numeric(All[, 2]), list(All[,1]), sd)
  StandardDeviation_corrected=aggregate(as.numeric(All[, 3]), list(All[,1]), sd)
  Gpdmean=aggregate(as.numeric(All[, 2]), list(All[,1]), mean)

  correlation1=cor(All1[,2],All2[,2]) #means of replicates
  correlation1

  correlation2=cor(Gpdmean[,2],All2[,2])
  correlation2
  
  Combined=cbind(All1[,1],Gpdmean[,2],All1[,2],All2[,2],StandardDeviation_gpdRes[,2],StandardDeviation_corrected[,2])
  colnames(Combined)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno")
  
  filename=paste("Correlation_Correctedphenotype_GBLUP_GPD_3Replicates",round,".txt",sep="")
  write.table(correlation1,filename,sep="\t",quote=F)
  
  filename=paste("Correlation_gpdResCor_GBLUP_GPD_3Replicates",round,".txt",sep="")
  write.table(correlation2,filename,sep="\t",quote=F)

  filename1=paste("Predictions_GBLUP_GPD_3Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
}  


#Apply so maximum of replicates is 2
{
  Only2reps=removereplicates(2,Only3reps)
  Only2reps_ready=Correcterfunction(Only2reps)
  if (nrow(Only2reps_ready)==length(unique(d6$Clovershort))*2){
    print("Number of replicates pr. genotype has been reduced to 2")
    GcloverReps=Make_GRM_matrices(Only2reps_ready)[[1]]
    CloverIndep=Make_GRM_matrices(Only2reps_ready)[[2]]
    tests=testpop_generator(Only2reps_ready)
    dataframe = Only2reps_ready
    print("Starting GBLUP prediction")
    results2=mclapply(tests, GP_GBLUP)
  }
}


# Summarize and make files with results
{
  first=results2[[1]]
  second=results2[[2]]
  third=results2[[3]]
  fourth=results2[[4]]
  fifth=results2[[5]]
  sixth=results2[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)
  
    
  All1=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  All2=aggregate(as.numeric(All[, 4]), list(All[,1]), mean)
  StandardDeviation_gpdRes=aggregate(as.numeric(All[, 2]), list(All[,1]), sd)
  StandardDeviation_corrected=aggregate(as.numeric(All[, 3]), list(All[,1]), sd)
  Gpdmean=aggregate(as.numeric(All[, 2]), list(All[,1]), mean)
  
  correlation1=cor(All1[,2],All2[,2]) #means of replicates
  correlation1

  correlation2=cor(Gpdmean[,2],All2[,2])
  correlation2
  
  Combined=cbind(All1[,1],Gpdmean[,2],All1[,2],All2[,2],StandardDeviation_gpdRes[,2],StandardDeviation_corrected[,2])
  colnames(Combined)=c("Individual","gpd_rescor","Observed_correctedForFixed","GEBVs","SD_gpdResCor","SD_CorrectedPheno")
  
  filename=paste("Correlation_Correctedphenotype_GBLUP_GPD_2Replicates",round,".txt",sep="")
  write.table(correlation1,filename,sep="\t",quote=F)
  
  filename=paste("Correlation_gpdResCor_GBLUP_GPD_2Replicates",round,".txt",sep="")
  write.table(correlation2,filename,sep="\t",quote=F)

  filename1=paste("Predictions_GBLUP_GPD_2Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
}  
