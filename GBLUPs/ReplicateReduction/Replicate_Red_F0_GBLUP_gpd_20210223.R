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

{
  d4$roundRep <- paste(d4$Round, d4$Replicate, sep='_')
  #d6=d4[-which(d4$roundRep=="1_2"),]
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


# Remove all genotypes that has <10 replicates
Genotypes=(unique(d6$Clovershort))

for (genotype in Genotypes){
  idx=which(d6$Clovershort==genotype)
  if (length(idx)<10){
    d6=d6[-idx,]
    print(paste(genotype,"removed",sep=" "))
    
    GRMidx = which(colnames(GRM1)==genotype)
    GRM1 = GRM1[-GRMidx,-GRMidx]
  }
}


# Clean up
{
  d6$Rhizobium=droplevels(d6$Rhizobium) # removing levels not used in actual data
  d6$Clover=droplevels(d6$Clover) # removing levels not used in actual data
  d6=d6[order(d6$Clovershort),] # make sure it is in alphabetic order like the GRM
}

# Divide into 6 populations for GP
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

grouping=list(testpop1,testpop2,testpop3,testpop4,testpop5,testpop6)
name=paste("grouping",round,".txt",sep="")
sink(name)
print(grouping)
sink()


############################################################

# Now remove replicates so each genotype has a maximum of of desired number (maxreplicates) and calculate mean phenotypes based on replicates left

removereplicates <- function(maxreplicates,dataframe){
  set.seed(NULL)
  for (genotype in Genotypes){
    replicateidx=which(dataframe$Clovershort==genotype)
    if (length(replicateidx)>maxreplicates){
      numbertoremove=length(replicateidx)-maxreplicates
      remove=sample(replicateidx,numbertoremove)
      dataframe=dataframe[-remove,]
    }
    print(paste("Number of replicates pr. genotype has been reduced to:",maxreplicates,sep=""))
    gpdmeans=aggregate(dataframe$growth_per_day, list(dataframe$Clovershort), mean) # calculate averages from reduced dataframe
    colnames(gpdmeans)=c("Clovershort","growth_per_day")
  }
  return(list(gpdmeans,dataframe))
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
  
  
  gpdmeans_training=dataframe[-testpop,] # limit the dataframe to only the individuals allowed for training the model
  gpdmeans_training_ready=na.omit(gpdmeans_training, cols = c("growth_per_day")) # remember that gpd na inidividuals should be removed whether or not they are in the training pop or not
  
  ind_not_in_train=dataframe$Clovershort[testpop]
  IndividualsToRemoveGRM=which(colnames(GRM1) %in% ind_not_in_train)
  
  GRM_trn = GRM1[-IndividualsToRemoveGRM,-IndividualsToRemoveGRM]
  
  # Run the GBLUP model on full training population to extract GEBVs
  yNA=gpdmeans_training_ready$growth_per_day
  ETA=list(list(K=GRM_trn,model="RKHS"))
  GBLUP=BGLR(y=yNA,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("GBLUP",round))
  matrix=cbind(as.character(gpdmeans_training_ready$Clovershort),as.numeric(gpdmeans_training_ready$growth_per_day),as.numeric(GBLUP$ETA[[1]]$u))
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
  matrix1=cbind(as.character(dataframe$Clovershort[testpop]),as.numeric(dataframe$growth_per_day[testpop]),as.numeric(as.character(GEBVpred)))
  colnames(matrix1)=c("ID", "Observed", "GEBV")
  return(matrix1)
}


#Apply so maximum of replicates is 10
{
  run10=removereplicates(10,d6)
  Only10reps_avg=run10[[1]]
  Only10reps=run10[[2]]
  head(Only10reps_avg)
  tests=testpop_generator(Only10reps_avg)
  dataframe = Only10reps_avg
  print("Starting GBLUP prediction")
  results10=mclapply(tests, GP_GBLUP)
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
  correlation=cor(as.numeric(as.character(All[,2])),as.numeric(as.character(All[,3]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_gpd_10Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_gpd_10Replicates",round,".txt",sep="")
  write.table(All,filename1,sep="\t",quote=F,row.names=F)
}  


#Apply so maximum of replicates is 9
{
  run9=removereplicates(9,Only10reps)
  Only9reps_avg=run9[[1]]
  Only9reps=run9[[2]]
  tests=testpop_generator(Only9reps_avg)
  dataframe = Only9reps_avg
  print("Starting GBLUP prediction")
  results9=mclapply(tests, GP_GBLUP)
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
  correlation=cor(as.numeric(as.character(All[,2])),as.numeric(as.character(All[,3]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_gpd_9Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_gpd_9Replicates",round,".txt",sep="")
  write.table(All,filename1,sep="\t",quote=F,row.names=F)
}  


#Apply so maximum of replicates is 8
{
  run8=removereplicates(8,Only9reps)
  Only8reps_avg=run8[[1]]
  Only8reps=run8[[2]]
  tests=testpop_generator(Only8reps_avg)
  dataframe = Only8reps_avg
  print("Starting GBLUP prediction")
  results8=mclapply(tests, GP_GBLUP)
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
  correlation=cor(as.numeric(as.character(All[,2])),as.numeric(as.character(All[,3]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_gpd_8Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_gpd_8Replicates",round,".txt",sep="")
  write.table(All,filename1,sep="\t",quote=F,row.names=F)
}  


#Apply so maximum of replicates is 7
{
  run7=removereplicates(7,Only8reps)
  Only7reps_avg=run7[[1]]
  Only7reps=run7[[2]]
  tests=testpop_generator(Only7reps_avg)
  dataframe = Only7reps_avg
  print("Starting GBLUP prediction")
  results7=mclapply(tests, GP_GBLUP)
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
  correlation=cor(as.numeric(as.character(All[,2])),as.numeric(as.character(All[,3]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_gpd_7Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_gpd_7Replicates",round,".txt",sep="")
  write.table(All,filename1,sep="\t",quote=F,row.names=F)
}  


#Apply so maximum of replicates is 6
{
  run6= removereplicates(6,Only7reps)
  Only6reps_avg=run6[[1]]
  Only6reps=run6[[2]]
  tests=testpop_generator(Only6reps_avg)
  dataframe = Only6reps_avg
  print("Starting GBLUP prediction")
  results6=mclapply(tests, GP_GBLUP)
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
  correlation=cor(as.numeric(as.character(All[,2])),as.numeric(as.character(All[,3]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_gpd_6Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_gpd_6Replicates",round,".txt",sep="")
  write.table(All,filename1,sep="\t",quote=F,row.names=F)
}  


#Apply so maximum of replicates is 5
{
  run5 = removereplicates(5,Only6reps)
  Only5reps_avg=run5[[1]]
  Only5reps=run5[[2]]
  tests=testpop_generator(Only5reps_avg)
  dataframe = Only5reps_avg
  print("Starting GBLUP prediction")
  results5=mclapply(tests, GP_GBLUP)
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
  correlation=cor(as.numeric(as.character(All[,2])),as.numeric(as.character(All[,3]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_gpd_5Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_gpd_5Replicates",round,".txt",sep="")
  write.table(All,filename1,sep="\t",quote=F,row.names=F)
}  

#Apply so maximum of replicates is 4
{
  run4= removereplicates(4,Only5reps)
  Only4reps_avg=run4[[1]]
  Only4reps=run4[[2]]
  tests=testpop_generator(Only4reps_avg)
  dataframe = Only4reps_avg
  print("Starting GBLUP prediction")
  results4=mclapply(tests, GP_GBLUP)
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
  correlation=cor(as.numeric(as.character(All[,2])),as.numeric(as.character(All[,3]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_gpd_4Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_gpd_4Replicates",round,".txt",sep="")
  write.table(All,filename1,sep="\t",quote=F,row.names=F)
}  




#Apply so maximum of replicates is 3
{
  run3= removereplicates(3,Only4reps)
  Only3reps_avg=run3[[1]]
  Only3reps=run3[[2]]
  tests=testpop_generator(Only3reps_avg)
  dataframe = Only3reps_avg
  print("Starting GBLUP prediction")
  results3=mclapply(tests, GP_GBLUP)
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
  correlation=cor(as.numeric(as.character(All[,2])),as.numeric(as.character(All[,3]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_gpd_3Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_gpd_3Replicates",round,".txt",sep="")
  write.table(All,filename1,sep="\t",quote=F,row.names=F)
}  


#Apply so maximum of replicates is 2
{
  run2=removereplicates(2,Only3reps)
  Only2reps_avg=run2[[1]]
  Only2reps=run2[[2]]
  tests=testpop_generator(Only2reps_avg)
  dataframe = Only2reps_avg
  print("Starting GBLUP prediction")
  results2=mclapply(tests, GP_GBLUP)
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
  correlation=cor(as.numeric(as.character(All[,2])),as.numeric(as.character(All[,3]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_gpd_2Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_gpd_2Replicates",round,".txt",sep="")
  write.table(All,filename1,sep="\t",quote=F,row.names=F)
}  


#Apply so maximum of replicates is 1
{
  run1=removereplicates(1,Only2reps)
  Only1reps_avg=run1[[1]]
  Only1reps=run1[[2]]
  tests=testpop_generator(Only1reps_avg)
  dataframe = Only1reps_avg
  print("Starting GBLUP prediction")
  results1=mclapply(tests, GP_GBLUP)
}

# Summarize and make files with results
{
  first=results1[[1]]
  second=results1[[2]]
  third=results1[[3]]
  fourth=results1[[4]]
  fifth=results1[[5]]
  sixth=results1[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)
  correlation=cor(as.numeric(as.character(All[,2])),as.numeric(as.character(All[,3]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_gpd_1Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_gpd_1Replicates",round,".txt",sep="")
  write.table(All,filename1,sep="\t",quote=F,row.names=F)
}  

