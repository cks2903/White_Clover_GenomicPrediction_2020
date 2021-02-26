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
  
  F1observations=read.table("F1results.csv",sep=",",header=T, stringsAsFactors = F)
  F1observations=as.data.frame(F1observations)
  
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

# Divide into F1 parental groups
{
  parentpop1names=c("Aoost_01","Aoost_08","Aoost_09","Banna_02","Banna_03","Banna_07")
  parentpop2names=c("Aearl_08","Ccyma_03","Llanc_06","Aaran_08")
  parentpop3names=c("Aearl_05","Clfin_02","Ctain_05","Mrida_04")
  parentpop4names=c("Clfin_03","Ctain_05","Volin_01","Aaran_04")
  parentpop5names=c("Aoost_02","Ilona_09","Llanc_09","Sster_01")
  parentpop6names=c("Ilona_05","Kdike_09","Llanc_09","Aalon_03")
  parentpop7names=c("Ancor_10","Borek_06","Ctain_09","Rbani_02")
  parentpop8names=c("Ancor_04","Aoost_10","Clfin_08","Kdike_08")
  parentpop9names=c("Aoost_01","Aoost_08","Banna_02","Rbani_02","Sster_01","Sster_06") 
}


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
    iSizemeans=aggregate(dataframe$InitialSize, list(dataframe$Clovershort), mean) # calculate averages from reduced dataframe
    colnames(iSizemeans)=c("Clovershort","InitialSize")
  }
  return(list(iSizemeans,dataframe))
}



  
testpop_generator<-function(dataframe){
  #Find indexes for test population  
  parentpop1names_idx=which(dataframe$Clovershort %in% parentpop1names)
  parentpop2names_idx=which(dataframe$Clovershort %in% parentpop2names)
  parentpop3names_idx=which(dataframe$Clovershort %in% parentpop3names)
  parentpop4names_idx=which(dataframe$Clovershort %in% parentpop4names)
  parentpop5names_idx=which(dataframe$Clovershort %in% parentpop5names)
  parentpop6names_idx=which(dataframe$Clovershort %in% parentpop6names)
  parentpop7names_idx=which(dataframe$Clovershort %in% parentpop7names)
  parentpop8names_idx=which(dataframe$Clovershort %in% parentpop8names)
  parentpop9names_idx=which(dataframe$Clovershort %in% parentpop9names)

  
  tests=list(parentpop1names_idx,parentpop2names_idx,parentpop3names_idx,parentpop4names_idx,parentpop5names_idx,parentpop6names_idx,parentpop7names_idx,parentpop8names_idx,parentpop9names_idx)
  
  return(tests)
}

GP_GBLUP<-function(testpop){
  
  ################  ################  ################  ################
  ##start by estimating GEBVs for training population individuals 
  ################  ################  ################  ################
  
  
  iSizemeans_training=dataframe[-testpop,] # limit the dataframe to only the individuals allowed for training the model
  iSizemeans_training_ready=na.omit(iSizemeans_training, cols = c("iSize")) # remember that gpd na inidividuals should be removed whether or not they are in the training pop or not
  
  ind_not_in_train=dataframe$Clovershort[testpop]
  IndividualsToRemoveGRM=which(colnames(GRM1) %in% ind_not_in_train)
  
  GRM_trn = GRM1[-IndividualsToRemoveGRM,-IndividualsToRemoveGRM]
  
  # Run the GBLUP model on full training population to extract GEBVs
  yNA=iSizemeans_training_ready$InitialSize
  ETA=list(list(K=GRM_trn,model="RKHS"))
  GBLUP=BGLR(y=yNA,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("GBLUP",round))
  matrix=cbind(as.character(iSizemeans_training_ready$Clovershort),as.numeric(iSizemeans_training_ready$InitialSize),as.numeric(GBLUP$ETA[[1]]$u))
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
  matrix1=cbind(as.character(dataframe$Clovershort[testpop]),as.numeric(dataframe$InitialSize[testpop]),as.numeric(as.character(GEBVpred)))
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
  GEBVsParentPop1=results10[[1]]
  GEBVsParentPop2=results10[[2]]
  GEBVsParentPop3=results10[[3]]
  GEBVsParentPop4=results10[[4]]
  GEBVsParentPop5=results10[[5]]
  GEBVsParentPop6=results10[[6]]
  GEBVsParentPop7=results10[[7]]
  GEBVsParentPop8=results10[[8]]
  GEBVsParentPop9=results10[[9]]

  GEBVsParentPop1_table1=  mean(as.numeric(as.character(GEBVsParentPop1[,2]))) # parental average phenotype
  GEBVsParentPop1_table2=  mean(as.numeric(as.character(GEBVsParentPop1[,3]))) # parental average GEBV
  GEBVsParentPop2_table1=  mean(as.numeric(as.character(GEBVsParentPop2[,2]))) # parental average phenotype
  GEBVsParentPop2_table2=  mean(as.numeric(as.character(GEBVsParentPop2[,3]))) # parental average GEBV
  GEBVsParentPop3_table1=  mean(as.numeric(as.character(GEBVsParentPop3[,2]))) # parental average phenotype
  GEBVsParentPop3_table2=  mean(as.numeric(as.character(GEBVsParentPop3[,3]))) # parental average GEBV
  GEBVsParentPop4_table1=  mean(as.numeric(as.character(GEBVsParentPop4[,2]))) # parental average phenotype
  GEBVsParentPop4_table2=  mean(as.numeric(as.character(GEBVsParentPop4[,3]))) # parental average GEBV
  GEBVsParentPop5_table1=  mean(as.numeric(as.character(GEBVsParentPop5[,2]))) # parental average phenotype
  GEBVsParentPop5_table2=  mean(as.numeric(as.character(GEBVsParentPop5[,3]))) # parental average GEBV
  GEBVsParentPop6_table1=  mean(as.numeric(as.character(GEBVsParentPop6[,2]))) # parental average phenotype
  GEBVsParentPop6_table2=  mean(as.numeric(as.character(GEBVsParentPop6[,3]))) # parental average GEBV
  GEBVsParentPop7_table1=  mean(as.numeric(as.character(GEBVsParentPop7[,2]))) # parental average phenotype
  GEBVsParentPop7_table2=  mean(as.numeric(as.character(GEBVsParentPop7[,3]))) # parental average GEBV
  GEBVsParentPop8_table1=  mean(as.numeric(as.character(GEBVsParentPop8[,2]))) # parental average phenotype
  GEBVsParentPop8_table2=  mean(as.numeric(as.character(GEBVsParentPop8[,3]))) # parental average GEBV
  GEBVsParentPop9_table1=  mean(as.numeric(as.character(GEBVsParentPop9[,2]))) # parental average phenotype
  GEBVsParentPop9_table2=  mean(as.numeric(as.character(GEBVsParentPop9[,3]))) # parental average GEBV
  
  All=rbind(GEBVsParentPop1_table2,GEBVsParentPop2_table2,GEBVsParentPop3_table2,GEBVsParentPop4_table2,GEBVsParentPop5_table2,GEBVsParentPop6_table2,GEBVsParentPop7_table2,GEBVsParentPop8_table2,GEBVsParentPop9_table2)
  colnames(All)="GEBV"
  All2=rbind(GEBVsParentPop1_table1,GEBVsParentPop2_table1,GEBVsParentPop3_table1,GEBVsParentPop4_table1,GEBVsParentPop5_table1,GEBVsParentPop6_table1,GEBVsParentPop7_table1,GEBVsParentPop8_table1,GEBVsParentPop9_table1)
  colnames(All2)="AvgParentalPheno"
  Merged= cbind(F1observations,All,All2)
  
  correlation=cor(as.numeric(as.character(Merged[,3])),as.numeric(as.character(Merged[,4]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_iSize_10Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_iSize_10Replicates",round,".txt",sep="")
  write.table(Merged,filename1,sep="\t",quote=F,row.names=F)
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
  GEBVsParentPop1=results9[[1]]
  GEBVsParentPop2=results9[[2]]
  GEBVsParentPop3=results9[[3]]
  GEBVsParentPop4=results9[[4]]
  GEBVsParentPop5=results9[[5]]
  GEBVsParentPop6=results9[[6]]
  GEBVsParentPop7=results9[[7]]
  GEBVsParentPop8=results9[[8]]
  GEBVsParentPop9=results9[[9]]
  
  GEBVsParentPop1_table1=  mean(as.numeric(as.character(GEBVsParentPop1[,2]))) # parental average phenotype
  GEBVsParentPop1_table2=  mean(as.numeric(as.character(GEBVsParentPop1[,3]))) # parental average GEBV
  GEBVsParentPop2_table1=  mean(as.numeric(as.character(GEBVsParentPop2[,2]))) # parental average phenotype
  GEBVsParentPop2_table2=  mean(as.numeric(as.character(GEBVsParentPop2[,3]))) # parental average GEBV
  GEBVsParentPop3_table1=  mean(as.numeric(as.character(GEBVsParentPop3[,2]))) # parental average phenotype
  GEBVsParentPop3_table2=  mean(as.numeric(as.character(GEBVsParentPop3[,3]))) # parental average GEBV
  GEBVsParentPop4_table1=  mean(as.numeric(as.character(GEBVsParentPop4[,2]))) # parental average phenotype
  GEBVsParentPop4_table2=  mean(as.numeric(as.character(GEBVsParentPop4[,3]))) # parental average GEBV
  GEBVsParentPop5_table1=  mean(as.numeric(as.character(GEBVsParentPop5[,2]))) # parental average phenotype
  GEBVsParentPop5_table2=  mean(as.numeric(as.character(GEBVsParentPop5[,3]))) # parental average GEBV
  GEBVsParentPop6_table1=  mean(as.numeric(as.character(GEBVsParentPop6[,2]))) # parental average phenotype
  GEBVsParentPop6_table2=  mean(as.numeric(as.character(GEBVsParentPop6[,3]))) # parental average GEBV
  GEBVsParentPop7_table1=  mean(as.numeric(as.character(GEBVsParentPop7[,2]))) # parental average phenotype
  GEBVsParentPop7_table2=  mean(as.numeric(as.character(GEBVsParentPop7[,3]))) # parental average GEBV
  GEBVsParentPop8_table1=  mean(as.numeric(as.character(GEBVsParentPop8[,2]))) # parental average phenotype
  GEBVsParentPop8_table2=  mean(as.numeric(as.character(GEBVsParentPop8[,3]))) # parental average GEBV
  GEBVsParentPop9_table1=  mean(as.numeric(as.character(GEBVsParentPop9[,2]))) # parental average phenotype
  GEBVsParentPop9_table2=  mean(as.numeric(as.character(GEBVsParentPop9[,3]))) # parental average GEBV
  
  All=rbind(GEBVsParentPop1_table2,GEBVsParentPop2_table2,GEBVsParentPop3_table2,GEBVsParentPop4_table2,GEBVsParentPop5_table2,GEBVsParentPop6_table2,GEBVsParentPop7_table2,GEBVsParentPop8_table2,GEBVsParentPop9_table2)
  colnames(All)="GEBV"
  All2=rbind(GEBVsParentPop1_table1,GEBVsParentPop2_table1,GEBVsParentPop3_table1,GEBVsParentPop4_table1,GEBVsParentPop5_table1,GEBVsParentPop6_table1,GEBVsParentPop7_table1,GEBVsParentPop8_table1,GEBVsParentPop9_table1)
  colnames(All2)="AvgParentalPheno"
  Merged= cbind(F1observations,All,All2)
  
  correlation=cor(as.numeric(as.character(Merged[,3])),as.numeric(as.character(Merged[,4]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_iSize_9Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_iSize_9Replicates",round,".txt",sep="")
  write.table(Merged,filename1,sep="\t",quote=F,row.names=F)
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
  GEBVsParentPop1=results8[[1]]
  GEBVsParentPop2=results8[[2]]
  GEBVsParentPop3=results8[[3]]
  GEBVsParentPop4=results8[[4]]
  GEBVsParentPop5=results8[[5]]
  GEBVsParentPop6=results8[[6]]
  GEBVsParentPop7=results8[[7]]
  GEBVsParentPop8=results8[[8]]
  GEBVsParentPop9=results8[[9]]
  
  GEBVsParentPop1_table1=  mean(as.numeric(as.character(GEBVsParentPop1[,2]))) # parental average phenotype
  GEBVsParentPop1_table2=  mean(as.numeric(as.character(GEBVsParentPop1[,3]))) # parental average GEBV
  GEBVsParentPop2_table1=  mean(as.numeric(as.character(GEBVsParentPop2[,2]))) # parental average phenotype
  GEBVsParentPop2_table2=  mean(as.numeric(as.character(GEBVsParentPop2[,3]))) # parental average GEBV
  GEBVsParentPop3_table1=  mean(as.numeric(as.character(GEBVsParentPop3[,2]))) # parental average phenotype
  GEBVsParentPop3_table2=  mean(as.numeric(as.character(GEBVsParentPop3[,3]))) # parental average GEBV
  GEBVsParentPop4_table1=  mean(as.numeric(as.character(GEBVsParentPop4[,2]))) # parental average phenotype
  GEBVsParentPop4_table2=  mean(as.numeric(as.character(GEBVsParentPop4[,3]))) # parental average GEBV
  GEBVsParentPop5_table1=  mean(as.numeric(as.character(GEBVsParentPop5[,2]))) # parental average phenotype
  GEBVsParentPop5_table2=  mean(as.numeric(as.character(GEBVsParentPop5[,3]))) # parental average GEBV
  GEBVsParentPop6_table1=  mean(as.numeric(as.character(GEBVsParentPop6[,2]))) # parental average phenotype
  GEBVsParentPop6_table2=  mean(as.numeric(as.character(GEBVsParentPop6[,3]))) # parental average GEBV
  GEBVsParentPop7_table1=  mean(as.numeric(as.character(GEBVsParentPop7[,2]))) # parental average phenotype
  GEBVsParentPop7_table2=  mean(as.numeric(as.character(GEBVsParentPop7[,3]))) # parental average GEBV
  GEBVsParentPop8_table1=  mean(as.numeric(as.character(GEBVsParentPop8[,2]))) # parental average phenotype
  GEBVsParentPop8_table2=  mean(as.numeric(as.character(GEBVsParentPop8[,3]))) # parental average GEBV
  GEBVsParentPop9_table1=  mean(as.numeric(as.character(GEBVsParentPop9[,2]))) # parental average phenotype
  GEBVsParentPop9_table2=  mean(as.numeric(as.character(GEBVsParentPop9[,3]))) # parental average GEBV
  
  All=rbind(GEBVsParentPop1_table2,GEBVsParentPop2_table2,GEBVsParentPop3_table2,GEBVsParentPop4_table2,GEBVsParentPop5_table2,GEBVsParentPop6_table2,GEBVsParentPop7_table2,GEBVsParentPop8_table2,GEBVsParentPop9_table2)
  colnames(All)="GEBV"
  All2=rbind(GEBVsParentPop1_table1,GEBVsParentPop2_table1,GEBVsParentPop3_table1,GEBVsParentPop4_table1,GEBVsParentPop5_table1,GEBVsParentPop6_table1,GEBVsParentPop7_table1,GEBVsParentPop8_table1,GEBVsParentPop9_table1)
  colnames(All2)="AvgParentalPheno"
  Merged= cbind(F1observations,All,All2)
  
  correlation=cor(as.numeric(as.character(Merged[,3])),as.numeric(as.character(Merged[,4]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_iSize_8Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_iSize_8Replicates",round,".txt",sep="")
  write.table(Merged,filename1,sep="\t",quote=F,row.names=F)
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
  GEBVsParentPop1=results7[[1]]
  GEBVsParentPop2=results7[[2]]
  GEBVsParentPop3=results7[[3]]
  GEBVsParentPop4=results7[[4]]
  GEBVsParentPop5=results7[[5]]
  GEBVsParentPop6=results7[[6]]
  GEBVsParentPop7=results7[[7]]
  GEBVsParentPop8=results7[[8]]
  GEBVsParentPop9=results7[[9]]
  
  GEBVsParentPop1_table1=  mean(as.numeric(as.character(GEBVsParentPop1[,2]))) # parental average phenotype
  GEBVsParentPop1_table2=  mean(as.numeric(as.character(GEBVsParentPop1[,3]))) # parental average GEBV
  GEBVsParentPop2_table1=  mean(as.numeric(as.character(GEBVsParentPop2[,2]))) # parental average phenotype
  GEBVsParentPop2_table2=  mean(as.numeric(as.character(GEBVsParentPop2[,3]))) # parental average GEBV
  GEBVsParentPop3_table1=  mean(as.numeric(as.character(GEBVsParentPop3[,2]))) # parental average phenotype
  GEBVsParentPop3_table2=  mean(as.numeric(as.character(GEBVsParentPop3[,3]))) # parental average GEBV
  GEBVsParentPop4_table1=  mean(as.numeric(as.character(GEBVsParentPop4[,2]))) # parental average phenotype
  GEBVsParentPop4_table2=  mean(as.numeric(as.character(GEBVsParentPop4[,3]))) # parental average GEBV
  GEBVsParentPop5_table1=  mean(as.numeric(as.character(GEBVsParentPop5[,2]))) # parental average phenotype
  GEBVsParentPop5_table2=  mean(as.numeric(as.character(GEBVsParentPop5[,3]))) # parental average GEBV
  GEBVsParentPop6_table1=  mean(as.numeric(as.character(GEBVsParentPop6[,2]))) # parental average phenotype
  GEBVsParentPop6_table2=  mean(as.numeric(as.character(GEBVsParentPop6[,3]))) # parental average GEBV
  GEBVsParentPop7_table1=  mean(as.numeric(as.character(GEBVsParentPop7[,2]))) # parental average phenotype
  GEBVsParentPop7_table2=  mean(as.numeric(as.character(GEBVsParentPop7[,3]))) # parental average GEBV
  GEBVsParentPop8_table1=  mean(as.numeric(as.character(GEBVsParentPop8[,2]))) # parental average phenotype
  GEBVsParentPop8_table2=  mean(as.numeric(as.character(GEBVsParentPop8[,3]))) # parental average GEBV
  GEBVsParentPop9_table1=  mean(as.numeric(as.character(GEBVsParentPop9[,2]))) # parental average phenotype
  GEBVsParentPop9_table2=  mean(as.numeric(as.character(GEBVsParentPop9[,3]))) # parental average GEBV
  
  All=rbind(GEBVsParentPop1_table2,GEBVsParentPop2_table2,GEBVsParentPop3_table2,GEBVsParentPop4_table2,GEBVsParentPop5_table2,GEBVsParentPop6_table2,GEBVsParentPop7_table2,GEBVsParentPop8_table2,GEBVsParentPop9_table2)
  colnames(All)="GEBV"
  All2=rbind(GEBVsParentPop1_table1,GEBVsParentPop2_table1,GEBVsParentPop3_table1,GEBVsParentPop4_table1,GEBVsParentPop5_table1,GEBVsParentPop6_table1,GEBVsParentPop7_table1,GEBVsParentPop8_table1,GEBVsParentPop9_table1)
  colnames(All2)="AvgParentalPheno"
  Merged= cbind(F1observations,All,All2)
  
  correlation=cor(as.numeric(as.character(Merged[,3])),as.numeric(as.character(Merged[,4]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_iSize_7Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_iSize_7Replicates",round,".txt",sep="")
  write.table(Merged,filename1,sep="\t",quote=F,row.names=F)
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
  GEBVsParentPop1=results6[[1]]
  GEBVsParentPop2=results6[[2]]
  GEBVsParentPop3=results6[[3]]
  GEBVsParentPop4=results6[[4]]
  GEBVsParentPop5=results6[[5]]
  GEBVsParentPop6=results6[[6]]
  GEBVsParentPop7=results6[[7]]
  GEBVsParentPop8=results6[[8]]
  GEBVsParentPop9=results6[[9]]
  
  GEBVsParentPop1_table1=  mean(as.numeric(as.character(GEBVsParentPop1[,2]))) # parental average phenotype
  GEBVsParentPop1_table2=  mean(as.numeric(as.character(GEBVsParentPop1[,3]))) # parental average GEBV
  GEBVsParentPop2_table1=  mean(as.numeric(as.character(GEBVsParentPop2[,2]))) # parental average phenotype
  GEBVsParentPop2_table2=  mean(as.numeric(as.character(GEBVsParentPop2[,3]))) # parental average GEBV
  GEBVsParentPop3_table1=  mean(as.numeric(as.character(GEBVsParentPop3[,2]))) # parental average phenotype
  GEBVsParentPop3_table2=  mean(as.numeric(as.character(GEBVsParentPop3[,3]))) # parental average GEBV
  GEBVsParentPop4_table1=  mean(as.numeric(as.character(GEBVsParentPop4[,2]))) # parental average phenotype
  GEBVsParentPop4_table2=  mean(as.numeric(as.character(GEBVsParentPop4[,3]))) # parental average GEBV
  GEBVsParentPop5_table1=  mean(as.numeric(as.character(GEBVsParentPop5[,2]))) # parental average phenotype
  GEBVsParentPop5_table2=  mean(as.numeric(as.character(GEBVsParentPop5[,3]))) # parental average GEBV
  GEBVsParentPop6_table1=  mean(as.numeric(as.character(GEBVsParentPop6[,2]))) # parental average phenotype
  GEBVsParentPop6_table2=  mean(as.numeric(as.character(GEBVsParentPop6[,3]))) # parental average GEBV
  GEBVsParentPop7_table1=  mean(as.numeric(as.character(GEBVsParentPop7[,2]))) # parental average phenotype
  GEBVsParentPop7_table2=  mean(as.numeric(as.character(GEBVsParentPop7[,3]))) # parental average GEBV
  GEBVsParentPop8_table1=  mean(as.numeric(as.character(GEBVsParentPop8[,2]))) # parental average phenotype
  GEBVsParentPop8_table2=  mean(as.numeric(as.character(GEBVsParentPop8[,3]))) # parental average GEBV
  GEBVsParentPop9_table1=  mean(as.numeric(as.character(GEBVsParentPop9[,2]))) # parental average phenotype
  GEBVsParentPop9_table2=  mean(as.numeric(as.character(GEBVsParentPop9[,3]))) # parental average GEBV
  
  All=rbind(GEBVsParentPop1_table2,GEBVsParentPop2_table2,GEBVsParentPop3_table2,GEBVsParentPop4_table2,GEBVsParentPop5_table2,GEBVsParentPop6_table2,GEBVsParentPop7_table2,GEBVsParentPop8_table2,GEBVsParentPop9_table2)
  colnames(All)="GEBV"
  All2=rbind(GEBVsParentPop1_table1,GEBVsParentPop2_table1,GEBVsParentPop3_table1,GEBVsParentPop4_table1,GEBVsParentPop5_table1,GEBVsParentPop6_table1,GEBVsParentPop7_table1,GEBVsParentPop8_table1,GEBVsParentPop9_table1)
  colnames(All2)="AvgParentalPheno"
  Merged= cbind(F1observations,All,All2)
  
  correlation=cor(as.numeric(as.character(Merged[,3])),as.numeric(as.character(Merged[,4]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_iSize_6Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_iSize_6Replicates",round,".txt",sep="")
  write.table(Merged,filename1,sep="\t",quote=F,row.names=F)
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
  GEBVsParentPop1=results5[[1]]
  GEBVsParentPop2=results5[[2]]
  GEBVsParentPop3=results5[[3]]
  GEBVsParentPop4=results5[[4]]
  GEBVsParentPop5=results5[[5]]
  GEBVsParentPop6=results5[[6]]
  GEBVsParentPop7=results5[[7]]
  GEBVsParentPop8=results5[[8]]
  GEBVsParentPop9=results5[[9]]
  
  GEBVsParentPop1_table1=  mean(as.numeric(as.character(GEBVsParentPop1[,2]))) # parental average phenotype
  GEBVsParentPop1_table2=  mean(as.numeric(as.character(GEBVsParentPop1[,3]))) # parental average GEBV
  GEBVsParentPop2_table1=  mean(as.numeric(as.character(GEBVsParentPop2[,2]))) # parental average phenotype
  GEBVsParentPop2_table2=  mean(as.numeric(as.character(GEBVsParentPop2[,3]))) # parental average GEBV
  GEBVsParentPop3_table1=  mean(as.numeric(as.character(GEBVsParentPop3[,2]))) # parental average phenotype
  GEBVsParentPop3_table2=  mean(as.numeric(as.character(GEBVsParentPop3[,3]))) # parental average GEBV
  GEBVsParentPop4_table1=  mean(as.numeric(as.character(GEBVsParentPop4[,2]))) # parental average phenotype
  GEBVsParentPop4_table2=  mean(as.numeric(as.character(GEBVsParentPop4[,3]))) # parental average GEBV
  GEBVsParentPop5_table1=  mean(as.numeric(as.character(GEBVsParentPop5[,2]))) # parental average phenotype
  GEBVsParentPop5_table2=  mean(as.numeric(as.character(GEBVsParentPop5[,3]))) # parental average GEBV
  GEBVsParentPop6_table1=  mean(as.numeric(as.character(GEBVsParentPop6[,2]))) # parental average phenotype
  GEBVsParentPop6_table2=  mean(as.numeric(as.character(GEBVsParentPop6[,3]))) # parental average GEBV
  GEBVsParentPop7_table1=  mean(as.numeric(as.character(GEBVsParentPop7[,2]))) # parental average phenotype
  GEBVsParentPop7_table2=  mean(as.numeric(as.character(GEBVsParentPop7[,3]))) # parental average GEBV
  GEBVsParentPop8_table1=  mean(as.numeric(as.character(GEBVsParentPop8[,2]))) # parental average phenotype
  GEBVsParentPop8_table2=  mean(as.numeric(as.character(GEBVsParentPop8[,3]))) # parental average GEBV
  GEBVsParentPop9_table1=  mean(as.numeric(as.character(GEBVsParentPop9[,2]))) # parental average phenotype
  GEBVsParentPop9_table2=  mean(as.numeric(as.character(GEBVsParentPop9[,3]))) # parental average GEBV
  
  All=rbind(GEBVsParentPop1_table2,GEBVsParentPop2_table2,GEBVsParentPop3_table2,GEBVsParentPop4_table2,GEBVsParentPop5_table2,GEBVsParentPop6_table2,GEBVsParentPop7_table2,GEBVsParentPop8_table2,GEBVsParentPop9_table2)
  colnames(All)="GEBV"
  All2=rbind(GEBVsParentPop1_table1,GEBVsParentPop2_table1,GEBVsParentPop3_table1,GEBVsParentPop4_table1,GEBVsParentPop5_table1,GEBVsParentPop6_table1,GEBVsParentPop7_table1,GEBVsParentPop8_table1,GEBVsParentPop9_table1)
  colnames(All2)="AvgParentalPheno"
  Merged= cbind(F1observations,All,All2)
  
  correlation=cor(as.numeric(as.character(Merged[,3])),as.numeric(as.character(Merged[,4]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_iSize_5Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_iSize_5Replicates",round,".txt",sep="")
  write.table(Merged,filename1,sep="\t",quote=F,row.names=F)
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
  GEBVsParentPop1=results4[[1]]
  GEBVsParentPop2=results4[[2]]
  GEBVsParentPop3=results4[[3]]
  GEBVsParentPop4=results4[[4]]
  GEBVsParentPop5=results4[[5]]
  GEBVsParentPop6=results4[[6]]
  GEBVsParentPop7=results4[[7]]
  GEBVsParentPop8=results4[[8]]
  GEBVsParentPop9=results4[[9]]
  
  GEBVsParentPop1_table1=  mean(as.numeric(as.character(GEBVsParentPop1[,2]))) # parental average phenotype
  GEBVsParentPop1_table2=  mean(as.numeric(as.character(GEBVsParentPop1[,3]))) # parental average GEBV
  GEBVsParentPop2_table1=  mean(as.numeric(as.character(GEBVsParentPop2[,2]))) # parental average phenotype
  GEBVsParentPop2_table2=  mean(as.numeric(as.character(GEBVsParentPop2[,3]))) # parental average GEBV
  GEBVsParentPop3_table1=  mean(as.numeric(as.character(GEBVsParentPop3[,2]))) # parental average phenotype
  GEBVsParentPop3_table2=  mean(as.numeric(as.character(GEBVsParentPop3[,3]))) # parental average GEBV
  GEBVsParentPop4_table1=  mean(as.numeric(as.character(GEBVsParentPop4[,2]))) # parental average phenotype
  GEBVsParentPop4_table2=  mean(as.numeric(as.character(GEBVsParentPop4[,3]))) # parental average GEBV
  GEBVsParentPop5_table1=  mean(as.numeric(as.character(GEBVsParentPop5[,2]))) # parental average phenotype
  GEBVsParentPop5_table2=  mean(as.numeric(as.character(GEBVsParentPop5[,3]))) # parental average GEBV
  GEBVsParentPop6_table1=  mean(as.numeric(as.character(GEBVsParentPop6[,2]))) # parental average phenotype
  GEBVsParentPop6_table2=  mean(as.numeric(as.character(GEBVsParentPop6[,3]))) # parental average GEBV
  GEBVsParentPop7_table1=  mean(as.numeric(as.character(GEBVsParentPop7[,2]))) # parental average phenotype
  GEBVsParentPop7_table2=  mean(as.numeric(as.character(GEBVsParentPop7[,3]))) # parental average GEBV
  GEBVsParentPop8_table1=  mean(as.numeric(as.character(GEBVsParentPop8[,2]))) # parental average phenotype
  GEBVsParentPop8_table2=  mean(as.numeric(as.character(GEBVsParentPop8[,3]))) # parental average GEBV
  GEBVsParentPop9_table1=  mean(as.numeric(as.character(GEBVsParentPop9[,2]))) # parental average phenotype
  GEBVsParentPop9_table2=  mean(as.numeric(as.character(GEBVsParentPop9[,3]))) # parental average GEBV
  
  All=rbind(GEBVsParentPop1_table2,GEBVsParentPop2_table2,GEBVsParentPop3_table2,GEBVsParentPop4_table2,GEBVsParentPop5_table2,GEBVsParentPop6_table2,GEBVsParentPop7_table2,GEBVsParentPop8_table2,GEBVsParentPop9_table2)
  colnames(All)="GEBV"
  All2=rbind(GEBVsParentPop1_table1,GEBVsParentPop2_table1,GEBVsParentPop3_table1,GEBVsParentPop4_table1,GEBVsParentPop5_table1,GEBVsParentPop6_table1,GEBVsParentPop7_table1,GEBVsParentPop8_table1,GEBVsParentPop9_table1)
  colnames(All2)="AvgParentalPheno"
  Merged= cbind(F1observations,All,All2)
  
  correlation=cor(as.numeric(as.character(Merged[,3])),as.numeric(as.character(Merged[,4]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_iSize_4Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_iSize_4Replicates",round,".txt",sep="")
  write.table(Merged,filename1,sep="\t",quote=F,row.names=F)
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
  GEBVsParentPop1=results3[[1]]
  GEBVsParentPop2=results3[[2]]
  GEBVsParentPop3=results3[[3]]
  GEBVsParentPop4=results3[[4]]
  GEBVsParentPop5=results3[[5]]
  GEBVsParentPop6=results3[[6]]
  GEBVsParentPop7=results3[[7]]
  GEBVsParentPop8=results3[[8]]
  GEBVsParentPop9=results3[[9]]
  
  GEBVsParentPop1_table1=  mean(as.numeric(as.character(GEBVsParentPop1[,2]))) # parental average phenotype
  GEBVsParentPop1_table2=  mean(as.numeric(as.character(GEBVsParentPop1[,3]))) # parental average GEBV
  GEBVsParentPop2_table1=  mean(as.numeric(as.character(GEBVsParentPop2[,2]))) # parental average phenotype
  GEBVsParentPop2_table2=  mean(as.numeric(as.character(GEBVsParentPop2[,3]))) # parental average GEBV
  GEBVsParentPop3_table1=  mean(as.numeric(as.character(GEBVsParentPop3[,2]))) # parental average phenotype
  GEBVsParentPop3_table2=  mean(as.numeric(as.character(GEBVsParentPop3[,3]))) # parental average GEBV
  GEBVsParentPop4_table1=  mean(as.numeric(as.character(GEBVsParentPop4[,2]))) # parental average phenotype
  GEBVsParentPop4_table2=  mean(as.numeric(as.character(GEBVsParentPop4[,3]))) # parental average GEBV
  GEBVsParentPop5_table1=  mean(as.numeric(as.character(GEBVsParentPop5[,2]))) # parental average phenotype
  GEBVsParentPop5_table2=  mean(as.numeric(as.character(GEBVsParentPop5[,3]))) # parental average GEBV
  GEBVsParentPop6_table1=  mean(as.numeric(as.character(GEBVsParentPop6[,2]))) # parental average phenotype
  GEBVsParentPop6_table2=  mean(as.numeric(as.character(GEBVsParentPop6[,3]))) # parental average GEBV
  GEBVsParentPop7_table1=  mean(as.numeric(as.character(GEBVsParentPop7[,2]))) # parental average phenotype
  GEBVsParentPop7_table2=  mean(as.numeric(as.character(GEBVsParentPop7[,3]))) # parental average GEBV
  GEBVsParentPop8_table1=  mean(as.numeric(as.character(GEBVsParentPop8[,2]))) # parental average phenotype
  GEBVsParentPop8_table2=  mean(as.numeric(as.character(GEBVsParentPop8[,3]))) # parental average GEBV
  GEBVsParentPop9_table1=  mean(as.numeric(as.character(GEBVsParentPop9[,2]))) # parental average phenotype
  GEBVsParentPop9_table2=  mean(as.numeric(as.character(GEBVsParentPop9[,3]))) # parental average GEBV
  
  All=rbind(GEBVsParentPop1_table2,GEBVsParentPop2_table2,GEBVsParentPop3_table2,GEBVsParentPop4_table2,GEBVsParentPop5_table2,GEBVsParentPop6_table2,GEBVsParentPop7_table2,GEBVsParentPop8_table2,GEBVsParentPop9_table2)
  colnames(All)="GEBV"
  All2=rbind(GEBVsParentPop1_table1,GEBVsParentPop2_table1,GEBVsParentPop3_table1,GEBVsParentPop4_table1,GEBVsParentPop5_table1,GEBVsParentPop6_table1,GEBVsParentPop7_table1,GEBVsParentPop8_table1,GEBVsParentPop9_table1)
  colnames(All2)="AvgParentalPheno"
  Merged= cbind(F1observations,All,All2)
  
  correlation=cor(as.numeric(as.character(Merged[,3])),as.numeric(as.character(Merged[,4]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_iSize_3Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_iSize_3Replicates",round,".txt",sep="")
  write.table(Merged,filename1,sep="\t",quote=F,row.names=F)
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
  GEBVsParentPop1=results2[[1]]
  GEBVsParentPop2=results2[[2]]
  GEBVsParentPop3=results2[[3]]
  GEBVsParentPop4=results2[[4]]
  GEBVsParentPop5=results2[[5]]
  GEBVsParentPop6=results2[[6]]
  GEBVsParentPop7=results2[[7]]
  GEBVsParentPop8=results2[[8]]
  GEBVsParentPop9=results2[[9]]
  
  GEBVsParentPop1_table1=  mean(as.numeric(as.character(GEBVsParentPop1[,2]))) # parental average phenotype
  GEBVsParentPop1_table2=  mean(as.numeric(as.character(GEBVsParentPop1[,3]))) # parental average GEBV
  GEBVsParentPop2_table1=  mean(as.numeric(as.character(GEBVsParentPop2[,2]))) # parental average phenotype
  GEBVsParentPop2_table2=  mean(as.numeric(as.character(GEBVsParentPop2[,3]))) # parental average GEBV
  GEBVsParentPop3_table1=  mean(as.numeric(as.character(GEBVsParentPop3[,2]))) # parental average phenotype
  GEBVsParentPop3_table2=  mean(as.numeric(as.character(GEBVsParentPop3[,3]))) # parental average GEBV
  GEBVsParentPop4_table1=  mean(as.numeric(as.character(GEBVsParentPop4[,2]))) # parental average phenotype
  GEBVsParentPop4_table2=  mean(as.numeric(as.character(GEBVsParentPop4[,3]))) # parental average GEBV
  GEBVsParentPop5_table1=  mean(as.numeric(as.character(GEBVsParentPop5[,2]))) # parental average phenotype
  GEBVsParentPop5_table2=  mean(as.numeric(as.character(GEBVsParentPop5[,3]))) # parental average GEBV
  GEBVsParentPop6_table1=  mean(as.numeric(as.character(GEBVsParentPop6[,2]))) # parental average phenotype
  GEBVsParentPop6_table2=  mean(as.numeric(as.character(GEBVsParentPop6[,3]))) # parental average GEBV
  GEBVsParentPop7_table1=  mean(as.numeric(as.character(GEBVsParentPop7[,2]))) # parental average phenotype
  GEBVsParentPop7_table2=  mean(as.numeric(as.character(GEBVsParentPop7[,3]))) # parental average GEBV
  GEBVsParentPop8_table1=  mean(as.numeric(as.character(GEBVsParentPop8[,2]))) # parental average phenotype
  GEBVsParentPop8_table2=  mean(as.numeric(as.character(GEBVsParentPop8[,3]))) # parental average GEBV
  GEBVsParentPop9_table1=  mean(as.numeric(as.character(GEBVsParentPop9[,2]))) # parental average phenotype
  GEBVsParentPop9_table2=  mean(as.numeric(as.character(GEBVsParentPop9[,3]))) # parental average GEBV
  
  All=rbind(GEBVsParentPop1_table2,GEBVsParentPop2_table2,GEBVsParentPop3_table2,GEBVsParentPop4_table2,GEBVsParentPop5_table2,GEBVsParentPop6_table2,GEBVsParentPop7_table2,GEBVsParentPop8_table2,GEBVsParentPop9_table2)
  colnames(All)="GEBV"
  All2=rbind(GEBVsParentPop1_table1,GEBVsParentPop2_table1,GEBVsParentPop3_table1,GEBVsParentPop4_table1,GEBVsParentPop5_table1,GEBVsParentPop6_table1,GEBVsParentPop7_table1,GEBVsParentPop8_table1,GEBVsParentPop9_table1)
  colnames(All2)="AvgParentalPheno"
  Merged= cbind(F1observations,All,All2)
  
  correlation=cor(as.numeric(as.character(Merged[,3])),as.numeric(as.character(Merged[,4]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_iSize_2Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_iSize_2Replicates",round,".txt",sep="")
  write.table(Merged,filename1,sep="\t",quote=F,row.names=F)
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
  GEBVsParentPop1=results1[[1]]
  GEBVsParentPop2=results1[[2]]
  GEBVsParentPop3=results1[[3]]
  GEBVsParentPop4=results1[[4]]
  GEBVsParentPop5=results1[[5]]
  GEBVsParentPop6=results1[[6]]
  GEBVsParentPop7=results1[[7]]
  GEBVsParentPop8=results1[[8]]
  GEBVsParentPop9=results1[[9]]
  
  GEBVsParentPop1_table1=  mean(as.numeric(as.character(GEBVsParentPop1[,2]))) # parental average phenotype
  GEBVsParentPop1_table2=  mean(as.numeric(as.character(GEBVsParentPop1[,3]))) # parental average GEBV
  GEBVsParentPop2_table1=  mean(as.numeric(as.character(GEBVsParentPop2[,2]))) # parental average phenotype
  GEBVsParentPop2_table2=  mean(as.numeric(as.character(GEBVsParentPop2[,3]))) # parental average GEBV
  GEBVsParentPop3_table1=  mean(as.numeric(as.character(GEBVsParentPop3[,2]))) # parental average phenotype
  GEBVsParentPop3_table2=  mean(as.numeric(as.character(GEBVsParentPop3[,3]))) # parental average GEBV
  GEBVsParentPop4_table1=  mean(as.numeric(as.character(GEBVsParentPop4[,2]))) # parental average phenotype
  GEBVsParentPop4_table2=  mean(as.numeric(as.character(GEBVsParentPop4[,3]))) # parental average GEBV
  GEBVsParentPop5_table1=  mean(as.numeric(as.character(GEBVsParentPop5[,2]))) # parental average phenotype
  GEBVsParentPop5_table2=  mean(as.numeric(as.character(GEBVsParentPop5[,3]))) # parental average GEBV
  GEBVsParentPop6_table1=  mean(as.numeric(as.character(GEBVsParentPop6[,2]))) # parental average phenotype
  GEBVsParentPop6_table2=  mean(as.numeric(as.character(GEBVsParentPop6[,3]))) # parental average GEBV
  GEBVsParentPop7_table1=  mean(as.numeric(as.character(GEBVsParentPop7[,2]))) # parental average phenotype
  GEBVsParentPop7_table2=  mean(as.numeric(as.character(GEBVsParentPop7[,3]))) # parental average GEBV
  GEBVsParentPop8_table1=  mean(as.numeric(as.character(GEBVsParentPop8[,2]))) # parental average phenotype
  GEBVsParentPop8_table2=  mean(as.numeric(as.character(GEBVsParentPop8[,3]))) # parental average GEBV
  GEBVsParentPop9_table1=  mean(as.numeric(as.character(GEBVsParentPop9[,2]))) # parental average phenotype
  GEBVsParentPop9_table2=  mean(as.numeric(as.character(GEBVsParentPop9[,3]))) # parental average GEBV
  
  All=rbind(GEBVsParentPop1_table2,GEBVsParentPop2_table2,GEBVsParentPop3_table2,GEBVsParentPop4_table2,GEBVsParentPop5_table2,GEBVsParentPop6_table2,GEBVsParentPop7_table2,GEBVsParentPop8_table2,GEBVsParentPop9_table2)
  colnames(All)="GEBV"
  All2=rbind(GEBVsParentPop1_table1,GEBVsParentPop2_table1,GEBVsParentPop3_table1,GEBVsParentPop4_table1,GEBVsParentPop5_table1,GEBVsParentPop6_table1,GEBVsParentPop7_table1,GEBVsParentPop8_table1,GEBVsParentPop9_table1)
  colnames(All2)="AvgParentalPheno"
  Merged= cbind(F1observations,All,All2)
  
  correlation=cor(as.numeric(as.character(Merged[,3])),as.numeric(as.character(Merged[,4]))) #means of replicates
  correlation
  
  filename=paste("Correlation_GBLUP_iSize_1Replicates",round,".txt",sep="")
  write.table(correlation,filename,sep="\t",quote=F,row.names=F,col.names=F)
  
  filename1=paste("Predictions_GBLUP_iSize_1Replicates",round,".txt",sep="")
  write.table(Merged,filename1,sep="\t",quote=F,row.names=F)
}  

