######################################
# This is a combined script of the F1
# prediction script and the replicate 
# reduction script, allowing us to 
# predict the effects replicates have
# on the ability to succesfully predict
# the performance of F1 populations
# based on their parents
######################################

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
  set.seed(NULL)
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

# Divide into F1 parental groups
{
  parentpop1names=c("Aoost_01","Aoost_08","Aoost_09","Banna_02","Banna_03","Banna_07")
  parentpop2names=c("Rbani_02","Rbani_04","Rbani_08","Sster_06","Sster_07","Sster_09")
  parentpop3names=c("Aoost_03","Aoost_04","Aoost_06","Banna_03","Banna_07","Aaran_05","Aaran_07")
  parentpop4names=c("Banna_02","Banna_06","Banna_09","Aaran_01","Aaran_04","Aaran_06")
  parentpop5names=c("Kdike_03","Kdike_04","Kdike_10","Rling_07","Rling_08","Rling_10")
  parentpop6names=c("Aearl_08","Ccyma_03","Llanc_06","Aaran_08")
  parentpop7names=c("Aearl_05","Clfin_02","Ctain_05","Mrida_04")
  parentpop8names=c("Clfin_03","Ctain_05","Volin_01","Aaran_04")
  parentpop9names=c("Aoost_02","Ilona_09","Llanc_09","Sster_01")
  parentpop10names=c("Ilona_05","Kdike_09","Llanc_09","Aalon_03")
  parentpop11names=c("Ancor_10","Borek_06","Ctain_09","Rbani_02")
  parentpop12names=c("Ancor_04","Aoost_10","Clfin_08","Kdike_08")
  parentpop13names=c("Aoost_01","Aoost_08","Banna_02","Rbani_02","Sster_01","Sster_06")
  parentpop14names=c("Banna_02","Banna_09","Rbani_02","Rbani_07","Rbani_08")
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
  
  # Calculate a phenotype corrected for all fixed effects
  {
    Correctedforallfixed <- lmer(growth_per_day ~ factor(Round) + factor(NS) + factor(EW) + factor(Rhizobium) + inoculation_date + (1|Clover), data=dataframe) 
    summary(Correctedforallfixed)
    matrixOfeffects=model.matrix( ~ factor(Round) + factor(NS) + factor(EW)  + factor(Rhizobium) + inoculation_date, data=dataframe)
    column_to_remove=which(colnames(matrixOfeffects) %in% names(fixef(Correctedforallfixed))==F)
    matrixOfeffects_new=matrixOfeffects[,-column_to_remove]
    
    ycorr <- dataframe$growth_per_day - matrixOfeffects_new %*% fixef(Correctedforallfixed)
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
  
  parentpop1names_idx=which(dataframe$Clovershort%in%parentpop1names)
  parentpop2names_idx=which(dataframe$Clovershort%in%parentpop2names)
  parentpop3names_idx=which(dataframe$Clovershort%in%parentpop3names)
  parentpop4names_idx=which(dataframe$Clovershort%in%parentpop4names)
  parentpop5names_idx=which(dataframe$Clovershort%in%parentpop5names)
  parentpop6names_idx=which(dataframe$Clovershort%in%parentpop6names)
  parentpop7names_idx=which(dataframe$Clovershort%in%parentpop7names)
  parentpop8names_idx=which(dataframe$Clovershort%in%parentpop8names)
  parentpop9names_idx=which(dataframe$Clovershort%in%parentpop9names)
  parentpop10names_idx=which(dataframe$Clovershort%in%parentpop10names)
  parentpop11names_idx=which(dataframe$Clovershort%in%parentpop11names)
  parentpop12names_idx=which(dataframe$Clovershort%in%parentpop12names)
  parentpop13names_idx=which(dataframe$Clovershort%in%parentpop13names)
  parentpop14names_idx=which(dataframe$Clovershort%in%parentpop14names)
  
  tests=list(parentpop1names_idx,parentpop2names_idx,parentpop3names_idx,parentpop4names_idx,parentpop5names_idx,parentpop6names_idx,parentpop7names_idx,parentpop8names_idx,parentpop9names_idx,parentpop10names_idx,parentpop11names_idx,parentpop12names_idx,parentpop13names_idx,parentpop14names_idx)

  return(tests)
}

GP_GBLUP<-function(testpop){
  
  ################  ################  ################  ################
  ##start by estimating GEBVs for training population individuals 
  ################  ################  ################  ################
  
  d6_training=dataframe[-testpop,] # limit the dataframe to only the individuals allowed for training the model
  d6_training_ready=na.omit(d6_training, cols = c("growth_per_day")) # remember that gpd na inidividuals should be removed whether or not they are in the training pop or not
  
  ind_not_in_train=dataframe$Clovershort[testpop]
  IndividualsToRemoveGRM=which(colnames(GRM1) %in% ind_not_in_train)
  
  GRM_trn = GRM1[-IndividualsToRemoveGRM,-IndividualsToRemoveGRM]
  
  GcloverReps_trn = GcloverReps[-testpop,-testpop] # a GRM for the individuals used to predict
  CloverIndep_trn= CloverIndep[-testpop,-testpop] # another GRM for the individuals used to predict
  
  # Run the GBLUP model on full training population to extract GEBVs
  yNA=d6_training_ready$growth_per_day
  
  fixedmod=model.matrix(~factor(d6_training_ready$Round)+factor(d6_training_ready$NS)+factor(d6_training_ready$EW)+factor(d6_training_ready$Rhizobium)+factor(d6_training_ready$inoculation_date))
  ETA=list(list(K=GcloverReps_trn,model="RKHS"),list(K=CloverIndep_trn,model="RKHS"),list(X=fixedmod,model="FIXED"))
  GBLUP=BGLR(y=yNA,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("GBLUP",round))
  matrix=cbind(as.character(d6_training_ready$Clovershort),as.numeric(as.character(d6_training_ready$growth_per_day)),as.numeric(d6_training_ready$CorrectedPheno),as.numeric(GBLUP$ETA[[1]]$u),as.numeric(GBLUP$ETA[[2]]$u),(as.numeric(GBLUP$ETA[[1]]$u)+as.numeric(GBLUP$ETA[[2]]$u)))
  colnames(matrix)=c("ID", "gpd_NoCor","Observed", "1st GEBV contribution","2nd GEBV contribution (clover independent)","totalGEBV")
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
  matrix1=cbind(as.character(dataframe$Clovershort[testpop]),as.numeric(dataframe$growth_per_day[testpop]),as.numeric(dataframe$CorrectedPheno[testpop]),as.numeric(as.character(GEBVpred)))
  colnames(matrix1)=c("ID", "gpd_NoCor","Observed", "GEBV")
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
  GEBVsParentPop1=results10[[1]]
  GEBVsParentPop2=results10[[2]]
  GEBVsParentPop3=results10[[3]]
  GEBVsParentPop4=results10[[4]]
  GEBVsParentPop5=results10[[5]]
  GEBVsParentPop6=results10[[6]]
  GEBVsParentPop7=results10[[7]]
  GEBVsParentPop8=results10[[8]]
  GEBVsParentPop9=results10[[9]]
  GEBVsParentPop10=results10[[10]]
  GEBVsParentPop11=results10[[11]]
  GEBVsParentPop12=results10[[12]]
  GEBVsParentPop13=results10[[13]]
  GEBVsParentPop14=results10[[14]]

  # pop 1
  GEBVsParentPop1_table1=aggregate(as.numeric(GEBVsParentPop1[, 3]), list(GEBVsParentPop1[,1]), mean) #observed a mean for each individual
  GEBVsParentPop1_table2=aggregate(as.numeric(GEBVsParentPop1[, 4]), list(GEBVsParentPop1[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined1=cbind(GEBVsParentPop1_table1[,1],GEBVsParentPop1_table1[,2],GEBVsParentPop1_table2[,2])
  Pre_Combined1=as.data.frame(Pre_Combined1)
  colnames(Pre_Combined1)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 2
  GEBVsParentPop2_table1=aggregate(as.numeric(GEBVsParentPop2[, 3]), list(GEBVsParentPop2[,1]), mean) #observed a mean for each individual
  GEBVsParentPop2_table2=aggregate(as.numeric(GEBVsParentPop2[,4]), list(GEBVsParentPop2[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined2=cbind(GEBVsParentPop2_table1[,1],GEBVsParentPop2_table1[,2],GEBVsParentPop2_table2[,2])
  Pre_Combined2=as.data.frame(Pre_Combined2)
  colnames(Pre_Combined2)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 3
  GEBVsParentPop3_table1=aggregate(as.numeric(GEBVsParentPop3[, 3]), list(GEBVsParentPop3[,1]), mean) #observed a mean for each individual
  GEBVsParentPop3_table2=aggregate(as.numeric(GEBVsParentPop3[, 4]), list(GEBVsParentPop3[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined3=cbind(GEBVsParentPop3_table1[,1],GEBVsParentPop3_table1[,2],GEBVsParentPop3_table2[,2])
  Pre_Combined3=as.data.frame(Pre_Combined3)
  colnames(Pre_Combined3)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  #  pop 4
  GEBVsParentPop4_table1=aggregate(as.numeric(GEBVsParentPop4[, 3]), list(GEBVsParentPop4[,1]), mean) #observed a mean for each individual
  GEBVsParentPop4_table2=aggregate(as.numeric(GEBVsParentPop4[, 4]), list(GEBVsParentPop4[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined4=cbind(GEBVsParentPop4_table1[,1],GEBVsParentPop4_table1[,2],GEBVsParentPop4_table2[,2])
  Pre_Combined4=as.data.frame(Pre_Combined4)
  colnames(Pre_Combined4)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 5
  GEBVsParentPop5_table1=aggregate(as.numeric(GEBVsParentPop5[, 3]), list(GEBVsParentPop5[,1]), mean) #observed a mean for each individual
  GEBVsParentPop5_table2=aggregate(as.numeric(GEBVsParentPop5[, 4]), list(GEBVsParentPop5[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined5=cbind(GEBVsParentPop5_table1[,1],GEBVsParentPop5_table1[,2],GEBVsParentPop5_table2[,2])
  Pre_Combined5=as.data.frame(Pre_Combined5)
  colnames(Pre_Combined5)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 6
  GEBVsParentPop6_table1=aggregate(as.numeric(GEBVsParentPop6[, 3]), list(GEBVsParentPop6[,1]), mean) #observed a mean for each individual
  GEBVsParentPop6_table2=aggregate(as.numeric(GEBVsParentPop6[, 4]), list(GEBVsParentPop6[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined6=cbind(GEBVsParentPop6_table1[,1],GEBVsParentPop6_table1[,2],GEBVsParentPop6_table2[,2])
  Pre_Combined6=as.data.frame(Pre_Combined6)
  colnames(Pre_Combined6)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 7
  GEBVsParentPop7_table1=aggregate(as.numeric(GEBVsParentPop7[, 3]), list(GEBVsParentPop7[,1]), mean) #observed a mean for each individual
  GEBVsParentPop7_table2=aggregate(as.numeric(GEBVsParentPop7[, 4]), list(GEBVsParentPop7[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined7=cbind(GEBVsParentPop7_table1[,1],GEBVsParentPop7_table1[,2],GEBVsParentPop7_table2[,2])
  Pre_Combined7=as.data.frame(Pre_Combined7)
  colnames(Pre_Combined7)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 8
  GEBVsParentPop8_table1=aggregate(as.numeric(GEBVsParentPop8[, 3]), list(GEBVsParentPop8[,1]), mean) #observed a mean for each individual
  GEBVsParentPop8_table2=aggregate(as.numeric(GEBVsParentPop8[, 4]), list(GEBVsParentPop8[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined8=cbind(GEBVsParentPop8_table1[,1],GEBVsParentPop8_table1[,2],GEBVsParentPop8_table2[,2])
  Pre_Combined8=as.data.frame(Pre_Combined8)
  colnames(Pre_Combined8)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 9
  GEBVsParentPop9_table1=aggregate(as.numeric(GEBVsParentPop9[, 3]), list(GEBVsParentPop9[,1]), mean) #observed a mean for each individual
  GEBVsParentPop9_table2=aggregate(as.numeric(GEBVsParentPop9[, 4]), list(GEBVsParentPop9[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined9=cbind(GEBVsParentPop9_table1[,1],GEBVsParentPop9_table1[,2],GEBVsParentPop9_table2[,2])
  Pre_Combined9=as.data.frame(Pre_Combined9)
  colnames(Pre_Combined9)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 10
  GEBVsParentPop10_table1=aggregate(as.numeric(GEBVsParentPop10[, 3]), list(GEBVsParentPop10[,1]), mean) #observed a mean for each individual
  GEBVsParentPop10_table2=aggregate(as.numeric(GEBVsParentPop10[, 4]), list(GEBVsParentPop10[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined10=cbind(GEBVsParentPop10_table1[,1],GEBVsParentPop10_table1[,2],GEBVsParentPop10_table2[,2])
  Pre_Combined10=as.data.frame(Pre_Combined10)
  colnames(Pre_Combined10)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 11
  GEBVsParentPop11_table1=aggregate(as.numeric(GEBVsParentPop11[, 3]), list(GEBVsParentPop11[,1]), mean) #observed a mean for each individual
  GEBVsParentPop11_table2=aggregate(as.numeric(GEBVsParentPop11[, 4]), list(GEBVsParentPop11[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined11=cbind(GEBVsParentPop11_table1[,1],GEBVsParentPop11_table1[,2],GEBVsParentPop11_table2[,2])
  Pre_Combined11=as.data.frame(Pre_Combined11)
  colnames(Pre_Combined11)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 12
  GEBVsParentPop12_table1=aggregate(as.numeric(GEBVsParentPop12[, 3]), list(GEBVsParentPop12[,1]), mean) #observed a mean for each individual
  GEBVsParentPop12_table2=aggregate(as.numeric(GEBVsParentPop12[, 4]), list(GEBVsParentPop12[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined12=cbind(GEBVsParentPop12_table1[,1],GEBVsParentPop12_table1[,2],GEBVsParentPop12_table2[,2])
  colnames(Pre_Combined12)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 13
  GEBVsParentPop13_table1=aggregate(as.numeric(GEBVsParentPop13[, 3]), list(GEBVsParentPop13[,1]), mean) #observed a mean for each individual
  GEBVsParentPop13_table2=aggregate(as.numeric(GEBVsParentPop13[, 4]), list(GEBVsParentPop13[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined13=cbind(GEBVsParentPop13_table1[,1],GEBVsParentPop13_table1[,2],GEBVsParentPop13_table2[,2])
  Pre_Combined13=as.data.frame(Pre_Combined13)
  colnames(Pre_Combined13)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 14
  GEBVsParentPop14_table1=aggregate(as.numeric(GEBVsParentPop14[, 3]), list(GEBVsParentPop14[,1]), mean) #observed a mean for each individual
  GEBVsParentPop14_table2=aggregate(as.numeric(GEBVsParentPop14[, 4]), list(GEBVsParentPop14[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined14=cbind(GEBVsParentPop14_table1[,1],GEBVsParentPop14_table1[,2],GEBVsParentPop14_table2[,2])
  colnames(Pre_Combined14)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  Combined=rbind(Pre_Combined1,Pre_Combined2,Pre_Combined3,Pre_Combined4,Pre_Combined5,Pre_Combined6,Pre_Combined7,Pre_Combined8,Pre_Combined9,Pre_Combined10,Pre_Combined11,Pre_Combined12,Pre_Combined13,Pre_Combined14)
  filename1=paste("Predictions_GBLUP_GPD_10Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
  
  # Calculating parental average GEBV
  MeanParentalGEBVs=c(mean(GEBVsParentPop1_table2[,2]),mean(GEBVsParentPop2_table2[,2]),mean(GEBVsParentPop3_table2[,2]),mean(GEBVsParentPop4_table2[,2]),mean(GEBVsParentPop5_table2[,2]),mean(GEBVsParentPop6_table2[,2]),mean(GEBVsParentPop7_table2[,2]),mean(GEBVsParentPop8_table2[,2]),mean(GEBVsParentPop9_table2[,2]),mean(GEBVsParentPop10_table2[,2]),mean(GEBVsParentPop11_table2[,2]),mean(GEBVsParentPop12_table2[,2]),mean(GEBVsParentPop13_table2[,2]),mean(GEBVsParentPop14_table2[,2]))
  F1_prediction_table=cbind(F1observations$Population,MeanParentalGEBVs,F1observations$X20200226_dryweight,F1observations$X20200212_freshweight)
  F1_prediction_table=as.data.frame(F1_prediction_table)
  colnames(F1_prediction_table)=c("F1 population","Mean Parental predicted GEBV","F1 pop mean dryweight","F1 pop mean freshweight")
  filename2=paste("Correlation_Correctedphenotype_GBLUP_GPD_10Replicates",round,".txt",sep="")
  write.table(F1_prediction_table,filename2,sep="\t",quote=F)
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
  GEBVsParentPop10=results9[[10]]
  GEBVsParentPop11=results9[[11]]
  GEBVsParentPop12=results9[[12]]
  GEBVsParentPop13=results9[[13]]
  GEBVsParentPop14=results9[[14]]
  
  # pop 1
  GEBVsParentPop1_table1=aggregate(as.numeric(GEBVsParentPop1[, 3]), list(GEBVsParentPop1[,1]), mean) #observed a mean for each individual
  GEBVsParentPop1_table2=aggregate(as.numeric(GEBVsParentPop1[, 4]), list(GEBVsParentPop1[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined1=cbind(GEBVsParentPop1_table1[,1],GEBVsParentPop1_table1[,2],GEBVsParentPop1_table2[,2])
  Pre_Combined1=as.data.frame(Pre_Combined1)
  colnames(Pre_Combined1)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 2
  GEBVsParentPop2_table1=aggregate(as.numeric(GEBVsParentPop2[, 3]), list(GEBVsParentPop2[,1]), mean) #observed a mean for each individual
  GEBVsParentPop2_table2=aggregate(as.numeric(GEBVsParentPop2[,4]), list(GEBVsParentPop2[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined2=cbind(GEBVsParentPop2_table1[,1],GEBVsParentPop2_table1[,2],GEBVsParentPop2_table2[,2])
  Pre_Combined2=as.data.frame(Pre_Combined2)
  colnames(Pre_Combined2)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 3
  GEBVsParentPop3_table1=aggregate(as.numeric(GEBVsParentPop3[, 3]), list(GEBVsParentPop3[,1]), mean) #observed a mean for each individual
  GEBVsParentPop3_table2=aggregate(as.numeric(GEBVsParentPop3[, 4]), list(GEBVsParentPop3[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined3=cbind(GEBVsParentPop3_table1[,1],GEBVsParentPop3_table1[,2],GEBVsParentPop3_table2[,2])
  Pre_Combined3=as.data.frame(Pre_Combined3)
  colnames(Pre_Combined3)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  #  pop 4
  GEBVsParentPop4_table1=aggregate(as.numeric(GEBVsParentPop4[, 3]), list(GEBVsParentPop4[,1]), mean) #observed a mean for each individual
  GEBVsParentPop4_table2=aggregate(as.numeric(GEBVsParentPop4[, 4]), list(GEBVsParentPop4[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined4=cbind(GEBVsParentPop4_table1[,1],GEBVsParentPop4_table1[,2],GEBVsParentPop4_table2[,2])
  Pre_Combined4=as.data.frame(Pre_Combined4)
  colnames(Pre_Combined4)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 5
  GEBVsParentPop5_table1=aggregate(as.numeric(GEBVsParentPop5[, 3]), list(GEBVsParentPop5[,1]), mean) #observed a mean for each individual
  GEBVsParentPop5_table2=aggregate(as.numeric(GEBVsParentPop5[, 4]), list(GEBVsParentPop5[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined5=cbind(GEBVsParentPop5_table1[,1],GEBVsParentPop5_table1[,2],GEBVsParentPop5_table2[,2])
  Pre_Combined5=as.data.frame(Pre_Combined5)
  colnames(Pre_Combined5)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 6
  GEBVsParentPop6_table1=aggregate(as.numeric(GEBVsParentPop6[, 3]), list(GEBVsParentPop6[,1]), mean) #observed a mean for each individual
  GEBVsParentPop6_table2=aggregate(as.numeric(GEBVsParentPop6[, 4]), list(GEBVsParentPop6[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined6=cbind(GEBVsParentPop6_table1[,1],GEBVsParentPop6_table1[,2],GEBVsParentPop6_table2[,2])
  Pre_Combined6=as.data.frame(Pre_Combined6)
  colnames(Pre_Combined6)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 7
  GEBVsParentPop7_table1=aggregate(as.numeric(GEBVsParentPop7[, 3]), list(GEBVsParentPop7[,1]), mean) #observed a mean for each individual
  GEBVsParentPop7_table2=aggregate(as.numeric(GEBVsParentPop7[, 4]), list(GEBVsParentPop7[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined7=cbind(GEBVsParentPop7_table1[,1],GEBVsParentPop7_table1[,2],GEBVsParentPop7_table2[,2])
  Pre_Combined7=as.data.frame(Pre_Combined7)
  colnames(Pre_Combined7)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 8
  GEBVsParentPop8_table1=aggregate(as.numeric(GEBVsParentPop8[, 3]), list(GEBVsParentPop8[,1]), mean) #observed a mean for each individual
  GEBVsParentPop8_table2=aggregate(as.numeric(GEBVsParentPop8[, 4]), list(GEBVsParentPop8[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined8=cbind(GEBVsParentPop8_table1[,1],GEBVsParentPop8_table1[,2],GEBVsParentPop8_table2[,2])
  Pre_Combined8=as.data.frame(Pre_Combined8)
  colnames(Pre_Combined8)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 9
  GEBVsParentPop9_table1=aggregate(as.numeric(GEBVsParentPop9[, 3]), list(GEBVsParentPop9[,1]), mean) #observed a mean for each individual
  GEBVsParentPop9_table2=aggregate(as.numeric(GEBVsParentPop9[, 4]), list(GEBVsParentPop9[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined9=cbind(GEBVsParentPop9_table1[,1],GEBVsParentPop9_table1[,2],GEBVsParentPop9_table2[,2])
  Pre_Combined9=as.data.frame(Pre_Combined9)
  colnames(Pre_Combined9)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 10
  GEBVsParentPop10_table1=aggregate(as.numeric(GEBVsParentPop10[, 3]), list(GEBVsParentPop10[,1]), mean) #observed a mean for each individual
  GEBVsParentPop10_table2=aggregate(as.numeric(GEBVsParentPop10[, 4]), list(GEBVsParentPop10[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined10=cbind(GEBVsParentPop10_table1[,1],GEBVsParentPop10_table1[,2],GEBVsParentPop10_table2[,2])
  Pre_Combined10=as.data.frame(Pre_Combined10)
  colnames(Pre_Combined10)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 11
  GEBVsParentPop11_table1=aggregate(as.numeric(GEBVsParentPop11[, 3]), list(GEBVsParentPop11[,1]), mean) #observed a mean for each individual
  GEBVsParentPop11_table2=aggregate(as.numeric(GEBVsParentPop11[, 4]), list(GEBVsParentPop11[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined11=cbind(GEBVsParentPop11_table1[,1],GEBVsParentPop11_table1[,2],GEBVsParentPop11_table2[,2])
  Pre_Combined11=as.data.frame(Pre_Combined11)
  colnames(Pre_Combined11)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 12
  GEBVsParentPop12_table1=aggregate(as.numeric(GEBVsParentPop12[, 3]), list(GEBVsParentPop12[,1]), mean) #observed a mean for each individual
  GEBVsParentPop12_table2=aggregate(as.numeric(GEBVsParentPop12[, 4]), list(GEBVsParentPop12[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined12=cbind(GEBVsParentPop12_table1[,1],GEBVsParentPop12_table1[,2],GEBVsParentPop12_table2[,2])
  colnames(Pre_Combined12)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 13
  GEBVsParentPop13_table1=aggregate(as.numeric(GEBVsParentPop13[, 3]), list(GEBVsParentPop13[,1]), mean) #observed a mean for each individual
  GEBVsParentPop13_table2=aggregate(as.numeric(GEBVsParentPop13[, 4]), list(GEBVsParentPop13[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined13=cbind(GEBVsParentPop13_table1[,1],GEBVsParentPop13_table1[,2],GEBVsParentPop13_table2[,2])
  Pre_Combined13=as.data.frame(Pre_Combined13)
  colnames(Pre_Combined13)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 14
  GEBVsParentPop14_table1=aggregate(as.numeric(GEBVsParentPop14[, 3]), list(GEBVsParentPop14[,1]), mean) #observed a mean for each individual
  GEBVsParentPop14_table2=aggregate(as.numeric(GEBVsParentPop14[, 4]), list(GEBVsParentPop14[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined14=cbind(GEBVsParentPop14_table1[,1],GEBVsParentPop14_table1[,2],GEBVsParentPop14_table2[,2])
  colnames(Pre_Combined14)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  Combined=rbind(Pre_Combined1,Pre_Combined2,Pre_Combined3,Pre_Combined4,Pre_Combined5,Pre_Combined6,Pre_Combined7,Pre_Combined8,Pre_Combined9,Pre_Combined10,Pre_Combined11,Pre_Combined12,Pre_Combined13,Pre_Combined14)
  filename1=paste("Predictions_GBLUP_GPD_9Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
  
  # Calculating parental average GEBV
  MeanParentalGEBVs=c(mean(GEBVsParentPop1_table2[,2]),mean(GEBVsParentPop2_table2[,2]),mean(GEBVsParentPop3_table2[,2]),mean(GEBVsParentPop4_table2[,2]),mean(GEBVsParentPop5_table2[,2]),mean(GEBVsParentPop6_table2[,2]),mean(GEBVsParentPop7_table2[,2]),mean(GEBVsParentPop8_table2[,2]),mean(GEBVsParentPop9_table2[,2]),mean(GEBVsParentPop10_table2[,2]),mean(GEBVsParentPop11_table2[,2]),mean(GEBVsParentPop12_table2[,2]),mean(GEBVsParentPop13_table2[,2]),mean(GEBVsParentPop14_table2[,2]))
  F1_prediction_table=cbind(F1observations$Population,MeanParentalGEBVs,F1observations$X20200226_dryweight,F1observations$X20200212_freshweight)
  F1_prediction_table=as.data.frame(F1_prediction_table)
  colnames(F1_prediction_table)=c("F1 population","Mean Parental predicted GEBV","F1 pop mean dryweight","F1 pop mean freshweight")
  filename2=paste("Correlation_Correctedphenotype_GBLUP_GPD_9Replicates",round,".txt",sep="")
  write.table(F1_prediction_table,filename2,sep="\t",quote=F)
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
  GEBVsParentPop1=results8[[1]]
  GEBVsParentPop2=results8[[2]]
  GEBVsParentPop3=results8[[3]]
  GEBVsParentPop4=results8[[4]]
  GEBVsParentPop5=results8[[5]]
  GEBVsParentPop6=results8[[6]]
  GEBVsParentPop7=results8[[7]]
  GEBVsParentPop8=results8[[8]]
  GEBVsParentPop9=results8[[9]]
  GEBVsParentPop10=results8[[10]]
  GEBVsParentPop11=results8[[11]]
  GEBVsParentPop12=results8[[12]]
  GEBVsParentPop13=results8[[13]]
  GEBVsParentPop14=results8[[14]]
  
  # pop 1
  GEBVsParentPop1_table1=aggregate(as.numeric(GEBVsParentPop1[, 3]), list(GEBVsParentPop1[,1]), mean) #observed a mean for each individual
  GEBVsParentPop1_table2=aggregate(as.numeric(GEBVsParentPop1[, 4]), list(GEBVsParentPop1[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined1=cbind(GEBVsParentPop1_table1[,1],GEBVsParentPop1_table1[,2],GEBVsParentPop1_table2[,2])
  Pre_Combined1=as.data.frame(Pre_Combined1)
  colnames(Pre_Combined1)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 2
  GEBVsParentPop2_table1=aggregate(as.numeric(GEBVsParentPop2[, 3]), list(GEBVsParentPop2[,1]), mean) #observed a mean for each individual
  GEBVsParentPop2_table2=aggregate(as.numeric(GEBVsParentPop2[,4]), list(GEBVsParentPop2[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined2=cbind(GEBVsParentPop2_table1[,1],GEBVsParentPop2_table1[,2],GEBVsParentPop2_table2[,2])
  Pre_Combined2=as.data.frame(Pre_Combined2)
  colnames(Pre_Combined2)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 3
  GEBVsParentPop3_table1=aggregate(as.numeric(GEBVsParentPop3[, 3]), list(GEBVsParentPop3[,1]), mean) #observed a mean for each individual
  GEBVsParentPop3_table2=aggregate(as.numeric(GEBVsParentPop3[, 4]), list(GEBVsParentPop3[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined3=cbind(GEBVsParentPop3_table1[,1],GEBVsParentPop3_table1[,2],GEBVsParentPop3_table2[,2])
  Pre_Combined3=as.data.frame(Pre_Combined3)
  colnames(Pre_Combined3)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  #  pop 4
  GEBVsParentPop4_table1=aggregate(as.numeric(GEBVsParentPop4[, 3]), list(GEBVsParentPop4[,1]), mean) #observed a mean for each individual
  GEBVsParentPop4_table2=aggregate(as.numeric(GEBVsParentPop4[, 4]), list(GEBVsParentPop4[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined4=cbind(GEBVsParentPop4_table1[,1],GEBVsParentPop4_table1[,2],GEBVsParentPop4_table2[,2])
  Pre_Combined4=as.data.frame(Pre_Combined4)
  colnames(Pre_Combined4)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 5
  GEBVsParentPop5_table1=aggregate(as.numeric(GEBVsParentPop5[, 3]), list(GEBVsParentPop5[,1]), mean) #observed a mean for each individual
  GEBVsParentPop5_table2=aggregate(as.numeric(GEBVsParentPop5[, 4]), list(GEBVsParentPop5[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined5=cbind(GEBVsParentPop5_table1[,1],GEBVsParentPop5_table1[,2],GEBVsParentPop5_table2[,2])
  Pre_Combined5=as.data.frame(Pre_Combined5)
  colnames(Pre_Combined5)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 6
  GEBVsParentPop6_table1=aggregate(as.numeric(GEBVsParentPop6[, 3]), list(GEBVsParentPop6[,1]), mean) #observed a mean for each individual
  GEBVsParentPop6_table2=aggregate(as.numeric(GEBVsParentPop6[, 4]), list(GEBVsParentPop6[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined6=cbind(GEBVsParentPop6_table1[,1],GEBVsParentPop6_table1[,2],GEBVsParentPop6_table2[,2])
  Pre_Combined6=as.data.frame(Pre_Combined6)
  colnames(Pre_Combined6)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 7
  GEBVsParentPop7_table1=aggregate(as.numeric(GEBVsParentPop7[, 3]), list(GEBVsParentPop7[,1]), mean) #observed a mean for each individual
  GEBVsParentPop7_table2=aggregate(as.numeric(GEBVsParentPop7[, 4]), list(GEBVsParentPop7[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined7=cbind(GEBVsParentPop7_table1[,1],GEBVsParentPop7_table1[,2],GEBVsParentPop7_table2[,2])
  Pre_Combined7=as.data.frame(Pre_Combined7)
  colnames(Pre_Combined7)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 8
  GEBVsParentPop8_table1=aggregate(as.numeric(GEBVsParentPop8[, 3]), list(GEBVsParentPop8[,1]), mean) #observed a mean for each individual
  GEBVsParentPop8_table2=aggregate(as.numeric(GEBVsParentPop8[, 4]), list(GEBVsParentPop8[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined8=cbind(GEBVsParentPop8_table1[,1],GEBVsParentPop8_table1[,2],GEBVsParentPop8_table2[,2])
  Pre_Combined8=as.data.frame(Pre_Combined8)
  colnames(Pre_Combined8)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 9
  GEBVsParentPop9_table1=aggregate(as.numeric(GEBVsParentPop9[, 3]), list(GEBVsParentPop9[,1]), mean) #observed a mean for each individual
  GEBVsParentPop9_table2=aggregate(as.numeric(GEBVsParentPop9[, 4]), list(GEBVsParentPop9[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined9=cbind(GEBVsParentPop9_table1[,1],GEBVsParentPop9_table1[,2],GEBVsParentPop9_table2[,2])
  Pre_Combined9=as.data.frame(Pre_Combined9)
  colnames(Pre_Combined9)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 10
  GEBVsParentPop10_table1=aggregate(as.numeric(GEBVsParentPop10[, 3]), list(GEBVsParentPop10[,1]), mean) #observed a mean for each individual
  GEBVsParentPop10_table2=aggregate(as.numeric(GEBVsParentPop10[, 4]), list(GEBVsParentPop10[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined10=cbind(GEBVsParentPop10_table1[,1],GEBVsParentPop10_table1[,2],GEBVsParentPop10_table2[,2])
  Pre_Combined10=as.data.frame(Pre_Combined10)
  colnames(Pre_Combined10)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 11
  GEBVsParentPop11_table1=aggregate(as.numeric(GEBVsParentPop11[, 3]), list(GEBVsParentPop11[,1]), mean) #observed a mean for each individual
  GEBVsParentPop11_table2=aggregate(as.numeric(GEBVsParentPop11[, 4]), list(GEBVsParentPop11[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined11=cbind(GEBVsParentPop11_table1[,1],GEBVsParentPop11_table1[,2],GEBVsParentPop11_table2[,2])
  Pre_Combined11=as.data.frame(Pre_Combined11)
  colnames(Pre_Combined11)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 12
  GEBVsParentPop12_table1=aggregate(as.numeric(GEBVsParentPop12[, 3]), list(GEBVsParentPop12[,1]), mean) #observed a mean for each individual
  GEBVsParentPop12_table2=aggregate(as.numeric(GEBVsParentPop12[, 4]), list(GEBVsParentPop12[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined12=cbind(GEBVsParentPop12_table1[,1],GEBVsParentPop12_table1[,2],GEBVsParentPop12_table2[,2])
  colnames(Pre_Combined12)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 13
  GEBVsParentPop13_table1=aggregate(as.numeric(GEBVsParentPop13[, 3]), list(GEBVsParentPop13[,1]), mean) #observed a mean for each individual
  GEBVsParentPop13_table2=aggregate(as.numeric(GEBVsParentPop13[, 4]), list(GEBVsParentPop13[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined13=cbind(GEBVsParentPop13_table1[,1],GEBVsParentPop13_table1[,2],GEBVsParentPop13_table2[,2])
  Pre_Combined13=as.data.frame(Pre_Combined13)
  colnames(Pre_Combined13)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 14
  GEBVsParentPop14_table1=aggregate(as.numeric(GEBVsParentPop14[, 3]), list(GEBVsParentPop14[,1]), mean) #observed a mean for each individual
  GEBVsParentPop14_table2=aggregate(as.numeric(GEBVsParentPop14[, 4]), list(GEBVsParentPop14[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined14=cbind(GEBVsParentPop14_table1[,1],GEBVsParentPop14_table1[,2],GEBVsParentPop14_table2[,2])
  colnames(Pre_Combined14)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  Combined=rbind(Pre_Combined1,Pre_Combined2,Pre_Combined3,Pre_Combined4,Pre_Combined5,Pre_Combined6,Pre_Combined7,Pre_Combined8,Pre_Combined9,Pre_Combined10,Pre_Combined11,Pre_Combined12,Pre_Combined13,Pre_Combined14)
  filename1=paste("Predictions_GBLUP_GPD_8Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
  
  # Calculating parental average GEBV
  MeanParentalGEBVs=c(mean(GEBVsParentPop1_table2[,2]),mean(GEBVsParentPop2_table2[,2]),mean(GEBVsParentPop3_table2[,2]),mean(GEBVsParentPop4_table2[,2]),mean(GEBVsParentPop5_table2[,2]),mean(GEBVsParentPop6_table2[,2]),mean(GEBVsParentPop7_table2[,2]),mean(GEBVsParentPop8_table2[,2]),mean(GEBVsParentPop9_table2[,2]),mean(GEBVsParentPop10_table2[,2]),mean(GEBVsParentPop11_table2[,2]),mean(GEBVsParentPop12_table2[,2]),mean(GEBVsParentPop13_table2[,2]),mean(GEBVsParentPop14_table2[,2]))
  F1_prediction_table=cbind(F1observations$Population,MeanParentalGEBVs,F1observations$X20200226_dryweight,F1observations$X20200212_freshweight)
  F1_prediction_table=as.data.frame(F1_prediction_table)
  colnames(F1_prediction_table)=c("F1 population","Mean Parental predicted GEBV","F1 pop mean dryweight","F1 pop mean freshweight")
  filename2=paste("Correlation_Correctedphenotype_GBLUP_GPD_8Replicates",round,".txt",sep="")
  write.table(F1_prediction_table,filename2,sep="\t",quote=F)
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
  GEBVsParentPop1=results7[[1]]
  GEBVsParentPop2=results7[[2]]
  GEBVsParentPop3=results7[[3]]
  GEBVsParentPop4=results7[[4]]
  GEBVsParentPop5=results7[[5]]
  GEBVsParentPop6=results7[[6]]
  GEBVsParentPop7=results7[[7]]
  GEBVsParentPop8=results7[[8]]
  GEBVsParentPop9=results7[[9]]
  GEBVsParentPop10=results7[[10]]
  GEBVsParentPop11=results7[[11]]
  GEBVsParentPop12=results7[[12]]
  GEBVsParentPop13=results7[[13]]
  GEBVsParentPop14=results7[[14]]
  
  # pop 1
  GEBVsParentPop1_table1=aggregate(as.numeric(GEBVsParentPop1[, 3]), list(GEBVsParentPop1[,1]), mean) #observed a mean for each individual
  GEBVsParentPop1_table2=aggregate(as.numeric(GEBVsParentPop1[, 4]), list(GEBVsParentPop1[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined1=cbind(GEBVsParentPop1_table1[,1],GEBVsParentPop1_table1[,2],GEBVsParentPop1_table2[,2])
  Pre_Combined1=as.data.frame(Pre_Combined1)
  colnames(Pre_Combined1)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 2
  GEBVsParentPop2_table1=aggregate(as.numeric(GEBVsParentPop2[, 3]), list(GEBVsParentPop2[,1]), mean) #observed a mean for each individual
  GEBVsParentPop2_table2=aggregate(as.numeric(GEBVsParentPop2[,4]), list(GEBVsParentPop2[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined2=cbind(GEBVsParentPop2_table1[,1],GEBVsParentPop2_table1[,2],GEBVsParentPop2_table2[,2])
  Pre_Combined2=as.data.frame(Pre_Combined2)
  colnames(Pre_Combined2)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 3
  GEBVsParentPop3_table1=aggregate(as.numeric(GEBVsParentPop3[, 3]), list(GEBVsParentPop3[,1]), mean) #observed a mean for each individual
  GEBVsParentPop3_table2=aggregate(as.numeric(GEBVsParentPop3[, 4]), list(GEBVsParentPop3[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined3=cbind(GEBVsParentPop3_table1[,1],GEBVsParentPop3_table1[,2],GEBVsParentPop3_table2[,2])
  Pre_Combined3=as.data.frame(Pre_Combined3)
  colnames(Pre_Combined3)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  #  pop 4
  GEBVsParentPop4_table1=aggregate(as.numeric(GEBVsParentPop4[, 3]), list(GEBVsParentPop4[,1]), mean) #observed a mean for each individual
  GEBVsParentPop4_table2=aggregate(as.numeric(GEBVsParentPop4[, 4]), list(GEBVsParentPop4[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined4=cbind(GEBVsParentPop4_table1[,1],GEBVsParentPop4_table1[,2],GEBVsParentPop4_table2[,2])
  Pre_Combined4=as.data.frame(Pre_Combined4)
  colnames(Pre_Combined4)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 5
  GEBVsParentPop5_table1=aggregate(as.numeric(GEBVsParentPop5[, 3]), list(GEBVsParentPop5[,1]), mean) #observed a mean for each individual
  GEBVsParentPop5_table2=aggregate(as.numeric(GEBVsParentPop5[, 4]), list(GEBVsParentPop5[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined5=cbind(GEBVsParentPop5_table1[,1],GEBVsParentPop5_table1[,2],GEBVsParentPop5_table2[,2])
  Pre_Combined5=as.data.frame(Pre_Combined5)
  colnames(Pre_Combined5)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 6
  GEBVsParentPop6_table1=aggregate(as.numeric(GEBVsParentPop6[, 3]), list(GEBVsParentPop6[,1]), mean) #observed a mean for each individual
  GEBVsParentPop6_table2=aggregate(as.numeric(GEBVsParentPop6[, 4]), list(GEBVsParentPop6[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined6=cbind(GEBVsParentPop6_table1[,1],GEBVsParentPop6_table1[,2],GEBVsParentPop6_table2[,2])
  Pre_Combined6=as.data.frame(Pre_Combined6)
  colnames(Pre_Combined6)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 7
  GEBVsParentPop7_table1=aggregate(as.numeric(GEBVsParentPop7[, 3]), list(GEBVsParentPop7[,1]), mean) #observed a mean for each individual
  GEBVsParentPop7_table2=aggregate(as.numeric(GEBVsParentPop7[, 4]), list(GEBVsParentPop7[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined7=cbind(GEBVsParentPop7_table1[,1],GEBVsParentPop7_table1[,2],GEBVsParentPop7_table2[,2])
  Pre_Combined7=as.data.frame(Pre_Combined7)
  colnames(Pre_Combined7)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 8
  GEBVsParentPop8_table1=aggregate(as.numeric(GEBVsParentPop8[, 3]), list(GEBVsParentPop8[,1]), mean) #observed a mean for each individual
  GEBVsParentPop8_table2=aggregate(as.numeric(GEBVsParentPop8[, 4]), list(GEBVsParentPop8[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined8=cbind(GEBVsParentPop8_table1[,1],GEBVsParentPop8_table1[,2],GEBVsParentPop8_table2[,2])
  Pre_Combined8=as.data.frame(Pre_Combined8)
  colnames(Pre_Combined8)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 9
  GEBVsParentPop9_table1=aggregate(as.numeric(GEBVsParentPop9[, 3]), list(GEBVsParentPop9[,1]), mean) #observed a mean for each individual
  GEBVsParentPop9_table2=aggregate(as.numeric(GEBVsParentPop9[, 4]), list(GEBVsParentPop9[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined9=cbind(GEBVsParentPop9_table1[,1],GEBVsParentPop9_table1[,2],GEBVsParentPop9_table2[,2])
  Pre_Combined9=as.data.frame(Pre_Combined9)
  colnames(Pre_Combined9)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 10
  GEBVsParentPop10_table1=aggregate(as.numeric(GEBVsParentPop10[, 3]), list(GEBVsParentPop10[,1]), mean) #observed a mean for each individual
  GEBVsParentPop10_table2=aggregate(as.numeric(GEBVsParentPop10[, 4]), list(GEBVsParentPop10[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined10=cbind(GEBVsParentPop10_table1[,1],GEBVsParentPop10_table1[,2],GEBVsParentPop10_table2[,2])
  Pre_Combined10=as.data.frame(Pre_Combined10)
  colnames(Pre_Combined10)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 11
  GEBVsParentPop11_table1=aggregate(as.numeric(GEBVsParentPop11[, 3]), list(GEBVsParentPop11[,1]), mean) #observed a mean for each individual
  GEBVsParentPop11_table2=aggregate(as.numeric(GEBVsParentPop11[, 4]), list(GEBVsParentPop11[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined11=cbind(GEBVsParentPop11_table1[,1],GEBVsParentPop11_table1[,2],GEBVsParentPop11_table2[,2])
  Pre_Combined11=as.data.frame(Pre_Combined11)
  colnames(Pre_Combined11)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 12
  GEBVsParentPop12_table1=aggregate(as.numeric(GEBVsParentPop12[, 3]), list(GEBVsParentPop12[,1]), mean) #observed a mean for each individual
  GEBVsParentPop12_table2=aggregate(as.numeric(GEBVsParentPop12[, 4]), list(GEBVsParentPop12[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined12=cbind(GEBVsParentPop12_table1[,1],GEBVsParentPop12_table1[,2],GEBVsParentPop12_table2[,2])
  colnames(Pre_Combined12)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 13
  GEBVsParentPop13_table1=aggregate(as.numeric(GEBVsParentPop13[, 3]), list(GEBVsParentPop13[,1]), mean) #observed a mean for each individual
  GEBVsParentPop13_table2=aggregate(as.numeric(GEBVsParentPop13[, 4]), list(GEBVsParentPop13[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined13=cbind(GEBVsParentPop13_table1[,1],GEBVsParentPop13_table1[,2],GEBVsParentPop13_table2[,2])
  Pre_Combined13=as.data.frame(Pre_Combined13)
  colnames(Pre_Combined13)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 14
  GEBVsParentPop14_table1=aggregate(as.numeric(GEBVsParentPop14[, 3]), list(GEBVsParentPop14[,1]), mean) #observed a mean for each individual
  GEBVsParentPop14_table2=aggregate(as.numeric(GEBVsParentPop14[, 4]), list(GEBVsParentPop14[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined14=cbind(GEBVsParentPop14_table1[,1],GEBVsParentPop14_table1[,2],GEBVsParentPop14_table2[,2])
  colnames(Pre_Combined14)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  Combined=rbind(Pre_Combined1,Pre_Combined2,Pre_Combined3,Pre_Combined4,Pre_Combined5,Pre_Combined6,Pre_Combined7,Pre_Combined8,Pre_Combined9,Pre_Combined10,Pre_Combined11,Pre_Combined12,Pre_Combined13,Pre_Combined14)
  filename1=paste("Predictions_GBLUP_GPD_7Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
  
  # Calculating parental average GEBV
  MeanParentalGEBVs=c(mean(GEBVsParentPop1_table2[,2]),mean(GEBVsParentPop2_table2[,2]),mean(GEBVsParentPop3_table2[,2]),mean(GEBVsParentPop4_table2[,2]),mean(GEBVsParentPop5_table2[,2]),mean(GEBVsParentPop6_table2[,2]),mean(GEBVsParentPop7_table2[,2]),mean(GEBVsParentPop8_table2[,2]),mean(GEBVsParentPop9_table2[,2]),mean(GEBVsParentPop10_table2[,2]),mean(GEBVsParentPop11_table2[,2]),mean(GEBVsParentPop12_table2[,2]),mean(GEBVsParentPop13_table2[,2]),mean(GEBVsParentPop14_table2[,2]))
  F1_prediction_table=cbind(F1observations$Population,MeanParentalGEBVs,F1observations$X20200226_dryweight,F1observations$X20200212_freshweight)
  F1_prediction_table=as.data.frame(F1_prediction_table)
  colnames(F1_prediction_table)=c("F1 population","Mean Parental predicted GEBV","F1 pop mean dryweight","F1 pop mean freshweight")
  filename2=paste("Correlation_Correctedphenotype_GBLUP_GPD_7Replicates",round,".txt",sep="")
  write.table(F1_prediction_table,filename2,sep="\t",quote=F)
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
  GEBVsParentPop1=results6[[1]]
  GEBVsParentPop2=results6[[2]]
  GEBVsParentPop3=results6[[3]]
  GEBVsParentPop4=results6[[4]]
  GEBVsParentPop5=results6[[5]]
  GEBVsParentPop6=results6[[6]]
  GEBVsParentPop7=results6[[7]]
  GEBVsParentPop8=results6[[8]]
  GEBVsParentPop9=results6[[9]]
  GEBVsParentPop10=results6[[10]]
  GEBVsParentPop11=results6[[11]]
  GEBVsParentPop12=results6[[12]]
  GEBVsParentPop13=results6[[13]]
  GEBVsParentPop14=results6[[14]]
  
  # pop 1
  GEBVsParentPop1_table1=aggregate(as.numeric(GEBVsParentPop1[, 3]), list(GEBVsParentPop1[,1]), mean) #observed a mean for each individual
  GEBVsParentPop1_table2=aggregate(as.numeric(GEBVsParentPop1[, 4]), list(GEBVsParentPop1[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined1=cbind(GEBVsParentPop1_table1[,1],GEBVsParentPop1_table1[,2],GEBVsParentPop1_table2[,2])
  Pre_Combined1=as.data.frame(Pre_Combined1)
  colnames(Pre_Combined1)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 2
  GEBVsParentPop2_table1=aggregate(as.numeric(GEBVsParentPop2[, 3]), list(GEBVsParentPop2[,1]), mean) #observed a mean for each individual
  GEBVsParentPop2_table2=aggregate(as.numeric(GEBVsParentPop2[,4]), list(GEBVsParentPop2[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined2=cbind(GEBVsParentPop2_table1[,1],GEBVsParentPop2_table1[,2],GEBVsParentPop2_table2[,2])
  Pre_Combined2=as.data.frame(Pre_Combined2)
  colnames(Pre_Combined2)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 3
  GEBVsParentPop3_table1=aggregate(as.numeric(GEBVsParentPop3[, 3]), list(GEBVsParentPop3[,1]), mean) #observed a mean for each individual
  GEBVsParentPop3_table2=aggregate(as.numeric(GEBVsParentPop3[, 4]), list(GEBVsParentPop3[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined3=cbind(GEBVsParentPop3_table1[,1],GEBVsParentPop3_table1[,2],GEBVsParentPop3_table2[,2])
  Pre_Combined3=as.data.frame(Pre_Combined3)
  colnames(Pre_Combined3)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  #  pop 4
  GEBVsParentPop4_table1=aggregate(as.numeric(GEBVsParentPop4[, 3]), list(GEBVsParentPop4[,1]), mean) #observed a mean for each individual
  GEBVsParentPop4_table2=aggregate(as.numeric(GEBVsParentPop4[, 4]), list(GEBVsParentPop4[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined4=cbind(GEBVsParentPop4_table1[,1],GEBVsParentPop4_table1[,2],GEBVsParentPop4_table2[,2])
  Pre_Combined4=as.data.frame(Pre_Combined4)
  colnames(Pre_Combined4)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 5
  GEBVsParentPop5_table1=aggregate(as.numeric(GEBVsParentPop5[, 3]), list(GEBVsParentPop5[,1]), mean) #observed a mean for each individual
  GEBVsParentPop5_table2=aggregate(as.numeric(GEBVsParentPop5[, 4]), list(GEBVsParentPop5[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined5=cbind(GEBVsParentPop5_table1[,1],GEBVsParentPop5_table1[,2],GEBVsParentPop5_table2[,2])
  Pre_Combined5=as.data.frame(Pre_Combined5)
  colnames(Pre_Combined5)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 6
  GEBVsParentPop6_table1=aggregate(as.numeric(GEBVsParentPop6[, 3]), list(GEBVsParentPop6[,1]), mean) #observed a mean for each individual
  GEBVsParentPop6_table2=aggregate(as.numeric(GEBVsParentPop6[, 4]), list(GEBVsParentPop6[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined6=cbind(GEBVsParentPop6_table1[,1],GEBVsParentPop6_table1[,2],GEBVsParentPop6_table2[,2])
  Pre_Combined6=as.data.frame(Pre_Combined6)
  colnames(Pre_Combined6)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 7
  GEBVsParentPop7_table1=aggregate(as.numeric(GEBVsParentPop7[, 3]), list(GEBVsParentPop7[,1]), mean) #observed a mean for each individual
  GEBVsParentPop7_table2=aggregate(as.numeric(GEBVsParentPop7[, 4]), list(GEBVsParentPop7[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined7=cbind(GEBVsParentPop7_table1[,1],GEBVsParentPop7_table1[,2],GEBVsParentPop7_table2[,2])
  Pre_Combined7=as.data.frame(Pre_Combined7)
  colnames(Pre_Combined7)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 8
  GEBVsParentPop8_table1=aggregate(as.numeric(GEBVsParentPop8[, 3]), list(GEBVsParentPop8[,1]), mean) #observed a mean for each individual
  GEBVsParentPop8_table2=aggregate(as.numeric(GEBVsParentPop8[, 4]), list(GEBVsParentPop8[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined8=cbind(GEBVsParentPop8_table1[,1],GEBVsParentPop8_table1[,2],GEBVsParentPop8_table2[,2])
  Pre_Combined8=as.data.frame(Pre_Combined8)
  colnames(Pre_Combined8)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 9
  GEBVsParentPop9_table1=aggregate(as.numeric(GEBVsParentPop9[, 3]), list(GEBVsParentPop9[,1]), mean) #observed a mean for each individual
  GEBVsParentPop9_table2=aggregate(as.numeric(GEBVsParentPop9[, 4]), list(GEBVsParentPop9[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined9=cbind(GEBVsParentPop9_table1[,1],GEBVsParentPop9_table1[,2],GEBVsParentPop9_table2[,2])
  Pre_Combined9=as.data.frame(Pre_Combined9)
  colnames(Pre_Combined9)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 10
  GEBVsParentPop10_table1=aggregate(as.numeric(GEBVsParentPop10[, 3]), list(GEBVsParentPop10[,1]), mean) #observed a mean for each individual
  GEBVsParentPop10_table2=aggregate(as.numeric(GEBVsParentPop10[, 4]), list(GEBVsParentPop10[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined10=cbind(GEBVsParentPop10_table1[,1],GEBVsParentPop10_table1[,2],GEBVsParentPop10_table2[,2])
  Pre_Combined10=as.data.frame(Pre_Combined10)
  colnames(Pre_Combined10)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 11
  GEBVsParentPop11_table1=aggregate(as.numeric(GEBVsParentPop11[, 3]), list(GEBVsParentPop11[,1]), mean) #observed a mean for each individual
  GEBVsParentPop11_table2=aggregate(as.numeric(GEBVsParentPop11[, 4]), list(GEBVsParentPop11[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined11=cbind(GEBVsParentPop11_table1[,1],GEBVsParentPop11_table1[,2],GEBVsParentPop11_table2[,2])
  Pre_Combined11=as.data.frame(Pre_Combined11)
  colnames(Pre_Combined11)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 12
  GEBVsParentPop12_table1=aggregate(as.numeric(GEBVsParentPop12[, 3]), list(GEBVsParentPop12[,1]), mean) #observed a mean for each individual
  GEBVsParentPop12_table2=aggregate(as.numeric(GEBVsParentPop12[, 4]), list(GEBVsParentPop12[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined12=cbind(GEBVsParentPop12_table1[,1],GEBVsParentPop12_table1[,2],GEBVsParentPop12_table2[,2])
  colnames(Pre_Combined12)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 13
  GEBVsParentPop13_table1=aggregate(as.numeric(GEBVsParentPop13[, 3]), list(GEBVsParentPop13[,1]), mean) #observed a mean for each individual
  GEBVsParentPop13_table2=aggregate(as.numeric(GEBVsParentPop13[, 4]), list(GEBVsParentPop13[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined13=cbind(GEBVsParentPop13_table1[,1],GEBVsParentPop13_table1[,2],GEBVsParentPop13_table2[,2])
  Pre_Combined13=as.data.frame(Pre_Combined13)
  colnames(Pre_Combined13)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 14
  GEBVsParentPop14_table1=aggregate(as.numeric(GEBVsParentPop14[, 3]), list(GEBVsParentPop14[,1]), mean) #observed a mean for each individual
  GEBVsParentPop14_table2=aggregate(as.numeric(GEBVsParentPop14[, 4]), list(GEBVsParentPop14[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined14=cbind(GEBVsParentPop14_table1[,1],GEBVsParentPop14_table1[,2],GEBVsParentPop14_table2[,2])
  colnames(Pre_Combined14)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  Combined=rbind(Pre_Combined1,Pre_Combined2,Pre_Combined3,Pre_Combined4,Pre_Combined5,Pre_Combined6,Pre_Combined7,Pre_Combined8,Pre_Combined9,Pre_Combined10,Pre_Combined11,Pre_Combined12,Pre_Combined13,Pre_Combined14)
  filename1=paste("Predictions_GBLUP_GPD_6Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
  
  # Calculating parental average GEBV
  MeanParentalGEBVs=c(mean(GEBVsParentPop1_table2[,2]),mean(GEBVsParentPop2_table2[,2]),mean(GEBVsParentPop3_table2[,2]),mean(GEBVsParentPop4_table2[,2]),mean(GEBVsParentPop5_table2[,2]),mean(GEBVsParentPop6_table2[,2]),mean(GEBVsParentPop7_table2[,2]),mean(GEBVsParentPop8_table2[,2]),mean(GEBVsParentPop9_table2[,2]),mean(GEBVsParentPop10_table2[,2]),mean(GEBVsParentPop11_table2[,2]),mean(GEBVsParentPop12_table2[,2]),mean(GEBVsParentPop13_table2[,2]),mean(GEBVsParentPop14_table2[,2]))
  F1_prediction_table=cbind(F1observations$Population,MeanParentalGEBVs,F1observations$X20200226_dryweight,F1observations$X20200212_freshweight)
  F1_prediction_table=as.data.frame(F1_prediction_table)
  colnames(F1_prediction_table)=c("F1 population","Mean Parental predicted GEBV","F1 pop mean dryweight","F1 pop mean freshweight")
  filename2=paste("Correlation_Correctedphenotype_GBLUP_GPD_6Replicates",round,".txt",sep="")
  write.table(F1_prediction_table,filename2,sep="\t",quote=F)
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
  GEBVsParentPop1=results5[[1]]
  GEBVsParentPop2=results5[[2]]
  GEBVsParentPop3=results5[[3]]
  GEBVsParentPop4=results5[[4]]
  GEBVsParentPop5=results5[[5]]
  GEBVsParentPop6=results5[[6]]
  GEBVsParentPop7=results5[[7]]
  GEBVsParentPop8=results5[[8]]
  GEBVsParentPop9=results5[[9]]
  GEBVsParentPop10=results5[[10]]
  GEBVsParentPop11=results5[[11]]
  GEBVsParentPop12=results5[[12]]
  GEBVsParentPop13=results5[[13]]
  GEBVsParentPop14=results5[[14]]
  
  # pop 1
  GEBVsParentPop1_table1=aggregate(as.numeric(GEBVsParentPop1[, 3]), list(GEBVsParentPop1[,1]), mean) #observed a mean for each individual
  GEBVsParentPop1_table2=aggregate(as.numeric(GEBVsParentPop1[, 4]), list(GEBVsParentPop1[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined1=cbind(GEBVsParentPop1_table1[,1],GEBVsParentPop1_table1[,2],GEBVsParentPop1_table2[,2])
  Pre_Combined1=as.data.frame(Pre_Combined1)
  colnames(Pre_Combined1)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 2
  GEBVsParentPop2_table1=aggregate(as.numeric(GEBVsParentPop2[, 3]), list(GEBVsParentPop2[,1]), mean) #observed a mean for each individual
  GEBVsParentPop2_table2=aggregate(as.numeric(GEBVsParentPop2[,4]), list(GEBVsParentPop2[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined2=cbind(GEBVsParentPop2_table1[,1],GEBVsParentPop2_table1[,2],GEBVsParentPop2_table2[,2])
  Pre_Combined2=as.data.frame(Pre_Combined2)
  colnames(Pre_Combined2)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 3
  GEBVsParentPop3_table1=aggregate(as.numeric(GEBVsParentPop3[, 3]), list(GEBVsParentPop3[,1]), mean) #observed a mean for each individual
  GEBVsParentPop3_table2=aggregate(as.numeric(GEBVsParentPop3[, 4]), list(GEBVsParentPop3[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined3=cbind(GEBVsParentPop3_table1[,1],GEBVsParentPop3_table1[,2],GEBVsParentPop3_table2[,2])
  Pre_Combined3=as.data.frame(Pre_Combined3)
  colnames(Pre_Combined3)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  #  pop 4
  GEBVsParentPop4_table1=aggregate(as.numeric(GEBVsParentPop4[, 3]), list(GEBVsParentPop4[,1]), mean) #observed a mean for each individual
  GEBVsParentPop4_table2=aggregate(as.numeric(GEBVsParentPop4[, 4]), list(GEBVsParentPop4[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined4=cbind(GEBVsParentPop4_table1[,1],GEBVsParentPop4_table1[,2],GEBVsParentPop4_table2[,2])
  Pre_Combined4=as.data.frame(Pre_Combined4)
  colnames(Pre_Combined4)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 5
  GEBVsParentPop5_table1=aggregate(as.numeric(GEBVsParentPop5[, 3]), list(GEBVsParentPop5[,1]), mean) #observed a mean for each individual
  GEBVsParentPop5_table2=aggregate(as.numeric(GEBVsParentPop5[, 4]), list(GEBVsParentPop5[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined5=cbind(GEBVsParentPop5_table1[,1],GEBVsParentPop5_table1[,2],GEBVsParentPop5_table2[,2])
  Pre_Combined5=as.data.frame(Pre_Combined5)
  colnames(Pre_Combined5)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 6
  GEBVsParentPop6_table1=aggregate(as.numeric(GEBVsParentPop6[, 3]), list(GEBVsParentPop6[,1]), mean) #observed a mean for each individual
  GEBVsParentPop6_table2=aggregate(as.numeric(GEBVsParentPop6[, 4]), list(GEBVsParentPop6[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined6=cbind(GEBVsParentPop6_table1[,1],GEBVsParentPop6_table1[,2],GEBVsParentPop6_table2[,2])
  Pre_Combined6=as.data.frame(Pre_Combined6)
  colnames(Pre_Combined6)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 7
  GEBVsParentPop7_table1=aggregate(as.numeric(GEBVsParentPop7[, 3]), list(GEBVsParentPop7[,1]), mean) #observed a mean for each individual
  GEBVsParentPop7_table2=aggregate(as.numeric(GEBVsParentPop7[, 4]), list(GEBVsParentPop7[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined7=cbind(GEBVsParentPop7_table1[,1],GEBVsParentPop7_table1[,2],GEBVsParentPop7_table2[,2])
  Pre_Combined7=as.data.frame(Pre_Combined7)
  colnames(Pre_Combined7)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 8
  GEBVsParentPop8_table1=aggregate(as.numeric(GEBVsParentPop8[, 3]), list(GEBVsParentPop8[,1]), mean) #observed a mean for each individual
  GEBVsParentPop8_table2=aggregate(as.numeric(GEBVsParentPop8[, 4]), list(GEBVsParentPop8[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined8=cbind(GEBVsParentPop8_table1[,1],GEBVsParentPop8_table1[,2],GEBVsParentPop8_table2[,2])
  Pre_Combined8=as.data.frame(Pre_Combined8)
  colnames(Pre_Combined8)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 9
  GEBVsParentPop9_table1=aggregate(as.numeric(GEBVsParentPop9[, 3]), list(GEBVsParentPop9[,1]), mean) #observed a mean for each individual
  GEBVsParentPop9_table2=aggregate(as.numeric(GEBVsParentPop9[, 4]), list(GEBVsParentPop9[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined9=cbind(GEBVsParentPop9_table1[,1],GEBVsParentPop9_table1[,2],GEBVsParentPop9_table2[,2])
  Pre_Combined9=as.data.frame(Pre_Combined9)
  colnames(Pre_Combined9)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 10
  GEBVsParentPop10_table1=aggregate(as.numeric(GEBVsParentPop10[, 3]), list(GEBVsParentPop10[,1]), mean) #observed a mean for each individual
  GEBVsParentPop10_table2=aggregate(as.numeric(GEBVsParentPop10[, 4]), list(GEBVsParentPop10[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined10=cbind(GEBVsParentPop10_table1[,1],GEBVsParentPop10_table1[,2],GEBVsParentPop10_table2[,2])
  Pre_Combined10=as.data.frame(Pre_Combined10)
  colnames(Pre_Combined10)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 11
  GEBVsParentPop11_table1=aggregate(as.numeric(GEBVsParentPop11[, 3]), list(GEBVsParentPop11[,1]), mean) #observed a mean for each individual
  GEBVsParentPop11_table2=aggregate(as.numeric(GEBVsParentPop11[, 4]), list(GEBVsParentPop11[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined11=cbind(GEBVsParentPop11_table1[,1],GEBVsParentPop11_table1[,2],GEBVsParentPop11_table2[,2])
  Pre_Combined11=as.data.frame(Pre_Combined11)
  colnames(Pre_Combined11)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 12
  GEBVsParentPop12_table1=aggregate(as.numeric(GEBVsParentPop12[, 3]), list(GEBVsParentPop12[,1]), mean) #observed a mean for each individual
  GEBVsParentPop12_table2=aggregate(as.numeric(GEBVsParentPop12[, 4]), list(GEBVsParentPop12[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined12=cbind(GEBVsParentPop12_table1[,1],GEBVsParentPop12_table1[,2],GEBVsParentPop12_table2[,2])
  colnames(Pre_Combined12)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 13
  GEBVsParentPop13_table1=aggregate(as.numeric(GEBVsParentPop13[, 3]), list(GEBVsParentPop13[,1]), mean) #observed a mean for each individual
  GEBVsParentPop13_table2=aggregate(as.numeric(GEBVsParentPop13[, 4]), list(GEBVsParentPop13[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined13=cbind(GEBVsParentPop13_table1[,1],GEBVsParentPop13_table1[,2],GEBVsParentPop13_table2[,2])
  Pre_Combined13=as.data.frame(Pre_Combined13)
  colnames(Pre_Combined13)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 14
  GEBVsParentPop14_table1=aggregate(as.numeric(GEBVsParentPop14[, 3]), list(GEBVsParentPop14[,1]), mean) #observed a mean for each individual
  GEBVsParentPop14_table2=aggregate(as.numeric(GEBVsParentPop14[, 4]), list(GEBVsParentPop14[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined14=cbind(GEBVsParentPop14_table1[,1],GEBVsParentPop14_table1[,2],GEBVsParentPop14_table2[,2])
  colnames(Pre_Combined14)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  Combined=rbind(Pre_Combined1,Pre_Combined2,Pre_Combined3,Pre_Combined4,Pre_Combined5,Pre_Combined6,Pre_Combined7,Pre_Combined8,Pre_Combined9,Pre_Combined10,Pre_Combined11,Pre_Combined12,Pre_Combined13,Pre_Combined14)
  filename1=paste("Predictions_GBLUP_GPD_5Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
  
  # Calculating parental average GEBV
  MeanParentalGEBVs=c(mean(GEBVsParentPop1_table2[,2]),mean(GEBVsParentPop2_table2[,2]),mean(GEBVsParentPop3_table2[,2]),mean(GEBVsParentPop4_table2[,2]),mean(GEBVsParentPop5_table2[,2]),mean(GEBVsParentPop6_table2[,2]),mean(GEBVsParentPop7_table2[,2]),mean(GEBVsParentPop8_table2[,2]),mean(GEBVsParentPop9_table2[,2]),mean(GEBVsParentPop10_table2[,2]),mean(GEBVsParentPop11_table2[,2]),mean(GEBVsParentPop12_table2[,2]),mean(GEBVsParentPop13_table2[,2]),mean(GEBVsParentPop14_table2[,2]))
  F1_prediction_table=cbind(F1observations$Population,MeanParentalGEBVs,F1observations$X20200226_dryweight,F1observations$X20200212_freshweight)
  F1_prediction_table=as.data.frame(F1_prediction_table)
  colnames(F1_prediction_table)=c("F1 population","Mean Parental predicted GEBV","F1 pop mean dryweight","F1 pop mean freshweight")
  filename2=paste("Correlation_Correctedphenotype_GBLUP_GPD_5Replicates",round,".txt",sep="")
  write.table(F1_prediction_table,filename2,sep="\t",quote=F)
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
  GEBVsParentPop1=results4[[1]]
  GEBVsParentPop2=results4[[2]]
  GEBVsParentPop3=results4[[3]]
  GEBVsParentPop4=results4[[4]]
  GEBVsParentPop5=results4[[5]]
  GEBVsParentPop6=results4[[6]]
  GEBVsParentPop7=results4[[7]]
  GEBVsParentPop8=results4[[8]]
  GEBVsParentPop9=results4[[9]]
  GEBVsParentPop10=results4[[10]]
  GEBVsParentPop11=results4[[11]]
  GEBVsParentPop12=results4[[12]]
  GEBVsParentPop13=results4[[13]]
  GEBVsParentPop14=results4[[14]]
  
  # pop 1
  GEBVsParentPop1_table1=aggregate(as.numeric(GEBVsParentPop1[, 3]), list(GEBVsParentPop1[,1]), mean) #observed a mean for each individual
  GEBVsParentPop1_table2=aggregate(as.numeric(GEBVsParentPop1[, 4]), list(GEBVsParentPop1[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined1=cbind(GEBVsParentPop1_table1[,1],GEBVsParentPop1_table1[,2],GEBVsParentPop1_table2[,2])
  Pre_Combined1=as.data.frame(Pre_Combined1)
  colnames(Pre_Combined1)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 2
  GEBVsParentPop2_table1=aggregate(as.numeric(GEBVsParentPop2[, 3]), list(GEBVsParentPop2[,1]), mean) #observed a mean for each individual
  GEBVsParentPop2_table2=aggregate(as.numeric(GEBVsParentPop2[,4]), list(GEBVsParentPop2[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined2=cbind(GEBVsParentPop2_table1[,1],GEBVsParentPop2_table1[,2],GEBVsParentPop2_table2[,2])
  Pre_Combined2=as.data.frame(Pre_Combined2)
  colnames(Pre_Combined2)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 3
  GEBVsParentPop3_table1=aggregate(as.numeric(GEBVsParentPop3[, 3]), list(GEBVsParentPop3[,1]), mean) #observed a mean for each individual
  GEBVsParentPop3_table2=aggregate(as.numeric(GEBVsParentPop3[, 4]), list(GEBVsParentPop3[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined3=cbind(GEBVsParentPop3_table1[,1],GEBVsParentPop3_table1[,2],GEBVsParentPop3_table2[,2])
  Pre_Combined3=as.data.frame(Pre_Combined3)
  colnames(Pre_Combined3)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  #  pop 4
  GEBVsParentPop4_table1=aggregate(as.numeric(GEBVsParentPop4[, 3]), list(GEBVsParentPop4[,1]), mean) #observed a mean for each individual
  GEBVsParentPop4_table2=aggregate(as.numeric(GEBVsParentPop4[, 4]), list(GEBVsParentPop4[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined4=cbind(GEBVsParentPop4_table1[,1],GEBVsParentPop4_table1[,2],GEBVsParentPop4_table2[,2])
  Pre_Combined4=as.data.frame(Pre_Combined4)
  colnames(Pre_Combined4)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 5
  GEBVsParentPop5_table1=aggregate(as.numeric(GEBVsParentPop5[, 3]), list(GEBVsParentPop5[,1]), mean) #observed a mean for each individual
  GEBVsParentPop5_table2=aggregate(as.numeric(GEBVsParentPop5[, 4]), list(GEBVsParentPop5[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined5=cbind(GEBVsParentPop5_table1[,1],GEBVsParentPop5_table1[,2],GEBVsParentPop5_table2[,2])
  Pre_Combined5=as.data.frame(Pre_Combined5)
  colnames(Pre_Combined5)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 6
  GEBVsParentPop6_table1=aggregate(as.numeric(GEBVsParentPop6[, 3]), list(GEBVsParentPop6[,1]), mean) #observed a mean for each individual
  GEBVsParentPop6_table2=aggregate(as.numeric(GEBVsParentPop6[, 4]), list(GEBVsParentPop6[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined6=cbind(GEBVsParentPop6_table1[,1],GEBVsParentPop6_table1[,2],GEBVsParentPop6_table2[,2])
  Pre_Combined6=as.data.frame(Pre_Combined6)
  colnames(Pre_Combined6)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 7
  GEBVsParentPop7_table1=aggregate(as.numeric(GEBVsParentPop7[, 3]), list(GEBVsParentPop7[,1]), mean) #observed a mean for each individual
  GEBVsParentPop7_table2=aggregate(as.numeric(GEBVsParentPop7[, 4]), list(GEBVsParentPop7[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined7=cbind(GEBVsParentPop7_table1[,1],GEBVsParentPop7_table1[,2],GEBVsParentPop7_table2[,2])
  Pre_Combined7=as.data.frame(Pre_Combined7)
  colnames(Pre_Combined7)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 8
  GEBVsParentPop8_table1=aggregate(as.numeric(GEBVsParentPop8[, 3]), list(GEBVsParentPop8[,1]), mean) #observed a mean for each individual
  GEBVsParentPop8_table2=aggregate(as.numeric(GEBVsParentPop8[, 4]), list(GEBVsParentPop8[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined8=cbind(GEBVsParentPop8_table1[,1],GEBVsParentPop8_table1[,2],GEBVsParentPop8_table2[,2])
  Pre_Combined8=as.data.frame(Pre_Combined8)
  colnames(Pre_Combined8)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 9
  GEBVsParentPop9_table1=aggregate(as.numeric(GEBVsParentPop9[, 3]), list(GEBVsParentPop9[,1]), mean) #observed a mean for each individual
  GEBVsParentPop9_table2=aggregate(as.numeric(GEBVsParentPop9[, 4]), list(GEBVsParentPop9[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined9=cbind(GEBVsParentPop9_table1[,1],GEBVsParentPop9_table1[,2],GEBVsParentPop9_table2[,2])
  Pre_Combined9=as.data.frame(Pre_Combined9)
  colnames(Pre_Combined9)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 10
  GEBVsParentPop10_table1=aggregate(as.numeric(GEBVsParentPop10[, 3]), list(GEBVsParentPop10[,1]), mean) #observed a mean for each individual
  GEBVsParentPop10_table2=aggregate(as.numeric(GEBVsParentPop10[, 4]), list(GEBVsParentPop10[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined10=cbind(GEBVsParentPop10_table1[,1],GEBVsParentPop10_table1[,2],GEBVsParentPop10_table2[,2])
  Pre_Combined10=as.data.frame(Pre_Combined10)
  colnames(Pre_Combined10)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 11
  GEBVsParentPop11_table1=aggregate(as.numeric(GEBVsParentPop11[, 3]), list(GEBVsParentPop11[,1]), mean) #observed a mean for each individual
  GEBVsParentPop11_table2=aggregate(as.numeric(GEBVsParentPop11[, 4]), list(GEBVsParentPop11[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined11=cbind(GEBVsParentPop11_table1[,1],GEBVsParentPop11_table1[,2],GEBVsParentPop11_table2[,2])
  Pre_Combined11=as.data.frame(Pre_Combined11)
  colnames(Pre_Combined11)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 12
  GEBVsParentPop12_table1=aggregate(as.numeric(GEBVsParentPop12[, 3]), list(GEBVsParentPop12[,1]), mean) #observed a mean for each individual
  GEBVsParentPop12_table2=aggregate(as.numeric(GEBVsParentPop12[, 4]), list(GEBVsParentPop12[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined12=cbind(GEBVsParentPop12_table1[,1],GEBVsParentPop12_table1[,2],GEBVsParentPop12_table2[,2])
  colnames(Pre_Combined12)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 13
  GEBVsParentPop13_table1=aggregate(as.numeric(GEBVsParentPop13[, 3]), list(GEBVsParentPop13[,1]), mean) #observed a mean for each individual
  GEBVsParentPop13_table2=aggregate(as.numeric(GEBVsParentPop13[, 4]), list(GEBVsParentPop13[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined13=cbind(GEBVsParentPop13_table1[,1],GEBVsParentPop13_table1[,2],GEBVsParentPop13_table2[,2])
  Pre_Combined13=as.data.frame(Pre_Combined13)
  colnames(Pre_Combined13)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 14
  GEBVsParentPop14_table1=aggregate(as.numeric(GEBVsParentPop14[, 3]), list(GEBVsParentPop14[,1]), mean) #observed a mean for each individual
  GEBVsParentPop14_table2=aggregate(as.numeric(GEBVsParentPop14[, 4]), list(GEBVsParentPop14[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined14=cbind(GEBVsParentPop14_table1[,1],GEBVsParentPop14_table1[,2],GEBVsParentPop14_table2[,2])
  colnames(Pre_Combined14)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  Combined=rbind(Pre_Combined1,Pre_Combined2,Pre_Combined3,Pre_Combined4,Pre_Combined5,Pre_Combined6,Pre_Combined7,Pre_Combined8,Pre_Combined9,Pre_Combined10,Pre_Combined11,Pre_Combined12,Pre_Combined13,Pre_Combined14)
  filename1=paste("Predictions_GBLUP_GPD_4Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
  
  # Calculating parental average GEBV
  MeanParentalGEBVs=c(mean(GEBVsParentPop1_table2[,2]),mean(GEBVsParentPop2_table2[,2]),mean(GEBVsParentPop3_table2[,2]),mean(GEBVsParentPop4_table2[,2]),mean(GEBVsParentPop5_table2[,2]),mean(GEBVsParentPop6_table2[,2]),mean(GEBVsParentPop7_table2[,2]),mean(GEBVsParentPop8_table2[,2]),mean(GEBVsParentPop9_table2[,2]),mean(GEBVsParentPop10_table2[,2]),mean(GEBVsParentPop11_table2[,2]),mean(GEBVsParentPop12_table2[,2]),mean(GEBVsParentPop13_table2[,2]),mean(GEBVsParentPop14_table2[,2]))
  F1_prediction_table=cbind(F1observations$Population,MeanParentalGEBVs,F1observations$X20200226_dryweight,F1observations$X20200212_freshweight)
  F1_prediction_table=as.data.frame(F1_prediction_table)
  colnames(F1_prediction_table)=c("F1 population","Mean Parental predicted GEBV","F1 pop mean dryweight","F1 pop mean freshweight")
  filename2=paste("Correlation_Correctedphenotype_GBLUP_GPD_4Replicates",round,".txt",sep="")
  write.table(F1_prediction_table,filename2,sep="\t",quote=F)
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
  GEBVsParentPop1=results3[[1]]
  GEBVsParentPop2=results3[[2]]
  GEBVsParentPop3=results3[[3]]
  GEBVsParentPop4=results3[[4]]
  GEBVsParentPop5=results3[[5]]
  GEBVsParentPop6=results3[[6]]
  GEBVsParentPop7=results3[[7]]
  GEBVsParentPop8=results3[[8]]
  GEBVsParentPop9=results3[[9]]
  GEBVsParentPop10=results3[[10]]
  GEBVsParentPop11=results3[[11]]
  GEBVsParentPop12=results3[[12]]
  GEBVsParentPop13=results3[[13]]
  GEBVsParentPop14=results3[[14]]
  
  # pop 1
  GEBVsParentPop1_table1=aggregate(as.numeric(GEBVsParentPop1[, 3]), list(GEBVsParentPop1[,1]), mean) #observed a mean for each individual
  GEBVsParentPop1_table2=aggregate(as.numeric(GEBVsParentPop1[, 4]), list(GEBVsParentPop1[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined1=cbind(GEBVsParentPop1_table1[,1],GEBVsParentPop1_table1[,2],GEBVsParentPop1_table2[,2])
  Pre_Combined1=as.data.frame(Pre_Combined1)
  colnames(Pre_Combined1)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 2
  GEBVsParentPop2_table1=aggregate(as.numeric(GEBVsParentPop2[, 3]), list(GEBVsParentPop2[,1]), mean) #observed a mean for each individual
  GEBVsParentPop2_table2=aggregate(as.numeric(GEBVsParentPop2[,4]), list(GEBVsParentPop2[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined2=cbind(GEBVsParentPop2_table1[,1],GEBVsParentPop2_table1[,2],GEBVsParentPop2_table2[,2])
  Pre_Combined2=as.data.frame(Pre_Combined2)
  colnames(Pre_Combined2)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 3
  GEBVsParentPop3_table1=aggregate(as.numeric(GEBVsParentPop3[, 3]), list(GEBVsParentPop3[,1]), mean) #observed a mean for each individual
  GEBVsParentPop3_table2=aggregate(as.numeric(GEBVsParentPop3[, 4]), list(GEBVsParentPop3[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined3=cbind(GEBVsParentPop3_table1[,1],GEBVsParentPop3_table1[,2],GEBVsParentPop3_table2[,2])
  Pre_Combined3=as.data.frame(Pre_Combined3)
  colnames(Pre_Combined3)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  #  pop 4
  GEBVsParentPop4_table1=aggregate(as.numeric(GEBVsParentPop4[, 3]), list(GEBVsParentPop4[,1]), mean) #observed a mean for each individual
  GEBVsParentPop4_table2=aggregate(as.numeric(GEBVsParentPop4[, 4]), list(GEBVsParentPop4[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined4=cbind(GEBVsParentPop4_table1[,1],GEBVsParentPop4_table1[,2],GEBVsParentPop4_table2[,2])
  Pre_Combined4=as.data.frame(Pre_Combined4)
  colnames(Pre_Combined4)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 5
  GEBVsParentPop5_table1=aggregate(as.numeric(GEBVsParentPop5[, 3]), list(GEBVsParentPop5[,1]), mean) #observed a mean for each individual
  GEBVsParentPop5_table2=aggregate(as.numeric(GEBVsParentPop5[, 4]), list(GEBVsParentPop5[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined5=cbind(GEBVsParentPop5_table1[,1],GEBVsParentPop5_table1[,2],GEBVsParentPop5_table2[,2])
  Pre_Combined5=as.data.frame(Pre_Combined5)
  colnames(Pre_Combined5)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 6
  GEBVsParentPop6_table1=aggregate(as.numeric(GEBVsParentPop6[, 3]), list(GEBVsParentPop6[,1]), mean) #observed a mean for each individual
  GEBVsParentPop6_table2=aggregate(as.numeric(GEBVsParentPop6[, 4]), list(GEBVsParentPop6[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined6=cbind(GEBVsParentPop6_table1[,1],GEBVsParentPop6_table1[,2],GEBVsParentPop6_table2[,2])
  Pre_Combined6=as.data.frame(Pre_Combined6)
  colnames(Pre_Combined6)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 7
  GEBVsParentPop7_table1=aggregate(as.numeric(GEBVsParentPop7[, 3]), list(GEBVsParentPop7[,1]), mean) #observed a mean for each individual
  GEBVsParentPop7_table2=aggregate(as.numeric(GEBVsParentPop7[, 4]), list(GEBVsParentPop7[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined7=cbind(GEBVsParentPop7_table1[,1],GEBVsParentPop7_table1[,2],GEBVsParentPop7_table2[,2])
  Pre_Combined7=as.data.frame(Pre_Combined7)
  colnames(Pre_Combined7)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 8
  GEBVsParentPop8_table1=aggregate(as.numeric(GEBVsParentPop8[, 3]), list(GEBVsParentPop8[,1]), mean) #observed a mean for each individual
  GEBVsParentPop8_table2=aggregate(as.numeric(GEBVsParentPop8[, 4]), list(GEBVsParentPop8[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined8=cbind(GEBVsParentPop8_table1[,1],GEBVsParentPop8_table1[,2],GEBVsParentPop8_table2[,2])
  Pre_Combined8=as.data.frame(Pre_Combined8)
  colnames(Pre_Combined8)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 9
  GEBVsParentPop9_table1=aggregate(as.numeric(GEBVsParentPop9[, 3]), list(GEBVsParentPop9[,1]), mean) #observed a mean for each individual
  GEBVsParentPop9_table2=aggregate(as.numeric(GEBVsParentPop9[, 4]), list(GEBVsParentPop9[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined9=cbind(GEBVsParentPop9_table1[,1],GEBVsParentPop9_table1[,2],GEBVsParentPop9_table2[,2])
  Pre_Combined9=as.data.frame(Pre_Combined9)
  colnames(Pre_Combined9)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 10
  GEBVsParentPop10_table1=aggregate(as.numeric(GEBVsParentPop10[, 3]), list(GEBVsParentPop10[,1]), mean) #observed a mean for each individual
  GEBVsParentPop10_table2=aggregate(as.numeric(GEBVsParentPop10[, 4]), list(GEBVsParentPop10[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined10=cbind(GEBVsParentPop10_table1[,1],GEBVsParentPop10_table1[,2],GEBVsParentPop10_table2[,2])
  Pre_Combined10=as.data.frame(Pre_Combined10)
  colnames(Pre_Combined10)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 11
  GEBVsParentPop11_table1=aggregate(as.numeric(GEBVsParentPop11[, 3]), list(GEBVsParentPop11[,1]), mean) #observed a mean for each individual
  GEBVsParentPop11_table2=aggregate(as.numeric(GEBVsParentPop11[, 4]), list(GEBVsParentPop11[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined11=cbind(GEBVsParentPop11_table1[,1],GEBVsParentPop11_table1[,2],GEBVsParentPop11_table2[,2])
  Pre_Combined11=as.data.frame(Pre_Combined11)
  colnames(Pre_Combined11)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 12
  GEBVsParentPop12_table1=aggregate(as.numeric(GEBVsParentPop12[, 3]), list(GEBVsParentPop12[,1]), mean) #observed a mean for each individual
  GEBVsParentPop12_table2=aggregate(as.numeric(GEBVsParentPop12[, 4]), list(GEBVsParentPop12[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined12=cbind(GEBVsParentPop12_table1[,1],GEBVsParentPop12_table1[,2],GEBVsParentPop12_table2[,2])
  colnames(Pre_Combined12)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 13
  GEBVsParentPop13_table1=aggregate(as.numeric(GEBVsParentPop13[, 3]), list(GEBVsParentPop13[,1]), mean) #observed a mean for each individual
  GEBVsParentPop13_table2=aggregate(as.numeric(GEBVsParentPop13[, 4]), list(GEBVsParentPop13[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined13=cbind(GEBVsParentPop13_table1[,1],GEBVsParentPop13_table1[,2],GEBVsParentPop13_table2[,2])
  Pre_Combined13=as.data.frame(Pre_Combined13)
  colnames(Pre_Combined13)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 14
  GEBVsParentPop14_table1=aggregate(as.numeric(GEBVsParentPop14[, 3]), list(GEBVsParentPop14[,1]), mean) #observed a mean for each individual
  GEBVsParentPop14_table2=aggregate(as.numeric(GEBVsParentPop14[, 4]), list(GEBVsParentPop14[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined14=cbind(GEBVsParentPop14_table1[,1],GEBVsParentPop14_table1[,2],GEBVsParentPop14_table2[,2])
  colnames(Pre_Combined14)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  Combined=rbind(Pre_Combined1,Pre_Combined2,Pre_Combined3,Pre_Combined4,Pre_Combined5,Pre_Combined6,Pre_Combined7,Pre_Combined8,Pre_Combined9,Pre_Combined10,Pre_Combined11,Pre_Combined12,Pre_Combined13,Pre_Combined14)
  filename1=paste("Predictions_GBLUP_GPD_3Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
  
  # Calculating parental average GEBV
  MeanParentalGEBVs=c(mean(GEBVsParentPop1_table2[,2]),mean(GEBVsParentPop2_table2[,2]),mean(GEBVsParentPop3_table2[,2]),mean(GEBVsParentPop4_table2[,2]),mean(GEBVsParentPop5_table2[,2]),mean(GEBVsParentPop6_table2[,2]),mean(GEBVsParentPop7_table2[,2]),mean(GEBVsParentPop8_table2[,2]),mean(GEBVsParentPop9_table2[,2]),mean(GEBVsParentPop10_table2[,2]),mean(GEBVsParentPop11_table2[,2]),mean(GEBVsParentPop12_table2[,2]),mean(GEBVsParentPop13_table2[,2]),mean(GEBVsParentPop14_table2[,2]))
  F1_prediction_table=cbind(F1observations$Population,MeanParentalGEBVs,F1observations$X20200226_dryweight,F1observations$X20200212_freshweight)
  F1_prediction_table=as.data.frame(F1_prediction_table)
  colnames(F1_prediction_table)=c("F1 population","Mean Parental predicted GEBV","F1 pop mean dryweight","F1 pop mean freshweight")
  filename2=paste("Correlation_Correctedphenotype_GBLUP_GPD_3Replicates",round,".txt",sep="")
  write.table(F1_prediction_table,filename2,sep="\t",quote=F)
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
  GEBVsParentPop1=results2[[1]]
  GEBVsParentPop2=results2[[2]]
  GEBVsParentPop3=results2[[3]]
  GEBVsParentPop4=results2[[4]]
  GEBVsParentPop5=results2[[5]]
  GEBVsParentPop6=results2[[6]]
  GEBVsParentPop7=results2[[7]]
  GEBVsParentPop8=results2[[8]]
  GEBVsParentPop9=results2[[9]]
  GEBVsParentPop10=results2[[10]]
  GEBVsParentPop11=results2[[11]]
  GEBVsParentPop12=results2[[12]]
  GEBVsParentPop13=results2[[13]]
  GEBVsParentPop14=results2[[14]]
  
  # pop 1
  GEBVsParentPop1_table1=aggregate(as.numeric(GEBVsParentPop1[, 3]), list(GEBVsParentPop1[,1]), mean) #observed a mean for each individual
  GEBVsParentPop1_table2=aggregate(as.numeric(GEBVsParentPop1[, 4]), list(GEBVsParentPop1[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined1=cbind(GEBVsParentPop1_table1[,1],GEBVsParentPop1_table1[,2],GEBVsParentPop1_table2[,2])
  Pre_Combined1=as.data.frame(Pre_Combined1)
  colnames(Pre_Combined1)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 2
  GEBVsParentPop2_table1=aggregate(as.numeric(GEBVsParentPop2[, 3]), list(GEBVsParentPop2[,1]), mean) #observed a mean for each individual
  GEBVsParentPop2_table2=aggregate(as.numeric(GEBVsParentPop2[,4]), list(GEBVsParentPop2[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined2=cbind(GEBVsParentPop2_table1[,1],GEBVsParentPop2_table1[,2],GEBVsParentPop2_table2[,2])
  Pre_Combined2=as.data.frame(Pre_Combined2)
  colnames(Pre_Combined2)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 3
  GEBVsParentPop3_table1=aggregate(as.numeric(GEBVsParentPop3[, 3]), list(GEBVsParentPop3[,1]), mean) #observed a mean for each individual
  GEBVsParentPop3_table2=aggregate(as.numeric(GEBVsParentPop3[, 4]), list(GEBVsParentPop3[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined3=cbind(GEBVsParentPop3_table1[,1],GEBVsParentPop3_table1[,2],GEBVsParentPop3_table2[,2])
  Pre_Combined3=as.data.frame(Pre_Combined3)
  colnames(Pre_Combined3)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  #  pop 4
  GEBVsParentPop4_table1=aggregate(as.numeric(GEBVsParentPop4[, 3]), list(GEBVsParentPop4[,1]), mean) #observed a mean for each individual
  GEBVsParentPop4_table2=aggregate(as.numeric(GEBVsParentPop4[, 4]), list(GEBVsParentPop4[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined4=cbind(GEBVsParentPop4_table1[,1],GEBVsParentPop4_table1[,2],GEBVsParentPop4_table2[,2])
  Pre_Combined4=as.data.frame(Pre_Combined4)
  colnames(Pre_Combined4)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 5
  GEBVsParentPop5_table1=aggregate(as.numeric(GEBVsParentPop5[, 3]), list(GEBVsParentPop5[,1]), mean) #observed a mean for each individual
  GEBVsParentPop5_table2=aggregate(as.numeric(GEBVsParentPop5[, 4]), list(GEBVsParentPop5[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined5=cbind(GEBVsParentPop5_table1[,1],GEBVsParentPop5_table1[,2],GEBVsParentPop5_table2[,2])
  Pre_Combined5=as.data.frame(Pre_Combined5)
  colnames(Pre_Combined5)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 6
  GEBVsParentPop6_table1=aggregate(as.numeric(GEBVsParentPop6[, 3]), list(GEBVsParentPop6[,1]), mean) #observed a mean for each individual
  GEBVsParentPop6_table2=aggregate(as.numeric(GEBVsParentPop6[, 4]), list(GEBVsParentPop6[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined6=cbind(GEBVsParentPop6_table1[,1],GEBVsParentPop6_table1[,2],GEBVsParentPop6_table2[,2])
  Pre_Combined6=as.data.frame(Pre_Combined6)
  colnames(Pre_Combined6)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 7
  GEBVsParentPop7_table1=aggregate(as.numeric(GEBVsParentPop7[, 3]), list(GEBVsParentPop7[,1]), mean) #observed a mean for each individual
  GEBVsParentPop7_table2=aggregate(as.numeric(GEBVsParentPop7[, 4]), list(GEBVsParentPop7[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined7=cbind(GEBVsParentPop7_table1[,1],GEBVsParentPop7_table1[,2],GEBVsParentPop7_table2[,2])
  Pre_Combined7=as.data.frame(Pre_Combined7)
  colnames(Pre_Combined7)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 8
  GEBVsParentPop8_table1=aggregate(as.numeric(GEBVsParentPop8[, 3]), list(GEBVsParentPop8[,1]), mean) #observed a mean for each individual
  GEBVsParentPop8_table2=aggregate(as.numeric(GEBVsParentPop8[, 4]), list(GEBVsParentPop8[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined8=cbind(GEBVsParentPop8_table1[,1],GEBVsParentPop8_table1[,2],GEBVsParentPop8_table2[,2])
  Pre_Combined8=as.data.frame(Pre_Combined8)
  colnames(Pre_Combined8)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 9
  GEBVsParentPop9_table1=aggregate(as.numeric(GEBVsParentPop9[, 3]), list(GEBVsParentPop9[,1]), mean) #observed a mean for each individual
  GEBVsParentPop9_table2=aggregate(as.numeric(GEBVsParentPop9[, 4]), list(GEBVsParentPop9[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined9=cbind(GEBVsParentPop9_table1[,1],GEBVsParentPop9_table1[,2],GEBVsParentPop9_table2[,2])
  Pre_Combined9=as.data.frame(Pre_Combined9)
  colnames(Pre_Combined9)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 10
  GEBVsParentPop10_table1=aggregate(as.numeric(GEBVsParentPop10[, 3]), list(GEBVsParentPop10[,1]), mean) #observed a mean for each individual
  GEBVsParentPop10_table2=aggregate(as.numeric(GEBVsParentPop10[, 4]), list(GEBVsParentPop10[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined10=cbind(GEBVsParentPop10_table1[,1],GEBVsParentPop10_table1[,2],GEBVsParentPop10_table2[,2])
  Pre_Combined10=as.data.frame(Pre_Combined10)
  colnames(Pre_Combined10)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 11
  GEBVsParentPop11_table1=aggregate(as.numeric(GEBVsParentPop11[, 3]), list(GEBVsParentPop11[,1]), mean) #observed a mean for each individual
  GEBVsParentPop11_table2=aggregate(as.numeric(GEBVsParentPop11[, 4]), list(GEBVsParentPop11[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined11=cbind(GEBVsParentPop11_table1[,1],GEBVsParentPop11_table1[,2],GEBVsParentPop11_table2[,2])
  Pre_Combined11=as.data.frame(Pre_Combined11)
  colnames(Pre_Combined11)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 12
  GEBVsParentPop12_table1=aggregate(as.numeric(GEBVsParentPop12[, 3]), list(GEBVsParentPop12[,1]), mean) #observed a mean for each individual
  GEBVsParentPop12_table2=aggregate(as.numeric(GEBVsParentPop12[, 4]), list(GEBVsParentPop12[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined12=cbind(GEBVsParentPop12_table1[,1],GEBVsParentPop12_table1[,2],GEBVsParentPop12_table2[,2])
  colnames(Pre_Combined12)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 13
  GEBVsParentPop13_table1=aggregate(as.numeric(GEBVsParentPop13[, 3]), list(GEBVsParentPop13[,1]), mean) #observed a mean for each individual
  GEBVsParentPop13_table2=aggregate(as.numeric(GEBVsParentPop13[, 4]), list(GEBVsParentPop13[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined13=cbind(GEBVsParentPop13_table1[,1],GEBVsParentPop13_table1[,2],GEBVsParentPop13_table2[,2])
  Pre_Combined13=as.data.frame(Pre_Combined13)
  colnames(Pre_Combined13)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  # pop 14
  GEBVsParentPop14_table1=aggregate(as.numeric(GEBVsParentPop14[, 3]), list(GEBVsParentPop14[,1]), mean) #observed a mean for each individual
  GEBVsParentPop14_table2=aggregate(as.numeric(GEBVsParentPop14[, 4]), list(GEBVsParentPop14[,1]), mean) #total GEBV, a mean for each individual
  
  Pre_Combined14=cbind(GEBVsParentPop14_table1[,1],GEBVsParentPop14_table1[,2],GEBVsParentPop14_table2[,2])
  colnames(Pre_Combined14)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
  
  Combined=rbind(Pre_Combined1,Pre_Combined2,Pre_Combined3,Pre_Combined4,Pre_Combined5,Pre_Combined6,Pre_Combined7,Pre_Combined8,Pre_Combined9,Pre_Combined10,Pre_Combined11,Pre_Combined12,Pre_Combined13,Pre_Combined14)
  filename1=paste("Predictions_GBLUP_GPD_2Replicates",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t",quote=F)
  
  # Calculating parental average GEBV
  MeanParentalGEBVs=c(mean(GEBVsParentPop1_table2[,2]),mean(GEBVsParentPop2_table2[,2]),mean(GEBVsParentPop3_table2[,2]),mean(GEBVsParentPop4_table2[,2]),mean(GEBVsParentPop5_table2[,2]),mean(GEBVsParentPop6_table2[,2]),mean(GEBVsParentPop7_table2[,2]),mean(GEBVsParentPop8_table2[,2]),mean(GEBVsParentPop9_table2[,2]),mean(GEBVsParentPop10_table2[,2]),mean(GEBVsParentPop11_table2[,2]),mean(GEBVsParentPop12_table2[,2]),mean(GEBVsParentPop13_table2[,2]),mean(GEBVsParentPop14_table2[,2]))
  F1_prediction_table=cbind(F1observations$Population,MeanParentalGEBVs,F1observations$X20200226_dryweight,F1observations$X20200212_freshweight)
  F1_prediction_table=as.data.frame(F1_prediction_table)
  colnames(F1_prediction_table)=c("F1 population","Mean Parental predicted GEBV","F1 pop mean dryweight","F1 pop mean freshweight")
  filename2=paste("Correlation_Correctedphenotype_GBLUP_GPD_2Replicates",round,".txt",sep="")
  write.table(F1_prediction_table,filename2,sep="\t",quote=F)
}  
