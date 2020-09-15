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
  
  d6_training=d6[-testpop,] # limit the dataframe to only the individuals allowed for training the model
  d6_training_ready=na.omit(d6_training, cols = c("growth_per_day")) # remember that gpd na inidividuals should be removed whether or not they are in the training pop or not
  
  ind_not_in_train=d6$Clovershort[testpop]
  IndividualsToRemoveGRM=which(colnames(GRM1) %in% ind_not_in_train)
  
  GRM_trn = GRM1[-IndividualsToRemoveGRM,-IndividualsToRemoveGRM]
  
  GcloverReps_trn = GcloverReps[-testpop,-testpop] # a GRM for the individuals used to predict
  CloverIndep_trn= CloverIndep[-testpop,-testpop] # another GRM for the individuals used to predict
  
  # Run the GBLUP model on full training population to extract GEBVs
  yNA=d6_training_ready$growth_per_day
  
  fixedmod=model.matrix(~factor(d6_training_ready$Round)+factor(d6_training_ready$NS)+factor(d6_training_ready$EW)+factor(d6_training_ready$Rhizobium)+factor(d6_training_ready$inoculation_date))
  ETA=list(list(K=GcloverReps_trn,model="RKHS"),list(K=CloverIndep_trn,model="RKHS"),list(X=fixedmod,model="FIXED"))
  GBLUP=BGLR(y=yNA,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("GBLUP",round))
  matrix=cbind(as.character(d6_training_ready$Clovershort),as.numeric(d6_training_ready$CorrectedPheno),as.numeric(GBLUP$ETA[[1]]$u),as.numeric(GBLUP$ETA[[2]]$u),(as.numeric(GBLUP$ETA[[1]]$u)+as.numeric(GBLUP$ETA[[2]]$u)))
  colnames(matrix)=c("ID", "Observed", "1st GEBV contribution","2nd GEBV contribution (clover independent)","totalGEBV")
  GEBV_contribution1data=as.numeric(as.character(matrix[,3]))
  GEBV_contribution1data=as.matrix(GEBV_contribution1data)
  GEBV_contribution2data=as.numeric(as.character(matrix[,4]))
  
  ################  ################ 
  ## Now predict testing population 
  ################  ################  
  
  GcloverReps_test = GcloverReps[testpop,testpop] # GRM for individuals that will be predicted
  GcloverReps_covar = GcloverReps[testpop,-testpop] # Covariance between training and testing pop.
  
  #GEBVpred_contr1 = GcloverReps_covar%*%solve(GcloverReps_trn) %*% GEBV_contribution1data 
  GEBVpred = GcloverReps_covar%*%ginv(GcloverReps_trn) %*% GEBV_contribution1data 
  #GEBVpred_contr1 = GcloverReps_covar%*%solve(GcloverReps_trn + diag(0.01, 1661, 1661)) %*% GEBV_contribution1data 
  
  # Output matrix with prediction results
  matrix1=cbind(as.character(d6$Clovershort[testpop]),as.numeric(d6$CorrectedPheno[testpop]),as.numeric(as.character(GEBVpred)))
  colnames(matrix1)=c("ID", "Observed", "GEBV")
  return(matrix1)
}
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

# Remove round 1 replicate 2
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




# Calculate a phenotype corrected for all fixed effects
{
  Correctedforallfixed <- lmer(growth_per_day ~ factor(Round) + factor(NS) + factor(EW) + factor(Rhizobium) + inoculation_date + (1|Clover), data=d6) 
  summary(Correctedforallfixed)
  matrixOfeffects=model.matrix( ~ factor(Round) + factor(NS) + factor(EW)  + factor(Rhizobium) + inoculation_date, data=d6)
  column_to_remove=which(colnames(matrixOfeffects) %in% names(fixef(Correctedforallfixed))==F)
  matrixOfeffects_new=matrixOfeffects[,-column_to_remove]

  ycorr <- d6$growth_per_day - matrixOfeffects_new %*% fixef(Correctedforallfixed)
  d6$CorrectedPheno <- ycorr 
}  

# Make the matrices that goes into the model
{
  CloverDesign <- model.matrix(~0+d6$Clovershort)
  GcloverReps <- CloverDesign %*% GRM1 %*% t(CloverDesign) 
  #GcloverReps is a G-matrix that will match the size and layout in the data and can be used in BGLR at the K=
  #Add a clover effect to capture non-additve variance that makes broad sense heritability.
  #CloverIndep catch clover-effects without relationships (clovers are independent), that is another matrix that can go into the model. 
  CloverIndep <- CloverDesign %*% t(CloverDesign)
  dim(CloverIndep)
}

# Setting up 6-fold CV system
{
# if groups are not provided 
# Make six groups for 6-fold CV

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
  
  testpop1_idx=which(d6$Clovershort %in% testpop1)
  testpop2_idx=which(d6$Clovershort %in% testpop2)
  testpop3_idx=which(d6$Clovershort %in% testpop3)
  testpop4_idx=which(d6$Clovershort %in% testpop4)
  testpop5_idx=which(d6$Clovershort %in% testpop5)
  testpop6_idx=which(d6$Clovershort %in% testpop6)
  
  tests=list(testpop1_idx,testpop2_idx,testpop3_idx,testpop4_idx,testpop5_idx,testpop6_idx)
}
}

#  Calculate heritability
{
  print("Calculating heritability")
  y=d6$growth_per_day
  fixedmod=model.matrix(~factor(d6$Round)+factor(d6$NS)+factor(d6$EW)+factor(d6$Rhizobium)+factor(d6$inoculation_date))
  ETA=list(list(K=GcloverReps,model="RKHS"),list(K=CloverIndep,model="RKHS"),list(X=fixedmod,model="FIXED")) 
  GBLUP=BGLR(y=y,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("heritabilityestimation",round))
  
  print(GBLUP$ETA[[1]]$varU)
  print(GBLUP$ETA[[2]]$varU)
  print(GBLUP$varE)
  h2=(GBLUP$ETA[[1]]$varU)/(GBLUP$varE+GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)
  h2_=(GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)/(GBLUP$varE+GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)
  write.table(h2_,paste("h2_includingNonAdditive",round,".txt",sep=""),sep="\t")
  write.table(h2,paste("h2",round,".txt",sep=""),sep="\t")
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
  
  All1=aggregate(as.numeric(All[, 2]), list(All[,1]), mean)
  All2=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  correlation=cor(All1[,2],All2[,2]) #means of replicates
  correlation
  
  Combined=cbind(All1[,1],All1[,2],All2[,2])
  colnames(Combined)=c("Individual","Observed_correctedForFixed","GEBVs")
  
  filename=paste("Correlation_GBLUP_GPD",round,".txt",sep="")
  write.table(correlation,filename,sep="\t")
  
  filename1=paste("Predictions_GBLUP_GPD",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t")
}  
