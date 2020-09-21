####################################################################################
# Genomic prediction fitting fixed effects on yield data with replicates.          #
####################################################################################
#                       GBLUP                                                      #
####################################################################################

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

# read in data
{
  d <- read.csv("/home/cks/NChain/faststorage/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/greenhouse_area.csv", header = TRUE, sep = ",")
  f=read.csv("/home/cks/NChain/faststorage/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/2018_weight.csv",header=T,sep=";")
  colnames(f)[1]="Barcode"
  df=merge(d,f,by="Barcode")
  d=df
  
}

#Calculate growth pr. day
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

# Fix data matrix so that names are the same as in GRM, that is 118 individuals in each
{
  GRM=read.table(args[1],sep=",",header=T)
  #GRM=read.table("/home/cks/NChain/faststorage/WHITE_CLOVER/GPJune2020_CKS/Preparation_of_genotypeFile/GRM_Clover_Fullfiltering_20200623.csv",sep=",",header=T)
  #GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/Nchain/Genomic_Prediction_of_Yield_redoneJune2020/GRM_Clover_Fullfiltering_20200623_noutliers.csv",sep=",",header=T)
  #GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/Nchain/Genomic_Prediction_of_Yield_redoneJune2020/GRM_Clover_Fullfiltering_20200623.csv",sep=",",header=T)
  #GRM=read.table("GRM_Clover_Fullfiltering_145NoRbani8ind_20200701.csv",sep=",",header=T)
  #GRM=read.table("GRM_Clover_Fullfiltering_145NoRbani8indNoAearl_07.csv",sep=",",header=T)
  dim(GRM)
  d2$Clovershort <- strtrim(d2$Clover,8)
  d3=d2[order(d2$Clovershort,decreasing=F),]
  length(unique(d3$Clovershort)) #149
  d4=d3[which(d3$Clovershort %in% colnames(GRM)),]
  length(unique(d4$Clovershort)) #147 unique genotypes with GPD data
  remove=GRM[which(colnames(GRM) %in% d4$Clovershort),which(colnames(GRM) %in% d4$Clovershort)]
  print(remove)
  IndWithGenoButNoPheno=colnames(GRM)[which(colnames(GRM) %in% d4$Clovershort==F)]
  which(IndWithGenoButNoPheno %in% c("Aoost_03","Aoost_06","Aaran_06","Aaran_07","Kdike_10"))
  toremove=which(colnames(GRM) %in% d4$Clovershort==F)[-which(IndWithGenoButNoPheno %in% c("Aoost_03","Aoost_06","Aaran_06","Aaran_07","Kdike_10"))]
  namestoremove=colnames(GRM)[toremove]
  idxtoremove=which(colnames(GRM) %in% namestoremove)
  
  GRM1=GRM[-idxtoremove,-idxtoremove]
  dim(GRM1)
  GRM1=data.matrix(GRM1)
}


#Aberpearl_07 contribute 700 of the datapoints and thus influence the variance a lot. Cut down Aberpearl_07 data so that we have only 6 different Rhizobia left like the other clovers
Aberpearl_07=which(d4$Clovershort=="Aearl_07")
Inocolums=unique(d4$Rhizobium[Aberpearl_07])
set.seed(15)
sample=sample(Inocolums,6)
#sample=c("MIX","SM22","SM25","SM149C","SM31","SM155A")
print(sample)
set.seed(NULL)


which(d4$Rhizobium[Aberpearl_07]==sample[1]) #4
which(d4$Rhizobium[Aberpearl_07]==sample[2]) #4
which(d4$Rhizobium[Aberpearl_07]==sample[3]) #4
which(d4$Rhizobium[Aberpearl_07]==sample[4]) #4
which(d4$Rhizobium[Aberpearl_07]==sample[5]) #4
which(d4$Rhizobium[Aberpearl_07]==sample[6]) #4

remove=which((d4$Rhizobium[Aberpearl_07] %in% sample)==FALSE)

d4=d4[-Aberpearl_07[remove],]
nrow(d4)

d6=d4
length(which(d6$Clovershort=="Aearl_07"))
d6$Rhizobium=droplevels(d6$Rhizobium) # removing levels not used in actual data
d6$Clover=droplevels(d6$Clover) # removing levels not used in actual data
nrow(d6)

# Correct GPD for the part of initial size that does not carry genetic information
{
  lm.fit <- lm(d6$InitialSize ~ d6$Clover)
  summary(lm.fit)
  d6$residuals <- round(lm.fit$residuals,2)
  
  fit <- lmer(growth_per_day ~ residuals + (1|Clover), data=d6) 
  ycorr <- d6$growth_per_day - model.matrix( ~ residuals, data=d6) %*% fixef(fit)
  d6$gpd_dryweight_cor <- ycorr #this is the new corrected dry weight
  
  print(cor(d6$gpd_dryweight_cor,d6$InitialSize))
  print(cor(d6$growth_per_day,d6$InitialSize))
}


# Calculate a phenotype corrected for all fixed effects
{
  Correctedforallfixed <- lmer(gpd_dryweight_cor ~ factor(Round) + factor(NS) + factor(EW) + factor(Rhizobium) + factor(inoculation_date) + (1|Clover), data=d6) 
  summary(Correctedforallfixed)
  matrixOfeffects=model.matrix( ~ factor(Round) + factor(NS) + factor(EW) + factor(Rhizobium) + factor(inoculation_date) , data=d6)
  column_to_remove=which(colnames(matrixOfeffects) %in% names(fixef(Correctedforallfixed))==F)
  matrixOfeffects_new=matrixOfeffects[,-column_to_remove]
  
  ycorr <- d6$gpd_dryweight_cor - matrixOfeffects_new %*% fixef(Correctedforallfixed)
  d6$CorrectedPheno <- ycorr 
}  


# make 4 fake rows in d6
fakedata1=c(rep("NA",7),"Aoost_03",rep("NA",18),"Aoost_03",rep("NA",4))
fakedata2=c(rep("NA",7),"Aoost_06",rep("NA",18),"Aoost_06",rep("NA",4))
fakedata3=c(rep("NA",7),"Aaran_06",rep("NA",18),"Aaran_06",rep("NA",4))
fakedata4=c(rep("NA",7),"Kdike_10",rep("NA",18),"Kdike_10",rep("NA",4))
d6=rbind(d6,fakedata1,fakedata2,fakedata3,fakedata4)
d6=d6[order(d6$Clovershort),] 
d6$gpd_dryweight_cor=as.numeric(as.character(d6$gpd_dryweight_cor))

d6$Round=as.numeric(as.character(d6$Round))
d6$NS=as.numeric(as.character(d6$NS))
d6$EW=as.numeric(as.character(d6$EW))
d6$Rhizobium=as.character(d6$Rhizobium)

# make a column that specifies variety
d6$Variety <- strtrim(d6$Clovershort,5)

# Divide into F1 parental groups
{
  Variety1_idx=which(d6$Variety==unique(d6$Variety)[1])
  Variety2_idx=which(d6$Variety==unique(d6$Variety)[2])
  Variety3_idx=which(d6$Variety==unique(d6$Variety)[3])
  Variety4_idx=which(d6$Variety==unique(d6$Variety)[4])
  Variety5_idx=which(d6$Variety==unique(d6$Variety)[5])
  Variety6_idx=which(d6$Variety==unique(d6$Variety)[6])
  Variety7_idx=which(d6$Variety==unique(d6$Variety)[7])
  Variety8_idx=which(d6$Variety==unique(d6$Variety)[8])
  Variety9_idx=which(d6$Variety==unique(d6$Variety)[9])
  Variety10_idx=which(d6$Variety==unique(d6$Variety)[10])
  Variety11_idx=which(d6$Variety==unique(d6$Variety)[11])
  Variety12_idx=which(d6$Variety==unique(d6$Variety)[12])
  Variety13_idx=which(d6$Variety==unique(d6$Variety)[13])
  Variety14_idx=which(d6$Variety==unique(d6$Variety)[14])
  Variety15_idx=which(d6$Variety==unique(d6$Variety)[15])
  Variety16_idx=which(d6$Variety==unique(d6$Variety)[16])
  Variety17_idx=which(d6$Variety==unique(d6$Variety)[17])
  Variety18_idx=which(d6$Variety==unique(d6$Variety)[18])
  Variety19_idx=which(d6$Variety==unique(d6$Variety)[19])
  Variety20_idx=which(d6$Variety==unique(d6$Variety)[20])
  
  tests=list(Variety1_idx,Variety2_idx,Variety3_idx,Variety4_idx,Variety5_idx,Variety6_idx,Variety7_idx,Variety8_idx,Variety9_idx,Variety10_idx,Variety11_idx,Variety12_idx,Variety13_idx,Variety14_idx,Variety15_idx,Variety16_idx,Variety17_idx,Variety18_idx,Variety19_idx,Variety20_idx)
}


# Make the matrices that goes into the model
{
  CloverDesign <- model.matrix(~0+d6$Clovershort)
  
  GcloverReps <- CloverDesign %*% GRM1 %*% t(CloverDesign) 
  #GcloverReps is a G-matrix that will match the size and lay-out in the data and can be used in BGLR at the K=
  #Add a clover effect to capture non-additve variance that makes broad sense heritability.
  #CloverIndep catch clover-effects without relationships (clovers are independent), that is another matrix that can go into the model. 
  CloverIndep <- CloverDesign %*% t(CloverDesign)
  dim(CloverIndep)
}








GP_GBLUP<-function(testpop){
  
  ################  ################  ################  ################
  ##start by estimating GEBVs for training population individuals 
  ################  ################  ################  ################
  
  d6_training=d6[-testpop,] # limit the dataframe to only the individuals allowed for training the model
  d6_training_ready=na.omit(d6_training, cols = c("gpd_dryweight_cor")) # remember that gpd na inidividuals should be removed whether or not they are in the training pop or not
  
  additional_individuals_to_remove=which(is.na(d6$gpd_dryweight_cor)==TRUE) #no genotypes
  IndividualsToRemove=c(testpop,additional_individuals_to_remove) # remember that gpd na inidividuals should be removed whether or not they are in the training pop or not
  
  ind_not_in_train=d6$Clovershort[IndividualsToRemove]
  IndividualsToRemoveGRM=which(colnames(GRM1) %in% ind_not_in_train)
  
  GRM_trn = GRM1[-IndividualsToRemoveGRM,-IndividualsToRemoveGRM]
  
  GcloverReps_trn = GcloverReps[-IndividualsToRemove,-IndividualsToRemove] # a GRM for the individuals used to predict
  CloverIndep_trn= CloverIndep[-IndividualsToRemove,-IndividualsToRemove] # another GRM for the individuals used to predict
  
  # Run the GBLUP model on full training population to extract GEBVs
  yNA=d6_training_ready$gpd_dryweight_cor
  
  fixedmod=model.matrix(~factor(d6_training_ready$Round)+factor(d6_training_ready$NS)+factor(d6_training_ready$EW)+factor(d6_training_ready$Rhizobium) + factor(d6_training_ready$inoculation_date))
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
  GcloverReps_covar = GcloverReps[testpop,-IndividualsToRemove] # Covariance between training and testing pop.
  
  #GEBVpred_contr1 = GcloverReps_covar%*%solve(GcloverReps_trn) %*% GEBV_contribution1data 
  GEBVpred = GcloverReps_covar%*%ginv(GcloverReps_trn) %*% GEBV_contribution1data 
  #GEBVpred_contr1 = GcloverReps_covar%*%solve(GcloverReps_trn + diag(0.01, 1661, 1661)) %*% GEBV_contribution1data 
  
  # Output matrix with prediction results
  matrix1=cbind(as.character(d6$Clovershort[testpop]),as.numeric(d6$CorrectedPheno[testpop]),as.numeric(as.character(GEBVpred)))
  colnames(matrix1)=c("ID", "Observed", "GEBV")
  return(matrix1)
}

# print the d6 table
#write.table(d6,"JustTheWhole_d6Table.txt", col.names=T,row.names=F, quote=F,sep="\t")


#  Calculate heritability
print("Calculating heritability")

d6_without_NA=na.omit(d6, cols = c("gpd_dryweight_cor")) # remove individuals with no gpd information from data
d6_without_NA$Clover=droplevels(d6_without_NA$Clover) # removing levels not used in actual data

y=d6_without_NA$gpd_dryweight_cor

remove=which(d6$Clovershort %in% c("Aoost_03","Aoost_06","Aaran_06","Kdike_10"))
GcloverReps_=GcloverReps[-remove,-remove]
CloverIndep_=CloverIndep[-remove,-remove]

fixedmod=model.matrix(~factor(d6_without_NA$Round)+factor(d6_without_NA$NS)+factor(d6_without_NA$EW)+factor(d6_without_NA$Rhizobium)+factor(d6_without_NA$inoculation_date))
ETA=list(list(K=GcloverReps_,model="RKHS"),list(K=CloverIndep_,model="RKHS"),list(X=fixedmod,model="FIXED")) 
GBLUP=BGLR(y=y,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("heritabilityestimation",round))


h2=(GBLUP$ETA[[1]]$varU)/(GBLUP$varE+GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)
h2_=(GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)/(GBLUP$varE+GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)
write.table(h2_,paste("h2_includingNonAdditive",round,".txt",sep=""),sep="\t")
write.table(h2,paste("h2",round,".txt",sep=""),sep="\t")


# Prediction
print("Starting GBLUP prediction")

results=mclapply(tests,GP_GBLUP)

GEBVsVariety1=results[[1]]
print(paste("variety 1 is",unique(d6$Variety)[1]))
GEBVsVariety2=results[[2]]
print(paste("variety 2 is",unique(d6$Variety)[2]))
GEBVsVariety3=results[[3]]
print(paste("variety 3 is",unique(d6$Variety)[3]))
GEBVsVariety4=results[[4]]
print(paste("variety 4 is",unique(d6$Variety)[4]))
GEBVsVariety5=results[[5]]
print(paste("variety 5 is",unique(d6$Variety)[5]))
GEBVsVariety6=results[[6]]
print(paste("variety 6 is",unique(d6$Variety)[6]))
GEBVsVariety7=results[[7]]
print(paste("variety 7 is",unique(d6$Variety)[7]))
GEBVsVariety8=results[[8]]
print(paste("variety 8 is",unique(d6$Variety)[8]))
GEBVsVariety9=results[[9]]
print(paste("variety 9 is",unique(d6$Variety)[9]))
GEBVsVariety10=results[[10]]
print(paste("variety 10 is",unique(d6$Variety)[10]))
GEBVsVariety11=results[[11]]
print(paste("variety 11 is",unique(d6$Variety)[11]))
GEBVsVariety12=results[[12]]
print(paste("variety 12 is",unique(d6$Variety)[12]))
GEBVsVariety13=results[[13]]
print(paste("variety 13 is",unique(d6$Variety)[13]))
GEBVsVariety14=results[[14]]
print(paste("variety 14 is",unique(d6$Variety)[14]))
GEBVsVariety15=results[[15]]
print(paste("variety 15 is",unique(d6$Variety)[15]))
GEBVsVariety16=results[[16]]
print(paste("variety 16 is",unique(d6$Variety)[16]))
GEBVsVariety17=results[[17]]
print(paste("variety 17 is",unique(d6$Variety)[17]))
GEBVsVariety18=results[[18]]
print(paste("variety 18 is",unique(d6$Variety)[18]))
GEBVsVariety19=results[[19]]
print(paste("variety 19 is",unique(d6$Variety)[19]))
GEBVsVariety20=results[[20]]
print(paste("variety 20 is",unique(d6$Variety)[20]))

# pop 1
GEBVsVariety1_table1=aggregate(as.numeric(GEBVsVariety1[, 2]), list(GEBVsVariety1[,1]), mean) #observed a mean for each individual
GEBVsVariety1_table2=aggregate(as.numeric(GEBVsVariety1[, 3]), list(GEBVsVariety1[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined1=cbind(GEBVsVariety1_table1[,1],GEBVsVariety1_table1[,2],GEBVsVariety1_table2[,2])
Pre_Combined1=as.data.frame(Pre_Combined1)
colnames(Pre_Combined1)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined1)

# pop 2
GEBVsVariety2_table1=aggregate(as.numeric(GEBVsVariety2[, 2]), list(GEBVsVariety2[,1]), mean) #observed a mean for each individual
GEBVsVariety2_table2=aggregate(as.numeric(GEBVsVariety2[, 3]), list(GEBVsVariety2[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined2=cbind(GEBVsVariety2_table1[,1],GEBVsVariety2_table1[,2],GEBVsVariety2_table2[,2])
Pre_Combined2=as.data.frame(Pre_Combined2)
colnames(Pre_Combined2)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined2)

# pop 3
GEBVsVariety3_table1=aggregate(as.numeric(GEBVsVariety3[, 2]), list(GEBVsVariety3[,1]), mean) #observed a mean for each individual
GEBVsVariety3_table2=aggregate(as.numeric(GEBVsVariety3[, 3]), list(GEBVsVariety3[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined3=cbind(GEBVsVariety3_table1[,1],GEBVsVariety3_table1[,2],GEBVsVariety3_table2[,2])
Pre_Combined3=as.data.frame(Pre_Combined3)
colnames(Pre_Combined3)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined3)

#  pop 4
GEBVsVariety4_table1=aggregate(as.numeric(GEBVsVariety4[, 2]), list(GEBVsVariety4[,1]), mean) #observed a mean for each individual
GEBVsVariety4_table2=aggregate(as.numeric(GEBVsVariety4[, 3]), list(GEBVsVariety4[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined4=cbind(GEBVsVariety4_table1[,1],GEBVsVariety4_table1[,2],GEBVsVariety4_table2[,2])
Pre_Combined4=as.data.frame(Pre_Combined4)
colnames(Pre_Combined4)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined4)

# pop 5
GEBVsVariety5_table1=aggregate(as.numeric(GEBVsVariety5[, 2]), list(GEBVsVariety5[,1]), mean) #observed a mean for each individual
GEBVsVariety5_table2=aggregate(as.numeric(GEBVsVariety5[, 3]), list(GEBVsVariety5[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined5=cbind(GEBVsVariety5_table1[,1],GEBVsVariety5_table1[,2],GEBVsVariety5_table2[,2])
Pre_Combined5=as.data.frame(Pre_Combined5)
colnames(Pre_Combined5)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined5)

# pop 6
GEBVsVariety6_table1=aggregate(as.numeric(GEBVsVariety6[, 2]), list(GEBVsVariety6[,1]), mean) #observed a mean for each individual
GEBVsVariety6_table2=aggregate(as.numeric(GEBVsVariety6[, 3]), list(GEBVsVariety6[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined6=cbind(GEBVsVariety6_table1[,1],GEBVsVariety6_table1[,2],GEBVsVariety6_table2[,2])
Pre_Combined6=as.data.frame(Pre_Combined6)
colnames(Pre_Combined6)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined6)

# pop 7
GEBVsVariety7_table1=aggregate(as.numeric(GEBVsVariety7[, 2]), list(GEBVsVariety7[,1]), mean) #observed a mean for each individual
GEBVsVariety7_table2=aggregate(as.numeric(GEBVsVariety7[, 3]), list(GEBVsVariety7[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined7=cbind(GEBVsVariety7_table1[,1],GEBVsVariety7_table1[,2],GEBVsVariety7_table2[,2])
Pre_Combined7=as.data.frame(Pre_Combined7)
colnames(Pre_Combined7)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined7)


# pop 8
GEBVsVariety8_table1=aggregate(as.numeric(GEBVsVariety8[, 2]), list(GEBVsVariety8[,1]), mean) #observed a mean for each individual
GEBVsVariety8_table2=aggregate(as.numeric(GEBVsVariety8[, 3]), list(GEBVsVariety8[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined8=cbind(GEBVsVariety8_table1[,1],GEBVsVariety8_table1[,2],GEBVsVariety8_table2[,2])
Pre_Combined8=as.data.frame(Pre_Combined8)
colnames(Pre_Combined8)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined8)

# pop 9
GEBVsVariety9_table1=aggregate(as.numeric(GEBVsVariety9[, 2]), list(GEBVsVariety9[,1]), mean) #observed a mean for each individual
GEBVsVariety9_table2=aggregate(as.numeric(GEBVsVariety9[, 3]), list(GEBVsVariety9[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined9=cbind(GEBVsVariety9_table1[,1],GEBVsVariety9_table1[,2],GEBVsVariety9_table2[,2])
Pre_Combined9=as.data.frame(Pre_Combined9)
colnames(Pre_Combined9)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined9)

# pop 10
GEBVsVariety10_table1=aggregate(as.numeric(GEBVsVariety10[, 2]), list(GEBVsVariety10[,1]), mean) #observed a mean for each individual
GEBVsVariety10_table2=aggregate(as.numeric(GEBVsVariety10[, 3]), list(GEBVsVariety10[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined10=cbind(GEBVsVariety10_table1[,1],GEBVsVariety10_table1[,2],GEBVsVariety10_table2[,2])
Pre_Combined10=as.data.frame(Pre_Combined10)
colnames(Pre_Combined10)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined10)

# pop 11
GEBVsVariety11_table1=aggregate(as.numeric(GEBVsVariety11[, 2]), list(GEBVsVariety11[,1]), mean) #observed a mean for each individual
GEBVsVariety11_table2=aggregate(as.numeric(GEBVsVariety11[, 3]), list(GEBVsVariety11[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined11=cbind(GEBVsVariety11_table1[,1],GEBVsVariety11_table1[,2],GEBVsVariety11_table2[,2])
Pre_Combined11=as.data.frame(Pre_Combined11)
colnames(Pre_Combined11)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined11)

# pop 12
GEBVsVariety12_table1=aggregate(as.numeric(GEBVsVariety12[, 2]), list(GEBVsVariety12[,1]), mean) #observed a mean for each individual
GEBVsVariety12_table2=aggregate(as.numeric(GEBVsVariety12[, 3]), list(GEBVsVariety12[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined12=cbind(GEBVsVariety12_table1[,1],GEBVsVariety12_table1[,2],GEBVsVariety12_table2[,2])
Pre_Combined12=as.data.frame(Pre_Combined12)
colnames(Pre_Combined12)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined12)

# pop 13
GEBVsVariety13_table1=aggregate(as.numeric(GEBVsVariety13[, 2]), list(GEBVsVariety13[,1]), mean) #observed a mean for each individual
GEBVsVariety13_table2=aggregate(as.numeric(GEBVsVariety13[, 3]), list(GEBVsVariety13[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined13=cbind(GEBVsVariety13_table1[,1],GEBVsVariety13_table1[,2],GEBVsVariety13_table2[,2])
Pre_Combined13=as.data.frame(Pre_Combined13)
colnames(Pre_Combined13)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined13)

# pop 14
GEBVsVariety14_table1=aggregate(as.numeric(GEBVsVariety14[, 2]), list(GEBVsVariety14[,1]), mean) #observed a mean for each individual
GEBVsVariety14_table2=aggregate(as.numeric(GEBVsVariety14[, 3]), list(GEBVsVariety14[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined14=cbind(GEBVsVariety14_table1[,1],GEBVsVariety14_table1[,2],GEBVsVariety14_table2[,2])
Pre_Combined14=as.data.frame(Pre_Combined14)
colnames(Pre_Combined14)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined14)

# pop 15
GEBVsVariety15_table1=aggregate(as.numeric(GEBVsVariety15[, 2]), list(GEBVsVariety15[,1]), mean) #observed a mean for each individual
GEBVsVariety15_table2=aggregate(as.numeric(GEBVsVariety15[, 3]), list(GEBVsVariety15[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined15=cbind(GEBVsVariety15_table1[,1],GEBVsVariety15_table1[,2],GEBVsVariety15_table2[,2])
Pre_Combined15=as.data.frame(Pre_Combined15)
colnames(Pre_Combined15)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined15)

# pop 16
GEBVsVariety16_table1=aggregate(as.numeric(GEBVsVariety16[, 2]), list(GEBVsVariety16[,1]), mean) #observed a mean for each individual
GEBVsVariety16_table2=aggregate(as.numeric(GEBVsVariety16[, 3]), list(GEBVsVariety16[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined16=cbind(GEBVsVariety16_table1[,1],GEBVsVariety16_table1[,2],GEBVsVariety16_table2[,2])
Pre_Combined16=as.data.frame(Pre_Combined16)
colnames(Pre_Combined16)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined16)

# pop 17
GEBVsVariety17_table1=aggregate(as.numeric(GEBVsVariety17[, 2]), list(GEBVsVariety17[,1]), mean) #observed a mean for each individual
GEBVsVariety17_table2=aggregate(as.numeric(GEBVsVariety17[, 3]), list(GEBVsVariety17[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined17=cbind(GEBVsVariety17_table1[,1],GEBVsVariety17_table1[,2],GEBVsVariety17_table2[,2])
Pre_Combined17=as.data.frame(Pre_Combined17)
colnames(Pre_Combined17)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined17)

# pop 18
GEBVsVariety18_table1=aggregate(as.numeric(GEBVsVariety18[, 2]), list(GEBVsVariety18[,1]), mean) #observed a mean for each individual
GEBVsVariety18_table2=aggregate(as.numeric(GEBVsVariety18[, 3]), list(GEBVsVariety18[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined18=cbind(GEBVsVariety18_table1[,1],GEBVsVariety18_table1[,2],GEBVsVariety18_table2[,2])
Pre_Combined18=as.data.frame(Pre_Combined18)
colnames(Pre_Combined18)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined18)

# pop 19
GEBVsVariety19_table1=aggregate(as.numeric(GEBVsVariety19[, 2]), list(GEBVsVariety19[,1]), mean) #observed a mean for each individual
GEBVsVariety19_table2=aggregate(as.numeric(GEBVsVariety19[, 3]), list(GEBVsVariety19[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined19=cbind(GEBVsVariety19_table1[,1],GEBVsVariety19_table1[,2],GEBVsVariety19_table2[,2])
Pre_Combined19=as.data.frame(Pre_Combined19)
colnames(Pre_Combined19)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined19)

# pop 20
GEBVsVariety20_table1=aggregate(as.numeric(GEBVsVariety20[, 2]), list(GEBVsVariety20[,1]), mean) #observed a mean for each individual
GEBVsVariety20_table2=aggregate(as.numeric(GEBVsVariety20[, 3]), list(GEBVsVariety20[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined20=cbind(GEBVsVariety20_table1[,1],GEBVsVariety20_table1[,2],GEBVsVariety20_table2[,2])
Pre_Combined20=as.data.frame(Pre_Combined20)
colnames(Pre_Combined20)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")
head(Pre_Combined20)


# combine
Combined=rbind(Pre_Combined1,Pre_Combined2,Pre_Combined3,Pre_Combined4,Pre_Combined5,Pre_Combined6,Pre_Combined7,Pre_Combined8,Pre_Combined9,Pre_Combined10,Pre_Combined11,Pre_Combined12,Pre_Combined13,Pre_Combined14,Pre_Combined15,Pre_Combined16,Pre_Combined17,Pre_Combined18,Pre_Combined19,Pre_Combined20)
filename1=paste("Predictions_GBLUP_GPD",round,".txt",sep="")
write.table(Combined,filename1,sep="\t",quote=F)

# correlations
Combined_withoutNAPheno=na.omit(Combined)

All1=aggregate(as.numeric(as.character(Combined_withoutNAPheno[, 2])), list(Combined_withoutNAPheno[,1]), mean) # mean observed
All2=aggregate(as.numeric(as.character(Combined_withoutNAPheno[, 3])), list(Combined_withoutNAPheno[,1]), mean) # mean GEBV
correlation=cor(All1[,2],All2[,2]) #means of replicates
correlation

filename=paste("Correlation_GBLUP_GPD",round,".txt",sep="")
write.table(correlation,filename,sep="\t")

