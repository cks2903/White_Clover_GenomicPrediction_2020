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
  
  parentpop1names_idx=which(d6$Clovershort%in%parentpop1names)
  parentpop2names_idx=which(d6$Clovershort%in%parentpop2names)
  parentpop3names_idx=which(d6$Clovershort%in%parentpop3names)
  parentpop4names_idx=which(d6$Clovershort%in%parentpop4names)
  parentpop5names_idx=which(d6$Clovershort%in%parentpop5names)
  parentpop6names_idx=which(d6$Clovershort%in%parentpop6names)
  parentpop7names_idx=which(d6$Clovershort%in%parentpop7names)
  parentpop8names_idx=which(d6$Clovershort%in%parentpop8names)
  parentpop9names_idx=which(d6$Clovershort%in%parentpop9names)
  parentpop10names_idx=which(d6$Clovershort%in%parentpop10names)
  parentpop11names_idx=which(d6$Clovershort%in%parentpop11names)
  parentpop12names_idx=which(d6$Clovershort%in%parentpop12names)
  parentpop13names_idx=which(d6$Clovershort%in%parentpop13names)
  parentpop14names_idx=which(d6$Clovershort%in%parentpop14names)
  
  tests=list(parentpop1names_idx,parentpop2names_idx,parentpop3names_idx,parentpop4names_idx,parentpop5names_idx,parentpop6names_idx,parentpop7names_idx,parentpop8names_idx,parentpop9names_idx,parentpop10names_idx,parentpop11names_idx,parentpop12names_idx,parentpop13names_idx,parentpop14names_idx)
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

GEBVsParentPop1=results[[1]]
GEBVsParentPop2=results[[2]]
GEBVsParentPop3=results[[3]]
GEBVsParentPop4=results[[4]]
GEBVsParentPop5=results[[5]]
GEBVsParentPop6=results[[6]]
GEBVsParentPop7=results[[7]]
GEBVsParentPop8=results[[8]]
GEBVsParentPop9=results[[9]]
GEBVsParentPop10=results[[10]]
GEBVsParentPop11=results[[11]]
GEBVsParentPop12=results[[12]]
GEBVsParentPop13=results[[13]]
GEBVsParentPop14=results[[14]]

# pop 1
GEBVsParentPop1_table1=aggregate(as.numeric(GEBVsParentPop1[, 2]), list(GEBVsParentPop1[,1]), mean) #observed a mean for each individual
GEBVsParentPop1_table2=aggregate(as.numeric(GEBVsParentPop1[, 3]), list(GEBVsParentPop1[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined1=cbind(GEBVsParentPop1_table1[,1],GEBVsParentPop1_table1[,2],GEBVsParentPop1_table2[,2])
Pre_Combined1=as.data.frame(Pre_Combined1)
colnames(Pre_Combined1)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

# pop 2
GEBVsParentPop2_table1=aggregate(as.numeric(GEBVsParentPop2[, 2]), list(GEBVsParentPop2[,1]), mean) #observed a mean for each individual
GEBVsParentPop2_table2=aggregate(as.numeric(GEBVsParentPop2[, 3]), list(GEBVsParentPop2[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined2=cbind(GEBVsParentPop2_table1[,1],GEBVsParentPop2_table1[,2],GEBVsParentPop2_table2[,2])
Pre_Combined2=as.data.frame(Pre_Combined2)
colnames(Pre_Combined2)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

# pop 3
GEBVsParentPop3_table1=aggregate(as.numeric(GEBVsParentPop3[, 2]), list(GEBVsParentPop3[,1]), mean) #observed a mean for each individual
GEBVsParentPop3_table2=aggregate(as.numeric(GEBVsParentPop3[, 3]), list(GEBVsParentPop3[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined3=cbind(GEBVsParentPop3_table1[,1],GEBVsParentPop3_table1[,2],GEBVsParentPop3_table2[,2])
Pre_Combined3=as.data.frame(Pre_Combined3)
colnames(Pre_Combined3)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

#  pop 4
GEBVsParentPop4_table1=aggregate(as.numeric(GEBVsParentPop4[, 2]), list(GEBVsParentPop4[,1]), mean) #observed a mean for each individual
GEBVsParentPop4_table2=aggregate(as.numeric(GEBVsParentPop4[, 3]), list(GEBVsParentPop4[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined4=cbind(GEBVsParentPop4_table1[,1],GEBVsParentPop4_table1[,2],GEBVsParentPop4_table2[,2])
Pre_Combined4=as.data.frame(Pre_Combined4)
colnames(Pre_Combined4)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

# pop 5
GEBVsParentPop5_table1=aggregate(as.numeric(GEBVsParentPop5[, 2]), list(GEBVsParentPop5[,1]), mean) #observed a mean for each individual
GEBVsParentPop5_table2=aggregate(as.numeric(GEBVsParentPop5[, 3]), list(GEBVsParentPop5[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined5=cbind(GEBVsParentPop5_table1[,1],GEBVsParentPop5_table1[,2],GEBVsParentPop5_table2[,2])
Pre_Combined5=as.data.frame(Pre_Combined5)
colnames(Pre_Combined5)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

# pop 6
GEBVsParentPop6_table1=aggregate(as.numeric(GEBVsParentPop6[, 2]), list(GEBVsParentPop6[,1]), mean) #observed a mean for each individual
GEBVsParentPop6_table2=aggregate(as.numeric(GEBVsParentPop6[, 3]), list(GEBVsParentPop6[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined6=cbind(GEBVsParentPop6_table1[,1],GEBVsParentPop6_table1[,2],GEBVsParentPop6_table2[,2])
Pre_Combined6=as.data.frame(Pre_Combined6)
colnames(Pre_Combined6)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

# pop 7
GEBVsParentPop7_table1=aggregate(as.numeric(GEBVsParentPop7[, 2]), list(GEBVsParentPop7[,1]), mean) #observed a mean for each individual
GEBVsParentPop7_table2=aggregate(as.numeric(GEBVsParentPop7[, 3]), list(GEBVsParentPop7[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined7=cbind(GEBVsParentPop7_table1[,1],GEBVsParentPop7_table1[,2],GEBVsParentPop7_table2[,2])
Pre_Combined7=as.data.frame(Pre_Combined7)
colnames(Pre_Combined7)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

# pop 8
GEBVsParentPop8_table1=aggregate(as.numeric(GEBVsParentPop8[, 2]), list(GEBVsParentPop8[,1]), mean) #observed a mean for each individual
GEBVsParentPop8_table2=aggregate(as.numeric(GEBVsParentPop8[, 3]), list(GEBVsParentPop8[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined8=cbind(GEBVsParentPop8_table1[,1],GEBVsParentPop8_table1[,2],GEBVsParentPop8_table2[,2])
Pre_Combined8=as.data.frame(Pre_Combined8)
colnames(Pre_Combined8)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

# pop 9
GEBVsParentPop9_table1=aggregate(as.numeric(GEBVsParentPop9[, 2]), list(GEBVsParentPop9[,1]), mean) #observed a mean for each individual
GEBVsParentPop9_table2=aggregate(as.numeric(GEBVsParentPop9[, 3]), list(GEBVsParentPop9[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined9=cbind(GEBVsParentPop9_table1[,1],GEBVsParentPop9_table1[,2],GEBVsParentPop9_table2[,2])
Pre_Combined9=as.data.frame(Pre_Combined9)
colnames(Pre_Combined9)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

# pop 10
GEBVsParentPop10_table1=aggregate(as.numeric(GEBVsParentPop10[, 2]), list(GEBVsParentPop10[,1]), mean) #observed a mean for each individual
GEBVsParentPop10_table2=aggregate(as.numeric(GEBVsParentPop10[, 3]), list(GEBVsParentPop10[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined10=cbind(GEBVsParentPop10_table1[,1],GEBVsParentPop10_table1[,2],GEBVsParentPop10_table2[,2])
Pre_Combined10=as.data.frame(Pre_Combined10)
colnames(Pre_Combined10)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

# pop 11
GEBVsParentPop11_table1=aggregate(as.numeric(GEBVsParentPop11[, 2]), list(GEBVsParentPop11[,1]), mean) #observed a mean for each individual
GEBVsParentPop11_table2=aggregate(as.numeric(GEBVsParentPop11[, 3]), list(GEBVsParentPop11[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined11=cbind(GEBVsParentPop11_table1[,1],GEBVsParentPop11_table1[,2],GEBVsParentPop11_table2[,2])
Pre_Combined11=as.data.frame(Pre_Combined11)
colnames(Pre_Combined11)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

# pop 12
GEBVsParentPop12_table1=aggregate(as.numeric(GEBVsParentPop12[, 2]), list(GEBVsParentPop12[,1]), mean) #observed a mean for each individual
GEBVsParentPop12_table2=aggregate(as.numeric(GEBVsParentPop12[, 3]), list(GEBVsParentPop12[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined12=cbind(GEBVsParentPop12_table1[,1],GEBVsParentPop12_table1[,2],GEBVsParentPop12_table2[,2])
colnames(Pre_Combined12)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

# pop 13
GEBVsParentPop13_table1=aggregate(as.numeric(GEBVsParentPop13[, 2]), list(GEBVsParentPop13[,1]), mean) #observed a mean for each individual
GEBVsParentPop13_table2=aggregate(as.numeric(GEBVsParentPop13[, 3]), list(GEBVsParentPop13[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined13=cbind(GEBVsParentPop13_table1[,1],GEBVsParentPop13_table1[,2],GEBVsParentPop13_table2[,2])
Pre_Combined13=as.data.frame(Pre_Combined13)
colnames(Pre_Combined13)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")

# pop 14
GEBVsParentPop14_table1=aggregate(as.numeric(GEBVsParentPop14[, 2]), list(GEBVsParentPop14[,1]), mean) #observed a mean for each individual
GEBVsParentPop14_table2=aggregate(as.numeric(GEBVsParentPop14[, 3]), list(GEBVsParentPop14[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined14=cbind(GEBVsParentPop14_table1[,1],GEBVsParentPop14_table1[,2],GEBVsParentPop14_table2[,2])
colnames(Pre_Combined14)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs")




# combine
Combined=rbind(Pre_Combined1,Pre_Combined2,Pre_Combined3,Pre_Combined4,Pre_Combined5,Pre_Combined6,Pre_Combined7,Pre_Combined8,Pre_Combined9,Pre_Combined10,Pre_Combined11,Pre_Combined12,Pre_Combined13,Pre_Combined14)
filename1=paste("Predictions_GBLUP_GPD",round,".txt",sep="")
write.table(Combined,filename1,sep="\t",quote=F)


# Calculating parental average GEBV
F1observations=read.table("F1results.csv",sep=",",header=T, stringsAsFactors = F)
F1observations=as.data.frame(F1observations)

MeanParentalGEBVs=c(mean(GEBVsParentPop1_table2[,2]),mean(GEBVsParentPop2_table2[,2]),mean(GEBVsParentPop3_table2[,2]),mean(GEBVsParentPop4_table2[,2]),mean(GEBVsParentPop5_table2[,2]),mean(GEBVsParentPop6_table2[,2]),mean(GEBVsParentPop7_table2[,2]),mean(GEBVsParentPop8_table2[,2]),mean(GEBVsParentPop9_table2[,2]),mean(GEBVsParentPop10_table2[,2]),mean(GEBVsParentPop11_table2[,2]),mean(GEBVsParentPop12_table2[,2]),mean(GEBVsParentPop13_table2[,2]),mean(GEBVsParentPop14_table2[,2]))
F1_prediction_table=cbind(F1observations$Population,MeanParentalGEBVs,F1observations$X20200226_dryweight,F1observations$X20200212_freshweight)
F1_prediction_table=as.data.frame(F1_prediction_table)
colnames(F1_prediction_table)=c("F1 population","Mean Parental predicted GEBV","F1 pop mean dryweight","F1 pop mean freshweight")
filename2=paste("Correlations_GBLUP_GPD",round,".txt",sep="")
write.table(F1_prediction_table,filename2,sep="\t",quote=F)
