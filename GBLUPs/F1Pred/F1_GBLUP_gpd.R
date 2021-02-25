####################################################################################
# Prediction of F1 clover generation.      gpd                                     #
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


# Fix data matrix so that names are the same as in GRM, that is 118 individuals in each
{
  GRM=read.table(args[1],sep=",",header=T)
  #GRM=read.table("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/Preparation_of_genotypeFile/GRM_Clover_Fullfiltering_20200728.csv",sep=",",header=T)
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

d6=d4
length(which(d6$Clovershort=="Aearl_07"))
d6$Rhizobium=droplevels(d6$Rhizobium) # removing levels not used in actual data
d6$Clover=droplevels(d6$Clover) # removing levels not used in actual data
nrow(d6)


# make 4 fake rows in d6
fakedata1=c(rep("NA",7),"Aoost_03",rep("NA",8),"Aoost_03",rep("NA",9),"Aoost_03")
fakedata2=c(rep("NA",7),"Aoost_06",rep("NA",8),"Aoost_06",rep("NA",9),"Aoost_06")
fakedata3=c(rep("NA",7),"Aaran_06",rep("NA",8),"Aaran_06",rep("NA",9),"Aaran_06")
fakedata4=c(rep("NA",7),"Kdike_10",rep("NA",8),"Kdike_10",rep("NA",9),"Kdike_10")
d6=rbind(d6,fakedata1,fakedata2,fakedata3,fakedata4)
d6=d6[order(d6$Clovershort),] 
d6$growth_per_day=as.numeric(as.character(d6$growth_per_day))


# Clean up
{
  d6$Clover=droplevels(d6$Clover) # removing levels not used in actual data
  d6=d6[order(d6$Clovershort),] # make sure it is in alphabetic order like the GRM
}

# calculate means of genotypes
{
  gpdmeans=aggregate(d6$growth_per_day, list(d6$Clovershort), mean)
  colnames(gpdmeans)=c("Individual","gpdNoCor")
}


# Divide into F1 parental groups
{
  parentpop1names=c("Aoost_01","Aoost_08","Aoost_09","Banna_02","Banna_03","Banna_07") #DLF1
  parentpop2names=c("Aearl_08","Ccyma_03","Llanc_06","Aaran_08") #LJ-L6
  parentpop3names=c("Aearl_05","Clfin_02","Ctain_05","Mrida_04") #LJ-L7
  parentpop4names=c("Clfin_03","Ctain_05","Volin_01","Aaran_04") #LJ-L9
  parentpop5names=c("Aoost_02","Ilona_09","Llanc_09","Sster_01") #LJ-H1
  parentpop6names=c("Ilona_05","Kdike_09","Llanc_09","Aalon_03") #LJ-H2
  parentpop7names=c("Ancor_10","Borek_06","Ctain_09","Rbani_02") #LJ-H3
  parentpop8names=c("Ancor_04","Aoost_10","Clfin_08","Kdike_08") #LJ-H5
  parentpop9names=c("Aoost_01","Aoost_08","Banna_02","Rbani_02","Sster_01","Sster_06") #SUA
  
  parentpop1names_idx=which(gpdmeans$Individual%in%parentpop1names)
  parentpop2names_idx=which(gpdmeans$Individual%in%parentpop2names)
  parentpop3names_idx=which(gpdmeans$Individual%in%parentpop3names)
  parentpop4names_idx=which(gpdmeans$Individual%in%parentpop4names)
  parentpop5names_idx=which(gpdmeans$Individual%in%parentpop5names)
  parentpop6names_idx=which(gpdmeans$Individual%in%parentpop6names)
  parentpop7names_idx=which(gpdmeans$Individual%in%parentpop7names)
  parentpop8names_idx=which(gpdmeans$Individual%in%parentpop8names)
  parentpop9names_idx=which(gpdmeans$Individual%in%parentpop9names)

  tests=list(parentpop1names_idx,parentpop2names_idx,parentpop3names_idx,parentpop4names_idx,parentpop5names_idx,parentpop6names_idx,parentpop7names_idx,parentpop8names_idx,parentpop9names_idx)
}



GP_GBLUP<-function(testpop){
  
  ################  ################  ################  ################
  ##start by estimating GEBVs for training population individuals 
  ################  ################  ################  ################
  
  gpdmeans_training=gpdmeans[-testpop,] # limit the dataframe to only the individuals allowed for training the model
  print(paste("this is the dimensions of training phenotypes",dim(gpdmeans_training),sep=" "))
  gpdmeans_training_ready_pre=na.omit(gpdmeans_training, cols = c("gpdNoCor")) # do not use individuals with no phenotypes for training
  additional_individuals_to_remove=which(gpdmeans_training_ready_pre$Individual %in% c("Aoost_03","Aoost_06","Aaran_06","Kdike_10")) # do not use individuals with no genotypes for training
  if (length(additional_individuals_to_remove)!=0){
    gpdmeans_training_ready = gpdmeans_training_ready_pre[-additional_individuals_to_remove,]
  }else{
    gpdmeans_training_ready = gpdmeans_training_ready_pre
  }
  print(paste("Using",nrow(gpdmeans_training_ready), "Individuals to train the model",sep=" "))
  
  ind_not_in_train=gpdmeans$Individual[testpop]
  IndividualsToRemoveGRM=which(colnames(GRM1) %in% ind_not_in_train)
  
  
  GRM_trn = GRM1[colnames(GRM1) %in% gpdmeans_training_ready$Individual,colnames(GRM1) %in% gpdmeans_training_ready$Individual]

  # Run the GBLUP model on full training population to extract GEBVs
  yNA=gpdmeans_training_ready$gpdNoCor
  ETA=list(list(K=GRM_trn,model="RKHS"))
  GBLUP=BGLR(y=yNA,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("GBLUP"))
  matrix=cbind(as.character(gpdmeans_training_ready$Individual),as.numeric(gpdmeans_training_ready$gpdNoCor),as.numeric(GBLUP$ETA[[1]]$u))
  colnames(matrix)=c("ID", "Observed", "GEBV")
  GEBV_contribution1data=as.numeric(as.character(matrix[,3]))
  
  ################  ################ 
  ## Now predict testing population 
  ################  ################  
  
  GRMforpred_test = GRM1[testpop,testpop] # GRM for individuals that will be predicted
  GRMforpred_covar = GRM1[testpop,colnames(GRM1) %in% gpdmeans_training_ready$Individual] # Covariance between training and testing pop.
  
  GEBVpred = GRMforpred_covar%*%ginv(GRM_trn) %*% GEBV_contribution1data 

  # Output matrix with prediction results
  matrix1=cbind(as.character(gpdmeans$Individual[testpop]),as.numeric(gpdmeans$gpdNoCor[testpop]),as.numeric(as.character(GEBVpred)))
  colnames(matrix1)=c("ID", "Observed", "GEBV")
  return(matrix1)
}


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

# pop 1
GEBVsParentPop1_table1=aggregate(as.numeric(GEBVsParentPop1[, 2]), list(GEBVsParentPop1[,1]), mean) #observed a mean for each individual
GEBVsParentPop1_table2=aggregate(as.numeric(GEBVsParentPop1[, 3]), list(GEBVsParentPop1[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined1=cbind(GEBVsParentPop1_table1[,1],GEBVsParentPop1_table1[,2],GEBVsParentPop1_table2[,2])
Pre_Combined1=as.data.frame(Pre_Combined1)
colnames(Pre_Combined1)=c("Individual","Observed","GEBVs")

# pop 2
GEBVsParentPop2_table1=aggregate(as.numeric(GEBVsParentPop2[, 2]), list(GEBVsParentPop2[,1]), mean) #observed a mean for each individual
GEBVsParentPop2_table2=aggregate(as.numeric(GEBVsParentPop2[, 3]), list(GEBVsParentPop2[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined2=cbind(GEBVsParentPop2_table1[,1],GEBVsParentPop2_table1[,2],GEBVsParentPop2_table2[,2])
Pre_Combined2=as.data.frame(Pre_Combined2)
colnames(Pre_Combined2)=c("Individual","Observed","GEBVs")

# pop 3
GEBVsParentPop3_table1=aggregate(as.numeric(GEBVsParentPop3[, 2]), list(GEBVsParentPop3[,1]), mean) #observed a mean for each individual
GEBVsParentPop3_table2=aggregate(as.numeric(GEBVsParentPop3[, 3]), list(GEBVsParentPop3[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined3=cbind(GEBVsParentPop3_table1[,1],GEBVsParentPop3_table1[,2],GEBVsParentPop3_table2[,2])
Pre_Combined3=as.data.frame(Pre_Combined3)
colnames(Pre_Combined3)=c("Individual","Observed","GEBVs")

#  pop 4
GEBVsParentPop4_table1=aggregate(as.numeric(GEBVsParentPop4[, 2]), list(GEBVsParentPop4[,1]), mean) #observed a mean for each individual
GEBVsParentPop4_table2=aggregate(as.numeric(GEBVsParentPop4[, 3]), list(GEBVsParentPop4[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined4=cbind(GEBVsParentPop4_table1[,1],GEBVsParentPop4_table1[,2],GEBVsParentPop4_table2[,2])
Pre_Combined4=as.data.frame(Pre_Combined4)
colnames(Pre_Combined4)=c("Individual","Observed","GEBVs")

# pop 5
GEBVsParentPop5_table1=aggregate(as.numeric(GEBVsParentPop5[, 2]), list(GEBVsParentPop5[,1]), mean) #observed a mean for each individual
GEBVsParentPop5_table2=aggregate(as.numeric(GEBVsParentPop5[, 3]), list(GEBVsParentPop5[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined5=cbind(GEBVsParentPop5_table1[,1],GEBVsParentPop5_table1[,2],GEBVsParentPop5_table2[,2])
Pre_Combined5=as.data.frame(Pre_Combined5)
colnames(Pre_Combined5)=c("Individual","Observed","GEBVs")

# pop 6
GEBVsParentPop6_table1=aggregate(as.numeric(GEBVsParentPop6[, 2]), list(GEBVsParentPop6[,1]), mean) #observed a mean for each individual
GEBVsParentPop6_table2=aggregate(as.numeric(GEBVsParentPop6[, 3]), list(GEBVsParentPop6[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined6=cbind(GEBVsParentPop6_table1[,1],GEBVsParentPop6_table1[,2],GEBVsParentPop6_table2[,2])
Pre_Combined6=as.data.frame(Pre_Combined6)
colnames(Pre_Combined6)=c("Individual","Observed","GEBVs")

# pop 7
GEBVsParentPop7_table1=aggregate(as.numeric(GEBVsParentPop7[, 2]), list(GEBVsParentPop7[,1]), mean) #observed a mean for each individual
GEBVsParentPop7_table2=aggregate(as.numeric(GEBVsParentPop7[, 3]), list(GEBVsParentPop7[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined7=cbind(GEBVsParentPop7_table1[,1],GEBVsParentPop7_table1[,2],GEBVsParentPop7_table2[,2])
Pre_Combined7=as.data.frame(Pre_Combined7)
colnames(Pre_Combined7)=c("Individual","Observed","GEBVs")

# pop 8
GEBVsParentPop8_table1=aggregate(as.numeric(GEBVsParentPop8[, 2]), list(GEBVsParentPop8[,1]), mean) #observed a mean for each individual
GEBVsParentPop8_table2=aggregate(as.numeric(GEBVsParentPop8[, 3]), list(GEBVsParentPop8[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined8=cbind(GEBVsParentPop8_table1[,1],GEBVsParentPop8_table1[,2],GEBVsParentPop8_table2[,2])
Pre_Combined8=as.data.frame(Pre_Combined8)
colnames(Pre_Combined8)=c("Individual","Observed","GEBVs")

# pop 9
GEBVsParentPop9_table1=aggregate(as.numeric(GEBVsParentPop9[, 2]), list(GEBVsParentPop9[,1]), mean) #observed a mean for each individual
GEBVsParentPop9_table2=aggregate(as.numeric(GEBVsParentPop9[, 3]), list(GEBVsParentPop9[,1]), mean) #total GEBV, a mean for each individual

Pre_Combined9=cbind(GEBVsParentPop9_table1[,1],GEBVsParentPop9_table1[,2],GEBVsParentPop9_table2[,2])
Pre_Combined9=as.data.frame(Pre_Combined9)
colnames(Pre_Combined9)=c("Individual","Observed","GEBVs")



# combine
Combined=rbind(Pre_Combined1,Pre_Combined2,Pre_Combined3,Pre_Combined4,Pre_Combined5,Pre_Combined6,Pre_Combined7,Pre_Combined8,Pre_Combined9)
filename1="Predictions_GBLUP_GPD.txt"
write.table(Combined,filename1,sep="\t",quote=F)


# Calculating parental average GEBV
F1observations=read.table("F1results.csv",sep=",",header=T, stringsAsFactors = F)
F1observations=as.data.frame(F1observations)

MeanParentalGEBVs=c(mean(GEBVsParentPop1_table2[,2]),mean(GEBVsParentPop2_table2[,2]),mean(GEBVsParentPop3_table2[,2]),mean(GEBVsParentPop4_table2[,2]),mean(GEBVsParentPop5_table2[,2]),mean(GEBVsParentPop6_table2[,2]),mean(GEBVsParentPop7_table2[,2]),mean(GEBVsParentPop8_table2[,2]),mean(GEBVsParentPop9_table2[,2]))
F1_prediction_table=cbind(F1observations$Population,MeanParentalGEBVs,F1observations$X20200226_dryweight,F1observations$X20200212_freshweight)
F1_prediction_table=as.data.frame(F1_prediction_table)
colnames(F1_prediction_table)=c("F1 population","Mean Parental predicted GEBV","F1 pop mean dryweight","F1 pop mean freshweight")
filename2="Correlations_GBLUP_GPD.txt"
write.table(F1_prediction_table,filename2,sep="\t",quote=F)