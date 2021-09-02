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


# Load data
{
  d6 <- read.csv("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/20210901_Genomic_Prediction_full_model/FullData.csv", header = TRUE, sep = ",")
}

# Sample down SM42 to include 5 clover plants (that is the absolut max seen for a single strain in the data)
SM42=which(d6$Rhizobium=="SM42")
Plants=unique(d6$Clover[SM42])
set.seed(15)

# how many plants is one inoculum normally associated with? 
sample=sample(Plants,5)
print(sample) # Aoost_02 Llanc_04 Aoost_04 Aalon_05 Llanc_05

remove=which((d6$Clover[SM42] %in% sample)==FALSE)
d7=d6[-SM42[remove],]
nrow(d7)


d7$CloverRhiz=paste(d7$Clover,d7$Rhizobium,sep=":")


# Clean up
{
  d7$Rhizobium=droplevels(d7$Rhizobium) # removing levels not used in actual data
  d7$Clover=droplevels(d7$Clover) # removing levels not used in actual data
  d7=d7[order(d7$Clover),] # make sure it is in alphabetic order like the GRM
}

# Load genomic relationship matrix and make sure clover genotypes match data
{
  GRM=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Heritabilities_20201217/GRM_Clover_Fullfiltering_20200728.csv",sep=",",header=T)
  dim(GRM)
  GRM1=GRM[which(colnames(GRM) %in% d7$Clover),which(colnames(GRM) %in% d7$Clover)]
  dim(GRM1)
  nrow(GRM1)==length(unique(d7$Clover)) #check
  GRM1=data.matrix(GRM1)
  all(colnames(GRM1)==unique(d7$Clover)) #check
}

# Make the matrices that goes into the model
{
  CloverDesign <- model.matrix(~0+d7$Clover)
  GcloverReps <- CloverDesign %*% GRM1 %*% t(CloverDesign) 
  colnames(GcloverReps)=d7$Clover
}



# Setting up 6-fold CV system
{
  set.seed(NULL)
  #Load groups from file
  f=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/20210901_Genomic_Prediction_full_model/grouping1.txt",fill = TRUE)
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
  
  testpop1_idx=which(d7$Clover %in% testpop1)
  testpop2_idx=which(d7$Clover %in% testpop2)
  testpop3_idx=which(d7$Clover %in% testpop3)
  testpop4_idx=which(d7$Clover %in% testpop4)
  testpop5_idx=which(d7$Clover %in% testpop5)
  testpop6_idx=which(d7$Clover %in% testpop6)
  
  tests=list(testpop1_idx,testpop2_idx,testpop3_idx,testpop4_idx,testpop5_idx,testpop6_idx)
}


# Get BLUE estimates for pseudo-phenotypes
{
  library(BayzR)
  d7$Clover=as.factor(d7$Clover)
  d7$NS=as.factor(d7$NS)
  d7$EW=as.factor(d7$EW)
  d7$Rhizobium=as.factor(d7$Rhizobium)
  d7$Inoculation.date=as.factor(d7$Inoculation.date)
  d7$CloverRhiz=as.factor(d7$CloverRhiz)
  
  fit = bayz(gpd ~ rn(Clover) + rn(Rhizobium) + rn(CloverRhiz)+ fixf(NS) + fixf(EW) + fixf(Inoculation.date), data=d7,chain=c(20000,5000,10))
  summary(fit)
  BLUPs= fit$Estimates[1851:1995,1]
  
  matrix=cbind(as.character(unique(d7$Clover)),BLUPs)
  colnames(matrix)=c("Clover", "BLUP")
  
  write.table(matrix,"/Volumes/NAT_MBG-PMg/Cathrine/Nchain/20210901_Genomic_Prediction_full_model/CloverBLUPs.txt",sep="\t",quote=F)
  

}