



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
    
    
    d7_training=d7[-testpop,] # limit the dataframe to only the individuals allowed for training the model

    ind_not_in_train=d7$Clover[testpop]
    IndividualsToRemoveGRM=which(colnames(GcloverReps) %in% ind_not_in_train)
    
    GcloverReps_trn = GcloverReps[-IndividualsToRemoveGRM,-IndividualsToRemoveGRM]
    
    # Run the GBLUP model on full training population to extract GEBVs
    yNA=d7_training$gpdCor
    fixedmod=model.matrix(~factor(d7_training$NS)+factor(d7_training$EW)+factor(d7_training$Rhizobium)+factor(d7_training$Inoculation.date))
    ETA=list(list(K=GcloverReps_trn,model="RKHS"),list(X=fixedmod,model="FIXED"))
    GBLUP=BGLR(y=yNA,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("GBLUP",round))
    matrix=cbind(as.character(d7_training$Clover),as.numeric(d7_training$gpdCor),as.numeric(GBLUP$ETA[[1]]$u))
    colnames(matrix)=c("ID", "Observed", "GEBV")
    GEBV=as.numeric(as.character(matrix[,3]))
    
    
    ################  ################ 
    ## Now predict testing population 
    ################  ################  
    
    GRMforpred_test = GcloverReps[testpop,testpop] # GRM for individuals that will be predicted
    GRMforpred_covar = GcloverReps[testpop,-testpop] # Covariance between training and testing pop.
    
    GEBVpred = GRMforpred_covar%*%ginv(GcloverReps_trn) %*% GEBV 

    # Output matrix with prediction results
    matrix1=cbind(as.character(d7$Clover[testpop]),as.numeric(d7$gpdCor[testpop]),as.numeric(as.character(GEBVpred)))
    colnames(matrix1)=c("ID", "Observed", "GEBV")
    return(matrix1)
  }
  
  GP_GBLUP_onBLUPs<-function(testpop){
    
    ################  ################  ################  ################
    ##start by estimating GEBVs for training population individuals 
    ################  ################  ################  ################
    
    
    PseudoPhenotypes_training=PseudoPhenotypes[-testpop,] # limit the dataframe to only the individuals allowed for training the model
    
    ind_not_in_train=PseudoPhenotypes$Clover[testpop]
    IndividualsToRemoveGRM=which(colnames(GRM1) %in% ind_not_in_train)
    
    GRM1_trn = GRM1[-IndividualsToRemoveGRM,-IndividualsToRemoveGRM]
    
    # Run the GBLUP model on full training population to extract GEBVs
    yNA=PseudoPhenotypes_training$BLUPs_gpdCor_AllGenetic
    ETA=list(list(K=GRM1_trn,model="RKHS"))
    GBLUP=BGLR(y=yNA,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("GBLUP",round))
    matrix=cbind(as.character(PseudoPhenotypes_training$Clover),as.numeric(PseudoPhenotypes_training$BLUPs_gpdCor_AllGenetic),as.numeric(GBLUP$ETA[[1]]$u))
    colnames(matrix)=c("ID", "Observed", "GEBV")
    GEBV=as.numeric(as.character(matrix[,3]))
    
    
    ################  ################ 
    ## Now predict testing population 
    ################  ################  
    
    GRMforpred_test = GRM1[testpop,testpop] # GRM for individuals that will be predicted
    GRMforpred_covar = GRM1[testpop,-testpop] # Covariance between training and testing pop.
    
    GEBVpred = GRMforpred_covar%*%ginv(GRM1_trn) %*% GEBV 
    
    # Output matrix with prediction results
    matrix1=cbind(as.character(PseudoPhenotypes$Clover[testpop]),as.numeric(PseudoPhenotypes$BLUPs_gpdCor_AllGenetic[testpop]),as.numeric(as.character(GEBVpred)))
    colnames(matrix1)=c("ID", "Observed", "GEBV")
    return(matrix1)
  }
}

# Load data
{
  d6 <- read.csv("FullData.csv", header = TRUE, sep = ",")
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
  GRM=read.table(args[1],sep=",",header=T)
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

#BLUPs 
{
  PseudoPhenotypes=read.table("Clover_BLUPs.txt",head=T)
}


# Setting up 6-fold CV system
{
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
    
    testpop1_idx=which(d7$Clover %in% testpop1)
    testpop2_idx=which(d7$Clover %in% testpop2)
    testpop3_idx=which(d7$Clover %in% testpop3)
    testpop4_idx=which(d7$Clover %in% testpop4)
    testpop5_idx=which(d7$Clover %in% testpop5)
    testpop6_idx=which(d7$Clover %in% testpop6)
    
    testpop1_idx_=which(PseudoPhenotypes$Clover %in% testpop1)
    testpop2_idx_=which(PseudoPhenotypes$Clover %in% testpop2)
    testpop3_idx_=which(PseudoPhenotypes$Clover %in% testpop3)
    testpop4_idx_=which(PseudoPhenotypes$Clover %in% testpop4)
    testpop5_idx_=which(PseudoPhenotypes$Clover %in% testpop5)
    testpop6_idx_=which(PseudoPhenotypes$Clover %in% testpop6)
    
    tests_pseudo=list(testpop1_idx_,testpop2_idx_,testpop3_idx_,testpop4_idx_,testpop5_idx_,testpop6_idx_)
    
    
    tests=list(testpop1_idx,testpop2_idx,testpop3_idx,testpop4_idx,testpop5_idx,testpop6_idx)
  }

# Prediction
{
  print("Starting GBLUP prediction")
  results=mclapply(tests,GP_GBLUP)  
  results_2=mclapply(tests_pseudo,GP_GBLUP_onBLUPs)  
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
  correlation=cor(as.numeric(as.character(All1[,2])),as.numeric(as.character(All2[,2]))) #means of replicates
  correlation
  Combined=cbind(All1[,1],All1[,2],All2[,2])
  colnames(Combined)=c("Individual","Observed_correctedForFixed","GEBVs")
  
  filename=paste("Correlation_GBLUPwithfixedeffects_gpdCor",round,".txt",sep="")
  write.table(correlation,filename,sep="\t")
  
  filename1=paste("Predictions_GBLUPwithfixedeffects_gpdCor",round,".txt",sep="")
  write.table(Combined,filename1,sep="\t")
  
} 

# Summarize and make files with results
{
  first=results_2[[1]]
  second=results_2[[2]]
  third=results_2[[3]]
  fourth=results_2[[4]]
  fifth=results_2[[5]]
  sixth=results_2[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)
  correlation2=cor(as.numeric(as.character(All[,2])),as.numeric(as.character(All[,3])))

  filename=paste("Correlation_GBLUPonBLUPs_gpdCor",round,".txt",sep="")
  write.table(correlation2,filename,sep="\t")
  
  filename1=paste("Predictions_GBLUPonBLUPss_gpdCor",round,".txt",sep="")
  write.table(All,filename1,sep="\t")
  
}  