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
  GRM=read.table("/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/Preparation_of_genotypeFile/GRM_Clover_Fullfiltering_20200728.csv",sep=",",header=T)
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


# Correct GPD for the part of initial size that does not carry genetic information
{
  lm.fit <- lm(d6$InitialSize ~ d6$Clover)
  summary(lm.fit)
  d6$residuals <- round(lm.fit$residuals,2)
  
  fit <- lm(growth_per_day ~ residuals, data=d6) 
  ycorr <- d6$growth_per_day - model.matrix( ~ residuals, data=d6) %*% summary(fit)$coefficients[,1]
  d6$gpd_dryweight_cor <- ycorr #this is the new corrected dry weight
  print(cor(d6$gpd_dryweight_cor,d6$InitialSize))
  print(cor(d6$growth_per_day,d6$InitialSize))

}


# Calculate a phenotype corrected for all fixed effects
{
  Correctedforallfixed <- lmer(gpd_dryweight_cor ~ factor(Round) + factor(NS) + factor(EW) + factor(Rhizobium) + inoculation_date + (1|Clover), data=d6) 
  summary(Correctedforallfixed)
  matrixOfeffects=model.matrix( ~ factor(Round) + factor(NS) + factor(EW)  + factor(Rhizobium) + inoculation_date, data=d6)
  column_to_remove=which(colnames(matrixOfeffects) %in% names(fixef(Correctedforallfixed))==F)
  matrixOfeffects_new=matrixOfeffects[,-column_to_remove]

  ycorr <- d6$gpd_dryweight_cor - matrixOfeffects_new %*% fixef(Correctedforallfixed)
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




#  Calculate clover heritability
{
  print("Calculating clover heritability by using individual data points")
  y=d6$gpd_dryweight_cor
  fixedmod=model.matrix(~factor(d6$Round)+factor(d6$NS)+factor(d6$EW)+factor(d6$Rhizobium)+factor(d6$inoculation_date))
  ETA=list(list(K=GcloverReps,model="RKHS"),list(K=CloverIndep,model="RKHS"),list(X=fixedmod,model="FIXED")) 
  avg_reps = nrow(d6)/length(unique(d6$Clovershort))
  print(paste("Average number of replicates pr. genotype is",avg_reps,sep=" "))
  m <- matrix(0, ncol = 6, nrow = 100)
  m=as.data.frame(m)
  colnames(m)=c("Round","additiveGRMvar","NonadditiveIndepVar","ResVar","h2_singleplot","h2_linemean")
  
  for (i in seq(1:100)){
    GBLUP=BGLR(y=y,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("heritabilityestimation"))
    h2_singleplot=(GBLUP$ETA[[1]]$varU)/(GBLUP$varE+GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)
    m$h2_singleplot[i]=h2_singleplot
    m$additiveGRMvar[i]=GBLUP$ETA[[1]]$varU
    m$NonadditiveIndepVar[i]=GBLUP$ETA[[2]]$varU
    m$ResVar[i]=GBLUP$varE
    h2_linemean=(GBLUP$ETA[[1]]$varU)/(GBLUP$varE/(avg_reps)+GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)
    m$h2_linemean[i]=h2_linemean
    m$Round[i]=i
    }
  
  write.table(m,paste("Clover_heritabilities_fulldata",".txt",sep=""),sep="\t")
}  

{
  print("Calculating clover heritability by using mean of lines")
  ypre=d6$CorrectedPheno #already corrected for all fixed effects
  LineMeans=aggregate(d6$CorrectedPheno, list(d6$Clovershort), mean)
  y=LineMeans[,2]
  ETA=list(list(K=GRM1,model="RKHS")) 
  m2 <- matrix(0, ncol = 4, nrow = 100)
  m2 =as.data.frame(m2)
  colnames(m2)=c("Round","additiveGRMvar","ResVar","h2_linemean")
  
  for (i in seq(1:100)){
    GBLUP=BGLR(y=y,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("heritabilityestimation"))
    h2_linemean=(GBLUP$ETA[[1]]$varU)/(GBLUP$varE+GBLUP$ETA[[1]]$varU)
    m2$h2_linemean[i]=h2_linemean
    m2$Round[i]=i
    m2$additiveGRMvar[i]=GBLUP$ETA[[1]]$varU
    m2$ResVar[i]=GBLUP$varE
  }
  write.table(m2,paste("Clover_h2_linemean_estimatingOnAverages",".txt",sep=""),sep="\t")
}  


# Upload rhiz. matrix
GRM_rhi=read.table("RhizobiumGRM.txt",sep="\t",head=T)
dim(GRM_rhi)
idx_to_keep=which(colnames(GRM_rhi) %in% unique(d6$rhizobium))
rhiz_with_no_pheno=colnames(GRM_rhi)[-idx_to_keep]
idx_to_keep_2=which(unique(d6$rhizobium) %in% colnames(GRM_rhi))
rhiz_with_no_geno=unique(d6$rhizobium)[-idx_to_keep_2]

Rhiz_to_remove=(c(list(rhiz_with_no_pheno), list(as.character(rhiz_with_no_geno))))
Rhiz_to_remove=unlist(Rhiz_to_remove)
d6_new=d6[-which(d6$Rhizobium %in% Rhiz_to_remove),]
d6_new$rhizobium=droplevels(d6_new$rhizobium) # removing levels not used in actual data
GRM_rhi_filt=GRM_rhi[-which(colnames(GRM_rhi) %in% Rhiz_to_remove),-which(colnames(GRM_rhi) %in% Rhiz_to_remove)]
GRM_rhi_filt=data.matrix(GRM_rhi_filt)
all(dim(GRM_rhi_filt)==length(unique(d6_new$rhizobium))) # Check

# Make the matrices that goes into the model
{
  RhizDesign <- model.matrix(~0+d6_new$rhizobium)
  GRhizReps <- RhizDesign %*% GRM_rhi_filt %*% t(RhizDesign) 
  #GcloverReps is a G-matrix that will match the size and layout in the data and can be used in BGLR at the K=
  #Add a clover effect to capture non-additve variance that makes broad sense heritability.
  #CloverIndep catch clover-effects without relationships (clovers are independent), that is another matrix that can go into the model. 
  RhizIndep <- RhizDesign %*% t(RhizDesign)
  dim(RhizIndep)
}




#  Calculate rhizobium heritability
{
  print("Calculating rhizobum heritability by using individual data points")
  y=d6_new$gpd_dryweight_cor
  length(y)
  fixedmod=model.matrix(~factor(d6_new$Round)+factor(d6_new$NS)+factor(d6_new$EW)+factor(d6_new$Clovershort)+factor(d6_new$inoculation_date))
  ETA=list(list(K=GRhizReps,model="RKHS"),list(K=RhizIndep,model="RKHS"),list(X=fixedmod,model="FIXED")) 
  avg_reps = nrow(d6_new)/length(unique(d6_new$rhizobium))
  print(paste("Average number of replicates pr. genotype is",avg_reps,sep=" "))
  m <- matrix(0, ncol = 6, nrow = 100)
  m=as.data.frame(m)
  colnames(m)=c("Round","additiveGRMvar","NonadditiveIndepVar","ResVar","h2_singleplot","h2_linemean")
  
  for (i in seq(1:100)){
    GBLUP=BGLR(y=y,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("heritabilityestimation"))
    h2_singleplot=(GBLUP$ETA[[1]]$varU)/(GBLUP$varE+GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)
    m$h2_singleplot[i]=h2_singleplot
    m$additiveGRMvar[i]=GBLUP$ETA[[1]]$varU
    m$NonadditiveIndepVar[i]=GBLUP$ETA[[2]]$varU
    m$ResVar[i]=GBLUP$varE
    h2_linemean=(GBLUP$ETA[[1]]$varU)/(GBLUP$varE/(avg_reps)+GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)
    m$h2_linemean[i]=h2_linemean
    m$Round[i]=i
  }
  
  write.table(m,paste("Rhiz_heritabilities_fulldata",".txt",sep=""),sep="\t")
}  


# Calculate a phenotype corrected for all fixed effects
{
  Correctedforallfixed_rhiz <- lmer(gpd_dryweight_cor ~ factor(Round) + factor(NS) + factor(EW) + (1|Rhizobium) + inoculation_date + factor(Clover), data=d6_new) 
  summary(Correctedforallfixed_rhiz)
  matrixOfeffects=model.matrix( ~ factor(Round) + factor(NS) + factor(EW)  + factor(Clover) + inoculation_date, data=d6_new)
  column_to_remove=which(colnames(matrixOfeffects) %in% names(fixef(Correctedforallfixed_rhiz))==F)
  matrixOfeffects_new=matrixOfeffects[,-column_to_remove]
  
  ycorr <- d6_new$gpd_dryweight_cor - matrixOfeffects_new %*% fixef(Correctedforallfixed_rhiz)
  d6_new$CorrectedPhenoRhiz <- ycorr 
}  


{
  print("Calculating rhizobium heritability by using mean of lines")
  ypre=d6_new$CorrectedPhenoRhiz #already corrected for all fixed effects
  LineMeans=aggregate(d6_new$CorrectedPhenoRhiz, list(d6_new$rhizobium), mean)
  y=LineMeans[,2]
  ETA=list(list(K=GRM_rhi_filt,model="RKHS")) 
  m2 <- matrix(0, ncol = 4, nrow = 100)
  m2 =as.data.frame(m2)
  colnames(m2)=c("Round","additiveGRMvar","ResVar","h2_linemean")
  
  for (i in seq(1:100)){
    GBLUP=BGLR(y=y,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("heritabilityestimation"))
    h2_linemean=(GBLUP$ETA[[1]]$varU)/(GBLUP$varE+GBLUP$ETA[[1]]$varU)
    m2$h2_linemean[i]=h2_linemean
    m2$Round[i]=i
    m2$additiveGRMvar[i]=GBLUP$ETA[[1]]$varU
    m2$ResVar[i]=GBLUP$varE
  }
  write.table(m2,paste("Rhiz_h2_linemean_estimatingOnAverages",".txt",sep=""),sep="\t")
}  