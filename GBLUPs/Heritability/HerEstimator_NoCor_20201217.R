###################################################################################################################
###################################################################################################################
### this is a script to run GBLUP on replicate data with predefined division into training and testing groups   ###
###################################################################################################################
###################################################################################################################
source("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Heritabilities_20201217/HPDbayz.R")

# Load libraries
{
  library(lme4)
  library("parallel")
  library("methods")
  library("Matrix")
  library("MASS")
  library("BayzR")
}


# Load data
{
  d <- read.csv("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Heritabilities_20201217/greenhouse_area.csv", header = TRUE, sep = ",")
  f=read.csv("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Heritabilities_20201217/2018_weight.csv",header=T,sep=";")
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
  remove = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Heritabilities_20201217/Barcodes_removed_based_on_single_Observations_2021-01-06.txt")
  removeidx = which(d$Barcode %in% remove)
  d005 = d[-removeidx,]
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
  GRM=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Heritabilities_20201217/GRM_Clover_Fullfiltering_20200728.csv",sep=",",header=T)
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
  print("Calculating broad sense clover heritability by using individual data points")
  
  Cloverfit <- bayz(growth_per_day ~ ranf(Clover) + fixf(NS) + fixf(EW) + fixf(Rhizobium) + fixf(inoculation_date) ,
               data = d6, chain=c(20000, 5000, 10))
  
  summary(Cloverfit) # 26% broad sense clover

  
  avg_reps = nrow(d6)/length(unique(d6$Clovershort))
  print(paste("Average number of replicates pr. genotype is",avg_reps,sep=" "))
  m <- matrix(0, ncol = 6, nrow = 1)
  m=as.data.frame(m)
  colnames(m)=c("CloverVar","CloverVar_SD","ResVar","ResVar_SD","H2_singleplot","H2_linemean")
  m$CloverVar=Cloverfit$Estimates$postMean[2]
  m$CloverVar_SD=Cloverfit$Estimates$postSD[2]
  m$ResVar=Cloverfit$Estimates$postMean[1]
  m$ResVar_SD=Cloverfit$Estimates$postSD[1]
  m$H2_singleplot= Cloverfit$Estimates$postMean[2]/(Cloverfit$Estimates$postMean[2]+Cloverfit$Estimates$postMean[1])
  m$H2_linemean = Cloverfit$Estimates$postMean[2]/(Cloverfit$Estimates$postMean[2]+Cloverfit$Estimates$postMean[1]/avg_reps)
  
  write.table(m,paste("Clover_heritabilities_fulldata",".txt",sep=""),col.names=T, row.names=F,quote=F,sep="\t")
}  



# Upload rhiz. matrix
GRM_rhi=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Heritabilities_20201217/RhizobiumGRM.txt",sep="\t",head=T)
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
  print("Calculating broad sense rhizobum heritability by using individual data points")
  Rhizfit <- bayz(growth_per_day ~ ranf(Rhizobium) + fixf(NS) + fixf(EW) + fixf(Clover) + fixf(inoculation_date) ,
                    data = d6, chain=c(20000, 5000, 10))
  
  summary(Rhizfit)

  avg_reps = nrow(d6)/length(unique(d6$Rhizobium))
  print(paste("Average number of replicates pr. genotype is",avg_reps,sep=" "))
  m <- matrix(0, ncol = 6, nrow = 1)
  m=as.data.frame(m)
  colnames(m)=c("RhizVar","RhizVar_SD","ResVar","ResVar_SD","H2_singleplot","H2_linemean")
  m$RhizVar=Rhizfit$Estimates$postMean[2]
  m$RhizVar_SD=Rhizfit$Estimates$postSD[2]
  m$ResVar=Rhizfit$Estimates$postMean[1]
  m$ResVar_SD=Rhizfit$Estimates$postSD[1]
  m$H2_singleplot= Rhizfit$Estimates$postMean[2]/(Rhizfit$Estimates$postMean[2]+Rhizfit$Estimates$postMean[1])
  m$H2_linemean = Rhizfit$Estimates$postMean[2]/(Rhizfit$Estimates$postMean[2]+Rhizfit$Estimates$postMean[1]/avg_reps)
  plot(Rhizfit$Samples)
  write.table(m,paste("Rhiz_heritabilities_fulldata",".txt",sep=""),col.names=T, row.names=F,quote=F,sep="\t")
}  



#  Calculate clover heritability
{
  print("Calculating broad sense and narrow sense clover heritability by using genotype averages")
  
  # calculate means of genotypes
  {
    gpdmeans=aggregate(d6$growth_per_day, list(d6$Clovershort), mean)
    colnames(gpdmeans)=c("Individual","gpdNoCor")
  }
  
  Cloverfit2 <- bayz(gpdNoCor ~ ranf(Individual,V=GRM1),
                    data = gpdmeans, chain=c(20000, 5000, 10))
  
  summary(Cloverfit2) 
  HPDbayz(Cloverfit2$Samples) # to get confindence intervals
  
  
  m <- matrix(0, ncol = 7, nrow = 1)
  m=as.data.frame(m)
  colnames(m)=c("CloverAdVar","CloverAdVar_SD","ResVar","ResVar_SD","h2","CloverAdVar_HPD","ResVar_HPD")
  m$CloverAdVar=Cloverfit2$Estimates$postMean[2]
  m$CloverAdVar_SD=Cloverfit2$Estimates$postSD[2]
  m$ResVar=Cloverfit2$Estimates$postMean[1]
  m$ResVar_SD=Cloverfit2$Estimates$postSD[1]
  m$h2 = (Cloverfit2$Estimates$postMean[2])/(Cloverfit2$Estimates$postMean[2]+Cloverfit2$Estimates$postMean[1])
  m$CloverAdVar_HPD=paste(HPDbayz(Cloverfit2$Samples)[2,1],HPDbayz(Cloverfit2$Samples)[2,2],sep=",")
  m$ResVar_HPD=paste(HPDbayz(Cloverfit2$Samples)[1,1],HPDbayz(Cloverfit2$Samples)[1,2],sep=",")
  
  
  write.table(m,paste("Clover_heritabilities_averages",".txt",sep=""),col.names=T, row.names=F,quote=F,sep="\t")
}  


#  Calculate rhiz heritability
{
  print("Calculating broad sense and narrow sense rhiz heritability by using genotype averages")
  
  # calculate means of genotypes
  {
    gpdmeans=aggregate(d6_new$growth_per_day, list(d6_new$Rhizobium), mean)
    colnames(gpdmeans)=c("Individual","gpdNoCor")
    gpdmeans$Individual=droplevels(gpdmeans$Individual) # removing levels not used in actual data
    
  }
  
  Rhizfit2 <- bayz(gpdNoCor ~ ranf(Individual,V=GRM_rhi_filt),
                     data = gpdmeans, chain=c(20000, 5000, 10))
  
  summary(Rhizfit2) 
  m <- matrix(0, ncol = 7, nrow = 1)
  m=as.data.frame(m)
  colnames(m)=c("RhizAdVar","RhizAdVar_SD","ResVar","ResVar_SD","h2","RhizAdVar_HPD","ResVar_HPD")
  m$RhizAdVar=Rhizfit2$Estimates$postMean[3]
  m$RhizAdVar_SD=Rhizfit2$Estimates$postSD[3]
  m$ResVar=Rhizfit2$Estimates$postMean[1]
  m$ResVar_SD=Rhizfit2$Estimates$postSD[1]
  m$h2 = (Rhizfit2$Estimates$postMean[2])/(Rhizfit2$Estimates$postMean[1]+Rhizfit2$Estimates$postMean[2])
  m$RhizAdVar_HPD=paste(HPDbayz(Rhizfit2$Samples)[2,1],HPDbayz(Rhizfit2$Samples)[2,2],sep=",")
  m$ResVar_HPD=paste(HPDbayz(Rhizfit2$Samples)[1,1],HPDbayz(Rhizfit2$Samples)[1,2],sep=",")
  
  write.table(m,paste("Rhiz_heritabilities_averages",".txt",sep=""),col.names=T, row.names=F,quote=F,sep="\t")
}  


FullModel <- bayz(growth_per_day ~ ranf(Clover) + fixf(NS) + fixf(EW) + ranf(Rhizobium) + fixf(inoculation_date) +ran2f(Rhizobium,Clover),
                  data = d6, chain=c(20000, 5000, 10))

summary(FullModel) 
HPDbayz(FullModel$Samples,bound="var")
