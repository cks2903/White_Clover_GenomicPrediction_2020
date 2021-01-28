
library(ggplot2)
library("BayzR")
library("wesanderson")


setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Heritabilities_20201217/")

# read in remaining info
greenhouseinfo=read.table("2018_weight.csv",sep=";",head=T)
colnames(greenhouseinfo)[1]="Barcode"


# gpd trait

gpd_data=read.table("d6gpd.txt",head=T,stringsAsFactors = F)
head(gpd_data)


gpd_data_rhiz = merge(gpd_data,greenhouseinfo,by="Barcode")

# average
# calculate means of genotypes
{
  gpdmeans=aggregate(gpd_data$growth_per_day, list(gpd_data$Clover), mean)
  colnames(gpdmeans)=c("Individual","gpd")
  
  # means of rhiz genotypes
  gpdmeans_rhiz=aggregate(gpd_data_rhiz$growth_per_day, list(gpd_data_rhiz$rhizobium), mean)
  
}

# Load genomic relationship matrix and make sure clover genotypes match data
{
  GRM=read.table("GRM_Clover_Fullfiltering_20200728.csv",sep=",",header=T)
  dim(GRM)
  gpdmeans$Clovershort <- strtrim(gpdmeans$Individual,8)
  remove=GRM[which(colnames(GRM) %in% gpdmeans$Clovershort), which(colnames(GRM) %in% gpdmeans$Clovershort)]
  print(remove)
  GRM1=GRM[which(colnames(GRM) %in% gpdmeans$Clovershort),which(colnames(GRM) %in% gpdmeans$Clovershort)]
  dim(GRM1)
  nrow(GRM1)==length(unique(gpdmeans$Clovershort))
  GRM1=data.matrix(GRM1)
  length(colnames(GRM1)==unique(gpdmeans$Clovershort))==nrow(GRM1) #check
}

# GRM for rhizobium
{
  GRMrhiz=read.table("RhizobiumGRM.txt",sep="\t",header=T)
  dim(GRMrhiz)
  GRMrhiz0=GRMrhiz[which(colnames(GRMrhiz) %in% gpdmeans_rhiz$Group.1), which(colnames(GRMrhiz) %in% gpdmeans_rhiz$Group.1)]
  remove=gpdmeans_rhiz$Group.1 [which(gpdmeans_rhiz$Group.1 %in% colnames(GRMrhiz0)==FALSE)]
  gpdmeans_rhiz_new=gpdmeans_rhiz[-which(gpdmeans_rhiz$Group.1 %in% remove),]
  dim(GRMrhiz0)
  nrow(GRMrhiz0)==length(unique(gpdmeans_rhiz_new$Group.1))
  colnames(gpdmeans_rhiz_new)=c("rhizobium","gpd")
}



# fit model
{
  Cloverfit_gpd <- bayz(gpd ~ ranf(Clovershort) + ranf(Clovershort,V=GRM1),
                        data = gpdmeans, chain=c(20000, 5000, 10))
  summary(Cloverfit_gpd)
  BLUPs_gpd_nonadditive=Cloverfit_gpd$Estimates[5:149,1]
  BLUPs_gpd_additive=Cloverfit_gpd$Estimates[150:294,1]
  length(BLUPs_gpd_nonadditive)==length(BLUPs_gpd_additive) #check
  
  BLUPs_gpd=BLUPs_gpd_nonadditive+BLUPs_gpd_additive # non-additive plus additive part
  BLUPs_gpd # sorted in alphabetic order 
  BLUPs_gpd=as.data.frame(BLUPs_gpd)
  BLUPs_gpd$Individual = gpdmeans$Individual
  colnames(BLUPs_gpd)[1]="blup"
  
  
  
  gpdmeans_rhiz_new$rhizobium=droplevels(gpdmeans_rhiz_new$rhizobium) # removing levels not used in actual data
  
  Rhizfit_gpd <- bayz(gpd ~ ranf(rhizobium) + ranf(rhizobium,V=GRMrhiz0),
                        data = gpdmeans_rhiz_new, chain=c(20000, 5000, 10))
  
  summary(Rhizfit_gpd)
  BLUPs_gpd_nonadditive_rhiz=Rhizfit_gpd$Estimates[5:169,1]
  BLUPs_gpd_additive_rhiz=Rhizfit_gpd$Estimates[170:334,1]
  length(BLUPs_gpd_additive_rhiz)==length(BLUPs_gpd_nonadditive_rhiz) #check
  
  BLUPs_gpd_rhiz=BLUPs_gpd_additive_rhiz+BLUPs_gpd_nonadditive_rhiz # non-additive plus additive part
  BLUPs_gpd_rhiz # sorted in alphabetic order 
  BLUPs_gpd_rhiz=as.data.frame(BLUPs_gpd_rhiz)
  BLUPs_gpd_rhiz$Individual = gpdmeans_rhiz_new$rhizobium
  colnames(BLUPs_gpd_rhiz)[1]="blup"
  
  }



# gpdfixcor

gpdfix_data=read.table("d6gpdFixcor.txt",head=T,stringsAsFactors = F)
head(gpdfix_data)

# average
# calculate means of genotypes
{
  gpdfixmeans=aggregate(gpdfix_data$gpd_dryweight_cor, list(gpdfix_data$Clover), mean)
  colnames(gpdfixmeans)=c("Individual","gpdfix")
  
  gpd_data_rhiz = merge(gpdfix_data,greenhouseinfo,by="Barcode")
  gpdmeans_rhiz=aggregate(gpd_data_rhiz$gpd_dryweight_cor, list(gpd_data_rhiz$rhizobium), mean)
  
}

# Load genomic relationship matrix and make sure clover genotypes match data
{
  gpdfixmeans$Clovershort <- strtrim(gpdfixmeans$Individual,8)
  remove=GRM[which(colnames(GRM) %in% gpdfixmeans$Clovershort), which(colnames(GRM) %in% gpdfixmeans$Clovershort)]
  print(remove)
  GRM1=GRM[which(colnames(GRM) %in% gpdfixmeans$Clovershort),which(colnames(GRM) %in% gpdfixmeans$Clovershort)]
  dim(GRM1)
  nrow(GRM1)==length(unique(gpdfixmeans$Clovershort))
  GRM1=data.matrix(GRM1)
  length(colnames(GRM1)==unique(gpdfixmeans$Clovershort))==nrow(GRM1) #check
  
  
  GRMrhiz0=GRMrhiz[which(colnames(GRMrhiz) %in% gpdmeans_rhiz$Group.1), which(colnames(GRMrhiz) %in% gpdmeans_rhiz$Group.1)]
  remove=gpdmeans_rhiz$Group.1 [which(gpdmeans_rhiz$Group.1 %in% colnames(GRMrhiz0)==FALSE)]
  gpdmeans_rhiz_new=gpdmeans_rhiz[-which(gpdmeans_rhiz$Group.1 %in% remove),]
  dim(GRMrhiz0)
  nrow(GRMrhiz0)==length(unique(gpdmeans_rhiz_new$Group.1))
  colnames(gpdmeans_rhiz_new)=c("rhizobium","gpdfix")
}



# fit model
{
  Cloverfit_gpd <- bayz(gpdfix ~ ranf(Clovershort) + ranf(Clovershort,V=GRM1),
                        data = gpdfixmeans, chain=c(20000, 5000, 10))
  summary(Cloverfit_gpd)
  BLUPs_gpd_nonadditive=Cloverfit_gpd$Estimates[5:149,1]
  BLUPs_gpd_additive=Cloverfit_gpd$Estimates[150:294,1]
  length(BLUPs_gpd_nonadditive)==length(BLUPs_gpd_additive) #check
  
  BLUPs_gpdFix=BLUPs_gpd_nonadditive+BLUPs_gpd_additive # non-additive plus additive part
  BLUPs_gpdFix # sorted in alphabetic order 
  BLUPs_gpdFix=as.data.frame(BLUPs_gpdFix)
  BLUPs_gpdFix$Individual = gpdfixmeans$Individual
  colnames(BLUPs_gpdFix)[1]="blup"
  
  
  
  gpdmeans_rhiz_new$rhizobium=droplevels(gpdmeans_rhiz_new$rhizobium) # removing levels not used in actual data
  
  Rhizfit_gpd <- bayz(gpdfix ~ ranf(rhizobium) + ranf(rhizobium,V=GRMrhiz0),
                      data = gpdmeans_rhiz_new, chain=c(20000, 5000, 10))
  
  summary(Rhizfit_gpd)
  BLUPs_gpd_nonadditive_rhiz=Rhizfit_gpd$Estimates[5:169,1]
  BLUPs_gpd_additive_rhiz=Rhizfit_gpd$Estimates[170:334,1]
  length(BLUPs_gpd_additive_rhiz)==length(BLUPs_gpd_nonadditive_rhiz) #check
  
  BLUPs_gpdFix_rhiz=BLUPs_gpd_additive_rhiz+BLUPs_gpd_nonadditive_rhiz # non-additive plus additive part
  BLUPs_gpdFix_rhiz # sorted in alphabetic order 
  BLUPs_gpdFix_rhiz=as.data.frame(BLUPs_gpdFix_rhiz)
  BLUPs_gpdFix_rhiz$Individual = gpdmeans_rhiz_new$rhizobium
  colnames(BLUPs_gpdFix_rhiz)[1]="blup"
}



# gpi

gpi_data=read.table("d6gpi.txt",head=T,stringsAsFactors = F)
head(gpi_data)

# average
# calculate means of genotypes
{
  gpi_datameans=aggregate(gpi_data$GPD_in_interval, list(gpi_data$Clover), mean)
  colnames(gpi_datameans)=c("Individual","gpi")
  
  gpd_data_rhiz = merge(gpi_data,greenhouseinfo,by="Barcode")
  gpdmeans_rhiz=aggregate(gpd_data_rhiz$GPD_in_interval, list(gpd_data_rhiz$rhizobium), mean)
}

# Load genomic relationship matrix and make sure clover genotypes match data
{
  gpi_datameans$Clovershort <- strtrim(gpi_datameans$Individual,8)
  remove=GRM[which(colnames(GRM) %in% gpi_datameans$Clovershort), which(colnames(GRM) %in% gpi_datameans$Clovershort)]
  print(remove)
  GRM1=GRM[which(colnames(GRM) %in% gpi_datameans$Clovershort),which(colnames(GRM) %in% gpi_datameans$Clovershort)]
  dim(GRM1)
  nrow(GRM1)==length(unique(gpi_datameans$Clovershort))
  GRM1=data.matrix(GRM1)
  length(colnames(GRM1)==unique(gpi_datameans$Clovershort))==nrow(GRM1) #check
  
  
  GRMrhiz0=GRMrhiz[which(colnames(GRMrhiz) %in% gpdmeans_rhiz$Group.1), which(colnames(GRMrhiz) %in% gpdmeans_rhiz$Group.1)]
  remove=gpdmeans_rhiz$Group.1 [which(gpdmeans_rhiz$Group.1 %in% colnames(GRMrhiz0)==FALSE)]
  gpdmeans_rhiz_new=gpdmeans_rhiz[-which(gpdmeans_rhiz$Group.1 %in% remove),]
  dim(GRMrhiz0)
  nrow(GRMrhiz0)==length(unique(gpdmeans_rhiz_new$Group.1))
  colnames(gpdmeans_rhiz_new)=c("rhizobium","gpi")
}



# fit model
{
  Cloverfit_gpd <- bayz(gpi ~ ranf(Clovershort) + ranf(Clovershort,V=GRM1),
                        data = gpi_datameans, chain=c(20000, 5000, 10))
  summary(Cloverfit_gpd)
  BLUPs_gpd_nonadditive=Cloverfit_gpd$Estimates[5:149,1]
  BLUPs_gpd_additive=Cloverfit_gpd$Estimates[150:294,1]
  length(BLUPs_gpd_nonadditive)==length(BLUPs_gpd_additive) #check
  
  BLUPs_gpi=BLUPs_gpd_nonadditive+BLUPs_gpd_additive # non-additive plus additive part
  BLUPs_gpi # sorted in alphabetic order 
  BLUPs_gpi=as.data.frame(BLUPs_gpi)
  BLUPs_gpi$Individual = gpi_datameans$Individual
  colnames(BLUPs_gpi)[1]="blup"
  
  
  gpdmeans_rhiz_new$rhizobium=droplevels(gpdmeans_rhiz_new$rhizobium) # removing levels not used in actual data
  
  Rhizfit_gpd <- bayz(gpi ~ ranf(rhizobium) + ranf(rhizobium,V=GRMrhiz0),
                      data = gpdmeans_rhiz_new, chain=c(20000, 5000, 10))
  
  summary(Rhizfit_gpd)
  BLUPs_gpd_nonadditive_rhiz=Rhizfit_gpd$Estimates[5:169,1]
  BLUPs_gpd_additive_rhiz=Rhizfit_gpd$Estimates[170:334,1]
  length(BLUPs_gpd_additive_rhiz)==length(BLUPs_gpd_nonadditive_rhiz) #check
  
  BLUPs_gpi_rhiz=BLUPs_gpd_additive_rhiz+BLUPs_gpd_nonadditive_rhiz # non-additive plus additive part
  BLUPs_gpi_rhiz # sorted in alphabetic order 
  BLUPs_gpi_rhiz=as.data.frame(BLUPs_gpi_rhiz)
  BLUPs_gpi_rhiz$Individual = gpdmeans_rhiz_new$rhizobium
  colnames(BLUPs_gpi_rhiz)[1]="blup"
}





# gpicor

gpiCor_data=read.table("d6gpiCor.txt",head=T,stringsAsFactors = F)
head(gpiCor_data)

# average
# calculate means of genotypes
{
  gpiCor_datameans=aggregate(gpiCor_data$gpi_cor, list(gpiCor_data$Clover), mean)
  colnames(gpiCor_datameans)=c("Individual","gpicor")
  
  
  gpd_data_rhiz = merge(gpiCor_data,greenhouseinfo,by="Barcode")
  gpdmeans_rhiz=aggregate(gpd_data_rhiz$gpi_cor, list(gpd_data_rhiz$rhizobium), mean)
}

# Load genomic relationship matrix and make sure clover genotypes match data
{
  gpiCor_datameans$Clovershort <- strtrim(gpiCor_datameans$Individual,8)
  remove=GRM[which(colnames(GRM) %in% gpiCor_datameans$Clovershort), which(colnames(GRM) %in% gpiCor_datameans$Clovershort)]
  print(remove)
  GRM1=GRM[which(colnames(GRM) %in% gpiCor_datameans$Clovershort),which(colnames(GRM) %in% gpiCor_datameans$Clovershort)]
  dim(GRM1)
  nrow(GRM1)==length(unique(gpiCor_datameans$Clovershort))
  GRM1=data.matrix(GRM1)
  length(colnames(GRM1)==unique(gpiCor_datameans$Clovershort))==nrow(GRM1) #check
  
  
  GRMrhiz0=GRMrhiz[which(colnames(GRMrhiz) %in% gpdmeans_rhiz$Group.1), which(colnames(GRMrhiz) %in% gpdmeans_rhiz$Group.1)]
  remove=gpdmeans_rhiz$Group.1 [which(gpdmeans_rhiz$Group.1 %in% colnames(GRMrhiz0)==FALSE)]
  gpdmeans_rhiz_new=gpdmeans_rhiz[-which(gpdmeans_rhiz$Group.1 %in% remove),]
  dim(GRMrhiz0)
  nrow(GRMrhiz0)==length(unique(gpdmeans_rhiz_new$Group.1))
  colnames(gpdmeans_rhiz_new)=c("rhizobium","gpicor")
}



# fit model
{
  Cloverfit_gpd <- bayz(gpicor ~ ranf(Clovershort) + ranf(Clovershort,V=GRM1),
                        data = gpiCor_datameans, chain=c(20000, 5000, 10))
  summary(Cloverfit_gpd)
  BLUPs_gpd_nonadditive=Cloverfit_gpd$Estimates[5:149,1]
  BLUPs_gpd_additive=Cloverfit_gpd$Estimates[150:294,1]
  length(BLUPs_gpd_nonadditive)==length(BLUPs_gpd_additive) #check
  
  BLUPs_gpiCor=BLUPs_gpd_nonadditive+BLUPs_gpd_additive # non-additive plus additive part
  BLUPs_gpiCor # sorted in alphabetic order 
  BLUPs_gpiCor=as.data.frame(BLUPs_gpiCor)
  BLUPs_gpiCor$Individual = gpiCor_datameans$Individual
  colnames(BLUPs_gpiCor)[1]="blup"
  
  gpdmeans_rhiz_new$rhizobium=droplevels(gpdmeans_rhiz_new$rhizobium) # removing levels not used in actual data
  
  
  Rhizfit_gpd <- bayz(gpicor ~ ranf(rhizobium) + ranf(rhizobium,V=GRMrhiz0),
                      data = gpdmeans_rhiz_new, chain=c(20000, 5000, 10))
  
  summary(Rhizfit_gpd)
  BLUPs_gpd_nonadditive_rhiz=Rhizfit_gpd$Estimates[5:169,1]
  BLUPs_gpd_additive_rhiz=Rhizfit_gpd$Estimates[170:334,1]
  length(BLUPs_gpd_additive_rhiz)==length(BLUPs_gpd_nonadditive_rhiz) #check
  
  BLUPs_gpiCor_rhiz=BLUPs_gpd_additive_rhiz+BLUPs_gpd_nonadditive_rhiz # non-additive plus additive part
  BLUPs_gpiCor_rhiz # sorted in alphabetic order 
  BLUPs_gpiCor_rhiz=as.data.frame(BLUPs_gpiCor_rhiz)
  BLUPs_gpiCor_rhiz$Individual = gpdmeans_rhiz_new$rhizobium
  colnames(BLUPs_gpiCor_rhiz)[1]="blup"
}



# iSize

iSize_data=read.table("d6iSize.txt",head=T,stringsAsFactors = F)
head(iSize_data)

# average
# calculate means of genotypes
{
  iSize_datameans=aggregate(iSize_data$InitialSize, list(iSize_data$Clover), mean)
  colnames(iSize_datameans)=c("Individual","InitialSize")
  
  gpd_data_rhiz = merge(iSize_data,greenhouseinfo,by="Barcode")
  gpdmeans_rhiz=aggregate(gpd_data_rhiz$InitialSize, list(gpd_data_rhiz$rhizobium), mean)
}

# Load genomic relationship matrix and make sure clover genotypes match data
{
  iSize_datameans$Clovershort <- strtrim(iSize_datameans$Individual,8)
  remove=GRM[which(colnames(GRM) %in% iSize_datameans$Clovershort), which(colnames(GRM) %in% iSize_datameans$Clovershort)]
  print(remove)
  GRM1=GRM[which(colnames(GRM) %in% iSize_datameans$Clovershort),which(colnames(GRM) %in% iSize_datameans$Clovershort)]
  dim(GRM1)
  nrow(GRM1)==length(unique(iSize_datameans$Clovershort))
  GRM1=data.matrix(GRM1)
  length(colnames(GRM1)==unique(iSize_datameans$Clovershort))==nrow(GRM1) #check
  
  
  GRMrhiz0=GRMrhiz[which(colnames(GRMrhiz) %in% gpdmeans_rhiz$Group.1), which(colnames(GRMrhiz) %in% gpdmeans_rhiz$Group.1)]
  remove=gpdmeans_rhiz$Group.1 [which(gpdmeans_rhiz$Group.1 %in% colnames(GRMrhiz0)==FALSE)]
  gpdmeans_rhiz_new=gpdmeans_rhiz[-which(gpdmeans_rhiz$Group.1 %in% remove),]
  dim(GRMrhiz0)
  nrow(GRMrhiz0)==length(unique(gpdmeans_rhiz_new$Group.1))
  colnames(gpdmeans_rhiz_new)=c("rhizobium","InitialSize")
}



# fit model
{
  Cloverfit_gpd <- bayz(InitialSize ~ ranf(Clovershort) + ranf(Clovershort,V=GRM1),
                        data = iSize_datameans, chain=c(20000, 5000, 10))
  summary(Cloverfit_gpd)
  BLUPs_gpd_nonadditive=Cloverfit_gpd$Estimates[5:149,1]
  BLUPs_gpd_additive=Cloverfit_gpd$Estimates[150:294,1]
  length(BLUPs_gpd_nonadditive)==length(BLUPs_gpd_additive) #check
  
  BLUPs_iSize=BLUPs_gpd_nonadditive+BLUPs_gpd_additive # non-additive plus additive part
  BLUPs_iSize # sorted in alphabetic order 
  BLUPs_iSize=as.data.frame(BLUPs_iSize)
  BLUPs_iSize$Individual = iSize_datameans$Individual
  colnames(BLUPs_iSize)[1]="blup"
  
  
  
  gpdmeans_rhiz_new$rhizobium=droplevels(gpdmeans_rhiz_new$rhizobium) # removing levels not used in actual data
  
  
  Rhizfit_gpd <- bayz(InitialSize ~ ranf(rhizobium) + ranf(rhizobium,V=GRMrhiz0),
                      data = gpdmeans_rhiz_new, chain=c(20000, 5000, 10))
  
  summary(Rhizfit_gpd)
  BLUPs_gpd_nonadditive_rhiz=Rhizfit_gpd$Estimates[5:169,1]
  BLUPs_gpd_additive_rhiz=Rhizfit_gpd$Estimates[170:334,1]
  length(BLUPs_gpd_additive_rhiz)==length(BLUPs_gpd_nonadditive_rhiz) #check
  
  BLUPs_iSize_rhiz=BLUPs_gpd_additive_rhiz+BLUPs_gpd_nonadditive_rhiz # non-additive plus additive part
  BLUPs_iSize_rhiz # sorted in alphabetic order 
  BLUPs_iSize_rhiz=as.data.frame(BLUPs_iSize_rhiz)
  BLUPs_iSize_rhiz$Individual = gpdmeans_rhiz_new$rhizobium
  colnames(BLUPs_iSize_rhiz)[1]="blup"
}



##################################
##################################
# PLOTTING
##################################
##################################

sp7=ggplot() +
  geom_point(data = BLUPs_gpd, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "#0B775E",cex=5,alpha=0)+
  geom_point(data = BLUPs_gpdFix, aes(x = Individual, y = blup), color= wes_palette("Rushmore1")[2],cex=5,alpha=0.8) +
  geom_point(data = BLUPs_gpi, aes(x = Individual, y = blup), color = "brown4",cex=5,alpha=0.8) +
  geom_point(data = BLUPs_gpiCor, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "cadetblue",cex=5,alpha=0.8) +
  geom_point(data = BLUPs_iSize, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = wes_palette("Rushmore1")[4],cex=5,alpha=0.8) +
  geom_point(data = BLUPs_gpd, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "#0B775E",cex=5,alpha=0.8)+
  
  expand_limits(x = 175) +
  #ylim(c(-0.3,0.3)) +
  xlab("Clover genotype") +
  ylab("BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp7



# plot individually 
sp=ggplot() +
  geom_point(data = BLUPs_gpd, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "#0B775E",cex=5,alpha=0)+
  geom_point(data = BLUPs_gpd, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "#0B775E",cex=3,alpha=0.8)+
  expand_limits(x = 175) +
  ylim(c(-0.3,0.38)) +
  xlab("Clover genotype") +
  ylab("BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/BLUPsClover_gpd",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)

spr=ggplot() +
  geom_point(data = BLUPs_gpd_rhiz, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "#0B775E",cex=5,alpha=0)+
  geom_point(data = BLUPs_gpd_rhiz, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "#0B775E",cex=3,alpha=0.8)+
  expand_limits(x = 175) +
  ylim(c(-0.3,0.38)) +
  xlab("Rhizobium genotype") +
  ylab("BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
spr
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/BLUPsRhiz_gpd",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)






# plot individually 
sp=ggplot() +
  geom_point(data = BLUPs_gpd, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "#0B775E",cex=5,alpha=0)+
  geom_point(data = BLUPs_gpdFix, aes(x = Individual, y = blup), color= wes_palette("Rushmore1")[2],cex=3,alpha=0.8) +
  expand_limits(x = 175) +
  ylim(c(-0.3,0.38)) +
  xlab("Clover genotype") +
  ylab("BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/BLUPsClover_gpdFixCor",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)

spr=ggplot() +
  geom_point(data = BLUPs_gpd_rhiz, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = wes_palette("Rushmore1")[2],cex=5,alpha=0)+
  geom_point(data = BLUPs_gpdFix_rhiz, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = wes_palette("Rushmore1")[2],cex=3,alpha=0.8)+
  expand_limits(x = 175) +
  ylim(c(-0.3,0.38)) +
  xlab("Rhizobium genotype") +
  ylab("BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
spr
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/BLUPsRhiz_gpdFixCor",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)





# plot individually 
sp=ggplot() +
  geom_point(data = BLUPs_gpd, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "#0B775E",cex=5,alpha=0)+
  geom_point(data = BLUPs_gpi, aes(x = Individual, y = blup), color = "brown4",cex=3,alpha=0.8) +
  expand_limits(x = 175) +
  ylim(c(-1300,1600)) +
  xlab("Clover genotype") +
  ylab("BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/BLUPsClover_gpi",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)

spr=ggplot() +
  geom_point(data = BLUPs_gpd_rhiz, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = wes_palette("Rushmore1")[2],cex=5,alpha=0)+
  geom_point(data = BLUPs_gpi_rhiz, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "brown4",cex=3,alpha=0.8)+
  expand_limits(x = 175) +
  ylim(c(-1300,1600)) +
  xlab("Rhizobium genotype") +
  ylab("BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
spr
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/BLUPsRhiz_gpi",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)





# plot individually 
sp=ggplot() +
  geom_point(data = BLUPs_gpd, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "#0B775E",cex=5,alpha=0)+
  geom_point(data = BLUPs_gpiCor, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "cadetblue",cex=3,alpha=0.8) +
  expand_limits(x = 175) +
  ylim(c(-1000,1100)) +
  xlab("Clover genotype") +
  ylab("BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/BLUPsClover_gpiCor",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)


spr=ggplot() +
  geom_point(data = BLUPs_gpd_rhiz, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = wes_palette("Rushmore1")[2],cex=5,alpha=0)+
  geom_point(data = BLUPs_gpiCor_rhiz, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "cadetblue",cex=3,alpha=0.8)+
  expand_limits(x = 175) +
  ylim(c(-1000,1100)) +
  xlab("Rhizobium genotype") +
  ylab("BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
spr
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/BLUPsRhiz_gpiCor",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)





# plot individually 
sp=ggplot() +
  geom_point(data = BLUPs_gpd, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "#0B775E",cex=5,alpha=0)+
  geom_point(data = BLUPs_iSize, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = wes_palette("Rushmore1")[4],cex=3,alpha=0.8) +
  expand_limits(x = 175) +
  ylim(c(-15000,15000)) +
  xlab("Clover genotype") +
  ylab("BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
sp
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/BLUPsClover_iSize",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)


spr=ggplot() +
  geom_point(data = BLUPs_gpd_rhiz, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = wes_palette("Rushmore1")[2],cex=5,alpha=0)+
  geom_point(data = BLUPs_iSize_rhiz, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = wes_palette("Rushmore1")[4],cex=3,alpha=0.8)+
  expand_limits(x = 175) +
  ylim(c(-15000,15000)) +
  xlab("Rhizobium genotype") +
  ylab("BLUPs") +
  theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
spr
ggsave(paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/BLUPsRhiz_iSize",Sys.Date(),".pdf",sep=""), width =60/5, height = 30/3, units = "cm",useDingbats=FALSE)












# can we scale?
#BLUPs_gpd_scaled=BLUPs_gpd
#BLUPs_gpd_scaledcol=scale(BLUPs_gpd$blup)
#BLUPs_gpd_scaled$blup=BLUPs_gpd_scaledcol

#BLUPs_gpdFix_scaled=BLUPs_gpdFix
#BLUPs_gpdFix_scaledcol=scale(BLUPs_gpdFix$blup)
#BLUPs_gpdFix_scaled$blup=BLUPs_gpdFix_scaledcol

#BLUPs_gpi_scaled=BLUPs_gpi
#BLUPs_gpi_scaledcol=scale(BLUPs_gpi$blup)
#BLUPs_gpi_scaled$blup=BLUPs_gpi_scaledcol

#BLUPs_gpiCor_scaled=BLUPs_gpiCor
#BLUPs_gpiCor_scaledcol=scale(BLUPs_gpiCor$blup)
#BLUPs_gpiCor_scaled$blup=BLUPs_gpiCor_scaledcol

#BLUPs_iSize_scaled=BLUPs_iSize
#BLUPs_iSize_scaledcol=scale(BLUPs_iSize$blup)
#BLUPs_iSize_scaled$blup=BLUPs_iSize_scaledcol



#sp8=ggplot() +
#  geom_point(data = BLUPs_gpd_scaled, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "#0B775E",cex=3.5,alpha=0)+
# geom_point(data = BLUPs_gpdFix_scaled, aes(x = Individual, y = blup), color= wes_palette("Rushmore1")[2],cex=3.5,alpha=0.8) +
#  geom_point(data = BLUPs_gpi_scaled, aes(x = Individual, y = blup), color = "brown4",cex=3.5,alpha=0.8) +
 # geom_point(data = BLUPs_gpiCor_scaled, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "cadetblue",cex=3.5,alpha=0.8) +
 # geom_point(data = BLUPs_iSize_scaled, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = wes_palette("Rushmore1")[4],cex=3.5,alpha=0.8) +
#  geom_point(data = BLUPs_gpd_scaled, aes(x = reorder(Individual, blup, FUN =median), y = blup), color = "#0B775E",cex=3.5,alpha=0.8)+
  
 # expand_limits(x = 175) +
  #ylim(c(-0.3,0.3)) +
 # xlab("Clover genotype") +
 # ylab("BLUPs") +
 # theme_classic() +theme(axis.text.x = element_text(angle = 90)) +theme(axis.title=element_text(size=20))
#sp8

