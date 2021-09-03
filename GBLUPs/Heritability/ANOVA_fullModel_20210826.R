#library(devtools)
#devtools::install_github("https://github.com/MarniTausen/BayzR/tree/9b02adcd7929a9428d851b75d9ff2dcc2cd56a45") #think this is the right one, 20/01/2021, problems with fit4 and fit5
#devtools::install_github("ljanss/BayzR")
library(BayzR) 
source("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/20210826_ANOVA_Model_InteractionTest_NarrowingDownSM42/HPDbayz.R")


setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/ANOVA_Model_InteractionTest")

#load data, filtered but single observations
df = read.table("SupplementaryFile2.xls - FullData_20210330.csv",sep=",",head=T)
dim(df)
head(df)

# Sample down SM42 to include 5 clover plants (that is the absolut max seen for a single strain in the data)
SM42=which(df$Rhizobium=="SM42")
Plants=unique(df$Clover[SM42])
set.seed(15)

# how many plants is one inoculum normally associated with? 
sample=sample(Plants,5)
print(sample) # Aoost_02 Llanc_04 Aoost_04 Aalon_05 Llanc_05

remove=which((df$Clover[SM42] %in% sample)==FALSE)
d7=df[-SM42[remove],]
nrow(d7)


d7$CloverRhiz=paste(d7$Clover,d7$Rhizobium,sep=":")
d7=d7[which(is.na(d7$gpiCor)==FALSE),]

# Make data ready for Bayz
d7$Clover=as.factor(d7$Clover)
d7$Rhizobium=as.factor(d7$Rhizobium)
d7$NS=as.factor(d7$NS)
d7$EW=as.factor(d7$EW)
d7$Inoculation.date=as.factor(d7$Inoculation.date)
d7$CloverRhiz=as.factor(d7$CloverRhiz)

# Clean up
{
  d7$Rhizobium=droplevels(d7$Rhizobium) # removing levels not used in actual data
  d7$Clover=droplevels(d7$Clover) # removing levels not used in actual data
  d7=d7[order(d7$Clover),] # make sure it is in alphabetic order like the GRM
}


# model
fit1 = bayz(gpd ~ rn(Clover) + fixf(NS) + fixf(EW) + rn(Rhizobium) + fixf(Inoculation.date) +rn(CloverRhiz), data=d7,chain=c(20000,5000,10))
summary(fit1)
HPDbayz(fit1$Samples,bound="var") # to get confindence intervals
varcolumns = ( substr(colnames(fit1$Samples),0,3)=="var" )
proportions = t(apply(fit1$Samples[,varcolumns],1,function(x){x/sum(x)}))
HPDbayz(proportions, bound="prob")


fit2 = bayz(iSize ~ rn(Clover) + fixf(NS) + fixf(EW) + rn(Rhizobium) + fixf(Inoculation.date) +rn(CloverRhiz), data=d7,chain=c(20000,5000,10))
summary(fit2)
HPDbayz(fit2$Samples,bound="var") # to get confindence intervals
varcolumns = ( substr(colnames(fit2$Samples),0,3)=="var" )
proportions = t(apply(fit2$Samples[,varcolumns],1,function(x){x/sum(x)}))
HPDbayz(proportions, bound="prob")

fit3 = bayz(gpdCor ~ rn(Clover) + fixf(NS) + fixf(EW) + rn(Rhizobium) + fixf(Inoculation.date) +rn(CloverRhiz), data=d7,chain=c(20000,5000,10))
summary(fit3)
HPDbayz(fit3$Samples,bound="var") # to get confindence intervals
varcolumns = ( substr(colnames(fit3$Samples),0,3)=="var" )
proportions = t(apply(fit3$Samples[,varcolumns],1,function(x){x/sum(x)}))
HPDbayz(proportions, bound="prob")

fit4 = bayz(gpi ~rn(Clover) + fixf(NS) + fixf(EW) + rn(Rhizobium) + fixf(Inoculation.date) +rn(CloverRhiz), data=d7,chain=c(20000,5000,10))
summary(fit4)
HPDbayz(fit4$Samples,bound="var") # to get confindence intervals
varcolumns = ( substr(colnames(fit4$Samples),0,3)=="var" )
proportions = t(apply(fit4$Samples[,varcolumns],1,function(x){x/sum(x)}))
HPDbayz(proportions, bound="prob")


fit5 = bayz(gpiCor ~ rn(Clover) + fixf(NS) + fixf(EW) + rn(Rhizobium) + fixf(Inoculation.date) +rn(CloverRhiz), data=d7,chain=c(20000,5000,10))
summary(fit5)
HPDbayz(fit5$Samples,bound="var") # to get confindence intervals
varcolumns = ( substr(colnames(fit5$Samples),0,3)=="var" )
proportions = t(apply(fit5$Samples[,varcolumns],1,function(x){x/sum(x)}))
HPDbayz(proportions, bound="prob")

