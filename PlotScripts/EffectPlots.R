#################################################
#################################################
#     Effect plots of top 25 SNPs in GWAS
#
#################################################
#################################################

# load libraries need
library(ggplot2)

# set wd
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/GWAS_all_data_averages_20210215")

# load data
data=read.table("_pid1_iSizeMean_emmax_none_t75.pvals",head=T,sep=",")
data=as.data.frame(data)
head(data)
data_sorted= data[order(data$scores),] 
head(data_sorted)

data1=read.table("_pid2_gpd_emmax_none_t75.pvals",head=T,sep=",")
data1=as.data.frame(data1)
head(data1)
data_sorted1= data1[order(data1$scores),] 
head(data_sorted1)

genotypes= read.table("Geno_LDfiltered0.5and1.00_20200728_noextracolumns.csv",head=T,sep=",")
genotypes=as.data.frame(genotypes)
head(genotypes)
colnames(genotypes)[1:2]=c("chromosomes","positions")

phenotypes =read.table("Phenotypes_20210215.csv",head=T,sep=",")
phenotypes=as.data.frame(phenotypes)


# Function doing effect plot for a given number of top SNPs
Effectplots <- function(data,m,trait){
  
  dataofinterest = data[1:m,]
  dataofinterest$identifier=paste(dataofinterest$chromosomes,dataofinterest$positions,sep="-")
  
  
  for (i in seq(1:m)){
    markerofinterest = dataofinterest$identifier[i]
    genotypes$identifier=paste(genotypes$chromosomes,genotypes$positions,sep="-")
    genotyperow = which(genotypes$identifier==markerofinterest)
    genotypesforeffectplot = genotypes[genotyperow,c(3:164)]
  
    ind_keep= which(colnames(genotypesforeffectplot) %in% phenotypes$Accession)
    genotypesforeffectplot_filt= genotypesforeffectplot[,ind_keep]
    dim(genotypesforeffectplot_filt)
  
    genotypesforeffectplot_filt_t =t(genotypesforeffectplot_filt)
    colnames(genotypesforeffectplot_filt_t)="Genotype"
    genotypesforeffectplot_filt_t=as.data.frame(genotypesforeffectplot_filt_t)
    genotypesforeffectplot_filt_t$Accession = rownames(genotypesforeffectplot_filt_t)
  
    merged = merge(genotypesforeffectplot_filt_t,phenotypes,by="Accession")
    merged$Genotype=as.factor(merged$Genotype)
  
    label1 = strsplit(markerofinterest,"-")
    Chromosome = label1[[1]][1]
    Position = label1[[1]][2]
    label2=paste("Chr",Chromosome,"_Pos",Position,sep="")
    
    if (trait=="iSize"){
      y = "iSizeMean"
      
      }else{
        y ="gpd"
      }
  
    p = ggplot(merged, aes(x = Genotype, y = merged[,which(colnames(merged)==y)])) + geom_boxplot(fill="grey",alpha=0.5)  + geom_jitter(position=position_jitter(0.2)) +
      ylab(trait) +
      ggtitle(label2) +
      theme_classic()
    p
  
    ggsave(paste('EffectPlot',label2,trait,Sys.Date(),'.png',sep="_"), plot = p, width = 15, height = 10, unit = 'cm')
    
  }
}


Effectplots(data_sorted,25,"iSize")

Effectplots(data_sorted,25,"gpd")
