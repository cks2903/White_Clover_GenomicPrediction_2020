
library(ggplot2)
library(multcompView)
library(forcats)


data=read.table("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/F1_prediction/F1Phenotypes_AllData.txt",sep="\t",head=T)
head(data)

#png(file="F1distribution.png",width=20,height=10,units="cm",res=300)
p=ggplot(data,aes(x=reorder(Cross,X20200226_dryweight,FUN=median), y=X20200226_dryweight)) + 
  geom_boxplot(fill="steelblue",alpha=0.7)+
  labs(title="Dry weight of F1 populations",x="F1 Population", y = "Dry weight")+
  theme_classic() +theme(axis.text.x = element_text(size = 10, angle = 90,hjust=1),axis.text.y = element_text(size = 10),axis.title=element_text(size=15))
p + geom_jitter(shape=16, position=position_jitter(0.2))


fit <- aov(X20200226_dryweight~ factor(Cross)  , data=data)

results <- TukeyHSD(fit, ordered=TRUE)
multcompLetters4(fit, results)
#dev.off()


dataNew=data[-which(data$Cross=="DLF2" | data$Cross=="DLF3" | data$Cross=="DLF4" | data$Cross=="DLF5" | data$Cross=="YT"),]
length(unique(dataNew$Cross))==9

#png(file="F1distribution_9pop.png",width=20,height=12,units="cm",res=300)
p=ggplot(dataNew,aes(x=reorder(Cross,X20200226_dryweight,FUN=median), y=X20200226_dryweight)) + 
  geom_boxplot(fill="#4A9667",alpha=0.9)+
  labs(title="Dry weight of F1 populations",x="F1 Population", y = "Dry weight")+
  theme_classic() +theme(axis.text.x = element_text(size = 10, angle = 90,hjust=1),axis.text.y = element_text(size = 10),axis.title=element_text(size=15))
p + geom_jitter(shape=16, position=position_jitter(0.2), size =5)


#png(file="F1distribution_9pop_ColoredByPop_dark.png",width=60,height=20,units="cm",res=300)
pops=c("DLF1","LJ L6","LJ L7","LJ L9","LJ H1","LJ H2","LJ H3","LJ H5","SUA ")
colors=c("#242c5e","#b01f35","#f6c29d","#f48988","#443659","#580b1d","#98d3ce","#71a188","#426468")
colormatrix=cbind(pops,colors)
colormatrix=as.data.frame(colormatrix)
colnames(colormatrix)=c("Cross","col")

dataNew_colMerged=merge(dataNew,colormatrix,by="Cross")
nrow(dataNew)==nrow(dataNew_colMerged)

colors=c("#242c5e","#443659","#580b1d","#98d3ce","#71a188","#b01f35","#f6c29d","#f48988","#426468")
p=ggplot(dataNew_colMerged,aes(x=reorder(Cross,X20200226_dryweight,FUN=median), y=X20200226_dryweight,fill=Cross)) + 
  geom_boxplot()+
  scale_fill_manual(values=colors)+
  labs(title="Dry weight of F1 populations",x="F1 Population", y = "Dry weight")+
  theme_classic() +theme(axis.text.x = element_text(size = 10, angle = 90,hjust=1),axis.text.y = element_text(size = 10),axis.title=element_text(size=15))
p + geom_jitter(shape=16, position=position_jitter(0.2), size = 5, alpha=0.9)


#dev.off()

fit <- aov(X20200226_dryweight~ factor(Cross)  , data=dataNew)
results <- TukeyHSD(fit, ordered=TRUE)
multcompLetters4(fit, results)
ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/F1distribution_9pop.pdf", width =60, height = 20,units = "cm",useDingbats=FALSE)



p=ggplot(dataNew_colMerged,aes(x=reorder(Cross,X20200226_dryweight,FUN=median), y=X20200226_dryweight,fill=Cross)) + 
  geom_boxplot()+
  scale_fill_manual(values=colors)+
  labs(title="Dry weight of F1 populations",x="F1 Population", y = "Dry weight")+
  theme_classic() +theme(axis.text.x = element_text(size = 10, angle = 90,hjust=1),axis.text.y = element_text(size = 10),axis.title=element_text(size=15))
p + geom_jitter(shape=16, position=position_jitter(0.2), size = 5, alpha=0.3)
ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/F1distribution_9pop_loweralpha.pdf", width =60, height = 20,units = "cm",useDingbats=FALSE)
