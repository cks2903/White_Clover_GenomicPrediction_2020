##############################################################################
#             Checking heritability of growth periods in Clover. 
##############################################################################

# Load in data
{
  setwd("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/New_fixation_Trait_20210106")
  library("lme4")
  #library("BGLR")
  library("BayzR")
  library("parallel")
  d=read.table("greenhouse_data_normalized_24092019.csv",header=T,sep=";",stringsAsFactors = FALSE)
  dim(d)
}

# Apply the same filtering and corrections as in the Genomic Prediction of yield
{
  before=dim(d)[1] # number of observations before filtering
  d=d[-which(d$Rhizobium=="SM73"),]
  d=d[-which(d$Rhizobium=="NO"),]
  d=d[-which(d$n.stolons<4),]
  
  d$Clover[which(d$Clover=="AAran_0104")]="Aaran_0104"
  d$Clover[which(d$Clover=="AAran_0206")]="Aaran_0206"
  
  d$roundRep <- paste(d$Round, d$Replicate, sep='_')
  
  after=dim(d)[1]
  before-after
  print(paste("removed",before-after,"observations,",after, "remains"),sep="")
}


# A function to return a matrix with all individuals and the growth pr. days during a certain period
growthprdaysinintervals<-function(start,end,data){
  matrix = matrix(,nrow = length(unique(data$Barcode)), ncol = 5)
  colnames(matrix)=c("Barcode","Start(days)","End(days)","number_of_days","GPD_in_interval")
  plants=unique(data$Barcode)
  for (plant in plants){
    name=paste(plant,"d",sep="")
    d2=subset(data,Barcode==plant)
    
    period_length=end-start+1
    
    days_matrix=rep(NA,100*2)
    days_matrix <- matrix(days_matrix,nrow=100,ncol=2)
    colnames(days_matrix)=c("day","index")
    for (day in seq(start,end,1)){
      dayidx=which(d2$NormTime==day)
      if (length(dayidx)!=0){
        days_matrix[day,]=c(day,dayidx) # this is so we don't have a problem if the a plant is not measured at the given date
      }
    }
    days_matrix=na.omit(days_matrix)
    
    NAS=period_length-nrow(days_matrix)    
    if (NAS<=period_length*0.25){ #if more than 25% is missing in given growth interval, do not calculate GPD_interval
      start_idx=days_matrix[1,2] 
      end_idx=days_matrix[nrow(days_matrix),2]
      reg_coef=lm(d2$Size[start_idx:end_idx]~(seq(start_idx:end_idx)))$coefficients[2]
      GPD_interval=reg_coef
      context=cbind(plant,start,end,period_length,GPD_interval)
      matrix[which(plants==plant),]=context
    }
  }
  return(matrix)
}

# A function to return a matrix with all individuals and the growth pr. days from a given day to the end of their period (is different for each plant)
OverallgrowthFromDayToEndPeriod<-function(start,data){
  matrix = matrix(,nrow = length(unique(data$Barcode)), ncol = 5)
  colnames(matrix)=c("Barcode","Start(days)","End(days)","number_of_days","GPD_in_interval")
  plants=unique(data$Barcode)
  
  for (plant in plants){
    name=paste(plant,"d",sep="")
    d2=subset(data,Barcode==plant)
    last_measurement = max(d2$NormTime)
    period_length=last_measurement-start+1
    
    days_matrix=rep(NA,100*2)
    days_matrix <- matrix(days_matrix,nrow=100,ncol=2)
    colnames(days_matrix)=c("day","index")
    
    for (day in seq(start,last_measurement,1)){
      dayidx=which(d2$NormTime==day)
      if (length(dayidx)!=0){
        days_matrix[day,]=c(day,dayidx)
      }
    }
    days_matrix=na.omit(days_matrix)
    start_idx=days_matrix[1,2]
    end_idx=days_matrix[nrow(days_matrix),2]
    reg_coef=lm(d2$Size[start_idx:end_idx]~(seq(start_idx:end_idx)))$coefficients[2]
    GPD_interval=reg_coef
    context=cbind(plant,start,last_measurement,period_length,GPD_interval)
    matrix[which(plants==plant),]=context
  }
  return(matrix)
}

# Remove plants that do not grow from day 10 to 20
# These are considered dead or contaminated if they start growing after day 20
first20days=growthprdaysinintervals(10,20,d)
first20days=na.omit(first20days)
nrow(first20days)
hist(first20days[,5],breaks=200) 
mean(first20days[,5])
Dead=(which(first20days[,5]<100)) #These plants die
DeadBarcodes=first20days[Dead,1]
length(DeadBarcodes) #150 plants

# Take a look at the Barcodes being removed
data_frame_of_dead_barcodes= d[which(d$Barcode %in% DeadBarcodes),]
df_list <- split(data_frame_of_dead_barcodes, as.factor(data_frame_of_dead_barcodes$Barcode)) # make a dataframe for each barcode being removed

for (dataframe in df_list){
  name_for_plot = paste("Barcode_",dataframe[1,1],"RemovedBecauseOfNowGrowthDuringDay10and20_",Sys.Date(),".png",sep="")
  png(filename=paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/New_fixation_Trait_20210106",name_for_plot,sep="/"))
  p = plot(dataframe$NormTime,dataframe$Size)
  dev.off()
}

remove=which(d$Barcode %in% DeadBarcodes)
d1=d[-remove,]




# Remove plants that showed an overall negative regression coefficient from 10 dpi to the remaining growth period
Fromday10=OverallgrowthFromDayToEndPeriod(10,d)
Fromday10=na.omit(Fromday10)
nrow(Fromday10)
hist(Fromday10[,5],breaks=200) 
mean(Fromday10[,5])
Weird=(which(Fromday10[,5]<0)) #These plants die
WeirdBarcodes=Fromday10[Weird,1]
length(WeirdBarcodes) #13

# Take a look at the Barcodes being removed
data_frame_of_negative_growth_barcodes= d1[which(d1$Barcode %in% WeirdBarcodes),]
df_list2 <- split(data_frame_of_negative_growth_barcodes, as.factor(data_frame_of_negative_growth_barcodes$Barcode)) # make a dataframe for each barcode being removed

for (dataframe in df_list2){
  name_for_plot = paste("Barcode_",dataframe[1,1],"RemovedBecauseOfNegativeGrowth_",Sys.Date(),".png",sep="")
  png(filename=paste("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/New_fixation_Trait_20210106",name_for_plot,sep="/"))
  p = plot(dataframe$NormTime,dataframe$Size)
  dev.off()
}

remove2=which(d1$Barcode %in% WeirdBarcodes)
d2=d1[-remove2,]

table_of_removed_barcodes = rbind(c(as.vector(DeadBarcodes),as.vector(WeirdBarcodes)))
write.table(table_of_removed_barcodes,paste("Barcodes_removed_based_on_single_Observations_",Sys.Date(),".txt",sep=""),row.names=FALSE,col.names=F,quote=F)


# A function to calculate broad sense heritability for a certain growth period
H2Calc_lme4<-function(dataframe){
  dataframe=as.data.frame(dataframe)
  dataframe =na.omit(dataframe)
  
  # Now connect Barcodes with Clover, Rhizobium, EW, ES, RoundRep and calculate heritability for 1-5 days
  # Initial Size is an average of the size of the 5 first days
  # Correct for initial size
  
  dataframe$Clover=as.character(rep(0,nrow(dataframe)))
  dataframe$Rhizobium=as.character(rep(0,nrow(dataframe)))
  dataframe$EW=as.character(rep(0,nrow(dataframe)))
  dataframe$NS=as.character(rep(0,nrow(dataframe)))
  dataframe$inoculation_date=as.character(rep(0,nrow(dataframe)))
  
  for (i in 1:nrow(dataframe)){
    dataframe$Clover[i]=d2$Clover[which(d2$Barcode==dataframe$Barcode[i])[1]]
    dataframe$Rhizobium[i]=d2$Rhizobium[which(d2$Barcode==dataframe$Barcode[i])[1]]
    dataframe$EW[i]=d2$EW[which(d2$Barcode==dataframe$Barcode[i])[1]]
    dataframe$NS[i]=d2$NS[which(d2$Barcode==dataframe$Barcode[i])[1]]
    dataframe$inoculation_date[i]=d2$inoculation_date[which(d2$Barcode==dataframe$Barcode[i])[1]]
  }
  dataframe=na.omit(dataframe)
  # Check heritabilities for all growth periods
  dataframe$Rhizobium=as.factor(dataframe$Rhizobium)
  dataframe$EW=as.factor(dataframe$EW)
  dataframe$NS=as.factor(dataframe$NS)
  dataframe$inoculation_date=as.factor(dataframe$inoculation_date)
  dataframe$Clover=as.factor(dataframe$Clover)
  
  Correctedforallfixed <- lmer(GPD_in_interval ~ factor(EW) + factor(NS) + Rhizobium + factor(inoculation_date) + (1|Clover), data=dataframe)
  summary(Correctedforallfixed)
  re_dat = as.data.frame(VarCorr(Correctedforallfixed))
  VarClover=re_dat[1,'vcov']
  VarResidual=re_dat[2,'vcov']
  H2=VarClover/(VarClover+VarResidual) 
  return(H2)
}

###############################################################################################
#
#           GPD in 5 days overlapping intervals and clover heritability calculation  
#
###############################################################################################

# Calculate GPD in 5 days overlapping intervals 

for (i in 1:52){
  start=i
  end=i+4
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d2))
}


# Heritability calculation; applying the function to different growth intervals

matrix1=rep(NA,52*3)
matrix1 <- matrix(matrix1,nrow=52,ncol=3)
colnames(matrix1)=c("start","end","H2")

for (i in 1:52){
  start=i
  end=i+4
  name=paste("day",i,"to",i+4,sep="")
  get(name)
  BroadSenseHer_lmer4=H2Calc_lme4(get(name))
  matrix1[i,]=cbind(c(start,end,BroadSenseHer_lmer4))
  
}

matrix1
matrix1=na.omit(matrix1)
plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1))
cor(matrix1[,1],matrix1[,3])

plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))


write.table(matrix1,paste("H2_5days_intervals",Sys.Date(),".csv",sep="_"),quote=F,row.names = F)


###############################################################################################
#
#           GPD in 10 days overlapping intervals and clover heritability calculation
#
###############################################################################################

# Calculate GPD in 10 days overlapping intervals 

for (i in 1:52){
  start=i
  end=i+9
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d2))
}


# Heritability calculation; applying the function to different growth intervals

matrix2=rep(NA,52*3)
matrix2 <- matrix(matrix2,nrow=52,ncol=3)
colnames(matrix2)=c("start","end","H2")

for (i in 1:52){
  start=i
  end=i+9
  name=paste("day",i,"to",i+9,sep="")
  get(name)
  Her=H2Calc_lme4(get(name))
  matrix2[i,]=cbind(c(start,end,Her))
  
}

matrix2
matrix2=na.omit(matrix2)
plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1))
cor(matrix2[,1],matrix2[,3])

plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

write.table(matrix2,paste("H2_10days_intervals",Sys.Date(),".csv",sep="_"),quote=F,row.names = F)

##############################################################################################
#
#           GPD in 15 days overlapping intervals and clover heritability calculation
#
###############################################################################################

# Calculate GPD in 15 days overlapping intervals 

for (i in 1:52){
  start=i
  end=i+14
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d2))
}


# Heritability calculation; applying the function to different growth intervals

matrix3=rep(NA,52*3)
matrix3 <- matrix(matrix3,nrow=52,ncol=3)
colnames(matrix3)=c("start","end","H2")

for (i in 1:52){
  start=i
  end=i+14
  name=paste("day",i,"to",i+14,sep="")
  get(name)
  Her=H2Calc_lme4(get(name))
  matrix3[i,]=cbind(c(start,end,Her))
  
}

matrix3
matrix3=na.omit(matrix3)
plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1))
cor(matrix3[,1],matrix3[,3])

plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))


write.table(matrix3,paste("H2_15days_intervals",Sys.Date(),".csv",sep="_"),quote=F,row.names = F)

##############################################################################################
#
#           GPD in 20 days overlapping intervals and clover heritability calculation
#
###############################################################################################


# Calculate GPD in 20 days overlapping intervals 

for (i in 1:52){
  start=i
  end=i+19
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d2))
}


# Heritability calculation; applying the function to different growth intervals

matrix4=rep(NA,52*3)
matrix4 <- matrix(matrix4,nrow=52,ncol=3)
colnames(matrix4)=c("start","end","H2")

for (i in 1:52){
  start=i
  end=i+19
  name=paste("day",i,"to",i+19,sep="")
  get(name)
  Her=H2Calc_lme4(get(name))
  matrix4[i,]=cbind(c(start,end,Her))
  
}

matrix4
matrix4=na.omit(matrix4)
plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1))
cor(matrix4[,1],matrix4[,3])

plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))


write.table(matrix4,paste("H2_20days_intervals",Sys.Date(),".csv",sep="_"),quote=F,row.names = F)

################################
# Get highly her phenotype     #
#day 11 to 25, median around day 18
head(day11to25)
phenotype=merge(d2,day11to25,by="Barcode")
length(unique(phenotype$Barcode))

phenotype1=phenotype[,c(1,5,7:12,15:17,22,26)]
nrow(phenotype1)

phenotype2=phenotype1[!duplicated(phenotype1$Barcode), ]       
nrow(phenotype2)

#write.table(phenotype2,"GPD_day11to25.csv",quote=F,row.names = F)


# Plot with all intervals
library(ggplot2)
library(reshape2)

matrix1_=data.frame(matrix1)
matrix2_=data.frame(matrix2)
matrix3_=data.frame(matrix3)
matrix4_=data.frame(matrix4)

matrix_all=merge(matrix1_,matrix2_,by= "start")
matrix_all=merge(matrix_all,matrix3_,by="start")
matrix_all=merge(matrix_all,matrix4_,by="start")
colnames(matrix_all)[3]="H2_5days"
colnames(matrix_all)[5]="H2_10days"
colnames(matrix_all)[7]="H2_15days"
colnames(matrix_all)[9]="H2_20days"
matrix_all[,2]=NULL
matrix_all[,3]=NULL
matrix_all[,4]=NULL
matrix_all[,5]=NULL
head(matrix_all)

premelted=matrix_all[1:30,]
melted=melt(premelted,id.var="start")


ggplot(data=melted,aes(x=start, y=value,group=variable,fill=variable)) +
  geom_line(aes(color=variable))+
  geom_point(aes(color=variable)) +
  scale_color_manual(values=c("#984EA3","#1B9E77","#386CB0","#FB8072")) +
  xlab("Start of interval") + ylab("Clover heritability") +
  #geom_hline(yintercept=0.45, linetype="dashed", color = "black", size=0.3) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Clover heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold"))


# This plot is fine but I want to plot the median of an interval instead
# this will make it easier to compare the maximum of each curve, because they actually have maximums around the same time




# Calculate median for each period
matrix1_=data.frame(matrix1)
matrix1_=matrix1[1:30,]
matrix2_=data.frame(matrix2)
matrix2_=matrix2[1:30,]
matrix3_=data.frame(matrix3)
matrix3_=matrix3[1:30,]
matrix4_=data.frame(matrix4)
matrix4_=matrix4[1:30,]


matrix_all=merge(matrix1_,matrix2_,by= "start")
matrix_all=merge(matrix_all,matrix3_,by="start")
matrix_all=merge(matrix_all,matrix4_,by="start")
colnames(matrix_all)[3]="H2_5days"
colnames(matrix_all)[5]="H2_10days"
colnames(matrix_all)[7]="H2_15days"
colnames(matrix_all)[9]="H2_20days"
colnames(matrix_all)[2]="end.5"
colnames(matrix_all)[4]="end.10"
colnames(matrix_all)[6]="end.15"
colnames(matrix_all)[8]="end.20"
median.5=apply(matrix_all[,c("start","end.5")],1,median)
median.10=apply(matrix_all[,c("start","end.10")],1,median)
median.15=apply(matrix_all[,c("start","end.15")],1,median)
median.20=apply(matrix_all[,c("start","end.20")],1,median)


newmatrix5days=cbind(median.5,matrix_all$H2_5days)
colnames(newmatrix5days)=c("median","H2_5days")
newmatrix10days=cbind(median.10,matrix_all$H2_10days)
colnames(newmatrix10days)=c("median","H2_10days")
newmatrix15days=cbind(median.15,matrix_all$H2_15days)
colnames(newmatrix15days)=c("median","H2_15days")
newmatrix20days=cbind(median.20,matrix_all$H2_20days)
colnames(newmatrix20days)=c("median","H2_20days")


matrix_all_=rbind(newmatrix5days,newmatrix10days,newmatrix15days,newmatrix20days)
matrix_all_=as.data.frame(matrix_all_)
matrix_all_$group[1:30]="5 days intervals"
matrix_all_$group[31:60]="10 days intervals"
matrix_all_$group[61:90]="15 days intervals"
matrix_all_$group[91:nrow(matrix_all_)]="20 days intervals"

matrix_all_=matrix_all_[order(matrix_all_[,1]),] 

#you can remove all median values above 30
matrix_all_=matrix_all_[-which(matrix_all_$median>30),]


ggplot(data=matrix_all_,aes(x=median, y=H2_5days,group=group,fill=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group)) +
  scale_color_manual(values=c("#1B9E77","#386CB0","#FB8072","#984EA3")) +
  xlab("Middle of interval") + ylab("Clover heritability") +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Clover heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold"))

# save as 7 x 10 pdf

# Now try to add error bars to the point
matrix_all_$xmin=NA
matrix_all_$xmax=NA



growthinterval1=which(matrix_all_$group=="5 days intervals")
matrix_all_$xmin[growthinterval1]=matrix_all_$median[growthinterval1]-2
matrix_all_$xmax[growthinterval1]=matrix_all_$median[growthinterval1]+2

growthinterval2=which(matrix_all_$group=="10 days intervals")
matrix_all_$xmin[growthinterval2]=matrix_all_$median[growthinterval2]-4.5
matrix_all_$xmax[growthinterval2]=matrix_all_$median[growthinterval2]+4.5

growthinterval3=which(matrix_all_$group=="15 days intervals")
matrix_all_$xmin[growthinterval3]=matrix_all_$median[growthinterval3]-7
matrix_all_$xmax[growthinterval3]=matrix_all_$median[growthinterval3]+7

growthinterval4=which(matrix_all_$group=="20 days intervals")
matrix_all_$xmin[growthinterval4]=matrix_all_$median[growthinterval4]-9.5
matrix_all_$xmax[growthinterval4]=matrix_all_$median[growthinterval4]+9.5


ggplot(data=matrix_all_,aes(x=median, y=H2_5days,group=group,fill=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group)) +
  geom_errorbarh(data=matrix_all_,aes(xmin = xmin,xmax = xmax),color="gray",alpha=0.4,height=0, size=0.5) +
  scale_color_manual(values=c("#1B9E77","#386CB0","#FB8072","#984EA3")) +
  xlab("Middle of interval") + ylab("Clover heritability") +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Clover heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold")) 


ggplot(data=matrix_all_,aes(x=median, y=H2_5days,group=group,fill=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group)) +
  geom_errorbarh(data=matrix_all_,aes(color=group,alpha=0.01,xmin = xmin,xmax = xmax),height=0, size=0.3) +
  scale_color_manual(values=c("#1B9E77","#386CB0","#FB8072","#984EA3")) +
  xlab("Middle of interval") + ylab("Clover heritability") +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Clover heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold")) 













########################################################################################################################################
#                                                                                                                                      #
#                   Now calculate rhizobium heritability instead of clover                                                             #
#                                                                                                                                      #
########################################################################################################################################

#Filtering as in Clover heritability calculations
H2Calc_Rhiz_lme4<-function(dataframe){
  dataframe=as.data.frame(dataframe)
  dataframe =na.omit(dataframe)
    
  # Now connect Barcodes with Clover, Rhizobium, EW, ES, RoundRep and calculate heritability for 1-5 days
  # Initial Size is an average of the size of the 5 first days
  # Correct for initial size
    
  dataframe$Clover=as.character(rep(0,nrow(dataframe)))
  dataframe$Rhizobium=as.character(rep(0,nrow(dataframe)))
  dataframe$EW=as.character(rep(0,nrow(dataframe)))
  dataframe$NS=as.character(rep(0,nrow(dataframe)))
  dataframe$inoculation_date=as.character(rep(0,nrow(dataframe)))
    
  for (i in 1:nrow(dataframe)){
    dataframe$Clover[i]=d2$Clover[which(d2$Barcode==dataframe$Barcode[i])[1]]
    dataframe$Rhizobium[i]=d2$Rhizobium[which(d2$Barcode==dataframe$Barcode[i])[1]]
    dataframe$EW[i]=d2$EW[which(d2$Barcode==dataframe$Barcode[i])[1]]
    dataframe$NS[i]=d2$NS[which(d2$Barcode==dataframe$Barcode[i])[1]]
    dataframe$inoculation_date[i]=d2$inoculation_date[which(d2$Barcode==dataframe$Barcode[i])[1]]
    }
  dataframe=na.omit(dataframe)
  
  # Check heritabilities for all growth periods
  dataframe$Rhizobium=as.factor(dataframe$Rhizobium)
  dataframe$EW=as.factor(dataframe$EW)
  dataframe$NS=as.factor(dataframe$NS)
  dataframe$inoculation_date=as.factor(dataframe$inoculation_date)
  dataframe$Clover=as.factor(dataframe$Clover)
    
  Correctedforallfixed <- lmer(GPD_in_interval ~ factor(EW) + factor(NS) + Clover + factor(inoculation_date) + (1|Rhizobium), data=dataframe)
  summary(Correctedforallfixed)
  re_dat = as.data.frame(VarCorr(Correctedforallfixed))
  VarRhiz=re_dat[1,'vcov']
  VarResidual=re_dat[2,'vcov']
  H2=VarRhiz/(VarRhiz+VarResidual) 
  return(H2)
  }
  
  
# Heritability calculation; applying the function to different growth intervals
  
  matrix1_rhiz=rep(NA,62*3)
  matrix1_rhiz <- matrix(matrix1_rhiz,nrow=62,ncol=3)
  colnames(matrix1_rhiz)=c("start","end","H2")
  
  for (i in 1:62){
    start=i
    end=i+4
    name=paste("day",i,"to",i+4,sep="")
    get(name)
    Her=H2Calc_Rhiz_lme4(get(name))
    matrix1_rhiz[i,]=cbind(c(start,end,Her))
    
  }
  
  matrix1_rhiz
  matrix1_rhiz=na.omit(matrix1_rhiz)
  plot(matrix1_rhiz[,1],matrix1_rhiz[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1))
  cor(matrix1_rhiz[,1],matrix1_rhiz[,3])
  
  plot(matrix1_rhiz[,1],matrix1_rhiz[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))
  
  write.table(matrix1_rhiz,paste("Rhiz_H2_5days_intervals",Sys.Date(),".csv"),quote=F,row.names = F)
  
  
  # Heritability calculation; applying the function to different growth intervals
  
  matrix2_rhiz=rep(NA,62*3)
  matrix2_rhiz <- matrix(matrix2_rhiz,nrow=62,ncol=3)
  colnames(matrix2_rhiz)=c("start","end","H2")
  
  for (i in 1:62){
    start=i
    end=i+9
    name=paste("day",i,"to",i+9,sep="")
    get(name)
    Her=H2Calc_Rhiz_lme4(get(name))
    matrix2_rhiz[i,]=cbind(c(start,end,Her))
    
  }
  
  matrix2_rhiz
  matrix2_rhiz=na.omit(matrix2_rhiz)
  plot(matrix2_rhiz[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1))
  cor(matrix2_rhiz[,1],matrix2[,3])
  
  plot(matrix2_rhiz[,1],matrix2_rhiz[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))
  
  write.table(matrix2_rhiz,paste("Rhiz_H2_10days_intervals",Sys.Date(),".csv"),quote=F,row.names = F)
  
  
  
  # Heritability calculation; applying the function to different growth intervals
  
  matrix3_rhiz=rep(NA,62*3)
  matrix3_rhiz <- matrix(matrix3_rhiz,nrow=62,ncol=3)
  colnames(matrix3_rhiz)=c("start","end","H2")
  
  for (i in 1:62){
    start=i
    end=i+14
    name=paste("day",i,"to",i+14,sep="")
    get(name)
    Her=H2Calc_Rhiz_lme4(get(name))
    matrix3_rhiz[i,]=cbind(c(start,end,Her))
    
  }
  
  matrix3_rhiz
  matrix3_rhiz=na.omit(matrix3_rhiz)
  plot(matrix3_rhiz[,1],matrix3_rhiz[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1))
  cor(matrix3_rhiz[,1],matrix3_rhiz[,3])
  
  plot(matrix3_rhiz[,1],matrix3_rhiz[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))
  
  write.table(matrix3_rhiz,paste("Rhiz_H2_15days_intervals",Sys.Date(),".csv"),quote=F,row.names = F)
  
  
  # Heritability calculation; applying the function to different growth intervals
  
  matrix4_rhiz=rep(NA,62*3)
  matrix4_rhiz <- matrix(matrix4_rhiz,nrow=62,ncol=3)
  colnames(matrix4_rhiz)=c("start","end","H2")
  
  for (i in 1:62){
    start=i
    end=i+19
    name=paste("day",i,"to",i+19,sep="")
    get(name)
    Her=H2Calc_Rhiz_lme4(get(name))
    matrix4_rhiz[i,]=cbind(c(start,end,Her))
    
  }
  
  matrix4_rhiz
  matrix4_rhiz=na.omit(matrix4_rhiz)
  plot(matrix4_rhiz[,1],matrix4_rhiz[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1))
  cor(matrix4_rhiz[,1],matrix4_rhiz[,3])
  
  plot(matrix4_rhiz[,1],matrix4_rhiz[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))
  write.table(matrix4_rhiz,paste("Rhiz_H2_20days_intervals",Sys.Date(),".csv"),quote=F,row.names = F)
  

  # Plot with all intervals
  library(ggplot2)
  library(reshape2)
  
  matrix1_rhiz_=data.frame(matrix1_rhiz)
  matrix2_rhiz_=data.frame(matrix2_rhiz)
  matrix3_rhiz_=data.frame(matrix3_rhiz)
  matrix4_rhiz_=data.frame(matrix4_rhiz)
  
  matrix_all=merge(matrix1_rhiz_,matrix2_rhiz_,by= "start")
  matrix_all=merge(matrix_all,matrix3_rhiz_,by="start")
  matrix_all=merge(matrix_all,matrix4_rhiz_,by="start")
  colnames(matrix_all)[3]="H2_5days"
  colnames(matrix_all)[5]="H2_10days"
  colnames(matrix_all)[7]="H2_15days"
  colnames(matrix_all)[9]="H2_20days"
  matrix_all[,2]=NULL
  matrix_all[,3]=NULL
  matrix_all[,4]=NULL
  matrix_all[,5]=NULL
  head(matrix_all)
  
  premelted=matrix_all[1:30,]
  melted=melt(premelted,id.var="start")
  
  
  ggplot(data=melted,aes(x=start, y=value,group=variable,fill=variable)) +
    geom_line(aes(color=variable))+
    geom_point(aes(color=variable)) +
    scale_color_manual(values=c("#984EA3","#1B9E77","#386CB0","#FB8072")) +
    xlab("Start of interval") + ylab("Rhizobium heritability") +
    geom_hline(yintercept=0.20, linetype="dashed", color = "white", size=0.3,alpha=0) +
    scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
    scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
    ggtitle("Rhizobium heritability for different growth intervals") +
    theme( panel.background = element_rect(fill = NA),
           panel.grid.major = element_line(colour = "grey50",size=0.05),
           plot.title=element_text(color="black", size=16, face="bold"))
  
  # Calculate median for each period
  matrix1_=data.frame(matrix1_rhiz)
  matrix1_=matrix1_[1:30,]
  matrix2_=data.frame(matrix2_rhiz)
  matrix2_=matrix2_[1:30,]
  matrix3_=data.frame(matrix3_rhiz)
  matrix3_=matrix3_[1:30,]
  matrix4_=data.frame(matrix4_rhiz)
  matrix4_=matrix4_[1:30,]
  
  
  matrix_all=merge(matrix1_,matrix2_,by= "start")
  matrix_all=merge(matrix_all,matrix3_,by="start")
  matrix_all=merge(matrix_all,matrix4_,by="start")
  colnames(matrix_all)[3]="H2_5days"
  colnames(matrix_all)[5]="H2_10days"
  colnames(matrix_all)[7]="H2_15days"
  colnames(matrix_all)[9]="H2_20days"
  colnames(matrix_all)[2]="end.5"
  colnames(matrix_all)[4]="end.10"
  colnames(matrix_all)[6]="end.15"
  colnames(matrix_all)[8]="end.20"
  median.5=apply(matrix_all[,c("start","end.5")],1,median)
  median.10=apply(matrix_all[,c("start","end.10")],1,median)
  median.15=apply(matrix_all[,c("start","end.15")],1,median)
  median.20=apply(matrix_all[,c("start","end.20")],1,median)
  
  
  newmatrix5days=cbind(median.5,matrix_all$H2_5days)
  colnames(newmatrix5days)=c("median","H2_5days")
  newmatrix10days=cbind(median.10,matrix_all$H2_10days)
  colnames(newmatrix10days)=c("median","H2_10days")
  newmatrix15days=cbind(median.15,matrix_all$H2_15days)
  colnames(newmatrix15days)=c("median","H2_15days")
  newmatrix20days=cbind(median.20,matrix_all$H2_20days)
  colnames(newmatrix20days)=c("median","H2_20days")
  
  
  matrix_all_=rbind(newmatrix5days,newmatrix10days,newmatrix15days,newmatrix20days)
  matrix_all_=as.data.frame(matrix_all_)
  matrix_all_$group[1:30]="5 days intervals"
  matrix_all_$group[31:60]="10 days intervals"
  matrix_all_$group[61:90]="15 days intervals"
  matrix_all_$group[91:nrow(matrix_all_)]="20 days intervals"
  
  matrix_all_=matrix_all_[order(matrix_all_[,1]),] 
  
  #you can remove all median values above 30
  matrix_all_=matrix_all_[-which(matrix_all_$median>30),]
  
  
  ggplot(data=matrix_all_,aes(x=median, y=H2_5days,group=group,fill=group)) +
    geom_line(aes(color=group))+
    geom_point(aes(color=group)) +
    scale_color_manual(values=c("#1B9E77","#386CB0","#FB8072","#984EA3")) +
    xlab("Middle of interval") + ylab("Rhizobium heritability") +
    geom_hline(yintercept=0.20, linetype="dashed", color = "white", size=0.3,alpha=0) +
    scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
    scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
    ggtitle("Rhizobium heritability for different growth intervals") +
    theme( panel.background = element_rect(fill = NA),
           panel.grid.major = element_line(colour = "grey50",size=0.05),
           plot.title=element_text(color="black", size=16, face="bold"))
  
  
  # Now try to add error bars to the point
  matrix_all_$xmin=NA
  matrix_all_$xmax=NA
  
  
  
  growthinterval1=which(matrix_all_$group=="5 days intervals")
  matrix_all_$xmin[growthinterval1]=matrix_all_$median[growthinterval1]-2
  matrix_all_$xmax[growthinterval1]=matrix_all_$median[growthinterval1]+2
  
  growthinterval2=which(matrix_all_$group=="10 days intervals")
  matrix_all_$xmin[growthinterval2]=matrix_all_$median[growthinterval2]-4.5
  matrix_all_$xmax[growthinterval2]=matrix_all_$median[growthinterval2]+4.5
  
  growthinterval3=which(matrix_all_$group=="15 days intervals")
  matrix_all_$xmin[growthinterval3]=matrix_all_$median[growthinterval3]-7
  matrix_all_$xmax[growthinterval3]=matrix_all_$median[growthinterval3]+7
  
  growthinterval4=which(matrix_all_$group=="20 days intervals")
  matrix_all_$xmin[growthinterval4]=matrix_all_$median[growthinterval4]-9.5
  matrix_all_$xmax[growthinterval4]=matrix_all_$median[growthinterval4]+9.5
  
  
  ggplot(data=matrix_all_,aes(x=median, y=H2_5days,group=group,fill=group)) +
    geom_line(aes(color=group))+
    geom_point(aes(color=group)) +
    geom_errorbarh(data=matrix_all_,aes(xmin = xmin,xmax = xmax),color="gray",alpha=0.4,height=0, size=0.5) +
    scale_color_manual(values=c("#1B9E77","#386CB0","#FB8072","#984EA3")) +
    xlab("Middle of interval") + ylab("Rhizobium heritability") +
    geom_hline(yintercept=0.20, linetype="dashed", color = "white", size=0.3,alpha=0) +
    scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
    scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
    ggtitle("Rhizobium heritability for different growth intervals") +
    theme( panel.background = element_rect(fill = NA),
           panel.grid.major = element_line(colour = "grey50",size=0.05),
           plot.title=element_text(color="black", size=16, face="bold")) 
  
  
