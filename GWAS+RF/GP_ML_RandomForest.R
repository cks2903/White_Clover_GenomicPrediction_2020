#########################################################################
#                 Using machine learning for genomic prediction         #
#                                                                       #
#########################################################################

library(caret)
args=commandArgs(trailingOnly = TRUE)

##An example comes with the mlbench package
{
library(mlbench)
#data("iris") #5 columns sepal length, sepal width. petal length, petal widt and species (classification)
#iris$Species = as.factor(iris$Species)
#To divide data into training and testing (25% of the data)
#inTraining <- createDataPartition(iris$Species, p = .75, list = FALSE)
#training <- iris[ inTraining,]
#testing  <- iris[-inTraining,]


#Create the Cross validation method 
#FitControl=trainControl(method="repeatedcv",number=10,repeats=10)		# A 10-fold CV

# Train a model
#set.seed(12345)


#model=train(training[,2:ncol(training)],as.factor(training[,1]),method= "ranger",importance=permutation) #fit training
#model

#gbmFit1 <- train(Species ~ ., data = training, 
                 #method = "gbm", 
                 #trControl = FitControl,
                 #verbose = FALSE)

#gbmFit1

#Accuracy is used by Caret to chose the final model. It is the overall agreement rate btw. the CV methods
#Kappat is another statistical method used for accessing models with categorical variables 
#CARET chose the first model with an interaction depth of 1, number of trees at 50, an accuracy of 97% and a Kappa of 95%.
#The models very in shrinkage 


#Finally we can use the training model to predict classifications and probabilities for the test data set
#predictions<-predict(object=model,testing[,2:ncol(testing)])
#predictions
}

####################### Random forest model on your data######################################

#Load phenotype data
{
pheno=read.table(args[1],sep=",",header=T)
print(head(pheno))
iSize=pheno[,2]
print(iSize)
}

# Divide into predefined training populations, the ones used for GWAS as well
{
f=read.table(args[2],fill = TRUE)

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

testpop1_idx=which(pheno$Accession %in% testpop1)
testpop2_idx=which(pheno$Accession %in% testpop2)
testpop3_idx=which(pheno$Accession %in% testpop3)
testpop4_idx=which(pheno$Accession %in% testpop4)
testpop5_idx=which(pheno$Accession %in% testpop5)
testpop6_idx=which(pheno$Accession %in% testpop6)

tests=list(testpop1_idx,testpop2_idx,testpop3_idx,testpop4_idx,testpop5_idx,testpop6_idx)
}


# Predict group 1
{
# Load genotypes for group 1 
{
  geno=read.table(args[3],sep=",",header=T)  
  dim(geno)
  print(head(geno))

  geno_shorter=geno[,3:(ncol(geno)-2)]
  geno_shorter_t=t(geno_shorter)
  rownames(geno_shorter_t)=colnames(geno_shorter)
}

#Combine genotypes and phenotypes into one large dataframe
{  
  cloverdata=cbind(as.numeric(as.character(iSize)),geno_shorter_t)
  colnames(cloverdata)=c("iSize",as.character(geno$identifier))
}

# Machine learning using the first group as testing and group2,3,4,5,6 as training
training <- cloverdata[-testpop1_idx,]
training_df=as.data.frame(training)
testing <- cloverdata[testpop1_idx,]
testing_df <- as.data.frame(testing)
print(head(training_df))
print(dim(training_df))


model=train(training_df[,2:ncol(training_df)],training_df$iSize,method= "ranger",importance = 'permutation') #fit training
model
varImp(model)

predictions_group1 =predict(object=model, testing[,2:ncol(testing)])
predictions_group1


Predicted_obs1=cbind(as.character(pheno[testpop1_idx,1]),predictions_group1,pheno[testpop1_idx,2])
colnames(Predicted_obs1)=c("Individual","Prediction","Observed")

VaiableImportance_obs1=varImp(model)
VaiableImportance_obs1_df=VaiableImportance_obs1$importance
VaiableImportance_obs1_df$SNP=rownames(VaiableImportance_obs1_df)
colnames(VaiableImportance_obs1_df)=c("SNP","Overall importance")
}





# Predict group 2
{
# Load genotypes for group 2
{
  geno=read.table(args[4],sep=",",header=T)  
  dim(geno)
  geno_shorter=geno[,3:(ncol(geno)-2)]
  geno_shorter_t=t(geno_shorter)
  rownames(geno_shorter_t)=colnames(geno_shorter)
}

#Combine genotypes and phenotypes into one large dataframe
{  
  if (all(rownames(geno_shorter_t)==pheno[,1])){
    cloverdata=cbind(as.numeric(as.character(iSize)),geno_shorter_t)
    colnames(cloverdata)=c("iSize",as.character(geno$identifier))
  }
}
# Machine learning using the second group as testing and group1,3,4,5,6 as training
training <- cloverdata[-testpop2_idx,]
training_df=as.data.frame(training)
testing <- cloverdata[testpop2_idx,]
testing_df <- as.data.frame(testing)


model=train(training_df[,2:ncol(training)],training_df$iSize,method= "ranger", importance = 'permutation') #fit training
model

predictions_group2 =predict(object=model, testing[,2:ncol(testing)])
predictions_group2

Predicted_obs2=cbind(as.character(pheno[testpop2_idx,1]),predictions_group2,pheno[testpop2_idx,2])
colnames(Predicted_obs2)=c("Individual","Prediction","Observed")

VaiableImportance_obs2=varImp(model)
VaiableImportance_obs2_df=VaiableImportance_obs2$importance
VaiableImportance_obs2_df$SNP=rownames(VaiableImportance_obs2_df)
colnames(VaiableImportance_obs2_df)=c("SNP","Overall importance")

}



# Predict group 3
{
# Load genotypes for group 3
{
  geno=read.table(args[5],sep=",",header=T)  
  dim(geno)
  geno_shorter=geno[,3:(ncol(geno)-2)]
  geno_shorter_t=t(geno_shorter)
  rownames(geno_shorter_t)=colnames(geno_shorter)
}

#Combine genotypes and phenotypes into one large dataframe
{  
  if (all(rownames(geno_shorter_t)==pheno[,1])){
    cloverdata=cbind(as.numeric(as.character(iSize)),geno_shorter_t)
    colnames(cloverdata)=c("iSize",as.character(geno$identifier))
  }
}
# Machine learning using the third group as testing and group1,2,4,5,6 as training
training <- cloverdata[-testpop3_idx,]
training_df=as.data.frame(training)
testing <- cloverdata[testpop3_idx,]
testing_df <- as.data.frame(testing)


model=train(training_df[,2:ncol(training)],training_df$iSize,method= "ranger",importance = 'permutation') #fit training
model

predictions_group3 =predict(object=model, testing[,2:ncol(testing)])
predictions_group3

Predicted_obs3=cbind(as.character(pheno[testpop3_idx,1]),predictions_group3,pheno[testpop3_idx,2])
colnames(Predicted_obs3)=c("Individual","Prediction","Observed")

VaiableImportance_obs3=varImp(model)
VaiableImportance_obs3_df=VaiableImportance_obs3$importance
VaiableImportance_obs3_df$SNP=rownames(VaiableImportance_obs3_df)
colnames(VaiableImportance_obs3_df)=c("SNP","Overall importance")

}




# Predict group 4
{
# Load genotypes for group 4
{
  geno=read.table(args[6],sep=",",header=T)  
  dim(geno)
  geno_shorter=geno[,3:(ncol(geno)-2)]
  geno_shorter_t=t(geno_shorter)
  rownames(geno_shorter_t)=colnames(geno_shorter)
}

#Combine genotypes and phenotypes into one large dataframe
{  
  if (all(rownames(geno_shorter_t)==pheno[,1])){
    cloverdata=cbind(as.numeric(as.character(iSize)),geno_shorter_t)
    colnames(cloverdata)=c("iSize",as.character(geno$identifier))
  }
}

# Machine learning using the fourth group as testing and group1,2,3,5,6 as training
training <- cloverdata[-testpop4_idx,]
training_df=as.data.frame(training)
testing <- cloverdata[testpop4_idx,]
testing_df <- as.data.frame(testing)


model=train(training_df[,2:ncol(training)],training_df$iSize,method= "ranger",importance = 'permutation') #fit training
model

predictions_group4 =predict(object=model, testing[,2:ncol(testing)])
predictions_group4

Predicted_obs4=cbind(as.character(pheno[testpop4_idx,1]),predictions_group4,pheno[testpop4_idx,2])
colnames(Predicted_obs4)=c("Individual","Prediction","Observed")

VaiableImportance_obs4=varImp(model)
VaiableImportance_obs4_df=VaiableImportance_obs4$importance
VaiableImportance_obs4_df$SNP=rownames(VaiableImportance_obs4_df)
colnames(VaiableImportance_obs4_df)=c("SNP","Overall importance")

}



# Predict group 5
{
# Load genotypes for group 5
{
  geno=read.table(args[7],sep=",",header=T)  
  dim(geno)
  geno_shorter=geno[,3:(ncol(geno)-2)]
  geno_shorter_t=t(geno_shorter)
  rownames(geno_shorter_t)=colnames(geno_shorter)
}

#Combine genotypes and phenotypes into one large dataframe
{  
  if (all(rownames(geno_shorter_t)==pheno[,1])){
    cloverdata=cbind(as.numeric(as.character(iSize)),geno_shorter_t)
    colnames(cloverdata)=c("iSize",as.character(geno$identifier))
  }
}


# Machine learning using the fifth group as testing and group1,2,3,4,6 as training
training <- cloverdata[-testpop5_idx,]
training_df=as.data.frame(training)
testing <- cloverdata[testpop5_idx,]
testing_df <- as.data.frame(testing)


model=train(training_df[,2:ncol(training)],training_df$iSize,method= "ranger",importance = 'permutation') #fit training
model

predictions_group5 =predict(object=model, testing[,2:ncol(testing)])
predictions_group5

Predicted_obs5=cbind(as.character(pheno[testpop5_idx,1]),predictions_group5,pheno[testpop5_idx,2])
colnames(Predicted_obs5)=c("Individual","Prediction","Observed")

VaiableImportance_obs5=varImp(model)
VaiableImportance_obs5_df=VaiableImportance_obs5$importance
VaiableImportance_obs5_df$SNP=rownames(VaiableImportance_obs5_df)
colnames(VaiableImportance_obs5_df)=c("SNP","Overall importance")


}

# Predict group 6
{
# Load genotypes for group 6
{
  geno=read.table(args[8],sep=",",header=T)  
  dim(geno)
  geno_shorter=geno[,3:(ncol(geno)-2)]
  geno_shorter_t=t(geno_shorter)
  rownames(geno_shorter_t)=colnames(geno_shorter)
}

#Combine genotypes and phenotypes into one large dataframe
{  
  if (all(rownames(geno_shorter_t)==pheno[,1])){
    cloverdata=cbind(as.numeric(as.character(iSize)),geno_shorter_t)
    colnames(cloverdata)=c("iSize",as.character(geno$identifier))
  }
}


# Machine learning using the sixth group as testing and group1,2,3,4,5 as training
training <- cloverdata[-testpop6_idx,]
training_df=as.data.frame(training)
testing <- cloverdata[testpop6_idx,]
testing_df <- as.data.frame(testing)


model=train(training_df[,2:ncol(training)],training_df$iSize,method= "ranger",importance = 'permutation') #fit training
model

predictions_group6 =predict(object=model, testing[,2:ncol(testing)])
predictions_group6

Predicted_obs6=cbind(as.character(pheno[testpop6_idx,1]),predictions_group6,pheno[testpop6_idx,2])
colnames(Predicted_obs6)=c("Individual","Prediction","Observed")

VaiableImportance_obs6=varImp(model)
VaiableImportance_obs6_df=VaiableImportance_obs6$importance
VaiableImportance_obs6_df$SNP=rownames(VaiableImportance_obs6_df)
colnames(VaiableImportance_obs6_df)=c("SNP","Overall importance")

}


#Combine
Combined=rbind(Predicted_obs1,Predicted_obs2,Predicted_obs3,Predicted_obs4,Predicted_obs5,Predicted_obs6)
Correlation=cor(as.numeric(as.character(Combined[,2])),as.numeric(as.character(Combined[,3])))
Combined_weights=rbind(VaiableImportance_obs1_df,VaiableImportance_obs2_df,VaiableImportance_obs3_df,VaiableImportance_obs4_df,VaiableImportance_obs5_df,VaiableImportance_obs6_df)


write.table(Combined,args[9],sep="\t",row.names=F,col.names=T,quote=F)
write.table(Correlation,args[10],sep="\t",row.names=F,col.names=T,quote=F)
write.table(Combined_weights,args[11],sep="\t",row.names=F,col.names=T,quote=F)