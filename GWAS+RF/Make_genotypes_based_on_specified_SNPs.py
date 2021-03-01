#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 10:59:46 2019

@author: CathrineKiel
"""

import os
import pandas as pd
import numpy as np
from sys import argv

file=pd.read_csv(argv[1],header=0)
#"_pid1_Yellow.tip_emmax_none_t75_round2.pvals"
# first filter out MAF<0.05

row,col=file.shape

remove=[]
for i in range(0,row):
    
    if file.mafs[i]<0.05:
        remove.append(i)

len(remove)
MAFfiltered=file.drop(file.index[[remove]])

len(MAFfiltered)

# Print a list of the 25 SNPs to use i genotype file
MAFfiltered.sort_values('scores',inplace=True)
top25=MAFfiltered.head(25)
headers=list(top25.columns)
headers=str(headers)

np.savetxt(argv[3],top25,fmt='%s',delimiter=',',header=headers)
#YellowTipRound2top200snps.txt"


# Print a list of the 200 SNPs to use i genotype file
MAFfiltered.sort_values('scores',inplace=True)
top200=MAFfiltered.head(200)
headers=list(top200.columns)
headers=str(headers)

np.savetxt(argv[4],top200,fmt='%s',delimiter=',',header=headers)
#YellowTipRound2top200snps.txt"


# Print a list of the 500 SNPs to use i genotype file
MAFfiltered.sort_values('scores',inplace=True)
top500=MAFfiltered.head(500)
headers=list(top500.columns)
headers=str(headers)

np.savetxt(argv[5],top500,fmt='%s',delimiter=',',header=headers)
#YellowTipRound2top200snps.txt"


# Print a list of 25 random SNPs to use i genotype file
random25=MAFfiltered.sample(n=25)
headers=list(random25.columns)
headers=str(headers)

np.savetxt(argv[6],random25,fmt='%s',delimiter=',',header=headers)
#YellowTipRound2top200snps.txt"



# Print a list of 200 random SNPs to use i genotype file
random200=MAFfiltered.sample(n=200)
headers=list(random200.columns)
headers=str(headers)

np.savetxt(argv[7],random200,fmt='%s',delimiter=',',header=headers)
#YellowTipRound2top200snps.txt"

# Print a list of 500 random SNPs to use i genotype file
random500=MAFfiltered.sample(n=500)
headers=list(random500.columns)
headers=str(headers)

np.savetxt(argv[8],random500,fmt='%s',delimiter=',',header=headers)
#YellowTipRound2top200snps.txt"



# Now load in genotype file and remove SNPs not used
# We need to have the genotype file with all individuals, not just training
geno = pd.read_csv(argv[2],header=0)
#"/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/Preparation_of_genotypeFile/Geno_LDfiltered0.5and1.00_20200728.csv"

# Remember to delete the individuals with no phenotype records and the 6 outlier individuals
dropcolumns=["Aaran_06","Aoost_03","Aoost_06","Banna_05","Ccyma_05","Ccyma_07","Ccyma_08","Ctain_08","Kdike_05","Kdike_10","Llanc_01","Llanc_02","Llanc_10","Mrida_03","Volin_06","Volin_08","Volin_10"]
geno=geno.drop(dropcolumns,axis=1)



looking_for=top25.ix[:,0:2]

geno_quick=np.array(geno) #converted genotype dataframe to numpy array for quicker looping.

indexes_needed_in_geno=[]

for i in range(0,25):
    current_suspect=looking_for.loc[looking_for.index[i],'positions']
    index=np.where(geno_quick == current_suspect)
    
    if len(index[0])==0:
        print(current_suspect,"not in genotype file")
    
    if len(index[0])==1:
        indexes_needed_in_geno.append(index[0][0])
        print(current_suspect,"found")
    
    if len(index[0])>1: #this is if the position ID is not unique
        chromosome=looking_for.loc[looking_for.index[i],'chromosomes']
        
        lst1=list(np.where(geno_quick == current_suspect))
        lst2=list(np.where(geno_quick == chromosome))
        
        c = [x for x in lst1[0] if x in lst2[0]]
        indexes_needed_in_geno.append(c[0])


indexes_needed_in_geno
len(indexes_needed_in_geno)
newgeno = geno_quick[indexes_needed_in_geno]
newgeno.shape
newgeno=np.array(newgeno)
newgeno_pd=pd.DataFrame(newgeno, columns = list(geno.columns.values))

# Save file
newgeno_pd.to_csv(argv[9], index = None, header= 1)






looking_for=top200.ix[:,0:2]

geno_quick=np.array(geno) #converted genotype dataframe to numpy array for quicker looping.

indexes_needed_in_geno=[]

for i in range(0,200):
    current_suspect=looking_for.loc[looking_for.index[i],'positions']
    index=np.where(geno_quick == current_suspect)
    
    if len(index[0])==0:
        print(current_suspect,"not in genotype file")
    
    if len(index[0])==1:
        indexes_needed_in_geno.append(index[0][0])
        print(current_suspect,"found")
    
    if len(index[0])>1: #this is if the position ID is not unique
        chromosome=looking_for.loc[looking_for.index[i],'chromosomes']
        
        lst1=list(np.where(geno_quick == current_suspect))
        lst2=list(np.where(geno_quick == chromosome))
        
        c = [x for x in lst1[0] if x in lst2[0]]
        indexes_needed_in_geno.append(c[0])


indexes_needed_in_geno
len(indexes_needed_in_geno)
newgeno = geno_quick[indexes_needed_in_geno]
newgeno.shape
newgeno=np.array(newgeno)
newgeno_pd=pd.DataFrame(newgeno, columns = list(geno.columns.values))

# Save file
newgeno_pd.to_csv(argv[10], index = None, header= 1)





# Now load in genotype file and remove SNPs not used
looking_for=top500.ix[:,0:2]

geno_quick=np.array(geno) #converted genotype dataframe to numpy array for quicker looping.

indexes_needed_in_geno=[]

for i in range(0,500):
    current_suspect=looking_for.loc[looking_for.index[i],'positions']
    index=np.where(geno_quick == current_suspect)
    
    if len(index[0])==0:
        print(current_suspect,"not in genotype file")
    
    if len(index[0])==1:
        indexes_needed_in_geno.append(index[0][0])
        print(current_suspect,"found")
    
    if len(index[0])>1: #this is if the position ID is not unique
        chromosome=looking_for.loc[looking_for.index[i],'chromosomes']
        
        lst1=list(np.where(geno_quick == current_suspect))
        lst2=list(np.where(geno_quick == chromosome))
        
        c = [x for x in lst1[0] if x in lst2[0]]
        indexes_needed_in_geno.append(c[0])


indexes_needed_in_geno
len(indexes_needed_in_geno)
newgeno = geno_quick[indexes_needed_in_geno]
newgeno.shape
newgeno=np.array(newgeno)
newgeno_pd=pd.DataFrame(newgeno, columns = list(geno.columns.values))

# Save file
newgeno_pd.to_csv(argv[11], index = None, header= 1)





# Now load in genotype file and remove SNPs not used
looking_for=random25.ix[:,0:2]

geno_quick=np.array(geno) #converted genotype dataframe to numpy array for quicker looping.

indexes_needed_in_geno=[]

for i in range(0,25):
    current_suspect=looking_for.loc[looking_for.index[i],'positions']
    index=np.where(geno_quick == current_suspect)
    
    if len(index[0])==0:
        print(current_suspect,"not in genotype file")
    
    if len(index[0])==1:
        indexes_needed_in_geno.append(index[0][0])
        print(current_suspect,"found")
    
    if len(index[0])>1: #this is if the position ID is not unique
        chromosome=looking_for.loc[looking_for.index[i],'chromosomes']
        
        lst1=list(np.where(geno_quick == current_suspect))
        lst2=list(np.where(geno_quick == chromosome))
        
        c = [x for x in lst1[0] if x in lst2[0]]
        indexes_needed_in_geno.append(c[0])


indexes_needed_in_geno
len(indexes_needed_in_geno)
newgeno = geno_quick[indexes_needed_in_geno]
newgeno.shape
newgeno=np.array(newgeno)
newgeno_pd=pd.DataFrame(newgeno, columns = list(geno.columns.values))

# Save file
newgeno_pd.to_csv(argv[12], index = None, header= 1)







looking_for=random200.ix[:,0:2]

geno_quick=np.array(geno) #converted genotype dataframe to numpy array for quicker looping.

indexes_needed_in_geno=[]

for i in range(0,200):
    current_suspect=looking_for.loc[looking_for.index[i],'positions']
    index=np.where(geno_quick == current_suspect)
    
    if len(index[0])==0:
        print(current_suspect,"not in genotype file")
    
    if len(index[0])==1:
        indexes_needed_in_geno.append(index[0][0])
        print(current_suspect,"found")
    
    if len(index[0])>1: #this is if the position ID is not unique
        chromosome=looking_for.loc[looking_for.index[i],'chromosomes']
        
        lst1=list(np.where(geno_quick == current_suspect))
        lst2=list(np.where(geno_quick == chromosome))
        
        c = [x for x in lst1[0] if x in lst2[0]]
        indexes_needed_in_geno.append(c[0])


indexes_needed_in_geno
len(indexes_needed_in_geno)
newgeno = geno_quick[indexes_needed_in_geno]
newgeno.shape
newgeno=np.array(newgeno)
newgeno_pd=pd.DataFrame(newgeno, columns = list(geno.columns.values))

# Save file
newgeno_pd.to_csv(argv[13], index = None, header= 1)



# Now load in genotype file and remove SNPs not used
looking_for=random500.ix[:,0:2]

geno_quick=np.array(geno) #converted genotype dataframe to numpy array for quicker looping.

indexes_needed_in_geno=[]

for i in range(0,500):
    current_suspect=looking_for.loc[looking_for.index[i],'positions']
    index=np.where(geno_quick == current_suspect)
    
    if len(index[0])==0:
        print(current_suspect,"not in genotype file")
    
    if len(index[0])==1:
        indexes_needed_in_geno.append(index[0][0])
        print(current_suspect,"found")
    
    if len(index[0])>1: #this is if the position ID is not unique
        chromosome=looking_for.loc[looking_for.index[i],'chromosomes']
        
        lst1=list(np.where(geno_quick == current_suspect))
        lst2=list(np.where(geno_quick == chromosome))
        
        c = [x for x in lst1[0] if x in lst2[0]]
        indexes_needed_in_geno.append(c[0])


indexes_needed_in_geno
len(indexes_needed_in_geno)
newgeno = geno_quick[indexes_needed_in_geno]
newgeno.shape
newgeno=np.array(newgeno)
newgeno_pd=pd.DataFrame(newgeno, columns = list(geno.columns.values))

# Save file
newgeno_pd.to_csv(argv[14], index = None, header= 1)


