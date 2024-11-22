#!/usr/bin/env python
# coding: utf-8

# In[1]:


from statannot import add_stat_annotation
from scipy import stats
import csv
import numpy as np
import pandas as pd
import os
from os import path
import re
import shutil
import zipfile
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
import itertools
from itertools import combinations
from itertools import (takewhile,repeat)


# In[8]:


path_to_files='/home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Replicas_aa/VDJ/'
izlist=['IGH','IGHA','IGHG','IGHM']
sort_type='bulk'
orglist=['LN','tum','PBMC']

topn=100

files=[f for f in os.listdir(path_to_files) if re.match('.*'+sort_type+'.*'+'IGH'+'.*'+'.txt', f)]
files.sort()


# In[9]:


data = pd.DataFrame()

for f in files:
    df=pd.read_csv(path_to_files+f,comment="#", sep="\t", header=0)
    df=df.head(topn).copy()
    df['filename']=f
    data=data.append(df)
    
data['isotype']=[x.split('.')[-2] for x in data['filename']]
data['patient']=[x.split('_')[1] for x in data['filename']]
data['cancer']=[x.split('_')[0] for x in data['filename']]
data['organ']=[x.split('_')[2] for x in data['filename']]
data['sample']=[x.split('_')[4]+x.split('_')[5] for x in data['filename']]
data['cdr3len']=[len(x) for x in data['cdr3aa']]
    


# In[14]:


ord=['col','lun','mel']
test='t-test_ind' #suspect normal distibution, but it needs to be checked


for iz in izlist:
    figure, axes = plt.subplots(nrows=1, ncols=3,figsize=(20, 6))
    figure.suptitle('aa-CDR3 length, top '+str(topn)+', '+str(iz)+', '+test+', Bonferroni corr.',size=20)
    for i in range(len(orglist)):
        org=orglist[i]
        df=data[(data['isotype']==iz) & (data['organ']==org)].copy()
        df=df.groupby(['cancer','patient']).agg({'cdr3len':'mean'})
        #df=df.reset_index()
        #df=df.groupby(['cancer']).agg({'cdr3len':'mean'})
        df=df.reset_index()
        
        print(org,iz,'Number of patients',len(df['patient'].unique()))
        
        sns.scatterplot(x="cancer", y="cdr3len",alpha=0.8,data=df,color='maroon',s=150,ax=axes[i]).set_title(str(orglist[i]),fontsize=15)
        sns.boxplot(x="cancer", y="cdr3len", color='white',data=df,ax=axes[i]).set_title(str(orglist[i]),fontsize=20)
        axes[i].set(xlabel='', ylabel='mean of CDR3 length')
        axes[i].set_xlabel('',fontsize=15)
        axes[i].set_ylabel('mean of CDR3 length',fontsize=15)
        axes[i].tick_params(labelsize=15)
        test_results = add_stat_annotation(ax=axes[i], data=df, x="cancer", y="cdr3len", order=ord, box_pairs=[('col','lun'),('mel','lun'),('col','mel')],test=test,text_format='simple',loc='inside',fontsize='x-large', verbose=2)
        
        plt.savefig("Pictures/CDR3_length_for_different_cancers_"+iz+".png")
        plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/CDR3_length_for_different_cancers_"+iz+".pdf")


# In[16]:


ord=['LN','PBMC','tum']
box_pairs=list(itertools.combinations(ord, 2))
test='t-test_ind' #suspect normal distibution, but it needs to be checked

for iz in izlist:
    df=data[data['isotype']==iz].copy()
    #df=df[df['cancer']=='lun']
    df=df.groupby(['cancer','patient','organ']).agg({'cdr3len':'mean'})
    df=df.reset_index()
    
    df=df[df['organ']!='norm']
    
    print('Number of patients',len(df['patient'].unique()))
    
    figure, axes = plt.subplots(nrows=1, ncols=1,figsize=(6, 6))
    
    sns.stripplot(x="organ", y="cdr3len",alpha=0.7,data=df,color='maroon',s=15).set_title(str(iz))
    sns.boxplot(x="organ", y="cdr3len", color='white',data=df).set_title('CDR3aa length, top '+str(topn)+' clonotypes, '+str(iz)+',\n'+test+', Bonferroni cor.',fontsize=15)
    sns.color_palette("husl", 9)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    axes.set(xlabel='', ylabel='mean of CDR3 length')
    test_results = add_stat_annotation(ax=axes, data=df, x="organ", y="cdr3len", order=ord, box_pairs=box_pairs,test=test,text_format='simple',loc='inside',fontsize='x-large', verbose=2)
    
    axes.set_xlabel('Tissue',fontsize=15)
    axes.set_ylabel('CDR3 length, amino acids',fontsize=15)
    axes.tick_params(labelsize=15)
    #axes.legend(fontsize='x-large')
    
    plt.savefig("Pictures/CDR3_length_for_different_tissues"+iz+".png")
    plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/CDR3_length_for_different_tissues"+iz+".pdf")


# In[6]:


df=data[data['isotype']=='IGH'].copy()
df=df[df['organ']=='LN']
dist=df.groupby(['patient','sample']).agg({'cdr3len':'mean'}).reset_index()



sample=np.array(dist.cdr3len)
F,p=stats.shapiro(sample)
print(F,p)

fig, axes=plt.subplots(nrows=1,ncols=1,figsize=(8,6))

plt.hist(dist.cdr3len,color='dimgray')
plt.title('Mean amino acid CDR3 length distribution',fontsize=20)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.savefig("Pictures/CDR3_length_distribution.png")


# In[6]:


#NOT SURE IF PLOTTING DIFFERENT ISOTYPES IS FAIR - THEIR CDR3 LENGTH MUST BE HIGHLY DEPENDANT
ord=['LN','PBMC','tum']
box_pairs=list(itertools.combinations(ord, 2))
test='t-test_ind'

df=data.copy()
    #df=df[df['cancer']=='lun']
df=df.groupby(['cancer','patient','organ','isotype']).agg({'cdr3len':'mean'})
df=df.reset_index()
    
df=df[df['organ']!='norm']
df=df[df['isotype']!='IGH']
    
figure, axes = plt.subplots(nrows=1, ncols=1,figsize=(6, 6))
    
sns.scatterplot(x="organ", y="cdr3len",alpha=0.8,data=df,hue='isotype',color='maroon',s=100).set_title(str(iz))
sns.boxplot(x="organ", y="cdr3len", color='white',data=df).set_title('CDR3aa length, top '+str(topn)+' clonotypes,\n'+test+', Bonferroni cor.',fontsize=15)
sns.color_palette("husl", 9)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
axes.set(xlabel='', ylabel='mean of CDR3 length')
test_results = add_stat_annotation(ax=axes, data=df, x="organ", y="cdr3len", order=ord, box_pairs=box_pairs,test=test,text_format='simple',loc='inside', verbose=2)
    
axes.set_xlabel('Tissue',fontsize=15)
axes.set_ylabel('CDR3 length',fontsize=15)
axes.tick_params(labelsize=15)
axes.legend(fontsize='x-large')
    
plt.savefig("Pictures/CDR3_length_for_different_tissues_all_isotypes_samples.png")


# In[18]:


ord=['LN','PBMC','tum','norm']

for iz in izlist:
    df=data[data['isotype']==iz].copy()
    #df=df.groupby(['cancer','patient','organ']).agg({'cdr3len':'min'})
    df=df.reset_index()
    
    #df=df[df['organ']!='norm']
    
    figure, axes = plt.subplots(nrows=1, ncols=1,figsize=(6, 8))
    
    #sns.scatterplot(x="organ", y="cdr3len",alpha=0.8,data=df,hue='cancer',color='maroon',s=100).set_title(str(iz))
    sns.violinplot(x="organ", y="cdr3len",data=df).set_title('CDR3aa length, top '+str(topn)+' clonotypes, '+str(iz),fontsize=15)
    sns.boxplot(x="organ", y="cdr3len",color='white',data=df).set_title('CDR3aa length, top '+str(topn)+' clonotypes, '+str(iz),fontsize=15)
    sns.color_palette("husl", 9)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    axes.set(xlabel='', ylabel='CDR3 lengths')
    test_results = add_stat_annotation(ax=axes, data=df, x="organ", y="cdr3len", order=ord, box_pairs=[('LN','PBMC'),('LN','tum'),('PBMC','tum'),('LN','norm'),('norm','tum'),('PBMC','norm')],test='t-test_ind',text_format='simple',loc='inside', verbose=2)
        
    plt.savefig("Pictures/CDR3_length_violinplots_with_norm"+iz+".png")


# In[ ]:




