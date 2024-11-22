#!/usr/bin/env python
# coding: utf-8

# In[11]:


#LIBRARIES
from statannot import add_stat_annotation
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import csv
import numpy as np
import pandas as pd
import os
from os import path
import re
import shutil
import zipfile
import math
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
import itertools
from itertools import combinations
from itertools import (takewhile,repeat)
import glob
import ternary

import warnings
warnings.filterwarnings("ignore") 

def union_dataframes(path_to_files,pattern):
    df=pd.DataFrame()
    files=[f for f in os.listdir(path_to_files) if re.match(pattern,f)]
    for f in files:
        data=pd.read_csv(path_to_files+f,sep='\t')
        data['filename']=f
        df=df.append(data)
    df['patient']=[x.split('_')[1] for x in df['1_sample_id']]
    df['sort_type']=[x.split('_')[3] for x in df['1_sample_id']]
    df['organ_1']=[x.split('_')[2] for x in df['1_sample_id']]
    df['organ_2']=[x.split('_')[2] for x in df['2_sample_id']]
    df['sample_1']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5] for x in df['1_sample_id']]
    df['sample_2']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5] for x in df['2_sample_id']]
    df['rep_1']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5]+x.split('_')[6].split('.')[0] for x in df['1_sample_id']]
    df['rep_2']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5]+x.split('_')[6].split('.')[0] for x in df['2_sample_id']]
    return df

def comparison_type(s1,s2,org1,org2):
    if s1==s2:
        ans='same'
    elif org1==org2:
        ans='different'
    else:
        ans='tissues'
    return ans

def tissue_pair(x,y):
    l=[x,y]
    l.sort()
    ans=l[0]+'/'+l[1]
    return ans

def segment_pair(x,y):
    l=[x,y]
    l.sort()
    ans=l[0]+'/'+l[1]
    return ans


# In[12]:


#DATA READING AND ALTERING

data=union_dataframes('/home/sofyakrasik/PwDistances/SeparateReplicas_with_singles/','.*_IGH.intersect.batch.aa.txt') #_with_singles
#data=union_dataframes('/home/sofyakrasik/PwDistances/SeparateReplicas/','.*_IGH.intersect.batch.aa.txt') # no singles
data['comparison_type']=[comparison_type(x,y,z,q) for x,y,z,q in zip(data['sample_1'],data['sample_2'],data['organ_1'],data['organ_2'])]
data['tissue_pair']=[tissue_pair(x,y) for x,y in zip(data['organ_1'],data['organ_2'])]
data['segment_pair']=[tissue_pair(x,y) for x,y in zip(data['sample_1'],data['sample_2'])]

patients=list(data['patient'].unique())
patients.sort()

#data[data['tissue_pair']=='LN/tum']


# In[13]:


#ALL SEGMENT OVERLAP WITH REPLICATES + ANOVA 
ts1='LN'
ts2='tum' #must be in ts1 ts2 alphabetical order

for pat in patients:
    df=data[data['patient']==pat].copy()
    df=df[df['comparison_type']=='tissues']
    df=df[df['tissue_pair']==ts1+'/'+ts2]
    
    if (len(df)>0) & (len(list(df['segment_pair'].unique()))>1):
        fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(10, 8))
        
        result = df.groupby('segment_pair')['F2'].apply(list)
        F, p = stats.f_oneway(*result)
        
        sns.scatterplot(x="segment_pair", y="F2",alpha=0.5,data=df,s=200, linewidth=0,hue='sample_1').set_title(pat+', ANOVA p-value='+str(round(p,5)),fontsize=15)
        
        axes.set_xlabel('',fontsize=15)
        axes.set_ylabel('',fontsize=15)
        axes.tick_params(labelsize=12,axis='x', rotation=45)
        axes.tick_params(labelsize=15,axis='y')
        
        axes.legend(fontsize='x-large', title_fontsize='40')


# In[14]:


#CHECK FOR ANOVA (SHAPIRO-WILK?) 
for pat in patients:
    df=data[data['patient']==pat].copy()
    df=df[df['comparison_type']=='tissues']
    df=df[df['tissue_pair']==ts1+'/'+ts2]
    
    df=df[df['patient']=='ccp3']
    
    for p in list(df['segment_pair'].unique()):
        sample=np.array(df[df['segment_pair']==p].F2)
        F,p=stats.shapiro(sample)
        print(F,p)
        
#P.S. not really helpful for sample size of 4


# In[15]:


#MAYBE BOOTSTRAP 100 TIMES OF 100 CLONOTYPES FROM EACH SEGMENT TO PROVE NORMALITY OF ITS DISTRIBUTION


# In[16]:


#ALL SEGMENTS VS ALL SEGMENTS: COMPARE MEANS FOR LN SEGMENTS AND FOR TUMOR (INTER-TISSUE HETEROGENEITY)
#lcp3 and mp9 are done with signles, otherwise tum 13 is missing and mp9 has no overlap 
#patients=['ccp2','ccp3','ccp5','ccp6','lcp1','lcp2','lcp3','lcp4','mp1','mp2','mp3','mp7','mp8','mp9']
patients=['ccp2','ccp3','ccp5','ccp6','lcp1','lcp2','lcp4','mp1','mp2','mp3','mp7','mp8'] #patietns done without singles
patients=['lcp3','mp9'] #patients done with singles, otherwise some samples are lost

ts1='LN'
ts2='tum' #must be in ts1 ts2 alphabetical order

category='sample_1'
hue='sample_2'

for pat in patients:
    df=data[data['patient']==pat].copy()
    df=df[df['comparison_type']=='tissues']
    df=df[df['tissue_pair']==ts1+'/'+ts2]
    
    #df=df[df['patient']=='mp9']
    
    if (len(df)>0) & (len(list(df[category].unique()))>1):
        fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(6, 6))
        
        result = df.groupby(category)['F2'].apply(list)
        F, p = stats.f_oneway(*result)
        
        ord=list(df[category].unique())
        box_pairs=list(itertools.combinations(ord, 2))
        test='Mann-Whitney'
        
        sns.stripplot(x=category, y="F2",alpha=0.5,data=df,s=15, linewidth=0,hue=hue,palette='bright').set_title(pat+', ANOVA p-value='+str(round(p,2)),fontsize=15)
        sns.boxplot(x=category, y="F2",data=df,color='white').set_title(pat+', LN fragments overlap with tum,\n'+test+', Bonferroni corr',fontsize=15)
        test_results = add_stat_annotation(ax=axes,data=df, x=category, y='F2', order=ord, box_pairs=box_pairs,test=test,text_format='simple',loc='inside',fontsize='x-large', verbose=2)
        
        axes.set_xlabel('',fontsize=15)
        axes.set_ylabel('F2',fontsize=15)
        axes.tick_params(labelsize=15,axis='x')
        axes.tick_params(labelsize=15,axis='y')
        axes.legend(fontsize='x-large', title_fontsize='40')
        
        plt.savefig('/home/sofyakrasik/Pipelines/Pictures/Heterogeneity_F2_'+pat+'_LN_segments_vs_tum.png',bbox_inches='tight') #tum segments vs LN should be done with singles, otherwise lcp3 tum3 is lost
        plt.savefig('/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/Heterogeneity_F2_'+pat+'_tum_segments_vs_LN.pdf',bbox_inches='tight') #_tum_segments_vs_LN _LN_segments_vs_tum


# In[49]:


#REPLICATES F2 vs ALL SEGMENTS PAIRS F2 (INTRA-TISSUE HETEROGENEITY) FOR LN AND TUM
#tum is done with singles
ts='LN'


for pat in patients:
    df=data[data['patient']==pat].copy()
    #df=df[df['comparison_type']]
    df=df[df['tissue_pair']==ts+'/'+ts]
    
    df=df.sort_values(['comparison_type'])

    
    if (len(df)>0) & (len(list(df['comparison_type'].unique()))>1):
        fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(6, 6))
        
        #result = df.groupby(category)['F2'].apply(list)
        #F, p = stats.f_oneway(*result)
        
        df=df.sort_values(['comparison_type'],ascending=False)
        
        ord=list(df['comparison_type'].unique())
        box_pairs=list(itertools.combinations(ord, 2))
        test='Mann-Whitney'
        
        sns.stripplot(x='comparison_type', y="F2",alpha=0.5,data=df,s=15, linewidth=0,hue='comparison_type',palette='Set1').set_title(pat+', ANOVA p-value='+str(p),fontsize=15)
        #sns.stripplot(x='comparison_type', y="F2",alpha=0.5,data=df,s=15, linewidth=0).set_title(pat+', ANOVA p-value='+str(p),fontsize=15)
        sns.boxplot(x='comparison_type', y="F2",data=df,color='white').set_title(ts+' samples repertoire overlap, '+pat,fontsize=15)
        test_results = add_stat_annotation(ax=axes,data=df, x='comparison_type', y='F2', order=ord, box_pairs=box_pairs,test=test,text_format='simple',loc='inside',fontsize='x-large', verbose=2)
        
        axes.set_xlabel('',fontsize=15)
        axes.set_ylabel('',fontsize=15)
        axes.tick_params(labelsize=15,axis='x')
        axes.tick_params(labelsize=15,axis='y')
        axes.legend(fontsize='x-large', title_fontsize='40')
        axes.legend([],[], frameon=False)
        
        plt.savefig('/home/sofyakrasik/Pipelines/Pictures/Heterogeneity_F2_'+pat+'_within_'+ts+'.png',bbox_inches='tight')


# In[177]:


df


# In[ ]:




