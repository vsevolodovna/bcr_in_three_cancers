#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[16]:


path_to_data="/home/sofyakrasik/Hypermutations/ChangeOut/"
files=glob.glob(path_to_data+"*_pooled_replicas_alldata_CDR3_clone-pass.tsv")

topn=30

data=pd.DataFrame()

for f in files:
    df=pd.read_csv(f,sep='\t')
    df['filename']=f.split('/')[-1]
    data=data.append(df)


# In[17]:


#data=data[data['allCHitsWithScore'].notna()]
data['cancer']=[x.split('_')[0] for x in data['basename']]
data['patient']=[x.split('_')[1] for x in data['basename']]
data['organ']=[x.split('_')[2] for x in data['basename']]
data['sample']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5] for x in data['basename']]
#data['sample']=[x.split('_')[2]+x.split('_')[4] for x in data['basename']]
#data['isotype']=[x[3] for x in data['allCHitsWithScore']]
data['mutnum']=[x.split(";")[0].count('S')+y.split(";")[0].count('S') for x,y in zip(data['allvalignments'],data['alljalignments'])]
data=data.sort_values(by=['patient'])

#data=data[data['patient']!='mp1']  #because only PBMC tissue

data['cluster_size']=data.groupby(['filename','clone_id'])['sequence'].transform('count')


# In[21]:


#ALL PATIENTS TOGETHER

topn=100

#filter out those without tree
treshold=1
data=data[data['cluster_size']>=treshold].copy()

groups='cancer'
isotype_list=['A','G','M','D','E']

for iz in isotype_list:
    df=data.copy()
    df=df[df['isotype'].isin([iz])]
    #df=df[df['organ'].isin(['tum','norm'])]
    #df=df[df['cancer']=='lun']
    df=df[df['organ']!='norm']
    #df=df[df['organ']=='tum']
    ord=list(df[groups].unique())
    box_pairs=list(itertools.combinations(ord, 2))
    
    df=df.sort_values([groups,'clonefraction'],ascending=False)
    df=df.groupby(['cancer','patient',groups,'sample']).head(topn)

    #df['numisotypes']=df.groupby(['patient','targetSequences'])['isotype'].transform('nunique')
    #df=df[df['numisotypes']>1]

    #df=df[df['mutnum']!=0]
    #df=df[df['organ']=='tum'].copy()
    df=df.groupby(['cancer','patient','organ','sample']).agg({'mutnum':'mean','sequence':'count'})
    df=df.reset_index()

    figure, axes = plt.subplots(nrows=1, ncols=1,figsize=(8, 6))
    sns.stripplot(x=groups, y='mutnum',alpha=0.5,data=df,s=15,hue='organ').set_title(groups)
    #sns.violinplot(x=groups, y='mutnum',alpha=0.3,data=df)
    sns.boxplot(x=groups, y='mutnum', color='white',data=df).set_title("Mean number of hypermutations, Ig"+iz+", Mann-Whitney test,\ncluster size threshold="+str(treshold)+', top '+str(topn),fontsize=15)
    test_results = add_stat_annotation(ax=axes, data=df, x=groups, y='mutnum', order=ord, box_pairs=box_pairs,test='Mann-Whitney',text_format='simple',loc='inside',fontsize='x-large', verbose=2)

    axes.set_xlabel(groups,fontsize=15)
    axes.set_ylabel(groups,fontsize=15)
    axes.tick_params(labelsize=15)
    
    df1=df.groupby([groups])['mutnum'].mean()
    df1=df1.to_frame().reset_index()
    df1['label']=['-' if x==0 else str(round(x,1)) for x in df1['mutnum']]
    df1=df1.sort_values(groups).copy()
    
    #QUANTILE PLOTTING
    #the_table=plt.table(cellText=[list(df1['label'])],loc='upper center',cellLoc='center',edges='open',fontsize=20)
    
    #the_table.set_fontsize(20)
    plt.xticks(fontsize=20)
    
    
    #plt.savefig("Pictures/Mutation_number_by_"+groups+"_"+iz+".png",bbox_inches='tight')
    plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/Mutation_number_by_"+groups+"_"+iz+".pdf",bbox_inches='tight')


# In[ ]:




