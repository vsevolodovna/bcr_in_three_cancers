#!/usr/bin/env python
# coding: utf-8

# In[1]:


from statannot import add_stat_annotation
from scipy import stats
from scipy.stats import pearsonr
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


# In[11]:


sort_type='bulk'
path_to_files='/home/sofyakrasik/BCR_CDR3_data/5RACE/All_data_filtered/'

izlist=['A','G','M']

org1="tum"
org2="LN"

top_min=100 #We do not use info from too small files

files=[f for f in os.listdir(path_to_files) if re.match('.*'+sort_type+'.*'+'IGH.txt', f)]
files.sort()


# In[12]:


cols=['cloneFraction','allCHitsWithScore','filename']

mixcr_data = pd.DataFrame(columns=cols)

for f in files:
    data=pd.read_csv(path_to_files+f,comment="#", sep="\t", header=0)
    if len(data)>=top_min:
        data['filename']=f
        data=data[['cloneFraction','allCHitsWithScore','filename']].copy()
        mixcr_data=mixcr_data.append(data)


# In[9]:


for iz in izlist:
    patdata=mixcr_data.dropna().copy()
    #FILTERING FOR ISOTYPE
    patdata=patdata[patdata['allCHitsWithScore'].str.contains('IGH'+iz)].copy()
    patdata=patdata.groupby('filename')['cloneFraction'].agg('sum')
    patdata=patdata.to_frame().reset_index()
    
    patdata['cancer']=[x.split('_')[0] for x in patdata['filename']]
    patdata['patient']=[x.split('_')[1] for x in patdata['filename']]
    patdata['organ']=[x.split('_')[2] for x in patdata['filename']]
    patdata['sample']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5] for x in patdata['filename']]
    patdata=patdata.groupby(['cancer','patient','organ','sample'])['cloneFraction'].agg('mean')
    patdata=patdata.to_frame().reset_index()
    patdata=patdata.groupby(['cancer','patient','organ'])['cloneFraction'].agg('mean')
    patdata=patdata.to_frame().reset_index()
    
    df_1=patdata[patdata['organ']==org1]
    df_2=patdata[patdata['organ']==org2]
    df=pd.merge(df_1,df_2,how='inner',on='patient')

    plt.figure(figsize=(5,5))
    p=sns.regplot(x="cloneFraction_x", y="cloneFraction_y",data=df,ci=95,color='grey')
    
    corr,pval=pearsonr(df['cloneFraction_x'],df['cloneFraction_y'])
    sns.scatterplot(x="cloneFraction_x", y="cloneFraction_y",hue='cancer_x',edgecolors="black",linewidth=1,alpha=0.9,
                    data=df,s=100).set_title(iz+",isotype fraction correlation,"+'\n'+"corr="+str(round(corr,2))+", p-value="+str(round(pval,2))+", 95% ci",fontsize=15)
    
    plt.xlabel(org1+' isotype fraction',fontsize=15)
    plt.ylabel(org2+' isotype fraction',fontsize=15)
    plt.tick_params(labelsize=15)
    plt.legend(fontsize='x-large', title_fontsize='40')
    
    plt.savefig('Pictures/Isotype_corr_'+org1+'_'+org2+'_'+iz+'.png',bbox_inches='tight')
    plt.savefig('/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/Isotype_corr_'+org1+'_'+org2+'_'+iz+'.pdf',bbox_inches='tight')


# In[13]:


#CLONOTYPES TOGETHER OR STATISTICAL TRICK WITH IGE IGD
patdata=mixcr_data.dropna().copy()
patdata['isotype']=[x[3:4] for x in patdata['allCHitsWithScore']]
patdata=patdata[patdata['isotype'].isin(izlist)]
patdata['isotype']=['Ig'+x for x in patdata['isotype']]
patdata=patdata.groupby(['filename','isotype'])['cloneFraction'].agg('sum')
patdata=patdata.to_frame().reset_index()

patdata['cancer']=[x.split('_')[0] for x in patdata['filename']]
patdata['patient']=[x.split('_')[1] for x in patdata['filename']]
patdata['organ']=[x.split('_')[2] for x in patdata['filename']]
patdata['sample']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5] for x in patdata['filename']]
patdata=patdata.groupby(['cancer','patient','organ','sample','isotype'])['cloneFraction'].agg('mean')
patdata=patdata.to_frame().reset_index()
patdata=patdata.groupby(['cancer','patient','organ','isotype'])['cloneFraction'].agg('mean')
patdata=patdata.to_frame().reset_index()
    
df_tum=patdata[patdata['organ']==org1]
df_PBMC=patdata[patdata['organ']==org2]
df=pd.merge(df_tum,df_PBMC,how='inner',on=['patient','isotype'])

plt.figure(figsize=(5,5))
p=sns.regplot(x="cloneFraction_x", y="cloneFraction_y",data=df,ci=95,color='grey')
    
corr,pval=pearsonr(df['cloneFraction_x'],df['cloneFraction_y'])
sns.scatterplot(x="cloneFraction_x", y="cloneFraction_y",hue='isotype',edgecolors="black",linewidth=1,alpha=0.9,
                data=df,s=100).set_title("Isotype fraction correlation,"+'\n'+"corr="+str(round(corr,2))+", p-value="+str(round(pval,6))+", 95% ci",fontsize=15)

plt.xlabel(org1+' isotype fraction',fontsize=15)
plt.ylabel(org2+' isotype fraction',fontsize=15)
plt.tick_params(labelsize=15)
plt.legend(fontsize='x-large', title_fontsize='40')
    
plt.savefig('Pictures/Isotype_corr_all_isotypes_'+org1+'_'+org2+'.png',bbox_inches='tight')
plt.savefig('/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/Isotype_corr_all_isotypes_'+org1+'_'+org2+'.pdf',bbox_inches='tight')

