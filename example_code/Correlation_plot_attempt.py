#!/usr/bin/env python
# coding: utf-8

# In[1]:


#add metadata input, isotypes
import csv
import shutil
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
from itertools import combinations

VDJTOOLS='/software/bin/vdjtools'


# In[2]:


#ДОБАВИТЬ ЧТЕНИЕ МЕТАДАТЫ
file_pairs=["/home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Samples_aa/mixcr/col_ccp2_tum_bulk_pooled.clonotypes.IGH.txt","/home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Samples_aa/mixcr/col_ccp2_tum_PCsort_pooled.clonotypes.IGH.txt"]

isotypes=['IGH','IGHA','IGHG','IGHM']

topn=400


# In[4]:


data1=pd.read_csv(file_pairs[0],comment="#", sep="\t", header=0)
data2=pd.read_csv(file_pairs[1],comment="#", sep="\t", header=0)

fname1=file_pairs[0].split("/")[-1]
fname2=file_pairs[1].split("/")[-1]

fmin=1
fmax=0

fig, axes = plt.subplots(nrows=2, ncols=2,figsize=(12, 12))
fig.suptitle(fname1.split('.')[0]+' vs '+fname2.split('.')[0]+', top '+str(topn),fontsize=20)


for i in range(len(isotypes)):
    iz=isotypes[i]
    df1=data1[data1['allCHitsWithScore'].str.contains(iz, na=False)].copy()
    df2=data2[data2['allCHitsWithScore'].str.contains(iz, na=False)].copy()
    
    df1=df1.sort_values(['cloneFraction'],ascending=False)
    df2=df2.sort_values(['cloneFraction'],ascending=False)
    
    df1['cloneFraction']=df1.groupby(['aaSeqCDR3'])['cloneFraction'].transform('sum')
    df1=df1.drop_duplicates(['aaSeqCDR3','cloneFraction'])
    df2['cloneFraction']=df2.groupby(['aaSeqCDR3'])['cloneFraction'].transform('sum')
    df2=df2.drop_duplicates(['aaSeqCDR3','cloneFraction'])
    
    c=min(len(df1),len(df2),topn)
    df1=df1.head(c)
    df2=df2.head(c)
    
    df1['cloneFraction']=df1['cloneFraction']/df1['cloneFraction'].sum()
    df2['cloneFraction']=df2['cloneFraction']/df2['cloneFraction'].sum()
    
    df1['isotype']=df1['cloneFraction']/df1['cloneFraction'].sum()
    
    df=pd.merge(df1,df2,on=['aaSeqCDR3'],how='outer')
    
    overlap_percent=round(len(df[(df['cloneFraction_x'].notna()) & (df['cloneFraction_y'].notna())])/c,2)
    overlap_percent=overlap_percent*100
    overlap_clones=df[(df['cloneFraction_x'].notna()) & (df['cloneFraction_y'].notna())]
    
    fmin1=df['cloneFraction_x'].min()
    fmax1=df['cloneFraction_x'].max()
    fmin2=df['cloneFraction_y'].min()
    fmax2=df['cloneFraction_y'].max()
    fmin=min(fmin1,fmin2,fmin)
    fmax=max(fmax1,fmax2,fmax)
    
    df=df.apply(lambda x: x.fillna(fmin*0.7))
    ax=axes[int(i//2)][int(i%2)]
    
    plt.figure(figsize=(5,5))
            
    ax.set_xlim(fmin/2,fmax)
    ax.set_ylim(fmin/2,fmax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    sns.scatterplot(x="cloneFraction_x", y="cloneFraction_y",alpha=0.3,data=df,s=100, linewidth=0,ax=ax).set_title(str(iz)+', '+str(overlap_percent)+'% identical',fontsize=20)
    
    ax.set_xlabel('', fontsize=15)
    ax.set_ylabel('', fontsize=15)
    

axes[0][0].tick_params(labelsize=15)
axes[1][0].tick_params(labelsize=15)
axes[0][1].tick_params(labelsize=15)
axes[1][1].tick_params(labelsize=15)

axes[0][0].set_xlim(fmin/2,fmax)
axes[0][0].set_ylim(fmin/2,fmax)

axes[1][0].set_xlim(fmin/2,fmax)
axes[1][0].set_ylim(fmin/2,fmax)

axes[0][1].set_xlim(fmin/2,fmax)
axes[0][1].set_ylim(fmin/2,fmax)

axes[1][1].set_xlim(fmin/2,fmax)
axes[1][1].set_ylim(fmin/2,fmax)

fig.text(0.5, 0.04, fname1.split("_")[1]+'_'+fname1.split("_")[2]+'_'+fname1.split("_")[3], ha='center',fontsize=15)
fig.text(0.04, 0.5, fname2.split("_")[1]+'_'+fname2.split("_")[2]+'_'+fname2.split("_")[3], va='center', rotation='vertical',fontsize=15)
#fig.savefig('Pictures/Correlation_plot_'+fname1.split('.')[0]+'_vs_'+fname2.split('.')[0]+'.png')


# In[ ]:


pip install --upgrade matplotlib


# In[ ]:


overlap_clones[['targetSequences_x','targetSequences_y','cloneFraction_x','cloneFraction_y']]


# In[ ]:


df1=df1[df1['allCHitsWithScore'].str.contains('IGHA', na=False)].copy()
df2=df2[df2['allCHitsWithScore'].str.contains('IGHA', na=False)].copy()

s1=df1['cloneFraction'].sum()
s2=df2['cloneFraction'].sum()
df1['cloneFraction']=df1['cloneFraction']/s1
df2['cloneFraction']=df2['cloneFraction']/s2

if path.exists("Sosulka/")==0:
    os.makedirs("Sosulka/")
if path.exists("Sosulka/DataByIsotypeNorm/")==0:
    os.makedirs("Sosulka/DataByIsotypeNorm/")
if path.exists("Sosulka/VDJdata/")==0:
    os.makedirs("Sosulka/VDJdata/")

#добавить выуживание названий файлов и их новую сборку
df1.to_csv("Sosulka/DataByIsotypeNorm/"+'lun_lcp3_norm_bulk_1_2_2.clonotypes.IGHA.txt',sep="\t",index=False)
df2.to_csv("Sosulka/DataByIsotypeNorm/"+'mel_mp7_LN_bulk_1_1_2.clonotypes.IGHA.txt',sep="\t",index=False)


# In[ ]:


get_ipython().system('$VDJTOOLS Convert -S mixcr Sosulka/DataByIsotypeNorm/lun_lcp3_norm_bulk_1_2_2.clonotypes.IGHA.txt Sosulka/DataByIsotypeNorm/mel_mp7_LN_bulk_1_1_2.clonotypes.IGHA.txt Sosulka/VDJdata/')


# In[ ]:


df1=pd.read_csv("Sosulka/VDJdata/lun_lcp3_norm_bulk_1_2_2.clonotypes.IGHA.txt",comment="#", sep="\t", header=0)
df2=pd.read_csv("Sosulka/VDJdata/mel_mp7_LN_bulk_1_1_2.clonotypes.IGHA.txt",comment="#", sep="\t", header=0)


# In[ ]:




