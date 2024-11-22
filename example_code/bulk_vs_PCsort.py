#!/usr/bin/env python
# coding: utf-8

# In[92]:


import csv
import numpy as np
import pandas as pd
import os
from os import path
import re
from shutil import copy
from itertools import (takewhile,repeat)
import glob

from statannot import add_stat_annotation
from scipy import stats
import csv
import zipfile
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns

def rawincount(filename):
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )

VDJTOOLS='/software/bin/vdjtools'

izlist=['IGH','IGHA','IGHG','IGHM']

def pair(x,y):
    l=[x,y]
    l.sort()
    ans=l[0]+"/"+l[1]
    return ans

def pairmain(x):
    if x.find('tum')!=-1 and x.find('LN')!=-1:
        ans='LN/tum'
    elif x.find('PBMC')!=-1 and x.find('LN')!=-1:
        ans='LN/PBMC'
    elif x.find('tum')!=-1 and x.find('PBMC')!=-1:
        ans='tum/PBMC'
    else: ans="other"
    return ans


# In[93]:


#Calculate Paiwise Distances. VDJ format is required
def pwd_calculate(path_to_files,path_to_save,izlist,sort_type):
    for f in glob.glob(path_to_save+'*'+sort_type+'*'):
        os.remove(f)
    files=[f for f in os.listdir(path_to_files) if re.match('.*'+sort_type+'.*'+'.clonotypes.IGH.*.txt', f)] 
    patients=list(set([x.split('_')[1] for x in files]))
    patients.sort()
    for iz in izlist:
        for p in patients:
            get_ipython().system('$VDJTOOLS CalcPairwiseDistances $path_to_files*$p*$sort_type*$iz".txt" $path_to_save"$p"_"$sort_type"_"$iz"')
        


# In[94]:


get_ipython().run_cell_magic('capture', '', 'pwd_calculate("/home/sofyakrasik/BCR_CDR3_data/5RACE/Top_for_patient/Pooled_Replicas/VDJ/","/home/sofyakrasik/PwDistances/",[\'IGH\',\'IGHA\',\'IGHG\',\'IGHM\'],\'PCsort\')\npwd_calculate("/home/sofyakrasik/BCR_CDR3_data/5RACE/Top_for_patient/Pooled_Replicas/VDJ/","/home/sofyakrasik/PwDistances/",[\'IGH\',\'IGHA\',\'IGHG\',\'IGHM\'],\'bulk\')')


# In[95]:


df=pd.DataFrame()

for f in glob.glob('/home/sofyakrasik/PwDistances/*ccp2*.txt'):
    data=pd.read_csv(f,sep='\t')
    df=df.append(data)
    
for f in glob.glob('/home/sofyakrasik/PwDistances/*ccp3*.txt'):
    data=pd.read_csv(f,sep='\t')
    df=df.append(data)
    
df.to_csv("PwD_bulk_vs_PCsort.txt",sep='\t')


# In[96]:


df=pd.read_csv("PwD_bulk_vs_PCsort.txt",comment="#", sep="\t", header=0)

df=df[['1_sample_id','2_sample_id','div1','div12','F2','R']].copy()

df['patient']=[x.split('_')[1] for x in df['1_sample_id']]
df['isotype']=[x.split('.')[2] for x in df['1_sample_id']]
df['org1']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5] for x in df['1_sample_id']]
df['org2']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5] for x in df['2_sample_id']]
df['sort']=[x.split('_')[3] for x in df['1_sample_id']]
df['pair']=[pair(x,y) for x,y in zip(df['org1'],df['org2'])]
df['pairmain']=[pairmain(x) for x in df['pair']]
df['overlap_share']=df['div12']/df['div1']
#df=df[ df['pair'].str.contains('LN') & df['pair'].str.contains('tum') ]
df=df.sort_values('sort')

#df=df[df['pair']=='LN21/tum12']
df


# In[6]:


for iz in izlist:
    df1=df[df['isotype']==iz].copy()
    figure, axes = plt.subplots(nrows=1, ncols=3,figsize=(17, 6))
    pairs=['LN/PBMC','LN/tum','tum/PBMC']
    ord=['bulk','PCsort']
    figure.suptitle('bulk vs PCsort overlap %, paired t-test, '+str(iz),size=20)

    for i in range(len(pairs)):
        df2=df1[df1['pairmain']==pairs[i]].copy()
        df2=df2.sort_values(["pair","sort","patient"])
        df2=df2.sort_values('sort')
        sns.scatterplot(x="sort", y="overlap_share",hue="pair", alpha=0.5,data=df2,s=200,ax=axes[i]).set_title(pairs[i],fontsize=15)
        axes[i].set(xlabel='', ylabel='% overlap')
        test_results = add_stat_annotation(ax=axes[i], data=df2, x="sort", y="overlap_share", order=ord, box_pairs=[("bulk", "PCsort")],
                                           test='Mann-Whitney',text_format='simple',fontsize='x-large',loc='inside', verbose=2)
        for x in list(df2['patient'].unique()):
            for y in list(df2['pair'].unique()):
                m=df2[(df2['pair']==y) & (df2['patient']==x)]
                axes[i].plot(m['sort'],m['overlap_share'],color='grey',alpha=1,linewidth=0.5)
        axes[i].legend(loc='lower center',fontsize='small')
        df3=df2.groupby(["sort"]).agg({'overlap_share':'mean'})
        df3=df3.reset_index()
        sns.scatterplot(x="sort", y="overlap_share",alpha=1,data=df3,s=500,ax=axes[i],color='black',marker='_')
        axes[i].legend(loc='lower center',fontsize='small')
        axes[i].set_xlabel('',fontsize=15)
        axes[i].set_ylabel('',fontsize=15)
        axes[i].tick_params(labelsize=15)
        axes[i].legend(loc='lower center',fontsize='x-large', title_fontsize='40')
        plt.savefig("/home/sofyakrasik/Pipelines/Pictures/bulk_vs_PCsort_overlap_share_"+iz+".png")
        plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/bulk_vs_PCsort_overlap_share_"+iz+".pdf",bbox_inches='tight')
        


# In[7]:


for iz in izlist:
    df1=df[df['isotype']==iz].copy()
    figure, axes = plt.subplots(nrows=1, ncols=3,figsize=(17, 6))
    pairs=['LN/PBMC','LN/tum','tum/PBMC']
    ord=['bulk','PCsort']
    figure.suptitle('bulk vs PCsort, R metric, paired t-test, '+str(iz),size=20)

    for i in range(len(pairs)):
        df2=df1[df1['pairmain']==pairs[i]].copy()
        df2=df2.sort_values(["pair","sort","patient"])
        df2=df2.sort_values('sort')
        sns.scatterplot(x="sort", y="R",hue="pair", alpha=0.5,data=df2,s=200,ax=axes[i]).set_title(pairs[i],fontsize=15)
        axes[i].set(xlabel='', ylabel='R')
        test_results = add_stat_annotation(ax=axes[i], data=df2, x="sort", y="R", order=ord, box_pairs=[("bulk", "PCsort")],
                                           test='Mann-Whitney',text_format='simple',fontsize='x-large',loc='inside', verbose=2)
        for x in list(df2['patient'].unique()):
            for y in list(df2['pair'].unique()):
                m=df2[(df2['pair']==y) & (df2['patient']==x)]
                axes[i].plot(m['sort'],m['R'],color='grey',alpha=1,linewidth=0.5)
        axes[i].legend(loc='lower center',fontsize='small')
        df3=df2.groupby(["sort"]).agg({'R':'mean'})
        df3=df3.reset_index()
        sns.scatterplot(x="sort", y="R",alpha=1,data=df3,s=500,ax=axes[i],color='black',marker='_')
        axes[i].legend(loc='lower center',fontsize='small')
        axes[i].set_xlabel('',fontsize=15)
        axes[i].set_ylabel('',fontsize=15)
        axes[i].tick_params(labelsize=15)
        axes[i].legend(loc='lower center',fontsize='x-large', title_fontsize='40')
        plt.savefig("/home/sofyakrasik/Pipelines/Pictures/bulk_vs_PCsort_R_"+iz+".png")
        plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/bulk_vs_PCsort_R_"+iz+".pdf",bbox_inches='tight')


# In[97]:


test='Wilcoxon'

for iz in izlist:
    df1=df[df['isotype']==iz].copy()
    figure, axes = plt.subplots(nrows=1, ncols=3,figsize=(17, 6))
    pairs=['LN/PBMC','LN/tum','tum/PBMC']
    ord=['bulk','PCsort']
    figure.suptitle('bulk vs PCsort overlap, metric F2, '+test+', '+str(iz),size=30, y=1.01)

    for i in range(len(pairs)):
        df2=df1[df1['pairmain']==pairs[i]].copy()
        df2=df2.sort_values(["pair","sort","patient"])
        df2=df2.sort_values('sort',ascending=False)
        sns.scatterplot(x="sort", y="F2",hue="pair", alpha=0.5,data=df2,s=200,ax=axes[i]).set_title(pairs[i],fontsize=20)
        axes[i].set(xlabel='', ylabel='F2')
        test_results = add_stat_annotation(ax=axes[i], data=df2, x="sort", y="F2", order=ord, box_pairs=[("bulk", "PCsort")],
                                           test=test,text_format='simple',loc='inside',fontsize='x-large', verbose=2)
        for x in list(df2['patient'].unique()):
            for y in list(df2['pair'].unique()):
                m=df2[(df2['pair']==y) & (df2['patient']==x)]
                axes[i].plot(m['sort'],m['F2'],color='grey',alpha=1,linewidth=0.5)
        df3=df2.groupby(["sort"]).agg({'F2':'mean'})
        df3=df3.reset_index()
        sns.scatterplot(x="sort", y="F2",alpha=1,data=df3,s=500,ax=axes[i],color='black',marker='_')
        axes[i].set_xlabel('',fontsize=15)
        axes[i].set_ylabel('',fontsize=15)
        axes[i].tick_params(labelsize=15)
        axes[i].legend(loc='lower center',fontsize='x-large', title_fontsize='40')
        plt.savefig("/home/sofyakrasik/Pipelines/Pictures/bulk_vs_PCsort_F2_"+iz+".png",bbox_inches='tight')
        plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/bulk_vs_PCsort_F2_"+iz+".pdf",bbox_inches='tight')


# In[82]:


#ISOTYPE COMPOSITION
def union_dataframes_isotype_composition(path_to_files,pattern):
    df=pd.DataFrame()
    files=[f for f in os.listdir(path_to_files) if re.match(pattern,f)]
    for f in files:
        data=pd.read_csv(path_to_files+f,sep='\t')
        df=df.append(data)
    df.to_csv("isotypecomposition.txt",sep='\t', index=False)
    
union_dataframes_isotype_composition('/home/sofyakrasik/BCR_CDR3_data/5RACE/IsotypeComposition/','.*IGH.txt')


# In[83]:


df=pd.read_csv('isotypecomposition.txt',sep='\t')

df.columns=['filename', 'isotype', 'total_frequency']

df['patient']=[x.split("_")[1] for x in df['filename']]
df['organ']=[x.split("_")[2] for x in df['filename']]
df['sort_type']=[x.split("_")[3] for x in df['filename']]
df['sample']=[x.split("_")[4]+x.split("_")[5] for x in df['filename']]
df['replicate']=[x.split("_")[6].split(".")[0] for x in df['filename']]


df=df[['patient','organ','isotype','sample','replicate','sort_type','total_frequency']].copy()

#pool replicas
df['total_frequency']=df.groupby(['patient','organ','isotype','sample','sort_type'])['total_frequency'].transform('mean')
df=df.drop_duplicates(['patient','organ','isotype','sample','sort_type','total_frequency'])


#make sure both sort types are available for the patient
df['count_sort_type']=df.groupby(['patient'])['sort_type'].transform('nunique')
df=df[df['patient'].isin(['ccp2','ccp3'])]

#organs=list(df['organ'].unique())
organs=['LN','tum','PBMC']
patients=list(df['patient'].unique())


# In[85]:


test='Mann-Whitney'
izlist=['IGHA','IGHG','IGHM']

for iz in izlist:
    figure,axes=plt.subplots(nrows=1,ncols=3,figsize=(18,6))
    figure.suptitle('Isotype composition, '+test+', '+iz,fontsize=25,y=1.01)
    
    df1=df[df['isotype']==iz].copy()
    for i in range(len(organs)):
        data=df1[df1['organ']==organs[i]].copy()
        data=data.sort_values(['sort_type','patient','isotype'],ascending=False)
        sns.stripplot(ax=axes[i],x="sort_type",y='total_frequency',data=data,hue='patient',size=15,alpha=0.8).set_title(organs[i])
        sns.boxplot(x="sort_type", y="total_frequency", color='white',data=data,ax=axes[i]).set_title(organs[i],fontsize=20)
        
        data=data.sort_values(['sort_type','patient','isotype'])
        test_results = add_stat_annotation(ax=axes[i], data=data, x="sort_type", y="total_frequency", order=ord, box_pairs=[("bulk", "PCsort")],test=test,text_format='simple',loc='inside',fontsize='x-large',verbose=2,comparisons_correction=None)
        axes[i].set_xlabel('',fontsize=15)
        axes[i].set_ylabel('',fontsize=15)
        axes[i].tick_params(labelsize=15)
        axes[i].legend(fontsize='x-large', title_fontsize='40')
    plt.savefig("/home/sofyakrasik/Pipelines/Pictures/bulk_vs_PCsort_isotypes_"+iz+".png",bbox_inches='tight')
    plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/bulk_vs_PCsort_isotypes_"+iz+".pdf",bbox_inches='tight')
    


# In[86]:


#CDR3 length
def union_dataframes(path_to_files,pattern):
    df=pd.DataFrame()
    files=[f for f in os.listdir(path_to_files) if re.match(pattern,f)]
    for f in files:
        data=pd.read_csv(path_to_files+f,sep='\t')
        data['filename']=f
        df=df.append(data)
    df['patient']=[x.split('_')[1] for x in df['filename']]
    df['organ']=[x.split('_')[2] for x in df['filename']]
    df['sort_type']=[x.split('_')[3] for x in df['filename']]
    df['sample']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5]+x.split('_')[6].split('.')[0] for x in df['filename']]
    df=df[df['allCHitsWithScore'].notna()]
    df['isotype']=[x[0:4] for x in df['allCHitsWithScore']]
    return df
    
df_PCsort=union_dataframes('/home/sofyakrasik/BCR_CDR3_data/5RACE/All_data_filtered/','.*ccp[23].*PCsort.*IGH.txt')
df_bulk=union_dataframes('/home/sofyakrasik/BCR_CDR3_data/5RACE/All_data_filtered/','.*ccp[23].*bulk.*IGH.txt')

df=df_PCsort.append(df_bulk)


# In[87]:


test='Mann-Whitney'
izlist=['IGH','IGHA','IGHG','IGHM']
parameter='cdr3_length'

for iz in izlist:
    figure,axes=plt.subplots(nrows=1,ncols=3,figsize=(18,6))
    figure.suptitle(parameter+', '+test+', '+iz,fontsize=25,y=1.01)
    
    df1=df[df['isotype'].str.contains(iz)].copy()
    top_min=min(df1.groupby(['filename'])['targetSequences'].transform('count').min(),50)
    df1=df1.groupby(['filename']).head(top_min)
    df1['cloneFraction']=df1.groupby(['filename'])['cloneFraction'].transform(lambda x: x/x.sum())

    df1['cdr3_length']=[len(x) for x in df1['aaSeqCDR3']]
    df1['cdr3_length_wt']=[len(x)*y for x,y in zip(df1['aaSeqCDR3'],df1['cloneFraction'])]
    df1=df1.groupby(['filename','patient','sort_type','organ']).agg({parameter:'mean'}).reset_index()
    for i in range(len(organs)):
        data=df1[df1['organ']==organs[i]].copy()
        data=data.sort_values(['sort_type','patient','organ'],ascending=False)
        sns.stripplot(ax=axes[i],x="sort_type",y=parameter,data=data,hue='patient',size=15,alpha=0.8,palette="Set2").set_title(organs[i])
        sns.boxplot(x="sort_type", y=parameter, color='white',data=data,ax=axes[i]).set_title(organs[i],fontsize=20)
        data=data.sort_values(['sort_type','patient','organ'])
        test_results = add_stat_annotation(ax=axes[i], data=data, x="sort_type", y=parameter, order=ord, box_pairs=[("bulk", "PCsort")],test=test,text_format='simple',loc='inside',fontsize='x-large',verbose=2,comparisons_correction=None)
        axes[i].set_xlabel('',fontsize=15)
        axes[i].set_ylabel('',fontsize=15)
        axes[i].tick_params(labelsize=15)
        axes[i].legend(fontsize='x-large', title_fontsize='40')
    plt.savefig("/home/sofyakrasik/Pipelines/Pictures/bulk_vs_PCsort_"+parameter+"_"+iz+".png",bbox_inches='tight')
    plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/bulk_vs_PCsort_"+parameter+"_"+iz+".pdf",bbox_inches='tight')


# In[73]:


#df=df[df['isotype']=='IGH'].copy()
#df['cdr3_length']=[len(x) for x in df['aaSeqCDR3']]
#df=df.groupby(['patient','organ','sample','sort_type']).agg({'cdr3_length':'mean'}).reset_index()
figure,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,6))
figure.suptitle('',fontsize=25)
sns.stripplot(ax=axes,x="sort_type",y='cdr3_length',data=df,hue='patient',size=15,alpha=0.8,palette="Set2")
sns.boxplot(ax=axes,x="sort_type",y='cdr3_length',data=df,color='white')
test_results = add_stat_annotation(ax=axes, data=df, x="sort_type", y='cdr3_length', order=["bulk", "PCsort"], box_pairs=[("bulk", "PCsort")],test=test,text_format='simple',loc='inside',fontsize='x-large',verbose=2,comparisons_correction=None)


# In[88]:


get_ipython().run_cell_magic('capture', '', '#CALCULATE DIVERSITY STATS\n#CLONALITY\nif path.exists("DiversityStats/")==0:\n    os.makedirs("DiversityStats/")\n    \n!rm DiversityStats/*\n\nfor iz in izlist:\n    !$VDJTOOLS CalcDiversityStats -i aa /home/sofyakrasik/BCR_CDR3_data/5RACE/Top_for_patient/Separate_Replicas/VDJ/*"ccp2"*"PCsort"*"$iz".txt ./DiversityStats/DiversityStats_"ccp2_PCsort"_$iz \n    !$VDJTOOLS CalcDiversityStats -i aa /home/sofyakrasik/BCR_CDR3_data/5RACE/Top_for_patient/Separate_Replicas/VDJ/*"ccp2"*"bulk"*"$iz".txt ./DiversityStats/DiversityStats_"ccp2_bulk"_$iz \n    !$VDJTOOLS CalcDiversityStats -i aa /home/sofyakrasik/BCR_CDR3_data/5RACE/Top_for_patient/Separate_Replicas/VDJ/*"ccp3"*"PCsort"*"$iz".txt ./DiversityStats/DiversityStats_"ccp3_PCsort"_$iz \n    !$VDJTOOLS CalcDiversityStats -i aa /home/sofyakrasik/BCR_CDR3_data/5RACE/Top_for_patient/Separate_Replicas/VDJ/*"ccp3"*"bulk"*"$iz".txt ./DiversityStats/DiversityStats_"ccp3_bulk"_$iz ')


# In[89]:


#CLONALITY
def union_dataframes(path_to_files,pattern):
    df=pd.DataFrame()
    files=[f for f in os.listdir(path_to_files) if re.match(pattern,f)]
    for f in files:
        data=pd.read_csv(path_to_files+f,sep='\t')
        df=df.append(data)
    df['patient']=[x.split('_')[1] for x in df['sample_id']]
    df['organ']=[x.split('_')[2] for x in df['sample_id']]
    df['sort_type']=[x.split('_')[3] for x in df['sample_id']]
    df['sample']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5]+x.split('_')[6].split('.')[0] for x in df['sample_id']]
    df['isotype']=[x.split('.')[-1] for x in df['sample_id']]
    return df

df=union_dataframes('/home/sofyakrasik/Pipelines/DiversityStats/','.*.diversity.aa.resampled.txt')
data=df.copy()
data['clonality']=[1-x for x in data['normalizedShannonWienerIndex_mean']]


# In[91]:


test='Mann-Whitney'
izlist=['IGH','IGHA','IGHG','IGHM']
organs=['LN','tum','PBMC']

for iz in izlist:
    df=data[data['isotype']==iz].copy()
    figure,axes=plt.subplots(nrows=1,ncols=3,figsize=(18,6))
    figure.suptitle('Clonality'+', '+test+', '+iz,fontsize=25,y=1.01)
    for i in range(len(organs)):
        df1=df[df['organ']==organs[i]].copy()
        
        df1=df1.sort_values(['sort_type','patient','organ'],ascending=False)
        sns.stripplot(ax=axes[i],x="sort_type",y='clonality',data=df1,hue='patient',size=15,alpha=0.8,palette="Set1").set_title(organs[i],fontsize=20)
        sns.boxplot(x="sort_type", y='clonality', color='white',data=df1,ax=axes[i]).set_title(organs[i],fontsize=20)
        df1=df1.sort_values(['sort_type','patient','organ'])
        test_results = add_stat_annotation(ax=axes[i], data=df1, x="sort_type", y='clonality', order=ord, box_pairs=[("bulk", "PCsort")],test=test,text_format='simple',loc='inside',fontsize='x-large',verbose=2,comparisons_correction=None)
        axes[i].set_xlabel('',fontsize=15)
        axes[i].set_ylabel('',fontsize=15)
        axes[i].tick_params(labelsize=15)
        axes[i].legend(fontsize='x-large', title_fontsize='40')
    plt.savefig("/home/sofyakrasik/Pipelines/Pictures/bulk_vs_PCsort_clonality_"+iz+".png",bbox_inches='tight')
    plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/bulk_vs_PCsort_clonality_"+iz+".pdf",bbox_inches='tight')


# In[ ]:




