#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#ИЗМЕНИТЬ НОРМИРОВКУ ПЕРЕД СЛИВОМ РЕПЛИК И ОБРАЗЦОВ, ВЫБРАСЫВАТЬ СИЛЬНО МАЛЕНЬКИЕ ОБРАЗЦЫ?


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

def rawincount(filename):
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )

VDJTOOLS='/software/bin/vdjtools'


# In[7]:


#pat='lcp3'
sort_type='bulk'
path_to_files='/home/sofyakrasik/BCR_CDR3_data/5RACE/All_data_filtered/'

izlist=['IGH','IGHA','IGHG','IGHM']

org1="tum"
org2="LN"

files=[f for f in os.listdir(path_to_files) if re.match('.*'+sort_type+'.*'+'IGH.txt', f)]
files.sort()
#files.remove('mel_mp2_LN_bulk_1_1_1.clonotypes.IGH.txt')


# In[8]:


if path.exists("DataByIsotype/")==0:
    os.makedirs("DataByIsotype/")

get_ipython().system('rm DataByIsotype/*')
    
cols=['cloneFraction','allCHitsWithScore','aaSeqCDR3','filename']

mixcr_data = pd.DataFrame(columns=cols)

for f in files:
    df=pd.read_csv(path_to_files+f,comment="#", sep="\t", header=0)
    df=df.dropna(subset=['allCHitsWithScore'])
    for iz in izlist:
        df_iz=df[df['allCHitsWithScore'].str.contains(iz)].copy()
        df_iz['cloneFraction']=df_iz['cloneFraction']/df_iz['cloneFraction'].sum()
        fname=f.split('.')[0]+'.'+f.split('.')[1]+'.'+iz+'.txt'
        df_iz.to_csv('DataByIsotype/'+fname,sep='\t',index=None)


# In[9]:


get_ipython().run_cell_magic('capture', '', '\nif path.exists("VDJdata/")==0:\n    os.makedirs("VDJdata/")\n\n!rm VDJdata/*\n\nfor iz in izlist:\n   !$VDJTOOLS Convert -S mixcr DataByIsotype/*$sort_type*${iz}.txt VDJdata/')


# In[10]:


#POOLING REPLICAS

if path.exists("PooledReplicas/")==0:
    os.makedirs("PooledReplicas/")
    
get_ipython().system('rm PooledReplicas/*')
    
for can in list(set([f.split('_')[0] for f in os.listdir('VDJdata/') if re.match('.*'+sort_type+'.*'+'IGH'+'.*', f)])):
    for pat in list(set([f.split('_')[1] for f in os.listdir('VDJdata/') if re.match(can+'.*'+sort_type+'.*'+'IGH'+'.*', f)])):
        for org in list(set([f.split('_')[2] for f in os.listdir('VDJdata/') if re.match(can+'.*'+pat+'.*'+sort_type+'.*'+'IGH'+'.*', f)])):
            for iz in list(set([f.split('.')[2] for f in os.listdir('VDJdata/') if re.match(can+'.*'+pat+'.*'+org+'.*'+sort_type+'.*'+'IGH'+'.*', f)])):
                for smpl in list(set(['_'+f.split('_')[4]+'_'+f.split('_')[5]+'_' for f in os.listdir('VDJdata/') if re.match(can+'.*'+pat+'.*'+org+'.*'+sort_type+'.*'+iz+'.*', f)])):
                    l=[f for f in os.listdir('VDJdata/') if re.match(can+'.*'+pat+'.*'+org+'.*'+sort_type+'.*'+smpl+'.*'+iz+'.txt', f)]
                    cnt=len(l)
                    if cnt==1:
                        #shutil.copy("VDJdata/"+l[0],"PooledReplicas/"+l[0])
                        df=pd.read_csv('VDJdata/'+l[0],comment="#", sep="\t", header=0)
                        df['freq']=df.groupby('cdr3aa')['freq'].transform('sum')
                        df=df.drop_duplicates('cdr3aa')
                        df.to_csv("PooledReplicas/"+l[0],sep='\t',index=None)
                    else:
                        top_max=rawincount('VDJdata/'+l[0])-1
                        
                        for f in l:
                            c=rawincount('VDJdata/'+f)
                            if c<top_max+1:
                                top_max=c-1
                        #print(pat,org,smpl,top_max)
                        
                        cols=['count', 'freq', 'cdr3nt', 'cdr3aa', 'v', 'd', 'j', 'VEnd', 'DStart','DEnd', 'JStart']
                        df = pd.DataFrame(columns=cols)

                        for f in l:
                            data=pd.read_csv('VDJdata/'+l[0],comment="#", sep="\t", header=0)
                            data=data.sort_values('freq',ascending=False)
                            #Uncomment to make normalization work
                            #data=data.head(top_max).copy() 
                            data['freq']=data['freq']/data['freq'].sum()
                            df=df.append(data)
                        df['freq']=df.groupby('cdr3aa')['freq'].transform('sum')/cnt
                        df=df.drop_duplicates('cdr3aa')
                        l.sort()
                        df.to_csv("PooledReplicas/"+l[0],sep='\t',index=None)


# In[11]:


#POOLING SAMPLES
if path.exists("PooledSamples/")==0:
    os.makedirs("PooledSamples/")

get_ipython().system('rm PooledSamples/*')
    
for can in list(set([f.split('_')[0] for f in os.listdir('PooledReplicas/') if re.match('.*'+sort_type+'.*'+'IGH'+'.*', f)])):
    for pat in list(set([f.split('_')[1] for f in os.listdir('PooledReplicas/') if re.match(can+'.*'+sort_type+'.*'+'IGH'+'.*', f)])):
        for org in list(set([f.split('_')[2] for f in os.listdir('PooledReplicas/') if re.match(can+'.*'+pat+'.*'+sort_type+'.*'+'IGH'+'.*', f)])):
            for iz in list(set([f.split('.')[2] for f in os.listdir('PooledReplicas/') if re.match(can+'.*'+pat+'.*'+org+'.*'+sort_type+'.*'+'IGH'+'.*', f)])):
                l=[f for f in os.listdir('PooledReplicas/') if re.match(can+'.*'+pat+'.*'+org+'.*'+sort_type+'.*'+iz+'.txt', f)]
                cnt=len(l)
                if cnt==1:
                    shutil.copy("PooledReplicas/"+l[0],"PooledSamples/"+l[0])
                else:
                    top_max=rawincount('PooledReplicas/'+l[0])-1
                        
                    for f in l:
                        c=rawincount('PooledReplicas/'+f)
                        if c<top_max+1:
                            top_max=c-1
                    #print(pat,org,smpl,top_max)
                        
                    cols=['count', 'freq', 'cdr3nt', 'cdr3aa', 'v', 'd', 'j', 'VEnd', 'DStart','DEnd', 'JStart']
                    df = pd.DataFrame(columns=cols)

                    for f in l:
                        data=pd.read_csv('PooledReplicas/'+l[0],comment="#", sep="\t", header=0)
                        data=data.sort_values('freq',ascending=False)
                        #Uncomment to make normalization work
                        #data=data.head(top_max).copy()
                        data['freq']=data['freq']/data['freq'].sum()
                        df=df.append(data)
                    df['freq']=df.groupby('cdr3aa')['freq'].transform('sum')/cnt
                    df=df.drop_duplicates('cdr3aa')
                    l.sort()
                    df.to_csv("PooledSamples/"+l[0],sep='\t',index=None)


# In[12]:


#DELETE POOLED SAMPLES WITH LESS THAN 50 CLONOTYPES
files=[f for f in os.listdir("PooledSamples/") if re.match('.*'+sort_type+'.*'+'IGH.*txt', f)]
files2=files.copy()

for f in files2:
    if f.find("norm")!=-1:
        files.remove(f)


for f in files:
    df=pd.read_csv("PooledSamples/"+f,sep='\t')
    if len(df)<51:
        os.remove(os.path.join("PooledSamples/", f))


# In[13]:


get_ipython().run_cell_magic('capture', '', '#CALCULATE DIVERSITY STATS\nif path.exists("DiversityStats/")==0:\n    os.makedirs("DiversityStats/")\n    \n!rm DiversityStats/*\n\nfor iz in izlist:\n    !$VDJTOOLS CalcDiversityStats -i aa PooledSamples/*$sort_type*"$iz".txt ./DiversityStats/DiversityStats_"$org1"_$iz \n    #!$VDJTOOLS CalcDiversityStats -i aa PooledSamples/*"tum"*$sort_type*"$iz".txt ./DiversityStats/DiversityStats_"$org2"_$iz\n    #!$VDJTOOLS CalcDiversityStats -i aa PooledSamples/*"LN"*$sort_type*"$iz".txt ./DiversityStats/DiversityStats_LN_$iz\n    #!$VDJTOOLS CalcDiversityStats -i aa PooledSamples/*"norm"*$sort_type*"$iz".txt ./DiversityStats/DiversityStats_norm_$iz')


# In[21]:


files=[f for f in os.listdir("DiversityStats/") if re.match('.*'+'.diversity.aa.resampled.txt', f)]
df = pd.DataFrame()

for f in files:
    data=pd.read_csv("DiversityStats/"+f,comment="#", sep="\t", header=0)
    df=df.append(data)
    
df['cancer']=[x.split('_')[0] for x in df['sample_id']]
df['patient']=[x.split('_')[1] for x in df['sample_id']]
df['organ']=[x.split('_')[2] for x in df['sample_id']]
df['isotype']=[x.split('.')[-1] for x in df['sample_id']]
df=df.sort_values(['cancer','patient','organ','isotype'],ascending=True)
df=df[df['organ']!='norm'].copy()
data=df.copy()
data['organ_cnt']=data.groupby(['patient','isotype'])['organ'].transform('count')
#data=data[data['organ_cnt']>=3]


# In[15]:


for iz in izlist:
    files=[f for f in os.listdir("DiversityStats/") if re.match('.*'+iz+'.diversity.aa.resampled.txt', f)]
    df_1=pd.read_csv('DiversityStats/'+files[0],comment="#", sep="\t", header=0)
    df_1=df_1[['sample_id','normalizedShannonWienerIndex_mean']]
    df_1['patient']=[x.split('_')[1] for x in df_1['sample_id']]
    files=[f for f in os.listdir("DiversityStats/") if re.match('.*'+iz+'.diversity.aa.resampled.txt', f)]
    df_2=pd.read_csv('DiversityStats/'+files[0],comment="#", sep="\t", header=0)
    df_2=df_2[['sample_id','normalizedShannonWienerIndex_mean']]
    df_2['patient']=[x.split('_')[1] for x in df_2['sample_id']]
    
    df_1['organ']=[x.split('_')[2] for x in df_1['sample_id']]
    df_1=df_1[df_1['organ']==org1].copy()
    
    df_2['organ']=[x.split('_')[2] for x in df_2['sample_id']]
    df_2=df_2[df_2['organ']==org2].copy()
    
    df=pd.merge(df_2,df_1,how='inner',on='patient')
    df['cancer']=[x.split('_')[0] for x in df['sample_id_x']]
    df=df[['patient','normalizedShannonWienerIndex_mean_x','normalizedShannonWienerIndex_mean_y','cancer']].copy()
    df.columns=['patient','nShW_2','nShW_1','cancer']
    #df=df[df['patient']!='lcp1']
    
    figure, axes = plt.subplots(nrows=1, ncols=1,figsize=(6, 6))
    p=sns.regplot(x="nShW_2", y="nShW_1",data=df,ci=95,color='grey')
    
    corr,pval=pearsonr(df['nShW_2'],df['nShW_1'])
    sns.scatterplot(x="nShW_2", y="nShW_1",hue='cancer',edgecolors="black",linewidth=1,alpha=0.9,data=df,s=100).set_title(iz+", Normalized Shannon-Wiener "+org1+" vs "+org2+","+'\n'
                                                                                                                              +"corr="+str(round(corr,2))+", p-value="+str(round(pval,2))+", 95% ci",fontsize=15)
    #plt.xlabel(org2,fontsize=15) 
    #plt.ylabel(org1,fontsize=15)
    
    axes.set_xlabel(org2,fontsize=15)
    axes.set_ylabel(org1,fontsize=15)
    axes.tick_params(labelsize=15)
    axes.legend(fontsize='x-large')
    
    plt.savefig("Pictures/nShW_"+org2+"_"+org1+"_correlation_"+iz+".png")
    plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/nShW_"+org2+"_"+org1+"_correlation_"+iz+".pdf",bbox_inches='tight')
    
    plt.show()


# In[16]:


files=[f for f in os.listdir("DiversityStats/") if re.match('.*.diversity.aa.resampled.txt', f)]
#df_1=pd.read_csv('DiversityStats/'+files[0],comment="#", sep="\t", header=0)
df_1=pd.DataFrame()
for f in files:
    df=pd.read_csv('DiversityStats/'+f,comment="#", sep="\t", header=0)
    df_1=df_1.append(df)
df_1=df_1[['sample_id','normalizedShannonWienerIndex_mean']]
df_1['patient']=[x.split('_')[1] for x in df_1['sample_id']]
#files=[f for f in os.listdir("DiversityStats/") if re.match('.*.diversity.aa.resampled.txt', f)]
#df_2=pd.read_csv('DiversityStats/'+files[0],comment="#", sep="\t", header=0)
df_2=df_1.copy()
df_2=df_2[['sample_id','normalizedShannonWienerIndex_mean']]
df_2['patient']=[x.split('_')[1] for x in df_2['sample_id']]
    
df_1['organ']=[x.split('_')[2] for x in df_1['sample_id']]
df_1['isotype']=[x.split('.')[2] for x in df_1['sample_id']]
df_1=df_1[df_1['organ']==org1].copy()
    
df_2['organ']=[x.split('_')[2] for x in df_2['sample_id']]
df_2['isotype']=[x.split('.')[2] for x in df_2['sample_id']]
df_2=df_2[df_2['organ']==org2].copy()
    
df=pd.merge(df_2,df_1,how='inner',on=['patient','isotype'])
df['cancer']=[x.split('_')[0] for x in df['sample_id_x']]
df=df[['patient','normalizedShannonWienerIndex_mean_x','normalizedShannonWienerIndex_mean_y','cancer','isotype']].copy()
df.columns=['patient','nShW_2','nShW_1','cancer','isotype']
    #df=df[df['patient']!='lcp1']
    
figure, axes = plt.subplots(nrows=1, ncols=1,figsize=(6, 6))
p=sns.regplot(x="nShW_2", y="nShW_1",data=df,ci=95,color='grey')
    
corr,pval=pearsonr(df['nShW_2'],df['nShW_1'])
sns.scatterplot(x="nShW_2", y="nShW_1",hue='isotype',edgecolors="black",linewidth=1,alpha=0.9,data=df,s=100).set_title(iz+", Normalized Shannon-Wiener "+org1+" vs "+org2+","+'\n'
                                                                                                                              +"corr="+str(round(corr,2))+", p-value="+str(round(pval,2))+", 95% ci",fontsize=15)
    #plt.xlabel(org2,fontsize=15) 
    #plt.ylabel(org1,fontsize=15)
    
axes.set_xlabel(org2,fontsize=15)
axes.set_ylabel(org1,fontsize=15)
axes.tick_params(labelsize=15)
axes.legend(fontsize='x-large')
    
plt.savefig("Pictures/nShW_"+org2+"_"+org1+"_correlation_all_isotypes.png")
plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/nShW_"+org2+"_"+org1+"_correlation_all_isotypes.pdf",bbox_inches='tight')

plt.show()


# In[22]:


ord=list(df["organ"].unique())
box_pairs=list(itertools.combinations(ord, 2))

df=data[data['isotype']!='IGH']
df['normalizedShannonWienerIndex_mean']=[1-x for x in df['normalizedShannonWienerIndex_mean']]  
test='Mann-Whitney'


figure, axes = plt.subplots(nrows=1, ncols=1,figsize=(6, 6))
sns.scatterplot(x="organ", y="normalizedShannonWienerIndex_mean",hue='isotype',edgecolors="black",linewidth=1,alpha=0.9,data=df,s=100).set_title(iz,fontsize=15)
sns.boxplot(x="organ", y="normalizedShannonWienerIndex_mean",color="white",linewidth=1,data=df).set_title('Clonality in different tissues,\n'+test+', Bonferroni cor.',fontsize=15)
    
test_results = add_stat_annotation(ax=axes, data=df, x="organ", y="normalizedShannonWienerIndex_mean", order=ord, box_pairs=box_pairs,test=test,text_format='simple',loc='inside',fontsize='x-large', verbose=2)
axes.set_xlabel('',fontsize=15)
axes.set_ylabel('Clonality index',fontsize=15)
axes.tick_params(labelsize=15)
axes.legend(fontsize='x-large')
axes.set_ylim(0,1.2)

plt.savefig("Pictures/Clonality_index_1_nShW_in_tissues_all_isotypes"".png",bbox_inches='tight')
plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/Clonality_index_1_nShW_in_tissues_all_isotypes"".pdf",bbox_inches='tight')


# In[23]:


ord=list(df["organ"].unique())
box_pairs=list(itertools.combinations(ord, 2))

test='Mann-Whitney'

for iz in izlist:
    df=data[data['isotype']==iz].copy()
    #df['normalizedShannonWienerIndex_mean']=[1-x+np.random.uniform(0,0.05) for x in df['normalizedShannonWienerIndex_mean']]
    df['normalizedShannonWienerIndex_mean']=[1-x for x in df['normalizedShannonWienerIndex_mean']]
    figure, axes = plt.subplots(nrows=1, ncols=1,figsize=(6, 6))
    
    sns.scatterplot(x="organ", y="normalizedShannonWienerIndex_mean",hue='cancer',edgecolors="black",linewidth=1,alpha=0.9,data=df,s=100).set_title(iz,fontsize=15)
    sns.boxplot(x="organ", y="normalizedShannonWienerIndex_mean",color="white",linewidth=1,data=df).set_title(iz+', Clonality in different tissues,\n'+test+', Bonferroni cor.',fontsize=15)
    axes.set_ylabel('Clonality index',fontsize=15)
    axes.set_xlabel('',fontsize=15)
    axes.tick_params(labelsize=15)
    axes.legend(fontsize='x-large')
    axes.set_ylim(0,1.2)
    
    test_results = add_stat_annotation(ax=axes, data=df, x="organ", y="normalizedShannonWienerIndex_mean", order=ord, box_pairs=box_pairs,test=test,text_format='simple',loc='inside',fontsize='x-large', verbose=2)
    plt.savefig("Pictures/Clonality_index_1_nShW_in_tissues_"+iz+".png",bbox_inches='tight')
    plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/Clonality_index_1_nShW_in_tissues_"+iz+".pdf",bbox_inches='tight')


# In[ ]:




