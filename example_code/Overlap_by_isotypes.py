#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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

def rawincount(filename):
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )

def org_pair(org1,org2):
    l=[org1,org2]
    l.sort()
    ans=l[0]+'/'+l[1]
    return ans

VDJTOOLS='/software/bin/vdjtools'


# In[5]:


sort_type='bulk'
path_to_files='/home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Replicas_aa/VDJ/'

izlist=['IGHA','IGHG']

top_max=100
top_min=101

files=[f for f in os.listdir(path_to_files) if re.match('.*'+sort_type+'.*'+'IGH'+'.*'+'.txt', f)]
files.sort()


# In[6]:


get_ipython().run_cell_magic('capture', '', '#NORMILIZING SAMPLES: SAME FOR ALL ISOTYPES\n\nif path.exists("PooledReplicas_Top/")==0:\n    os.makedirs("PooledReplicas_Top/")\n\nfor f in files:\n    c=rawincount(path_to_files+f)\n    if top_min<=c<top_max+1:\n        top_max=c-1\n\nprint(top_max)\n\n!rm PooledReplicas_Top/*\n        \nfor f in files:\n    c=rawincount(path_to_files+f)\n    if c>=top_max:\n        !$VDJTOOLS SelectTop  -x $top_max $path_to_files$f PooledReplicas_Top/')


# In[7]:


get_ipython().run_cell_magic('capture', '', '#CALCULATING PAIRWISE DISTANCES (EACH ISOTYPE, CANCER, PATIENT, SEPARATELY)\nif path.exists("PwDistances/")==0:\n    os.makedirs("PwDistances/")\n\n!rm PwDistances/*\n    \nfor can in list(set([f.split(\'_\')[0] for f in os.listdir(\'PooledReplicas_Top/\') if re.match(\'.*\'+sort_type+\'.*\'+\'IGH\'+\'.*\', f)])):\n    for pat in list(set([f.split(\'_\')[1] for f in os.listdir(\'PooledReplicas_Top/\') if re.match(can+\'.*\'+sort_type+\'.*\'+\'IGH\'+\'.*\', f)])):\n            for iz in list(set([f.split(\'.\')[2] for f in os.listdir(\'PooledReplicas_Top/\') if re.match(can+\'.*\'+pat+\'.*\'+sort_type+\'.*\'+\'IGH\'+\'.*\', f)])):\n                files=[f for f in os.listdir(\'PooledReplicas_Top/\') if re.match(can+\'.*\'+pat+\'.*\'+sort_type+\'.*\'+iz+\'.txt\', f)]\n                if len(files)>=3:\n                    !$VDJTOOLS CalcPairwiseDistances PooledReplicas_Top/$can*$pat*$sort_type*"$iz".txt PwDistances/"$can"_"$pat"_"$sort_type"_"$iz".txt')


# In[8]:


get_ipython().run_cell_magic('capture', '', '#TEMPORARY CELL\n#CALCULATING PAIRWISE DISTANCES (EACH ISOTYPE, CANCER, PATIENT, SEPARATELY)\nif path.exists("PwDistances/")==0:\n    os.makedirs("PwDistances/")\n\n!rm PwDistances/*\n    \nfor can in list(set([f.split(\'_\')[0] for f in os.listdir(\'PooledSamples_Top/\') if re.match(\'.*\'+sort_type+\'.*\'+\'IGH\'+\'.*\', f)])):\n    for pat in list(set([f.split(\'_\')[1] for f in os.listdir(\'PooledSamples_Top/\') if re.match(can+\'.*\'+sort_type+\'.*\'+\'IGH\'+\'.*\', f)])):\n            for iz in list(set([f.split(\'.\')[2] for f in os.listdir(\'PooledSamples_Top/\') if re.match(can+\'.*\'+pat+\'.*\'+sort_type+\'.*\'+\'IGH\'+\'.*\', f)])):\n                files=[f for f in os.listdir(\'PooledSamples_Top/\') if re.match(can+\'.*\'+pat+\'.*\'+sort_type+\'.*\'+iz+\'.txt\', f)]\n                if len(files)>=2:\n                    !$VDJTOOLS CalcPairwiseDistances PooledSamples_Top/$can*$pat*$sort_type*"$iz".txt PwDistances/"$can"_"$pat"_"$sort_type"_"$iz".txt')


# In[9]:


get_ipython().run_cell_magic('capture', '', '#TEMPORARY CELL 2\n#CALCULATING PAIRWISE DISTANCES (EACH ISOTYPE SEPARATELY)\nif path.exists("PwDistances/")==0:\n    os.makedirs("PwDistances/")\n\n!rm PwDistances/*\n    \nfor iz in list(set([f.split(\'.\')[2] for f in os.listdir(\'PooledSamples_Top/\') if re.match(\'.*\'+sort_type+\'.*\'+\'IGH\'+\'.*\', f)])):\n    files=[f for f in os.listdir(\'PooledSamples_Top/\') if re.match(\'.*\'+sort_type+\'.*\'+iz+\'.txt\', f)]\n    if len(files)>=2:\n        !$VDJTOOLS CalcPairwiseDistances PooledSamples_Top/*$sort_type*"$iz".txt PwDistances/"$sort_type"_"$iz".txt')


# In[10]:


cols=['1_sample_id', '2_sample_id', 'div1', 'div2', 'div12', 'div21',
       'count1', 'count2', 'count12', 'count21', 'freq1', 'freq2', 'freq12',
       'freq21', 'R', 'Rs', 'D', 'F', 'F2']

data = pd.DataFrame(columns=cols)

files=[f for f in os.listdir('PwDistances/') if re.match('.*'+sort_type+'.*'+'IGH'+'.*', f)]

for f in files:
    df=pd.read_csv('PwDistances/'+f,comment="#", sep="\t", header=0)
    df=df[['1_sample_id', '2_sample_id', 'div1', 'div2', 'div12', 'div21','count1', 'count2', 'count12', 'count21', 'freq1', 'freq2', 'freq12','freq21', 'R', 'Rs', 'D', 'F', 'F2']].copy()
    data=data.append(df)


# In[11]:


#TEMPORARY CELL 2
cols=['1_sample_id', '2_sample_id', 'div1', 'div2', 'div12', 'div21',
       'count1', 'count2', 'count12', 'count21', 'freq1', 'freq2', 'freq12',
       'freq21', 'R', 'Rs', 'D', 'F', 'F2']

data = pd.DataFrame(columns=cols)

files=[f for f in os.listdir('PwDistances/') if re.match('.*'+sort_type+'.*'+'IGH'+'.*', f)]

for f in files:
    df=pd.read_csv('PwDistances/'+f,comment="#", sep="\t", header=0)
    df=df[['1_sample_id', '2_sample_id', 'div1', 'div2', 'div12', 'div21','count1', 'count2', 'count12', 'count21', 'freq1', 'freq2', 'freq12','freq21', 'R', 'Rs', 'D', 'F', 'F2']].copy()
    data=data.append(df)
    
data['patient_1']=[x.split('_')[1] for x in data['1_sample_id']]
data['patient_2']=[x.split('_')[1] for x in data['2_sample_id']]
data=data[data['patient_1']==data['patient_2']]


# In[12]:


data['cancer']=[x.split('_')[0] for x in data['1_sample_id']]
data['patient']=[x.split('_')[1] for x in data['1_sample_id']]
data['org1']=[x.split('_')[2] for x in data['1_sample_id']]
data['org2']=[x.split('_')[2] for x in data['2_sample_id']]
data['pair']=[org_pair(x,y) for x,y in zip(data['org1'],data['org2'])]
pairlist=['LN/tum','PBMC/tum','LN/PBMC']
data=data[data['pair'].isin(pairlist)]
data['isotype']=[x.split('.')[2] for x in data['1_sample_id']]
data=data[data['isotype'].isin(izlist)]
data['isotype']=['Ig'+x[-1] for x in data['isotype']]
data=data.sort_values(['patient','pair','isotype'])
data['count']=data.groupby(['patient','isotype'])['pair'].transform('nunique').copy()
#data=data[data['count']==3]


# In[16]:


test='Mann-Whitney'

data=data.sort_values(['isotype','cancer'])

for pr in pairlist:
    df=data[data['pair']==pr].copy()
    
    print('Number of patient 1',len(df['patient_2'].unique()))
    
    figure, axes = plt.subplots(nrows=1, ncols=1,figsize=(8, 6))
    sns.stripplot(x="isotype", y='F2',hue='cancer',edgecolors="black",linewidth=1,alpha=0.7,data=df,s=15)
    sns.boxplot(x="isotype", y='F2',color='white',data=df).set_title(pr+" overlap, "+test+", top "+str(top_max),fontsize=15)
    
    ord=list(data['isotype'].unique())
    box_pairs=itertools.combinations(ord,2)
    test_results = add_stat_annotation(ax=axes,data=df, x="isotype", y='F2', order=ord, box_pairs=box_pairs,test=test,text_format='simple',loc='inside',fontsize='x-large', verbose=2)
    
    axes.set_xlabel('',fontsize=15)
    axes.set_ylabel('F2 metric',fontsize=15)
    axes.tick_params(labelsize=15)
    axes.legend(fontsize='x-large')
    
    #plt.savefig("Pictures/Overlap_by_isotype_"+pr.split('/')[0]+"_"+pr.split('/')[1]+"_top"+str(top_max)+"_F2_distinct_segments.png",bbox_inches='tight')
    plt.savefig("Pictures/Overlap_by_isotype_"+pr.split('/')[0]+"_"+pr.split('/')[1]+"_top"+str(top_max)+"_F2_pooled_samples.png",bbox_inches='tight')
    
    #plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/Overlap_by_isotype_"+pr.split('/')[0]+"_"+pr.split('/')[1]+"_top"+str(top_max)+"_F2_distinct_segments.pdf",bbox_inches='tight')
    plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/Overlap_by_isotype_"+pr.split('/')[0]+"_"+pr.split('/')[1]+"_top"+str(top_max)+"_F2_pooled_samples.pdf",bbox_inches='tight')


# In[19]:


data[data['patient_1']==data['patient_2']]


# In[ ]:




