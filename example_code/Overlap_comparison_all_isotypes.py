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


# In[2]:


sort_type='bulk'
path_to_files='/home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Samples_aa/VDJ/'

izlist=['IGHA','IGHG','IGHM']

top_max=1000
top_min=101

files=[f for f in os.listdir(path_to_files) if re.match('.*'+sort_type+'.*'+'IGH'+'.*'+'.txt', f)]
files.sort()


# In[3]:


get_ipython().run_cell_magic('capture', '', '#NORMILIZING POOLED SAMPLES: SAME FOR ALL ISOTYPES\n\nif path.exists("PooledSamples_Top/")==0:\n    os.makedirs("PooledSamples_Top/")\n\nfor f in files:\n    c=rawincount(path_to_files+f)\n    if top_min<=c<top_max+1:\n        top_max=c-1\n\nprint(top_max)\n\n!rm PooledSamples_Top/*\n        \nfor f in files:\n    c=rawincount(path_to_files+f)\n    if c>=top_max:\n        !$VDJTOOLS SelectTop  -x $top_max $path_to_files$f PooledSamples_Top/')


# In[4]:


top_max


# In[5]:


get_ipython().run_cell_magic('capture', '', '#CALCULATING PAIRWISE DISTANCES (EACH ISOTYPE, CANCER, PATIENT, SEPARATELY)\nif path.exists("PwDistances/")==0:\n    os.makedirs("PwDistances/")\n\n!rm PwDistances/*\n    \nfor can in list(set([f.split(\'_\')[0] for f in os.listdir(\'PooledSamples_Top/\') if re.match(\'.*\'+sort_type+\'.*\'+\'IGH\'+\'.*\', f)])):\n    for pat in list(set([f.split(\'_\')[1] for f in os.listdir(\'PooledSamples_Top/\') if re.match(can+\'.*\'+sort_type+\'.*\'+\'IGH\'+\'.*\', f)])):\n            for iz in list(set([f.split(\'.\')[2] for f in os.listdir(\'PooledSamples_Top/\') if re.match(can+\'.*\'+pat+\'.*\'+sort_type+\'.*\'+\'IGH\'+\'.*\', f)])):\n                files=[f for f in os.listdir(\'PooledSamples_Top/\') if re.match(can+\'.*\'+pat+\'.*\'+sort_type+\'.*\'+iz+\'.txt\', f)]\n                if len(files)>=3:\n                    !$VDJTOOLS CalcPairwiseDistances PooledSamples_Top/$can*$pat*$sort_type*"$iz".txt PwDistances/"$can"_"$pat"_"$sort_type"_"$iz".txt')


# In[6]:


cols=['1_sample_id', '2_sample_id', 'div1', 'div2', 'div12', 'div21',
       'count1', 'count2', 'count12', 'count21', 'freq1', 'freq2', 'freq12',
       'freq21', 'R', 'Rs', 'D', 'F', 'F2']

data = pd.DataFrame(columns=cols)

files=[f for f in os.listdir('PwDistances/') if re.match('.*'+sort_type+'.*'+'IGH'+'.*', f)]

for f in files:
    df=pd.read_csv('PwDistances/'+f,comment="#", sep="\t", header=0)
    df=df[['1_sample_id', '2_sample_id', 'div1', 'div2', 'div12', 'div21','count1', 'count2', 'count12', 'count21', 'freq1', 'freq2', 'freq12','freq21', 'R', 'Rs', 'D', 'F', 'F2']].copy()
    data=data.append(df)


# In[7]:


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
data=data[data['count']==3]

#data=data[data['cancer']=='mel']


# In[8]:


ord=['LN/PBMC','LN/tum','PBMC/tum']
metrics=['D','R','F2']
test='Wilcoxon' #because samples are not normal and dependant

for m in metrics:
    figure, axes = plt.subplots(nrows=1, ncols=1,figsize=(8, 6))
    sns.scatterplot(x="pair", y=m,hue='isotype',edgecolors="black",linewidth=1,alpha=0.7,data=data,s=200).set_title(m+" metric, "+test+", Bonferroni corr, isotypes A, G, M, top "+str(top_max),fontsize=15)
    test_results = add_stat_annotation(ax=axes,data=data, x="pair", y=m, order=ord, box_pairs=[("LN/PBMC", "LN/tum"), ("LN/PBMC", "PBMC/tum"), ("PBMC/tum", "LN/tum")],test=test,text_format='simple',loc='inside',fontsize='x-large', verbose=2)
    plt.legend(loc='upper left')
    sns.pointplot(x="pair", y=m,color="grey",edgecolors="black",linewidth=1,alpha=0.3,data=data,markers="D",s=100,errwidth=1,ci=0,scale = 0.5,linestyles="--",ax=axes)
    axes.set_xlabel('',fontsize=15)
    axes.set_ylabel(m+' metric',fontsize=15)
    axes.tick_params(labelsize=15)
    axes.legend(fontsize='x-large')
    
    
    plt.savefig("Pictures/Overlaps_comparison_all_isotypes_top"+str(top_max)+'_'+m+".png",bbox_inches='tight')
    
    plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/Overlaps_comparison_all_isotypes_top"+str(top_max)+'_'+m+".pdf",bbox_inches='tight')


# In[10]:


len(data['patient'].unique())


# In[40]:


df1=data[data['pair']=='PBMC/tum']
df2=data[data['pair']=='LN/tum']
mt=pd.merge(df1,df2,how='outer',on=['patient','isotype'])

stats.wilcoxon(mt['F2_x'],mt['F2_y'])


# In[12]:


data['org2'].unique()


# In[23]:


data[data['pair']=='LN/PBMC'].sort_values('1_sample_id')


# In[24]:


data[data['pair']=='LN/tum'].sort_values('1_sample_id')


# In[ ]:




