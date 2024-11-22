#!/usr/bin/env python
# coding: utf-8

# In[2]:


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
from matplotlib.backends.backend_pdf import PdfPages
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
import itertools
from itertools import combinations
from itertools import (takewhile,repeat)
import glob

def rawincount(filename):
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )

VDJTOOLS='/software/bin/vdjtools'


# In[3]:


top_max=355
org_list=['tum','LN']
izlist=['A','G','M']

def build_scatterplots(path_to_data,pattern,sort_type,norm_type,group_column,compare_column):
    files=glob.glob(path_to_data+pattern+'*'+sort_type+'*.clonotypes.IGH.txt')
    
    #read data
    data=pd.DataFrame()

    for f in files:
        df=pd.read_csv(f,sep='\t')
        #df=df[['sequence','clonefraction','clone_id', 'filename','cancer', 'patient', 'organ','isotype']]
        df['filename']=f.split('/')[-1]
        data=data.append(df)
        
    data['cancer']=[x.split('_')[0] for x in data['filename']]
    data['patient']=[x.split('_')[1] for x in data['filename']]
    data['organ']=[x.split('_')[2] for x in data['filename']]
    data['sample']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5] if len(x.split('_'))>=6 else '-' for x in data['filename']]
    data['replicate']=[x.split('_')[2]+x.split('_')[4]+x.split('_')[5]+x.split('_')[6].split('.')[0] if len(x.split('_'))==7 else '-' for x in data['filename']]
    data=data[data['allCHitsWithScore'].notna()]
    data['isotype']=[x[3:4] for x in data['allCHitsWithScore']]
    data=data[data['isotype'].isin(izlist)]
    data['v_call']=[x.split('*')[0] for x in data['allVHitsWithScore']]
    data['j_call']=[x.split('*')[0] for x in data['allJHitsWithScore']]
    
    data['cloneFraction']=data.groupby(['filename','aaSeqCDR3','isotype','v_call','j_call'])['cloneFraction'].transform('sum')
    data=data.drop_duplicates(['filename','aaSeqCDR3','isotype','v_call','j_call','cloneFraction'])
    
    data=data.sort_values('cloneFraction',ascending=False)
    #data=data.groupby('filename').head(top_max)
    
    data=data[data['organ'].isin(org_list)].copy()
    
    for pat in list(data['patient'].unique()):
        #just to define common fmin and fmax
        df_temp=data[data['patient']==pat].copy()
        df_temp=df_temp.sort_values('cloneFraction',ascending=False)
        df_temp=df_temp.groupby(['filename','replicate']).head(min(df_temp.groupby(compare_column)['aaSeqCDR3'].transform('count')))
        df_temp['cloneFraction']=df_temp.groupby(['filename','replicate'])['cloneFraction'].transform(lambda x: x/x.sum())
        fmin=df_temp['cloneFraction'].min()
        fmax=df_temp['cloneFraction'].max()
        for smp in list(data[group_column].unique()):
            print(smp)
            df=data[(data['patient']==pat)&(data[group_column]==smp)].copy()
            df['files_per_group']=df.groupby(['patient',group_column,compare_column])['aaSeqCDR3'].transform('nunique')
            df['cloneFraction']=df.groupby(['patient',group_column,compare_column,'aaSeqCDR3'])['cloneFraction'].transform('sum')
            df['cloneFraction']=df['cloneFraction']/df['files_per_group']
            df=df.drop_duplicates(['patient',group_column,compare_column,'aaSeqCDR3'])                                                                                                                
            
            if len(list(df[compare_column].unique()))>1:
                compare_list=list(df[compare_column].unique())
                combinations=itertools.combinations(compare_list,2)
                for comb in combinations:
                    comb=list(comb)
                    print(comb)
                    comb.sort()
                    r1=comb[0]
                    r2=comb[1]
                
                    df1=df[df[compare_column]==r1].copy()
                    df2=df[df[compare_column]==r2].copy()
                
                    top_head=min(len(df1),len(df2),top_max)
                    df1=df1.head(top_head)
                    df2=df2.head(top_head)
                
                    if norm_type=='normalized':
                        df1['cloneFraction']=df1['cloneFraction']/df1['cloneFraction'].sum()
                        df2['cloneFraction']=df2['cloneFraction']/df2['cloneFraction'].sum()
                
                    fmin=min(df1['cloneFraction'].min(),df2['cloneFraction'].min(),fmin)
                    fmax=max(df1['cloneFraction'].max(),df2['cloneFraction'].max(),fmax)
                    
                    plot_data=pd.merge(df1,df2,on=['aaSeqCDR3','isotype','v_call','j_call'],how='outer')
                    plot_data=plot_data[['aaSeqCDR3','isotype','v_call','j_call','cloneFraction_x','cloneFraction_y']]
                    plot_data= plot_data.apply(lambda x: x.fillna(fmin*0.7))
                    plot_data=plot_data.sort_values('isotype')
                    
                    plt.figure(figsize=(5,5))
            
                    plt.xlim(fmin/2,fmax)
                    plt.ylim(fmin/2,fmax)
                    plt.xscale('log')
                    plt.yscale('log')
                    
                    #sns.scatterplot(x="cloneFraction_x", y="cloneFraction_y",alpha=0.5,data=plot_data,s=100,hue='isotype',linewidth=0)
                    sns.scatterplot(x="cloneFraction_x", y="cloneFraction_y",alpha=0.5,data=plot_data,s=100,linewidth=0,color='g')
                    plt.xlabel(pat+'_'+r1.split('.')[0], fontsize=15)
                    plt.ylabel(pat+'_'+r2.split('.')[0], fontsize=15)
                    plt.title('top='+f'{top_head:,}'.replace(',',' ')+', '+norm_type, fontsize=20)
                    #plt.legend(fontsize='x-large')
                    plt.tick_params(labelsize=15)
                    
                    plt.savefig('/home/sofyakrasik/Pipelines/Pictures/Scatterplots/Correlation_scatterplot_'+pat+'_'+compare_column+'s_'+sort_type+'_'+norm_type+'_'+r1+'_'+r2+'_top_'+str(top_head)+'.png',bbox_inches='tight')
                    plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/Scatterplots/Correlation_scatterplot_"+pat+'_'+compare_column+'s_'+sort_type+'_'+norm_type+'_'+r1+'_'+r2+'_top_'+str(top_head)+".pdf",bbox_inches='tight')
                    plt.show()
    
    #return plot_data
    


# In[14]:


pwd


# In[93]:


path_to_data='/home/sofyakrasik/BCR_CDR3_data/5RACE/All_data_filtered/'
pattern='*lcp1*'
sort_type='bulk'
norm_type='normalized'
group_column='sample'
compare_column='replicate'
build_scatterplots(path_to_data,pattern,sort_type,norm_type,group_column,compare_column)


# In[95]:


path_to_data='/home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Replicas/mixcr/'
pattern='*ccp6*'
sort_type='bulk'
norm_type='normalized'
group_column='organ'
compare_column='sample'
build_scatterplots(path_to_data,pattern,sort_type,norm_type,group_column,compare_column)


# In[26]:


path_to_data='/home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Samples/mixcr/'
pattern='*'
sort_type='bulk'
norm_type='normalized'
group_column='patient'
compare_column='organ'
build_scatterplots(path_to_data,pattern,sort_type,norm_type,group_column,compare_column)


# In[ ]:




