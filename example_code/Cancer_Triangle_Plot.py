#!/usr/bin/env python
# coding: utf-8

# In[1]:


from statannot import add_stat_annotation
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import ks_2samp
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


# In[2]:


#POOLED SAMPLES
normalization_type='alldata' # 'alldata', 'top' or 'quantile' #ADD NORMALIZATION FOR ALL SAMPLES< SO WE WOULD HAVE MORE DATA

path_to_data="/home/sofyakrasik/Hypermutations/ChangeOut/"
files=glob.glob(path_to_data+"*_pooled_samples_"+normalization_type+"_CDR3_with_singles_clone-pass.tsv")
files


# In[3]:


#READ DATA
data=pd.DataFrame()

for f in files:
    df=pd.read_csv(f,sep='\t')
    df=df[['sequence','clonefraction','clone_id', 'filename','cancer', 'patient', 'organ','isotype']]
    data=data.append(df)


# In[4]:


#DEFINE PARAMETERS
color_dict={"A":"blue","E":"brown","M":"green","G":"red","D":"yellow","no":"grey"}
isotype_list=['A','G','M']

def color_by_isotype(x):
    return color_dict[x]

l1=-0.1
l2=+0.1
dominance_treshold=0.6

size_threshold=4

data=data[data['patient'].notna()]
patients=list(data['patient'].unique())
patients.sort()

def chi_goodness_of_fit(coordinates_list,colors,sizes):
    datapoints=list()
    for i in range(len(colors)):
        datapoints.append([xi for xi in coordinates_list[i]]) #Don't forget sizes[i]/5 divide by 5 because sizes are multiplied by 5 for proper visualization colors[i]*
    y=np.array([xi for xi in datapoints])
    y=y.transpose()
    f_obs=[sum(rw) for rw in y]
    f_exp=[sum(f_obs)/len(f_obs)]*len(f_obs)
    res=stats.chisquare(f_obs,f_exp)
    #print(res.pvalue)
    mass_center=[x/sum(f_obs) for x in f_obs]
    return mass_center,res.pvalue


# In[5]:


#CALCULATE AND PLOT
#patients=['ccp3']
data_with_top=pd.DataFrame() #normalized data that will be used further

for pat in patients:
    df=data[data['patient']==pat].copy()
    df=df[df['organ']!='norm']
    if df['organ'].nunique()>=3:
        #Consider all trios of tissues (e.g. if norm is present)
        organ_list=['tum','LN','PBMC'] #list(df['organ'].unique())
        organ_combinations=itertools.combinations(organ_list, 3)
        
        df_save=df.copy() #save dataframe to iterate further
        
        idx=0
        for com in organ_combinations:
            idx=idx+1
            df=df_save[df_save['organ'].isin(com)]
            
            #BEFORE WE START CALCULATIONS, NORMALIZE EACH ISOTYPE EQUALLY AMONG TISSUES
            df=df.sort_values(['clonefraction'],ascending=False)
            
            for iz in isotype_list:
                df1=df[df['isotype']==iz].copy()
                df1=df1.groupby(['organ']).head(df1.groupby(['organ'])['sequence'].transform('count').min())
                data_with_top=data_with_top.append(df1)

cancer_list=list(data_with_top['cancer'].unique())                
for can in cancer_list:
    df=data_with_top[data_with_top['cancer']==can].copy()
    
    df['cluster_size']=df.groupby(['clone_id'])['sequence'].transform('count')
    df=df[df['cluster_size']>=size_threshold]
            
    if df['organ'].nunique()==3: #check that there are still 3 tissues left
        df['isotype_share']=df.groupby(['clone_id','isotype'])['sequence'].transform('count')/df['cluster_size']
        #df=df.groupby(['clone_id']).filter(lambda x: (x['isotype_share']>=dominance_treshold).any())
                
        df['dominant_isotype']=df.assign(dominant_isotype=df['isotype'].where(df['isotype_share']>=dominance_treshold)).fillna('').groupby('clone_id')['dominant_isotype'].transform(lambda x: x.max() if x.any() else 'no') #СЛОЖНА.. но вроде правильно
                
        df=df.groupby(['clone_id','cluster_size','organ','dominant_isotype']).agg({'sequence':'count'})
        df=df.reset_index()
        df=df.sort_values('organ')
        df=df.pivot(index=['clone_id','cluster_size','dominant_isotype'],columns='organ',values=['sequence'])
        df.columns = df.columns.droplevel(0)
        df = df.reset_index()
        df=df.fillna(0)
        df['coordinates_fixed']=[[x/(x+y+z),y/(x+y+z),z/(x+y+z)] for x,y,z in zip(df[df.columns[-3]],df[df.columns[-2]],df[df.columns[-1]])]
        df['coordinates']=[[(x+np.random.uniform(l1, l2))/(x+y+z),(y+np.random.uniform(l1, l2))/(x+y+z),(z+np.random.uniform(l1, l2))/(x+y+z)] for x,y,z in zip(df[df.columns[-4]],df[df.columns[-3]],df[df.columns[-2]])]
        
        df['color_by_dominant_isotype']=[color_by_isotype(x) for x in df['dominant_isotype']]
        df=df.sort_values('dominant_isotype')
        
        scale = 1
        fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6, 5))
        ax.axis("off")
        figure, tax = ternary.figure(ax=ax,scale=scale)
        tax.gridlines(multiple=0.25, color="blue")
        tax.right_corner_label(df.columns[-6], fontsize=20)
        tax.top_corner_label(df.columns[-5], fontsize=20)
        tax.left_corner_label(df.columns[-4], fontsize=20)
                
        stat='significance: '
        avg_tum_share='tum share: '
                
        for iz in list(df['dominant_isotype'].unique()):
            df2=df[df['dominant_isotype']==iz].copy()
            points=list(df2['coordinates'])
            points_fixed=list(df2['coordinates_fixed'])
            colors=list(df2['color_by_dominant_isotype'])
            sizes=list(df2['cluster_size']*5)
            tax.scatter(points,s=sizes,c=color_dict[iz],alpha=0.5,label=iz)
            mass_center,pval=chi_goodness_of_fit(points_fixed,colors,sizes)
            if iz!='no':
                stat=stat+iz+' '+str(round(pval,6))+', '
                avg_tum_share=avg_tum_share+iz+' '+str(round([sum(sub_list) / len(sub_list) for sub_list in zip(*points_fixed)][2],2))+', '
            tax.scatter([mass_center],marker="*",s=300,color=color_by_isotype(iz),edgecolor='black')
                    
            tax.legend(['blue','brown','red','green','yellow',"grey"], ["IGHA","IGHE","IGHM","IGHG","IGHD","no"])
                    
            tax.legend(fontsize='x-large')
            plt.xlim(-0.15,1.15)
            plt.ylim(-0.1,0.95)
            tax.set_title(can+", trees by tissue, pvalues:\n"+stat+'\n'+avg_tum_share+'\nmass center weighted by coordinates', fontsize=15, pad=30)
                    
            plt.savefig("Pictures/Cancer_traingle_plot_"+can+"_"+str(df.columns[-6])+'_'+str(df.columns[-5])+'_'+str(df.columns[-4])+"_by_isotype_"+normalization_type+"_normalized.png",bbox_inches='tight')
            plt.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/Cancer_traingle_plot_"+can+"_"+str(df.columns[-6])+'_'+str(df.columns[-5])+'_'+str(df.columns[-4])+"_by_isotype_"+normalization_type+"_normalized.pdf",bbox_inches='tight')
        
            


# In[33]:


column_average = [sum(sub_list) / len(sub_list) for sub_list in zip(*points_fixed)]
column_average


# In[ ]:




