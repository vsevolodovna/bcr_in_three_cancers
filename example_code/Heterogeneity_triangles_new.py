#!/usr/bin/env python
# coding: utf-8

# In[6]:


#NEWEST VERSION
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


# In[7]:


#POOLED REPLICAS
normalization_type='alldata' # 'alldata', 'top' or 'quantile' #ADD NORMALIZATION FOR ALL SAMPLES< SO WE WOULD HAVE MORE DATA

org1='tum' #on the angles by sample
org2='tum'

path_to_data="/home/sofyakrasik/Hypermutations/ChangeOut/"
folders=glob.glob(path_to_data+"*_pooled_replicas_"+normalization_type+"_CDR3_with_singles/")
files=glob.glob(path_to_data+"*_pooled_replicas_"+normalization_type+"_CDR3_with_singles_clone-pass.tsv")
files


# In[8]:


#READ DATA
data=pd.DataFrame()

for f in files:
    df=pd.read_csv(f,sep='\t')
    df=df[['sequence','clonefraction','clone_id', 'filename','cancer', 'patient', 'organ', 'sample']]
    data=data.append(df)


# In[9]:


#DEFINE PARAMETERS

l1=-0.1
l2=+0.1

size_threshold=5
org1_threshold=2
org2_threshold=1

def chi_goodness_of_fit(coordinates_list,colors,sizes):
    datapoints=list()
    for i in range(len(colors)):
        datapoints.append([xi*colors[i]*sizes[i]/5 for xi in coordinates_list[i]]) #Don't forget sizes[i]/5 divide by 5 because sizes are multiplied by 5 for proper visualization colors[i]*
    y=np.array([xi for xi in datapoints])
    y=y.transpose()
    f_obs=[sum(rw) for rw in y]
    f_exp=[sum(f_obs)/len(f_obs)]*len(f_obs)
    print(f_obs,f_exp)
    res=stats.chisquare(f_obs,f_exp)
    print(res.pvalue)
    mass_center=[x/sum(f_obs) for x in f_obs]
    return mass_center,res.pvalue

data=data[data['patient'].notna()]
patients=list(data['patient'].unique())
patients.sort()


# In[10]:


#patients=['ccp3']
organ_list=list(data['organ'].unique())

for pat in patients:
    print(pat)
    df=data[data['patient']==pat].copy()
    
    #Take only clonal lineages larger than threshold
    df['cluster_size']=df.groupby(['patient','clone_id'])['sequence'].transform('count')
    df=df[df['cluster_size']>=size_threshold]
    
    df=df[df['organ'].isin([org1,org2])].copy()
    if (org2 in list(df['organ'].unique())) and (df[df['organ']==org1]['sample'].nunique()>=3):
        
        #Taking the same top from each of org1 samples
        data_with_top=pd.DataFrame()
        for org in organ_list:
            df1=df[df['organ']==org].copy()
            if org==org1:
                df1=df1.sort_values('clonefraction',ascending=False)
                df1=df1.groupby(['sample']).head(df1.groupby(['sample'])['sequence'].transform('count').min())
            data_with_top=data_with_top.append(df1)
        df=data_with_top.copy()
        
        #If number of org1 samples is bigger than 3, take all possible combinations:
        org1_sample_list=list(df[df['organ']==org1]['sample'].unique())
        org1_samples_combinations=itertools.combinations(org1_sample_list, 3)
        df_save=df.copy() #save dataframe to iterate further
        idx=0
        for comb in org1_samples_combinations:
            idx=idx+1
            df=df_save[(df_save['organ']!=org1) | (df_save['sample'].isin(list(comb)))].copy()
            #df['cluster_size']=df.groupby(['patient','clone_id'])['sequence'].transform('count')
            #df=df[df['cluster_size']>=size_threshold]
            
            #Define lineages that contain both org1 and org2
            grouped=df.groupby(['clone_id'])
            df=grouped.filter(lambda x: (x['organ'].str.contains(org1).any())&(x['organ'].str.contains(org2).any()))
            
        
            #Define org2 share in each group by number of unique sequences
            if org1==org2:
                df['org2_share']=df.groupby(['clone_id'])['organ'].transform(lambda x: (x==org2).sum()/len(x))
                tshd_str=", size tshd="+str(size_threshold)
                triangle_wd=6.5
                title_str=pat+", "+org2+" vs "+org1+" samples,\nChi2 pvalue=" #add tissue threshold?
            else:
                df['org2_share']=df.groupby(['clone_id'])['organ'].transform(lambda x: (x==org2).sum()/((x==org2).sum()+(x==org1).sum()))
                tshd_str=", size tshd="+str(size_threshold)+", "+org1+" tshd="+str(org1_threshold)+", "+org2+" tshd="+str(org2_threshold)
                triangle_wd=8
                title_str=pat+", "+org2+" vs "+org1+" samples,\nmass center of coordinates weighted by "+org2+" share,\nChi2 pvalue="
            
            #Optional: condition on minimal number of clonotypes of each organ in the lineage
            df=df.groupby(['clone_id']).filter(lambda x: len(x[x['organ']==org1])>=org1_threshold)
            df=df.groupby(['clone_id']).filter(lambda x: len(x[x['organ']==org2])>=org2_threshold)
        
            #org1 samples
            org1_samples=list(df[df['organ']==org1]['sample'].unique())
            org1_samples.sort()
        
            #Create pivot table
            df=df.groupby(['clone_id','cluster_size','organ','sample','org2_share']).agg({'sequence':'count'})
            df=df.reset_index()
            df=df.pivot(index=['clone_id','cluster_size','org2_share'],columns='sample',values=['sequence'])
            df.columns = df.columns.droplevel(0)
            df = df.reset_index()
            df=df[['clone_id','cluster_size','org2_share']+org1_samples]
            df[org1_samples]=df[org1_samples].fillna(0)
        
            #Calculating coordinates on triangle
            df['coordinates_fixed']=[[x/(x+y+z),y/(x+y+z),z/(x+y+z)] for x,y,z in zip(df[org1_samples[0]],df[org1_samples[1]],df[org1_samples[2]])]
            df['coordinates']=[[(x+np.random.uniform(l1, l2))/(x+y+z),(y+np.random.uniform(l1, l2))/(x+y+z),(z+np.random.uniform(l1, l2))/(x+y+z)] for x,y,z in zip(df[org1_samples[0]],df[org1_samples[1]],df[org1_samples[2]])]
            points=list(df['coordinates'])
            points_fixed=list(df['coordinates_fixed'])
            sizes=list(df['cluster_size']*5)
            colors=list(df['org2_share'])
        
            scale = 1
            fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(triangle_wd, 5))
            ax.axis("off")
            figure, tax = ternary.figure(ax=ax,scale=scale)

            tax.gridlines(multiple=0.25, color="blue")
            if org1=='LN' and org2=='tum':
                tax.scatter(points,s=sizes,vmax=max(colors), colormap=plt.cm.coolwarm, colorbar=True, c=colors, cmap=plt.cm.coolwarm,cbarlabel=org2+' prevalence ->',alpha=0.7) #for tum vs LN samples
            elif org1=='tum' and org2=='LN':
                tax.scatter(points,s=sizes,vmax=max(colors), colormap=plt.cm.YlGn, colorbar=True, c=colors, cmap=plt.cm.YlGn,cbarlabel=org2+' prevalence ->',alpha=0.7) #for LN vs tum samples
                #figure.colorbar(points,location='right', shrink=0.6)
            elif org1=='tum' and org2=='tum':
                tax.scatter(points,s=sizes,vmax=max(colors), colormap=plt.cm.Reds, c=colors, cmap=plt.cm.Reds,cbarlabel=org2+' prevalence ->',alpha=0.3) #for tum vs tum add colorbar=True for legend
            else: #org1=='LN' and org2=='tum'
                tax.scatter(points,s=sizes,vmax=max(colors), colormap=plt.cm.Greens, c=colors, cmap=plt.cm.Greens,cbarlabel=org2+' prevalence ->',alpha=0.3) #for LN vs LN samples
            plt.xlim(-0.15,1.15)
            plt.ylim(-0.1,0.95)
            tax.set_title(pat+", "+org2+" vs "+org1+" samples", fontsize=15, pad=30)
            tax.right_corner_label(org1_samples[0], fontsize=13)
            tax.top_corner_label(org1_samples[1], fontsize=13)
            tax.left_corner_label(org1_samples[2], fontsize=13)
            #plt.axis('off')
            
            mass_center,pval=chi_goodness_of_fit(points_fixed,colors,sizes)
            tax.scatter([mass_center],marker="*",s=100,color='black').set_title(title_str+str(round(pval,6))+tshd_str, fontsize=15, pad=30)
            
            #tax.savefig("Pictures/Heterogeneity_triangle_"+pat+"_"+org2+"_vs_"+org1+"_samples_all_data_normalized_"+str(idx)+".png",bbox_inches='tight')
            #tax.savefig("/home/sofyakrasik/Pictures/BCR_Manuscript_Pictures_20220331/Heterogeneity_triangle_"+pat+"_"+org2+"_vs_"+org1+"_samples_all_data_normalized_"+str(idx)+".pdf",bbox_inches='tight')


# In[ ]:




