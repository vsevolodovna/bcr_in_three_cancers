#!/usr/bin/env python
# coding: utf-8

# In[2]:


import csv
import numpy as np
import pandas as pd
import os
from os import path
import re
from shutil import copy
from itertools import (takewhile,repeat)
import glob

def rawincount(filename):
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )

VDJTOOLS='/software/bin/vdjtools'


# In[3]:


path_to_data='/home/sofyakrasik/BCR_CDR3_data/5RACE/All_data_filtered/' #data in mixcr format
top_max=100

izlist=['IGH','IGHA','IGHM','IGHG','IGHD','IGHE']


path_to_mixcr_by_isotype="/home/sofyakrasik/BCR_CDR3_data/5RACE/All_data_filtered_byIsotype/"
path_to_VDJconverted="/home/sofyakrasik/BCR_CDR3_data/5RACE/VDJconverted/"
path_to_pooled_replicas_aa="/home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Replicas_aa/"
path_to_pooled_samples_aa="/home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Samples_aa/"
path_to_top="/home/sofyakrasik/BCR_CDR3_data/5RACE/Top"+str(top_max)+"_or_max/"
path_to_patient_top="/home/sofyakrasik/BCR_CDR3_data/5RACE/Top_for_patient/"
path_too_patient_quantile="/home/sofyakrasik/BCR_CDR3_data/5RACE/Top_Quantile/"
path_to_isotype_composition="/home/sofyakrasik/BCR_CDR3_data/5RACE/IsotypeComposition/"
path_to_vdj_with_aa_properties="/home/sofyakrasik/BCR_CDR3_data/5RACE/VDJ_with_aa_properties/"
path_to_aa_properties_stats="/home/sofyakrasik/BCR_CDR3_data/5RACE/AminoAcidPropertiesStats/"
path_to_pooled_replicas_nt="/home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Replicas_nt/"
path_to_pooled_samples_nt="/home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Samples_nt/"
path_to_patient_top_alternative="/home/sofyakrasik/BCR_CDR3_data/5RACE/Top_for_patient_Alternative/"

if path.exists(path_to_data)==0:
    os.makedirs(path_to_data)

if path.exists(path_to_mixcr_by_isotype)==0:
    os.makedirs(path_to_mixcr_by_isotype)

if path.exists(path_to_VDJconverted)==0:
    os.makedirs(path_to_VDJconverted)
    
if path.exists(path_to_pooled_replicas_aa)==0:
    os.makedirs(path_to_pooled_replicas_aa)
    os.makedirs(path_to_pooled_replicas_aa+"mixcr/")
    os.makedirs(path_to_pooled_replicas_aa+"VDJ/")

if path.exists(path_to_pooled_replicas_nt)==0:
    os.makedirs(path_to_pooled_replicas_nt)
    os.makedirs(path_to_pooled_replicas_nt+"mixcr/")
    os.makedirs(path_to_pooled_replicas_nt+"VDJ/")
    
if path.exists(path_to_pooled_samples_aa)==0:
    os.makedirs(path_to_pooled_samples_aa)
    os.makedirs(path_to_pooled_samples_aa+"mixcr/")
    os.makedirs(path_to_pooled_samples_aa+"VDJ/")
    
if path.exists(path_to_pooled_samples_nt)==0:
    os.makedirs(path_to_pooled_samples_nt)
    os.makedirs(path_to_pooled_samples_nt+"mixcr/")
    os.makedirs(path_to_pooled_samples_nt+"VDJ/")
    
if path.exists(path_to_top)==0:
    os.makedirs(path_to_top)
    os.makedirs(path_to_top+"Pooled_Replicas/")
    os.makedirs(path_to_top+"Pooled_Replicas/mixcr/")
    os.makedirs(path_to_top+"Pooled_Replicas/VDJ/")
    os.makedirs(path_to_top+"Pooled_Samples/")
    os.makedirs(path_to_top+"Pooled_Samples/mixcr/")
    os.makedirs(path_to_top+"Pooled_Samples/VDJ/")
    os.makedirs(path_to_top+"Separate_Replicas/")
    os.makedirs(path_to_top+"Separate_Replicas/mixcr/")
    os.makedirs(path_to_top+"Separate_Replicas/VDJ/")

if path.exists(path_to_patient_top)==0:
    os.makedirs(path_to_patient_top)
    os.makedirs(path_to_patient_top+"Pooled_Replicas/")
    os.makedirs(path_to_patient_top+"Pooled_Replicas/mixcr/")
    os.makedirs(path_to_patient_top+"Pooled_Replicas/VDJ/")
    os.makedirs(path_to_patient_top+"Pooled_Samples/")
    os.makedirs(path_to_patient_top+"Pooled_Samples/mixcr/")
    os.makedirs(path_to_patient_top+"Pooled_Samples/VDJ/")
    os.makedirs(path_to_patient_top+"Separate_Replicas/")
    os.makedirs(path_to_patient_top+"Separate_Replicas/mixcr/")
    os.makedirs(path_to_patient_top+"Separate_Replicas/VDJ/")
    
    
if path.exists(path_too_patient_quantile)==0:
    os.makedirs(path_too_patient_quantile)
    os.makedirs(path_too_patient_quantile+"Pooled_Replicas/")
    os.makedirs(path_too_patient_quantile+"Pooled_Replicas/mixcr/")
    os.makedirs(path_too_patient_quantile+"Pooled_Replicas/VDJ/")
    os.makedirs(path_too_patient_quantile+"Pooled_Samples/")
    os.makedirs(path_too_patient_quantile+"Pooled_Samples/mixcr/")
    os.makedirs(path_too_patient_quantile+"Pooled_Samples/VDJ/")
    os.makedirs(path_too_patient_quantile+"Separate_Replicas/")
    os.makedirs(path_too_patient_quantile+"Separate_Replicas/mixcr/")
    os.makedirs(path_too_patient_quantile+"Separate_Replicas/VDJ/")
    
if path.exists(path_to_isotype_composition)==0:
    os.makedirs(path_to_isotype_composition)
    
if path.exists(path_to_vdj_with_aa_properties)==0:
    os.makedirs(path_to_vdj_with_aa_properties)
    os.makedirs(path_to_vdj_with_aa_properties+"Pooled_Replicas/")
    os.makedirs(path_to_vdj_with_aa_properties+"Pooled_Samples/")
    os.makedirs(path_to_vdj_with_aa_properties+"Separate_Replicas/")  
    
if path.exists(path_to_aa_properties_stats)==0:
    os.makedirs(path_to_aa_properties_stats)
    os.makedirs(path_to_aa_properties_stats+"Pooled_Replicas/")
    os.makedirs(path_to_aa_properties_stats+"Pooled_Samples/")
    os.makedirs(path_to_aa_properties_stats+"Separate_Replicas/") 

if path.exists(path_to_patient_top_alternative)==0:
    os.makedirs(path_to_patient_top_alternative)
    os.makedirs(path_to_patient_top_alternative+"Pooled_Replicas/")
    os.makedirs(path_to_patient_top_alternative+"Pooled_Replicas/mixcr/")
    os.makedirs(path_to_patient_top_alternative+"Pooled_Replicas/VDJ/")
    os.makedirs(path_to_patient_top_alternative+"Pooled_Samples/")
    os.makedirs(path_to_patient_top_alternative+"Pooled_Samples/mixcr/")
    os.makedirs(path_to_patient_top_alternative+"Pooled_Samples/VDJ/")
    os.makedirs(path_to_patient_top_alternative+"Separate_Replicas/")
    os.makedirs(path_to_patient_top_alternative+"Separate_Replicas/mixcr/")
    os.makedirs(path_to_patient_top_alternative+"Separate_Replicas/VDJ/")


# In[55]:


#Dividing MIXCR by Isotype
files=[f for f in os.listdir(path_to_data) if re.match('.*'+'clonotypes.IGH.*.txt', f)]

for f in files:
    #copy(path_to_data+"/"+f, "/home/sofyakrasik/All_data_filtered_byIsotype/")
    basename=f.split('.')[0]
    df=pd.read_csv(path_to_data+"/"+f,comment="#", sep="\t", header=0)
    df=df[~df['allCHitsWithScore'].isnull()]
    for iz in izlist:
        df_is=df[df['allCHitsWithScore'].str.contains(iz)].copy()
        df_is['cloneFraction']=df_is['cloneFraction']/df_is['cloneFraction'].sum()
        df_is.to_csv(path_to_mixcr_by_isotype+basename+'.clonotypes.'+iz+'.txt',index=False,sep='\t')


# In[56]:


get_ipython().run_cell_magic('capture', '', "#Converting MIXCR to VDJ\nfiles=[f for f in os.listdir(path_to_mixcr_by_isotype) if re.match('.*'+'clonotypes.IGH.*.txt', f)]\n\nfor f in files:\n    !$VDJTOOLS Convert -S mixcr $path_to_mixcr_by_isotype*$f $path_to_VDJconverted")


# In[57]:


#Pooling replicas
def pool_replicas(path_to_files,path_to_save,columns_to_groupby,freq_column,count_column,izlist):
    files=[f for f in os.listdir(path_to_files) if re.match('.*'+'.clonotypes.IGH.*.txt', f)] 
    basenames=list(set([x.split('.')[0][:-2] for x in files]))
    basenames.sort() #why no
    for bname in basenames:
        for iz in izlist:
            replist=[f for f in os.listdir(path_to_files) if re.match(bname+'.*.clonotypes.'+iz+'.txt', f)]
            repl_cnt=len(replist)
            df=pd.DataFrame()
            for f in replist:
                data=pd.read_csv(path_to_files+f,comment="#", sep="\t", header=0)
                df=df.append(data)
            df[freq_column]=df.groupby(columns_to_groupby)[freq_column].transform('sum')/repl_cnt
            df[count_column]=df[count_column].sum()*df[freq_column]
            df[count_column]=df[count_column].round()
            df[count_column]=df[count_column].apply(np.ceil).astype(int, errors='ignore')
            df=df.drop_duplicates(columns_to_groupby).copy()
            
            df.to_csv(path_to_save+bname+'_pooled.clonotypes.'+iz+'.txt',sep='\t', index=False)
    #return df

#POOL BY AA (ONLY MAKES SENSE FOR CDR3 SEQUENCES)
pool_replicas(path_to_data,path_to_pooled_replicas_aa+"mixcr/",['allCHitsWithScore','aaSeqCDR3'],'cloneFraction','cloneCount',['IGH']) 
pool_replicas(path_to_VDJconverted,path_to_pooled_replicas_aa+"VDJ/",['cdr3aa'],'freq','count',['IGH','IGHA','IGHG','IGHM']) 

#POOL BY NT
pool_replicas(path_to_data,path_to_pooled_replicas_nt+"mixcr/",['allCHitsWithScore','targetSequences'],'cloneFraction','cloneCount',['IGH']) 
pool_replicas(path_to_VDJconverted,path_to_pooled_replicas_nt+"VDJ/",['cdr3nt'],'freq','count',['IGH','IGHA','IGHG','IGHM']) 


# In[58]:


#Pool samples
def pool_samples(path_to_files,path_to_save,columns_to_groupby,freq_column,count_column,izlist):
    files=[f for f in os.listdir(path_to_files) if re.match('.*'+'_pooled.clonotypes.IGH.*.txt', f)] 
    basenames=list(set([x.split('.')[0][:-11] for x in files]))
    basenames.sort() #why not
    for bname in basenames:
        for iz in izlist:
            samplelist=[f for f in os.listdir(path_to_files) if re.match(bname+'.*_pooled.clonotypes.'+iz+'.txt', f)]
            sample_cnt=len(samplelist)
            df=pd.DataFrame()
            for f in samplelist:
                data=pd.read_csv(path_to_files+f,comment="#", sep="\t", header=0)
                df=df.append(data)
            df[freq_column]=df.groupby(columns_to_groupby)[freq_column].transform('sum')/sample_cnt
            df[count_column]=df[count_column].sum()*df[freq_column]
            df[count_column]=df[count_column].round()
            df[count_column]=df[count_column].apply(np.ceil).astype(int, errors='ignore')
            df=df.drop_duplicates(columns_to_groupby).copy()
            
            df.to_csv(path_to_save+bname+'_pooled.clonotypes.'+iz+'.txt',sep='\t', index=False)
    #return basenames
#POOL BY AA
pool_samples(path_to_pooled_replicas_aa+"mixcr/",path_to_pooled_samples_aa+"mixcr/",['allCHitsWithScore','aaSeqCDR3'],'cloneFraction','cloneCount',['IGH'])
pool_samples(path_to_pooled_replicas_aa+"VDJ/",path_to_pooled_samples_aa+"VDJ/",['cdr3aa'],'freq','count',['IGH','IGHA','IGHG','IGHM'])

#POOL BY NT
pool_samples(path_to_pooled_replicas_nt+"mixcr/",path_to_pooled_samples_nt+"mixcr/",['allCHitsWithScore','targetSequences'],'cloneFraction','cloneCount',['IGH'])
pool_samples(path_to_pooled_replicas_nt+"VDJ/",path_to_pooled_samples_nt+"VDJ/",['cdr3nt'],'freq','count',['IGH','IGHA','IGHG','IGHM'])


# In[78]:


#Take Top
def take_top(path_to_files,path_to_save,freq_column,izlist,sort_type):
    print(path_to_files,path_to_save)
    for f in glob.glob(path_to_save+'*'+sort_type+'*'):
        os.remove(f)
    files=[f for f in os.listdir(path_to_files) if re.match('.*'+sort_type+'.*'+'.clonotypes.IGH.*.txt', f)] 
    patients=list(set([x.split('_')[1] for x in files]))
    patients.sort()    
    for iz in izlist:
        for p in patients:
            sample_list=[f for f in os.listdir(path_to_files) if re.match('.*'+p+'.*'+sort_type+'.*.clonotypes.'+iz+'.txt', f)]
            topn=top_max
            for f in sample_list:
                cnt=rawincount(path_to_files+f)
                if cnt<topn+1:
                    if cnt<50: #If file is too small, we do not use it
                        sample_list.remove(f)
                    else:
                        topn=cnt-1
            for f in sample_list:
                #bname=f.split('.')[0].split('_')[0]+'_'+f.split('.')[0].split('_')[1]+'_'+f.split('.')[0].split('_')[2]+'_'+f.split('.')[0].split('_')[3]
                bname=f.split('_pooled')[0].split('.clonotypes')[0]
                df=pd.read_csv(path_to_files+f,comment="#", sep="\t", header=0)
                df = df.sample(frac=1) #shuffle tail of counts
                df=df.sort_values(freq_column,ascending=False)
                df=df.head(topn).copy()
                df[freq_column]=df[freq_column]/df[freq_column].sum()
                df.to_csv(path_to_save+bname+'_top'+str(topn)+'.clonotypes.'+iz+'.txt',sep='\t', index=False)
    #return sample_list


# In[79]:


#If need to distinguish between bulk and PC_sort
take_top(path_to_pooled_replicas_nt+"mixcr/",path_to_top+"Pooled_Replicas/mixcr/",'cloneFraction',['IGH'],'PCsort')
take_top(path_to_pooled_replicas_nt+"mixcr/",path_to_top+"Pooled_Replicas/mixcr/",'cloneFraction',['IGH'],'bulk')
take_top(path_to_pooled_replicas_nt+"VDJ/",path_to_top+"Pooled_Replicas/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'PCsort')
take_top(path_to_pooled_replicas_nt+"VDJ/",path_to_top+"Pooled_Replicas/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'bulk')
take_top(path_to_pooled_samples_nt+"mixcr/",path_to_top+"Pooled_Samples/mixcr/",'cloneFraction',['IGH'],'PCsort')
take_top(path_to_pooled_samples_nt+"mixcr/",path_to_top+"Pooled_Samples/mixcr/",'cloneFraction',['IGH'],'bulk')
take_top(path_to_pooled_samples_nt+"VDJ/",path_to_top+"Pooled_Samples/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'PCsort')
take_top(path_to_pooled_samples_nt+"VDJ/",path_to_top+"Pooled_Samples/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'bulk')
take_top(path_to_data,path_to_top+"Separate_Replicas/mixcr/",'cloneFraction',['IGH'],'PCsort')
take_top(path_to_data,path_to_top+"Separate_Replicas/mixcr/",'cloneFraction',['IGH'],'bulk')
take_top(path_to_VDJconverted,path_to_top+"Separate_Replicas/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'PCsort')
take_top(path_to_VDJconverted,path_to_top+"Separate_Replicas/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'bulk')


# In[ ]:


#If need same top for bulk and PC_sort
take_top(path_to_pooled_replicas_nt+"mixcr/",path_to_top+"/Pooled_Replicas/mixcr/",'cloneFraction',['IGH'],'')
take_top(path_to_pooled_replicas_nt+"VDJ/",path_to_top+"/Pooled_Replicas/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'')
take_top(path_to_pooled_samples_nt+"mixcr/",path_to_top+"/Pooled_Samples/mixcr/",'cloneFraction',['IGH'],'')
take_top(path_to_pooled_samples_nt+"VDJ/",path_to_top+"/Pooled_Samples/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'')
take_top(path_to_data,path_to_top+"Separate_Replicas/mixcr/",'cloneFraction',['IGH'],'')
take_top(path_to_VDJconverted,path_to_top+"Separate_Replicas/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'')


# In[ ]:


#If need for expanded
#take_top("/home/sofyakrasik/Pooled_Samples/VDJ/","/home/sofyakrasik/ExpandedClonotypes/Data/",'freq',['IGH','IGHA','IGHG','IGHM'],'bulk')
take_top(path_to_pooled_replicas+"VDJ/","/home/sofyakrasik/ExpandedClonotypes/Data/",'freq',['IGH','IGHA','IGHG','IGHM'],'bulk')


# In[62]:


#Isotype composition, MIXCR format
def isotype_composition(path_to_files,path_to_save):
    files=[f for f in os.listdir(path_to_files) if re.match('.*IGH.txt',f)]
    for f in files:
        df=pd.read_csv(path_to_files+f,comment="#", sep="\t", header=0)
        if len(df)>=100: #Only trust isotype composition if sample is big enough
            df['allCHitsWithScore']=[str(x)[0:4] for x in df['allCHitsWithScore']]
            df=df[df['allCHitsWithScore'].str.contains('IGH')].copy()
            df['filename']=f
            df=df.groupby(['filename','allCHitsWithScore']).agg({'cloneFraction':'sum'})
            df=df.reset_index()
            df.to_csv(path_to_save+f,sep='\t', index=False)
    #return df

isotype_composition(path_to_data,path_to_isotype_composition)


# In[63]:


#Take Patient Top, unlimited
def take_patient_top(path_to_files,path_to_save,freq_column,izlist,sort_type):
    for f in glob.glob(path_to_save+'*'+sort_type+'*'):
        os.remove(f)
    files=[f for f in os.listdir(path_to_files) if re.match('.*'+sort_type+'.*'+'.clonotypes.IGH.*.txt', f)] 
    patients=list(set([x.split('_')[1] for x in files]))
    patients.sort()    
    for iz in izlist:
        for p in patients:
            sample_list=[f for f in os.listdir(path_to_files) if re.match('.*'+p+'.*'+sort_type+'.*.clonotypes.'+iz+'.txt', f)]
            topn=rawincount(path_to_files+sample_list[0])
            for f in sample_list:
                cnt=rawincount(path_to_files+f)
                if cnt<topn+1:
                    if cnt<50:
                        sample_list.remove(f)
                    else:
                        topn=cnt-1
            topn=round(topn*0.9)
            for f in sample_list:
                #bname=f.split('.')[0].split('_')[0]+'_'+f.split('.')[0].split('_')[1]+'_'+f.split('.')[0].split('_')[2]+'_'+f.split('.')[0].split('_')[3]
                bname=f.split('_pooled')[0].split('.clonotypes')[0]
                df=pd.read_csv(path_to_files+f,comment="#", sep="\t", header=0)
                df = df.sample(frac=1) #shuffle tail of counts
                df=df.sort_values(freq_column,ascending=False)
                df=df.head(topn).copy()
                df[freq_column]=df[freq_column]/df[freq_column].sum()
                df.to_csv(path_to_save+bname+'_top'+str(topn)+'.clonotypes.'+iz+'.txt',sep='\t', index=False)
    #return sample_list


# In[41]:


#If need to distinguish between bulk and PC_sort
take_patient_top(path_to_pooled_replicas_nt+"mixcr/",path_to_patient_top+"Pooled_Replicas/mixcr/",'cloneFraction',['IGH'],'PCsort')
take_patient_top(path_to_pooled_replicas_nt+"mixcr/",path_to_patient_top+"Pooled_Replicas/mixcr/",'cloneFraction',['IGH'],'bulk')
take_patient_top(path_to_pooled_replicas_nt+"VDJ/",path_to_patient_top+"Pooled_Replicas/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'PCsort')
take_patient_top(path_to_pooled_replicas_nt+"VDJ/",path_to_patient_top+"Pooled_Replicas/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'bulk')
take_patient_top(path_to_pooled_samples_nt+"mixcr/",path_to_patient_top+"Pooled_Samples/mixcr/",'cloneFraction',['IGH'],'PCsort')
take_patient_top(path_to_pooled_samples_nt+"mixcr/",path_to_patient_top+"Pooled_Samples/mixcr/",'cloneFraction',['IGH'],'bulk')
take_patient_top(path_to_pooled_samples_nt+"VDJ/",path_to_patient_top+"Pooled_Samples/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'PCsort')
take_patient_top(path_to_pooled_samples_nt+"VDJ/",path_to_patient_top+"Pooled_Samples/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'bulk')
take_patient_top(path_to_data,path_to_patient_top+"Separate_Replicas/mixcr/",'cloneFraction',['IGH'],'PCsort')
take_patient_top(path_to_data,path_to_patient_top+"Separate_Replicas/mixcr/",'cloneFraction',['IGH'],'bulk')
take_patient_top(path_to_VDJconverted,path_to_patient_top+"Separate_Replicas/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'PCsort')
take_patient_top(path_to_VDJconverted,path_to_patient_top+"Separate_Replicas/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'bulk')


# In[64]:


#If need same top for bulk and PC_sort
take_patient_top(path_to_pooled_replicas_nt+"mixcr/",path_to_patient_top+"Pooled_Replicas/mixcr/",'cloneFraction',['IGH'],'')
take_patient_top(path_to_pooled_replicas_nt+"VDJ/",path_to_patient_top+"Pooled_Replicas/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'')
take_patient_top(path_to_pooled_samples_nt+"mixcr/",path_to_patient_top+"Pooled_Samples/mixcr/",'cloneFraction',['IGH'],'')
take_patient_top(path_to_pooled_samples_nt+"VDJ/",path_to_patient_top+"Pooled_Samples/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'')
take_patient_top(path_to_data,path_to_patient_top+"Separate_Replicas/mixcr/",'cloneFraction',['IGH'],'')
take_patient_top(path_to_VDJconverted,path_to_patient_top+"Separate_Replicas/VDJ/",'freq',['IGH','IGHA','IGHG','IGHM'],'')


# In[4]:


#Take Patient Top ALTERNATIVE, unlimited
def take_patient_top_alternative(path_to_files,path_to_save,freq_column,izlist,sort_type):
    for f in glob.glob(path_to_save+'*'+sort_type+'*'):
        os.remove(f)
    files=[f for f in os.listdir(path_to_files) if re.match('.*'+sort_type+'.*'+'.clonotypes.IGH.*.txt', f)] 
    patients=list(set([x.split('_')[1] for x in files]))
    patients.sort()    
    for iz in izlist:
        for p in patients:
            sample_list=[f for f in os.listdir(path_to_files) if re.match('.*'+p+'.*'+sort_type+'.*.clonotypes.'+iz+'.txt', f)]
            topn=rawincount(path_to_files+sample_list[0])
            for f in sample_list:
                cnt=rawincount(path_to_files+f)
                if cnt<topn+1:
                    if cnt<50:
                        sample_list.remove(f)
                    else:
                        topn=cnt-1
            topn=round(topn*0.9)
            for f in sample_list:
                #bname=f.split('.')[0].split('_')[0]+'_'+f.split('.')[0].split('_')[1]+'_'+f.split('.')[0].split('_')[2]+'_'+f.split('.')[0].split('_')[3]
                bname=f.split('_pooled')[0].split('.clonotypes')[0]
                df=pd.read_csv(path_to_files+f,comment="#", sep="\t", header=0)
                df = df.sample(frac=1) #shuffle tail of counts
                df=df.sort_values(freq_column,ascending=False)
                df=df.head(topn).copy()
                df[freq_column]=df[freq_column]/df[freq_column].sum()
                df.to_csv(path_to_save+bname+'_top'+str(topn)+'.clonotypes.'+iz+'.txt',sep='\t', index=False)
    #return sample_list


# In[5]:


#If need to distinguish between bulk and PC_sort
take_patient_top_alternative("/home/sofyakrasik/BCR_CDR3_data/5RACE/Alternative_Replicas_Pooling_aa/",path_to_patient_top_alternative+"Pooled_Replicas/mixcr/",'cloneFraction',['IGH'],'PCsort')
take_patient_top_alternative("/home/sofyakrasik/BCR_CDR3_data/5RACE/Alternative_Replicas_Pooling_aa/",path_to_patient_top_alternative+"Pooled_Replicas/mixcr/",'cloneFraction',['IGH'],'bulk')


# In[65]:


get_ipython().run_cell_magic('capture', '', '#CALCULATE PAIRWISE DISTANCES\ndef pwd_calculation(path_to_files,path_to_save,izlist,sort_type):\n    if path.exists(path_to_save)==0:\n        os.makedirs(path_to_save)\n    files=[f for f in os.listdir(path_to_files) if re.match(\'.*\'+sort_type+\'.*\'+\'.clonotypes.IGH.*.txt\', f)] \n    patients=list(set([x.split(\'_\')[1] for x in files]))\n    patients.sort()    \n    for iz in izlist:\n        for p in patients:\n            print(p,iz)\n            #sample_list=[f for f in os.listdir(path_to_files) if re.match(\'.*\'+p+\'.*\'+sort_type+\'.*.clonotypes.\'+iz+\'.txt\', f)]\n            !$VDJTOOLS CalcPairwiseDistances $path_to_files*$p*$sort_type*$iz".txt" $path_to_save$p"_"$sort_type"_"$iz\n    return files\n            \npwd_calculation(path_to_patient_top+\'Separate_Replicas/VDJ/\',\'/home/sofyakrasik/PwDistances/SeparateReplicas_with_singles/\',[\'IGH\',\'IGHA\',\'IGHG\',\'IGHM\'],\'\')')


# In[66]:


#TAKE TOP QUANTILE
def take_quantile(path_to_files,path_to_save,freq_column,Qn):
    for f in glob.glob(path_to_save+'*'):
        os.remove(f)
    files=[f for f in os.listdir(path_to_files) if re.match('.*'+'.clonotypes.IGH.*.txt', f)] 
    for f in files:
        bname=f.split('.clonotypes')[0]
        endname=f.split('.clonotypes')[1]
        df=pd.read_csv(path_to_files+f,sep='\t')
        df=df.sort_values(freq_column,ascending=False)
        cum_freq=0
        n=0
        for freq in df[freq_column]:
            if cum_freq>Qn:
                break
            n=n+1
            cum_freq=cum_freq+freq
        df=df.head(n)
        df[freq_column]=df[freq_column]/df[freq_column].sum()
        df.to_csv(path_to_save+bname+'_top_'+str(round(Qn*100))+'_quantile.clonotypes'+endname,sep='\t', index=False)


# In[68]:


take_quantile(path_to_pooled_replicas_nt+"mixcr/",path_too_patient_quantile+"Pooled_Replicas/mixcr/",'cloneFraction',0.4)
take_quantile(path_to_pooled_replicas_nt+"VDJ/",path_too_patient_quantile+"Pooled_Replicas/VDJ/",'freq',0.4)
take_quantile(path_to_pooled_samples_nt+"mixcr/",path_too_patient_quantile+"Pooled_Samples/mixcr/",'cloneFraction',0.4)
take_quantile(path_to_pooled_samples_nt+"VDJ/",path_too_patient_quantile+"Pooled_Samples/VDJ/",'freq',0.4)
take_quantile(path_to_data,path_too_patient_quantile+"Separate_Replicas/mixcr/",'cloneFraction',0.4)
take_quantile(path_to_VDJconverted,path_too_patient_quantile+"Separate_Replicas/VDJ/",'freq',0.4)


# In[19]:


get_ipython().run_cell_magic('capture', '', "##VDJ FILES WITH AMINO ACID PROPERTIES: CHOOSE 7 FROM THE CENTER\ndef calculate_aa_properties(path_to_files,path_to_save):\n    for f in glob.glob(path_to_save+'*'):\n        os.remove(f)\n    files=[f for f in os.listdir(path_to_files) if re.match('.*'+'clonotypes.IGH.*.txt', f)]\n\n    for f in files:\n        df=pd.read_csv(path_to_files+f,sep='\\t')\n        df['cdr3aa']=[x[int((len(x)-5)/2):int((len(x)-5)/2)+5] for x in df['cdr3aa']]\n        df.to_csv(path_to_save+f,sep='\\t',index=None)\n        !$VDJTOOLS Annotate -a  cdr3length,strength,ndnsize,insertsize,hydropathy,charge,polarity,strength,cdr3contact,disorder,volume,mjenergy,kf1,kf2,kf3,kf4,kf5,kf6,kf7,kf8,kf9,kf10 $path_to_save*$f $path_to_save\n        \n#cdr3length,strength,ndnsize,insertsize,hydropathy,charge,polarity,strength,cdr3contact,kf1,kf2,kf3,kf4,kf5,kf6,kf7,kf8,kf9,kf10")


# In[20]:


get_ipython().run_cell_magic('capture', '', "calculate_aa_properties(path_to_pooled_replicas_nt+'VDJ/',path_to_vdj_with_aa_properties+'Pooled_Replicas/')\ncalculate_aa_properties(path_to_VDJconverted,path_to_vdj_with_aa_properties+'Separate_Replicas/')\ncalculate_aa_properties(path_to_pooled_samples_nt+'VDJ/',path_to_vdj_with_aa_properties+'Pooled_Samples/')")


# In[11]:


get_ipython().run_cell_magic('capture', '', "#AMINO ACID PROPERTIES STATS CDR3 CENTER 5, FOR TOP N CLONOTYPES\ndef aa_properties_stats(path_to_files,path_to_save):\n    for f in glob.glob(path_to_save+'*'):\n        os.remove(f)\n        \n    !$VDJTOOLS CalcCdrAaStats -r cdr3-center-5 -a hydropathy,strength,charge,disorder,volume,mjenergy,count,kf1,kf2,kf3,kf4,kf5,kf6,kf7,kf8,kf9,kf10 $path_to_files*IGH*txt $path_to_save\n    !$VDJTOOLS CalcCdrAaStats -r cdr3-center-5 -w -a hydropathy,strength,charge,disorder,volume,mjenergy,count,kf1,kf2,kf3,kf4,kf5,kf6,kf7,kf8,kf9,kf10 $path_to_files*IGH*txt $path_to_save\n    !$VDJTOOLS CalcCdrAaStats -r cdr3-center-5 -n -a hydropathy,strength,charge,disorder,volume,mjenergy,count,kf1,kf2,kf3,kf4,kf5,kf6,kf7,kf8,kf9,kf10 $path_to_files*IGH*txt $path_to_save\n    !$VDJTOOLS CalcCdrAaStats -r cdr3-center-5 -w -n -a hydropathy,strength,charge,disorder,volume,mjenergy,count,kf1,kf2,kf3,kf4,kf5,kf6,kf7,kf8,kf9,kf10 $path_to_files*IGH*txt $path_to_save\n    !$VDJTOOLS CalcCdrAaStats -r CDR3-full -a count $path_to_files*IGH*txt $path_to_save'cdr3_length'\n\naa_properties_stats(path_to_top+'Pooled_Samples/VDJ/',path_to_aa_properties_stats+'Pooled_Samples/')\naa_properties_stats(path_to_top+'Pooled_Replicas/VDJ/',path_to_aa_properties_stats+'Pooled_Replicas/')\naa_properties_stats(path_to_top+'Separate_Replicas/VDJ/',path_to_aa_properties_stats+'Separate_Replicas/')")


# In[ ]:




