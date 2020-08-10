# -*- coding: utf-8 -*-
"""
Code developed in the laboratory of Cancer Cell Biology,
Engelhardt Institute of Molecular Biology, Moscow, Russia.
Not for commercial use.
If you have any questions about the code or its use contact:
lebedevtd@gmai.com
vr.elmira@gmail.com 
"""

import pandas as pd
from scipy.stats import mannwhitneyu
from itertools import chain
import numpy as np
import os


def alt_names(input_drugs): #extract alternative drug names and pertubation ids from L1000 database
    df_meta= pd.read_csv('Datasets\\L1000FWD_drugs_metadata.csv')
    pert_lower=df_meta['pert_id'].str.lower()
    df_meta['pert_id']=pert_lower
    alt_lower=df_meta['alt_name'].str.lower()
    df_meta['alt_name']=alt_lower
    iname_lower=df_meta['pert_iname'].str.lower()
    df_meta['pert_iname']=iname_lower
    meta_names=df_meta[{'pert_id', 'alt_name', 'pert_iname'}] 
    output_names=set()
    for drug in input_drugs:
        name_id= meta_names[(meta_names['pert_id']==drug) | 
                (meta_names['alt_name']==drug) | 
                (meta_names['pert_iname']==drug)].index
        drug_names=set()
        for ind in name_id:
            names_list= list(dict.fromkeys(meta_names.loc[ind].dropna()))
            names_separated= set(list(chain.from_iterable(el.split("|") for el in names_list)))
            drug_names= drug_names | names_separated
        output_names= output_names | drug_names 
    output_drugs= list(set(input_drugs) | output_names)
    return output_drugs

def cmap_drug_id(input_drugs): #create DataFrame with pertubation ids and drug names from L1000 database
    df_meta= pd.read_csv('Datasets\\L1000FWD_drugs_metadata.csv')
    pert_lower=df_meta['pert_id'].str.lower()
    df_meta['pert_id']=pert_lower
    alt_lower=df_meta['alt_name'].str.lower()
    df_meta['alt_name']=alt_lower
    iname_lower=df_meta['pert_iname'].str.lower()
    df_meta['pert_iname']=iname_lower
    meta_names=df_meta[{'pert_id', 'alt_name', 'pert_iname'}] 
    drug_meta=pd.DataFrame()
    output_drugs= alt_names(input_drugs)
    for drug in output_drugs:
        name_id= meta_names[(meta_names['pert_id']==drug) | 
                 (meta_names['alt_name']==drug) | 
                 (meta_names['pert_iname']==drug)]
        drug_meta=pd.concat([drug_meta, name_id[{'pert_id', 'alt_name', 'pert_iname'}]])
    drug_meta.drop_duplicates(inplace=True)
    return drug_meta

def FScore(gene_list): #calculate fibrosis scores
    df1_s= ipf.loc[gene_list].transpose()
    df2_s= norm.loc[gene_list].transpose()
    df1_mean= df1_s.mean(axis=0)
    df2_mean= df2_s.mean(axis=0)
    score=0
    n_genes=len(gene_list)
    score_list= list()
    for gene in gene_list:
        if gene in ipf.index:
            stat, p= mannwhitneyu(df1_s[gene], df2_s[gene], alternative='two-sided')
            if p<0.05:
                if df1_mean.loc[gene]>df2_mean.loc[gene]:
                    score_list.append(1)
                    score+=1
                else:
                    score-=1
                    score_list.append(-1)
            else:
                score_list.append(0)
        else:
            score_list.append(0)
    df_out= pd.DataFrame({'Gene': gene_list, 'FScore': score_list})
    mean_score= score/n_genes
    return mean_score, df_out

def VScore(gene_list): #calculate viral scores
    score=0
    mean_score=0
    n_genes=len(gene_list)
    score_list= list()
    for gene in gene_list:
        if gene in vir_up:
            score+=1
            score_list.append(1)
        elif gene in vir_down:
            score-=1
            score_list.append(-1)
        else:
            score_list.append(0)
    df_out= pd.DataFrame({'Gene': gene_list, 'VScore': score_list})
    mean_score= score/n_genes
    return mean_score, df_out

def drug_signature(drugs, pert): #record viral and fibrosis drug signatures
    pert_id_list= list()
    df= pd.DataFrame()
    for drug in drugs:
        df_slice= df_A549[df_A549[pert]==drug]
        df_slice.set_index('sig_id', inplace=True)
        sig_ids= list(dict.fromkeys(df_slice.index))
        up_genes= list()
        down_genes= list()
        for ind in sig_ids:
            up_genes= up_genes+df_bsig.at[ind, 'up_genes']
            down_genes= down_genes+df_bsig.at[ind, 'down_genes']
            pert_id= df_slice.at[ind, 'pert_id']
        pert_id_list.append(pert_id)
        upregulated= list(set(up_genes))
        downregulated= list(set(down_genes))
        df_temp=pd.DataFrame({drug+'_up': pd.Series(upregulated), drug+'_down': pd.Series(downregulated)})
        df= pd.concat([df, df_temp], axis=1)
    return df, pert_id_list

def ViFi_score (drugs, df): #calculate Vi-Fi scores for list of drugs based on provided drug-induced signature (df)
    FS_list= list()
    VS_list= list()
    Sig_dict= {}
    for drug in drugs:
        drug_up= list(dict.fromkeys(df[drug+'_up'].dropna()))
        drug_down= list(dict.fromkeys(df[drug+'_down'].dropna()))
        FS_up, df_Fup= FScore(drug_up)
        FS_down, df_Fdown= FScore(drug_down)
        FS= FS_up- FS_down
        FS_list.append(FS)
        df_Fup['change']='up'
        df_Fdown['change']='down'
        df_Fsig=pd.concat([df_Fup, df_Fdown])
        VS_up, df_Vup= VScore(drug_up)
        VS_down, df_Vdown= VScore(drug_down)
        VS= VS_up- VS_down
        VS_list.append(VS)
        df_Vup['change']='up'
        df_Vdown['change']='down'
        df_Vsig=pd.concat([df_Vup, df_Vdown])
        Sig_dict[drug+'VSig']= df_Vsig
        Sig_dict[drug+'FSig']= df_Fsig
    return VS_list, FS_list , Sig_dict

"""define working directory"""
os.chdir('C:\\Vi-Fi scoring')
#Load DSigDB databases and user drug lists
input_file=input('Drug list file:')
df_user= pd.read_excel(input_file+'.xlsx')
df_bsig= pd.read_json('Datasets\\CD_signatures_binary_42809.json')
df_bsig.set_index('sig_id', inplace= True)
df_sigid= pd.read_csv('Datasets\\CD_signature_metadata.csv')
#Extract drug names and select signature ids for A549 from c-map 
user_drugs= list(dict.fromkeys(df_user['Drug'].str.lower()))
df_drugs_id= cmap_drug_id(user_drugs)
df_A549= df_sigid[df_sigid['cell_id']=='A549']
A549_idlower=df_A549['pert_id'].str.lower()
df_A549['pert_id']=A549_idlower
A549_namelower=df_A549['pert_desc'].str.lower()
df_A549['pert_desc']=A549_namelower
cmap_names= set(dict.fromkeys(df_A549['pert_desc']))
drug_names= list(set(user_drugs) & cmap_names)
cmapA549_ids= set(dict.fromkeys(df_A549['pert_id']))
user_ids= set(dict.fromkeys(df_drugs_id['pert_id'])) #select drugs present in c-map
drugs_id= list(user_ids & cmapA549_ids)
#Extract drug signatures from c-map into df_DS
df_DS_names, name_ids= drug_signature(drug_names, 'pert_desc')
df_DS_ids, pert_ids= drug_signature(drugs_id, 'pert_id')

#Load upregulated/downregulated genes used for fibrosis scoring 
ipf= pd.read_excel('Datasets\\IPF lung gene expression.xlsx')
norm= pd.read_excel('Datasets\\Healthy lung gene expression.xlsx')
norm.set_index('NAME', inplace=True)
ipf.set_index('NAME', inplace=True)
#Load upregulated/downregulated genes used for virus scoring 
df_A549vir= pd.read_excel('Datasets\\A549_vir_genes.xlsx')
df_NBHEvir= pd.read_excel('Datasets\\NBHE_vir_genes.xlsx')
A549_up= list(dict.fromkeys(df_A549vir['up'].dropna()))
A549_down= list(dict.fromkeys(df_A549vir['down'].dropna()))
NBHE_up= df_NBHEvir['up'].dropna()
NBHE_down= df_NBHEvir['down'].dropna()
up_set= set(A549_up) | set(NBHE_up)
down_set= set(A549_down) | set(NBHE_down)
vir_up= list(up_set-(up_set & down_set))
vir_down= list(down_set-(up_set & down_set))
#Calculate fibrosis and virus score and record to df_score
VS_names, FS_names, Sig_names= ViFi_score (drug_names, df_DS_names)
df_score_name= pd.DataFrame({'Drug_id': name_ids, 'Name': drug_names, 'alt_name': np.nan,
                        'Viral Score': VS_names, 'Fibrotic Score': FS_names, 'Source': 'pert_desc'})
VS_id, FS_id, Sig_ids= ViFi_score (drugs_id, df_DS_ids)
df_drugs_id.set_index('pert_id', inplace=True)
df_id_names= df_drugs_id.loc[drugs_id]
df_score_id= pd.DataFrame({'Drug_id': drugs_id, 'Name': df_id_names['pert_iname'], 'alt_name': df_id_names['alt_name'],
                        'Viral Score': VS_id, 'Fibrotic Score': FS_id, 'Source': 'pert_id'})
df_score=pd.concat([df_score_name, df_score_id])
#sort signatures based on highest difference and remove duplicates
df_score['score']= df_score['Viral Score']**2 + df_score['Fibrotic Score']**2
df_score.sort_values('score', ascending=False, inplace=True)
df_score.drop_duplicates(subset='Drug_id', inplace=True)
df_score.drop_duplicates(subset='Name', inplace=True)
df_score.set_index('Drug_id', inplace=True)
#extract viral and fibrosis drug signatures to Vi-Fi Signatures folder
df_score1= df_score[df_score['Source']=='pert_desc'] # for drugs identified by name
for drug in list(dict.fromkeys(df_score1['Name'].dropna())):
    df_Vtemp=Sig_names[drug+'VSig']
    df_Ftemp=Sig_names[drug+'FSig']
    df_Vtemp.to_excel('Vi-Fi Signatures\\'+drug+'_Viral_sig.xlsx')
    df_Ftemp.to_excel('Vi-Fi Signatures\\'+drug+'_Fibrosis_sig.xlsx')
df_score2= df_score[df_score['Source']=='pert_id'] # for drugs identified by pert_id
for drug in list(dict.fromkeys(df_score2.index)):
    df_Vtemp=Sig_ids[drug+'VSig']
    df_Ftemp=Sig_ids[drug+'FSig']
    drug_name= df_score2.at[drug, 'Name']
    df_Vtemp.to_excel('Vi-Fi Signatures\\'+drug_name+'_Viral_sig.xlsx')
    df_Ftemp.to_excel('Vi-Fi Signatures\\'+drug_name+'_Fibrosis_sig.xlsx')
df_out=df_score.drop(['score', 'Source'], axis=1)
df_out.to_excel(input_file+'_ViFi_scores.xlsx')
