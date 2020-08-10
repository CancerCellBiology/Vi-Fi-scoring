# -*- coding: utf-8 -*-
"""
Code developed in the laboratory of Cancer Cell Biology,
Engelhardt Institute of Molecular Biology, Moscow, Russia.
Not for commercial use.
If you have any questions about the code or its use contact:
lebedevtd@gmai.com
vr.elmira@gmail.com 
"""

from chembl_webresource_client.new_client import new_client
import pandas as pd
import numpy as np
import urllib.parse
import urllib.request
import os

def to_chembl(df_genes): #converts UNIPROT ids to ChEMBL ids
    uniprot_list = list(dict.fromkeys(df_genes['Uniprot ID']))
    query = ' '.join([str(elem) for elem in uniprot_list]) 
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
    'from': 'ACC+ID',
    'to': 'CHEMBL_ID',
    'format': 'tab',
    'query': query
    }
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
       response = f.read()
    result= response.decode('utf-8')
    stringList = result.split('\n')
    df= pd.DataFrame([x.split('\t') for x in stringList]).drop(0)
    df.reset_index(drop=True, inplace=True)
    df.dropna(inplace=True)
    df.columns=['Uniprot_id', 'ChEMBLid']
    validgenes=list()
    df1= df_genes.set_index('Uniprot ID')
    for ind in df['Uniprot_id']:
        validgenes.append(df1.at[ind, 'Gene'])
    df['Gene']= validgenes
    return df

def count_targets(df): #counts target genes/proteins and creates a list of targets for each drug
    targets_df={}
    counts_df={}
    for target in list(dict.fromkeys(df.columns)):
        drugs= df[target].dropna().tolist()
        for drug in drugs:
            if drug not in targets_df.keys():
                targets_df[drug]=[target]
                counts_df[drug]=1
            else:
                targets_df[drug].append(target)
                counts_df[drug]+=1
    df_out=pd.DataFrame({'Drug': list(targets_df.keys()), 'Number_of_targets': list(counts_df.values()), 'Targets': list(targets_df.values())})
    return df_out

def name_ids(df): #adds drug names and max clinical trial phases for ChEMBL drug ids 
    name_list=list()
    phase_list=list()
    df.rename(columns={'Drug':'Drug_id'}, inplace=True)
    df.set_index('Drug_id', inplace=True)
    for ind in df.index: #records names and max_phases for each CHEMBL id 
        drug_stage= new_client.molecule.filter(molecule_chembl_id__exact=ind).only(['max_phase', 'pref_name'])
        for act in drug_stage: #cycles through one Django element to get single values from dict
            name_list.append(act['pref_name'])
            phase_list.append(act['max_phase'])
    df.insert(0,'Drug',name_list)
    df.insert(2,'max_phase',phase_list)
    return df

"""define working directory"""
os.chdir('C:\\Vi-Fi scoring')
#Part1: INPUT gene lists
file_name=input('Gene list file:')
df_input= pd.read_excel(file_name+'.xlsx')
df_genes= to_chembl(df_input)
gene_list = list(dict.fromkeys(df_genes['Gene']))
chEMBL_id = list(dict.fromkeys(df_genes['ChEMBLid']))
#Part2: create DataFrame (df_ch) with CHEMBL drug ids for each input gene 
df_ch= pd.DataFrame()
for id in chEMBL_id:
    drugs = new_client.mechanism.filter(target_chembl_id__exact=id).only(['molecule_chembl_id'])
    drug_list= list()
    for act in drugs:
        if act['molecule_chembl_id'] in drug_list: #check for duplication
            continue
        else:
            drug_list.append(act['molecule_chembl_id'])
    if len(drugs)==0: #record NaN values for genes with no drugs, to keep same df size
        df_temp=pd.DataFrame([np.NaN])
    else:
        df_temp=pd.DataFrame(drug_list)
    df_ch= pd.concat([df_ch, df_temp], axis=1)
df_ch.columns=gene_list
#Part3: create output DataFrame with drug targets and drug names
df_targets= count_targets(df_ch)
df_names= name_ids(df_targets)
df_names.to_excel(file_name+'_ChEMBL_Results.xlsx')