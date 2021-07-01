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

def find_chembl(df):
    """Finds drugs for targets from df in ChEMBL database"""
    df_genes= to_chembl(df)
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
    return df_names

def find_dsigdb(df_input, DB): 
    """Finds drugs that target listed genes in DSigDB databases:
        kinases- data for kinase assays
        mining- data for text mining"""
    if DB=='kinases':
        df=pd.read_excel('DSigDB_kinases.xlsx')
    elif DB=='mining':
        df=pd.read_excel('DSigDB_mining.xlsx')
    else:
        return print('invalid DB')
    input_genes = set(dict.fromkeys(df_input['Gene']))
    DB_genes= set(dict.fromkeys(df['Gene']))
    genes= input_genes & DB_genes
    df_out= pd.DataFrame()
    for gene in genes:
        df1=df[df['Gene']==gene]
        drugs= set(dict.fromkeys(df1['Drug']))
        df_temp= pd.DataFrame(drugs)
        df_out= pd.concat([df_out, df_temp], axis=1)
    df_out.columns= genes
    return df_out
 
def FDA_drugs(df): #finds FDA approved drugs using from FDA_list_ChEMBL_lower.xlsx 
    df_FDA=pd.read_excel('FDA_list_ChEMBL.xlsx')
    approved= list()
    FDA_drugs= list(dict.fromkeys(df_FDA['drug_name']))
    for drug in list(dict.fromkeys(df['Drug'])):
        if drug in FDA_drugs:
            approved.append('yes')
        else:
            approved.append('no')
    df.insert(2,'FDA_approved', approved)
    return df

"""define working directory"""
#os.chdir('C:\\Vi-Fi scoring')
os.chdir('C:\\Lab\\Python\\Virus Paper\\Github\\')
#Part1: INPUT gene list of potential drug targets
file_name=input('Enter gene list file name in excel (for exapmle "Test_targets"):')
df_in= pd.read_excel(file_name+'.xlsx')
out_dir=os.path.join(os.getcwd(), 'Results')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
#find drugs in ChEMBL
df_chembl= find_chembl(df_in)
df_chembl.to_excel(os.path.join(out_dir, file_name+'_ChEMBL_Results.xlsx'))
#find drug targets in kinase inhibitors DSigDB databases and export to xlsx file
kinases_drugs= find_dsigdb(df_in, 'kinases')
kinases_counts= count_targets(kinases_drugs)
kinases_result= FDA_drugs(kinases_counts)
kinases_result.to_excel(os.path.join(out_dir,file_name+'_DSigDB_kinases.xlsx'))
#find drug targets in text mininng and computational DSigDB databases and export to xlsx file
mining_drugs= find_dsigdb(df_in, 'mining')
mining_counts= count_targets(mining_drugs)
mining_result= FDA_drugs(mining_counts)
mining_result.to_excel(os.path.join(out_dir,file_name+'_DSigDB_mining.xlsx'))
os.chdir('C:\\Lab\\Python\\Virus Paper\\Github\\')


