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
import os


def Find_drugs(df_input, DB): #finds drugs that target listed genes in DSigDB databases
    if DB=='kinases':
        df=pd.read_excel('Datasets\\DSigDB_kinases.xlsx')
    elif DB=='mining':
        df=pd.read_excel('Datasets\\DSigDB_mining.xlsx')
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
    
def FDA_drugs(df): #finds FDA approved drugs using from FDA_list_ChEMBL_lower.xlsx in Dataset folder  
    df_FDA=pd.read_excel('D:\\Python\\Datasets\\FDA_list_ChEMBL_lower.xlsx')
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
os.chdir('C:\\Vi-Fi scoring')
file_name=input('Gene list file:')
df_genes= pd.read_excel(file_name+'.xlsx')
#find drug targets in kinase inhibitors DSigDB databases and export to xlsx file
kinases_drugs= Find_drugs(df_genes, 'kinases')
kinases_counts= count_targets(kinases_drugs)
kinases_result= FDA_drugs(kinases_counts)
kinases_result.to_excel(file_name+'_DSigDB_kinases.xlsx')
#find drug targets in text mininng and computational DSigDB databases and export to xlsx file
mining_drugs= Find_drugs (df_genes, 'mining')
mining_counts= count_targets(mining_drugs)
mining_result= FDA_drugs(mining_counts)
mining_result.to_excel(file_name+'_DSigDB_mining.xlsx')